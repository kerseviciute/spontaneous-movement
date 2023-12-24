def read_pickle(filename):
    import pickle

    with open(filename, 'rb') as file:
        data = pickle.load(file)

    return data


def save_pickle(data, filename):
    import pickle

    with open(filename, 'wb') as file:
        pickle.dump(data, file)


def tkeo(signal):
    return signal[1:-1] ** 2 - signal[:-2] * signal[2:]


def calculate_tkeo(data, h_freq = 20):
    import numpy as np
    import mne

    channels = data.ch_names
    sfreq = data.info['sfreq']
    tke = []

    for i, channel in enumerate(channels):
        signal = data.get_data(picks = channel)[0]
        tke.append(tkeo(signal))

    tke = np.array(tke)

    ch_types = ['emg' for i in range(0, len(channels))]
    info = mne.create_info(ch_names = channels, sfreq = sfreq, ch_types = ch_types)
    tke = mne.io.RawArray(tke, info, verbose = False)

    tke.filter(l_freq = None, h_freq = h_freq, fir_design = 'firwin', picks = 'emg', verbose = False)

    return tke


def extract_events(tkeo, sampling_rate, threshold = 0.02, min_event_length = 0.1, min_break = 0.05,
                   expand_by: float = 0.0, extract_movement = True):
    """
    Extracts events by analysing a single channel TKEO.

    Parameters:
    - tkeo (array): single channel TKEO signal.
    - sampling_rate (int): data sampling rate.
    - threshold (float): signal threshold for movement detection.
    - min_event_length (float): minimum event length, in seconds.
    - min_break (float): minimum time period between two movement events, in seconds. If the movements are separated by a shorter calm phase, they will be joined into a single event.
    - expand_by (float): movement start and end times will be expanded by specified time period, in seconds. Note, that overlapping events after time expansion are not merged.
    - extract_movement (bool): whether movement or no movement events should be extracted.

    Returns:
    A pandas DataFrame with the movement periods.
    """
    import numpy as np
    import pandas as pd

    # Analyse where the TKEO signal is above specified threshold
    signal_over_threshold = tkeo >= threshold
    change_indices = np.where(np.diff(signal_over_threshold))[0]

    movement_data = pd.DataFrame({
        'EventStart': np.insert(change_indices + 1, 0, 0),
        'EventEnd': np.append(change_indices, len(signal_over_threshold))
    })

    movement_data['Movement'] = signal_over_threshold[movement_data['EventStart']]

    # Calculate the movement time in points and merge events
    movement_data['EventLength'] = list(movement_data['EventEnd'] - movement_data['EventStart'])
    movement_data = merge_short(movement_data, min_break = min_break * sampling_rate)

    # Filter by minimum time
    min_points = sampling_rate * min_event_length
    movement_data = movement_data[movement_data['EventLength'] >= min_points]

    # Select only movement data
    if extract_movement:
        movement_data = movement_data[movement_data['Movement']]
    else:
        movement_data = movement_data[~movement_data['Movement']]

    # Add expansion
    movement_data['EventStart'] = movement_data['EventStart'] - sampling_rate * expand_by
    movement_data['EventStart'] = movement_data['EventStart'].apply(lambda x: 0 if x < 0 else x)

    movement_data['EventEnd'] = movement_data['EventEnd'] + sampling_rate * expand_by
    movement_data['EventEnd'] = movement_data['EventEnd'].apply(lambda x: x if x < len(tkeo) else len(tkeo))

    movement_data['EventLength'] = movement_data['EventEnd'] - movement_data['EventStart']

    # Convert to seconds
    movement_data['Start'] = movement_data['EventStart'] / sampling_rate
    movement_data['End'] = movement_data['EventEnd'] / sampling_rate
    movement_data['Length'] = movement_data['EventLength'] / sampling_rate

    return movement_data


def merge_short(data, min_break = 200):
    i = 1
    while i < data.shape[0] - 1:
        if data.iloc[i]['EventLength'] >= min_break:
            i = i + 1
            continue

        previous_el = data.iloc[i - 1]
        next_el = data.iloc[i + 1]
        current = data.iloc[i]

        data.at[i, 'EventLength'] = previous_el['EventLength'] + next_el['EventLength'] + current['EventLength']
        data.at[i, 'EventStart'] = previous_el['EventStart']
        data.at[i, 'EventEnd'] = next_el['EventEnd']
        data.at[i, 'Movement'] = previous_el['Movement']

        data = data.drop(i - 1)
        data = data.drop(i + 1)
        data = data.reset_index(drop = True)

    return data


def determine_exact_time(event, data, rest_amplitude, amplitude_diff = 1.5, min_move_amplitude = 0.25,
                         min_event_length = 0.1, window_size = 200, resample_factor = 4, verbose = False):
    import numpy as np
    import pandas as pd

    sfreq = data.info['sfreq']

    # Get event data
    channel = event['Channel']
    event_start = np.max([event['Start'] - 0.25, 0])
    event_end = np.min([event['End'] + 0.25, np.max(data.times)])

    event_data = data.copy().pick(picks = [channel]).crop(tmin = event_start, tmax = event_end)
    y = event_data.get_data()[0]

    # Resample the data
    y = y[::resample_factor]
    new_sfreq = int(sfreq / resample_factor)

    if verbose: print(f'Obtained sampling rate: {new_sfreq} Hz')

    # Smooth the signal by rolling average
    smooth_y = np.convolve(np.abs(y), np.ones(window_size) / window_size, mode = 'same')
    smooth_y = smooth_y - amplitude_diff * rest_amplitude

    # Find where the movement happens
    signal_over_threshold = smooth_y >= 0
    change_indices = np.where(np.diff(signal_over_threshold))[0]

    movement_data = pd.DataFrame({
        'EventStart': np.insert(change_indices + 1, 0, 0),
        'EventEnd': np.append(change_indices, len(signal_over_threshold))
    })

    movement_data['Movement'] = signal_over_threshold[movement_data['EventStart']]

    # Calculate the movement time in points
    movement_data['EventLength'] = list(movement_data['EventEnd'] - movement_data['EventStart'])
    movement_data = movement_data[movement_data['Movement']]

    # Filter events by minimum time length
    min_event_points = new_sfreq * min_event_length
    movement_data = movement_data[movement_data['EventLength'] > min_event_points]
    movement_data = movement_data.reset_index(drop = True)

    if verbose: print(f'Number of estimated movements: {len(movement_data)}')

    # Calculate the signal amplitude
    if verbose: print('Calculating the average signal amplitude during estimated movements')
    for index, movement in movement_data.iterrows():
        start = movement['EventStart']
        end = movement['EventEnd']
        dt = smooth_y[start:end]
        movement_data.at[index, 'Amplitude'] = np.sqrt(np.mean(np.square(dt)))

    # Filter by amplitude
    movement_data = movement_data[movement_data['Amplitude'] > min_move_amplitude]
    movement_data = movement_data.reset_index(drop = True)

    # Recalculate the event start and end times as they are relative to the event start
    movement_data['Start'] = movement_data['EventStart'] / new_sfreq + event_start
    movement_data['Start'] = movement_data['Start'].apply(lambda x: 0 if x < 0 else x)

    movement_data['End'] = movement_data['EventEnd'] / new_sfreq + event_start
    movement_data['End'] = movement_data['End'].apply(lambda x: x if x < np.max(data.times) else np.max(data.times))

    movement_data['Length'] = list(movement_data['End'] - movement_data['Start'])

    movement_data['EventStart'] = movement_data['Start'] * sfreq
    movement_data['EventEnd'] = movement_data['End'] * sfreq
    movement_data['EventLength'] = list(movement_data['EventEnd'] - movement_data['EventStart'])

    if verbose: print(f'Number of estimated movements after amplitude filtering: {len(movement_data)}')

    return movement_data


def signal_amplitude(signal):
    import numpy as np
    return np.sqrt(np.mean(np.square(signal)))


def event_amplitude(event, data):
    start = event['Start']
    end = event['End']
    channel = event['Channel']

    # Calculate signal amplitude during the event
    event_data = data.copy().pick(picks = [channel]).crop(tmin = start, tmax = end)
    return signal_amplitude(event_data.get_data()[0])


def amplitude_all_events(events, data):
    amplitudes = []
    for i, event in events.iterrows():
        amplitudes.append(event_amplitude(event, data))

    return amplitudes
