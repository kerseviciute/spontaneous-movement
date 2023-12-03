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
                   expand_by: float = 0.0):
    """
    Extracts events by analysing a single channel TKEO.

    Parameters:
    - tkeo (array): single channel TKEO signal.
    - sampling_rate (int): data sampling rate.
    - threshold (float): signal threshold for movement detection.
    - min_event_length (float): minimum event length, in seconds.
    - min_break (float): minimum time period between two movement events, in seconds. If the movements are separated by a shorter calm phase, they will be joined into a single event.
    - expand_by (float): movement start and end times will be expanded by specified time period, in seconds. Note, that overlapping events after time expansion are not merged.

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
    movement_data = movement_data[movement_data['Movement']]

    # Add expansion
    movement_data['EventStart'] = movement_data['EventStart'] - sampling_rate * expand_by
    movement_data['EventEnd'] = movement_data['EventEnd'] + sampling_rate * expand_by
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
