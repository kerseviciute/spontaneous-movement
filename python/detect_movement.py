import pickle
import mne
import pandas as pd
import numpy as np

with open('.detect_movement.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)

# with open('.detect_movement.py.pkl', 'rb') as file:
#     snakemake = pickle.load(file)

with open(f'{snakemake.scriptdir}/methods.py', 'r') as file:
    exec(file.read())

data = read_pickle(snakemake.input['filter'])

tkeo_max_freq = snakemake.params['tkeoMaxFreq']

print('Calculating Teager-Kaiser Energy Operator for all available channels')
print(f'TKEO high frequency cutoff at {tkeo_max_freq} Hz')

tke = calculate_tkeo(data, h_freq = tkeo_max_freq)

threshold = snakemake.params['threshold']
min_length = snakemake.params['minLength']
min_break = snakemake.params['minBreak']
expand_by = snakemake.params['expandBy']

threshold_no_move = snakemake.params['threshold_no_move']
min_length_no_move = snakemake.params['minLength_no_move']
min_break_no_move = snakemake.params['minBreak_no_move']
expand_by_no_move = snakemake.params['expandBy_no_move']

print('Extracting events in all channels')
print(f'(Movement) TKEO threshold: {threshold}')
print(f'(Movement) Min event length: {min_length} seconds')
print(f'(Movement) Min break between events: {min_break} seconds')
print(f'(Movement) Events expanded by: {expand_by} seconds')

print(f'(No movement) TKEO threshold: {threshold_no_move}')
print(f'(No movement) Min event length: {min_length_no_move} seconds')
print(f'(No movement) Min break between events: {min_break_no_move} seconds')
print(f'(No movement) Events expanded by: {expand_by_no_move} seconds')

events_move = []
events_no_move = []

for channel in tke.ch_names:
    channel_data = tke.get_data(picks = [channel])[0]
    channel_events_move = extract_events(
        tkeo = channel_data,
        sampling_rate = tke.info['sfreq'],
        threshold = threshold,
        min_event_length = min_length,
        min_break = min_break,
        expand_by = expand_by,
        extract_movement = True
    )

    channel_events_move['Channel'] = channel
    events_move.append(channel_events_move)

    channel_events_no_move = extract_events(
        tkeo = channel_data,
        sampling_rate = tke.info['sfreq'],
        threshold = threshold_no_move,
        min_event_length = min_length_no_move,
        min_break = min_break_no_move,
        expand_by = expand_by_no_move,
        extract_movement = False
    )

    channel_events_no_move['Channel'] = channel
    events_no_move.append(channel_events_no_move)

#
# Preprocess no movement events
#

events_no_move = pd.concat(events_no_move, ignore_index = True)
events_no_move['EventId'] = [f'Event{x}' for x in range(len(events_no_move))]

print(f'Detected {len(events_no_move)} no movement events through all the channels')

for i, event in events_no_move.iterrows():
    start = event['Start']
    end = event['End']
    channel = event['Channel']

    # Calculate signal amplitude during the event
    event_data = data.copy().pick(picks = [channel]).crop(tmin = start, tmax = end)
    events_no_move.at[i, 'Amplitude'] = np.sqrt(np.mean(np.square(event_data.get_data()[0])))

save_pickle(events_no_move, snakemake.output['all_no_movement'])

#
# Preprocess movement events
#

events_move = pd.concat(events_move, ignore_index = True)
events_move['EventId'] = [f'Event{x}' for x in range(len(events_move))]

print(f'Detected {len(events_move)} movement events through all the channels')

min_movement_amplitude = snakemake.params['minMovementAmplitude']
min_amplitude_difference = snakemake.params['minAmplitudeDifference']
calm_time = snakemake.params['calmTime']

print(f'Filtering detected events')
print(f'Minimum required movement signal amplitude: {min_movement_amplitude}')
print(f'Minimum difference between movement and calm period: {min_amplitude_difference}')
print(f'Calm period length: {calm_time} seconds')

for i, event in events_move.iterrows():
    start = event['Start']
    end = event['End']
    channel = event['Channel']

    # Calculate signal amplitude during the event
    event_data = data.copy().pick(picks = [channel]).crop(tmin = start, tmax = end)
    events_move.at[i, 'Amplitude'] = np.sqrt(np.mean(np.square(event_data.get_data()[0])))

    # Calculate signal amplitude before the event
    before_start = start - calm_time if start - calm_time >= 0 else 0
    if before_start == start:
        events_move.at[i, 'AmplitudeBefore'] = None
    else:
        before_event = data.copy().pick(picks = [channel]).crop(tmin = before_start, tmax = start)
        events_move.at[i, 'AmplitudeBefore'] = np.sqrt(np.mean(np.square(before_event.get_data()[0])))

    # Calculate signal amplitude after the event
    after_end = end + calm_time if end + calm_time <= np.max(data.times) else np.max(data.times)
    if after_end == end:
        events_move.at[i, 'AmplitudeAfter'] = None
    else:
        after_event = data.copy().pick(picks = [channel]).crop(tmin = end, tmax = after_end)
        events_move.at[i, 'AmplitudeAfter'] = np.sqrt(np.mean(np.square(after_event.get_data()[0])))


def mean_amplitude_diff(row):
    if np.isnan(row['AmplitudeBefore']):
        return row['Amplitude'] - row['AmplitudeAfter']

    if np.isnan(row['AmplitudeAfter']):
        return row['Amplitude'] - row['AmplitudeBefore']

    return (2 * row['Amplitude'] - row['AmplitudeAfter'] - row['AmplitudeBefore']) / 2


events_move['MeanAmplitudeDiff'] = events_move.apply(mean_amplitude_diff, axis = 1)

save_pickle(events_move, snakemake.output['all_movement'])

# Filter based on amplitude
events_move = events_move[events_move['Amplitude'] >= min_movement_amplitude]
events_move = events_move[events_move['MeanAmplitudeDiff'] >= min_amplitude_difference]

# Drop events that start at the very beginning or end at the very end of the recording
events_move = events_move[events_move['Start'] != 0.0]
events_move = events_move[events_move['End'] != np.max(data.times)]

events_move = events_move.reset_index(drop = True)

print(f'Events after filtering: {len(events_move)}')

save_pickle(events_move, snakemake.output['quality_movement'])

amplitude_diff = snakemake.config['movement']['minAmplitudeDiff']
min_move_amplitude = snakemake.config['movement']['minMovementAmplitude']
min_event_length = snakemake.config['movement']['minEventLength']
window_size = snakemake.config['movement']['smoothWindowSize']
resample_factor = snakemake.config['movement']['resampleFactor']

rest_amplitude = np.mean(events_no_move['Amplitude'])

exact_time_events = []

for i, event in events_move.iterrows():
    exact_event = determine_exact_time(
        event = event,
        data = data,
        rest_amplitude = rest_amplitude,
        amplitude_diff = amplitude_diff,
        min_move_amplitude = min_move_amplitude,
        min_event_length = min_event_length,
        window_size = window_size,
        resample_factor = resample_factor,
        verbose = False
    )

    exact_event['Channel'] = event['Channel']
    exact_time_events.append(exact_event)

exact_time_events = pd.concat(exact_time_events, ignore_index = True)

for i, event in exact_time_events.iterrows():
    start = event['Start']
    end = event['End']
    channel = event['Channel']

    # Calculate signal amplitude during the event
    event_data = data.copy().pick(picks = [channel]).crop(tmin = start, tmax = end)
    exact_time_events.at[i, 'Amplitude'] = np.sqrt(np.mean(np.square(event_data.get_data()[0])))

print(f'Events after determining exact onset/offset times: {len(exact_time_events)}')

save_pickle(exact_time_events, snakemake.output['final_movement'])
