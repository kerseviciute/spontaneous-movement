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

print('Extracting events in all channels')
print(f'TKEO threshold: {threshold}')
print(f'Min event length: {min_length} seconds')
print(f'Min break between events: {min_break} seconds')
print(f'Events expanded by: {expand_by} seconds')

events = []

for channel in tke.ch_names:
    channel_events = extract_events(
        tkeo = tke.get_data(picks = [channel])[0],
        sampling_rate = tke.info['sfreq'],
        threshold = threshold,
        min_event_length = min_length,
        min_break = min_break,
        expand_by = expand_by
    )

    channel_events['Channel'] = channel

    events.append(channel_events)

events = pd.concat(events, ignore_index = True)

print(f'Detected {len(events)} movements through all the channels')

min_movement_amplitude = snakemake.params['minMovementAmplitude']
min_amplitude_difference = snakemake.params['minAmplitudeDifference']
calm_time = snakemake.params['calmTime']

print(f'Filtering detected events')
print(f'Minimum required movement signal amplitude: {min_movement_amplitude}')
print(f'Minimum difference between movement and calm period: {min_amplitude_difference}')
print(f'Calm period length: {calm_time} seconds')

for i, event in events.iterrows():
    start = event['Start']
    end = event['End']
    channel = event['Channel']

    # Calculate signal amplitude during the event
    event_data = data.copy().pick(picks = [channel]).crop(tmin = start, tmax = end)
    events.at[i, 'Amplitude'] = np.sqrt(np.mean(np.square(event_data.get_data()[0])))

    # Calculate signal amplitude before the event
    before_start = start - calm_time if start - calm_time >= 0 else 0
    if before_start == start:
        events.at[i, 'AmplitudeBefore'] = None
    else:
        before_event = data.copy().pick(picks = [channel]).crop(tmin = before_start, tmax = start)
        events.at[i, 'AmplitudeBefore'] = np.sqrt(np.mean(np.square(before_event.get_data()[0])))

    # Calculate signal amplitude after the event
    after_end = end + calm_time if end + calm_time <= np.max(data.times) else np.max(data.times)
    if after_end == end:
        events.at[i, 'AmplitudeAfter'] = None
    else:
        after_event = data.copy().pick(picks = [channel]).crop(tmin = end, tmax = after_end)
        events.at[i, 'AmplitudeAfter'] = np.sqrt(np.mean(np.square(after_event.get_data()[0])))


def mean_amplitude_diff(row):
    if np.isnan(row['AmplitudeBefore']):
        return row['Amplitude'] - row['AmplitudeAfter']

    if np.isnan(row['AmplitudeAfter']):
        return row['Amplitude'] - row['AmplitudeBefore']

    return (2 * row['Amplitude'] - row['AmplitudeAfter'] - row['AmplitudeBefore']) / 2


events['MeanAmplitudeDiff'] = events.apply(mean_amplitude_diff, axis = 1)

save_pickle(events, snakemake.output['all'])

# Filter based on amplitude
events = events[events['Amplitude'] >= min_movement_amplitude]
events = events[events['MeanAmplitudeDiff'] >= min_amplitude_difference]

# Drop events that start at the very beginning or end at the very end of the recording
events = events[events['Start'] != 0.0]
events = events[events['End'] != np.max(data.times)]

events = events.reset_index(drop = True)

print(f'Events after filtering: {len(events)}')

# TODO: afterwards, cut the expanded times for more precise movement detection

save_pickle(events, snakemake.output['movement'])
