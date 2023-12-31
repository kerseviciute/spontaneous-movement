import mne
import pandas as pd

pd.to_pickle(snakemake, '.extract_events.py.pkl')
# snakemake = pd.read_pickle('.extract_events.py.pkl')

with open(f'{snakemake.scriptdir}/methods.py', 'r') as file:
    exec(file.read())

data = pd.read_pickle(snakemake.input['filter'])

tkeo_max_freq = snakemake.params['tkeoMaxFreq']

print('Calculating Teager-Kaiser Energy Operator for all available channels')
print(f'TKEO high frequency cutoff at {tkeo_max_freq} Hz')

tke = calculate_tkeo(data, h_freq = tkeo_max_freq)

threshold = snakemake.params['threshold']
min_length = snakemake.params['minLength']
min_break = snakemake.params['minBreak']
expand_by = snakemake.params['expandBy']
detect_movement = snakemake.params['detectMovement']

print('Extracting events in all channels')
print(f'Detecting movement events: {detect_movement}')
print(f'TKEO threshold: {threshold}')
print(f'Min event length: {min_length} seconds')
print(f'Min break between events: {min_break} seconds')
print(f'Events expanded by: {expand_by} seconds')

events = []

for i, channel in enumerate(tke.ch_names):
    channel_data = tke.get_data(picks = [channel])[0]
    channel_events = extract_events(
        tkeo = channel_data,
        sampling_rate = tke.info['sfreq'],
        threshold = threshold,
        min_event_length = min_length,
        min_break = min_break,
        expand_by = expand_by,
        extract_movement = detect_movement
    )

    channel_events['Channel'] = channel
    channel_events['ChannelID'] = i
    events.append(channel_events)

events = pd.concat(events, ignore_index = True)
events['EventId'] = [f'Event{x}' for x in range(len(events))]
events['Amplitude'] = amplitude_all_events(events, data)

print(f'Detected {len(events)} events through all the channels')

events.to_csv(snakemake.output['events'])
