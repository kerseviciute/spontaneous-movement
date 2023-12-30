import pickle
import mne
import pandas as pd
import numpy as np

with open('.filter_events.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)

# with open('.filter_events.py.pkl', 'rb') as file:
#     snakemake = pickle.load(file)

with open(f'{snakemake.scriptdir}/methods.py', 'r') as file:
    exec(file.read())

data = read_pickle(snakemake.input['filter'])
movement = read_pickle(snakemake.input['movement'])
no_movement = read_pickle(snakemake.input['no_movement'])


print('1. Filter by amplitude (compare with maximum no movement amplitude)')

amplitude_diff_times = snakemake.params['amplitude_diff_times']

amplitude_cutoff = np.max(no_movement['Amplitude']) * amplitude_diff_times
print(f'   Number of events before filtering: {len(movement)}')
movement = movement[ movement['Amplitude'] >= amplitude_cutoff ]
print(f'   Number of events after filtering: {len(movement)}')


# print('2. Drop events that are at the very beginning / very end of the recording')

calm_time = snakemake.params['calm_time']

# print(f'   Number of events before filtering: {len(movement)}')
# movement = movement[ movement['Start'] > calm_time ]
# movement = movement[ movement['End'] < np.max(data.times) - calm_time ]
# print(f'   Number of events after filtering: {len(movement)}')


print('3. Determine exact movement onset/offset time')

window_size = snakemake.params['window_size']
resample_factor = snakemake.params['resample_factor']

rest_amplitude = np.mean(no_movement['Amplitude'])
print(f'   Mean resting amplitude {rest_amplitude}')
# Note: using mean resting amplitude for EVENT DETECTION (less strict)
#       and maximum resting amplitude for EVENT FILTERING (more strict)

exact_time_events = []
for i, event in movement.iterrows():
    exact_event = determine_exact_time(
        event = event,
        data = data,
        rest_amplitude = rest_amplitude,
        expand_by = calm_time,
        window_size = window_size,
        resample_factor = resample_factor
    )

    exact_event['Channel'] = event['Channel']
    exact_event['ChannelId'] = event['ChannelId']
    exact_time_events.append(exact_event)

exact_time_events = pd.concat(exact_time_events, ignore_index = True)

print(f'   Detected events: {len(exact_time_events)}')


print('4. Filtering newly detected events')

max_merge_time = snakemake.params['max_merge_time']
calm_amplitude_difference = snakemake.params['calm_amplitude_difference']
min_movement_time = snakemake.params['min_movement_time']

print('   Filter by amplitude')

exact_time_events['Amplitude'] = amplitude_all_events(exact_time_events, data)

print(f'       Number of events before filtering: {len(exact_time_events)}')
exact_time_events = exact_time_events[exact_time_events['Amplitude'] > amplitude_cutoff]
print(f'       Number of events after filtering: {len(exact_time_events)}')


print(f'   Merge events that are less than {max_merge_time} seconds apart')

# TODO: move this into a function

final_events = []

for channel in np.unique(exact_time_events['Channel']):
    dt = exact_time_events[exact_time_events['Channel'] == channel]
    dt = dt.sort_values(by = 'Start')
    dt = dt.reset_index(drop = True)

    i = 1
    while i < len(dt):
        previous_element = dt.iloc[i - 1]
        current_element = dt.iloc[i]

        if previous_element['Channel'] != current_element['Channel']:
            i = i + 1
            continue

        if previous_element['End'] + max_merge_time < current_element['Start']:
            i = i + 1
            continue

        dt.at[i, 'EventStart'] = previous_element['EventStart']
        dt.at[i, 'EventEnd'] = current_element['EventEnd']

        dt = dt.drop(i - 1)
        dt = dt.reset_index(drop = True)

    final_events.append(dt)

final_events = pd.concat(final_events)
final_events = final_events.reset_index(drop = True)

# Recalculate event parameters
final_events['EventLength'] = list(final_events['EventEnd'] - final_events['EventStart'])
final_events['Start'] = final_events['EventStart'] / data.info['sfreq']
final_events['End'] = final_events['EventEnd'] / data.info['sfreq']
final_events['Length'] = list(final_events['End'] - final_events['Start'])
final_events['Amplitude'] = amplitude_all_events(final_events, data)

print(f'       Number of events after merging: {len(final_events)}')

print(f'   Calculating average amplitude before and after the events')
final_events['AmplitudeBefore'] = around_amplitude_all_events(final_events, data, calm_time = calm_time, before_event = True)
final_events['AmplitudeAfter'] = around_amplitude_all_events(final_events, data, calm_time = calm_time, before_event = False)

cutoff = np.mean(no_movement['Amplitude']) * calm_amplitude_difference
print(f'       Using before/after event amplitude cutoff of {round(cutoff, 3)} ({calm_amplitude_difference} * mean resting amplitude)')
final_events = final_events[ final_events['AmplitudeBefore'] < cutoff ]
final_events = final_events[ final_events['AmplitudeAfter'] < cutoff ]

print(f'       Number of events after filtering: {len(final_events)}')

print(f'   Filtering events by minimum time: {min_movement_time} seconds')
final_events = final_events[ final_events['Length'] >= min_movement_time ]

print(f'       Number of events after filtering: {len(final_events)}')

print('2. Drop events that are at the very beginning of the recording')

print(f'   Number of events before filtering: {len(movement)}')
final_events = final_events[ final_events['Start'] > calm_time ]
print(f'   Number of events after filtering: {len(movement)}')

print(f'Final number of events after onset/offset time correction and filtering: {len(final_events)}')

final_events = final_events.reset_index(drop = True)

save_pickle(final_events, snakemake.output['movement'])
