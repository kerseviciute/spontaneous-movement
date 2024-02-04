import numpy as np
import pandas as pd

pd.to_pickle(snakemake, '.ap_time_start.py.pkl')
# snakemake = pd.read_pickle('.ap_time_start.py.pkl')


def find_over_threshold(signal, ap_threshold = -10):
    over_threshold = np.where(signal > ap_threshold)[0]

    if len(over_threshold) == 0:
        return []

    idx_diff = np.where(np.diff(over_threshold) != 1)[0]
    ap_count = np.append(over_threshold[0], over_threshold[idx_diff + 1])

    return ap_count


def find_ap(signal, time, channel, sid, ap_threshold = -20):
    differential = np.diff(signal)
    aps = find_over_threshold(signal, ap_threshold = ap_threshold)

    if len(aps) == 0:
        return pd.DataFrame()

    action_potentials = []

    for ap_index in aps:
        start_idx = np.max([ap_index - 50, 0])
        end_idx = np.min([ap_index + 50, len(differential)])
        diff_over_threshold = np.where(differential[start_idx:end_idx] > 1)[0]

        if len(diff_over_threshold) == 0:
            max_val = np.max(differential[start_idx:end_idx])
            print(f'Skipping AP on {time[ap_index]} with voltage {signal[ap_index]} with max diff {max_val}')
            continue

        exact_onset = diff_over_threshold[0] + ap_index - 50

        action_potentials.append(pd.DataFrame({
            'Trial': [f'{channel}_{sid}'],
            'Time': [time[exact_onset]],
        }))

    action_potentials = pd.concat(action_potentials)
    action_potentials = action_potentials.reset_index(drop = True)

    return action_potentials


samples = pd.read_csv(snakemake.input['samples'])

action_potentials = []

for index, sample in samples.iterrows():
    print(f'Processing sample {sample["SID"]}')
    animal_id = sample['AnimalID']
    cell_name = sample['CellName']

    prefix = f'{snakemake.params["prefix"]}/{animal_id}/{cell_name}'

    movement_events = pd.read_csv(f'{prefix}/movement_start_final.csv')
    data = pd.read_pickle(f'{prefix}/vm/filter.pkl')

    for i, event in movement_events.iterrows():
        event_data = data.copy().pick(picks = event['ChannelID']).crop(tmin = event['Start'] - 0.4, tmax = event['Start'] + 0.4)

        time = event_data.times
        signal = event_data.get_data()[0, :]

        event_aps = find_ap(signal, time, channel = i, sid = sample['SID'], ap_threshold = -15)
        action_potentials.append(event_aps)

action_potentials = pd.concat(action_potentials)
action_potentials = action_potentials.reset_index(drop = True)

action_potentials.to_csv(snakemake.output['data'])
