import numpy as np
import pandas as pd


def find_ap(signal, time, ap_threshold = -10):
    over_threshold = np.where(signal > ap_threshold)[0]

    if len(over_threshold) == 0:
        return []

    idx_diff = np.where(np.diff(over_threshold) != 1)[0]
    ap_count = np.append(over_threshold[0], over_threshold[idx_diff + 1])

    return time[ap_count]


samples = pd.read_csv(snakemake.input['samples'])

ap_data = []

for index, sample in samples.iterrows():
    print(f'Processing sample {sample["SID"]}')
    animal_id = sample['AnimalID']
    cell_name = sample['CellName']
    prefix = f'{snakemake.params["prefix"]}/{animal_id}/{cell_name}'

    sample_data = pd.read_pickle(f'{prefix}/vm/filter.pkl')

    sample_events = pd.read_csv(f'{prefix}/movement_final.csv')
    event_data = []
    for i, event in sample_events.iterrows():
        start = event['Start']
        end = event['End']
        channel = event['ChannelID']

        channel_data = sample_data.copy().crop(tmin = start, tmax = end)

        event_data.append(pd.DataFrame({
            'EventID': [i],
            'CountAP': len(
                find_ap(channel_data.get_data(picks = [channel]).flatten(), channel_data.times, ap_threshold = -20)),
            'Region': sample['Region'],
            'SID': sample['SID'],
            'Type': 'Movement'
        }))

    ap_data.append(pd.concat(event_data))

    sample_events = pd.read_csv(f'{prefix}/no_movement_final.csv')
    event_data = []
    for i, event in sample_events.iterrows():
        start = event['Start']
        end = event['End']
        channel = event['ChannelID']

        channel_data = sample_data.copy().crop(tmin = start, tmax = end)

        event_data.append(pd.DataFrame({
            'EventID': [i],
            'CountAP': len(
                find_ap(channel_data.get_data(picks = [channel]).flatten(), channel_data.times, ap_threshold = -20)),
            'Region': sample['Region'],
            'SID': sample['SID'],
            'Type': 'No movement'
        }))

    ap_data.append(pd.concat(event_data))

ap_data = pd.concat(ap_data, ignore_index = True)

ap_data.to_csv(snakemake.output['data'])
