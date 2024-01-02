import pandas as pd

samples = pd.read_csv(snakemake.input['samples'])

final_data = []
for index, sample in samples.iterrows():
    print(f'Processing sample {sample["SID"]}')
    animal_id = sample['AnimalID']
    cell_name = sample['CellName']
    prefix = f'{snakemake.params["prefix"]}/{animal_id}/{cell_name}'

    sample_data = pd.read_pickle(f'{prefix}/vm/filter.pkl')

    sample_events = pd.read_csv(f'{prefix}/events2.csv')
    event_data = []
    for i, event in sample_events.iterrows():
        start = event['Start']
        end = event['End']
        channel = event['ChannelID']

        channel_data = sample_data.copy().crop(tmin = start, tmax = end)

        event_data.append(pd.DataFrame({
            'EventID': i,
            'Time': channel_data.times,
            'Vm': channel_data.get_data(picks = [channel]).flatten(),
            'Region': sample['Region'],
            'SID': sample['SID'],
            'Type': 'Movement'
        }))

    final_data.append(pd.concat(event_data))

    sample_events = pd.read_csv(f'{prefix}/no_movement_events.csv')
    event_data = []
    for i, event in sample_events.iterrows():
        start = event['Start']
        end = event['End']
        channel = event['ChannelID']

        channel_data = sample_data.copy().crop(tmin = start, tmax = end)

        event_data.append(pd.DataFrame({
            'EventID': i,
            'Time': channel_data.times,
            'Vm': channel_data.get_data(picks = [channel]).flatten(),
            'Region': sample['Region'],
            'SID': sample['SID'],
            'Type': 'No movement'
        }))

    final_data.append(pd.concat(event_data))

final_data = pd.concat(final_data, ignore_index = True)
final_data.to_csv(snakemake.output['data'])
