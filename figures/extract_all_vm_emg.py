import pandas as pd
import re

pd.to_pickle(snakemake, '.extract_all_vm_emg.py.pkl')
# snakemake = pd.read_pickle('.extract_average_vm_emg.py.pkl')

with open(f'{snakemake.scriptdir}/../python/methods.py', 'r') as file:
    exec(file.read())

samples = pd.read_csv(snakemake.input['samples'])
event_data = []

for sid in samples['SID']:
    print(f'Processing {sid}')
    animal_id = re.match('(W[0-9]+).*', sid).group(1)
    cell_name = re.match('.*(C[0-9]+)', sid).group(1)

    prefix = f'{snakemake.params["prefix"]}/{animal_id}/{cell_name}'

    movement_events = pd.read_pickle(f'{prefix}/emg/filtered_movement_events.pkl')
    data = pd.read_pickle(f'{prefix}/emg/filter.pkl')
    vm = pd.read_pickle(f'{prefix}/vm/filter.pkl')
    tke = calculate_tkeo(data, h_freq = 20)

    for index, event in movement_events.iterrows():
        emg_data = data.copy().crop(tmin = event['Start'] - 0.4, tmax = event['Start'] + 0.4)
        x = emg_data.times
        emg_y = emg_data.get_data(picks = [event['Channel']])[0]

        tke_data = tke.copy().crop(tmin = event['Start'] - 0.4, tmax = event['Start'] + 0.4)
        tke_y = tke_data.get_data(picks = [event['Channel']])[0]

        vm_data = vm.copy().crop(tmin = event['Start'] - 0.4, tmax = event['Start'] + 0.4)
        vm_y = vm_data.get_data(picks = [event['ChannelId']])[0]

        region = samples[samples['SID'] == sid].iloc[0]['Region']

        event_data.append(pd.DataFrame({
            'Time': x,
            'EMG': emg_y,
            'VM': vm_y,
            'TKEO': tke_y,
            'SID': sid,
            'Region': region,
            'EventId': index,
            'EventStart': 0.4
        }))

event_data = pd.concat(event_data)

event_data.to_csv(snakemake.output['data'], index = False)
