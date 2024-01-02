import pandas as pd
import numpy as np
import mne

pd.to_pickle(snakemake, '.extract_data.py.pkl')
# snakemake = pd.read_pickle('.extract_data.py.pkl')

with open(f'{snakemake.scriptdir}/../python/methods.py', 'r') as file:
    exec(file.read())

data = pd.read_pickle(snakemake.input['data'])
vm = pd.read_pickle(snakemake.input['vm'])
tke = calculate_tkeo(data, h_freq = 20)

info = []

for index, channel in enumerate(data.ch_names):
    ch_data = data.get_data(picks = [channel])[0]
    vm_data = vm.get_data(picks = [index])[0]
    tkeo_data = tke.get_data(picks = [index])[0]
    tkeo_data = np.insert(tkeo_data, [0 , len(tkeo_data)], [tkeo_data[0], tkeo_data[-1]])

    df = pd.DataFrame({
        'EMG': ch_data,
        'Time': data.times,
        'Channel': channel
    })

    info.append(pd.DataFrame({
        'VM': vm_data,
        'EMG': ch_data,
        'TKEO': tkeo_data,
        'Time': data.times,
        'Channel': channel,
        'ChannelID': index
    }))

info = pd.concat(info)

info.to_csv(snakemake.output['data'])
