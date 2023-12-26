import pickle
import mne
import numpy as np

with open('.filter.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)

# with open('.filter.py.pkl', 'rb') as file:
#     snakemake = pickle.load(file)

with open(f'{snakemake.scriptdir}/methods.py', 'r') as file:
    exec(file.read())

raw = read_pickle(snakemake.input['raw'])

l_freq = snakemake.params['drop_below']
h_freq = snakemake.params['drop_above']

ch_type = snakemake.params['ch_type']
raw.filter(l_freq = l_freq, h_freq = h_freq, fir_design = 'firwin', picks = ch_type)

if snakemake.params['scale']:
    print('Scaling the channel data')
    data = raw.get_data().transpose()
    mean = np.mean(data)
    std = np.std(data)
    data = (data - mean) / std
    raw._data = data.transpose()
else:
    print('Channels were not scaled')

save_pickle(raw, snakemake.output['filter'])
