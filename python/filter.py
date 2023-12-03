import pickle
import mne

with open('.filter.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)

# with open('.filter.py.pkl', 'rb') as file:
#     snakemake = pickle.load(file)

with open(f'{snakemake.scriptdir}/methods.py', 'r') as file:
    exec(file.read())

raw = read_pickle(snakemake.input['raw'])

l_freq = snakemake.params['drop_below']
h_freq = snakemake.params['drop_above']

raw.filter(l_freq = l_freq, h_freq = h_freq, fir_design = 'firwin', picks = 'emg')

save_pickle(raw, snakemake.output['filter'])