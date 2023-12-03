import pickle
import numpy as np
import csv
import mne

with open('.create_mne.py.pkl', 'wb') as file:
    pickle.dump(snakemake, file)

# with open('.create_mne.py.pkl', 'rb') as file:
#     snakemake = pickle.load(file)

with open(f'{snakemake.scriptdir}/methods.py', 'r') as file:
    exec(file.read())

# Read raw data
with open(snakemake.input['raw'], 'r') as file:
    data = np.array(list(csv.reader(file, delimiter = '\t')))

sfreq = snakemake.params['samplingRate']

# Create MNE object
channels = list(data[0, :])
channels = [channel.replace("'", "") for channel in channels]
ch_types = ['emg' for i in range(0, len(channels))]

data = data[1:, :].transpose()

info = mne.create_info(ch_names = channels, sfreq = sfreq, ch_types = ch_types)
raw = mne.io.RawArray(data, info)

save_pickle(raw, snakemake.output['raw'])
