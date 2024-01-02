import mne
import pandas as pd

pd.to_pickle(snakemake, '.create_mne.py.pkl')
# snakemake = pd.read_pickle('.create_mne.py.pkl')

with open(f'{snakemake.scriptdir}/methods.py', 'r') as file:
    exec(file.read())

# Read raw data
data = pd.read_csv(snakemake.input['raw'], sep = '\t', index_col = False)

sfreq = snakemake.params['samplingRate']

# Create MNE object
channels = data.columns
channels = [channel.replace("'", "") for channel in channels]

ch_type = snakemake.params['ch_type']
ch_types = [ch_type for i in range(0, len(channels))]

data = data.transpose()

info = mne.create_info(ch_names = channels, sfreq = sfreq, ch_types = ch_types)
raw = mne.io.RawArray(data, info)

pd.to_pickle(raw, snakemake.output['raw'])
