
configfile: 'config.yml'

rule all:
  input:
    expand('output/{project}/W1/C1/emg/filter.pkl', project = config['project'])

#
# Convert the raw EMG .txt files to MNE objects for further processing
#
rule emg_mne:
  input:
    raw = 'raw/{sid}_{cell}_EMG_rest.txt'
  output:
    raw = 'output/{project}/{sid}/{cell}/emg/raw.pkl'
  params:
    samplingRate = config['samplingRate']
  conda: 'env/mne.yml'
  script: 'python/create_mne.py'

#
# Filter the EMG data
#
rule emg_filter:
  input:
    raw = 'output/{project}/{sid}/{cell}/emg/raw.pkl'
  output:
    filter = 'output/{project}/{sid}/{cell}/emg/filter.pkl'
  params:
    drop_below = config['filter']['drop_below'],
    drop_above = config['filter']['drop_above']
  conda: 'env/mne.yml'
  script: 'python/filter.py'
