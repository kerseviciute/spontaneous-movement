
configfile: 'config.yml'

rule all:
  input:
    expand('output/{project}/W1/C1/emg/filter.pkl', project = config['project'])

#
# Convert the raw EMG .txt files to MNE objects for further processing.
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
# Filter the EMG data.
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

#
# Analyse filtered EMG data and detect movement episodes.
#
rule emg_detectMovement:
  input:
    filter = 'output/{project}/{sid}/{cell}/emg/filter.pkl'
  output:
    quality_movement = 'output/{project}/{sid}/{cell}/emg/final_movement.pkl',
    all_movement = 'output/{project}/{sid}/{cell}/emg/all_movement.pkl',
    all_no_movement = 'output/{project}/{sid}/{cell}/emg/all_no_movement.pkl'
  params:
    tkeoMaxFreq = config['movement']['tkeoMaxFreq'],
    threshold = config['movement']['tkeoThreshold'],
    minLength = config['movement']['minEventLength'],
    minBreak = config['movement']['minEventBreak'],
    expandBy = config['movement']['expandBy'],
    minMovementAmplitude = config['movement']['minMovementAmplitude'],
    minAmplitudeDifference = config['movement']['minAmplitudeDifference'],
    calmTime = config['movement']['calmTime'],
    threshold_no_move = config['no_movement']['tkeoThreshold'],
    minLength_no_move = config['no_movement']['minEventLength'],
    minBreak_no_move = config['no_movement']['minEventBreak'],
    expandBy_no_move = config['no_movement']['expandBy']
  conda: 'env/mne.yml'
  script: 'python/detect_movement.py'
