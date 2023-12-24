
configfile: 'config.yml'

rule all:
  input:
    expand('output/{project}/W1/C3/emg/filter.pkl', project = config['project']),
    expand('output/{project}/W1/C3/vm/filter.pkl', project = config['project']),
    expand('output/{project}/W1/C3/emg/movement_events.pkl', project = config['project']),
    expand('output/{project}/W1/C3/emg/no_movement_events.pkl', project = config['project'])

#
# Convert the raw EMG .txt files to MNE objects for further processing.
#
rule emg_mne:
  input:
    raw = 'raw/{sid}_{cell}_EMG_rest.txt'
  output:
    raw = 'output/{project}/{sid}/{cell}/emg/raw.pkl'
  params:
    samplingRate = config['samplingRate'],
    ch_type = 'emg'
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
    drop_below = config['filter']['emg']['drop_below'],
    drop_above = config['filter']['emg']['drop_above'],
    ch_type = 'emg'
  conda: 'env/mne.yml'
  script: 'python/filter.py'

#
# Analyse filtered EMG data and detect movement episodes.
#
# TODO: split movement detection and movement filtering into separate files
# rule emg_detectMovement:
#   input:
#     filter = 'output/{project}/{sid}/{cell}/emg/filter.pkl'
#   output:
#     quality_movement = 'output/{project}/{sid}/{cell}/emg/quality_movement.pkl',
#     all_movement = 'output/{project}/{sid}/{cell}/emg/all_movement.pkl',
#     all_no_movement = 'output/{project}/{sid}/{cell}/emg/all_no_movement.pkl',
#     final_movement = 'output/{project}/{sid}/{cell}/emg/final_movement.pkl'
#   params:
#     tkeoMaxFreq = config['movement']['tkeoMaxFreq'],
#     threshold = config['movement']['tkeoThreshold'],
#     minLength = config['movement']['minEventLength'],
#     minBreak = config['movement']['minEventBreak'],
#     expandBy = config['movement']['expandBy'],
#     minMovementAmplitude = config['movement']['minMovementAmplitude'],
#     minAmplitudeDifference = config['movement']['minAmplitudeDifference'],
#     calmTime = config['movement']['calmTime'],
#     threshold_no_move = config['no_movement']['tkeoThreshold'],
#     minLength_no_move = config['no_movement']['minEventLength'],
#     minBreak_no_move = config['no_movement']['minEventBreak'],
#     expandBy_no_move = config['no_movement']['expandBy']
#   conda: 'env/mne.yml'
#   script: 'python/detect_movement.py'

rule extract_events_move:
  input:
    filter = 'output/{project}/{sid}/{cell}/emg/filter.pkl'
  output:
    events = 'output/{project}/{sid}/{cell}/emg/movement_events.pkl'
  params:
    tkeoMaxFreq = config['movement']['tkeoMaxFreq'],
    threshold = config['movement']['tkeoThreshold'],
    minLength = config['movement']['minEventLength'],
    minBreak = config['movement']['minEventBreak'],
    expandBy = config['movement']['expandBy'],
    detectMovement = True
  conda: 'env/mne.yml'
  script: 'python/extract_events.py'

rule filter_events_move:
  input:
    filter = 'output/{project}/{sid}/{cell}/emg/filter.pkl',
    events = 'output/{project}/{sid}/{cell}/emg/movement_events.pkl'
  output:
    events = 'output/{project}/{sid}/{cell}/emg/filtered_movement_events.pkl'
  params:
    calmTime = config['movement']['calmTime']
  conda: 'env/mne.yml'
  script: 'python/extract_events.py'

rule extract_events_no_move:
  input:
    filter = 'output/{project}/{sid}/{cell}/emg/filter.pkl'
  output:
    events = 'output/{project}/{sid}/{cell}/emg/no_movement_events.pkl'
  params:
    tkeoMaxFreq = config['movement']['tkeoMaxFreq'],
    threshold = config['no_movement']['tkeoThreshold'],
    minLength = config['no_movement']['minEventLength'],
    minBreak = config['no_movement']['minEventBreak'],
    expandBy = config['no_movement']['expandBy'],
    detectMovement = False
  conda: 'env/mne.yml'
  script: 'python/extract_events.py'

#
# Convert the raw VM .txt files to MNE objects for further processing.
#
rule vm_mne:
  input:
    raw = 'raw/{sid}_{cell}_Vm_rest.txt'
  output:
    raw = 'output/{project}/{sid}/{cell}/vm/raw.pkl'
  params:
    samplingRate = config['samplingRate'],
    ch_type = 'bio'
  conda: 'env/mne.yml'
  script: 'python/create_mne.py'

#
# Filter the VM data.
#
rule vm_filter:
  input:
    raw = 'output/{project}/{sid}/{cell}/vm/raw.pkl'
  output:
    filter = 'output/{project}/{sid}/{cell}/vm/filter.pkl'
  params:
    drop_below = config['filter']['vm']['drop_below'],
    drop_above = config['filter']['vm']['drop_above'],
    ch_type = 'bio'
  conda: 'env/mne.yml'
  script: 'python/filter.py'
