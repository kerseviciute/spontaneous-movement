
configfile: 'config.yml'

rule all:
  input:
    expand('output/{project}/W1/C3/emg/filter.pkl', project = config['project']),
    expand('output/{project}/W1/C3/vm/filter.pkl', project = config['project']),
    expand('output/{project}/W1/C3/emg/filtered_movement_events.pkl', project = config['project']),
    expand('output/{project}/W1/C3/emg/no_movement_events.pkl', project = config['project'])

rule sample_sheet:
  input:
    data = 'raw/cells patched_naive.xlsx'
  output:
    sample_sheet = 'sample_sheet.csv'
  conda: 'env/r.yml'
  script: 'R/sample_sheet.R'

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

#
# Filter the movement episodes to keep only good quality data.
#
rule filter_events_move:
  input:
    filter = 'output/{project}/{sid}/{cell}/emg/filter.pkl',
    movement = 'output/{project}/{sid}/{cell}/emg/movement_events.pkl',
    no_movement = 'output/{project}/{sid}/{cell}/emg/no_movement_events.pkl'
  output:
    movement = 'output/{project}/{sid}/{cell}/emg/filtered_movement_events.pkl'
  params:
    amplitude_diff_times = config['movement']['amplitudeDifference'],
    calm_time = config['movement']['calmTime'],
    window_size = config['movement']['windowSize'],
    resample_factor = config['movement']['resampleFactor'],
    max_merge_time = config['movement']['maxMergeTime'],
    calm_amplitude_difference = config['movement']['calmAmplitudeDifference'],
    min_movement_time = config['movement']['minMovementTime']
  conda: 'env/mne.yml'
  script: 'python/filter_events.py'

#
# Analyse filtered EMG data and detect no movement episodes.
#
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

rule report_emg:
  input:
    raw = 'output/{project}/{sid}/{cell}/emg/raw.pkl',
    data = 'output/{project}/{sid}/{cell}/emg/filter.pkl',
    no_movement_filtered = 'output/{project}/{sid}/{cell}/emg/no_movement_events.pkl',
    movement = 'output/{project}/{sid}/{cell}/emg/movement_events.pkl',
    movement_filtered = 'output/{project}/{sid}/{cell}/emg/filtered_movement_events.pkl'
  output:
    report = 'output/{project}/report/{sid}_{cell}_emg.html'
  params:
    script = 'emg.Rmd'
  conda: 'env/r.yml'
  script: 'R/render.R'

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
