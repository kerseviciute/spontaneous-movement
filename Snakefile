import pandas as pd

configfile: 'config.yml'

samples = pd.read_csv(config['sample_sheet'])

rule all:
  input:
    expand('output/{project}/report/{page}.html',
      project = config['project'],
      page = config['report']['pages']
    )

rule report_index:
  output:
    report = 'output/{project}/report/index.html'
  params:
    script = 'reports/index.Rmd'
  conda: 'env/r.yml'
  script: 'R/render.R'

rule report_emg_index:
  input:
    reports = expand('output/{{project}}/report/{sid}_emg.html', sid = samples['SID'])
  output:
    report = 'output/{project}/report/index_emg.html'
  params:
    script = 'reports/emg/emg_index.Rmd'
  conda: 'env/r.yml'
  script: 'R/render.R'

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
    script = 'reports/emg/emg.Rmd'
  priority: 10
  threads: 4 # best not to run in parallel (figures get mixed up)
  conda: 'env/r.yml'
  script: 'R/render.R'

rule report_vm_index:
  input:
    reports = expand('output/{{project}}/report/{sid}_vm.html', sid = samples['SID'])
  output:
    report = 'output/{project}/report/index_vm.html'
  params:
    script = 'reports/vm/vm_index.Rmd'
  conda: 'env/r.yml'
  script: 'R/render.R'

rule report_vm:
  input:
    data = 'output/{project}/{sid}/{cell}/vm/filter.pkl',
    no_movement_filtered = 'output/{project}/{sid}/{cell}/emg/no_movement_events.pkl',
    movement_filtered = 'output/{project}/{sid}/{cell}/emg/filtered_movement_events.pkl'
  output:
    report = 'output/{project}/report/{sid}_{cell}_vm.html'
  params:
    script = 'reports/vm/vm.Rmd'
  priority: 10
  threads: 4 # best not to run in parallel (figures get mixed up)
  conda: 'env/r.yml'
  script: 'R/render.R'

rule report_vm_in_events:
  input:
    samples = config['sample_sheet'],
    movement = expand('output/{{project}}/{sid}/emg/filtered_movement_events.pkl', sid = samples['Location']),
    vm = expand('output/{{project}}/{sid}/vm/filter.pkl', sid = samples['Location'])
  output:
    report = 'output/{project}/report/{region}_events.html'
  params:
    script = 'reports/vm_in_events.Rmd'
  threads: 4
  conda: 'env/r.yml'
  script: 'R/render.R'

#
# Generate a sample sheet from the excel data.
#
rule sample_sheet:
  input:
    data = 'raw/cells patched_naive.xlsx'
  output:
    sample_sheet = config['sample_sheet'],
    samples = expand('output/{project}/samples.RDS', project = config['project'])
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
    scale = config['filter']['emg']['scale'],
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
    scale = config['filter']['vm']['scale'],
    ch_type = 'bio'
  conda: 'env/mne.yml'
  script: 'python/filter.py'
