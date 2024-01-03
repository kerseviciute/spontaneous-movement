import pandas as pd

configfile: 'config.yml'

samples = pd.read_csv(config['sample_sheet'])

rule all:
  input:
    expand('{project}/{page}.html',
      project = config['project'],
      page = config['report']['pages']
    ),
    expand('{project}/www/{supplementary}',
      project = config['project'],
      supplementary = config['report']['supplementary']
    )

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

rule extract_emg_data:
  input:
    data = 'output/{project}/{animal_id}/{cell_name}/emg/filter.pkl',
  output:
    data = 'output/{project}/{animal_id}/{cell_name}/emg_data.csv'
  params:
    extractTKEO = True,
    tkeoThreshold = config['movement']['tkeoMaxFreq'],
    dataType = 'EMG'
  conda: 'env/mne.yml'
  script: 'python/extract_data.py'

rule extract_vm_data:
  input:
    data = 'output/{project}/{animal_id}/{cell_name}/vm/filter.pkl',
  output:
    data = 'output/{project}/{animal_id}/{cell_name}/vm_data.csv'
  params:
    extractTKEO = False,
    dataType = 'VM'
  conda: 'env/mne.yml'
  script: 'python/extract_data.py'

rule extract_events_move:
  input:
    data = 'output/{project}/{animal_id}/{cell_name}/emg_data.csv'
  output:
    all_events = 'output/{project}/{animal_id}/{cell_name}/events1.csv',
    events = 'output/{project}/{animal_id}/{cell_name}/events2.csv'
  params:
    tkeoThreshold = config['movement']['tkeoThreshold'],
    maxTimeApart = config['movement']['maxTimeApart'],
    maxLength = config['movement']['maxLength'],
    minLength = config['movement']['minLength'],
    calmBeforeEvent = config['movement']['calmBeforeEvent'],
    minAmplitude = config['movement']['minAmplitude']
  conda: 'env/r.yml'
  script: 'R/extract_events.R'

#
# Analyse filtered EMG data and detect no movement episodes.
#
# TODO: use the same R script as for movement detection
rule extract_events_no_move:
  input:
    filter = 'output/{project}/{sid}/{cell}/emg/filter.pkl'
  output:
    events = 'output/{project}/{sid}/{cell}/no_movement_events.csv'
  params:
    tkeoMaxFreq = config['movement']['tkeoMaxFreq'],
    threshold = config['no_movement']['tkeoThreshold'],
    minLength = config['no_movement']['minLength'],
    minBreak = config['no_movement']['maxTimeApart'],
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


# ######################################################################
# # Reports
# ######################################################################

rule report_summary:
  output:
    report = '{project}/index.html'
  params:
    script = 'reports/index.Rmd'
  conda: 'env/r.yml'
  script: 'R/render.R'

# May take a while to run!
rule report_emg:
  input:
    samples = config['sample_sheet'],
    raw = expand('output/{{project}}/{sid}/emg/raw.pkl', sid = samples['Location']),
    data = expand('output/{{project}}/{sid}/emg/filter.pkl', sid = samples['Location']),
    no_movement_filtered = expand('output/{{project}}/{sid}/no_movement_events.csv', sid = samples['Location']),
    movement = expand('output/{{project}}/{sid}/events1.csv', sid = samples['Location']),
    movement_filtered = expand('output/{{project}}/{sid}/events2.csv', sid = samples['Location'])
  output:
    report = '{project}/emg.html'
  params:
    script = 'reports/emg.Rmd',
    prefix = 'output/{project}'
  conda: 'env/r.yml'
  script: 'R/render.R'

# May take a while to run!
rule report_vm:
  input:
    samples = config['sample_sheet'],
    data = expand('output/{{project}}/{sid}/vm/filter.pkl', sid = samples['Location']),
    movement_filtered = expand('output/{{project}}/{sid}/events2.csv', sid = samples['Location'])
  output:
    report = '{project}/vm.html'
  params:
    script = 'reports/vm.Rmd',
    prefix = 'output/{project}'
  conda: 'env/r.yml'
  script: 'R/render.R'

# ######################################################################
# # Figure data
# ######################################################################

rule combine_vm:
  input:
    samples = config['sample_sheet'],
    vm = expand('output/{{project}}/{sid}/vm/filter.pkl', sid = samples['Location']),
    rest = expand('output/{{project}}/{sid}/no_movement_events.csv', sid = samples['Location']),
    move = expand('output/{{project}}/{sid}/events2.csv', sid = samples['Location'])
  output:
    data = 'output/{project}/figures/combined_vm_data.csv'
  params:
    prefix = 'output/{project}'
  conda: 'env/mne.yml'
  script: 'python/figures/combine_vm.py'

rule combine_rest_info:
  input:
    samples = config['sample_sheet'],
    data = expand('output/{{project}}/{sid}/no_movement_events.csv', sid = samples['Location'])
  output:
    data = 'output/{project}/figures/combined_rest_information.csv'
  params:
    prefix = 'output/{project}'
  conda: 'env/mne.yml'
  script: 'python/figures/combine_rest.py'

rule ap_count:
  input:
    samples = config['sample_sheet'],
    rest = expand('output/{{project}}/{sid}/no_movement_events.csv', sid = samples['Location']),
    move = expand('output/{{project}}/{sid}/events2.csv', sid = samples['Location']),
    vm = expand('output/{{project}}/{sid}/vm/filter.pkl', sid = samples['Location'])
  output:
    data = 'output/{project}/figures/ap_count.csv'
  params:
    prefix = 'output/{project}'
  conda: 'env/mne.yml'
  script: 'python/figures/ap_count.py'

rule event_data:
  input:
    samples = config['sample_sheet'],
    emg = expand('output/{{project}}/{sid}/emg/filter.pkl', sid = samples['Location']),
    movement = expand('output/{{project}}/{sid}/events2.csv', sid = samples['Location']),
    vm = expand('output/{{project}}/{sid}/vm/filter.pkl', sid = samples['Location'])
  output:
    data = 'output/{project}/figures/event_data.csv'
  params:
    prefix = 'output/{project}'
  conda: 'env/mne.yml'
  script: 'python/figures/event_data.py'

# Takes painfully long :)
rule ap_time:
  input:
    all_events = 'output/{project}/figures/event_data.csv'
  output:
    data = 'output/{project}/figures/ap_time.csv'
  conda: 'env/mne.yml'
  script: 'python/figures/ap_time.py'

######################################################################
# Figures
######################################################################

rule figure1:
  input:
    emg_data = 'output/{project}/W1/C2/emg_data.csv',
    vm_data = 'output/{project}/W1/C2/vm_data.csv',
    movement_information = 'output/{project}/W1/C2/events2.csv',
    no_movement_information = 'output/{project}/W1/C2/no_movement_events.csv'
  output:
    png = '{project}/www/figure1.png'
  conda: 'env/r.yml'
  script: 'R/figures/figure1.R'

rule figure2:
  input:
    sample_sheet = config['sample_sheet'],
    movement_information = expand('output/{{project}}/{sid}/events2.csv', sid = samples['Location']),
    rest_information = 'output/{project}/figures/combined_rest_information.csv',
    vm = 'output/{project}/figures/combined_vm_data.csv',
    ap_count = 'output/{project}/figures/ap_count.csv'
  output:
    png = '{project}/www/figure2.png'
  conda: 'env/r.yml'
  script: 'R/figures/figure2.R'

rule figure3:
  input:
    event_data = 'output/{project}/figures/event_data.csv',
    ap_time = 'output/{project}/figures/ap_time.csv'
  output:
    png = '{project}/www/figure3.png'
  conda: 'env/r.yml'
  script: 'R/figures/figure3.R'

rule figure4:
  input:
    event_data = 'output/{project}/figures/event_data.csv'
  output:
    png = '{project}/www/figure4.png'
  conda: 'env/r.yml'
  script: 'R/figures/figure4.R'
