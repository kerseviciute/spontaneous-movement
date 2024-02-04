import pandas as pd
from snakemake.io import expand

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
    data = 'output/{project}/{animal_id}/{cell_name}/emg/filter.pkl'
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
    data = 'output/{project}/{animal_id}/{cell_name}/vm/filter.pkl'
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
    all_events = 'output/{project}/{animal_id}/{cell_name}/movement_all.csv',
    # Events filtered by previous event ending some time before the start (good for analysing event onset)
    events_start = 'output/{project}/{animal_id}/{cell_name}/movement_start_final.csv',
    # Events filtered by next event starting some time after the end (good for analysing event offset)
    events_end = 'output/{project}/{animal_id}/{cell_name}/movement_end_final.csv'
  params:
    tkeoThreshold = config['movement']['tkeoThreshold'],
    maxTimeApart = config['movement']['maxTimeApart'],
    maxLength = config['movement']['maxLength'],
    minLength = config['movement']['minLength'],
    calmBeforeEvent = config['movement']['calmBeforeEvent'],
    minAmplitude = config['movement']['minAmplitude']
  conda: 'env/r.yml'
  script: 'R/extract_movement.R'

#
# Analyse filtered EMG data and detect no movement episodes.
#
rule extract_events_no_move:
  input:
    data = 'output/{project}/{animal_id}/{cell_name}/emg_data.csv'
  output:
    all_events = 'output/{project}/{animal_id}/{cell_name}/no_movement_all.csv',
    events = 'output/{project}/{animal_id}/{cell_name}/no_movement_final.csv'
  params:
    tkeoThreshold = config['no_movement']['tkeoThreshold'],
    maxTimeApart = config['no_movement']['maxTimeApart'],
    minLength = config['no_movement']['minLength'],
    expandBy = config['no_movement']['expandBy']
  conda: 'env/r.yml'
  script: 'R/extract_no_movement.R'

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
    no_movement_filtered = expand('output/{{project}}/{sid}/no_movement_final.csv', sid = samples['Location']),
    movement = expand('output/{{project}}/{sid}/movement_all.csv', sid = samples['Location']),
    movement_filtered = expand('output/{{project}}/{sid}/movement_start_final.csv', sid = samples['Location'])
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
    movement_filtered = expand('output/{{project}}/{sid}/movement_start_final.csv', sid = samples['Location'])
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
    rest = expand('output/{{project}}/{sid}/no_movement_final.csv', sid = samples['Location']),
    move = expand('output/{{project}}/{sid}/movement_start_final.csv', sid = samples['Location'])
  output:
    data = 'output/{project}/figures/combined_vm_data.csv'
  params:
    prefix = 'output/{project}'
  conda: 'env/mne.yml'
  script: 'python/figures/combine_vm.py'

rule combine_rest_info:
  input:
    samples = config['sample_sheet'],
    data = expand('output/{{project}}/{sid}/no_movement_final.csv', sid = samples['Location'])
  output:
    data = 'output/{project}/figures/combined_rest_information.csv'
  params:
    prefix = 'output/{project}'
  conda: 'env/mne.yml'
  script: 'python/figures/combine_rest.py'

rule ap_count:
  input:
    samples = config['sample_sheet'],
    rest = expand('output/{{project}}/{sid}/no_movement_final.csv', sid = samples['Location']),
    move = expand('output/{{project}}/{sid}/movement_start_final.csv', sid = samples['Location']),
    vm = expand('output/{{project}}/{sid}/vm/filter.pkl', sid = samples['Location'])
  output:
    data = 'output/{project}/figures/ap_count.csv'
  params:
    prefix = 'output/{project}'
  conda: 'env/mne.yml'
  script: 'python/figures/ap_count.py'

rule event_data_start:
  input:
    samples = config['sample_sheet'],
    emg = expand('output/{{project}}/{sid}/emg/filter.pkl', sid = samples['Location']),
    movement = expand('output/{{project}}/{sid}/movement_start_final.csv', sid = samples['Location']),
    vm = expand('output/{{project}}/{sid}/vm/filter.pkl', sid = samples['Location'])
  output:
    data = 'output/{project}/figures/event_data_start.csv'
  params:
    prefix = 'output/{project}'
  conda: 'env/mne.yml'
  script: 'python/figures/event_data_start.py'

rule event_data_end:
  input:
    samples = config['sample_sheet'],
    emg = expand('output/{{project}}/{sid}/emg/filter.pkl', sid = samples['Location']),
    movement = expand('output/{{project}}/{sid}/movement_end_final.csv', sid = samples['Location']),
    vm = expand('output/{{project}}/{sid}/vm/filter.pkl', sid = samples['Location'])
  output:
    data = 'output/{project}/figures/event_data_end.csv'
  params:
    prefix = 'output/{project}'
  conda: 'env/mne.yml'
  script: 'python/figures/event_data_end.py'

rule ap_time_start:
  input:
    samples = config['sample_sheet'],
    movement = expand('output/{{project}}/{sid}/movement_start_final.csv', sid = samples['Location']),
    vm = expand('output/{{project}}/{sid}/vm/filter.pkl', sid = samples['Location'])
  output:
    data = 'output/{project}/figures/ap_time_start.csv'
  params:
    prefix = 'output/{project}'
  conda: 'env/mne.yml'
  script: 'python/figures/ap_time_start.py'

rule ap_time_end:
  input:
    samples = config['sample_sheet'],
    movement = expand('output/{{project}}/{sid}/movement_end_final.csv', sid = samples['Location']),
    vm = expand('output/{{project}}/{sid}/vm/filter.pkl', sid = samples['Location'])
  output:
    data = 'output/{project}/figures/ap_time_end.csv'
  params:
    prefix = 'output/{project}'
  conda: 'env/mne.yml'
  script: 'python/figures/ap_time_end.py'

######################################################################
# Figures
######################################################################

rule figure1:
  input:
    emg_data = 'output/{project}/W1/C2/emg_data.csv',
    vm_data = 'output/{project}/W1/C2/vm_data.csv',
    movement_information = 'output/{project}/W1/C2/movement_start_final.csv',
    no_movement_information = 'output/{project}/W1/C2/no_movement_final.csv'
  output:
    png = '{project}/www/figure1.png'
  conda: 'env/r.yml'
  script: 'R/figures/figure1.R'

# TODO: vector memory exhausted (limit reached?)
rule figure2:
  input:
    sample_sheet = config['sample_sheet'],
    movement_information = expand('output/{{project}}/{sid}/movement_start_final.csv', sid = samples['Location']),
    rest_information = 'output/{project}/figures/combined_rest_information.csv',
    vm = 'output/{project}/figures/combined_vm_data.csv',
    ap_count = 'output/{project}/figures/ap_count.csv'
  output:
    png = '{project}/www/figure2.png'
  conda: 'env/r.yml'
  script: 'R/figures/figure2.R'

rule figure3:
  input:
    event_data = 'output/{project}/figures/event_data_start.csv',
    ap_time = 'output/{project}/figures/ap_time_start.csv'
  output:
    png = '{project}/www/figure3.png'
  conda: 'env/r.yml'
  script: 'R/figures/figure3.R'

rule figure4:
  input:
    event_data = 'output/{project}/figures/event_data_start.csv'
  output:
    png = '{project}/www/figure4.png'
  conda: 'env/r.yml'
  script: 'R/figures/figure4.R'

rule figure5:
  input:
    event_data = 'output/{project}/figures/event_data_end.csv',
    ap_time = 'output/{project}/figures/ap_time_end.csv'
  output:
    png = '{project}/www/figure5.png'
  conda: 'env/r.yml'
  script: 'R/figures/figure5.R'

rule figure6:
  input:
    event_data = 'output/{project}/figures/event_data_end.csv'
  output:
    png = '{project}/www/figure6.png'
  conda: 'env/r.yml'
  script: 'R/figures/figure6.R'
