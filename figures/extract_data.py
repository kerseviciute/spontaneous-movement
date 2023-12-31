import pandas as pd
import numpy as np
import mne

pd.to_pickle(snakemake, '.extract_data.py.pkl')
# snakemake = pd.read_pickle('.extract_data.py.pkl')

data = pd.read_pickle(snakemake.input['data'])
vm = pd.read_pickle(snakemake.input['vm'])
movement = pd.read_pickle(snakemake.input['movement'])
no_movement = pd.read_pickle(snakemake.input['no_movement'])

info = []

for index, channel in enumerate(data.ch_names):
    ch_data = data.get_data(picks = [channel])[0]
    vm_data = vm.get_data(picks = [index])[0]

    df = pd.DataFrame({
        'EMG': ch_data,
        'Time': data.times,
        'Channel': channel
    })

    ch_movement = movement[ movement['Channel'] == channel ]
    ch_no_movement = no_movement[ no_movement['Channel'] == channel ]

    movement_true = []
    for index, event in ch_movement.iterrows():
        idx = np.logical_and(df['Time'] >= event['Start'], df['Time'] <= event['End'])
        movement_true.append(idx)

    movement_true = np.array(movement_true)
    movement_true = np.logical_xor.reduce(movement_true, 0).astype(np.int32)

    no_movement_true = []
    for index, event in ch_no_movement.iterrows():
        idx = np.logical_and(df['Time'] >= event['Start'], df['Time'] <= event['End'])
        no_movement_true.append(idx)

    no_movement_true = np.array(no_movement_true)
    no_movement_true = np.logical_xor.reduce(no_movement_true, 0).astype(np.int32)

    info.append(pd.DataFrame({
        'VM': vm_data,
        'EMG': ch_data,
        'Time': data.times,
        'Channel': channel,
        'Movement': movement_true == 1,
        'NoMovement': no_movement_true == 1
    }))

info = pd.concat(info)

info.to_csv(snakemake.output['data'])
