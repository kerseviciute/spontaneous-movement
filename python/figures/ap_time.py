import pandas as pd
import numpy as np

full = pd.read_csv(snakemake.input['all_events'])


def create_event_id(row):
    return f'{row["EventId"]}_{row["SID"]}'


def find_ap(signal, time, ap_threshold = -10):
    over_threshold = np.where(signal > ap_threshold)[0]

    if len(over_threshold) == 0:
        return []

    idx_diff = np.where(np.diff(over_threshold) != 1)[0]
    ap_count = np.append(over_threshold[0], over_threshold[idx_diff + 1])

    return time[ap_count]


full['Event'] = full.apply(create_event_id, axis = 1)
unique_events = np.unique(full['Event'])

ap_threshold = -20
aps = []
for index, event in enumerate(unique_events):
    print(index)
    signal = full[full['Event'] == event]['VM'].values
    time = full[full['Event'] == event]['Time'].values
    signal_ap = find_ap(signal, time, ap_threshold = ap_threshold)

    aps.append(pd.DataFrame({
        'Trial': event,
        'Time': signal_ap
    }))

aps = pd.concat(aps)

aps.to_csv(snakemake.output['data'])
