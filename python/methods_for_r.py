def read_pickle(filename):
    import pickle

    with open(filename, 'rb') as file:
        data = pickle.load(file)

    return data


def plot_all_events(data, channels, movement = None, no_movement = None, limit = 4, channel_col = 'Channel'):
    from matplotlib import pyplot as plt
    import numpy as np

    figure, axes = plt.subplots(nrows = len(data.ch_names), figsize = (20, 0.5 * len(data.ch_names)))
    plt.subplots_adjust(
        left = 0.1,
        bottom = 0.1,
        right = 0.9,
        top = 0.9,
        wspace = 0,
        hspace = 0
    )

    x = data.times

    for i, channel in enumerate(channels):
        y = data.get_data(picks = [channel])[0]

        axes[i].plot(x, y, linewidth = 0.5, color = 'black')

        if limit is not None:
            axes[i].set_ylim(-limit, limit)

        axes[i].get_xaxis().set_ticks([])
        axes[i].get_yaxis().set_ticks([])
        axes[i].set_ylabel(channel, rotation = 0, labelpad = 60, loc = 'center')
        axes[i].set_xlim(0, np.max(data.times))

        if movement is not None:
            events = movement[movement[channel_col] == channel]
            for index, row in events.iterrows():
                start = int(row['EventStart'])
                end = int(row['EventEnd'])

                axes[i].plot(x[start:end], y[start:end], linewidth = 1, color = 'red')

        if no_movement is not None:
            events = no_movement[no_movement[channel_col] == channel]
            for index, row in events.iterrows():
                start = int(row['EventStart'])
                end = int(row['EventEnd'])

                axes[i].plot(x[start:end], y[start:end], linewidth = 1, color = 'blue')

        if i != 0:
            axes[i].spines['top'].set_visible(False)

        if i != len(channels) - 1:
            axes[i].spines['bottom'].set_visible(False)

    plt.show()


def read_sample_event_data(sample_name, movement = True, time_shift = 0, event_length = 0.25):
    """
    Reads sample event data.

    :param sample_name: Sample identifier.
    :param movement: Whether to read movement or no movement data.
    :param time_shift: How much time to record before the event (negative shift) or after the event
                       (positive shift) (in seconds). If the shift is not possible, will return
                       zero-padded values.
    :param event_length: How much of the event time to return (in seconds). None corresponds to the
                         complete event.
    :return: np.array (number of events x time points)
    """
    import re
    import numpy as np

    animal_id = re.match('(W[0-9]+).*', sample_name).group(1)
    cell_name = re.match('.*(C[0-9]+)', sample_name).group(1)

    prefix = f'output/spontaneous-movement/{animal_id}/{cell_name}/'

    sample_data = read_pickle(f'{prefix}/vm/filter.pkl')

    if movement:
        sample_events = read_pickle(f'{prefix}/emg/filtered_movement_events.pkl')
    else:
        sample_events = read_pickle(f'{prefix}/emg/no_movement_events.pkl')

    before_event = time_shift < 0
    time_shift = np.abs(time_shift)

    n_expected_points = sample_data.info['sfreq'] * (time_shift + event_length) + 1

    event_data = []
    for i, event in sample_events.iterrows():
        if before_event:
            start = event['Start'] - time_shift
            start = np.max([start, np.min(sample_data.times)])
            end = event['Start'] + event_length
        else:
            start = event['End'] - event_length
            end = event['End'] + time_shift
            end = np.min([end, np.max(sample_data.times)])

        channel = event['ChannelId']

        channel_data = sample_data.copy().crop(tmin = start, tmax = end).get_data(picks = [channel]).flatten()

        # If the time shift was not possible, will have to append a few points here and there
        if len(channel_data) < n_expected_points:
            print('Too short')
            average = np.mean(channel_data)
            append = np.ones(int(n_expected_points - len(channel_data)))
            append = append * average

            assert len(append) + len(channel_data) == n_expected_points

            if before_event:
                channel_data = np.concatenate([append, channel_data])
            else:
                channel_data = np.concatenate([channel_data, append])

        event_data.append(channel_data)

    event_data = np.array(event_data)

    return event_data


def read_event_data(sample_ids, movement = True, time_shift = 0, event_length = 0.25):
    import numpy as np

    event_data = []
    for sample_id in sample_ids:
        event_data.append(
            read_sample_event_data(
                sample_id,
                movement = movement,
                time_shift = time_shift,
                event_length = event_length
            ))

    event_data = np.concatenate(event_data)
    return event_data

