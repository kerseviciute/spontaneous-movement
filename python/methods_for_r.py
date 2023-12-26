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
