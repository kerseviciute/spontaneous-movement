def read_pickle(filename):
    import pickle

    with open(filename, 'rb') as file:
        data = pickle.load(file)

    return data


def save_pickle(data, filename):
    import pickle

    with open(filename, 'wb') as file:
        pickle.dump(data, file)
