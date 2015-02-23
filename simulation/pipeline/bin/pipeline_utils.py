import yaml


def write_data(file_name, data):
    """
    Simple YAML output format for scoring or other pipeline associated data.

    :param file_name: the file to write output
    :param data: the data object
    """
    with open(file_name, 'w') as h_out:
        yaml.dump(data, h_out, default_flow_style=False)


def read_data(file_name):
    """
    Read simple YAML format file of pipeline associated data.

    :param file_name: the file from which to read input
    :return: loaded data object
    """
    with open(file_name, 'r') as h_in:
        return yaml.load(h_in)
