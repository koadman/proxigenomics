import yaml
import sys
import os.path
import glob
import fnmatch


def read(file_name):
    """
    Read the application configuration file (YAML format), exits
    on IOError.

    :param file_name: configuration file path
    :return: yaml configuration object
    :
    """
    try:
        with open(file_name, 'r') as h_cf:
            return yaml.load(h_cf)
    except IOError as e:
        print e
        sys.exit(e.errno)


def find_files(path, pattern, remove_top=True):
    matches = []
    for root, dirs, files in os.walk(path):
        for f in fnmatch.filter(files, pattern):
            if remove_top:
                matches.append(
                    os.path.join(
                        '/'.join(root.split('/')[1:]), f))
            else:
                matches.append(os.path.join(root, f))
    return matches


def get_files(path, suffix):
    return [os.path.abspath(f) for f in glob.glob(os.path.join(path, '*.{0}'.format(suffix))) if os.path.isfile(f)]


def get_folders(path):
    return [os.path.abspath(dn) for dn in glob.glob(os.path.join(path, '*')) if os.path.isdir(dn)]


def get_communities(config):
    """
    Return the list of community folders
    :param config: application config object
    :return: list of community folders
    """
    return get_folders(config['community']['folder'])


def get_wgs_reads(path, config):
    """
    Return the pair of generated WGS files for a given community
    :param path: containing path of WGS reads
    :param config: application config object
    :return: list of read files (R1, R2)
    """
    return [os.path.join(path, '{0}{1}.fq'.format(config['wgs_base'], n)) for n in range(1, 3)]