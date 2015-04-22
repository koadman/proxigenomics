import yaml
import sys
import os.path
import glob
import fnmatch
import types


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


def prepend_paths(path, fnames):
    if isinstance(fnames, types.StringTypes):
        fnames = [fnames]
    return [os.path.join(path, fn) for fn in fnames]


def find_files(path, pattern, remove_top=True):
    """
    Search the given path, top down for the file pattern.
    :param path: top of directory tree to search
    :param pattern: the file name pattern to match
    :param remove_top: optionally remove the top path element
    :return: list of file paths
    """
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


def search_up(path, filename):
    """
    Search up a directory path until a file by the given name is
    found. The method stops as the first occurence, that is the
    lowest point in the path containing a file by the given name.
    :param path: the path to search upwards
    :param filename: the file to find
    :return: the containing path with the filename as a string
    :raises: RuntimeError if a folder contains more than one file of the same name. This should never happen.
    """
    while path != '':
        filepath = [fn for fn in os.listdir(path) if fn == filename]
        nfiles = len(filepath)
        if nfiles > 1:
            raise RuntimeError('duplicate file names in same directory {0}'.format(p))
        elif len(filepath) == 1:
            return os.path.join(path, filepath[0])
        path = os.path.dirname(path)
    return None


def get_precedents(root, filename, strip_file=True, prepend_root=True):
    """
    Search a given directory tree top to bottom and return file resources matching
    the given name. Optionally, return the root path and file name itself.
    target actions.
    :param root: the root path to begin the search
    :param filename: the file name to find
    :param strip_file: remove the filename from found path elements
    :param prepend_root: include the root element in the return paths
    :return: the list of paths
    """

    # build a function to use in preparing elements to match
    # requested options.
    if strip_file:
        f1 = lambda x: os.path.dirname(x)
    else:
        f1 = lambda x: x

    if prepend_root:
        f2 = lambda x: os.path.join(root, f1(x))
    else:
        f2 = lambda x: f1(x)

    paths = [f2(pn) for pn in find_files(root, filename)]
    return paths


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