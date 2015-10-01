#!/usr/bin/env python

import inspect
import os
import logging
import gridjug
import jug

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('log.txt')
fh.setLevel(logging.DEBUG)
logger.addHandler(fh)


# http://stackoverflow.com/a/50905/2366781
THIS_DIR = os.path.dirname(
    os.path.abspath(inspect.getfile(inspect.currentframe()))
)

print THIS_DIR

PRIMES_JUGFILE = os.path.join(THIS_DIR, 'primes.py')

PRIMES_JUGDIR = os.path.join(THIS_DIR, 'primes.jugdata')

FAILING_JUGFILE = os.path.join(THIS_DIR, 'failing.py')

TEMP_DIR = os.path.join(THIS_DIR, 'tmp')
if not os.path.isdir(TEMP_DIR):
    os.mkdir(TEMP_DIR)

GRIDMAP_PARAMS = {
    'local': False,
    'require_cluster': True,
    'temp_dir': TEMP_DIR,
    'quiet': False,
    'queue': 'all.q',
    'interpreting_shell': '/bin/bash',
    'copy_env': False,
}


gridjug.grid_jug(jugfile=PRIMES_JUGFILE, jugdir=PRIMES_JUGDIR, **GRIDMAP_PARAMS)

_, jugspace = jug.init(jugfile=PRIMES_JUGFILE, jugdir=jugdir)

assert jug.value(jugspace['primes10']) == [True, True, False, True, False, True, False, False, False]
