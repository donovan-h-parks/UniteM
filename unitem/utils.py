###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import sys
import re
import logging
import errno


def check_file_exists(input_file):
    """Check if file exists."""

    if not os.path.exists(input_file) or not os.path.isfile(input_file):
        logger = logging.getLogger('timestamp')
        logger.error('Input file does not exists: ' + input_file + '\n')
        sys.exit(1)


def check_dir_exists(input_dir):
    """Check if directory exists."""

    if not os.path.exists(input_dir) or not os.path.isdir(input_dir):
        logger = logging.getLogger('timestamp')
        logger.error('Input directory does not exists: ' + input_dir + '\n')
        sys.exit(1)


def make_sure_path_exists(path):
    """Create directory if it does not exist."""

    if not path:
        # lack of a path qualifier is acceptable as this
        # simply specifies the current directory
        return

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logger = logging.getLogger('timestamp')
            logger.error('Specified path could not be created: ' + path + '\n')
            sys.exit(1)


def is_executable(fpath):
    """Check if file is executable.
    This is a Python implementation of the linux
    command 'which'.

    Parameters
    ----------
    fpath : str
        Path to file.

    Returns
    -------
    boolean
        True if executable, else False.
    """

    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    """Return path to program.

    This is a Python implementation of the linux
    command 'which'.
    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    Parameters
    ----------
    program : str
        Name of executable for program.

    Returns
    -------
    str
        Path to executable, or None if it isn't on the path.
    """

    fpath, _fname = os.path.split(program)
    if fpath:
        if is_executable(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_executable(exe_file):
                return exe_file

    return None


def check_on_path(program, exit_on_fail=True):
    """Check if program is on the system path.

    Parameters
    ----------
    program : str
        Name of executable for program.
    exit_on_fail : boolean
        Exit program with error code -1 if program in not on path.

    Returns
    -------
    boolean
        True if program is on path, else False.
    """

    if which(program):
        return True

    if exit_on_fail:
        print('%s is not on the system path.' % program)
        sys.exit(1)

    return False


def check_dependencies(programs, exit_on_fail=True):
    """Check if all required programs are on the system path.

    Parameters
    ----------
    programs : iterable
        Names of executable programs.
    exit_on_fail : boolean
        Exit program with error code -1 if any program in not on path.

    Returns
    -------
    boolean
        True if all programs are on path, else False.
    """

    for program in programs:
        if not check_on_path(program, exit_on_fail):
            return False

    return True


def alphanumeric_sort(l):
    """Sorts the given iterable alphanumerically.

    http://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python

    Parameters
    ----------
    l : iterable
        The iterable to be sorted.

    Returns
    -------
    iterable
        Iterable sorted alphanumerically.
    """
    def convert(text): return int(text) if text.isdigit() else text

    def alphanum_key(key): return [convert(c)
                                   for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)
