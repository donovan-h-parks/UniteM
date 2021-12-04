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
import ntpath
import argparse
import errno



class CustomHelpFormatter(argparse.HelpFormatter):
    """Provide a customized format for help output.

    http://stackoverflow.com/questions/9642692/argparse-help-without-duplicate-allcaps
    """

    def _get_help_string(self, action):
        """Place default value in help string."""
        h = action.help
        if '%(default)' not in action.help:
            if action.default != '' and action.default != [] and action.default is not None and type(action.default) != bool:
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [
                        argparse.OPTIONAL, argparse.ZERO_OR_MORE]

                    if action.option_strings or action.nargs in defaulting_nargs:
                        if '\n' in h:
                            lines = h.splitlines()
                            lines[0] += ' (default: %(default)s)'
                            h = '\n'.join(lines)
                        else:
                            h += ' (default: %(default)s)'
            return h

    def _format_action_invocation(self, action):
        """Removes duplicate ALLCAPS with positional arguments."""
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar

        else:
            parts = []

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append(option_string)

                return '%s %s' % (', '.join(parts), args_string)

            return ', '.join(parts)

    def _get_default_metavar_for_optional(self, action):
        return action.dest.upper()

    def _get_default_metavar_for_positional(self, action):
        return action.dest


def logger_setup(log_dir, log_file, program_name, version, silent):
    """Setup loggers.

    Two logger are setup which both print to the stdout and a 
    log file when the log_dir is not None. The first logger is
    named 'timestamp' and provides a timestamp with each call,
    while the other is named 'no_timestamp' and does not prepend
    any information. The attribution 'is_silent' is also added
    to each logger to indicate if the silent flag is thrown.

    Parameters
    ----------
    log_dir : str
        Output directory for log file.
    log_file : str
        Desired name of log file.
    program_name : str
        Name of program.
    version : str
        Program version number.
    silent : boolean
        Flag indicating if output to stdout should be suppressed.
    """

    # setup general properties of loggers
    timestamp_logger = logging.getLogger('timestamp')
    timestamp_logger.setLevel(logging.DEBUG)
    log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                   datefmt="%Y-%m-%d %H:%M:%S")

    no_timestamp_logger = logging.getLogger('no_timestamp')
    no_timestamp_logger.setLevel(logging.DEBUG)

    # setup logging to console
    timestamp_stream_logger = logging.StreamHandler(sys.stdout)
    timestamp_stream_logger.setFormatter(log_format)
    timestamp_logger.addHandler(timestamp_stream_logger)

    no_timestamp_stream_logger = logging.StreamHandler(sys.stdout)
    no_timestamp_stream_logger.setFormatter(None)
    no_timestamp_logger.addHandler(no_timestamp_stream_logger)

    timestamp_logger.is_silent = False
    no_timestamp_stream_logger.is_silent = False
    if silent:
        timestamp_logger.is_silent = True
        timestamp_stream_logger.setLevel(logging.ERROR)
        no_timestamp_stream_logger.is_silent = True

    if log_dir:
        make_sure_path_exists(log_dir)
        timestamp_file_logger = logging.FileHandler(
            os.path.join(log_dir, log_file), 'a')
        timestamp_file_logger.setFormatter(log_format)
        timestamp_logger.addHandler(timestamp_file_logger)

        no_timestamp_file_logger = logging.FileHandler(
            os.path.join(log_dir, log_file), 'a')
        no_timestamp_file_logger.setFormatter(None)
        no_timestamp_logger.addHandler(no_timestamp_file_logger)

    timestamp_logger.info(f'{program_name} v{version}')
    timestamp_logger.info(ntpath.basename(
        sys.argv[0]) + ' ' + ' '.join(sys.argv[1:]))


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
