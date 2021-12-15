#!/usr/bin/env python3

from setuptools import setup, find_packages

import os


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    version_file = open(os.path.join(setup_dir, 'unitem', 'VERSION'))
    return version_file.readline().strip()


if __name__ == '__main__':

    dirName = os.path.dirname(__file__)
    if dirName and os.getcwd() != dirName:
        os.chdir(dirName)

    setup(
        name='unitem',
        version=version(),
        author='Donovan Parks',
        author_email='donovan.parks@gmail.com',
        maintainer='Donovan Parks',
        maintainer_email='donovan.parks@gmail.com',
        packages=find_packages(),
        scripts=['bin/unitem'],
        package_data={'unitem': ['VERSION', '../markers/*']},
        url='http://pypi.python.org/pypi/unitem/',
        license='GPLv3',
        description='Ensemble binning strategies for combining the output of multiple binning methods.',
        classifiers=[
            'Development Status :: 4 - Beta',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=["numpy>=1.0", "svgwrite>=1.1.9"],
    )
