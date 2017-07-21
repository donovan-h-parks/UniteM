#!/usr/bin/env python

from setuptools import setup

import os


def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'unitem', 'VERSION'))
    return versionFile.read().strip()

if __name__ == '__main__':

    dirName = os.path.dirname(__file__)
    if dirName and os.getcwd() != dirName:
        os.chdir(dirName)

    setup(
        name='unitem',
        version=version(),
        author='Donovan Parks',
        author_email='donovan.parks@gmail.com',
        packages=['unitem'],
        scripts=['bin/unitem'],
        package_data={'unitem' : ['VERSION', './checkm_ms/*.ms']},
        url='http://pypi.python.org/pypi/unitem/',
        license='GPL3',
        description='Combines bins produced by independent binning methods.',
        install_requires=[
            "numpy>=1.9.0",
            "biolib>=0.0.19"],
    )
