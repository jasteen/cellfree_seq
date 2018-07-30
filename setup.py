#!/usr/bin/env python

from setuptools import setup

setup(
    name='cellfree_seq',
    version='0.1',
    author='Jason Steen',
    author_email='jason.steen@monash.edu',
    packages=['src'],
    entry_points={
        'console_scripts': ['cellfree_seq = src.main:main']
    },
    url='https://github.com/SoutheyLab/hiplexpipe',
    license='LICENSE.txt',
    description='this program will run a cellfree pipeline on any fastqs',
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus == 2.6.3",
        "drmaa == 0.7.6",
        "PyYAML == 3.11"
    ],
)
