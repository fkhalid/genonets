#!/usr/bin/env python

from setuptools import setup

# Read the contents of the README file
import os
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

setup(
    name='genonets',
    version='1.1.10',
    description='Framework for creating and analyzing genotype networks from '
                'data.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Fahad Khalid',
    author_email='fahad.khalid@ieu.uzh.ch',
    url='https://github.com/fkhalid/genonets',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Operating System :: OS Independent',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='genonets genotype phenotype network',
    packages=['genonets'],
    package_data={
        'genonets': [
            'sample/*.py',
            'sample/data/*.txt'
        ]
    },
    install_requires=[
        'python-igraph==0.7.1.post6',
        'numpy>=1.8.2'
    ]
)
