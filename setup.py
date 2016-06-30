#!/usr/bin/env python

from setuptools import setup

setup(
    name='genonets',
    version='1.1.0',
    description='Framework for creating and analyzing genotype networks from data.',
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
        'python-igraph>=0.6',
        'numpy>=1.8.2'
    ]
)
