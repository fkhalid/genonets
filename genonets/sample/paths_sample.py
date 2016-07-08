#!/usr/bin/env python

"""
    Use the following command to run the script:
    'python -u paths_sample.py DNA true data/genonets_sample_input.txt 0.35 results_paths -rc -v'

    Output files will be generated in 'results_simple/'

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

from genonets.cmdl_handler import CmdParser  # For parsing command line arguments
from genonets.genonets_interface import Genonets  # Interface to Genonets API


def process(args):
    # Create the Genonets object. This will load the input file into
    # memory.
    gn = Genonets(args)

    # Use 'gn' to create genotype networks for all genotype sets.
    gn.create()

    # Get a list of all accessible mutational paths between 'CCGTCTGC' and 'GCAGCTGC'
    paths_result = gn.accessible_paths_between('Ascl2', 'CCGTCTGC', 'GCAGCTGC')

    # Print paths as lists of lists of vertex IDs
    print('\nPaths as vertex IDs')
    print('Accessible mutational paths: ' + str(paths_result.accessible))
    print('All shortest paths: ' + str(paths_result.all))

    # Get the giant
    giant = gn.dominant_network('Ascl2')

    # Create lists of genotypes from the lists of vertex IDs
    paths_accessible = [
        giant.vs[path]["sequences"]
        for path in paths_result.accessible
    ]
    paths_all = [
        giant.vs[path]["sequences"]
        for path in paths_result.all
    ]

    # Print paths as lists of lists of genotypes
    print('\nPaths as genotypes')
    print('Accessible mutational paths: ' + str(paths_accessible))
    print('All shortest paths: ' + str(paths_all))

    # Scores
    print('\nScores: ')

    vertex_id = gn.vertex_id('Ascl2', 'CCGTCTGC')
    old_score = giant.vs[vertex_id]['escores']
    new_score = giant.vs[vertex_id]['escores'] = 100

    print('Old score for CCGTCTGC: ' + str(old_score))
    print('New score for CCGTCTGC: ' + str(new_score))


if __name__ == "__main__":
    # Parse the command line arguments using the Genonets command line handler, and
    # pass the list of arguments to 'process()'.
    process(CmdParser().getArgs())

    # Print message to indicate processing is done.
    print("\nDone.\n")
