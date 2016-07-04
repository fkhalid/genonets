#!/usr/bin/env python

"""
    genonets_exmpl_simple
    ~~~~~~~~~~~~~~~~~~~~~

    Demonstrates the steps required to create genotype networks, perform analyses, 
    and write results to files using Genonets.

    Use the following command to run the script:
    'python genonets_exmpl_simple.py DNA true data/genonets_sample_input.txt 0.35 results_simple'

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

    # Perform all available analyses on all genotype networks.
    gn.analyze()

    # Write all genotype networks to files in GML format. For a genotype network
    # with two or more components, two files are generated: One corresponds to the
    # entire network with all components, and the other corresponds to the dominant
    # component only.
    gn.save()

    # Save all genotype network level measures to 'Genotype_set_measures.txt'.
    gn.save_network_results()

    # Save all genotype level measures to '<genotypeSetName>_genotype_measures.txt'
    # files. One file per genotype set is generated.
    gn.save_genotype_results()


if __name__ == "__main__":
    # Parse the command line arguments using the Genonets command line handler, and
    # pass the list of arguments to 'process()'.
    process(CmdParser().getArgs())

    # Print message to indicate processing is done.
    print("\nDone.\n")
