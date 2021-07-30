#!/usr/bin/env python

"""
    Demonstrates the minimal code required to create genotype networks,
    perform analyses, and write results to files using Genonets with default
    settings.

    Use the following command to run the script:

        python -u -m genonets.sample.minimal --alphabet=DNA --include-indels \
            --input-file=input.csv --tau=0.35 --output-path=results_minimal \
            --verbose

    Output files will be generated in 'results_minimal/'

    :license: MIT, see LICENSE for more details.

"""

# Command-line parser
from genonets.cmdl_handler import CmdParser

# Interface to Genonets API
from genonets.interface import Genonets


if __name__ == "__main__":
    # This single line of code,
    #   1. Parses the command line arguments.
    #   2. Loads the input file into memory.
    #   3. Creates genotype networks corresponding to all genotype sets
    #      available in the input data set.
    #   4. Performs all available analyses on all genotype networks.
    #   5. Writes all genotype networks to files in GML format.
    #   6. Writes all genotype network level results to file.
    #   7. Writes all genotype level results to files.
    #   8. Creates the phenotype network.
    #   9. Writes the phenotype network to file in GML format.
    Genonets(CmdParser().get_args(), process=True)

    # Print message to indicate processing is done.
    print('\nDone.\n')
