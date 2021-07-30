#!/usr/bin/env python

"""
    Demonstrates the steps required to create genotype networks, perform
    analyses, and write results to files, for a selected subset of genotype
    sets from the input file using Genonets.

    Use the following command to run the script:

        python -u -m genonets.sample.selective --alphabet=DNA --include-indels \
            --input-file=input.csv --tau=0.35 --output-path=results_selective \
            --verbose

    Output files will be generated in 'results_selective/'

    :license: MIT, see LICENSE for more details.
"""

# Command-line parser
from genonets.cmdl_handler import CmdParser

# Interface to Genonets API
from genonets.interface import Genonets

# For analysis type constants
from genonets.constants import AnalysisConstants as Ac


def process(args):
    # Create the Genonets object. This will load the input file into
    # memory.
    gn = Genonets(args)

    # Use 'gn' to create genotype networks for all genotype sets.
    gn.create()

    # Perform only 'Robustness' and 'Evolvability' analyses on just two of
    # the genotype sets available in the input file, i.e., 'Foxa2' and 'Bbx'.
    gn.analyze(["Foxa2", "Bbx"], analyses=[Ac.ROBUSTNESS, Ac.EVOLVABILITY])

    # Write the given genotype networks to files in GML format. For a genotype
    # network with two or more components, two files are generated: One
    # corresponds to the entire network with all components, and the other
    # corresponds to the dominant component only.
    gn.save(["Foxa2", "Bbx"])

    # Save genotype network level measures for the given genotype sets to
    # 'Genotype_set_measures.csv'.
    gn.save_network_results(["Foxa2", "Bbx"])

    # Save all genotype level measures for the given genotype sets to
    # 'Foxa2_genotype_measures.csv' and 'Bbx_genotype_measures.csv' files.
    gn.save_genotype_results(["Foxa2", "Bbx"])


if __name__ == "__main__":
    # Parse the command line arguments using the Genonets command line handler,
    # and pass the list of arguments to 'process()'.
    process(CmdParser().get_args())

    # Print message to indicate processing is done.
    print("\nDone.\n")
