#!/usr/bin/env python

"""
    genonets_exmpl_selective
    ~~~~~~~~~~~~~~~~~~~~~~~~

    Demonstrates the steps required to create genotype networks, perform analyses, and
    write results to files, for a selected subset of genotype sets from the input file
    using Genonets.

    Use the following command to run the script:
    'python genonets_exmpl_selective.py DNA true data/genonets_sample_input.txt 0.35 results_selective'

    Output files will be generated in 'results_selective/'

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

from genonets.cmdl_handler import CmdParser  # For parsing command line arguments
from genonets.genonets_interface import Genonets  # Interface to Genonets API
from genonets.genonets_constants import AnalysisConstants as ac  # For analysis type constants


def process(args):
    # Create the Genonets object. This will load the input file into
    # memory.
    gn = Genonets(args)

    # Use 'gn' to create genotype networks for all genotype sets.
    gn.create()

    # Perform only 'Robustness' and 'Evolvability' analyses on just two of
    # the genotype sets available in the input file, i.e., 'Foxa2' and 'Bbx'.
    gn.analyze(["Foxa2", "Bbx"], analyses=[ac.ROBUSTNESS, ac.EVOLVABILITY])

    # Write the given genotype networks to files in GML format.
    # For a genotype network with two or more components, two files are generated:
    # One corresponds to the entire network with all components, and the other
    # corresponds to the dominant component only.
    gn.save(["Foxa2", "Bbx"])

    # Save genotype network level measures for the given genotype sets to
    # 'Genotype_set_measures.txt'.
    gn.save_network_results(["Foxa2", "Bbx"])

    # Save all genotype level measures for the given genotype sets to
    # 'Foxa2_genotype_measures.txt' and 'Bbx_genotype_measures.txt' files.
    gn.save_genotype_results(["Foxa2", "Bbx"])


if __name__ == "__main__":
    # Parse the command line arguments using the Genonets command line handler, and
    # pass the list of arguments to 'process()'.
    process(CmdParser().getArgs())

    # Print message to indicate processing is done.
    print("\nDone.\n")
