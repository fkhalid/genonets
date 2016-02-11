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

from genonets import cmdl_handler  # For parsing command line arguments
from genonets import genonets_constants  # For analysis type constants
from genonets import genonets_interface as gn_if  # Interface to get the Genonets object


def process(args):
    # Get a reference to the analysis constants. These constants are used
    # to specify analysis types.
    ac = genonets_constants.AnalysisConstants

    # Create the Genonets object. This will load the input file into
    # memory.
    gn = gn_if.Genonets(args)

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
    gn.saveNetResults(["Foxa2", "Bbx"])

    # Save all genotype level measures for the given genotype sets to
    # 'Foxa2_genotype_measures.txt' and 'Bbx_genotype_measures.txt' files.
    gn.saveGenotypeResults(["Foxa2", "Bbx"])


if __name__ == "__main__":
    # Parse the command line arguments using the Genonets command line handler, and
    # pass the list of arguments to 'process()'.
    process(cmdl_handler.CmdParser().getArgs())

    # Print message to indicate processing is done.
    print("Done.\n")
