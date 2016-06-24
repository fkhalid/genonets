#!/usr/bin/env python

"""
    genonets_exmpl_custom
    ~~~~~~~~~~~~~~~~~~~~~

    Illustrates the process of customizing the output by adding information to the
    result file that would not be added by Genonets by default.

    The 'Peaks' analysis is used as an example. By default, the 'Peaks' analysis 
    stores results in a dictionary of the form:
    {key=peakID : value=[genotypes in the peak]}. The sample code in this file
    customizes this dictionary by adding the score value corresponding to each
    genotype in the list. The resulting dictionary is of the format:
    {key=peakID : value=[(genotype1, score1), ..., (genotypeN, scoreN)]},
    i.e., it is a list of tuples.

    Use the following command to run the script:
    'python genonets_exmpl_custom.py DNA true data/genonets_sample_input.txt 0.35 results_custom'

    Output files will be generated in 'results_custom/'

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

    # Perform 'Peaks' analysis
    gn.analyze(analyses=[ac.PEAKS])

    # At this point, the analysis is done. We now need to extract the
    # data we need, i.e., the peaks dictionary.

    # For each genotype set,
    for genotypeSet in gn.genotype_sets():
        # Get the igraph object for the giant
        giant = gn.dominant_network(genotypeSet)

        # Get the dict of peaks {key=peakId : value=[list of sequences in the peak]}
        peaks = giant["Peaks"]

        # Update the dict of peaks, so that instead of just a list of
        # genotypes, we have a list of tuples (sequence, escore).
        newPeaks = {}

        # For each peak,
        for peak in peaks:
            # Initialize the list of tuples
            seqScrTuples = []

            # For each sequence in this peak,
            for sequence in peaks[peak]:
                # Find the corresponding vertex in the giant
                try:
                    vertex = giant.vs.find(sequences=sequence)
                except ValueError:
                    print("Oops! can't find " + sequence + " in giant.")

                # Get the escore
                score = giant.vs[vertex.index]["escores"]

                # Add the tuple to the list of tuples
                seqScrTuples.append((sequence, score))

            # Add the peak and the corresponding list of tuples to the
            # new peaks dict
            newPeaks[peak] = seqScrTuples

        # Replace the peaks dict in giant with the new dict. This is
        # useful because now if you use the genonets functions below to
        # save the network and results, the updated peaks dict with tuples
        # will be written automatically to file.
        giant["Peaks"] = newPeaks

    # Save networks to file in GML format
    gn.save()

    # Save the results to file from network level analysis
    gn.save_network_results()


if __name__ == "__main__":
    # Parse the command line arguments using the Genonets command line handler, and
    # pass the list of arguments to 'process()'.
    process(CmdParser().getArgs())

    # Print message to indicate processing is done.
    print("\nDone.\n")
