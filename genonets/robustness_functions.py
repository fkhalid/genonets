"""
    robustness_functions
    ~~~~~~~~~~~~~~~~~~~~

    Contains functions used for genotype and phenotype robustness computations.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import numpy as np
from tqdm import tqdm


class RobustnessAnalyzer:
    # Constructor
    def __init__(self, network, netBuilder, is_double_stranded, verbose):
        # Reference to the network on which to perform this
        # analysis
        self.network = network

        # Get a reference to the NetworkUtils object
        self.netBuilder = netBuilder

        # Get a reference to the BitSeqManipulator in use
        self.bitManip = netBuilder.bitManip

        # Flag indicating whether the genotypes correspond to double
        # stranded molecules
        self.is_double_stranded = is_double_stranded

        # Reference to the list of tuples, (seq, robustness)
        # populated when required
        self.robustnessVals = None

        self._verbose = verbose

    # Returns the arithmetic mean of the all values in the
    # robustness list.
    def getAvgRobustness(self):
        # Calculate the mean of all robustness values for the network
        return np.mean(self.getRobustnessAll())

    # Returnes a list of tuples of the form (sequence, robustness) that
    # consist of robustness values for all sequences represented by
    # vertices of the given network.
    def getRobustnessAll(self, recompute=False):
        if self._verbose:
            print('\n')

        # If either the robustness values have not been computed already,
        # or the caller has explicitly requested re-computing,
        if not self.robustnessVals or recompute:
            if self._verbose:
                iterable = tqdm(self.network.vs["sequences"])
            else:
                iterable = self.network.vs["sequences"]

            # Get robustness values for all sequences in the network
            self.robustnessVals = [
                self.getGenotypeRobustness(seq)
                for seq in iterable
            ]

        return self.robustnessVals

    # Calculates sequence robustness, which is defined as Genotype
    # robustness in,
    # Ref: "Robustness and evolvability: a paradox resolved". Andreas
    # Wagner, Proc. R. Soc. B 2008 275 91-100;
    # DOI: 10.1098/rspb.2007.1137. Published 7 January 2008
    def getGenotypeRobustness(self, sequence):
        # Degree of the vertex corresponding to the sequence
        degree = self.network.degree(sequence)

        # Get a list of all possible 1-neighbors for this sequence
        allNeighbors = self.bitManip.generateNeighbors(
            self.bitManip.seqToBits(sequence))

        # If reverse complements should be considered,
        if self.is_double_stranded:
            # Reverse complement of the focal genotype
            rc = self.bitManip.get_reverse_complement(
                self.bitManip.seqToBits(sequence))

            # If the reverse complement is also a 1-mutant,
            if rc in allNeighbors:
                # Do not consider it as one of all possible 1-mutants
                allNeighbors.remove(rc)

            # Make sure a genotype and its reverse complement are not both
            # considered separate neighbors
            for neighbor in list(allNeighbors):
                # Reverse complement of the neighbor
                rc_n = self.bitManip.get_reverse_complement(neighbor)

                # If the reverse complement is not the same as the neighbor,
                # and the reverse complement is already in the list of
                # external neighbors, remove the neighbor itself from the
                # list.
                if rc_n != neighbor and rc_n in allNeighbors:
                    allNeighbors.remove(neighbor)

        # Count the number of all possible 1-neighbors of this sequence
        numNeighbors = len(allNeighbors)

        try:
            return float(degree) / float(numNeighbors)
        except ZeroDivisionError:
            return 0
