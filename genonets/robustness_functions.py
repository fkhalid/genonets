
"""
    robustness_functions
    ~~~~~~~~~~~~~~~~~~~~

    Contains functions used for genotype and phenotype robustness computations.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import numpy as np


class RobustnessAnalyzer :
    # Constructor
    def __init__(self, network, netBuilder) :
        # Reference to the network on which to perform this
        # analysis
        self.network = network

        # Get a reference to the NetworkUtils object
        self.netBuilder = netBuilder

        # Get a reference to the BitSeqManipulator in use
        self.bitManip = netBuilder.bitManip

        # Refernce to the list of tuples, (seq, robustness)
        # poupulated when required
        self.robustnessVals = None

    # Returns the arithmetic mean of the all values in the
    # robustness list.
    def getAvgRobustness(self) :
        # Calculate the mean of all robustness values for the network
        return np.mean(self.getRobustnessAll())

    # Returnes a list of tuples of the form (sequence, robustness) that
    # consist of robustness values for all sequences represented by
    # vertices of the given network.
    def getRobustnessAll(self, recompute=False) :
        # If either the robustness values have not been computed already,
        # or the caller has explicitly requested re-computing,
        if not self.robustnessVals or recompute :
            # Get robustness values for all sequences in the network
            self.robustnessVals = [self.getGenotypeRobustness(seq) \
                                    for seq in self.network.vs["sequences"]]

        return self.robustnessVals

    # Calculates sequence robustness, which is defined as Genotype
    # robustness in,
    # Ref: "Robustness and evolvability: a paradox resolved". Andreas
    # Wagner, Proc. R. Soc. B 2008 275 91-100;
    # DOI: 10.1098/rspb.2007.1137. Published 7 January 2008
    def getGenotypeRobustness(self, sequence) :
        # Get the vertex that corresponds to this sequence
        vertex = self.netBuilder.getVertex(sequence, self.network)

        # Get vertex degree
        degree = self.network.degree(vertex)

        # Get a list of all possible 1-neighbors for this sequence
        allNeighbors = self.bitManip.generateNeighbors(
            self.bitManip.seqToBits(sequence))

        # Count the number of all possible 1-neighbors of this sequence
        numNeighbors = len(allNeighbors)

        try :
            return float(degree) / float(numNeighbors)
        except ZeroDivisionError :
            return 0
