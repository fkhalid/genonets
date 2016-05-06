
"""
    landscape_functions
    ~~~~~~~~~~~~~~~~~~~

    Wrapper for all landscape analysis functions.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import igraph

from peak_functions import PeakAnalyzer
from path_functions import PathAnalyzer
from epistasis_functions import EpistasisAnalyzer
from genonets_utils import Utils


class Landscape :
    # Constructor
    def __init__(self, network, netUtils, seqToEscrDict, delta, bitManip) :
        # Store reference to the network object
        self.network = network

        # Store reference to the netUtils object
        self.netUtils = netUtils

        # Create the peak analyzer
        self.peakAnalyzer = PeakAnalyzer(network, netUtils, delta)

        # Create the path analyzer
        self.pathAnalyzer = PathAnalyzer(network, netUtils, delta)

        # Create the epistasis analyzer
        self.epiAnalyzer = EpistasisAnalyzer(network, netUtils,
            seqToEscrDict, delta, bitManip)

        # Get a reference to the BitSeqManipulator in use
        self.bitManip = self.netUtils.bitManip

    # ----------------------------------------------------------------
    # Peak analysis methods
    # ----------------------------------------------------------------

    def getPeaks(self, recompute) :
        return self.peakAnalyzer.getPeaks(recompute)

    # ----------------------------------------------------------------
    # Path analysis methods
    # ----------------------------------------------------------------

    def getAccessiblePaths(self, pathLength=0):
        return self.pathAnalyzer.getAccessiblePaths(pathLength)

    # ----------------------------------------------------------------
    # Epistasis analysis methods
    # ----------------------------------------------------------------

    def getEpistasis(self) :
        return self.epiAnalyzer.getEpiAll()

    # ----------------------------------------------------------------
    # Calculate mutational distances from summit for all vertices
    # ----------------------------------------------------------------

    def populateDistsToSummit(self) :
        # Convenient handle for bitManip
        bm = self.bitManip

        # Get the summit sequence
        summit = Utils.getSeqWithMaxScore(self.network,
            self.bitManip.seqLength)

        # Get vertex that represents summit
        trgtVrtx = self.netUtils.getVertex(summit, self.network)

        # Reference to the list of sequences in the network
        vertices = [ self.netUtils.getVertex(seq, self.network) \
                     for seq in self.network.vs["sequences"] ]

        self.network.vs["Distance from Summit"] = \
            [ 	len(self.network.get_shortest_paths(srcVrtx, to=trgtVrtx, \
                        weights=None, mode=igraph.OUT, output="epath")[0]) \
                for srcVrtx in vertices ]
