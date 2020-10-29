
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


class Landscape:
    # Constructor
    def __init__(self, network, netUtils, seqToEscrDict, delta, bitManip):
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

    def getPeaks(self, recompute):
        return self.peakAnalyzer.get_all_peaks(recompute)

    # ----------------------------------------------------------------
    # Path analysis methods
    # ----------------------------------------------------------------

    def getAccessiblePaths(self, pathLength=0):
        return self.pathAnalyzer.getAccessiblePaths(pathLength)

    # ----------------------------------------------------------------
    # Epistasis analysis methods
    # ----------------------------------------------------------------

    def getEpistasis(self):
        return self.epiAnalyzer.getEpiAll()

    # ----------------------------------------------------------------
    # Calculate mutational distances from summit for all vertices
    # ----------------------------------------------------------------

    def populateDistsToSummit(self):
        # Get the summit sequence
        summit = Utils.getSeqWithMaxScore(
            self.network, self.bitManip.seqLength)

        # Get vertex that represents summit
        trgtVrtx = self.netUtils.getVertex(summit, self.network)

        self.network.vs["Distance from Summit"] = [
            len(
                self.network.get_shortest_paths(
                    srcVrtx, to=trgtVrtx, weights=None,
                    mode=igraph.OUT, output="epath")[0]
            )
            for srcVrtx in self.network.vs
        ]
