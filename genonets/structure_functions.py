"""
    structure_functions
    ~~~~~~~~~~~~~~~~~~~

    Encapsulates functions for structural analyses.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import igraph
from tqdm import tqdm


class StructureAnalyzer:
    # Constructor
    def __init__(self, network, netBuilder, verbose):
        # Reference to the network on which to perform this
        # analysis
        self.network = network

        # Get a reference to the NetworkBuilder object
        self.netBuilder = netBuilder

        self._verbose = verbose

        # Reference to giant
        self.giant = self.netBuilder.getGiantComponent(self.network)

        # Diameter value
        self._diameter = None

    # ----------------------------------------------------------------
    #	Network level properties
    # ----------------------------------------------------------------

    def getComponentSizes(self):
        return self.netBuilder.getComponents(self.network)

    def getNumComponents(self):
        return len(self.netBuilder.getComponents(self.network))

    def getDominantSize(self):
        return self.giant.vcount()

    def getPercentDominantSize(self):
        return float(self.giant.vcount()) / float(self.network.vcount())

    def getEdgeDensity(self):
        return self.giant.density()

    def getAvgClstrCoeff(self):
        return self.giant.transitivity_avglocal_undirected()

    def getAssortativity(self):
        return self.giant.assortativity_degree()

    # Note: Calculation of diameter can be very expensive for large
    #       networks. Therefore, it should be calculated only once.
    # FIXME: Can be very expensive for large networks ...
    def getDiameter(self):
        if self._diameter is None:
            self._diameter = self.giant.diameter()

        return self._diameter

    def getDiameterPath(self):
        if self._verbose:
            print('\nDiameter path:')

        # Get network diameter
        # FIXME: Can be very expensive for large networks ...
        diameter = self.getDiameter()

        # Get the vertex sequence
        vSeq = self.giant.vs

        # Calculate shortest paths between each pair of vertices in the
        # network until one is found of size 'diameter'.
        for i in tqdm(range(len(vSeq) - 1)):
            for j in range(i+1, len(vSeq)):
                # Calculate shortest path between vertx 'i' and vertex 'j'
                sPath = self.giant.get_shortest_paths(vSeq[i], vSeq[j],
                    mode=igraph.OUT, output="vpath")

                # Return the path is the No. of vertices in the path
                # is diameter + 1, i.e, the No. of edges = diameter.
                if len(sPath[0]) == (diameter + 1):
                    return sPath[0]

        return []

    # ----------------------------------------------------------------
    #	Vertex level properties
    # ----------------------------------------------------------------

    def getCoreness(self):
        return self.giant.coreness()

    # Returns the clustering coefficients for all sequences in the giant
    # component.
    def getClusteringCoefficients(self):
        # Get all sequences
        sequences = self.giant.vs

        # Get the list of clustering coefficients for all vertices
        # in giant. The list of sequences is passed as an argument
        # just to make sure that the order of coefficients in the
        # list corresponds to the order of passed sequences.
        clstrCoeffs = self.giant.transitivity_local_undirected(sequences)

        return clstrCoeffs
