"""
    path_functions
    ~~~~~~~~~~~~~~

    Contains functions used for computation of accessible mutational paths.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import igraph

# from memory_profiler import profile

from genonets_utils import Utils


class PathAnalyzer:
    # Constructor
    def __init__(self, network, netUtils, delta):
        # Reference to the network on which to perform this
        # analysis
        self.network = network

        # Get a reference to the NetworkUtils object
        self.netUtils = netUtils

        # Get a reference to the BitSeqManipulator in use
        self.bitManip = netUtils.bitManip

        # Keep a copy of the delta value
        self.delta = delta

        # Summit vertex, to be populated later.
        self.summitId = None

        # Longest mutational path in the network (among all shortest paths)
        self.max_path_length = 0

        # Dict to store all accessible paths. {vertexId : [paths]}
        self.allPathsToPeak = self.initPathsToPeak()

    def initPathsToPeak(self):
        return {vId: [] for vId in range(len(self.network.vs["sequences"]))}

    def getSummitId(self):
        return self.summitId

    def getAllPathsToPeak(self):
        return self.allPathsToPeak

    def getPathsThruVtxs(self):
        # List: Each element is the No. of paths through the vertex Id
        # 		corresponding to the index
        pathsThruVtx = [0 for i in range(self.network.vcount())]

        # Go through accessible paths corresponding to each vertex
        for vtxId in range(self.network.vcount()):
            # Get the list of accessible paths for this vertex
            vtxPaths = self.network.vs[vtxId]["pathsToSummit"]

            # For each path,
            for path in vtxPaths:
                # For each vertex in the path,
                for vtx in path:
                    # Increment the No. of paths that go through
                    # this vertex
                    pathsThruVtx[vtx] += 1

        return pathsThruVtx

    # Get the ratio of accessible paths to all paths of the
    # given length.
    def getAccessiblePaths(self, pathLength=0):
        # Stats for the entire network
        totalPaths = 0  # Total mutational paths (only shortest)
        allAccPaths = 0  # No. of accessible paths (only shortest)

        # Get the sequence with the highest score. This sequence
        # will represent the global peak. All paths use this
        # sequence as the target.
        summit = Utils.getSeqWithMaxScore(self.network, self.bitManip.seqLength)

        # Get the vertex object that represents the summit
        trgtVrtx = self.netUtils.getVertex(summit, self.network)

        # Store a copy of the summit
        self.summitId = trgtVrtx.index

        # Get a list of all sequences in the network
        sequences = self.network.vs["sequences"]
        # Remove the target itself
        sequences.remove(summit)

        # For each sequence in the network
        for source in sequences:
            # If we only need to calculate all accessible paths, regardless
            # of the path length,
            if pathLength == 0:
                # Get all shortest paths as well as all accessible
                # shortest paths from source to summit
                self.getShortestAccPaths(source, trgtVrtx, pathLength)
            else:   # Ratios should be calculated
                shrtPaths, accPaths = self.getShortestAccPaths(source,
                                                               trgtVrtx,
                                                               pathLength + 1)

                # If at least one shortest path of length == 'pathLength' was
                # found,
                if shrtPaths:
                    # Increment the counts
                    totalPaths += float(len(shrtPaths))
                    allAccPaths += float(len(accPaths))

        try:
            return float(allAccPaths) / float(totalPaths)
        except ZeroDivisionError:
            return 0

    # From within all shortest paths between source sequence and the
    # given target vertex, computes all accessible paths.
    # For computation, returns only those paths that are of the given
    # length. However, all accessible paths are stored regardless of
    # size. This is then used in the visualization.
    # @profile
    def getShortestAccPaths(self, source, trgtVrtx, pathLength):
        # Get the source and target vertices
        srcVrtx = self.netUtils.getVertex(source, self.network)

        # Get all shortest paths between source and target
        allShrtPaths = self.network.get_all_shortest_paths(srcVrtx,
                                                           trgtVrtx,
                                                           mode=igraph.OUT)

        # If we only need to calculate all accessible paths, regardless
        # of the path length,
        if pathLength == 0:
            # Get all shortest accessible paths
            shrtAccPaths = [
                self.network.vs[path].indices
                for path in allShrtPaths
                if self.isAccessible(path)
            ]

            # Store the paths
            self.allPathsToPeak[srcVrtx.index].extend(shrtAccPaths)

            # Update the value of the longest path
            if len(allShrtPaths[0]) > self.max_path_length:
                self.max_path_length = len(allShrtPaths[0])

            return None, None
        else:   # If only paths of 'path length' are required,
            # If the shortest path length is the same as the required
            # length,
            if len(allShrtPaths[0]) == pathLength:
                allShrtAccPaths = [
                    self.network.vs[path]["sequences"]
                    for path in allShrtPaths
                    if self.isAccessible(path)
                ]
            else:
                # We did not find a path of the required length
                allShrtPaths = None
                allShrtAccPaths = None

            return allShrtPaths, allShrtAccPaths

    # Determines whether the given path is an accessible path, i.e.,
    # scores on this path increase monotonously.
    def isAccessible(self, path):
        isAcc = True

        # Get a list of escores for the sequences in the path
        escores = self.network.vs[path]["escores"]

        # Place holder to keep track of the highest score
        # encountered inside the loop
        maxYet = -0.5
        # For each escore in the list
        for i in range(len(escores) - 1):
            # If the current score is higher than the max so far,
            if escores[i] > maxYet:
                # Assign current escore to max
                maxYet = escores[i]

            # If the next sequence in the path is in a lower bin,
            if maxYet - self.delta > escores[i + 1]:
                # The increase in e-score is not monotonous. Hence,
                # the path is not accessible.
                isAcc = False
                break

        return isAcc
