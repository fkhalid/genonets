"""
    evolvability_functions
    ~~~~~~~~~~~~~~~~~~~~~~

    Contains functions used for epistasis computation.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

from genonets_constants import EpistasisConstants as epi


class EpistasisAnalyzer:
    # Constructor
    def __init__(self, network, netUtils, seqToEscrDict, delta, bitManip):
        # Reference to the network on which to perform this
        # analysis
        self.network = network

        # Get reference to the NetworkUtils object
        self.netUtils = netUtils

        # Get reference to dict: key=sequence, value=escore
        self.seqToEscrDict = seqToEscrDict

        # Keep a copy of the delta value
        self.delta = delta

        # Reference to the bit manipulator object
        self.bitManip = bitManip

        # Dict {sequence : vertex ID}
        self.seqToVidDict = self.buildSeqToVidDict()

        # To store the computed squares
        self.squares = None

        # To store the computed epistasis class per square
        self.sqrEpi = None

    # Computes epistasis for the given network
    def getEpiAll(self):
        # Dictionary to keep count of the number of squares that
        # show each class of epistasis.
        epistasis = {epi.MAGNITUDE: 0, epi.SIGN: 0,
                     epi.RECIPROCAL_SIGN: 0, epi.NO_EPISTASIS: 0}

        # Get all squares
        squares = self.getSquares()

        # List of epistasis class corresponding to each square
        self.sqrEpi = []

        # For each square
        for square in squares:
            # Determine to which class of epistasis this square belongs
            epiClass = self.getEpistasis(square)

            # Store epistasis type for this square
            self.sqrEpi.append(epiClass)

            # Increment the count for the epistasis class to which the
            # square belongs
            epistasis[epiClass] += 1

        # There's no need to keep count of no epistasis
        del epistasis[epi.NO_EPISTASIS]

        # Convert counts into ratios
        for epiClass in epistasis.keys():
            try:
                epistasis[epiClass] = float(epistasis[epiClass]) / float(len(squares))
            except ZeroDivisionError:
                epistasis[epiClass] = 0

        return epistasis

    def getEpistasis(self, square):
        epsilon = 0
        epiClass = epi.NO_EPISTASIS

        # Get e-scores for all corners of the square
        esrc_AB = self.seqToEscrDict[square[3]]
        esrc_ab = self.seqToEscrDict[square[0]]
        esrc_Ab = self.seqToEscrDict[square[2]]
        esrc_aB = self.seqToEscrDict[square[1]]

        # Calculate the magnitude of epistasis
        epsilon = esrc_AB + esrc_ab - esrc_Ab - esrc_aB

        # If epsilon is greater than or equal to the noise threshold,
        # there is epistasis
        if abs(epsilon) >= self.delta:
            # Classify the type of epistasis
            epiClass = self.getEpiClass(esrc_ab, esrc_aB, esrc_Ab, esrc_AB)

            # Check if we have a case where all mutational effects are 0
            if epiClass == epi.NO_EPISTASIS:
                # No epistasis
                epsilon = 0
        else:
            # No epistasis
            epsilon = 0

        return epiClass

    def getEpiClass(self, esrc_ab, esrc_aB, esrc_Ab, esrc_AB):
        epiClass = 0

        # Calculate mutational effects
        dE_ab_Ab = esrc_Ab - esrc_ab
        dE_aB_AB = esrc_AB - esrc_aB
        dE_ab_aB = esrc_aB - esrc_ab
        dE_Ab_AB = esrc_AB - esrc_Ab

        # A mutational effect should only be considered if it is larger
        # than delta
        dE_ab_Ab = dE_ab_Ab if abs(dE_ab_Ab) >= self.delta else 0
        dE_aB_AB = dE_aB_AB if abs(dE_aB_AB) >= self.delta else 0
        dE_ab_aB = dE_ab_aB if abs(dE_ab_aB) >= self.delta else 0
        dE_Ab_AB = dE_Ab_AB if abs(dE_Ab_AB) >= self.delta else 0

        # If all mutational effects are 0, there is no epistasis
        if dE_ab_Ab == 0 and dE_aB_AB == 0 and dE_ab_aB == 0 and dE_Ab_AB == 0:
            return epi.NO_EPISTASIS

        # Conditions used to determine epistasis type
        condition1 = abs(dE_ab_Ab + dE_aB_AB) == abs(dE_ab_Ab) + abs(dE_aB_AB)
        condition2 = abs(dE_ab_aB + dE_Ab_AB) < abs(dE_ab_aB) + abs(dE_Ab_AB)

        # Evaluate the the conditions to determine epistasis type
        if condition1 and not condition2:
            # if condition 1 holds but not condition 2, there is
            # magnitude epistasis
            epiClass = epi.MAGNITUDE
        elif condition2 and not condition1:
            # if condition 2 holds but not condition 1, there is
            # reciprocal sign epistasis
            epiClass = epi.RECIPROCAL_SIGN
        else:
            # Otherwise, we have sign epistasis
            epiClass = epi.SIGN

        return epiClass

    # ----------------------------
    # Squares
    # ----------------------------

    def buildSeqToVidDict(self):
        # Reference to the list of sequences
        sequences = self.network.vs["sequences"]

        return {sequences[i]: i for i in range(len(sequences))}

    # Return all squares found in the genotype network
    # TODO: Look into performance optimization ...
    def getSquares(self, recompute=False):
        # If squares computation has been done already, and the
        # caller has not explicitly asked for re-running the
        # algorithms,
        if self.squares and not recompute:
            # Return the pre-computed results
            return self.squares

        # List of squares to be populated
        squares = []

        # Set of unique squares (permutation agnostic) with sequences
        # in bit format
        bitSqrs = set()

        # Get the list of all sequences
        sequences = self.network.vs["sequences"]

        # For each sequence in the network
        for sequence in sequences:
            # Get all 1-neighbor sequences
            neighbors = [self.network.vs[vid]["sequences"]
                         for vid in
                         self.network.neighbors(
                             self.netUtils.getVertex(sequence, self.network))]

            # If the number of neighbors is less than two, there's no point
            # in continuing any further
            if len(neighbors) < 2:
                # Move on to the next sequence
                continue

            # Construct all possible pairs of neighbors, where symmetric pairs
            # are considered only once. Also, pairs that neighbor each other
            # are not considered.
            pairs = [(neighbors[i], neighbors[j])
                     for i in range(len(neighbors) - 1)
                     for j in range(i + 1, len(neighbors))
                     if not self.netUtils.areConnected(
                    neighbors[i], neighbors[j])]

            # For each pair of neighbors
            for pair in pairs:
                # Get neighbors common to both sequences in the pair, but
                # without the sequence itself
                commonNeighs = self.getCommonNeighbors(pair, sequence)

                # For each common neighbor
                for node in commonNeighs:
                    # For a node that is an immediate neighbor of the parent,
                    # it does not make sense to compute epistasis, since a
                    # single mutation from parent is sufficient.
                    if node not in neighbors:
                        # Construct the square
                        square = [sequence, pair[0], pair[1], node]
                        # Square as a set of bitseqs
                        bitSqr = frozenset([self.bitManip.seqToBits(sequence),
                                            self.bitManip.seqToBits(pair[0]),
                                            self.bitManip.seqToBits(pair[1]),
                                            self.bitManip.seqToBits(node)])

                        # If the square has not already been found.
                        # Note: This is a set operation, which means due to
                        # the lack or ordering, all permutations will be tested
                        # with just this one condition.
                        if not bitSqr in bitSqrs:
                            squares.append(square)
                            bitSqrs.add(bitSqr)

        # Keep a reference for re-use
        self.squares = squares

        return squares

    # Returns a list of neighbors common to both elements in the pair. The
    # parent is not considered as a common neighbor.
    def getCommonNeighbors(self, pair, parent):
        # Get neighbors for the first sequence in the pair
        neighbors1 = [self.network.vs[vid]["sequences"] for vid in \
                      self.network.neighbors(self.seqToVidDict[pair[0]])]

        # Get neighbors for the second sequence in the pair
        neighbors2 = [self.network.vs[vid]["sequences"] for vid in \
                      self.network.neighbors(self.seqToVidDict[pair[1]])]

        # Get a list of common neighbors
        commonNeighbors = list(set(neighbors1) & set(neighbors2))

        # Remove the parent
        commonNeighbors.remove(parent)

        return commonNeighbors

    # Return a dict - {verterx id : [square ids in which this vertex appears]}
    def getVertexToSquaresDict(self):
        # Make sure there is at least one square
        if not self.squares or len(self.squares) < 1:
            return {}, []

        # Get the list of all sequences
        sequences = self.network.vs["sequences"]

        # Initialize the dict: Keys = vertex ids corresponding to sequence indices
        # in the list of sequences; Values = empty lists to be populated later with
        # square ids
        vtxToSqrs = {vId: [] for vId in range(len(sequences))}

        # List of squares with vertex ids instead of sequences as square elements
        squares = []

        # For each square,
        for sqrIndx in range(len(self.squares)):
            sqr = []
            # For each sequence in the square,
            for sequence in self.squares[sqrIndx]:
                # Get the vertex id corresponding to the sequence
                seqIndx = sequences.index(sequence)

                # Add the square Id as value to the vertex Id
                vtxToSqrs[seqIndx].append(sqrIndx)

                # Add vertex id to sqr
                sqr.append(seqIndx)

            # Append the square to the list of squares
            squares.append(sqr)

        return vtxToSqrs, squares

    def getSqrEpi(self):
        return self.sqrEpi
