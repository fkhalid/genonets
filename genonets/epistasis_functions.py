"""
    evolvability_functions
    ~~~~~~~~~~~~~~~~~~~~~~

    Contains functions used for epistasis computation.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import random

from tqdm import tqdm

from genonets_constants import EpistasisConstants as epi


class EpistasisAnalyzer:
    # Constructor
    def __init__(self, network, netUtils, seqToEscrDict, delta, bitManip, sample_size, verbose):
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

        # Sample size to use for epistasis
        self._sample_size = sample_size

        # Flag to indicate whether output should be verbose
        self._verbose = verbose

        # To store the computed squares, where each element in a
        # square is a vertex id
        self.squares = []

        # Dict {verterx ID: [square ids in which this vertex appears]}
        self.vtxToSqrs = {v.index: [] for v in network.vs}

        # To store the computed epistasis class per square
        self.sqrEpi = None

        if self._verbose:
            print('\nConstructing neighbor map ...')

        # Neighbor map
        self._neighbor_map = {
            s: set(self.netUtils.getNeighborSequences(s, self.network))
            for s in tqdm(network.vs['sequences'])
        }

    # Computes epistasis for the given network
    def getEpiAll(self):
        # Dictionary to keep count of the number of squares that
        # show each class of epistasis.
        epistasis = {
            epi.MAGNITUDE: 0,
            epi.SIGN: 0,
            epi.RECIPROCAL_SIGN: 0,
            epi.NO_EPISTASIS: 0
        }

        # Get all squares
        self.getSquares()

        # List of epistasis class corresponding to each square
        self.sqrEpi = []

        print('Calculating epistasis ...')

        # For each square
        for square in tqdm(self.squares):
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
        for epiClass in epistasis:
            try:
                epistasis[epiClass] = float(epistasis[epiClass]) / float(len(self.squares))
            except ZeroDivisionError:
                epistasis[epiClass] = 0

        return epistasis

    def getEpistasis(self, square):
        epiClass = epi.NO_EPISTASIS

        # Get e-scores for all corners of the square
        esrc_AB = self.seqToEscrDict[self.network.vs[square[3]]['name']]
        esrc_ab = self.seqToEscrDict[self.network.vs[square[0]]['name']]
        esrc_Ab = self.seqToEscrDict[self.network.vs[square[2]]['name']]
        esrc_aB = self.seqToEscrDict[self.network.vs[square[1]]['name']]

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
                pass
        else:
            # No epistasis
            pass

        return epiClass

    def getEpiClass(self, esrc_ab, esrc_aB, esrc_Ab, esrc_AB):
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

    # Return all squares found in the genotype network
    def getSquares(self, recompute=False):
        print('Constructing squares ...')

        # If squares computation has been done already, and the
        # caller has not explicitly asked for re-running the
        # algorithms,
        if self.squares and not recompute:
            # Return the pre-computed results
            return self.squares

        # Set of unique squares (permutation agnostic) with sequences
        sqrs_set = set()

        # Get the list of all sequences
        sequences = self.network.vs["sequences"]

        # For each sequence in the network
        for sequence in tqdm(sequences):
            # Get all 1-neighbor sequences
            neighbors = list(self._neighbor_map[sequence])

            # If the number of neighbors is less than two, there's no point
            # in continuing any further
            if len(neighbors) < 2:
                # Move on to the next sequence
                continue

            # Construct all possible pairs of neighbors, where symmetric pairs
            # are considered only once. Also, pairs that neighbor each other
            # are not considered.
            if self._sample_size == 0:
                pairs = [
                    (neighbors[i], neighbors[j])
                    for i in xrange(len(neighbors) - 1)
                    for j in xrange(i + 1, len(neighbors))
                    if neighbors[j] not in self._neighbor_map[neighbors[i]]
                    and neighbors[i] not in self._neighbor_map[neighbors[j]]
                ]
            else:
                pairs = [
                    (neighbors[i], neighbors[j])
                    for i in random.sample(xrange(len(neighbors) - 1), min(self._sample_size, len(neighbors)))
                    for j in random.sample(xrange(i + 1, len(neighbors)), min(self._sample_size, len(neighbors) - (i + 1)))
                    if neighbors[j] not in self._neighbor_map[neighbors[i]]
                       and neighbors[i] not in self._neighbor_map[neighbors[j]]
                ]

            # For each pair of neighbors
            while pairs:
                pair = pairs.pop()
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
                        square = [
                            self.network.vs.find(sequence).index,
                            self.network.vs.find(pair[0]).index,
                            self.network.vs.find(pair[1]).index,
                            self.network.vs.find(node).index
                        ]
                        square_set = frozenset(square)

                        # If the square has not already been found.
                        # Note: This is a set operation, which means due to
                        # the lack of ordering, all permutations will be tested
                        # with just this one condition.
                        if square_set not in sqrs_set:
                            sqrs_set.add(square_set)
                            self.squares.append(square)

                            for v_id in square:
                                self.vtxToSqrs[v_id].append(len(self.squares) - 1)

                del commonNeighs

            del pairs

        if self._verbose:
            print('No. of squares: \n' + str(len(self.squares)))

        return self.squares

    # Returns a list of neighbors common to both elements in the pair. The
    # parent is not considered as a common neighbor.
    def getCommonNeighbors(self, pair, parent):
        # Get neighbors for the first sequence in the pair
        neighbors1 = self._neighbor_map[pair[0]]

        # Get neighbors for the second sequence in the pair
        neighbors2 = self._neighbor_map[pair[1]]

        # Get a list of common neighbors
        commonNeighbors = neighbors1 & neighbors2

        # Remove the parent
        commonNeighbors.discard(parent)

        return commonNeighbors

    # Return a dict - {verterx id : [square ids in which this vertex appears]}
    def getVertexToSquaresDict(self):
        return self.vtxToSqrs, self.squares

    def getSqrEpi(self):
        return self.sqrEpi
