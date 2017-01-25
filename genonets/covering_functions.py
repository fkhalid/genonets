"""
    covering_functions
    ~~~~~~~~~~~~~~~~~~

    Contains functions used for computation of 'phenotype space covering'.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import igraph
import collections


class CoveringAnalyzer:
    def __init__(self, network, net_builder, evo_analyzer, sequence_length, total_phenotypes, is_double_stranded):
        # Reference to the network on which to perform this analysis
        self.network = network

        # Reference to the Network Builder
        self.netBuilder = net_builder

        # EvolvabilityAnalyzer object
        self.evoAnalyzer = evo_analyzer

        # Sequence length for each genotype
        self.seqLength = sequence_length

        # Total number of phenotypes other the focal genotype set
        self.NUM_PHENOTYPES = total_phenotypes - 1

        # Reference to Bit Manipulator
        self.bm = self.netBuilder.bitManip

        # Flag to indicate whether or not reverse complements should
        # be considered
        self.is_double_stranded = is_double_stranded

        # Set - All genotypes within the network in bit format
        self.bit_genotypes = frozenset(
            self.bm.seqToBits(seq)
            for seq in self.network.vs["sequences"])

        # FIXME: for debugging purposes only ...
        self.counter = 0

    def covering_unique(self, radius):
        # TODO: Add check for radius > diameter ...

        # Get the list of all genotypes for this network
        genotypes = self.network.vs["sequences"]

        return [
            self.covering_unique_all(genotype, radius)
            for genotype in genotypes
        ]

    def covering_unique_all(self, genotype, radius):
        return [
            self.covering_unique_genotype(genotype, r)
            for r in xrange(1, radius + 1)
        ]

    def covering_unique_genotype(self, focal_genotype, radius):
        # List of ratios
        ratios = []

        # Set of all genotypes except 'focal_genotype'
        genotypes = set(self.network.vs["sequences"]) - {focal_genotype}

        # If there are no genotypes to process,
        if not genotypes:
            return []

        # Set of phenotypes found in the 1-mutant neighborhood of the focal
        # genotype
        src_evo_trgts = set(
            self.evoAnalyzer.getSeqEvo(focal_genotype)
            .target_reps
            .keys()
        )

        # If the focal genotype does not have any targets,
        if not src_evo_trgts:
            return []

        # Vertex ID corresponding to the focal genotype
        src_vrtx = self.netBuilder.getVertex(focal_genotype, self.network)

        # For each genotype in 'genotypes',
        for genotype in genotypes:
            # Vertex corresponding to 'genotype'
            trgt_vrtx = self.netBuilder.getVertex(genotype, self.network)

            # Measure the mutational distance to the focal genotype
            distance = self.network.shortest_paths_dijkstra(
                src_vrtx.index, trgt_vrtx.index, mode=igraph.OUT)[0][0]

            # If the genotype is 'radius' away from 'focal_genotype',
            if distance == radius:
                # Set of phenotypes found in the 1-mutant neighborhood of the
                # genotype
                trgt_evo_trgts = set(
                    self.evoAnalyzer.getSeqEvo(genotype)
                    .target_reps
                    .keys()
                )

                # Ratio of target phenotypes unique to either the focal
                # or the current genotype.
                ratio = 1 - (
                    float(len(src_evo_trgts & trgt_evo_trgts)) /
                    len(src_evo_trgts))

                # Append to the list of results
                ratios.append(ratio)

        return ratios

    def covering_all(self, radius):
        # FIXME: for debugging purposes only ...
        # print

        # Radius value must be between 1, and sequence length
        if radius not in xrange(1, self.seqLength + 1):
            return []
        
        # For each sequence in the network, compute covering. This
        # results in a list of lists.
        return [
            self.covering(bit_genotype, radius)
            for bit_genotype in self.bit_genotypes
        ]

    def covering(self, bit_sequence, radius):
        # FIXME: for debugging purposes only ...
        # print '{0}\r'.format(self.counter),
        self.counter += 1

        # Construct the result as a named tuple
        result = collections.namedtuple("Result", ["covering", "radius", "targets"])

        # If the focal genotype set is the only genotype set available,
        if self.NUM_PHENOTYPES < 1:
            # There's no point on doing the calculations
            return result(
                [],
                0,
                []
            )

        # List to hold the covering value for each value of 'r';
        # initialization here with all '0's.
        covering = [float("NaN") for _ in xrange(radius)]

        # Set of all target repertoires covered within the radius
        target_reps = set()

        # Convert the input sequence value into the required format. It
        # needs to be an iterable for later convenience.
        genotypes = {bit_sequence}

        # Set of sequences that have already been considered in the
        # previous value of 'r'. These are initialized with all
        # genotypes within the focal GN, since these must always be
        # ignored while calculating overlap with other phenotypes.
        parents = set(self.bit_genotypes)

        # If reverse complements should be considered,
        if self.is_double_stranded:
            # Add the reverse complements to the set
            parents |= {
                self.bm.getReverseComplement(g)
                for g in self.bit_genotypes
            }

        # For radius starting with 1,
        for r in xrange(1, radius + 1):
            # If 'genotypes' is empty, there is no use in proceeding
            # any further.
            if not genotypes:
                break

            # Set to hold the genotypes that have already been computed
            # within the current value of 'r'.
            siblings = set()

            # For each genotype computed during the previous iteration,
            for genotype in genotypes:
                # Compute all possible 1-neighbors
                neighbors = set(self.bm.generateNeighbors(genotype))

                # If reverse complements should be considered,
                if self.is_double_stranded:
                    # Remove redundant genotypes from the list of neighbors
                    neighbors = self.remove_redundant_neighbors(genotype, neighbors)

                # Keep only those neighbors that are not part of the
                # 'parents' set.
                neighbors = neighbors - parents

                # Union of the two sets
                siblings |= neighbors

            # Compute overlap of all elements in 'siblings', with all
            # genotypes in other genotype sets (dominants only).
            overlap_result = self.overlap(siblings)

            # Store the percentage of covered phenotypes
            covering[r-1] = overlap_result.covering

            # Union of the all target phenotypes found so far
            target_reps |= set(overlap_result.targets.keys())

            # If 100% coverage has been achieved, no need to proceed
            # any further.
            if sum(covering[0:r]) >= 1:
                break

            # Current siblings are added to the parents, since we do
            # not want to allow backwards mutations.
            # FIXME: This may result in memory consumption issues ...
            parents |= siblings.copy()

            # Current siblings need to be used as 'genotypes' for the
            # next iteration.
            genotypes = siblings.copy()

        return result(
            covering=covering,
            radius=r,
            targets=target_reps
        )

    def overlap(self, genotypes):
        # Use the 'EvolvabilityAnalyzer' to get a dictionary with
        # {key=genotype set, value=[overlapping genotypes]}
        targets = self.evoAnalyzer.getEvoTargetReps(genotypes)

        # Calculate the value of 'covering'
        covering_val = float(len(targets)) / float(self.NUM_PHENOTYPES)

        # Construct the result as a named tuple
        result = collections.namedtuple("Result", ["covering", "targets"])

        return result(
            covering=covering_val,
            targets=targets
        )

    def remove_redundant_neighbors(self, genotype, neighbors):
        # Reverse complement of the focal genotype
        rc = self.bm.getReverseComplement(genotype)

        # If the reverse complement is also a 1-mutant,
        if rc in neighbors:
            # Do not consider it as one of all possible 1-mutants
            neighbors.remove(rc)

        # Make sure a genotype and its reverse complement are not both
        # considered separate neighbors
        for neighbor in list(neighbors):
            # Reverse complement of the neighbor
            rc_n = self.bm.getReverseComplement(neighbor)

            # If the reverse complement is not the same as the neighbor,
            # and the reverse complement is already in the list of
            # neighbors, remove the neighbor itself from the
            # list.
            if rc_n != neighbor and rc_n in neighbors:
                neighbors.remove(neighbor)

        return neighbors
