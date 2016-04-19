"""
    covering_functions
    ~~~~~~~~~~~~~~~~~~

    Contains functions used for computation of 'phenotype space covering'.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""


class CoveringAnalyzer:
    def __init__(self, network, net_builder, evo_analyzer, sequence_length, total_phenotypes):
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

        # Set - All genotypes within the network in bit format
        self.bit_genotypes = frozenset(
            self.bm.seqToBits(seq)
            for seq in self.network.vs["sequences"])

        # FIXME: for debugging purposes only ...
        self.counter = 0

    def covering_all(self, radius):
        # FIXME: for debugging purposes only ...
        print

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
        print '{0}\r'.format(self.counter),
        self.counter += 1

        # If the focal genotype set is the only genotype set available,
        if self.NUM_PHENOTYPES < 1:
            # There's no point on doing the calculations
            return 0

        # List to hold the covering value for each value of 'r';
        # initialization here with all '0's.
        covering = [0 for i in xrange(radius)]

        # Convert the input sequence value into the required format. It
        # needs to be an iterable for later convenience.
        genotypes = {bit_sequence}

        # Set of sequences that have already been considered in the
        # previous value of 'r'. These are initialized with all
        # genotypes within the focal GN, since these must always be
        # ignored while calculating overlap with other phenotypes.
        parents = set(self.bit_genotypes)

        # For radius starting with 1,
        for r in xrange(1, radius + 1):
            # If 'genotypes' is empty, there is not use in proceeding
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

                # Keep only those neighbors that are not part of the
                # 'parents' set.
                neighbors = neighbors - parents

                # Union of the two sets
                siblings |= neighbors

            # Compute overlap of all elements in 'siblings', with all
            # genotypes in other genotype sets (dominants only).
            covering[r-1] = self.overlap(siblings)

            # Current siblings are added to the parents, since we do
            # not want to allow backwards mutations.
            # FIXME: This may result in memory consumption issues ...
            parents |= siblings.copy()

            # Current siblings need to be used as 'genotypes' for the
            # next iteration.
            genotypes = siblings.copy()

        return covering

    def overlap(self, genotypes):
        # Use the 'EvolvabilityAnalyzer' to get a dictionary with
        # {key=genotype set, value=[overlapping genotypes]}
        targets = self.evoAnalyzer.getEvoTargetReps(genotypes)

        return float(len(targets)) / float(self.NUM_PHENOTYPES)


