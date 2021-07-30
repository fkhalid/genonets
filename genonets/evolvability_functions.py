"""
    evolvability_functions
    ~~~~~~~~~~~~~~~~~~~~~~

    Contains functions used for genotype and phenotype evolvability computations.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import sys

from tqdm import tqdm


class EvolvabilityAnalyzer:
    # Constructor
    def __init__(self, network, dataDict, seqToRepDict, repToGiantDict,
                 rcToSeqDict, bitsToSeqDict, netBuilder,
                 isDoubleStranded, verbose):
        # Reference to the network on which to perform this
        # analysis
        self.network = network

        # Reference to dict: key=rep, value=dict{seq: score}
        self.dataDict = dataDict

        # Copy dict: key=sequence, value=[repertoires]
        self.seqToRepDict = seqToRepDict

        # Copy of dict: key=repertoire, value=giant
        self.repToGiantDict = repToGiantDict

        # Dict {bit reverse complement: string seq} - Holds reverse
        # complement (bit format) to original sequence (string)
        self.rcToSeqDict = rcToSeqDict

        # Dict {bitseq : seq}
        self.bitsToSeqDict = bitsToSeqDict

        # Get a reference to the NetworkUtils object
        self.netBuilder = netBuilder

        # Flag indicating whether the genotypes correspond to double
        # stranded molecules
        self.isDoubleStranded = isDoubleStranded

        # Reference to the bit manipulator object
        self.bm = self.netBuilder.bitManip

        # Reference to the list of evolvability values for all sequences
        self.evoTuples = None

        self._verbose = verbose

        self._counter = 0

    @staticmethod
    def updateSeqToRepDict(seqToRepDict_original, repToGiantDict):
        # Make a deep copy of the dict
        seqToRepDict = seqToRepDict_original.copy()

        # For each sequence key in the dict,
        for seq in seqToRepDict:
            # Update the list of repertoires corresponding to this
            # sequence with only those repertoires for which this
            # sequence exists in the giant.
            seqToRepDict[seq] = [
                rep for rep in seqToRepDict[seq]
                if EvolvabilityAnalyzer.is_sequence_in_network(seq, repToGiantDict[rep])
            ]

        return seqToRepDict

    @staticmethod
    def is_sequence_in_network(sequence, network):
        try:
            network.vs.find(sequence)
            vertex_exists = True
        except ValueError:
            vertex_exists = False

        return vertex_exists

    @staticmethod
    def buildRcToSeqDict(seqToRepDict, bm):
        # Initialize the dict
        rcToSeqDict = dict()

        # For each sequence key in the dict,
        for seq in seqToRepDict.keys():
            # Compute the reverse complement
            rcBitSeq = bm.get_reverse_complement(
                bm.seqToBits(seq))

            # With the reverse complement bit sequence as the key,
            # add the original sequence (string format) as the value.
            rcToSeqDict[rcBitSeq] = seq

        return rcToSeqDict

    @staticmethod
    def buildBitsToSeqDict(seqToRepDict, rcToSeqDict, bm, isDoubleStranded):
        # Create the dictionary - {bit sequence: sequence}
        bitsToSeqDict = {
            bm.seqToBits(seq): seq for seq in seqToRepDict.keys()
        }

        # If we need to consider reverse complements,
        if isDoubleStranded:
            # Use the 'reverse complement to sequence dict' to get the
            # the reverse complements corresponding to all the
            # sequences.
            for rc in rcToSeqDict.keys():
                bitsToSeqDict[rc] = bm.bitsToSeq(rc)

        return bitsToSeqDict

    # Returns repertoire evolvability
    def getReportoireEvo(self):
        # Get evo tuples for all sequences in the network
        evoTuples = self.getEvoAll()

        # Get repertoire lists from the tuples
        repLists = [
            evoTuples[1].keys()
            for evoTuples in evoTuples if evoTuples[1]
        ]

        # Combine all repLists into one list
        targets = []
        for repList in repLists:
            targets.extend([target for target in repList])

        # Remove redundant repertoires
        targets = list(set(targets))

        try:
            evolvability = float(len(targets)) / float(len(self.dataDict) - 1)

            return evolvability, targets
        except ZeroDivisionError:
            return 0, targets

    def getEvoAll(self, recompute=False):
        if self._verbose:
            iterable = tqdm(self.network.vs['sequences'])
        else:
            iterable = self.network.vs['sequences']

        # If either the evolvability values have not been computed already,
        # or the caller has explicitly requested re-computing,
        if not self.evoTuples or recompute:
            # Get tuples: (evo, targetRepsDict) for all sequences in the network
            self.evoTuples = [self.getSeqEvo(seq) for seq in iterable]

        self._counter = 0

        return self.evoTuples

    # Returns evolvability score for the given sequence
    def getSeqEvo(self, sequence):
        # if self._verbose:
        #     sys.stdout.write(str(self._counter) + '\r')
        #     sys.stdout.flush()

        self._counter += 1

        # Get a list of sequences that are 1-neighbors, but not part of the
        # genotype network
        externNeighbors = self.netBuilder.getExternalNeighbors(sequence, self.network)

        # Get dictionary: key=repertoire, value=extNeighs. The keys represents the
        # list of repertoires to which the given sequence could join by undergoing
        # a single mutation.
        targetReps = self.getEvoTargetReps(externNeighbors)

        # Compute and return the evolvability score, i.e., No. of repertoire
        # targets found for this sequence / No. of repertoires available in
        # the dataset.
        # FIXME: ZeroDivisionError can only occur if the denominator is zero.
        # This already be checked at the very beginning of the function to
        # avoid the rest of the computation ...
        try:
            evolvability = float(len(targetReps)) / float(len(self.dataDict) - 1)
        except ZeroDivisionError:
            evolvability = 0

        return (evolvability, targetReps)

    # For the given sequence, list of external neighbors, and the sequence
    # to repertoire dictionary, return a list of repertoires for which the given
    # sequence could develop score >= tau by a single mutation.
    def getEvoTargetReps(self, extNeighs):
        targetReps = dict()

        # For each external neighbor,
        for extNeigh in extNeighs:
            try:
                # Get the sequence corresponding to the bitseq
                extNeighSeq = self.bitsToSeqDict[extNeigh]
            except KeyError:
                # If the bitseq corresponding to sequence could not be found,
                # it means that it does not exist in any other repertoire. We
                # can therefore safely skip the current iteration
                continue

            # If the external neighbor is in any repertoire,
            if extNeighSeq in self.seqToRepDict:
                # Append the corresponding repertoire(s) to the
                # list of target repertoires for this sequence.
                self.appendToTargets(extNeighSeq, targetReps)
            elif self.isDoubleStranded:
                # Get the bit representation for the extNeighSeq
                extNeighBits = self.bm.seqToBits(extNeighSeq)

                # If the external neighbor is amongst the reverse complements,
                if extNeighBits in self.rcToSeqDict:
                    # Get the string representation of the reverse complement
                    # of the external neighbor sequence.
                    # Note: The reason we are using the original sequence rather
                    # than the reverse complement is that if we use the
                    # reverse complement, the visualization in the Genonets
                    # Server would not be able to highlight the target node.
                    strSeq = self.rcToSeqDict[extNeighBits]

                    self.appendToTargets(strSeq, targetReps)

        # Sort the list of sequences in each target repository
        for r in targetReps:
            targetReps[r].sort()

        return targetReps

    def appendToTargets(self, seq, targetReps):
        # For each repertoire that contains the external neighbor,
        for repertoire in self.seqToRepDict[seq]:
            # The focal repertoire should not be considered
            if self.network["name"] in [repertoire, repertoire + "_dominant"]:
                continue

            # We have found a repertoire with a genotype network that
            # has a 1-neighbor to which  the sequence can evolve.
            if repertoire not in targetReps:
                targetReps[repertoire] = []

            targetReps[repertoire].append(seq)

    def compute_external_mutation_types(
            self,
            target_repertoires,
            targets_per_genotype):

        """
        For each genotype set, i.e., repertoire, goes through each genotype
        within the genotype set. If the genotype has one or more external
        neighbors, computes and records the mutation type against each external
        neighbor.

        :param target_repertoires: (list of strings) Names of all genotype sets
                in which there is at least one external neighbor.
        :param targets_per_genotype: (list) Each element is a dict, where the
                key is the name of a target genotype set and the value is a
                list of genotype targets within the target set. E.g.,
                    [{'Set2': ['AAG']}, {'Set2': ['AAG', 'ACC']}]

        :return: (dict) Keys are names of target genotypes sets with at least
                one external neighbor. Each value is a dict, where keys are
                mutation types and values are counts of the types.

        """

        # Dict to return: Mutation statistics per repertoire
        mutation_statistics = {
            repertoire: {} for repertoire in target_repertoires
        }

        # For each genotype in the network,
        for i in range(self.network.vcount()):
            # If the corresponding dict is not empty, i.e., at least one
            # external neighbor exists,
            if targets_per_genotype[i]:
                # Get the bit sequence for the source genotype
                bit_source_genotype = self.bm.seqToBits(
                    self.network.vs[i]['sequences']
                )

                self._update_mutation_statistics_for(
                    bit_source_genotype=bit_source_genotype,
                    genotype_targets=targets_per_genotype[i],
                    mutation_statistics=mutation_statistics
                )

        # Sort the dict so that the target repertoires are in lexicographical
        # order
        mutation_statistics = {
            k: v for k, v in sorted(
                mutation_statistics.items(), key=lambda item: item[0]
            )
        }

        # Sort the mutation types in lexicographical order
        for repertoire in mutation_statistics:
            mutation_statistics[repertoire] = {
                k: v for k, v in sorted(
                    mutation_statistics[repertoire].items(),
                    key=lambda item: item[0]
                )
            }

        return mutation_statistics

    def _update_mutation_statistics_for(
            self,
            bit_source_genotype,
            genotype_targets,
            mutation_statistics):

        """
        Updates mutation_statistics by computing the statistics for the given
        genotype.

        :param bit_source_genotype: (int) Bit value corresponding to the source
                genotype.
        :param genotype_targets: (list) Each element is a dict, where the
                key is the name of a target genotype set and the value is a
                list of genotype targets within the target set. E.g.,
                    [{'Set2': ['AAG']}, {'Set2': ['AAG', 'ACC']}]
        :param mutation_statistics: (dict) The dict with statistics to be
                updated. (dict) Keys are names of target genotypes sets with at
                least one external neighbor. Each value is a dict, where keys
                are mutation types and values are counts of the types.

        :return: (dict) Updated mutation_statistics.

        """

        # For each target repertoire to which this genotype can evolve,
        for repertoire in genotype_targets:
            # If the target repertoire has at least one genotype to
            # which the source genotype can evolve,
            if repertoire:
                # For each genotype in the current target repertoire,
                for target_genotype in genotype_targets[repertoire]:
                    # Determine the mutation type
                    mutation_type = self.bm.get_mutation_type(
                        bit_source_genotype,
                        self.bm.seqToBits(target_genotype)
                    )

                    # Update the statistics
                    if mutation_type in mutation_statistics[repertoire]:
                        mutation_statistics[repertoire][mutation_type] += 1
                    else:
                        mutation_statistics[repertoire][mutation_type] = 1
