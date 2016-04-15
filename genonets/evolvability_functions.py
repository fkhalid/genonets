"""
    evolvability_functions
    ~~~~~~~~~~~~~~~~~~~~~~

    Contains functions used for genotype and phenotype evolvability computations.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""


class EvolvabilityAnalyzer:
    # Constructor
    def __init__(self, network, dataDict, seqToRepDict, repToGiantDict,
                 rcToSeqDict, bitsToSeqDict, netBuilder,
                 isDoubleStranded):
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

    @staticmethod
    def updateSeqToRepDict(seqToRepDict_original, repToGiantDict):
        # Make a deep copy of the dict
        seqToRepDict = seqToRepDict_original.copy()

        # For each sequence key in the dict,
        for seq in seqToRepDict.keys():
            # Update the list of repertoires corresponding to this
            # sequence with only those repertoires for which this
            # sequence exists in the giant.
            seqToRepDict[seq] = [
                rep for rep in seqToRepDict[seq]
                if seq in repToGiantDict[rep].vs["sequences"]
            ]

        return seqToRepDict

    @staticmethod
    def buildRcToSeqDict(seqToRepDict, bm):
        # Initialize the dict
        rcToSeqDict = dict()

        # For each sequence key in the dict,
        for seq in seqToRepDict.keys():
            # Compute the reverse complement
            rcBitSeq = bm.getReverseComplement(
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
        # If either the evolvability values have not been computed already,
        # or the caller has explicitly requested re-computing,
        if not self.evoTuples or recompute:
            # Get tuples: (evo, targetRepsDict) for all sequences in the network
            self.evoTuples = [
                self.getSeqEvo(seq)
                for seq in self.network.vs["sequences"]
            ]

        return self.evoTuples

    # Returns evolvability score for the given sequence
    def getSeqEvo(self, sequence):
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
