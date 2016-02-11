
"""
    evolvability_functions
    ~~~~~~~~~~~~~~~~~~~~~~

    Contains functions used for genotype and phenotype evolvability computations.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

class EvolvabilityAnalyzer :
    # Constructor
    def __init__(self, network, dataDict, seqToRepDict, repToGiantDict, netBuilder, bitsToSeqDict) :
        # Reference to the network on which to perform this
        # analysis
        self.network = network

        # Reference to dict: key=rep, value=dict{seq: score}
        self.dataDict = dataDict

        # Copy dict: key=sequence, value=[repertoires]
        self.seqToRepDict = seqToRepDict.copy()

        # Copy of dict: key=repertoire, value=giant
        self.repToGiantDict = repToGiantDict

        # Remove all seq-rep mappings, where seq is not within the
        # giant
        self.updateSeqToRepDict()

        # Get a reference to the NetworkUtils object
        self.netBuilder = netBuilder

        # Dict {bitseq : seq}
        self.bitsToSeqDict = bitsToSeqDict

        # Reference to the list of evlvability values for all sequences
        self.evoTuples = None

    def updateSeqToRepDict(self) :
        # For each sequence key in the dict,
        for seq in self.seqToRepDict.keys() :
            # For each repertoire of which the sequence is a member,
            for rep in self.seqToRepDict[seq] :
                # If the sequence is not within the giant of the
                # repertoire,
                if seq not in self.repToGiantDict[rep].vs["sequences"] :
                    # Remove the repertoire from the list
                    self.seqToRepDict[seq].remove(rep)

    # Returns repertorire evolvability
    def getReportoireEvo(self) :
        # Get evo tuples for all sequences in the network
        evoTuples = self.getEvoAll()

        # Get repertoire lists from the tuples
        repLists = [evoTuples[1].keys() \
                        for evoTuples in evoTuples if evoTuples[1]]

        # Combine all repLists into one list
        targets = []
        for repList in repLists :
            targets.extend([target for target in repList])

        # Remove redundant repertoires
        targets = list(set(targets))

        try :
            evolvability = float(len(targets)) / float(len(self.dataDict) - 1)

            return evolvability, targets
        except ZeroDivisionError :
            return 0, targets

    def getEvoAll(self, recompute=False) :
        # If either the evolvability values have not been computed already,
        # or the caller has explicitly requested re-computing,
        if not self.evoTuples or recompute :
            # Get tuples: (evo, targetRepsDict) for all sequences in the network
            self.evoTuples = [self.getSeqEvo(seq) \
                                for seq in self.network.vs["sequences"]]

        return self.evoTuples

    # Returns evolvability score for the given sequence
    def getSeqEvo(self, sequence) :
        # Get a list of sequences that are 1-neighbors, but not part of the
        # genotype network
        externNeighbors = self.netBuilder.getExternalNeighbors(sequence, self.network)

        # Get dictionary: key=repertoire, value=extNeighs. The keys represents the
        # list of repertoires to which the given sequence could join by undergoing
        # a single mutation.
        targetReps = self.getEvoTargetReps(sequence, externNeighbors)

        # Compute and return the evolvability score, i.e., No. of repertoire
        # targets found for this sequence / No. of repertoires available in
        # the dataset.
        # FIXME: ZeroDivisionError can only occur if the denominator is zero. This
        #		 already be checked at the very beginning of the function to avoid
        #		 the rest of the computation ...
        try :
            evolvability = float(len(targetReps)) / float(len(self.dataDict) - 1)
        except ZeroDivisionError :
            evolvability = 0

        return (evolvability, targetReps)

    # For the given sequence, list of external neighbors, and the sequence
    # to repertoire dictionary, return a list of repertoires for which the given
    # sequence could develop score >= tau by a single mutation.
    def getEvoTargetReps(self, sequence, extNeighs) :
        targetReps = {}

        # For each external neighbor,
        for extNeigh in extNeighs :
            try :
                # Get the sequence corresponding to the bitseq
                extNeighSeq = self.bitsToSeqDict[extNeigh]
            except KeyError :
                # If the bitseq corresponding to sequence could not be found,
                # it means that it does not exist in any other repertoire. We
                # can therefore safely skip the current iteration
                continue

            # If the external neighbor is in any repertoires,
            if extNeighSeq in self.seqToRepDict :
                # For each repertoire that contains the external neighbor,
                for repertoire in self.seqToRepDict[extNeighSeq] :
                    # We have found a repertoire with a genotype network that
                    # has a 1-neighbor to which  the sequence can evolve.
                    if repertoire not in targetReps :
                        targetReps[repertoire] = []

                    targetReps[repertoire].append(extNeighSeq)

        return targetReps
