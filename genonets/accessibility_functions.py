"""
    accessibility_functions
    ~~~~~~~~~~~~~~~~~~~~~~~

    Exposes functions for accessibility, neighbor abundance, and diversity index computation.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""


class AccessibilityAnalyzer:
    # Constructor
    def __init__(self, repertoire, network, repToGiantDict, dataDict,
                 netBuilder, bitManip, isDoubleStranded):
        # Repertoire name
        self.repertoire = repertoire

        # Reference to the network on which to perform this
        # analysis
        self.network = network

        # Reference to dict: key=repertoire, value=network
        # FIXME: perhaps this dict should be replaced by the genonet
        #		 method getGraphFor() ...
        self.repToGiantDict = repToGiantDict

        # Reference to dict: key=rep, value=dict{seq: score}
        self.dataDict = dataDict

        # Get a reference to the NetworkBuilder object
        self.netBuilder = netBuilder

        # Reference to the bit manipulator object
        self.bitManip = bitManip

        # Flag to indicate whether reverse complements should be
        # considered for genotypes.
        self.isDoubelStranded = isDoubleStranded

    # Computes the accessibility value for the given repertoire.
    # Ref: Cowperthwaite et al. (2008) PLoS Comp. Biol.
    def getAccessibility(self):
        # Value to be calculated and returned
        accessibility = 0

        # Get the list of targets
        targets = self.dataDict.keys()

        # Remove the repertoire under consideration from the list of targets
        targets.remove(self.repertoire)

        # Get all sequences (in bit format) for the given repertoire
        sequences = [
            self.bitManip.seqToBits(seq) for seq in
            self.network.vs["sequences"]
        ]

        # If reverse complements should be considered,
        if self.isDoubelStranded:
            # Extend the list of sequences with their reverse complements
            sequences.extend([
                self.bitManip.getReverseComplement(seq)
                for seq in sequences
            ])

        # For each target, compute F(j,i)
        for target in targets:
            # Get the genotype network for this target
            targetNet = self.repToGiantDict[target]

            # Get all external neighbors for this target
            extNeighbors = self.netBuilder.getAllExtNeighbors(targetNet)

            # Find sequences that exist in both extNeighbors and
            # sequences, i.e., intersection of both sets
            commonSeqs = set(extNeighbors) & set(sequences)

            # Calculate what fraction of sequences in extNeighbors
            # is in commonSeqs
            try:
                fraction = float(len(commonSeqs)) / float(len(extNeighbors))
            except ZeroDivisionError:
                fraction = 0

            # Add the fraction to the total accessibility value
            accessibility += fraction

        return accessibility

    # Computes the neighbor abundance for the given repertoire.
    # Ref: Cowperthwaite et al. (2008) PLoS Comp. Biol.
    def getNeighborAbundance(self):
        # Value to be calculated and returned
        abundance = 0

        # Get all external neighbors for the given network
        extNeighbors = self.netBuilder.getAllExtNeighbors(self.network)

        # Get the list of targets
        targets = self.dataDict.keys()

        # Remove the repertoire under consideration from the list of targets
        targets.remove(self.repertoire)

        # For each target, compute F(i, j)
        for target in targets:
            # Get target sequences in bit format
            targetSeqs = [
                self.bitManip.seqToBits(seq) for seq in
                self.repToGiantDict[target].vs["sequences"]
            ]

            # If reverse complements should be considered,
            if self.isDoubelStranded:
                # Extend the list of sequences with their reverse complements
                targetSeqs.extend([
                    self.bitManip.getReverseComplement(seq)
                    for seq in targetSeqs
                ])

            # Find sequences that exist in both extNeighbors and
            # targetSeqs, i.e., intersection of both sets
            commonSeqs = set(extNeighbors) & set(targetSeqs)

            # Calculate what fraction of targetSeqs in extNeighbors
            # is in commonSeqs
            try:
                fraction = float(len(commonSeqs)) / float(len(extNeighbors))
            except ZeroDivisionError:
                fraction = 0

            # Multiply the fraction with the target's abundance, i.e., the size
            # of the target genotype network
            abundanceFrac = float(fraction) * float(len(targetSeqs))

            # Add the fraction to the total abundance value
            abundance += abundanceFrac

        return abundance

    # Computes the E-statistic value for the given repertoire.
    # Ref: Cowperthwaite et al. (2008) PLoS Comp. Biol.
    def getPhenotypicDivesity(self):
        # Value to be calculated and returned
        diversity = 0

        # Get all external neighbors for the given network
        extNeighbors = self.netBuilder.getAllExtNeighbors(self.network)

        # Get the list of targets
        targets = self.dataDict.keys()

        # Remove the repertoire under consideration from the list of targets
        targets.remove(self.repertoire)

        # For each target, compute F(i, j)
        for target in targets:
            # Get target sequences in bit format
            targetSeqs = [
                self.bitManip.seqToBits(seq) for seq in
                self.repToGiantDict[target].vs["sequences"]
            ]

            # If reverse complements should be considered,
            if self.isDoubelStranded:
                # Extend the list of sequences with their reverse complements
                targetSeqs.extend([
                    self.bitManip.getReverseComplement(seq)
                    for seq in targetSeqs
                ])

            # Find sequences that exist in both extNeighbors and
            # targetSeqs, i.e., intersection of both sets
            commonSeqs = set(extNeighbors) & set(targetSeqs)

            # Calculate what fraction of targetSeqs in extNeighbors
            # is in commonSeqs
            try:
                fraction = float(len(commonSeqs)) / float(len(extNeighbors))
            except ZeroDivisionError:
                fraction = 0

            # Add the fraction to the total diversity value
            diversity += (fraction * fraction)

        return diversity
