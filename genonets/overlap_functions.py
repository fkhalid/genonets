"""
    overlap_functions
    ~~~~~~~~~~~~~~~~~

    Contains functions used for computation of overlap.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""


class OverlapAnalyzer:
    # Constructor
    def __init__(self, repToGiantDict, repertoires):
        # Copy of dict: key=repertoire, value=giant
        self.repToGiantDict = repToGiantDict

        # List of available repertoires
        self.repertoires = repertoires

    def getOverlapData(self):
        if len(self.repertoires) < 2:
            print("Overlap computation triggered with only one repertoire!")
            print("Overlap can only be calculated with 2 or more repertoires.")

            return None, None

        # Overlap dict for all repertoires. Dict{rep : {seq : [target reps]}}
        allOverlap = {
            rep: {
                seq: []
                for seq in self.repToGiantDict[rep].vs["sequences"]
            }
            for rep in self.repertoires
        }

        # Initialize the overlap matrix with zeros
        overlapMat = [
            [0 for x in range(len(self.repertoires))]
            for x in range(len(self.repertoires))
        ]

        # For each repertoire,
        for i in range(len(self.repertoires) - 1):
            # Get giant for this repertoire
            giant_i = self.repToGiantDict[self.repertoires[i]]

            # Get sequence list from giant
            seqs_i = giant_i.vs["sequences"]

            # For the rest of the repertoires,
            for j in range(i + 1, len(self.repertoires)):
                # Get giant for this repertoire
                giant_j = self.repToGiantDict[self.repertoires[j]]

                # Get sequence list from giant
                seqs_j = giant_j.vs["sequences"]

                # Get the list of sequences that are common to both repertoires
                overlapList = self.getOverlapList(seqs_i, seqs_j)

                # Get the number of sequences that overlap between the
                # two repertoires, and set it in the matrix
                overlapMat[i][j] = overlapMat[j][i] = len(overlapList)

                # If there is at least one shared sequence,
                if len(overlapList) > 0:
                    # For each sequence in the list,
                    for sequence in overlapList:
                        # Append the target to the target list for this sequence
                        allOverlap[self.repertoires[i]][sequence].append(self.repertoires[j])
                        allOverlap[self.repertoires[j]][sequence].append(self.repertoires[i])

        return overlapMat, self.repertoires, allOverlap

    # Returns a list of sequences that are common to the two input sequence
    # lists
    def getOverlapList(self, seqs1, seqs2):
        return list(set(seqs1) & set(seqs2))
