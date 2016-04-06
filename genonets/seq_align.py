import numpy as np


class Aligner:
    def __init__(self, alphabet, seqlen, match, mismatch, gap):
        # Alphabet (or molecule type)
        self.alphabet = alphabet

        # Length of the sequences to align
        self.seqlen = seqlen

        # Scores
        self.MATCH = match
        self.MISMATCH = mismatch
        self.GAP = gap

        # Sequences to be aligned; populated later
        self.seq1 = None
        self.seq2 = None

        # Matrix: Initialized later
        self.mat = None

        # Map: holds information for backtracking
        self.btMap = None

        # Constants used in backtracking
        self.INDEL1 = 1
        self.INDEL2 = 2
        self.NOINDEL_MATCH = 3
        self.NOINDEL_MISMATCH = 4

        # Aligned sequences
        self.seq1_aligned = None
        self.seq2_aligned = None

        # Mismatch count
        self.numMismatches = None

    def align(self, seq1, seq2):
        # Store sequences as class members
        self.seq1 = seq1
        self.seq2 = seq2

        # Initialize the No. of mismatches
        self.numMismatches = 0

        # --------- Populate the matrix --------- #

        # Initialize info map for backtracking
        self.btMap = dict()

        # Initialize with zeroes
        self.mat = np.zeros(
            shape=(self.seqlen + 1, self.seqlen + 1)
        )

        # Fill in the first row
        for i in xrange(1, self.seqlen + 1):
            self.mat[0, i] = self.mat[0, i - 1] + self.GAP

            # # Store op in map
            # self.btMap[(0, i)] = (self.INDEL2, (0, i-1))

        # Fill in the first column
        for i in xrange(1, self.seqlen + 1):
            self.mat[i, 0] = self.mat[i - 1, 0] + self.GAP

            # # Store op in map
            # self.btMap[(i, 0)] = (self.INDEL1, (i-1, 0))

        # Fill in the rest of the matrix
        for i in xrange(1, self.seqlen + 1):
            for j in xrange(1, self.seqlen + 1):
                self.mat[i, j] = self.score(i, j)

        # --------- Backtrack the optimal path --------- #

        self.backtrack()

    def score(self, i, j):
        # Check if the letters match at the current position
        if self.seq1[i - 1] == self.seq2[j - 1]:
            match_score = self.MATCH
        else:
            match_score = self.MISMATCH

        # Case 1: S_i-1,j-1 + s(a_i, b_i)
        s1 = self.mat[i - 1, j - 1] + match_score

        if i in [1, self.seqlen] and j in [1, self.seqlen]:
            # Case 2: S_i-1,j - delta
            s2 = self.mat[i - 1, j] + self.GAP

            # Case 3: S_i,j-1 - delta
            s3 = self.mat[i, j - 1] + self.GAP
        else:
            s3 = s2 = s1 - 1

        # Compute the max
        maxval, indel, parent = self.max(s1, s2, s3, match_score, i=i, j=j)

        # Store op in map
        self.btMap[(i, j)] = (indel, parent)

        return maxval

    def max(self, s1, s2, s3, match_status, i, j):
        if s1 >= s2 and s1 >= s3:
            return s1, \
                   self.NOINDEL_MATCH if match_status == self.MATCH else self.NOINDEL_MISMATCH, \
                   (i - 1, j - 1)
        elif s2 >= s1 and s2 >= s3:
            return s2, self.INDEL1, (i - 1, j)
        else:
            return s3, self.INDEL2, (i, j - 1)

    def backtrack(self):
        # Aligned sequences
        self.seq1_aligned = ""
        self.seq2_aligned = ""

        # Starting cell index for backtracking
        i = self.seqlen
        j = self.seqlen

        try:
            while True:
                # Get information for this cell
                cellinfo = self.btMap[(i, j)]

                # If this is an indel,
                if cellinfo[0] in [self.NOINDEL_MATCH, self.NOINDEL_MISMATCH]:
                    self.seq1_aligned += self.seq1[i - 1]
                    self.seq2_aligned += self.seq2[j - 1]
                    i -= 1
                    j -= 1

                    if cellinfo[0] == self.NOINDEL_MISMATCH:
                        self.numMismatches += 1
                elif cellinfo[0] == self.INDEL1:
                    self.seq1_aligned += self.seq1[i - 1]
                    self.seq2_aligned += "-"
                    i -= 1
                    self.numMismatches += 1
                else:
                    self.seq1_aligned += "-"
                    self.seq2_aligned += self.seq2[j - 1]
                    j -= 1
                    self.numMismatches += 1
        except KeyError:
            # Done backtracking
            pass

    def revserse_complement(self, seq):
        # Dictionary that maps each base to its complement
        baseToComplementDict = {
            'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C'
        }

        reverse_seq = seq[::-1]

        rc = [
            baseToComplementDict[reverse_seq[i]]
            for i in xrange(len(reverse_seq))
        ]

        return ''.join(rc)

    def areNeighbors(self, seq1, seq2, bitmanip):
        areneighbors = False
        # error = False
        # complement = False

        self.align(seq1, seq2)

        if self.numMismatches == 1:
            areneighbors = True
        else:
            # complement = True
            self.align(self.seq1,
                       self.revserse_complement(seq2))

            if self.numMismatches == 1:
                areneighbors = True

        # if areneighbors:
        #     if not bitmanip.areNeighbors(bitmanip.seqToBits(seq1), bitmanip.seqToBits(seq2)):
        #         error = True
        # else:
        #     if bitmanip.areNeighbors(bitmanip.seqToBits(seq1), bitmanip.seqToBits(seq2)):
        #         error = True
        # if error:
        #         print("Mismatch:")
        #         print("Are neighbors: " + str(areneighbors))
        #         print("Reverse complement: " + str(complement))
        #         print("Sequence 1: " + seq1)
        #         print("Sequence 2: " + self.seq2)
        #         print("Sequence 2: " + self.seq2)
        #         print("RC Sequence 2: " + self.revserse_complement(self.seq2))
        #         print("Alignment: ")
        #         print("Sequence 1: " + self.seq1_aligned)
        #         print("Sequence 2: " + self.seq2_aligned)

        return areneighbors
