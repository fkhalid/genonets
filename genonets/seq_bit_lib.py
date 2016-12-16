"""
    seq_bit_lib
    ~~~~~~~~~~~

    Library of functions for manipulation of genotypes as bits.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

from array import array


# Abstract base class that provides all methods for manipulation of arbitrary
# character sequences as bits. Only concrete classes derived from this class
# are supposed to be instantiated for use.
class AbstractBitSeqManipulator:
    # Constructor
    def __init__(self, seqLength, useIndels):
        # Set the sequence length
        self.seqLength = seqLength

        # Flag to indicate whether shift mutations should be considered
        self.useIndels = useIndels

        # Compute and store the start indices, where each index
        # represents the start bit of an 'n'-bit code
        self.startIndices = self.computeStartIndices()

        # Generate and store masks to be used in distance
        # computation
        self.masks = self.getMasks()

    # Basic popcount
    def popcount(self, x):
        return bin(x).count('1')

    # Returns the bit sequence corresponding to the given letter
    def letterToBits(self, letter):
        return self.letterToBitDict[letter]

    # Return the letter corresponding to the given bit sequence
    def bitsToLetter(self, bits):
        return self.bitToLetterDict[bits]

    # Returns 'True' if the given number is even and 'False'
    # otherwise.
    def isEven(self, x):
        return (x % 2) == 0

    # Returns a list of indices, such that each index in the list
    # is the index of the start of bit of the n-bit code
    # representing a single letter in the sequence.
    def computeStartIndices(self):
        # Code length
        n = self.bitCodeLen
        # Total No. of bits in the sequence
        l = n * self.seqLength

        # Create a list of start indices based on the sequence
        # length
        indices = [x for x in range(n - 1, l, n)]

        return indices

    # Returns 'n' masks, where 'n = bitCodeLen'
    def getMasks(self):
        masks = []

        for i in range(self.bitCodeLen):
            if i == 0:
                # Zeros mask is no use here. We start with a mask
                # with the first bit set in each letter.
                masks.append(self.generateMask(1))
            else:
                # Mask with bit i-1 set in each letter. This is
                # done by taking the i-1 power of 2.
                masks.append(self.generateMask(2 << (i - 1)))

        return masks

    # Generates and returns a mask where only every pth bit is set,
    # where 'p' is determined by 'maskVal'.
    def generateMask(self, maskVal):
        # Total No. of bits in the sequence
        l = self.bitCodeLen * self.seqLength

        # Initialize mask to be constructed later
        mask = 0

        # Starting from the most significant bit, set maskval
        # in each nth position, where n=self.bitCodeLen
        for i in reversed(range(0, l, self.bitCodeLen)):
            mask |= maskVal << i

        return mask

    # Return True if the two sequences are 1-neighbors.
    # The 1-neighborhood includes all single point mutations, as well as all
    # left shift and right shift mutations.
    def areNeighbors(self, seq1, seq2):
        # If the distance between seq1 and seq2 is a point
        # mutation, they are neighbors; no need to check further.
        if self.distanceBwSeqs(seq1, seq2) == 1:
            return True

        # If shift mutations should not be considered,
        if not self.useIndels:
            return False

        # Since seq2 is not separated by a point mutation, test for
        # a left shift mutation

        # Get the right most letter of seq2
        rmn = self.getLetterAtIndex(seq2, self.seqLength - 1)

        # Left shift seq1 by one letter (n bits)
        temp = seq1 << self.bitCodeLen

        # Mutate the right most nucleotide in temp and see if that results
        # in seq2. If so, seq1 and seq2 are 1-neighbors separated by
        # a single left shift mutation.
        if seq2 == self.mutAftrLftShift(temp, self.seqLength - 1, rmn):
            return True

        # Since seq2 is neither separated by a point mutation nor by
        # a left shit, test for a right shift mutation

        # Get the left most nucleotide of seq2
        lmn = self.getLetterAtIndex(seq2, 0)

        # Right shift seq1 by one nucleotide (2 bits)
        temp = seq1 >> self.bitCodeLen

        # Mutate the right most nucleotide temp and see if that results
        # in seq2. If so, seq1 and seq2 are 1-neighbors separated by
        # a single right shift mutation. If not, we have exhausted all
        # possibilities, which means seq1 and seq2 are not 1-neighbors.
        if seq2 == self.mutateLetter(temp, 0, lmn):
            return True
        else:
            return False

    # Accepts a sequence in bit format and generates all possible 1-neighbors.
    # Returns a list of neighbors.
    def generateNeighbors(self, sequence):
        # No. of nucleotides in the sequence
        k = self.seqLength

        # Generate all possible single letter mutants
        neighbors = [
            self.mutateLetter(sequence, i, target)
            for i in range(k)
            for target in self.bitToLetterDict.keys()
            if target != self.getLetterAtIndex(sequence, i)
        ]

        # If shift mutations should be considered,
        if self.useIndels:
            # Generate all possible left shift mutations

            # Left shift sequence by one letter, i.e., n bits
            temp = sequence << self.bitCodeLen

            # Generate all possible left shift mutants by mutating the right most
            # nucleotide 4 times.
            lsNeighbors = self.getShiftMutants("left", temp)

            # Append the list of left shift mutants to the neighborhood. Avoid duplication,
            # which is possible for sequences where the same letter occupies all positions.
            # Also, do not include the mutant that results in the source sequence.
            neighbors.extend(self.getUniqueNeighbors(lsNeighbors, neighbors, sequence))

            # Generate all right shift mutants

            # Right shift sequence by one nucleotide (2 bits)
            temp = sequence >> self.bitCodeLen

            # Generate all possible right shift mutations by mutating the right most
            # nucleotide 4 times. However, do not include the mutant that results
            # in the source sequence.
            rsNeighbors = self.getShiftMutants("right", temp)

            # Append the list of right shift mutants to the neighborhood. Avoid duplication,
            # which is possible for sequences where the same letter occupies all positions.
            # Also, do not include the mutant that results in the source sequence.
            neighbors.extend(self.getUniqueNeighbors(rsNeighbors, neighbors, sequence))

        # # TODO:
        # #   The sequences for which the reverse complement is the same as the sequence,
        # #   a shift mutation can result in the source sequence being added to the list
        # #   of neighbors. However, these are already eliminated in the previous code
        # #   block. Then, is the following condition really necessary??
        # if self.useRC:
        #     # If the reverse complement of the given sequence has been
        #     # added as a neighbor, remove it
        #     try:
        #         neighbors.remove(self.getReverseComplement(sequence))
        #     except ValueError:
        #         pass

        return neighbors

    # Depending on the shiftType, i.e., "left" or "right", returns the
    # corresponding shift mutants.
    def getShiftMutants(self, shiftType, source):
        # Check if this is a left shift
        if shiftType == "left":
            # The letter index at which to perform the mutations is
            # 'seqLength - 1' for a left shift
            index = self.seqLength - 1

            # Get all left shift mutants
            mutants = [
                self.mutAftrLftShift(source, index, target)
                for target in self.bitToLetterDict.keys()
            ]
        else:  # Right shift
            # The letter index at which to perform the mutations is
            # '0' for a left shift
            index = 0

            # Get all right shift mutants
            mutants = [
                self.mutateLetter(source, index, target)
                for target in self.bitToLetterDict.keys()
            ]

        return mutants

    # Return a list of neighbors from 'newNeighbors' that are not already in
    # 'neighbors', and none of them is the same as the 'source' sequence for
    # which the neighborhood is being generated.
    def getUniqueNeighbors(self, newNeighbors, neighbors, source):
        return [
            neighbor
            for neighbor in newNeighbors
            if neighbor not in neighbors and
            neighbor != source
        ]

    # Computes the distance between two sequences (in bit format) in terms
    # of the total number of letters that differ in the two sequences.
    # FIXME: Should be extend to include distance for shift mutations ...
    def distanceBwSeqs(self, seq1, seq2):
        # Code length
        n = self.bitCodeLen

        # Get the list of numbers that each represent diff bits for a
        # particular position between 0 and n-1.
        diffList = self.buildDiffList(seq1, seq2)

        # Since each letter is represented as an n-bit sequence, we need
        # to make sure we count the difference in a single letter only once.
        # In order to do this, within a single letter, all bits other than
        # the left most bit are shifted to the position of the left most
        # bit. This means that any differences that existed only in the
        # other bits are now available in the left most positions as well.
        allOred = 0
        # For each number in the list
        for i in range(len(diffList)):
            # Move all differences to the left most position
            allOred |= diffList[i] << (n - i)

        # Popcount gives the distance
        distance = self.popcount(allOred)

        return distance

    # Builds a list of numbers, where each number has bits set only at
    # specific indices, and represents difference in bit value between
    # seq1 and seq2 at that index. These indices appear at intervals of
    # length 'n', where 'n = bitCodeLen'. Each number represents the bits
    # for each letter that differs on a certain index within the bit code.
    # E.g., one of the numbers contains only the difference bits on the first
    # index of each letter.
    def buildDiffList(self, seq1, seq2):
        diffList = []

        # Get the difference values for all bits. XOR gives us only those
        # bits that have different values in the two sequences.
        diffBits = seq1 ^ seq2  # Bitwise XOR

        # Populate the diffLists corresponding to each position
        for i in range(self.bitCodeLen):
            diffList.append(diffBits & self.masks[i])

        return diffList

    # Assuming the requested single letter mutation is being performed
    # on a sequence that has been shifted left, returns a mutant with the
    # letter at index 'position' in 'sequence' replaced with 'target'.
    # Note: This function is needed because for sequences that start with
    # 		a letter other than '0', a left shift can leave one or more
    #		'1' bits on the left. These must be unset so that integer
    #		comparisons return correct results.
    def mutAftrLftShift(self, sequence, position, target):
        # Mutate the required nucleotide
        mutant = self.mutateLetter(sequence, position, target)

        # Unset any possible set bits on the left at position -1
        mutant = self.mutateLetter(mutant, -1, 0)

        return mutant

    # Mutates a single letter in the bit sequence at index=position, to
    # value=target. Position is expected to have been determined using
    # zero-based indexing.
    def mutateLetter(self, sequence, position, target):
        # Code length
        n = self.bitCodeLen
        # Total No. of bits in the sequence
        l = n * self.seqLength

        # Get a mask that is all ones of length=n
        ones = (1 << n) - 1

        # Get index of the right bit in the n-bit sequence that
        # corresponds to the given position.
        bitIndex = (l - n) - (n * position)

        # Create an n-bit mask, i.e., all ones ('11...1') at the
        # index computed above
        mask = ones << bitIndex

        # Perform the mutation
        # TODO: Add description here ...
        mutatedSeq = sequence & (~mask) | (target << bitIndex)

        return mutatedSeq

    # Return the n-bit code at the given index, where index is in the
    # range '0' to 'k - 1', i.e, it is the index corresponding to the
    # letter index in the string sequence (increasing from left to right).
    def getLetterAtIndex(self, sequence, index):
        # Total No. of bits in the sequence
        l = self.bitCodeLen * self.seqLength

        # Calculate the index that corresponds to the to the start bit
        # for the letter on the given index. E.g., if k=8, n=5, l=40,
        # index=7: l-1=39, index*n=35, => i=4. Then the letter is
        # represented by bits on indices 4,3,2,1,0.
        i = (l - 1) - (index * self.bitCodeLen)

        # Get the letter on index 'i'
        letter = self.getbitsOnIndex(i, sequence)

        return letter

    # Convert a sequence into a bit representation
    def seqToBits(self, sequence):
        # Initializing the to-be-constructed bit sequence
        bitSequence = 0

        # For each letter in the sequence
        for i in range(0, self.seqLength):
            # Get the corresponding bit code
            bitValue = self.letterToBitDict[sequence[i]]

            # Add the bit code to the bit sequence
            bitSequence |= bitValue

            # If this is not the last letter in the sequence
            if i != (self.seqLength - 1):
                # Left shift by 'bitCodeLen' bits to make room
                # for the bit code corresponding to the next
                # letter in the sequence
                bitSequence <<= self.bitCodeLen

        return bitSequence

    # Convert a bit sequence into letter sequence
    def bitsToSeq(self, bitSequence):
        # k-mer length
        k = self.seqLength

        # String to store the sequence. The initialization here has two
        # purposes: 1) To get the desired size, 2) To ensure there are
        # no ones, since left shifting ones at the left most index can
        # cause problems.
        outputSequence = array('c', ['0' for _ in range(k)])

        # For each n-bit sequence (in a k-mer occupying n*k bits in an
        # integer),
        j = k - 1
        for i in self.startIndices:
            # Get the letter on index 'i'
            letter = self.getbitsOnIndex(i, bitSequence)

            # Convert the n-bit sequence into a letter and store
            # it in the correct index of the output sequence
            outputSequence[j] = self.bitToLetterDict[letter]

            j -= 1

        return outputSequence.tostring()

    # Extracts 'n' bits from 'bitSeq' starting at index 'i'.
    # Returns an integer that is a copy of only these bits.
    def getbitsOnIndex(self, i, bitSeq):
        # Code length
        n = self.bitCodeLen

        # Extract the ith bit
        letter = ((bitSeq & (1 << i)) >> i)

        # For the rest of the bits used to represent this letter
        for p in range(1, n):
            # Shift the bit by 1 bit to the left to make room
            letter <<= 1

            # Extract pth bit for the current letter
            letter |= ((bitSeq & (1 << (i - p))) >> (i - p))

        return letter
