"""
    seq_bit_impl
    ~~~~~~~~~~~~

    Contains classes that define alphabet types that can be used with
    the bit manipulation library.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import math

from genonets_utils import Utils
from seq_bit_lib import AbstractBitSeqManipulator


class CustomBitSeqManipulator(AbstractBitSeqManipulator):
    # TODO: In order to support reverse complements, we need a mapping
    #       from each letter to its complement, which should be provided
    #       by the user ...
    def __init__(self, seqLength, useIndels, letter_to_neighbors):
        # No. of unique letters in the alphabet
        self._alphabet_size = len(set(letter_to_neighbors.keys()))

        # Mapping from letters to all 1-neighbors
        self._letter_to_neighbors = letter_to_neighbors

        # No. of bits required to encode the alphabet
        self.bitCodeLen = int(math.ceil(math.log(self._alphabet_size, 2)))

        # Mapping from letters to the corresponding bit representation
        self.letterToBitDict = dict(
            (letter, bit) for bit, letter in enumerate(letter_to_neighbors)
        )

        # Mapping from bit representation to the corresponding letter
        self.bitToLetterDict = Utils.reverseDict(self.letterToBitDict)

        # Call the base initializer
        AbstractBitSeqManipulator.__init__(self, seqLength, useIndels)

        # Mapping from letters to the bit representation of 1-neighbors
        self._letter_to_neighbor_bits = {}

        # For each letter,
        for letter in letter_to_neighbors:
            # For each letter in the letter to bit map,
            if letter in self.letterToBitDict:
                # Set the key to the bit representation of the letter, and
                # value to the set of bit representations of the 1-neighbor
                # letters.
                self._letter_to_neighbor_bits[self.letterToBitDict[letter]] = {
                    self.letterToBitDict[n] for n in letter_to_neighbors[letter]
                    if n in self.letterToBitDict    # Just so * is ignored
                }

    def _are_neighbors(self, seq1, seq2):
        # If the distance between seq1 and seq2 is a point
        # mutation, they are neighbors; no need to check further.
        if self.distanceBwSeqs(seq1, seq2) == 1:
            return 'point'

        # If shift mutations should not be considered,
        if not self.useIndels:
            return None

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
            return 'left'

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
            return 'right'
        else:
            return None

    def areNeighbors(self, seq1, seq2):
        # Check if seq1 and seq2 are 1-neighbors
        result = self._are_neighbors(seq1, seq2)

        if result == 'point':
            # Full string representation of the bit string where only the
            # bits that differ between seq1 and seq2 are set. These
            # would correspond to the letter at which the two
            # sequences differ.
            diff = "{0:{fill}{code_length}b}".format(
                seq1 ^ seq2,
                fill='0',
                code_length=self.seqLength * self.bitCodeLen
            )

            # Index of the first set bit. Only the first is enough to
            # determine which letter it is.
            index = diff.index('1')

            # Index of the letter at which the two sequences differ
            i = max(0, index / self.bitCodeLen)

            # Get the bit values for letter at index 'i' both the source
            # and mutated sequences
            bit_seq_1 = self.getLetterAtIndex(seq1, i)
            bit_seq_2 = self.getLetterAtIndex(seq2, i)

            if bit_seq_2 not in self._letter_to_neighbor_bits[bit_seq_1]:
                return False

            return True
        elif result in {'right', 'left'}:
            return True

        return False

    def generateNeighbors(self, sequence):
        # No. of nucleotides in the sequence
        k = self.seqLength

        # Generate all possible single letter mutants
        # TODO: Modify the following comprehension such that letter
        #  being mutated can only be mutated to one of the other
        #  letters allowed by the genetic code
        neighbors = [
            self.mutateLetter(sequence, i, target)
            for i in range(k)
            for target in self._letter_to_neighbor_bits[self.getLetterAtIndex(sequence, i)]
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

        return neighbors


class BinaryBitSeqManipulator(AbstractBitSeqManipulator):
    # No. of bits required to encode a single bit
    bitCodeLen = 1

    # Dictionary: Maps single bit letter to bit value
    # (represented as integers here)
    letterToBitDict = {'0': 0, '1': 1}

    # Dictionary: Maps bit values to bit letters
    bitToLetterDict = Utils.reverseDict(letterToBitDict)


class ProteinBitSeqManipulator(AbstractBitSeqManipulator):
    # No. of bits required to encode a single amino acid
    bitCodeLen = 5

    # Dictionary: Maps single letter amino acid codes to bit
    # values (represented as integers here)
    letterToBitDict = {
        'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4,
        'E': 5, 'Q': 6, 'G': 7, 'H': 8, 'I': 9,
        'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
        'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19
    }

    # Dictionary: Maps bit values to single letter amino acid codes
    bitToLetterDict = Utils.reverseDict(letterToBitDict)

    # # Constructor
    # def __init__(self, seqLength, useIndels, letter_to_neighbors):
    #     AbstractBitSeqManipulator.__init__(self, seqLength, useIndels)
    #
    #     self.letter_to_neighbors = letter_to_neighbors
    #
    #     self.letter_to_neighbor_bits = {}
    #
    #     for letter in letter_to_neighbors:
    #         if letter in self.letterToBitDict:
    #             self.letter_to_neighbor_bits[self.letterToBitDict[letter]] = {
    #                 self.letterToBitDict[n] for n in letter_to_neighbors[letter]
    #                 if n in self.letterToBitDict    # Just so * is ignored
    #             }

    # def _are_neighbors(self, seq1, seq2):
    #     # If the distance between seq1 and seq2 is a point
    #     # mutation, they are neighbors; no need to check further.
    #     if self.distanceBwSeqs(seq1, seq2) == 1:
    #         return 'point'
    #
    #     # If shift mutations should not be considered,
    #     if not self.useIndels:
    #         return None
    #
    #     # Since seq2 is not separated by a point mutation, test for
    #     # a left shift mutation
    #
    #     # Get the right most letter of seq2
    #     rmn = self.getLetterAtIndex(seq2, self.seqLength - 1)
    #
    #     # Left shift seq1 by one letter (n bits)
    #     temp = seq1 << self.bitCodeLen
    #
    #     # Mutate the right most nucleotide in temp and see if that results
    #     # in seq2. If so, seq1 and seq2 are 1-neighbors separated by
    #     # a single left shift mutation.
    #     if seq2 == self.mutAftrLftShift(temp, self.seqLength - 1, rmn):
    #         return 'left'
    #
    #     # Since seq2 is neither separated by a point mutation nor by
    #     # a left shit, test for a right shift mutation
    #
    #     # Get the left most nucleotide of seq2
    #     lmn = self.getLetterAtIndex(seq2, 0)
    #
    #     # Right shift seq1 by one nucleotide (2 bits)
    #     temp = seq1 >> self.bitCodeLen
    #
    #     # Mutate the right most nucleotide temp and see if that results
    #     # in seq2. If so, seq1 and seq2 are 1-neighbors separated by
    #     # a single right shift mutation. If not, we have exhausted all
    #     # possibilities, which means seq1 and seq2 are not 1-neighbors.
    #     if seq2 == self.mutateLetter(temp, 0, lmn):
    #         return 'right'
    #     else:
    #         return None

    # def _mutation_is_allowed(self, seq1, seq2, i):
    #     # Get the letter at index 'i'
    #     bitValue1 = self.getLetterAtIndex(seq1, i)
    #     bitValue2 = self.getLetterAtIndex(seq2, i)
    #
    #     # If there was a mutation here,
    #     if bitValue1 != bitValue2:
    #         letter1 = self.bitToLetterDict[bitValue1]
    #         letter2 = self.bitToLetterDict[bitValue2]
    #
    #         if letter2 not in self.letter_to_neighbors[letter1]:
    #             return False
    #
    #     return True

    # def areNeighbors(self, seq1, seq2):
    #     # Check if seq1 and seq2 are 1-neighbors
    #     result = self._are_neighbors(seq1, seq2)
    #
    #     if self.letter_to_neighbors and result == 'point':
    #         diff = "{0:{fill}{code_length}b}".format(
    #             seq1 ^ seq2, fill='0', code_length=self.seqLength * self.bitCodeLen)
    #         index = diff.index('1')
    #         i = max(0, index / self.bitCodeLen)
    #
    #         # Get the letter at index 'i'
    #         bitValue1 = self.getLetterAtIndex(seq1, i)
    #         bitValue2 = self.getLetterAtIndex(seq2, i)
    #
    #         # If there was a mutation here,
    #         if bitValue1 != bitValue2:
    #             letter1 = self.bitToLetterDict[bitValue1]
    #             letter2 = self.bitToLetterDict[bitValue2]
    #
    #             if letter2 not in self.letter_to_neighbors[letter1]:
    #                 return False
    #
    #         return True
    #         # Check if the mutation is allowed by the genetic code
    #         # for i in range(self.seqLength):
    #         #     # if not self._mutation_is_allowed(seq1, seq2, i):
    #         #     #     return False
    #         #     # Get the letter at index 'i'
    #         #     bitValue1 = self.getLetterAtIndex(seq1, i)
    #         #     bitValue2 = self.getLetterAtIndex(seq2, i)
    #         #
    #         #     # If there was a mutation here,
    #         #     if bitValue1 != bitValue2:
    #         #         letter1 = self.bitToLetterDict[bitValue1]
    #         #         letter2 = self.bitToLetterDict[bitValue2]
    #         #
    #         #         if letter2 not in self.letter_to_neighbors[letter1]:
    #         #             return False
    #     elif result in {'point', 'right', 'left'}:
    #         return True
    #
    #     return False
    #
    # def generateNeighbors(self, sequence):
    #     # No. of nucleotides in the sequence
    #     k = self.seqLength
    #
    #     # Generate all possible single letter mutants
    #     # TODO: Modify the following comprehension such that letter
    #     #  being mutated can only be mutated to one of the other
    #     #  letters allowed by the genetic code
    #     neighbors = [
    #         self.mutateLetter(sequence, i, target)
    #         for i in range(k)
    #         for target in self.letter_to_neighbor_bits[self.getLetterAtIndex(sequence, i)]
    #     ]
    #
    #     # If shift mutations should be considered,
    #     if self.useIndels:
    #         # Generate all possible left shift mutations
    #
    #         # Left shift sequence by one letter, i.e., n bits
    #         temp = sequence << self.bitCodeLen
    #
    #         # Generate all possible left shift mutants by mutating the right most
    #         # nucleotide 4 times.
    #         lsNeighbors = self.getShiftMutants("left", temp)
    #
    #         # Append the list of left shift mutants to the neighborhood. Avoid duplication,
    #         # which is possible for sequences where the same letter occupies all positions.
    #         # Also, do not include the mutant that results in the source sequence.
    #         neighbors.extend(self.getUniqueNeighbors(lsNeighbors, neighbors, sequence))
    #
    #         # Generate all right shift mutants
    #
    #         # Right shift sequence by one nucleotide (2 bits)
    #         temp = sequence >> self.bitCodeLen
    #
    #         # Generate all possible right shift mutations by mutating the right most
    #         # nucleotide 4 times. However, do not include the mutant that results
    #         # in the source sequence.
    #         rsNeighbors = self.getShiftMutants("right", temp)
    #
    #         # Append the list of right shift mutants to the neighborhood. Avoid duplication,
    #         # which is possible for sequences where the same letter occupies all positions.
    #         # Also, do not include the mutant that results in the source sequence.
    #         neighbors.extend(self.getUniqueNeighbors(rsNeighbors, neighbors, sequence))
    #
    #     return neighbors

# Abstract base class the provides common attributes for RNA
# DNA sequence manipulation as bits. Only derived concrete
# classes are supposed to be instantiated.
class AbstractXnaBitSeqManipulator(AbstractBitSeqManipulator):
    # No. of bits required to encode a single nucleotide
    bitCodeLen = 2


# Concrete class for manipulation of arbitrary length RNA
# sequences as bits.
class RNABitSeqManipulator(AbstractXnaBitSeqManipulator):
    # Dictionary that maps nucleotide letters to bit values
    # (represented as integers here)
    letterToBitDict = {'A': 0, 'U': 1, 'C': 2, 'G': 3}

    # Dictionary that maps bit values to nucleotide letters
    bitToLetterDict = {0: 'A', 1: 'U', 2: 'C', 3: 'G'}


# Concrete class for manipulation of arbitrary length DNA
# sequences as bits.
class DNABitSeqManipulator(AbstractXnaBitSeqManipulator):
    # Dictionary that maps nucleotide letters to bit values
    # (represented as integers here)
    letterToBitDict = {'A': 0, 'T': 1, 'C': 2, 'G': 3}

    # Dictionary that maps bit values to nucleotide letters
    bitToLetterDict = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}

    # Dictionary that maps each base to its complement
    baseToComplementDict = {0: 1, 1: 0, 2: 3, 3: 2}

    # Constructor
    def __init__(self, seqLength, useIndels, use_reverse_complements):
        AbstractXnaBitSeqManipulator.__init__(self, seqLength, useIndels)

        self.useRC = use_reverse_complements

    # Computes and returns the bit representation of the reverse
    # complement of the given sequence
    def getReverseComplement(self, sequence):
        # Variable to hold the reverse complement
        revCompl = 0

        # Iterate over nucleotides in the given sequence
        # in reverse order
        for i in reversed(range(self.seqLength)):
            # Get nucleotide at index 'i'
            bitValue = self.getLetterAtIndex(sequence, i)

            # Add the complement base to the reverse complement
            revCompl |= self.baseToComplementDict[bitValue]

            # If this is not the final nucleotide to be processed,
            # i.e., the first nucleotide in the given sequence
            if i != 0:
                # Left shift by 2 bits to make room for the next
                # 2-bit sequence
                revCompl <<= 2

        return revCompl

    # Overrides the base implementation since for DNA, both
    # the sequence and its reverse complement have to be
    # considered.
    def areNeighbors(self, seq1, seq2):
        # Check if seq1 and seq2 are 1-neighbors
        if AbstractXnaBitSeqManipulator.areNeighbors(self, seq1, seq2):
            # If the two sequences are 1-neighbors, no need to
            # process any further
            return True

        # If reverse complements should be considered,
        if self.useRC:
            # Get the reverse complement of seq2
            revComplSeq2 = self.getReverseComplement(seq2)

            # Check if seq1 and reverse complement of seq2 are
            # 1-neighbors
            if AbstractXnaBitSeqManipulator.areNeighbors(self, seq1, revComplSeq2):
                return True

        # At this point, neither possibility is true, so we return False.
        return False


# ------------------------------------------------------------
# Factory method
# ------------------------------------------------------------

# TODO: Could probably do something similar to AnalysisHandler here, i.e.,
#		use a more elegant method for determining the BitSeqManipulator
#		corresponding to the given molecule type. E.g., using a dict,
#		{moleculeType : constructor} ...
class BitManipFactory:
    # Supported molecule types
    moleculeTypes = ["RNA", "DNA", "Protein", "Binary"]

    # Returns the list of supported molecule types
    @staticmethod
    def getMoleculeTypes():
        return BitManipFactory.moleculeTypes

    # Based on the given molecule type, returns the corresponding
    # bit seq manipulator object.
    @staticmethod
    def getBitSeqManip(moleculeType, seqLength, useIndels, letter_to_neighbors, useReverseComplements=False):
        if letter_to_neighbors:
            return CustomBitSeqManipulator(seqLength, useIndels, letter_to_neighbors)

        if moleculeType == "RNA":
            return RNABitSeqManipulator(seqLength, useIndels)
        elif moleculeType == "DNA":
            return DNABitSeqManipulator(seqLength, useIndels,  useReverseComplements)
        elif moleculeType == "Protein":
            # return ProteinBitSeqManipulator(seqLength, useIndels, letter_to_neighbors)
            return ProteinBitSeqManipulator(seqLength, useIndels)
        elif moleculeType == "Binary":
            return BinaryBitSeqManipulator(seqLength, useIndels)
        else:
            print("Unsupported moleculeType: " + str(moleculeType))
            print("Exiting program ...")

            exit()
