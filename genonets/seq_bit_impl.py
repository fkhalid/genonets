"""
    seq_bit_impl
    ~~~~~~~~~~~~

    Contains classes that define alphabet types that can be used with
    the bit manipulation library.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import math

from genonets.utils import Utils
from genonets.seq_bit_lib import AbstractBitSeqManipulator
from genonets.errors import InternalError


class GenericBitSeqManipulator(AbstractBitSeqManipulator):
    """
    Note: Custom genetic code handling is only permissible if the molecule
    type is Protein.

    """
    def __init__(
            self,
            sequence_length,
            use_indels,
            molecule_type,
            alphabet=None,
            letter_to_neighbors=None,
            use_reverse_complements=False):

        if alphabet is None and letter_to_neighbors is None:
            raise InternalError(
                'Creation of a Bit-sequence manipulator was attempted without '
                'providing the alphabet.'
            )

        # Mapping from letters to all 1-neighbors
        self._letter_to_neighbors = letter_to_neighbors

        if self._letter_to_neighbors is not None:
            self._alphabet = list(self._letter_to_neighbors.keys())
        else:
            self._alphabet = alphabet

        self._alphabet_size = len(self._alphabet)

        # No. of bits required to encode the alphabet
        self.bitCodeLen = int(math.ceil(math.log(self._alphabet_size, 2)))

        # Mapping from letters to the corresponding bit representation
        self.letterToBitDict = dict(
            (letter, bit)
            for bit, letter in enumerate(self._alphabet)
        )

        # Mapping from bit representation to the corresponding letter
        self.bitToLetterDict = Utils.reverseDict(self.letterToBitDict)

        if molecule_type == 'DNA':
            # Dictionary that maps each base to its complement
            self._base_to_complement_map = {0: 1, 1: 0, 2: 3, 3: 2}

            # Flag to indicate whether or not to use reverse complements
            self._use_reverse_complements = use_reverse_complements
        else:
            # Reverse complements are only available for DNA
            self._use_reverse_complements = False

        # Call the base initializer
        AbstractBitSeqManipulator.__init__(self, sequence_length, use_indels)

        # Mapping from letters to the bit representation of 1-neighbors
        self._letter_to_neighbor_bits = None

        if self._letter_to_neighbors is not None:
            self._letter_to_neighbor_bits = {}

            # For each letter,
            for letter in self._letter_to_neighbors:
                # Letter in numeric bit representation
                bit_letter = self.letterToBitDict[letter]

                # Set the key to the bit representation of the letter, and
                # value to the set of bit representations of the 1-neighbor
                # letters.
                self._letter_to_neighbor_bits[bit_letter] = {
                    self.letterToBitDict[n]
                    for n in letter_to_neighbors[letter]
                }

    # Computes and returns the bit representation of the reverse complement of
    # the given bit sequence
    def get_reverse_complement(self, sequence):
        # Variable to hold the reverse complement
        reverse_complement = 0

        # Iterate over nucleotides in the given sequence in reverse order
        for i in reversed(range(self.seqLength)):
            # Get nucleotide at index 'i'
            bit_value = self.getLetterAtIndex(sequence, i)

            # Add the complement base to the reverse complement
            reverse_complement |= self._base_to_complement_map[bit_value]

            # If this is not the final nucleotide to be processed, i.e., the
            # first nucleotide in the given sequence
            if i != 0:
                # Left shift by 2 bits to make room for the next 2-bit sequence
                reverse_complement <<= 2

        return reverse_complement

    def _are_neighbors(self, seq1, seq2):
        if self.distanceBwSeqs(seq1, seq2) == 1:
            return 'point'
        elif self.useIndels:
            # Get the right most letter of seq2
            rmn = self.getLetterAtIndex(seq2, self.seqLength - 1)

            # Left shift seq1 by one letter (n bits)
            temp = seq1 << self.bitCodeLen

            # Mutate the right most letter in temp and see if that results in
            # seq2. If so, seq1 and seq2 are 1-neighbors separated by a single
            # left shift mutation.
            if seq2 == self.mutAftrLftShift(temp, self.seqLength - 1, rmn):
                return 'left'

            # Since seq2 is not separated by a left shit, test for a right shift
            # mutation

            # Get the left most letter of seq2
            lmn = self.getLetterAtIndex(seq2, 0)

            # Right shift seq1 by one letter
            temp = seq1 >> self.bitCodeLen

            # Mutate the right most letter temp and see if that results in seq2.
            # If so, seq1 and seq2 are 1-neighbors separated by a single right
            # shift mutation.
            if seq2 == self.mutateLetter(temp, 0, lmn):
                return 'right'
        else:
            return None

    def areNeighbors(self, seq1, seq2):
        # Check if seq1 and seq2 are 1-neighbors
        result = self._are_neighbors(seq1, seq2)

        if self._letter_to_neighbor_bits is not None:
            if result == 'point':
                # Full string representation of the bit string where only the
                # bits that differ between seq1 and seq2 are set. These would
                # correspond to the letter at which the two sequences differ.
                diff = "{0:{fill}{code_length}b}".format(
                    seq1 ^ seq2,
                    fill='0',
                    code_length=self.seqLength * self.bitCodeLen
                )

                # Index of the first set bit. Only the first is enough to
                # determine which letter it is
                index = diff.index('1')

                # Index of the letter at which the two sequences differ
                i = max(0, int(index / self.bitCodeLen))

                # Get the bit values for letter at index 'i' in both the source
                # and mutated sequences
                bit_seq_1 = self.getLetterAtIndex(seq1, i)
                bit_seq_2 = self.getLetterAtIndex(seq2, i)

                if bit_seq_2 not in self._letter_to_neighbor_bits[bit_seq_1]:
                    return False

                return True
            elif result in {'right', 'left'}:
                return True
            else:
                return False
        else:
            if result is not None:
                return True
            else:
                if self._use_reverse_complements:
                    rev_compl_seq2 = self.get_reverse_complement(seq2)

                    result = self._are_neighbors(seq1, rev_compl_seq2)

                    return True if result is not None else False

        return False

    def generateNeighbors(self, sequence):
        # No. of letters in the sequence
        k = self.seqLength

        # Generate all possible single letter mutants
        if self._letter_to_neighbor_bits is not None:
            # TODO: Modify the following comprehension such that the letter
            #  being mutated can only be mutated to one of the other letters
            #  allowed by the genetic code
            neighbors = [
                self.mutateLetter(sequence, i, target)
                for i in range(k)
                for target in self._letter_to_neighbor_bits[
                    self.getLetterAtIndex(sequence, i)
                ]
            ]
        else:
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

            # Generate all possible left shift mutants by mutating the right
            # most nucleotide 4 times.
            lsNeighbors = self.getShiftMutants("left", temp)

            # Append the list of left shift mutants to the neighborhood. Avoid
            # duplication, which is possible for sequences where the same letter
            # occupies all positions. Also, do not include the mutant that
            # results in the source sequence.
            neighbors.extend(
                self.getUniqueNeighbors(lsNeighbors, neighbors, sequence)
            )

            # Generate all right shift mutants

            # Right shift sequence by one nucleotide (2 bits)
            temp = sequence >> self.bitCodeLen

            # Generate all possible right shift mutations by mutating the right
            # most nucleotide 4 times. However, do not include the mutant that
            # results in the source sequence.
            rsNeighbors = self.getShiftMutants("right", temp)

            # Append the list of right shift mutants to the neighborhood. Avoid
            # duplication, which is possible for sequences where the same letter
            # occupies all positions. Also, do not include the mutant that
            # results in the source sequence.
            neighbors.extend(
                self.getUniqueNeighbors(rsNeighbors, neighbors, sequence)
            )

        return neighbors

    # TODO: The code in this function is mostly shared with areNeighbors().
    #   Therefore, the two functions should be refactored such that the
    #   common code is factored out ...
    def get_mutation_type(self, bit_seq_1, bit_seq_2):
        # Check if seq1 and seq2 are 1-neighbors
        result = self._are_neighbors(bit_seq_1, bit_seq_2)

        if result == 'right' or result == 'left':
            return 'indel'
        elif result == 'point':
            # Full string representation of the bit string where only the bits
            # that differ between seq1 and seq2 are set. These would correspond
            # to the letter at which the two sequences differ.
            diff = "{0:{fill}{code_length}b}".format(
                bit_seq_1 ^ bit_seq_2,
                fill='0',
                code_length=self.seqLength * self.bitCodeLen
            )

            # Index of the first set bit. Only the first is enough to
            # determine which letter it is.
            index = diff.index('1')

            # Index of the letter at which the two sequences differ
            i = max(0, int(index / self.bitCodeLen))

            # Get the bit values for letter at index 'i' in both the source and
            # mutated sequences
            bit_letter_1 = self.getLetterAtIndex(bit_seq_1, i)
            bit_letter_2 = self.getLetterAtIndex(bit_seq_2, i)

            return f'{self.bitToLetterDict[bit_letter_1]}' \
                   f'{self.bitToLetterDict[bit_letter_2]}'
        else:
            if self._use_reverse_complements:
                rev_compl_seq2 = self.get_reverse_complement(bit_seq_2)

                result = self._are_neighbors(bit_seq_1, rev_compl_seq2)

                return result
            else:
                return None


class BitManipFactory:
    # Supported molecule types
    moleculeTypes = ["RNA", "DNA", "Protein", "Binary"]

    # Supported alphabet corresponding to each support molecule type
    alphabet = {
        'Binary': ['0', '1'],
        'RNA': ['A', 'U', 'C', 'G'],
        'DNA': ['A', 'T', 'C', 'G'],
        'Protein': [
            'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
        ]
    }

    @staticmethod
    def get_alphabet(molecule_type):
        return BitManipFactory.alphabet[molecule_type]

    # Returns the list of supported molecule types
    @staticmethod
    def getMoleculeTypes():
        return BitManipFactory.moleculeTypes

    # Based on the given molecule type, returns the corresponding
    # bit seq manipulator object.
    @staticmethod
    def getBitSeqManip(
            moleculeType,
            seqLength,
            useIndels,
            letter_to_neighbors,
            alphabet=None,
            useReverseComplements=False):

        return GenericBitSeqManipulator(
            sequence_length=seqLength,
            use_indels=useIndels,
            molecule_type=moleculeType,
            alphabet=alphabet,
            letter_to_neighbors=letter_to_neighbors,
            use_reverse_complements=useReverseComplements
        )
