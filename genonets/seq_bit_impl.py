"""
    seq_bit_impl
    ~~~~~~~~~~~~

    Contains classes that define alphabet types that can be used with
    the bit manipulation library.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

from seq_bit_lib import AbstractBitSeqManipulator
from genonets_utils import Utils


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
    def getBitSeqManip(moleculeType, seqLength, useIndels, useReverseComplements=False):
        if moleculeType == "RNA":
            return RNABitSeqManipulator(seqLength, useIndels)
        elif moleculeType == "DNA":
            return DNABitSeqManipulator(seqLength, useIndels,  useReverseComplements)
        elif moleculeType == "Protein":
            return ProteinBitSeqManipulator(seqLength, useIndels)
        elif moleculeType == "Binary":
            return BinaryBitSeqManipulator(seqLength, useIndels)
        else:
            print("Unsupported moleculeType: " + str(moleculeType))
            print("Exiting program ...")

            exit()
