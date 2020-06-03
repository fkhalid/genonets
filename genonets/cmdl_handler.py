"""
    Handles command line arguments.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import argparse

from genonets_writer import Writer
from seq_bit_impl import BitManipFactory
from genonets_constants import ErrorCodes
from genonets_exceptions import GenonetsError


# Parses the command line arguments using python's 'argparse' module
class CmdParser:
    # Constructor
    def __init__(self, arguments=None):
        # Initialize the parser object
        parser = argparse.ArgumentParser()

        # ---------------- Mandatory arguments ------------------ #

        # Add 'alphabetType' as an argument
        parser.add_argument("alphabetType",
                            help="Each genotype in the input file must be a string " +
                                 "of letters from the alphabet selected here.",
                            choices=BitManipFactory.getMoleculeTypes())

        # Add 'includeIndels' as an argument
        parser.add_argument("includeIndels",
                            help="Should be set to 'True' if mutations that shift " +
                                 "the entire genotype sequence by one letter should " +
                                 "also be considered. Only single point mutations are " +
                                 "considered when this parameter is set to 'False'.",
                            choices=["True", "true", "False", "false"])

        # Add 'inFilePath' as an argument
        parser.add_argument("inFilePath",
                            help="Complete path to the input file, relative to the " +
                                 "current directory.")

        # Add 'tau' as an argument
        parser.add_argument("tau",
                            help="Minimum score threshold: all genotypes in the input " +
                                 "file with 'score' values below 'Tau' will be ignored." +
                                 "Tau must be a number.",
                            type=float)

        # Add 'outPath' as an argument
        parser.add_argument("outPath",
                            help="Path to the directory where output files should be " +
                                 "generated. The directory will be created if it does " +
                                 "not already exist.")

        # ---------------- Optional arguments ------------------ #

        parser.add_argument('--genetic-code-file',
                            dest='genetic_code_file',
                            help='Path to the file which contains the mapping '
                                 'from the genetic code to each letter in the '
                                 'selected alphabet.')

        parser.add_argument('--codon-alphabet',
                            dest='codon_alphabet',
                            choices=['DNA', 'RNA'],
                            default='DNA',
                            help='Alphabet to use for the codons in '
                                 'the genetic code file.')

        # Add 'use_reverse_complements' as an argument
        parser.add_argument("-rc", "--use_reverse_complements",
                            dest="use_reverse_complements",
                            action="store_const", const=True,
                            help="When specified, reverse complements are considered during creation " +
                                 "and analysis of genotype networks for alphabet type DNA. This option is " +
                                 "not valid for other alphabet types.")

        # Add 'num_processes' as an argument
        parser.add_argument("-np", "--num_processes", dest="num_procs", action="store",
                            type=int, default="4",
                            help="No. of processes to be used in parallel processing")

        # Add 'verbose' as an argument
        parser.add_argument("-v", "--verbose", dest="verbose", action="store_const",
                            const=True,
                            help="Processing steps are printed to the screen during " +
                                 "program execution.")

        # Keep an object level copy of the arguments dict
        if arguments:
            # Parse the string of arguments received
            self.args = parser.parse_args(arguments)
        else:
            # Parse sys.argv
            self.args = parser.parse_args()

    # Returns parsed arguments
    def getArgs(self):
        return self.args


# Stores command line arguments
class CmdArgs:
    # Constructor. Accepts a list of arguments.
    def __init__(self, arguments):
        # Molecule type: RNA, DNA, Protein, etc.
        self.moleculeType = arguments.alphabetType

        # 'Use reverse complements' flag
        self.use_reverse_complements = True if arguments.use_reverse_complements else False

        # Report exception if 'use_reverse_complements' has been passed as an argument with
        # alphabet type other than DNA
        if self.use_reverse_complements and self.moleculeType != "DNA":
            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.RC_ALPHABET_MISMATCH))

            raise GenonetsError(ErrorCodes.RC_ALPHABET_MISMATCH)

        # Flag to indicate whether shift mutations should
        # be considered
        if arguments.includeIndels.lower() == "true":
            self.useIndels = True
        else:
            self.useIndels = False

        # Path to the input file
        self.inFilePath = arguments.inFilePath

        # Lower bound on fitness values to be used.
        self.tau = arguments.tau

        # Path to the output folder
        self.outPath = arguments.outPath

        # Make sure the path ends with "/", since this is needed
        # in the file writing routines
        if not self.outPath.endswith("/"):
            self.outPath += "/"

        # Optional file with the codon-to-letter mapping
        self.genetic_code_file = arguments.genetic_code_file

        # Alphabet to use for codons in the codon-to-letter mapping
        self.codon_alphabet = arguments.codon_alphabet

        # Maximum number of parallel processes to be used
        self.num_procs = arguments.num_procs

        # Verbose flag
        self.verbose = True if arguments.verbose else False

        # Create a dictionary of parameters
        paramsDict = {
            "alphabetType": self.moleculeType,
            "includeIndels": str(self.useIndels),
            "inFilePath": self.inFilePath,
            "tau": str(self.tau),
            "outPath": self.outPath,
            "genetic_code_file": self.genetic_code_file,
            "codon_alphabet": self.codon_alphabet,
            "useReverseComplements": str(self.use_reverse_complements),
            "num_procs": str(self.num_procs),
            "verbose": str(self.verbose)
        }

        # Print the parsed parameter values
        self.printInParams(paramsDict)

        # Write input parameters to file
        Writer.writeInParamsToFile(paramsDict, self.outPath)

    # Print parsed input parameters
    def printInParams(self, paramsDict):
        print("\nParsed input parameter values:")
        print("------------------------------")

        # For each parameter,
        for param in paramsDict:
            # Print the 'parameter : value' pair
            print(param + ": " + paramsDict[param])

        print("------------------------------\n")
