"""
    cmdl_handler
    ~~~~~~~~~~~~

    Handles command line arguments.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import argparse

from genonets_writer import Writer
from seq_bit_impl import BitManipFactory


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

        # Add 'num_processes' as an argument
        parser.add_argument("-np", "--num_processes", dest="num_procs", action="store",
                            type=int, default="4",
                            help="No. of processes to be used in parallel processing")

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

        # Maximum number of parallel processes to be used
        self.num_procs = arguments.num_procs

        # Create a dictionary of parameters
        paramsDict = {
            "alphabetType": self.moleculeType,
            "includeIndels": str(self.useIndels),
            "inFilePath": self.inFilePath,
            "tau": str(self.tau),
            "outPath": self.outPath,
            "num_procs": str(self.num_procs)
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
        for param in paramsDict.keys():
            # Print the 'parameter : value' pair
            print(param + ": " + paramsDict[param])

        print("------------------------------\n")
