"""
    genonets_reader
    ~~~~~~~~~~~~~~~

    Parser for the genonets input file format.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import csv

from genonets_exceptions import GenonetsError
from genonets_constants import ErrorCodes, SupportedAlphabet


class InReader:
    # Required input file column headers
    COL_HEADERS = ["Genotypeset", "Delta", "Genotype", "Score"]

    # Uses the given DictReader to construct a dictionary from file with,
    # key=Genotypeset, value=dict{key=sequence, value=score}
    # Returns the dictionary
    # TODO: This function can be extend to also return a reverse dictionary
    # with key=score, value=sequence ...
    @staticmethod
    def buildDataDicts(inFilePath, tau, alphabetType):
        # Dicts to be returned
        dataDict = {}
        deltaDict = {}
        sequences = []

        # Genotype length to be determined
        seqLength = 0

        # Get handles to the input file and a DictReader for the file
        reader, inFile = InReader.getDictReaderForFile(inFilePath)

        # Check if all the required column headers are available in the file
        if not InReader.reqHdrsArePresent(reader.fieldnames):
            inFile.close()

            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.INCONSISTENT_HEADER))

            raise GenonetsError(ErrorCodes.INCONSISTENT_HEADER)

        # For each data row in the file
        for row in reader:
            # Check for missing values in this row
            if any(row[col] in (None, "") for col in row.keys()):
                inFile.close()

                lineNum = str(int(reader.line_num))

                print("Error: " +
                      ErrorCodes.getErrDescription(ErrorCodes.MISSING_VALUE) +
                      ": Line No. " + lineNum)

                raise GenonetsError(ErrorCodes.MISSING_VALUE,
                                    "Line No. " + lineNum)

            # Get the fitness score
            try:
                score = float(row["Score"])
            except:
                inFile.close()

                lineNum = str(int(reader.line_num))

                print("Error: " +
                      ErrorCodes.getErrDescription(ErrorCodes.BAD_SCORE_FORMAT) +
                      ": Line No. " + lineNum)

                raise GenonetsError(ErrorCodes.BAD_SCORE_FORMAT,
                                    "Line No. " + lineNum)

            # If the score for this sequence is greater than or equal to
            # the given threshold,
            if score >= tau:
                # If the current genotypeset has not already been added
                if row["Genotypeset"] not in dataDict:
                    # Initialize dict for the genotypeset
                    dataDict[row["Genotypeset"]] = {}

                    # Get the delta value
                    try:
                        delta = float(row["Delta"])
                    except:
                        inFile.close()

                        lineNum = str(int(reader.line_num))

                        print("Error: " +
                              ErrorCodes.getErrDescription(ErrorCodes.BAD_DELTA_FORMAT) +
                              ": Line No. " + lineNum)

                        raise GenonetsError(
                            ErrorCodes.BAD_DELTA_FORMAT,
                            "Line No. " + lineNum)

                    deltaDict[row["Genotypeset"]] = delta

                # Get the genotype sequence
                sequence = row["Genotype"]

                # If sequence length has not been initialized yet, i.e., this is
                # the first row,
                if seqLength == 0:
                    # Set length of the current genotype as the sequence length
                    # for the entire dataset
                    seqLength = len(sequence)

                try:
                    InReader.verifyGenotype(sequence, seqLength, alphabetType,
                                            str(int(reader.line_num)))
                except Exception as e:
                    inFile.close()
                    raise e

                # Add sequence as key and score as value to the current
                # genotypeset
                dataDict[row["Genotypeset"]][sequence] = score

                # If the sequence has not already been read for any other
                # genotypeset,
                if sequence not in sequences:
                    # Add it to the list of unique sequences found in the
                    # input file
                    sequences.append(sequence)

        inFile.close()

        # If no genotypes were found with score >= tau,
        if not dataDict:
            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.NO_USABLE_SCORES)
                  + ": Tau=" + str(tau))

            raise GenonetsError(
                ErrorCodes.NO_USABLE_SCORES,
                "Tau=" + str(tau))

        # Dictionary: Key=Sequence, Value=[Genotypesets]. Reverse dictionary
        # that is used in functions like evolvability.
        seqToRepDict = InReader.getSeqToRepDict(sequences, dataDict)

        return dataDict, deltaDict, seqToRepDict, seqLength

    # Builds Dictionary: Key=Sequence, Value=[Genotypesets in which this sequence
    # exists]
    @staticmethod
    def getSeqToRepDict(sequences, dataDict):
        seqToRepDict = {}

        # For each sequence in the given list of sequences,
        for seq in sequences:
            # Initialize the list genotypesets
            seqToRepDict[seq] = []

            # For each genotypeset in the given data dictionary,
            for rep in dataDict.keys():
                # If the sequence exists in the genotypeset,
                if seq in dataDict[rep]:
                    # Append the genotypeset name to the list of values
                    # for this sequence
                    seqToRepDict[seq].append(rep)

        return seqToRepDict

    # If any of the required column headers are not in the file, returns
    # False. Otherwise return True.
    @staticmethod
    def reqHdrsArePresent(colNames):
        colsInFile = [col for col in colNames]

        if any(header not in colsInFile for header in InReader.COL_HEADERS):
            return False

        return True

    @staticmethod
    def verifyGenotype(genotype, seqLen, alphabetType, lineNum):
        # Verify the length
        if len(genotype) != seqLen:
            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.INCONSISTENT_SEQ_LEN) +
                  ": Line No. " + lineNum)

            raise GenonetsError(
                ErrorCodes.INCONSISTENT_SEQ_LEN,
                "Line No. " + lineNum)

        # Get the alphabet corresponding to the type received as
        # argument
        alphabet = SupportedAlphabet.getAlphabet(alphabetType)

        # Verify alphabet
        if any(letter not in alphabet for letter in genotype):
            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.ALPHABET_TYPE_MISMATCH) +
                  ": Line No. " + lineNum)

            raise GenonetsError(
                ErrorCodes.ALPHABET_TYPE_MISMATCH,
                "Line No. " + lineNum)

    # Open the given CSV file, and return a handle to the dict reader
    # for the file, as well as the file handle.
    @staticmethod
    def getDictReaderForFile(fileName):
        # Open file
        try:
            dataFile = open(fileName, 'rU')
        except Exception as e:
            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.UNKNOWN_PARSING_ERROR))

            raise GenonetsError(ErrorCodes.UNKNOWN_PARSING_ERROR)

        # Read the file into a dictionary
        reader = csv.DictReader(dataFile, delimiter="\t")

        return reader, dataFile
