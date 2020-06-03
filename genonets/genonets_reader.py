"""
    genonets_reader
    ~~~~~~~~~~~~~~~

    Parser for the genonets input file format.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

import csv

from genonets_exceptions import GenonetsError
from genonets_constants import ErrorCodes
from genonets_constants import SupportedAlphabet


class GeneticCodeReader:
    def __init__(self):
        pass

    COL_HEADERS = {'Codon', 'Letter'}

    @staticmethod
    def load_codon_to_letter_map(filename, alphabet_type):
        codon_length = 0
        codon_to_letter = {}

        # Get the alphabet corresponding to the type received as
        # argument
        alphabet = set(SupportedAlphabet.getAlphabet(alphabet_type))

        # Load the file
        try:
            mapping_file = open(filename, 'rU')
        except Exception as e:
            print('Error: Could not open the codon-to-letter '
                  'mapping file: ' + filename)

            raise GenonetsError(ErrorCodes.UNKNOWN_PARSING_ERROR)

        # Read the file into a dictionary
        reader = csv.DictReader(mapping_file, delimiter='\t')

        if set(reader.fieldnames) ^ GeneticCodeReader.COL_HEADERS:
            mapping_file.close()

            print('Error: Unsupported headers in the codon-to-letter '
                  'mapping file: ' + filename)

            raise GenonetsError(ErrorCodes.INCONSISTENT_HEADER)

        # For each row,
        for row in reader:
            # Check for missing values in this row
            if any(row[col] in (None, "") for col in row.keys()):
                mapping_file.close()

                line_number = str(int(reader.line_num))

                print('Error: Missing value on line No. ' + line_number)

                raise GenonetsError(ErrorCodes.MISSING_VALUE,
                                    "Line No. " + line_number)

            # Get the codon
            codon = row['Codon'].strip(' ')

            if codon_length == 0:
                codon_length = len(codon)
            elif len(codon) != codon_length:
                print('Error: Inconsisten codon length on line No. '
                      + str(int(reader.line_num)))

                raise GenonetsError(
                    ErrorCodes.INCONSISTENT_SEQ_LEN,
                    "Line No. " + str(int(reader.line_num)))

            # Verify alphabet
            if set(codon) - alphabet:
                print('Error: Alphabet type mismatch on line No. '
                      + str(int(reader.line_num)))

                raise GenonetsError(
                    ErrorCodes.ALPHABET_TYPE_MISMATCH,
                    "Line No. " + str(int(reader.line_num)))

            # Load the letter
            letter = row['Letter'].strip(' ')

            if len(letter) != 1:
                print('Error: Length of a letter cannot be more than '
                      'one. Line No. ' + str(int(reader.line_num)))

                raise GenonetsError('Unsupported letter length.')

            # Save the mapping
            codon_to_letter[codon] = letter

        mapping_file.close()

        return codon_to_letter, codon_length


class InReader:
    def __init__(self):
        pass

    # Required input file column headers
    COL_HEADERS = ["Genotypeset", "Delta", "Genotype", "Score"]

    # Uses the given DictReader to construct a dictionary from file with,
    # key=Genotypeset, value=dict{key=genotype, value=score}
    # Returns the dictionary
    # TODO: This function can be extend to also return a reverse dictionary
    # with key=score, value=sequence ...
    @staticmethod
    def build_data_dicts(in_file_path, tau, alphabet_type):
        # Data structures to be returned
        data_dict = {}
        delta_dict = {}
        genotypes = []      # List of unique genotypes across all genotype sets
        genotype_sets = []  # List of genotype sets in the order in which they are read from file

        # Genotype length to be determined
        genotype_length = 0

        # Get handle to the input file and a DictReader for the file
        reader, in_file = InReader.dict_reader_for_file(in_file_path)

        # Check if all the required column headers are available in the file
        if not InReader.req_hdrs_are_present(reader.fieldnames):
            in_file.close()

            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.INCONSISTENT_HEADER))

            raise GenonetsError(ErrorCodes.INCONSISTENT_HEADER)

        # For each data row in the file,
        for row in reader:
            # Check for missing values in this row
            if any(row[col] in (None, "") for col in row.keys()):
                in_file.close()

                line_number = str(int(reader.line_num))

                print("Error: " +
                      ErrorCodes.getErrDescription(ErrorCodes.MISSING_VALUE) +
                      ": Line No. " + line_number)

                raise GenonetsError(ErrorCodes.MISSING_VALUE,
                                    "Line No. " + line_number)

            # Get the fitness score
            try:
                score = float(row["Score"])
            except:
                in_file.close()

                line_number = str(int(reader.line_num))

                print("Error: " +
                      ErrorCodes.getErrDescription(ErrorCodes.BAD_SCORE_FORMAT) +
                      ": Line No. " + line_number)

                raise GenonetsError(ErrorCodes.BAD_SCORE_FORMAT,
                                    "Line No. " + line_number)

            # If the score for this genotype is greater than or equal to
            # the given threshold,
            if score >= tau:
                # If the current genotype set has not already been added,
                if row["Genotypeset"] not in data_dict:
                    # Initialize dict for the genotype set
                    data_dict[row["Genotypeset"]] = {}

                    # Add the name of this genotype set to the ordered list of
                    # genotype set names
                    genotype_sets.append(row["Genotypeset"])

                    # Get the delta value
                    try:
                        delta = float(row["Delta"])
                    except:
                        in_file.close()

                        line_number = str(int(reader.line_num))

                        print("Error: " +
                              ErrorCodes.getErrDescription(ErrorCodes.BAD_DELTA_FORMAT) +
                              ": Line No. " + line_number)

                        raise GenonetsError(
                            ErrorCodes.BAD_DELTA_FORMAT,
                            "Line No. " + line_number)

                    delta_dict[row["Genotypeset"]] = delta

                # Get the genotype sequence
                genotype = row["Genotype"]

                # If genotype length has not been initialized yet, i.e., this is
                # the first row,
                if genotype_length == 0:
                    # Set length of the current genotype as the genotype length
                    # for the entire dataset
                    genotype_length = len(genotype)

                try:
                    InReader.verify_genotype(genotype, genotype_length, alphabet_type,
                                             str(int(reader.line_num)))
                except Exception as e:
                    in_file.close()
                    raise e

                # Add genotype as key and score as value to the current
                # genotype set
                data_dict[row["Genotypeset"]][genotype] = score

                # If the genotype has not already been read for any other
                # genotype set,
                if genotype not in genotypes:
                    # Add it to the list of unique sequences found in the
                    # input file
                    genotypes.append(genotype)

        in_file.close()

        # If no genotypes were found with score >= tau,
        if not data_dict:
            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.NO_USABLE_SCORES) +
                  ": Tau=" + str(tau))

            raise GenonetsError(
                ErrorCodes.NO_USABLE_SCORES,
                "Tau=" + str(tau))

        # Dictionary: Key=Sequence, Value=[Genotype sets]. Reverse dictionary
        # that is used in functions like evolvability.
        genotype_to_set_dict = InReader.build_genotype_to_set_dict(genotypes, data_dict)

        return data_dict, delta_dict, genotype_to_set_dict, genotype_length, genotype_sets

    # Builds Dictionary: Key=Sequence, Value=[Genotype sets in which this sequence
    # exists]
    @staticmethod
    def build_genotype_to_set_dict(genotypes, data_dict):
        genotype_to_set_dict = {}

        # For each sequence in the given list of genotypes,
        for seq in genotypes:
            # Initialize the list genotype sets
            genotype_to_set_dict[seq] = []

            # For each genotype set in the given data dictionary,
            for rep in data_dict.keys():
                # If the sequence exists in the genotype set,
                if seq in data_dict[rep]:
                    # Append the genotype set name to the list of values
                    # for this sequence
                    genotype_to_set_dict[seq].append(rep)

        return genotype_to_set_dict

    # If any of the required column headers are not in the file, returns
    # False. Otherwise return True.
    @staticmethod
    def req_hdrs_are_present(col_names):
        cols_in_file = [col for col in col_names]

        if any(header not in cols_in_file for header in InReader.COL_HEADERS):
            return False

        return True

    @staticmethod
    def verify_genotype(genotype, genotype_length, alphabet_type, line_number):
        # Verify the length
        if len(genotype) != genotype_length:
            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.INCONSISTENT_SEQ_LEN) +
                  ": Line No. " + line_number)

            raise GenonetsError(
                ErrorCodes.INCONSISTENT_SEQ_LEN,
                "Line No. " + line_number)

        # Get the alphabet corresponding to the type received as
        # argument
        alphabet = SupportedAlphabet.getAlphabet(alphabet_type)

        # Verify alphabet
        if any(letter not in alphabet for letter in genotype):
            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.ALPHABET_TYPE_MISMATCH) +
                  ": Line No. " + line_number)

            raise GenonetsError(
                ErrorCodes.ALPHABET_TYPE_MISMATCH,
                "Line No. " + line_number)

    # Open the given CSV file, and return a handle to the dict reader
    # for the file, as well as the file handle.
    @staticmethod
    def dict_reader_for_file(file_name):
        # Open file
        try:
            data_file = open(file_name, 'rU')
        except Exception as e:
            print("Error: " +
                  ErrorCodes.getErrDescription(ErrorCodes.UNKNOWN_PARSING_ERROR))

            raise GenonetsError(ErrorCodes.UNKNOWN_PARSING_ERROR)

        # Read the file into a dictionary
        reader = csv.DictReader(data_file, delimiter="\t")

        return reader, data_file
