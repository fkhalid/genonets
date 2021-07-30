"""
    This module is responsible for parsing and storing application
    configuration.

    :license: MIT, see LICENSE for more details.

"""

import argparse
import configparser

from genonets.writer import Writer
from genonets.errors import FileError
from genonets.errors import MissingMandatoryArgumentError


# TODO: Improvement is required for handling of ignored options, i.e., if
#   one or more options are not supported, these should perhaps be explicitly
#   set to 'None' or the appropriate value used to signal these are not being
#   used. E.g., the [Genetic code] section used for a molecule type other than
#   Protein ...
class ConfigParser:
    """
    Parses and stores the application configuration file.

    """

    def __init__(self, args: argparse.Namespace):
        """
        Initializes with the given arguments.

        :param args: Object of type argparse.Namespace.

        """

        print('')

        self._tau = None
        self._verbose = None
        self._out_path = None
        self._num_procs = None
        self._use_indels = None
        self._in_file_path = None
        self._save_squares = None
        self._molecule_type = None
        self._codon_alphabet = None
        self._genetic_code_file = None
        self._use_rc_for_codons = None
        self._use_all_components = None
        self._use_reverse_complements = None
        self._include_indels_for_codons = None

        self._config = None
        self._cmd_args = args
        self._config_file = args.config_file

        if self._config_file is not None:
            self._config = configparser.ConfigParser()
            files_read = self._config.read(self._config_file)

            if not files_read:
                raise FileError(
                    f'Failed to read the provided configuration file: '
                    f'{self._config_file}'
                )
        else:
            print(
                'Warning: No configuration file specified. All '
                'parameters will be read from the command-line.'
            )

        self._set_all()

        # Create a dictionary of parameters
        params_dict = self._create_params_dict()

        # Print the parsed parameter values
        self.print_params(params_dict)

        # Write input parameters to file
        Writer.writeInParamsToFile(params_dict, self._out_path)

    def _set_all(self) -> None:
        self._set_general()
        self._set_genetic_code()

    def _set_general(self) -> None:
        section = 'General'

        # Verbose flag
        if self._cmd_args.verbose is not None:
            self._verbose = self._cmd_args.verbose
        else:
            if self._config is not None:
                self._verbose = self._config.getboolean(
                    section=section,
                    option='verbose',
                    fallback=False
                )
            else:
                self._verbose = False

        # --store-epistasis-squares
        if self._cmd_args.save_squares is not None:
            self._save_squares = self._cmd_args.save_squares
        else:
            if self._config is not None:
                self._save_squares = self._config.getboolean(
                    section=section,
                    option='store-epistasis-squares',
                    fallback=False
                )
            else:
                self._save_squares = False

        # Path to the output folder
        if self._cmd_args.out_path is not None:
            self._out_path = self._cmd_args.out_path
        else:
            if self._config is not None:
                self._out_path = self._config.get(
                    section=section,
                    option='output-path',
                    fallback=None
                )
            else:
                raise MissingMandatoryArgumentError(
                    f'Path to the output directory is a mandatory '
                    f'argument that must be specified.'
                )

        # Make sure the path ends with "/", since this is needed
        # in the file writing routines
        if not self._out_path.endswith("/"):
            self._out_path += "/"

        # Path to the input file
        if self._cmd_args.in_file_path is not None:
            self._in_file_path = self._cmd_args.in_file_path
        else:
            if self._config is not None:
                self._in_file_path = self._config.get(
                    section=section,
                    option='input-file',
                    fallback=None
                )
            else:
                raise MissingMandatoryArgumentError(
                    f'Path to the input directory is a mandatory '
                    f'argument that must be specified.'
                )

        # Tau
        if self._cmd_args.tau is not None:
            self._tau = self._cmd_args.tau
        else:
            if self._config is not None:
                self._tau = self._config.getfloat(
                    section=section,
                    option='tau',
                    fallback=None
                )
            else:
                raise MissingMandatoryArgumentError(
                    f'Tau is a mandatory argument that must be '
                    f'specified.'
                )

        # # Molecule type: RNA, DNA, Protein, etc.
        if self._cmd_args.alphabet_type is not None:
            self._molecule_type = self._cmd_args.alphabet_type
        else:
            if self._config is not None:
                self._molecule_type = self._config.get(
                    section=section,
                    option='alphabet',
                    fallback=None
                )
            else:
                raise MissingMandatoryArgumentError(
                    f'Alphabet type is a mandatory argument that must '
                    f'be specified.'
                )

        # Flag to indicate whether shift mutations should be considered
        if self._cmd_args.include_indels is not None:
            self._use_indels = self._cmd_args.include_indels
        else:
            if self._config is not None:
                self._use_indels = self._config.getboolean(
                    section=section,
                    option='include-indels',
                    fallback=False
                )
            else:
                self._use_indels = False

        # 'Use reverse complements' flag
        if self._cmd_args.use_reverse_complements is not None:
            self._use_reverse_complements = \
                self._cmd_args.use_reverse_complements
        else:
            if self._config is not None:
                self._use_reverse_complements = self._config.getboolean(
                    section=section,
                    option='use-reverse-complements',
                    fallback=False
                )
            else:
                self._use_reverse_complements = False

        # Issue a warning if 'use_reverse_complements' has been passed
        # as an argument with alphabet type other than DNA
        if self._use_reverse_complements and self._molecule_type != 'DNA':
            self._use_reverse_complements = False
            print(
                'Warning: Reverse complements will be ignored, because '
                'the alphabet type is not DNA.'
            )

        # No. of processes to use
        if self._cmd_args.num_procs is not None:
            self._num_procs = self._cmd_args.num_procs
        else:
            if self._config is not None:
                self._num_procs = self._config.getint(
                    section=section,
                    option='num-processes',
                    fallback=1
                )
            else:
                self._num_procs = 1

        if self._cmd_args.use_all_components is not None:
            self._use_all_components = self._cmd_args.use_all_components
        else:
            if self._config is not None:
                self._use_all_components = self._config.getboolean(
                    section=section,
                    option='use-all-components',
                    fallback=False
                )
            else:
                self._use_all_components = False

    def _set_genetic_code(self) -> None:
        section = 'Genetic Code'

        # Genetic code file
        if self._cmd_args.genetic_code_file is not None:
            self._genetic_code_file = self._cmd_args.genetic_code_file
        else:
            if self._config is not None:
                self._genetic_code_file = self._config.get(
                    section=section,
                    option='genetic-code-file',
                    fallback=None
                )
            else:
                self._genetic_code_file = None

        if self._genetic_code_file is not None and self._molecule_type != 'Protein':
            self._genetic_code_file = None
            print(
                f'Ignoring all [Genetic Code] options : These options are only '
                f'supported for molecule type Protein; these are not supported '
                f'for molecule type {self._molecule_type}.'
            )

        # Codon alphabet
        if self._cmd_args.codon_alphabet is not None:
            self._codon_alphabet = self._cmd_args.codon_alphabet
        else:
            if self._config is not None:
                self._codon_alphabet = self._config.get(
                    section=section,
                    option='codon-alphabet',
                    fallback=None
                )
            else:
                self._codon_alphabet = None

        # Include indels for codons
        if self._cmd_args.include_indels_for_codons is not None:
            self._include_indels_for_codons = \
                self._cmd_args.include_indels_for_codons
        else:
            if self._config is not None:
                self._include_indels_for_codons = self._config.getboolean(
                    section=section,
                    option='include-indels-for-codons',
                    fallback=False
                )
            else:
                self._include_indels_for_codons = False

        # Use reverse complements for codons
        if self._cmd_args.use_rc_for_codons is not None:
            self._use_rc_for_codons = \
                self._cmd_args.use_rc_for_codons
        else:
            if self._config is not None:
                self._use_rc_for_codons = self._config.getboolean(
                    section=section,
                    option='use-rc-for-codons',
                    fallback=False
                )
            else:
                self._use_rc_for_codons = False

        if self._codon_alphabet != 'DNA' and self._use_rc_for_codons:
            self._use_rc_for_codons = False
            print(
                f'Ignoring option --codon-use-rc: This option can only '
                f'be used when the codon alphabet is DNA; it is not '
                f'supported for codon alphabet {self._codon_alphabet}.'
            )

    def _create_params_dict(self) -> dict:
        """
        Create and returns a mapping such that each key is a description of an
        argument, and the value of the corresponding value.

        :return: The mapping {Descriptive name: value}

        """

        return {
            'Genonets input file': self._in_file_path,
            'Output path': self._out_path,
            'Genotype alphabet': self._molecule_type,
            'Tau': self.tau,
            'Include indels for Genotypes': self._use_indels,
            'Use reverse complements for genotypes':
                self._use_reverse_complements,
            'Store epistasis squares': self._save_squares,
            'Use all components for evolvability': self._use_all_components,
            'No. of parallel processes': self._num_procs,
            'Verbose': self._verbose,
            'Genetic code input file': self._genetic_code_file,
            'Codon alphabet': self._codon_alphabet,
            'Include indels for Codons': self._include_indels_for_codons,
            'Use reverse complements for codons': self._use_rc_for_codons,
        }

    @staticmethod
    def print_params(params: dict) -> None:
        """
        Write the key-value pairs from the received dict to stdout.

        :param params: Parameter names and values as key-value pairs.

        """

        print("\nParsed input parameter values:")
        print("------------------------------")

        # For each parameter,
        for param in params:
            # Print the 'parameter : value' pair
            print(f'{param}: {params[param]}')

        print("------------------------------\n")

    @property
    def tau(self) -> float:
        """
        Returns the score threshold.

        :return: The configured threshold value.

        """

        return self._tau

    @property
    def verbose(self) -> bool:
        """
        Returns the flag to indicate whether or not output should be verbose.

        :return: The configured flag value.

        """

        return self._verbose

    @property
    def save_squares(self) -> bool:
        """
        Returns the flag to indicate whether or not epistasis squares should
        be written to disk.

        :return: The configured flag value.

        """

        return self._save_squares

    @property
    def output_path(self) -> str:
        """
        Returns the path to the directory in which to generate the output files.

        :return: The configured path to the output directory.

        """

        return self._out_path

    @property
    def num_processes(self) -> int:
        """
        Returns the No. of parallel processes to use.

        :return: No. of configured parallel processes.

        """

        return self._num_procs

    @property
    def use_indels(self) -> bool:
        """
        Returns the flag that indicates whether or not indels should be
        considered during the construction and analysis of genotype networks.

        :return: The configured flag value.

        """

        return self._use_indels

    @property
    def in_file_path(self) -> str:
        """
        Returns the path to the directory that contains the input file.

        :return: The configured path.

        """

        return self._in_file_path

    @property
    def molecule_type(self) -> str:
        """
        Returns the molecule type that corresponding the genotype alphabet.

        :return: The configure molecule type.

        """

        return self._molecule_type

    @property
    def codon_alphabet(self) -> str:
        """
        Return the alphabet to use for the codons in the genetic code file.

        :return: The configure codon alphabet.

        """

        return self._codon_alphabet

    @property
    def genetic_code_file(self) -> str:
        """
        Returns the path to the genetic code file, including the file name.

        :return: The configured file name.

        """

        return self._genetic_code_file

    @property
    def use_rc_for_codons(self) -> bool:
        """
        Returns the flag the indicates whether or not reverse complements
        should be considered for codons.

        :return: The configured flag.

        """

        return self._use_rc_for_codons

    @property
    def use_reverse_complements(self) -> bool:
        """
        Returns the flag that indicates whether reverse complements should be
        used.

        :return: The configured flag.

        """

        return self._use_reverse_complements

    @property
    def include_indels_for_codons(self) -> bool:
        """
        Returns the flag that indicates whether or not indels should be
        considered for codons.

        :return: The configured flag.

        """

        return self._include_indels_for_codons

    @property
    def use_all_components(self):
        """
        Returns the flag that indicates whether or not all connected components
        should be considered for analysis.

        :return: The configured flag.

        """

        return self._use_all_components
