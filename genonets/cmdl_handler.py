"""
    Handles command line arguments.

    :license: MIT, see LICENSE for more details.

"""

import argparse

from genonets.seq_bit_impl import BitManipFactory


class CmdParser:
    # pylint: disable=too-few-public-methods
    """
    Parses the command line arguments using python's 'argparse' module.

    """

    def __init__(self, arguments: list = None):
        """
        Does the parsing.

        :param arguments: Optional list of command-line arguments that can
                be used instead of reading parameters from the command-line.
                This is particularly useful for testing and embedding the
                package in an external program.

        """

        # Initialize the parser object
        parser = argparse.ArgumentParser()

        parser.add_argument(
            '-c', '--config-file',
            dest='config_file', default=None,
            help='When specified, all other arguments are treated as optional. '
                 'Any other arguments provided in addition to this will '
                 'override the values provided in the configuration file. It '
                 'is sufficient to provide all other arguments in the config '
                 'file, significantly simplifying the command-line.'
        )

        # Add 'inFilePath' as an argument
        parser.add_argument(
            '--input-file', default=None,
            dest='in_file_path',
            help='Complete path to the input file, relative to the '
                 'current directory.'
        )

        # Add 'outPath' as an argument
        parser.add_argument(
            '--output-path',
            dest='out_path', default=None,
            help='Path to the directory where output files should be '
                 'generated. The directory will be created if it does '
                 'not already exist.'
        )

        # Add 'alphabetType' as an argument
        parser.add_argument(
            '--alphabet',
            dest='alphabet_type', default=None,
            help='Each genotype in the input file must be a string of '
                 'letters from the alphabet selected here.',
            choices=BitManipFactory.getMoleculeTypes()
        )

        # Add 'tau' as an argument
        parser.add_argument(
            '--tau',
            dest='tau', default=None,
            help='Minimum score threshold: all genotypes in the input '
                 'file with "score" values below "Tau" will be '
                 'ignored. Tau must be a number.',
            type=float
        )

        # Add 'includeIndels' as an argument
        parser.add_argument(
            '--include-indels',
            dest='include_indels',
            action='store_const', const=True, default=None,
            help='If specified, mutations that shift the entire genotype '
                 'sequence by one letter should also be considered. Only '
                 'single point mutations are considered otherwise',
        )

        # Add 'use_reverse_complements' as an argument
        parser.add_argument(
            '-rc', '--use-reverse-complements',
            dest='use_reverse_complements',
            action='store_const', const=True,
            help='When specified, reverse complements are considered '
                 'during creation and analysis of genotype networks '
                 'for alphabet type DNA. This option is not valid for '
                 'other alphabet types.'
        )

        parser.add_argument(
            '--store-epistasis-squares',
            dest='save_squares', action='store_const', const=True,
            help='If this argument is specified and epistasis analysis is '
                 'used, all squares are stored on disk.'
        )

        parser.add_argument(
            '--use-all-components',
            dest='use_all_components',
            action='store_const', const=True, default=None,
            help='If specified, evolvability analysis considers all connected '
                 'components, not just the dominant network.',
        )

        # Add 'num_processes' as an argument
        parser.add_argument(
            '-np', '--num-processes',
            dest='num_procs', action='store',
            type=int,  default=None,
            help='No. of processes to be used in parallel processing'
        )

        # Add 'verbose' as an argument
        parser.add_argument(
            '-v', '--verbose',
            dest='verbose', action='store_const', const=True,
            help='Processing steps are printed to the screen during '
                 'program execution.'
        )

        parser.add_argument(
            '--genetic-code-file',
            dest='genetic_code_file', default=None,
            help='Path to the file which contains the mapping from the genetic '
                 'code to each letter in the selected alphabet.'
        )

        parser.add_argument(
            '--codon-alphabet',
            dest='codon_alphabet',
            choices=['DNA', 'RNA'], default=None,
            help='Alphabet to use for the codons in the genetic code file.'
        )

        parser.add_argument(
            '--include-indels-for-codons',
            dest='include_indels_for_codons', action='store_const',
            const=True, default=None,
            help='When specified, indels are considered when checking for '
                 'mutations in codons.'
        )

        parser.add_argument(
            '--use-rc-for-codons',
            dest='use_rc_for_codons', action='store_const',
            const=True, default=None,
            help='When specified, reverse complements are considered when '
                 'checking for mutations in codons. This is only applicable '
                 'when the alphabet for codons is DNA.'
        )

        # Keep an object level copy of the arguments dict
        if arguments:
            # Parse the string of arguments received
            self._args = parser.parse_args(arguments)
        else:
            # Parse sys.argv
            self._args = parser.parse_args()

    # Returns parsed arguments
    def get_args(self) -> argparse.Namespace:
        """
        Returns the parsed command-line arguments.

        :return: The populated Namespace object with all arguments.

        """

        return self._args
