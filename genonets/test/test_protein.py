
import tempfile

import genonets.test.utils as utils
import genonets.test.compare_result_files as comparator

from genonets.cmdl_handler import CmdParser
from genonets.interface import Genonets


class TestProtein:
    @staticmethod
    def run_test(cmd_args, ground_truth_dir, data_dir):
        args = CmdParser(arguments=cmd_args).get_args()

        gn = Genonets(args)
        gn.create()
        gn.analyze()
        gn.save_network_results()
        gn.save_genotype_results()

        assert utils.num_files_matches(ground_truth_dir, data_dir)

        assert comparator.compare_genotype_set_measures(
            ground_truth_dir, data_dir
        )

        assert comparator.compare_genotype_measures(
            ground_truth_dir, data_dir
        )

    @staticmethod
    def test_no_indels():
        ground_truth_dir = 'genonets/test/data/ground_truth/protein/h3nb/no_indels'

        with tempfile.TemporaryDirectory(prefix='test_protein_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/protein/input_h3nb.tsv',
                f'--output-path={data_dir}'
            ]

            TestProtein.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_with_indels():
        ground_truth_dir = 'genonets/test/data/ground_truth/protein/h3nb/with_indels'

        with tempfile.TemporaryDirectory(prefix='test_protein_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--include-indels',
                '--input-file=genonets/test/data/inputs/protein/input_h3nb.tsv',
                f'--output-path={data_dir}'
            ]

            TestProtein.run_test(cmd_args, ground_truth_dir, data_dir)
