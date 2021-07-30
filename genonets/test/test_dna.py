
import tempfile

import genonets.test.utils as utils
import genonets.test.compare_result_files as comparator

from genonets.cmdl_handler import CmdParser
from genonets.interface import Genonets


class TestDNA:
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

        assert comparator.compare_overlap_results(
            ground_truth_dir, data_dir
        )

    @staticmethod
    def test_no_indels_no_rc():
        ground_truth_dir = 'genonets/test/data/ground_truth/dna/mus/no_indels_no_rc'

        with tempfile.TemporaryDirectory(prefix='test_dna_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.35',
                '--input-file=genonets/test/data/inputs/dna/input_sample_dna-mus.tsv',
                f'--output-path={data_dir}'
            ]

            TestDNA.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_no_indels_with_rc():
        ground_truth_dir = 'genonets/test/data/ground_truth/dna/mus/no_indels_with_rc'

        with tempfile.TemporaryDirectory(prefix='test_dna_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.35',
                '--input-file=genonets/test/data/inputs/dna/input_sample_dna-mus.tsv',
                '--use-reverse-complements',
                f'--output-path={data_dir}'
            ]

            TestDNA.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_with_indels_no_rc():
        ground_truth_dir = 'genonets/test/data/ground_truth/dna/mus/with_indels_no_rc'

        with tempfile.TemporaryDirectory(prefix='test_dna_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.35',
                '--input-file=genonets/test/data/inputs/dna/input_sample_dna-mus.tsv',
                '--include-indels',
                f'--output-path={data_dir}'
            ]

            TestDNA.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_with_indels_with_rc():
        ground_truth_dir = 'genonets/test/data/ground_truth/dna/mus/with_indels_with_rc'

        with tempfile.TemporaryDirectory(prefix='test_dna_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.35',
                '--input-file=genonets/test/data/inputs/dna/input_sample_dna-mus.tsv',
                '--include-indels',
                '--use-reverse-complements',
                f'--output-path={data_dir}'
            ]

            TestDNA.run_test(cmd_args, ground_truth_dir, data_dir)
