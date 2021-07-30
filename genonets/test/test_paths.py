
import tempfile

import genonets.test.utils as utils
import genonets.test.compare_result_files as comparator

from genonets.cmdl_handler import CmdParser
from genonets.interface import Genonets
from genonets.constants import AnalysisConstants as Ac


class TestPaths:
    @staticmethod
    def run_test(cmd_args, ground_truth_dir, data_dir):
        args = CmdParser(arguments=cmd_args).get_args()

        gn = Genonets(args)
        gn.create()
        gn.analyze(analyses=[Ac.PATHS])
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
    def test_1():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_1'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test1_input.tsv',
                '--codon-alphabet=RNA',
                '--genetic-code-file=genonets/test/data/inputs/paths/code_standard.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_2():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_2'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test2_input.tsv',
                '--codon-alphabet=RNA',
                '--genetic-code-file=genonets/test/data/inputs/paths/code_standard.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_3():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_3'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test3_input.tsv',
                '--codon-alphabet=RNA',
                '--genetic-code-file=genonets/test/data/inputs/paths/code_standard.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_4():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_4'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test4_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_5():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_5'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test5_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_6():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_6'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test6_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_7():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_7'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test7_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_8():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_8'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test8_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_9():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_9'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test9_input.tsv',
                '--codon-alphabet=RNA',
                '--genetic-code-file=genonets/test/data/inputs/paths/code_standard.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_10():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_10'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test10_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_11():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_11'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test11_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_12():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_12'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test12_input.tsv',
                '--codon-alphabet=RNA',
                '--genetic-code-file=genonets/test/data/inputs/paths/code_standard.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_13():
        ground_truth_dir = 'genonets/test/data/ground_truth/paths/test_13'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/paths/test13_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestPaths.run_test(cmd_args, ground_truth_dir, data_dir)
