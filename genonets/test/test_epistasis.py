
import tempfile

import genonets.test.utils as utils
import genonets.test.compare_result_files as comparator

from genonets.cmdl_handler import CmdParser
from genonets.interface import Genonets
from genonets.constants import AnalysisConstants as Ac


class TestEpistasis:
    @staticmethod
    def run_test(cmd_args, ground_truth_dir, data_dir):
        args = CmdParser(arguments=cmd_args).get_args()

        gn = Genonets(args)
        gn.create()
        gn.analyze(analyses=[Ac.EPISTASIS])
        gn.save_network_results()
        gn.save_genotype_results()

        assert utils.num_files_matches(ground_truth_dir, data_dir)

        assert comparator.compare_genotype_set_measures(
            ground_truth_dir, data_dir
        )

    @staticmethod
    def test_1():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_1'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test1_input.tsv',
                '--codon-alphabet=RNA',
                '--genetic-code-file=genonets/test/data/inputs/epistasis/code_standard.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_2():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_2'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test2_input.tsv',
                '--codon-alphabet=RNA',
                '--genetic-code-file=genonets/test/data/inputs/epistasis/code_standard.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_3():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_3'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test3_input.tsv',
                '--codon-alphabet=RNA',
                '--genetic-code-file=genonets/test/data/inputs/epistasis/code_standard.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_4():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_4'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test4_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_5():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_5'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test5_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_6():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_6'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test6_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_7():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_7'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test7_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_8():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_8'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test8_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_9():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_9'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test9_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_10():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_10'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test10_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_11():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_11'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test11_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_12():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_12'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test12_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_13():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_13'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test13_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_14():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_14'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test14_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_15():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_15'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test15_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_16():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_16'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test16_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_17():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_17'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test17_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_18():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_18'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test18_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_19():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_19'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test19_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_20():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_20'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test20_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_21():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_21'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test21_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)

    @staticmethod
    def test_22():
        ground_truth_dir = 'genonets/test/data/ground_truth/epistasis/test_22'

        with tempfile.TemporaryDirectory(prefix='test_epistasis_') as data_dir:
            cmd_args = [
                '--alphabet=Protein',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/epistasis/test22_input.tsv',
                f'--output-path={data_dir}'
            ]

            TestEpistasis.run_test(cmd_args, ground_truth_dir, data_dir)
