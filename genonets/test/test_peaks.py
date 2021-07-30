
import tempfile

import genonets.test.utils as utils
import genonets.test.compare_result_files as comparator

from genonets.cmdl_handler import CmdParser
from genonets.interface import Genonets
from genonets.constants import AnalysisConstants as Ac


class TestPeaks:
    @staticmethod
    def test_1():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_1'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test1_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_2():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_2'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test2_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_3():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_3'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test3_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_4():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_4'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test4_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_5():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_5'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test5_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_6():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_6'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test6_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_7():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_7'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test7_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_8():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_8'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test8_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_9():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_9'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test9_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_10():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_10'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test10_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_11():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_11'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test11_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_12():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_12'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test12_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

    @staticmethod
    def test_13():
        ground_truth_dir = 'genonets/test/data/ground_truth/peaks/test_13'

        with tempfile.TemporaryDirectory(prefix='test_peaks_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/peaks/test13_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.PEAKS])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )
