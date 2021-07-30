
import tempfile

import genonets.test.utils as utils
import genonets.test.compare_result_files as comparator

from genonets.cmdl_handler import CmdParser
from genonets.interface import Genonets
from genonets.constants import AnalysisConstants as Ac


class TestInterfaceEdgesForDominant:
    @staticmethod
    def test_1():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'dominant_only/test_1'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/dominant_only/test1_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
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
    def test_2():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'dominant_only/test_2'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/dominant_only/test2_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
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
    def test_3():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'dominant_only/test_3'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=RNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/dominant_only/test3_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
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
    def test_4():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'dominant_only/test_4'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/dominant_only/test4_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
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
    def test_5():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'dominant_only/test_5'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--include-indels',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/dominant_only/test5_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

            assert comparator.compare_genotype_measures(
                ground_truth_dir, data_dir
            )


class TestInterfaceEdgesForAll:
    @staticmethod
    def test_1():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'all_components/test_1'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--use-all-components',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/all_components/test1_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
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
    def test_2():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'all_components/test_2'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--use-all-components',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/all_components/test2_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
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
    def test_3():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'all_components/test_3'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=RNA',
                '--tau=0.0',
                '--use-all-components',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/all_components/test3_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
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
    def test_4():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'all_components/test_4'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--use-all-components',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/all_components/test4_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
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
    def test_5():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'all_components/test_5'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                '--include-indels',
                '--use-all-components',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/all_components/test5_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
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
    def test_6():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'all_components/test_6'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=RNA',
                '--tau=0.0',
                # '--use-all-components',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/all_components/test6_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
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
    def test_7():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'all_components/test_7'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=DNA',
                '--tau=0.0',
                # '--use-all-components',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/all_components/test7_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
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
    def test_8():
        ground_truth_dir = 'genonets/test/data/ground_truth/interface_edges/' \
                           'all_components/test_8'

        with tempfile.TemporaryDirectory(prefix='test_interface_edges_') as data_dir:
            cmd_args = [
                '--alphabet=RNA',
                '--tau=0.0',
                # '--use-all-components',
                '--input-file=genonets/test/data/inputs/interface_edges'
                '/all_components/test8_input.tsv',
                f'--output-path={data_dir}'
            ]

            args = CmdParser(arguments=cmd_args).get_args()

            gn = Genonets(args)
            gn.create()
            gn.analyze(analyses=[Ac.EVOLVABILITY])
            gn.save_network_results()
            gn.save_genotype_results()

            assert utils.num_files_matches(ground_truth_dir, data_dir)

            assert comparator.compare_genotype_set_measures(
                ground_truth_dir, data_dir
            )

            assert comparator.compare_genotype_measures(
                ground_truth_dir, data_dir
            )
