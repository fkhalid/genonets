
import os
import sys
import csv
import argparse

import numpy as np
import pandas as pd

# Enable support for handling large files
csv.field_size_limit(sys.maxsize)

attr_to_type_map = {
    'Peaks': dict,
    'Evolvability_targets': list,
    'Genotype_network_sizes': list,
    'Evolves_to_genotypes_in': list,
    'Overlaps_with_genotypes_in': list,
    'Ratio_of_accessible_mutational_paths': dict,
    'Epistasis squares': dict,
    'Interface_edges': dict,
}

genotype_attr_to_action_map = {
    'Sequence': str,
    'Robustness': float,
    'Evolvability': float,
    'Evolvability_targets': eval,
    'Evolves_to_genotypes_in': eval,
    'Overlaps_with_genotypes_in': eval,
    'Coreness': float,
    'Clustering_coefficient': float,
    'Distance from Summit': int,
    'Accessible_paths_through': int,
    'Accessible_paths_from': int,
    'Shortest_accessible_path_length': int,
    'Shortest_path_length': int,
}

genotype_set_attr_to_action_map = {
    'Genotype_set': str,
    'Robustness': float,
    'Evolvability': float,
    'Evolvability_targets': eval,
    'Accessibility': float,
    'Neighbor_abundance': float,
    'Diversity_index': float,
    'Diameter': float,
    'Number_of_genotype_networks': int,
    'Genotype_network_sizes': eval,
    'Size_of_dominant_genotype_network': int,
    'Proportional_size_of_dominant_genotype_network': float,
    'Edge_density': float,
    'Assortativity': float,
    'Average_clustering_coefficient_of_dominant_genotype_network': float,
    'Number_of_peaks': int,
    'Peaks': eval,
    'Number_of_squares': int,
    'Magnitude_epistasis': float,
    'Simple_sign_epistasis': float,
    'Reciprocal_sign_epistasis': float,
    'Ratio_of_overlapping_genotype_sets': float,
    'Overlapping_genotype_sets': eval,
    'Ratio_of_accessible_mutational_paths': eval,
    'Epistasis squares': eval,
    'Interface_edges': eval,
    'Summit': eval,
}


class FileComparator:
    @staticmethod
    def reader_from_file(filename):
        data_file = open(filename, 'r')

        return csv.DictReader(data_file, delimiter='\t')

    @staticmethod
    def get_filenames_with_substring(data_dir, substring):
        return [
            f for f in os.listdir(data_dir)
            if substring in f
        ]

    @staticmethod
    def field_names_from_reader(r):
        return list(set(r.fieldnames) - {''})

    @staticmethod
    def get_value(attr, value, attr_map):
        try:
            return attr_map[attr](value)
        except NameError:
            if attr_to_type_map[attr] == list:
                return []
            elif attr_to_type_map[attr] == dict:
                return {}

    @staticmethod
    def floats_match(v1, v2):
        return np.isclose(
            [v1], [v2], rtol=1e-09, atol=1e-09, equal_nan=True
        )

    @staticmethod
    def are_floats(v1, v2):
        return isinstance(v1, float) or isinstance(v2, float)

    @staticmethod
    def data_from_reader(rd, attr_map):
        field_names = FileComparator.field_names_from_reader(rd)

        data = {f: [] for f in field_names}

        for row in rd:
            for f in field_names:
                data[f].append(
                    FileComparator.get_value(f, row[f], attr_map)
                )

        return data

    @staticmethod
    def compare_dataframes(df_1, df_2):
        dataframes_match = True

        # Check if the number of columns is the same
        if not df_1.columns.sort_values().equals(df_2.columns.sort_values()):
            print('Mismatch in columns.')

            c_1 = set(df_1.columns.sort_values().values)
            c_2 = set(df_2.columns.sort_values().values)

            print(f'Columns in source not in target: {c_1 - c_2}')
            print(f'Columns in target not in source: {c_2 - c_1}')

            # Columns in source that are not in target should be
            # removed from source so that we do not try to read a
            # non-existent columns from target
            if c_1 - c_2:
                df_1 = df_1.drop(columns=list(c_1 - c_2))

            dataframes_match = False

        # Check if index column has the same values
        if not df_1.index.sort_values().equals(df_2.index.sort_values()):
            print('Mismatch in index values.')

            c_1 = set(df_1.index.sort_values().values)
            c_2 = set(df_2.index.sort_values().values)

            print(f'Index values in source not in target: {c_1 - c_2}')
            print(f'Index values in target not in source: {c_2 - c_1}')

            # Remove any index values from df_1 that do not exist in
            # df_2
            if c_1 - c_2:
                df_1 = df_1.drop(list(c_1 - c_2), axis=0)

            dataframes_match = False

        for k in df_1.index.values:
            for c in df_1.columns.values:
                if c == 'Epistasis squares':
                    continue
                if FileComparator.are_floats(df_1[c][k], df_2[c][k]):
                    if not FileComparator.floats_match(df_1[c][k], df_2[c][k]):
                        print(
                            f'Mismatch '
                            f'-- Key: {k} '
                            f'-- Column: {c} '
                            f'-- Source value: {df_1[c][k]} '
                            f'-- Target value: {df_2[c][k]}'
                        )
                        dataframes_match = False
                elif isinstance(df_1[c][k], list):
                    l_1 = set(df_1[c][k])
                    l_2 = set(df_2[c][k])

                    if l_1 ^ l_2:
                        print(
                            f'Mismatch '
                            f'-- Key: {k} '
                            f'-- Column: {c} '
                            f'-- Source value: {df_1[c][k]} '
                            f'-- Target value: {df_2[c][k]}'
                        )
                        dataframes_match = False
                elif isinstance(df_1[c][k], dict):
                    if df_1[c][k].keys() != df_2[c][k].keys():
                        print(
                            f'Mismatch '
                            f'-- Key: {k} '
                            f'-- Column: {c} '
                            f'-- Source value: {df_1[c][k].keys()} '
                            f'-- Target value: {df_2[c][k].keys()} '
                            f'-- Difference: '
                            f'{set(df_1[c][k].keys()) ^ set(df_2[c][k].keys())}'
                        )
                        dataframes_match = False

                    diff = set(df_1[c][k].keys()) ^ set(df_2[c][k].keys())

                    keys = set(df_1[c][k].keys()) - diff

                    for ky in keys:
                        v_1 = df_1[c][k][ky]
                        v_2 = df_2[c][k][ky]

                        if isinstance(v_1, list):
                            if set(v_1) ^ set(v_2):
                                print(
                                    f'Mismatch '
                                    f'-- Key: {k} '
                                    f'-- Column: {c} '
                                    f'-- Source value: {df_1[c][k]} '
                                    f'-- Target value: {df_2[c][k]}'
                                )
                                dataframes_match = False
                        elif isinstance(v_1, float) or isinstance(v_2, float):
                            if not FileComparator.are_floats(v_1, v_2):
                                print(
                                    f'Mismatch '
                                    f'-- Key: {k} '
                                    f'-- Column: {c} '
                                    f'-- Source value: {df_1[c][k]} '
                                    f'-- Target value: {df_2[c][k]}'
                                )
                                dataframes_match = False
                        else:
                            if v_1 != v_2:
                                print(
                                    f'Mismatch '
                                    f'-- Key: {k} '
                                    f'-- Column: {c} '
                                    f'-- Source value: {df_1[c][k]} '
                                    f'-- Target value: {df_2[c][k]}'
                                )
                                dataframes_match = False
                else:
                    if df_1[c][k] != df_2[c][k]:
                        print(
                            f'Mismatch '
                            f'-- Key: {k} '
                            f'-- Column: {c} '
                            f'-- Source value: {df_1[c][k]} '
                            f'-- Target value: {df_2[c][k]}'
                        )
                        dataframes_match = False

        return dataframes_match

    @staticmethod
    def compare_files(f1, f2, attr_map, index_column):
        # Create a reader for each file
        r1 = FileComparator.reader_from_file(f1)
        r2 = FileComparator.reader_from_file(f2)

        # Load data from readers into dicts
        d1 = FileComparator.data_from_reader(r1, attr_map)
        d2 = FileComparator.data_from_reader(r2, attr_map)

        # Construct Pandas data frames from dicts
        df1 = pd.DataFrame(d1)
        df2 = pd.DataFrame(d2)

        # Set the index column
        df1 = df1.set_index(index_column)
        df2 = df2.set_index(index_column)

        return FileComparator.compare_dataframes(df1, df2)


def compare_genotype_set_measures(source_dir, target_dir):
    filename = 'Genotype_set_measures.csv'

    f1 = os.path.join(source_dir, filename)
    f2 = os.path.join(target_dir, filename)

    print(f'{filename}: Starting ...')

    files_match = FileComparator.compare_files(
        f1, f2,
        attr_map=genotype_set_attr_to_action_map,
        index_column='Genotype_set'
    )

    print(f'{filename}: Done.')

    return files_match


def compare_genotype_measures(source_dir, target_dir):
    files_match = True

    # Pattern with which to identify the required files
    file_id_string = '_genotype_measures'

    source_filenames = FileComparator.get_filenames_with_substring(
        source_dir, file_id_string)
    target_filenames = FileComparator.get_filenames_with_substring(
        target_dir, file_id_string)

    print(f'Genotype measure files in source but not in target: '
          f'{set(source_filenames) - set(target_filenames)}')
    print(f'Genotype measure files in target but not in source: '
          f'{set(target_filenames) - set(source_filenames)}')

    # Only those files should be considered which are present in both
    # the source and target directories
    files_to_compare = set(source_filenames) & set(target_filenames)

    for f in files_to_compare:
        print(f'{f}: Starting ...')

        f1 = os.path.join(source_dir, f)
        f2 = os.path.join(target_dir, f)

        result = FileComparator.compare_files(
            f1, f2,
            attr_map=genotype_attr_to_action_map,
            index_column='Sequence'
        )

        if not result:
            files_match = False

        print(f'{f}: Done.')

    return files_match


def compare_overlap_results(source_dir, target_dir):
    files_match = True
    filename = 'Genotype_set_overlap.csv'

    if filename not in os.listdir(source_dir):
        print('Overlap measures not found in the source directory.')
        return

    print(f'{filename}: Starting ...')

    f1 = os.path.join(source_dir, filename)
    f2 = os.path.join(target_dir, filename)

    df1 = pd.read_csv(f1, sep='\t')
    df2 = pd.read_csv(f2, sep='\t')

    df1 = df1.rename(columns={'Unnamed: 0': 'Genotype_set'})
    df2 = df2.rename(columns={'Unnamed: 0': 'Genotype_set'})

    df1 = df1.set_index('Genotype_set')
    df2 = df2.set_index('Genotype_set')

    result = FileComparator.compare_dataframes(df1, df2)

    if not result:
        files_match = False

    print(f'{filename}: Done.')

    return files_match


def parse_command_line_args():
    # Initialize the parser object
    parser = argparse.ArgumentParser()

    # ---------------- Mandatory arguments ------------------ #

    parser.add_argument(
        dest='source_dir', default=None,
        help='Directory from which to load ground truth files.'
    )

    parser.add_argument(
        dest='target_dir', default=None,
        help='Directory from which to load the files which need to be '
             'compared with the ground truth files.'
    )

    return parser.parse_args()


if __name__ == '__main__':
    cmd_args = parse_command_line_args()

    compare_genotype_set_measures(cmd_args.source_dir, cmd_args.target_dir)

    compare_genotype_measures(cmd_args.source_dir, cmd_args.target_dir)

    compare_overlap_results(cmd_args.source_dir, cmd_args.target_dir)
