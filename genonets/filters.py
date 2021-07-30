"""
    genonets_filters
    ~~~~~~~~~~~~~~~~

    Contains filters used throughout the package.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""


class WriterFilter:
    def __init__(self):
        pass

    # Ordered list of genotype sets, initialized by the Genonets object
    ORDERED_GENOTYPE_SETS = None

    # Network level attribute to order map
    ordered_net_attrib = {
        "Robustness": 1,
        "Evolvability": 2,
        "Evolvability_targets": 3,
        "Interface_edges": 4,
        "Accessibility": 5,
        "Neighbor_abundance": 6,
        "Diversity_index": 7,
        "Number_of_genotype_networks": 8,
        "Genotype_network_sizes": 9,
        "Size_of_dominant_genotype_network": 10,
        "Proportional_size_of_dominant_genotype_network": 11,
        "Diameter": 12,
        "Edge_density": 13,
        "Assortativity": 14,
        "Average_clustering_coefficient_of_dominant_genotype_network": 15,
        "Number_of_peaks": 16,
        "Peaks": 17,
        "Number_of_squares": 18,
        "Magnitude_epistasis": 19,
        "Simple_sign_epistasis": 20,
        "Reciprocal_sign_epistasis": 21,
        "Summit": 22,
        "Ratio_of_accessible_mutational_paths": 23,
        "Ratio_of_overlapping_genotype_sets": 24,
        "Overlapping_genotype_sets": 25,
    }

    # Genotype level attribute to order map
    ordered_seq_attrib = {
        "Robustness": 1,
        "Evolvability": 2,
        "Evolvability_targets": 3,
        "Evolves_to_genotypes_in": 4,
        "Overlaps_with_genotypes_in": 5,
        "Coreness": 6,
        "Clustering_coefficient": 7,
        "Distance from Summit": 8,
        "Shortest_path_length": 9,
        "Shortest_accessible_path_length": 10,
        "Accessible_paths_through": 11,
        "Accessible_paths_from": 12
    }

    @staticmethod
    def gmlAttribsToIgnore(level):
        if level == "network":
            attrs = [
                "Evolvability_targets",
                "SqrEpi_list",
                "diameterPath_list",
                "Squares_list"
            ]
        elif level == "vertex":
            attrs = [
                "label",
                "pathsToSummit",
                "VtxToSqrs"
            ]

        return attrs

    @staticmethod
    def netAttribsToIgnore():
        return [
            "name",
            "SqrEpi_list",
            "diameterPath_list",
            "Squares_list",
        ]

    @staticmethod
    def seqAttribsToIgnore():
        return [
            "name",
            "sequences",
            "label",
            "escores",
            "pathsToSummit",
            "VtxToSqrs"
        ]

    @staticmethod
    def net_attribute_to_order(attribute):
        try:
            return WriterFilter.ordered_net_attrib[attribute]
        except KeyError:
            # Any custom attribute added by the user will be assigned
            # this arbitrary high value, so that the user-defined
            # attributes are placed at the very end.
            return 1000

    @staticmethod
    def seq_attribute_to_order(attribute):
        try:
            return WriterFilter.ordered_seq_attrib[attribute]
        except KeyError:
            # Any custom attribute added by the user will be assigned
            # this arbitrary high value, so that the user-defined
            # attributes are placed at the very end.
            return 1000

    @staticmethod
    def genotype_set_to_order(genotype_set):
        try:
            return WriterFilter.ORDERED_GENOTYPE_SETS.index(genotype_set)
        except ValueError:
            # If the requested genotype set was not found in the list, just
            # return an arbitrarily high value.
            return 9999999
