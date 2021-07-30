"""
    This module implements the functions required for computation of
    epistasis in a given network (when viewed as a fitness landscape).

    :license: MIT, see LICENSE for more details.
"""

import os
import math
import shutil
import collections

from tqdm import tqdm

from genonets.errors import FileError
from genonets.constants import EpistasisConstants as Constants


AlternateNetworkData = collections.namedtuple(
    'AlternateNetworkData',
    ['network', 'bit_manip', 'net_utils']
)
AlternateNetworkData.__doc__ = """
    Contains the alternate network, the corresponding bit manipulator and 
    network utils.
"""


class Buffer:
    """
    A simple buffer that flushes its content to a new file and  self
    re-initializes whenever it is completely filled, i.e, whenever
    the No. of elements in the buffer reaches the size provided. The
    context manager protocol is implemented.

    """

    def __init__(self, size, out_dir, genotype_to_score_map, verbose=False):
        """
        Initializer.

        :param size: (int) Total No. of squares to hold in the buffer before
                    flushing the content to the disk.
        :param out_dir: (str) The path to the output directory.
        :param genotype_to_score_map: (dict) {genotype: score}, for all
                    genotypes in the network.
        :param verbose: (bool) Flag to indicate whether or not verbose
                    output is required.

        """

        self._size = size
        self._out_dir = out_dir
        self._verbose = verbose
        self._squares_dir = os.path.join(out_dir, 'epistasis_squares')
        self._genotype_to_score_map = genotype_to_score_map

        self._create_dir(self._squares_dir)

        # Data structure to hold squares
        self._no_epistasis_data = []
        self._sign_epistasis_data = []
        self._magnitude_epistasis_data = []
        self._reciprocal_epistasis_data = []

        # File names
        no_epistasis_filename = os.path.join(
            self._squares_dir, 'no_epistasis.txt')
        sign_epistasis_filename = os.path.join(
            self._squares_dir, 'simple_sign_epistasis.txt')
        magnitude_epistasis_filename = os.path.join(
            self._squares_dir, 'magnitude_sign_epistasis.txt')
        reciprocal_epistasis_filename = os.path.join(
            self._squares_dir, 'reciprocal_sign_epistasis.txt')

        # File handles
        self._f_no_epistasis = open(no_epistasis_filename, 'w')
        self._f_sign_epistasis = open(sign_epistasis_filename, 'w')
        self._f_magnitude_epistasis = open(magnitude_epistasis_filename, 'w')
        self._f_reciprocal_epistasis = open(reciprocal_epistasis_filename, 'w')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _create_dir(self, path):
        """
        Creates the directory specified by the given path.

        :param path: (str) Path referring to the directory to create.

        """

        if self._verbose:
            print(f'Creating directory: {path}', end=' ... ', flush=True)

        if not os.path.exists(path):
            try:
                os.makedirs(path)
            except os.error:
                raise FileError(
                    f'Error: Failed to create directory ({path})'
                )
        else:
            if self._verbose:
                print('already exists', end=' ... ', flush=True)

        if self._verbose:
            print('Done.', flush=True)

    def push(self, square: list, epi_type: Constants):
        """
        Appends the given square to the buffer. If this results in a full
        buffer, it is flushed.

        :param square: (list) The square to append.
        :param epi_type: (EpistasisConstants) The type of epistasis exhibited
                    by the square.

        """

        # Reformat the square such that each element is a tuple of the form
        # (genotype, score)
        square_with_scores = [
            (genotype, self._genotype_to_score_map[genotype])
            for genotype in square
        ]

        # Swap the last two elements so that each genotype is a 1-neighbor of
        # the preceding genotype
        square_with_scores[-2], square_with_scores[-1] = \
            square_with_scores[-1], square_with_scores[-2]

        if epi_type == Constants.NO_EPISTASIS:
            self._no_epistasis_data.append(square_with_scores)
        elif epi_type == Constants.SIGN:
            self._sign_epistasis_data.append(square_with_scores)
        elif epi_type == Constants.MAGNITUDE:
            self._magnitude_epistasis_data.append(square_with_scores)
        else:
            self._reciprocal_epistasis_data.append(square_with_scores)

        if self._should_flush():
            self._flush()

    def close(self):
        """
        Serves as the destructor.

        """

        self._flush()

        self._f_no_epistasis.close()
        self._f_sign_epistasis.close()
        self._f_magnitude_epistasis.close()
        self._f_reciprocal_epistasis.close()

        shutil.make_archive(
            os.path.join(self._out_dir, 'epistasis_squares'),
            'zip',
            self._squares_dir
        )

        shutil.rmtree(self._squares_dir)

    def _should_flush(self):
        """
        Determines whether or not to flush the buffer contents to disk.

        :return: (bool) True if the contents should be flushed; False
                otherwise.

        """

        num_total_elements = len(self._no_epistasis_data) \
                             + len(self._sign_epistasis_data) \
                             + len(self._magnitude_epistasis_data) \
                             + len(self._reciprocal_epistasis_data)

        return True if num_total_elements >= self._size else False

    def _flush(self):
        """
        Flushes all contents of the buffer to disk.

        """

        self._save()

        del self._no_epistasis_data
        del self._sign_epistasis_data
        del self._magnitude_epistasis_data
        del self._reciprocal_epistasis_data

        self._no_epistasis_data = []
        self._sign_epistasis_data = []
        self._magnitude_epistasis_data = []
        self._reciprocal_epistasis_data = []

    @staticmethod
    def _save_to_file(data, file):
        """
        Writes the given data list to disk using the given file handle.

        :param data: (list) List of squares.
        :param file: (file object) File handle for the file in which to write
                    write the squares.

        """

        if data:
            for square in data:
                file.write(f'{str(square)}\n')

    def _save(self):
        """
        Write all buffered squares to file.

        """

        self._save_to_file(self._no_epistasis_data, self._f_no_epistasis)
        self._save_to_file(self._sign_epistasis_data, self._f_sign_epistasis)
        self._save_to_file(self._magnitude_epistasis_data, self._f_magnitude_epistasis)
        self._save_to_file(self._reciprocal_epistasis_data, self._f_reciprocal_epistasis)


class EpistasisAnalyzer:
    """
    Encapsulates the computation required for identification of
    squares and computation of epistasis. It is assumed that the given
    network is a single connected component; typically the dominant
    component of a genotype network.

    """

    def __init__(self, network, net_utils, genotype_to_score_map,
                 delta, alternate_net_data, verbose,
                 save_squares=False, out_dir=None):
        """
        Initializes the object.

        :param network: (igraph.Graph) The network within which to
                    identify peaks. The network must be connected, i.e.,
                    it must not consist of more than one connected
                    component.
        :param net_utils: (genonets.graph_utils.NetworkBuilder)
                    Reference to the object providing the required
                    network manipulation functions.
        :param genotype_to_score_map: (dict) {genotype: score}, for all
                    genotypes in the network.
        :param delta: (float) The experimental noise value to consider
                    when evaluating each genotype for optimality.
        :param verbose: (bool) Flag to indicate whether or not verbose
                    output is required.
        :param save_squares: (bool) Flag to indicate whether or not to write
                    squares to disk.
        :param out_dir: (str) Path to the directory in which to store squares.

        """

        self._delta = delta
        self._network = network
        self._verbose = verbose
        self._out_dir = out_dir
        self._net_utils = net_utils
        self._save_squares = save_squares
        self._bit_manip = net_utils.bitManip
        self._genotype_to_score_map = genotype_to_score_map

        # Data structures required when an alternate network has to be
        # considered
        self._alternate_network = None
        self._alternate_bit_manip = None
        self._alternate_net_utils = None
        self._alternate_neighbor_map = None

        # Initialize map {Epistasis type: square count}
        self._epi_type_to_count_map = {
            Constants.MAGNITUDE: 0,
            Constants.SIGN: 0,
            Constants.RECIPROCAL_SIGN: 0,
            Constants.NO_EPISTASIS: 0
        }

        if self._verbose:
            print('\nConstructing neighbor map ...')

        if self._verbose:
            iterable = tqdm(network.vs['sequences'])
        else:
            iterable = network.vs['sequences']

        # Neighbor map
        self._neighbor_map = {
            s: set(self._net_utils.getNeighborSequences(s, self._network))
            for s in iterable
        }

        # If the alternate network is to be consulted,
        if alternate_net_data is not None:
            self._alternate_network = alternate_net_data.network
            self._alternate_bit_manip = alternate_net_data.bit_manip
            self._alternate_net_utils = alternate_net_data.net_utils

            if self._verbose:
                print('Constructing alternate neighbor map ...')

            if self._verbose:
                iterable = tqdm(self._alternate_network.vs['sequences'])
            else:
                iterable = self._alternate_network.vs['sequences']

            self._alternate_neighbor_map = {
                s: set(self._alternate_net_utils.getNeighborSequences(s, self._alternate_network))
                for s in iterable
            }

    @property
    def num_squares(self) -> int:
        """
        Return the total No. of squares identified in the network.

        :return: No. of squares found.

        """

        return sum(self._epi_type_to_count_map.values())

    def compute_epistasis_ratios(self) -> dict:
        """
        Processes all squares and computes the ratio of each epistasis
        type in the network.

        :return: dict of the form {
                    epistasis type: ratio of squares of this type
                }

        """

        # Identify all squares and populate epistasis counts
        if self._save_squares:
            with Buffer(5_00_000, self._out_dir, self._genotype_to_score_map) as b:
                self._compute_and_analyze_squares(b)
        else:
            self._compute_and_analyze_squares(None)

        # Map {epistasis type: ratio}
        epi_type_to_ratio_map = {}

        # Convert counts into ratios
        for epi_type in self._epi_type_to_count_map:
            try:
                epi_type_to_ratio_map[epi_type] = \
                    float(self._epi_type_to_count_map[epi_type]) / float(self.num_squares)
            except ZeroDivisionError:
                epi_type_to_ratio_map[epi_type] = 0

        if self._verbose:
            print(f'No. of squares: {self.num_squares}')

        return epi_type_to_ratio_map

    def _epistasis_for(self, square: list) -> Constants:
        """
        Computes epistasis for the given square.

        :param square: Vertex Ids of the genotypes that constitute the
                    square.

        :return: Type of epistasis corresponding to the given square.

        """

        # Scores for all corners of the square
        score_AB = self._genotype_to_score_map[
            self._network.vs[square[0]]['name']]
        score_Ab = self._genotype_to_score_map[
            self._network.vs[square[1]]['name']]
        score_aB = self._genotype_to_score_map[
            self._network.vs[square[2]]['name']]
        score_ab = self._genotype_to_score_map[
            self._network.vs[square[3]]['name']]

        # Magnitude of epistasis
        epsilon = math.fsum([score_AB, score_ab, -score_Ab, -score_aB])

        # Determine the type of epistasis
        if math.fabs(epsilon) <= (2.0 * self._delta):
            return Constants.NO_EPISTASIS
        else:
            dE_ab_Ab = math.fsum([score_Ab, -score_ab])
            dE_ab_aB = math.fsum([score_aB, -score_ab])

            if dE_ab_Ab >= -self._delta and dE_ab_aB >= -self._delta:
                return Constants.MAGNITUDE
            elif dE_ab_Ab >= -self._delta or dE_ab_aB >= -self._delta:
                return Constants.SIGN
            else:
                return Constants.RECIPROCAL_SIGN

    def _get_sorted_genotypes(self) -> list:
        """
        Creates a list of all genotypes in the network; sorted in descending
        order of scores.

        :return: List of all genotypes sorted in descending order of scores.

        """

        # Zip genotypes and scores into tuples
        genotype_score_pairs = list(zip(
            self._network.vs['sequences'],
            self._network.vs['escores']
        ))

        # Sort the list in descending order of scores
        genotype_score_pairs.sort(key=lambda x: x[1], reverse=True)

        # Extract the list of sorted genotypes; discard the scores
        genotypes, _ = zip(*genotype_score_pairs)

        return genotypes

    def _compute_and_analyze_squares(self, squares_buffer):
        """
        Identifies all squares, and for each square, determines the
        type of epistasis.

        Populates the epistasis-to-squares-count map.

        """

        if self._verbose:
            print('Constructing and analyzing squares ...')

        # Set of genotypes that have already been processed
        processed_genotypes = set()

        # List of all genotypes
        genotypes = self._get_sorted_genotypes()

        iterable = tqdm(genotypes) if self._verbose else genotypes

        # For each genotype,
        for genotype in iterable:
            # All 1-neighbors
            neighbors = list(
                self._neighbor_map[genotype] - processed_genotypes
            )

            # If the number of neighbors is less than two, there's no
            # point in continuing any further
            if len(neighbors) < 2:
                # Move on to the next genotype
                continue

            # TODO: Figure out if we need to check whether or not the two
            #   genotypes in the pair are also not 1-neighbors in the
            #   alternate network ...

            # Construct all possible pairs of neighbors, where
            # symmetric pairs are considered only once. Also, pairs
            # that neighbor each other are not considered.
            pairs = [
                (neighbors[i], neighbors[j])
                for i in range(len(neighbors) - 1)
                for j in range(i + 1, len(neighbors))
                if neighbors[j] not in self._neighbor_map[neighbors[i]]
                and neighbors[i] not in self._neighbor_map[neighbors[j]]
            ]

            # Process all pairs
            while pairs:
                pair = pairs.pop()

                # Set of genotypes that are 1-neighbors of both the
                # genotypes in the pair, less the focal genotype and
                # all the genotypes that have already been processed
                # as focal genotypes.
                common_neighbors = ((self._neighbor_map[pair[0]] &
                                     self._neighbor_map[pair[1]]) -
                                    {genotype}) - processed_genotypes

                # If an alternate neighbor map is available,
                if self._alternate_neighbor_map is not None:
                    # Exclude 1-neighbors in the alternate network
                    common_neighbors -= self._alternate_neighbor_map[genotype]

                # For each common neighbor,
                for node in common_neighbors:
                    # For a node that is a 1-neighbor of the focal
                    # genotype, it does not make sense to compute
                    # epistasis, since a single mutation from parent is
                    # sufficient.
                    if node not in neighbors:
                        # Construct the square
                        square = [
                            self._network.vs.find(genotype).index,
                            self._network.vs.find(pair[0]).index,
                            self._network.vs.find(pair[1]).index,
                            self._network.vs.find(node).index
                        ]

                        # Determine the type of epistasis represented
                        # in the square
                        epi_type = self._epistasis_for(square)

                        # Increment the corresponding count
                        self._epi_type_to_count_map[epi_type] += 1

                        # Store the square
                        if squares_buffer is not None:
                            squares_buffer.push(
                                [genotype, pair[0], pair[1], node], epi_type
                            )

            # Mark the focal genotype as processed
            processed_genotypes.add(genotype)
