"""
    path_functions
    ~~~~~~~~~~~~~~

    Contains functions used for computation of accessible mutational paths.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

from tqdm import tqdm

from genonets.utils import Utils


class PathAnalyzer:
    # Constructor
    def __init__(self, network, netUtils, delta, verbose, genotype_to_score_map):
        # Reference to the network on which to perform this
        # analysis
        self._network = network

        # Get a reference to the NetworkUtils object
        self._net_utils = netUtils

        # Get a reference to the BitSeqManipulator in use
        self._bit_manip = self._net_utils.bitManip

        # Keep a copy of the delta value
        self._delta = delta

        # Flag to indicate whether output should be verbose
        self._verbose = verbose

        # Map: {key=genotype: value=escore}
        self._genotype_to_score_map = genotype_to_score_map

        # Summit vertex, to be populated later.
        self._summit_sequences = self._compute_summits()

        # TODO: Factor this out from all analyses that require it, so that is
        #   computed only once per gentoype set ...
        # Neighbor map
        if self._verbose:
            print('\nConstructing neighbor map ...')
            iterable = tqdm(self._network.vs['sequences'])
        else:
            iterable = self._network.vs['sequences']
        self._neighbor_map = {
            s: set(self._net_utils.getNeighborSequences(s, self._network))
            for s in iterable
        }

    @property
    def summit_sequences(self):
        return self._summit_sequences

    def _compute_summits(self):
        sorted_tuples = Utils.getSortedSeqEscArr(
            self._network,
            self._bit_manip.seqLength,
            sortOrder="ascending"
        )

        # Highest score
        max_score = sorted_tuples[-1]['escore']

        # Find the minimal index i such that sortedArr[i]['escore'] = maxScore
        i = -1
        while True:
            if sorted_tuples[i-1]['escore'] == max_score:
                i -= 1
            else:
                break

        return [x['sequence'] for x in sorted_tuples[i:]]

    def compute_shortest_paths(self):
        queue = []

        visited = {genotype: False for genotype in self._neighbor_map}
        shortest_path_length = {genotype: 0 for genotype in self._neighbor_map}
        min_scores_on_acc_paths = {genotype: {} for genotype in self._neighbor_map}
        num_accessible_paths = {genotype: 0 for genotype in self._neighbor_map}
        num_inaccessible_paths = {genotype: 0 for genotype in self._neighbor_map}

        # the global peak(s)
        for summit in self._summit_sequences:
            visited[summit] = True

        # score of the global peak(s)
        global_max = self._genotype_to_score_map[self._summit_sequences[0]]

        # add all neighbours of summit to the queue, update the helper arrays
        for summit in self._summit_sequences:
            neighbours = self._neighbor_map[summit]
            for N in neighbours:
                if not visited[N]:
                    queue.append(N)
                    shortest_path_length[N] = 1
                    min_scores_on_acc_paths[N] = {global_max: 1}
                    visited[N] = True

        # process the rest of the network using BFS
        while queue:
            # get the first element
            V = queue.pop(0)
            V_score = self._genotype_to_score_map[V]

            # update the number of accessible and inaccessible paths
            # and create a list with minimal scores (including V) of all
            # accessible paths from summit to V
            for min_score_on_path, count in min_scores_on_acc_paths[V].items():
                if V_score <= min_score_on_path + self._delta:
                    num_accessible_paths[V] += count
                else:
                    num_inaccessible_paths[V] += count


            # add all not yet discovered neighbours of V to the queue, update
            # the helper arrays
            neighbours = self._neighbor_map[V]
            for N in neighbours:
                # N discovered for the first time
                if not visited[N]:
                    queue.append(N)
                    shortest_path_length[N] = shortest_path_length[V] + 1
                    visited[N] = True
                # if N->V is part of some shortest paths N->summit
                if shortest_path_length[N] == shortest_path_length[V] + 1:
                    # update the minimal scores seen on all paths going from
                    # summit to N
                    for min_score_on_path, count in min_scores_on_acc_paths[V].items():
                        if V_score <= min_score_on_path + self._delta:
                            tmp_min_score = min(min_score_on_path, V_score)
                            if tmp_min_score in min_scores_on_acc_paths[N]:
                                min_scores_on_acc_paths[N][tmp_min_score] += count
                            else:
                                min_scores_on_acc_paths[N][tmp_min_score] = count
                    # update the number of inaccessible paths N->summit
                    num_inaccessible_paths[N] += num_inaccessible_paths[V]

        # compute the ratios of accessible paths
        max_path_length = max(shortest_path_length.values())
        ratio_accessible = {}
        for i in range(2, max_path_length + 1):
            num_accessible = sum([
                num_accessible_paths[k]
                for k in num_accessible_paths
                if shortest_path_length[k] == i
            ])
            num_inaccessible = sum([
                num_inaccessible_paths[k]
                for k in num_inaccessible_paths
                if shortest_path_length[k] == i
            ])
            ratio_accessible[i] = num_accessible / (num_accessible + num_inaccessible)

        results = {
            "Ratio_of_accessible_mutational_paths": ratio_accessible,
            "Accessible_paths_from": list(num_accessible_paths.values()),
            "Shortest_path_length": list(shortest_path_length.values())
        }

        return results

    def compute_indirect_paths(self):
        visited = {genotype: False for genotype in self._neighbor_map}
        shortest_path_length = {genotype: -1 for genotype in self._neighbor_map}
        maximal_min_score_on_acc_path = {genotype: [] for genotype in self._neighbor_map}
        queue = []

        # the summit sequence(s)
        for summit in self._summit_sequences:
            visited[summit] = True
            shortest_path_length[summit] = 0

        # score of the global peak
        global_max = self._genotype_to_score_map[self._summit_sequences[0]]
        for summit in self._summit_sequences:
            maximal_min_score_on_acc_path[summit] = global_max

        # add the summit(s) to the queue
        for summit in self._summit_sequences:
            queue.append([summit, 0, global_max])

        # process the rest of the network
        while queue:
            # get the first element
            V, path_length, min_score = queue.pop(0)

            # get all neighbors of V
            neighbours = self._neighbor_map[V]
            for N in neighbours:
                N_score = self._genotype_to_score_map[N]

                # we add N to the queue if:
                # 1. N not yet visited and the path after adding the (V, N)
                # edge is accessible
                # 2. N already visited but the min_score on the current path is
                # larger (i.e., we might be able to
                #   discover some new accessible paths)
                cond1 = not visited[N] and min_score >= N_score - self._delta
                cond2 = visited[N] and min(min_score, N_score) > maximal_min_score_on_acc_path[N]
                if cond1 or cond2:
                    queue.append([N, path_length + 1, min(min_score, N_score)])
                    maximal_min_score_on_acc_path[N] = min(min_score, N_score)
                    if cond1:
                        visited[N] = True
                        shortest_path_length[N] = path_length + 1

        return list(shortest_path_length.values())
