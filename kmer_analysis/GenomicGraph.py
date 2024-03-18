import networkx as nx
import numpy as np
import pickle as pkl
from GfaReader import GFA
import heapq


class GenomicGraph:
    def __init__(self, adj_path, cov_path, gfa):
        """
        graph: NetworkX graph object representing the genomic assembly
        coverage_data: Dictionary mapping each node to its coverage value
        gfa: GFA object from GfaReader
        """
        with open(adj_path, 'rb') as f:
            adj_matrix = pkl.load(f)

        with open(cov_path, 'rb') as f:
            self.coverage_data = pkl.load(f)
        self.graph = nx.from_numpy_array(adj_matrix, parallel_edges=False, create_using=nx.Graph)
        self.gfa = gfa

    def find_largest_connected_component(self):
        # Find the connected component in the graph with the maximum number of nodes
        largest_component = max(nx.connected_components(self.graph), key=len)
        return self.graph.subgraph(largest_component)

    def find_significant_cycle(self):
        def find_cycle(start_node):
            visited = set()
            stack = [(start_node, [start_node])]

            while stack:
                node, path = stack.pop()
                visited.add(node)

                for neighbor in self.graph.neighbors(node):
                    if neighbor == start_node and len(path) > 2:
                        return path
                    elif neighbor not in visited:
                        stack.append((neighbor, path + [neighbor]))

            return []

        # Find nodes that are connectors (nodes that connect two or more distinct paths)
        candidate_nodes = [node for node in self.graph.nodes if len(list(self.graph.neighbors(node))) > 1]

        # Explore cycles starting from candidate nodes
        largest_cycle = []
        for node in candidate_nodes:
            cycle = find_cycle(node)
            if len(cycle) > len(largest_cycle):
                largest_cycle = cycle

        # Convert node indices in the cycle back to segment names
        largest_cycle_contigs = [self.gfa.get_segment_name(node) for node in largest_cycle]

        return largest_cycle_contigs

    def find_highest_coverage_node(self, subgraph):
        highest_coverage_node = max(subgraph, key=lambda node: self.coverage_data.get(node, 0))
        return highest_coverage_node

    def a_star_search(self, start_node, largest_cluster_nodes):
        def calculate_average_coverage(largest_cluster_nodes):
            total_coverage = sum(
                self.coverage_data.get(self.gfa.get_segment_name(node), 0) for node in largest_cluster_nodes)
            return total_coverage / len(largest_cluster_nodes)

        def heuristic(node):
            contig_name = self.gfa.get_segment_name(node)
            return -self.coverage_data.get(contig_name, 0)

        def heuristic_version_2(node, path, largest_cluster_nodes):
            coverage = self.coverage_data.get(self.gfa.get_segment_name(node), 0)
            unvisited_nodes = len(largest_cluster_nodes) - len(set(path))
            return -(coverage + unvisited_nodes)

        def balanced_heuristic(node, path, average_coverage):
            node_coverage = self.coverage_data.get(self.gfa.get_segment_name(node), 0)
            coverage_deviation = abs(node_coverage - average_coverage)
            path_length_factor = len(path)
            return coverage_deviation + path_length_factor  # Seek to minimize deviation while considering path length

        # open_heap = [(heuristic(start_node), 0, start_node, [start_node])]
        open_heap = [(heuristic_version_2(start_node, [start_node], largest_cluster_nodes), 0, start_node, [start_node])]

        visited = set([start_node])

        while open_heap:
            _, g_score_current, current, path = heapq.heappop(open_heap)
            print(f"path: {[self.gfa.get_segment_name(index) for index in path]}")

            print(f"path: {path}")

            if current == start_node and len(path) > 1 and path[-2] != start_node:
                return path

            neighbors = list(self.graph.neighbors(current))
            for neighbor in neighbors:
                # Skip any nodes that would cause an immediate loop
                if neighbor == start_node or neighbor in path:
                    continue

                new_g_score = g_score_current + self.coverage_data.get(self.gfa.get_segment_name(neighbor), 0)
                new_path = path + [neighbor]
                f_score = new_g_score + heuristic_version_2(neighbor, new_path, largest_cluster_nodes)
                # f_score = new_g_score + heuristic(neighbor, new_path, largest_cluster_nodes)

                # If this neighbor leads to a path we've not visited, or it's part of the cluster and not visited yet, add to the heap
                if neighbor in largest_cluster_nodes and tuple(new_path) not in visited:
                    heapq.heappush(open_heap, (f_score, new_g_score, neighbor, new_path))
                    visited.add(tuple(new_path))

            # If we are at a leaf node (no unvisited neighbors), backtrack
            if not any(neighbor not in visited for neighbor in neighbors):
                for back_idx in range(len(path) - 2, -1, -1):
                    back_node = path[back_idx]
                    back_neighbors = list(self.graph.neighbors(back_node))
                    # Check if there are any unvisited branches from this point
                    unvisited_branches = [n for n in back_neighbors if tuple(path[:back_idx + 1] + [n]) not in visited]
                    if unvisited_branches:
                        # Backtrack to this branching point and explore the unvisited branches
                        for branch in unvisited_branches:
                            new_g_score = sum(
                                self.coverage_data.get(self.gfa.get_segment_name(p), 0) for p in path[:back_idx + 1])
                            new_path = path[:back_idx + 1] + [branch]
                            # f_score = new_g_score + heuristic(branch)
                            f_score = new_g_score + heuristic_version_2(branch, new_path, largest_cluster_nodes)
                            heapq.heappush(open_heap, (f_score, new_g_score, branch, new_path))
                            visited.add(tuple(new_path))
                        break  # Only backtrack to the first valid branching point

        return None  # No circular path found


    def reconstruct_path(self, came_from, current):
        path = [current]
        while current in came_from:
            print(f'recont path: {current}')
            current = came_from[current]
            path.append(current)
        path.reverse()
        return [self.gfa.get_segment_name(node) for node in path]


def main():
    gfa = GFA('assembly_2graph_combined.gfa')
    # print(f'edges size: {len(edges)}')
    GG = GenomicGraph('adj.pkl', 'cov_data.pkl', gfa)

    sc = GG.find_significant_cycle()
    print(sc)
    largest_component = GG.find_largest_connected_component()
    largest_contig_comp = [gfa.get_segment_name(index) for index in largest_component]

    print("Nodes in the largest connected component:", largest_component.nodes)
    print(len(largest_contig_comp))
    print(f'Contigs in the largest component: {largest_contig_comp}')
    largest_cluster = GG.find_largest_connected_component()

    start_contig = GG.find_highest_coverage_node(largest_contig_comp)
    print(f'start contig: {start_contig}')
    start_node = gfa.get_index(start_contig)
    print(f'start node: {start_node}')

    print(f"coverages: {GG.coverage_data}")
    start_score = GG.coverage_data[start_contig]
    print(f"start sc: {start_score}")

    x = [9, 1, 17, 7, 198, 197, 13, 22, 43, 68]
    contig_path = [gfa.get_segment_name(index) for index in x]
    print(contig_path)

    circular_path = GG.a_star_search(start_node, largest_cluster)
    print(f'circlular path: {circular_path}')
    contig_path = [gfa.get_segment_name(index) for index in circular_path]
    print(f'contig path: {contig_path}')

if __name__ == '__main__':
    main()
