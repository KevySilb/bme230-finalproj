import gfapy
import numpy as np
import pickle


class GFA:
    def __init__(self, file_path):
        self.file_path = file_path
        self.gfa = gfapy.Gfa.from_file(file_path)
        self.coverage_data = self.extract_coverage_data()
        self.adjacency_matrix = self.create_adj_matrix()
        print(self.coverage_data)

    def extract_coverage_data(self):
        coverage_data = {}
        for segment in self.gfa.segments:
            segment_name = segment.name

            if 'dp' in segment.tagnames:
                depth = segment.dp
            else:
                depth = 0  # Default to 0 if the tag doesn't exist

            coverage_data[segment_name] = depth

        return coverage_data

    def create_adj_matrix(self):
        # Create a list of segment names to index them
        segments = [segment.name for segment in self.gfa.segments]
        self.segment_indices = {name: idx for idx, name in enumerate(segments)}
        self.index_to_segment = {idx: name for name, idx in self.segment_indices.items()}  # Add this line

        N = len(segments)
        adjacency_matrix = np.zeros((N, N))

        for line in self.gfa.lines:
            if line.record_type == 'L':
                from_idx = self.segment_indices[line.from_segment.name]
                to_idx = self.segment_indices[line.to_segment.name]

                # Assume the adjacency matrix is for an undirected graph
                # and set the value to the average coverage of the two segments
                from_coverage = self.coverage_data.get(line.from_segment.name, 0)
                to_coverage = self.coverage_data.get(line.to_segment.name, 0)

                avg_coverage = (from_coverage + to_coverage) / 2
                adjacency_matrix[from_idx, to_idx] = avg_coverage
                adjacency_matrix[to_idx, from_idx] = avg_coverage  # for undirected graph

        return adjacency_matrix

    def dump_adj(self, file_name='adj.pkl'):
        with open(file_name, 'wb') as f:
            pickle.dump(self.adjacency_matrix, f)
            print(f'Dumped File')

    def get_segment_name(self, index):
        return self.index_to_segment[index]

    def get_index(self, contig_name):
        return self.segment_indices[contig_name]
    def dump_coverage_data(self, file_name='cov_data.pkl'):
        with open(file_name, 'wb') as f:
            pickle.dump(self.coverage_data, f)
            print('Dumped coverage')


if __name__ == "__main__":
    gfa_file_path = 'assembly_2graph_combined.gfa'
    gfa = GFA(gfa_file_path)
    gfa.create_adj_matrix()
    print(gfa.adjacency_matrix)

    gfa.dump_adj()
    x = np.max(gfa.adjacency_matrix)
    gfa.dump_coverage_data()

    print(x)
