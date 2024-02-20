import numpy as np
import pandas as pd
from commandLine import commandLine
import pickle
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster

class Cluster:

    """
    Arguments:
        df (DataFrame) : a pandas dataframe that is a joint probability matrix.
            Rows == contiguous sequences
            Columns == probability of specific 4 letter word in the nucleotide alphabet appearing in sequence.

        t = clustering threshold 
            distance : Forms flat clusters so that the original observations in each flat cluster have no greater a cophenetic distance than t.
    """
    
    def __init__(self, df, t) -> None:
        self.dist = squareform(pdist(df.values, metric='euclidean'))
        self.linkage = linkage(self.dist, method = 'ward')
        self.clusters = fcluster(self.linkage, t, criterion= 'distance')

"""
        self.cluster_assignments = fcluster(self.linkage_matrix, t, criterion='distance')
        self.cluster_sequences, self.mean_cluster_dists = self._cluster_sequences_and_distances()

    def _cluster_sequences_and_distances(self):
        cluster_sequences = {}
        mean_cluster_dists = {}

        for i, cluster_id in enumerate(self.cluster_assignments):
            if cluster_id not in cluster_sequences:
                cluster_sequences[cluster_id] = []
            cluster_sequences[cluster_id].append(self.df.index[i])

        for cluster_id, sequences in cluster_sequences.items():
            if len(sequences) > 1:
                distances = self.dist_matrix[np.ix_(sequences, sequences)]
                mean_distance = np.sum(distances) / (len(sequences) * (len(sequences) - 1))
                mean_cluster_dists[cluster_id] = mean_distance
            else:
                mean_cluster_dists[cluster_id] = 0  # Single-element clusters have 0 intra-distance

        return cluster_sequences, mean_cluster_dists

    def get_most_related_cluster(self):

        if not self.mean_cluster_dists:
            return None
        min_cluster_id = min(self.mean_cluster_dists, key=self.mean_cluster_dists.get)
        return self.cluster_sequences[min_cluster_id]
"""

def main():
    cmd = commandLine()
    with open(cmd.args.inFile, "rb") as file:
        df = pickle.load(file)
    t = 0.4
    clt = Cluster(df)
    

if __name__ == "__main__":
    main()
    