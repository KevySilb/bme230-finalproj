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

def main():
    cmd = commandLine()
    with open(cmd.args.inFile, "rb") as file:
        df = pickle.load(file)
    t = 0.4
    clt = Cluster(df)
    

if __name__ == "__main__":
    main()
    