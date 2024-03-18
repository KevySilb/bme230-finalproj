from fastaReader import FastAreader
from copy import deepcopy
import itertools as it
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os
import pickle as pkl
import random 
# from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cophenet
# from scipy.spatial.distance import pdist
import seaborn as sb
import sklearn
from sklearn.cluster import AgglomerativeClustering, KMeans
# from sklearn.preprocessing import StandardScaler
# import sklearn.metrics as sm
from collections import Counter
from sklearn.metrics import silhouette_score

class BuildProfile:

    def __init__(self, objList:list, templateKmerDict) -> None:
        rows_list = []
        self.col_index = [key for key in templateKmerDict.keys()]
        for obj in objList:
            row_dict = {key: obj.kmer_dict[key] for key in self.col_index}
            rows_list.append(row_dict)
        self.df = pd.DataFrame(rows_list, index = [obj.head for obj in objList], columns = self.col_index)
        self.length = len(obj.seq)
        self.k = obj.k

class Contig:
      
    def __init__(self, head, seq, templateDict, kmerMap) -> None:
        self.head = head
        self.seq = seq
        self.kmer_dict = deepcopy(templateDict)
        self.k = len(next(iter(kmerMap)))
        # count the occurance of each 4-mer in the sequence
        for i in range(len(self.seq) - self.k + 1):
            kmer = self.seq[i:i + self.k]
            if kmer not in kmerMap:
                continue
            key = tuple(sorted([kmer, kmerMap[kmer]]))
            if key in self.kmer_dict:
                self.kmer_dict[key] += 1.0
        # sum all the counts and normalize the counts into probabilities
        value_sum = np.longdouble(0.0)
        for value in self.kmer_dict.values():
            if value > 0:
                value_sum += value
        for key in self.kmer_dict.keys():
            self.kmer_dict[key] /= value_sum

def make_temp_dict(k):
    alphabet = ['A', 'T', 'G', 'C']
    nucMap = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    kmerPermutations = [''.join(p) for p in it.product(alphabet, repeat=k)]
    templateDict = dict()
    kmerMap = dict()
    # make kmerMap
    for kmer in kmerPermutations:
        revKmer = rev_compliment(kmer, nucMap)
        kmerMap[kmer] = revKmer
        key = tuple(sorted([kmer, revKmer]))
        if key not in templateDict:
            templateDict[key] = np.longdouble(k*2)
    return templateDict, kmerMap

def rev_compliment(kmer, nucMap):
    string = []
    for char in kmer:
        string.append(nucMap[char])
    return ''.join(string[::-1])
# rewrite this function to generate fasta files where the contig size varies, but the place in the sequence is random 
def make_subsequences(inFile, dirname):
    myReader = FastAreader(inFile)
    sequence = ""
    seq_length = 0
    for head, seq in myReader.readFasta():
        sequence = seq
        seq_length = len(seq)
        break

    def func(x):
        return 3*x
    
    basepairlengths = []
    val = 400
    for i in range(7):
        basepairlengths.append(val)
        val = func(val)

    if not os.path.exists(dirname):
        os.makedirs(dirname)
    
    for bplen in basepairlengths:
        with open(os.path.join(dirname, f'subsequences_{bplen}bp.fa'), 'w') as f:
            for i in range(0, seq_length - bplen + 1, bplen):
                subsequence = sequence[i:i+bplen]
                header = f'>CONTIG{i}:{i+bplen}:{bplen}\n'
                f.write(header)
                f.write(subsequence + '\n')

def get_profiles(dir_path = 'Escherichia_coli_subsequences'):
    templateKmerDict, kmerMap = make_temp_dict(4)
    profiles = dict()
    for filename in os.listdir(dir_path):
        contigList = []
        profile = None
        if filename.endswith('.fa'):
            file_path = os.path.join(dir_path, filename)
            myReader = FastAreader(file_path)
            for head, seq in myReader.readFasta():
                contig_obj = Contig(head, seq, templateKmerDict, kmerMap)
                contigList.append(contig_obj)
            profile = BuildProfile(contigList, templateKmerDict)
            profiles[filename] = profile
    with open(os.path.join(dir_path, 'profiles.pkl'), 'wb') as file:
        pkl.dump(profiles, file, protocol=pkl.HIGHEST_PROTOCOL)
        
def get_original_profile(inFile, k):
    outFile = f"{inFile.split('.fna')[0]}_{k}.pkl"
    templateKmerDict, kmerMap = make_temp_dict(k)
    myReader = FastAreader(inFile)
    for head, seq in myReader.readFasta():
        profile = BuildProfile([Contig(head, seq, templateKmerDict, kmerMap)], templateKmerDict)
    
    with open(outFile, 'wb') as file:
        pkl.dump(profile, file, protocol=pkl.HIGHEST_PROTOCOL)

def get_random_sample():
    filenames = ['Staphylococcus_aureus.fna', 'Pyrococcus_furiosus.fna', 'Escherichia_coli_complete.fna', 'Zika_virus_complete.fna', 'VibrioPhage_pVa-1.fna', 'Mycoplasma_canadense.fna', 'WestNile_complete.fna']
    readers = [FastAreader(file) for file in filenames]
    N = 20000
    with open('random_sample_7_species.fasta', 'w') as f:
        for _ in range(N):
            reader = np.random.choice(readers)
            contig_len = np.random.randint(2000, 6000)
            for head, seq in reader.readFasta():
                i_max = len(seq)
                name = head.split(' ')[1]
                i = random.randint(0, i_max - contig_len)
                sub_seq = seq[i:i + contig_len]
                f.write(f'>{name}_{i}:{i+contig_len}\n')
                f.write(f'{sub_seq}\n')

def kmeans_cluster(k, edges):
    templateKmerDict, kmerMap = make_temp_dict(k)
    objList = []
    myReader = FastAreader('assembly_2.fasta')
    for head, seq in myReader.readFasta():
        obj = Contig(head, seq, templateKmerDict, kmerMap)
        objList.append(obj)
    profile = BuildProfile(objList, templateKmerDict)
    X = profile.df.values
    # ---------------------------------------------------------------------------
    # wcss = []
    # k_values = range(2, 11)  # You can adjust the range of k values as needed
    # for j in k_values:
    #     kmeans = KMeans(n_clusters=j, random_state=42)
    #     kmeans.fit(X)
    #     wcss.append(kmeans.inertia_)  # .inertia_ gives the WCSS value for the fitted model

    # # Plotting the elbow plot
    # plt.figure(figsize=(8, 6))
    # plt.plot(k_values, wcss, marker='o')
    # plt.title('Elbow Plot for assembly.fasta')
    # plt.xlabel('Number of Clusters')
    # plt.ylabel('WCSS')
    # plt.xticks(k_values)
    # plt.grid(True)
    # plt.show()
    # -----------------------------------------------------------------------------
    # k_values = range(2, 11)
    # silhouette_scores = []
    # for k in k_values:
    #     kmeans = KMeans(n_clusters=k, random_state=42)
    #     cluster_labels = kmeans.fit_predict(X)
    #     silhouette_avg = silhouette_score(X, cluster_labels)
    #     silhouette_scores.append(silhouette_avg)
    # plt.figure(figsize=(8, 6))
    # plt.plot(k_values, silhouette_scores, marker='o')
    # plt.title('Silhouette scores for assembly.fasta')
    # plt.xlabel('Number of Clusters')
    # plt.ylabel('Silhouette Score')
    # plt.xticks(k_values)
    # plt.grid(True)
    # plt.show()
    # -----------------------------------------------------------------------------
    kmeans = KMeans(n_clusters=7)
    kmeans.fit(X)
    profile.df['cluster'] = kmeans.labels_
    clusters = profile.df.groupby('cluster')
    for number, cluster in clusters:
        cluster_set = set()
        for contig in cluster.index:
            cluster_set.add(int(contig.split('_')[1]))
        print(number, cluster_set & edges)

def main():
    edges = {2, 69, 44, 195, 176, 23, 14, 198, 199, 15, 12, 10, 11, 8, 161, 4, 18, 20}
    # myReader = FastAreader('assembly_2.fasta')
    # with open('majority_bandage.fasta', 'w') as file:
    #     for head, seq in myReader.readFasta():
    #         edge = int(head.split('_')[1])
    #         if edge not in edges:
    #             continue
    #         else:
    #             file.write(f'>{head}\n{seq}\n')
    # myReader = FastAreader('majority_bandage.fasta')
    # edges2 = set()
    # for head, seq in myReader.readFasta():
    #     edge = int(head.split('_')[1])
    #     edges2.add(edge)
    # print(edges2 ^ edges)

    kmeans_cluster(5, edges)



        

if __name__ == "__main__":
    main()