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
    # plt.title('assembly_2.fasta')
    # plt.xlabel('Number of Clusters (k)')
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
        print(cluster_set)

    # cluster_dicts = []
    # for _, cluster in clusters:
    #     names = cluster.index.str.split('_').str[0]
    #     count = Counter(names)
    #     # total = sum(count.values())
    #     proportions = {key: value for key, value in count.items()}
    #     cluster_dicts.append(proportions)
    # indices = np.arange(0, len(cluster_dicts))
    # all_species = set()
    # for cluster_dict in cluster_dicts:
    #     all_species.update(cluster_dict.keys())
    
    # fig, ax = plt.subplots()

    # bar_width = 0.1

    # for i, species in enumerate(all_species):
    #     # Extracting accuracy values for each cluster for the current bacterium
    #     accuracies = [cluster_dict.get(species, 0) for cluster_dict in cluster_dicts]
        
    #     # Plotting
    #     ax.bar(indices + i * bar_width, accuracies, width=bar_width, label=species)

    # # Setting the x-axis ticks and labels
    # ax.set_xticks(indices + (len(all_species) * bar_width / 2))
    # ax.set_xticklabels([f'Cluster {i+1}' for i in range(len(cluster_dicts))])

    # # Adding labels and title
    # ax.set_ylabel('Composition')
    # ax.set_title(f'K={k}, assembly_2.fasta')

    # # Adding a legend
    # ax.legend()

    # # Display the plot
    # plt.show()

def main():
    temp_dict, kmer_dict = make_temp_dict(5)
    myReader = FastAreader('assembly_2.fasta')
    objList = []
    for head, seq in myReader.readFasta():
        objList.append(Contig(head, seq, templateDict=temp_dict, kmerMap= kmer_dict))
    profile = BuildProfile(objList=objList, templateKmerDict=temp_dict)
    with open(os.path.join('assembly_2_profiles.pkl'), 'wb') as file:
        pkl.dump(profile.df, file, protocol=pkl.HIGHEST_PROTOCOL)

    # filenames = ['Mycobacteriumphage.fna', 'VibrioPhage_pVa-1.fna', 'Staphylococcus_aureus.fna', 'Escherichia_coli_complete.fna']
    # for file in filenames:
    #     for i in range(1, 11):
    #         get_original_profile(file, k = i)
    # distance = pdist(profile.df) #default is euclidean distance
    # linkage_matrix = linkage(distance, method='ward') # can be centroid or median for Euclidean distance 
    # The next step is to perform AgglomerativeClustering and then see how hierarchial clustering performs 
    # also play around with dendrograms at small sample sizes

    # also you need to actually work on bernicks samples 
    # edges = {69, 44, 195, 176, 23, 14, 198, 199, 15, 12, 10, 11, 8, 161, 4, 18, 20}
    # kmeans_cluster(5, edges)



if __name__ == "__main__":
    main()