import os
import pickle as pkl
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from profile_builder import Contig, BuildProfile
from scipy.stats import norm
import seaborn as sns

def get_profiles(dir_path):
    for filename in os.listdir(dir_path):
        if filename.endswith('.pkl'):
            file_path = os.path.join(dir_path, filename)
            with open(file_path, 'rb') as file:
                profiles = pkl.load(file)
    return profiles

def KL_divergence(p_i, q_i):
    return p_i * np.log2(p_i / q_i)

def plot_normal_dist(dir_path, inFile):
    def extract_contig_length(filename):
        return int(filename.split('_')[1].split('bp.fa')[0])
    profiles = get_profiles(dir_path)
    with open(inFile, 'rb') as file:
        original_genome = pkl.load(file)
    p = original_genome.values[0]
    for filename, obj in profiles.items():
        KL_matrix = np.zeros(obj.df.values.shape, dtype = np.longdouble)
        for i, q in enumerate(obj.df.values):
            assert len(q) == len(p)
            for j in range(len(p)):
                KL_matrix[i, j] = KL_divergence(p[j], q[j])
        data = np.sum(KL_matrix, axis = 1)
        data = data*136
        mu, std = norm.fit(data)
        xmin, xmax = min(data), max(data)
        x = np.linspace(xmin - 5, xmax + 5, 100)
        normal_dist = norm.pdf(x, mu, std)
        plt.plot(x, normal_dist, label = f'{extract_contig_length(filename)}')
        plt.scatter(mu, 0, label = f'{mu:.2g} bits')

    plt.xlabel('Kullback-Leibler divergence (bits)')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.show()

def contig_length_KL_div(dir_name, inFile):
    def extract_contig_length(filename):
        return int(filename.split('_')[1].split('bp.fa')[0])
    profiles = get_profiles(dir_name)
    with open(inFile, 'rb') as file:
        original_genome = pkl.load(file)
    p = original_genome.values[0]
    x = []
    y = []
    for filename, obj in profiles.items():
        KL_matrix = np.zeros(obj.df.values.shape, dtype = np.longdouble)
        for i, q in enumerate(obj.df.values):
            assert len(q) == len(p)
            for j in range(len(p)):
                KL_matrix[i, j] = KL_divergence(p[j], q[j])
        data = np.sum(KL_matrix, axis = 1)
        data = data * len(p)
        mu, _ = norm.fit(data)
        x.append(extract_contig_length(filename))
        y.append(mu)
    
    # now add the comparison between original genome and itself as a control. 
    KL_matrix = np.zeros_like(p)
    for i in range(len(p)):
        KL_matrix[i] = KL_divergence(p[i], p[i])
    data = np.sum(KL_matrix)
    data = data * len(p)
    mu, _ = norm.fit(data)
    plt.scatter(5498578, mu, c = 'r', label = 'control')

    combined = sorted(zip(x,y), key = lambda x: x[0])
    x, y = zip(*combined)
    plt.scatter(x, y, label = r'mean of KLD')
    plt.xlabel('contig length')
    plt.ylabel('bits')
    plt.xscale('log')
    plt.title('E.coli vs subsequences of E.coli')
    plt.legend()
    plt.show()

def compare_KLs(filenames:list):
    dfs_dict = dict()
    for filename in filenames:
        with open(filename, 'rb') as file:
            profile = pkl.load(file)
            dfs_dict[profile.length] = profile.df
    
    N = len(dfs_dict)

    lengths = [key for key in dfs_dict.keys()]
    KL_matrix = np.zeros((N, N), dtype = np.longdouble)
    row_index = []
    for i in range(KL_matrix.shape[0]):
        row = dfs_dict[lengths[i]]
        row_index.append(row.index[0].split(' ')[1] + f' {lengths[i]}')
        for j in range(KL_matrix.shape[1]):
            col = dfs_dict[lengths[j]]
            KLD = np.longdouble(0.0)
            for k in range(row.shape[1]):
                KLD += KL_divergence(row.iloc[0, k], col.iloc[0, k])
            KL_matrix[i, j] = KLD
    col_index = row_index
    KLD_DF = pd.DataFrame(KL_matrix, index = row_index, columns = col_index)
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(KLD_DF, annot = True, fmt = ".2f", cmap = "viridis", xticklabels = KLD_DF.columns, yticklabels=KLD_DF.index)
    ax.xaxis.tick_top()
    plt.title("")
    plt.show()

def compare_KL_kmer_size(filenames, delete = False):
    # filenames = []
    # for i in range(1, 11):
    #     filenames.append(f'VibrioPhage_pVa-1_{i}.pkl')
    #     filenames.append(f'Mycobacteriumphage_{i}.pkl')
    # compare_KL_kmer_size(filenames, delete = True)
    dfs_dict = dict()
    names = []
    kmersizes = []
    genomesizes = []
    for i in range(len(filenames)):
        filename = filenames[i]
        kmersize = filename.split('_')[-1].split('.pkl')[0]
        name = filename.split('_')[0]
        filetag = name + kmersize
        names.append(name)
        if int(kmersize) not in kmersizes:
            kmersizes.append(int(kmersize))
        with open(filename, 'rb') as file:
            profile = pkl.load(file)
            genomesizes.append(profile.length)
            dfs_dict[filetag] = profile.df

    dfs = [df.values for df in dfs_dict.values()]

    KLD_between_kmersize = []
    for i in range(0, len(dfs_dict) - 1, 2):
        KLD = np.longdouble(0.0)
        df1 = dfs[i][0]
        df2 = dfs[i+1][0]
        assert len(df1) == len(df2)
        for j in range(len(df1)):
            KLD += KL_divergence(df1[j], df2[j])
        KLD_between_kmersize.append(KLD)
    assert len(kmersizes) == len(KLD_between_kmersize)
    plt.scatter(kmersizes, KLD_between_kmersize)
    plt.xlabel('Kmer size')
    plt.ylabel('Kullbach-Leibler Divergence')
    plt.title('P = {} {:.3f}Kb, Q = {} {:.3f}Kb'.format(names[0], genomesizes[0] / 1e3, names[1], genomesizes[1] / 1e3))
    plt.xticks(kmersizes)
    plt.show()
    if delete:
        for filename in filenames:
            try:
                os.remove(filename)
                print(f"Deleted {filename}")
            except FileNotFoundError:
                print(f"{filename} not found, skipping.")

def main():
    pass

if __name__ == "__main__":
    main()