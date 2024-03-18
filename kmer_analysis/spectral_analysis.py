from sklearn.manifold import SpectralEmbedding
from sklearn.cluster import SpectralClustering
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import seaborn as sns
import pickle as pkl


class Spectral:
    def __init__(self):
        self.kl_normalized = None
        self.spectral_embedding = None
        self.df = pd.read_pickle('df2.pkl')

        # Convert the DataFrame to a numpy array
        self.array_data = self.df.values.astype(np.float64)
        print(f'Data Shape: {self.array_data.shape}')

    def compute_kl(self, mat):
        n = mat.shape[0]
        kl_scores = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):  # Start from i+1 to avoid redundant calculations
                kl_sum = 0
                print(f'i: {i}, j: {j}')
                for k in range(mat.shape[1]):
                    kl_sum += self.KL_divergence(mat[i, k], mat[j, k])
                kl_scores[i, j] = kl_sum
                kl_scores[j, i] = kl_sum  # Fill the symmetric entry

        # Normalize the KL divergence scores, if necessary
        # For example, using Min-Max normalization
        kl_scores = np.exp(-kl_scores)

        kl_min, kl_max = kl_scores.min(), kl_scores.max()
        self.kl_normalized = (kl_scores - kl_min + 1) / (kl_max - kl_min + 1)
        print(f'KL_shape: {self.kl_normalized.shape}')

        with open(os.path.join('sim.pkl'), 'wb') as f:
            pkl.dump(self.kl_normalized, f)
            print('dumped file')
        sns.heatmap(self.kl_normalized)

    def spectral(self, n_components, n_clusters, mat):
        # Perform spectral embedding on the KL divergence matrix
        self.embedding = SpectralEmbedding(n_components=n_components, affinity='precomputed')
        self.spectral_embedding = self.embedding.fit_transform(mat)
        self.clustering = SpectralClustering(n_clusters=n_clusters, affinity='precomputed')
        self.labels = self.clustering.fit_predict(mat)

        print(self.spectral_embedding.shape[0])

        # return spectral_embedding, clustering

    def KL_divergence(self, p, q):
        p = np.where(p == 0, 1e-9, p)  # Avoid division by zero
        q = np.where(q == 0, 1e-9, q)
        return np.sum(p * np.log2(p / q))

    def plot_embeddings(self, vec_tuple):
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.title("Spectral Embedding")
        plt.scatter(self.spectral_embedding[:, vec_tuple[0]],
                    self.spectral_embedding[:, vec_tuple[1]], c=self.labels, cmap='viridis')
        plt.colorbar(label='Cluster')
        plt.xlabel(f'Embedding Dimension {vec_tuple[0]}')
        plt.ylabel(f'Embedding Dimension {vec_tuple[1]}')

        plt.subplot(1, 2, 2)
        plt.title("Original Data In Feature Space")
        plt.scatter(self.array_data[:, vec_tuple[0]],
                    self.array_data[:, vec_tuple[1]], c=labels,
                    cmap='viridis')
        plt.colorbar(label='Cluster')
        plt.xlabel(f'Feature {vec_tuple[0]}')
        plt.ylabel(f'Feature {vec_tuple[1]}')

        plt.tight_layout()

        # plt.savefig(f'emb_1.png')
        plt.savefig(f'emb_{vec_tuple[0]}_{vec_tuple[1]}.png')

        plt.show()


if __name__ == '__main__':

    Spec = Spectral()
    # Spec.compute_kl(Spec.array_data)
    with open('sim.pkl', 'rb') as f:
        data = np.load(f, allow_pickle=True)

    from scipy.sparse.csgraph import laplacian
    from scipy.linalg import eigh

    laplacian_mat = laplacian(data, normed=True)
    eigval, eigvec = eigh(laplacian_mat)

    zero_thres = 1e-5
    zeo_eigs_count = np.sum(eigval < zero_thres)
    print(f'Num zero: {zeo_eigs_count}')
    Spec.spectral(4, 4, data)

    print(f"num vecs: {Spec.spectral_embedding.shape}")
    labels = Spec.labels
    df = Spec.df

    Spec.plot_embeddings((0, 1))
    Spec.plot_embeddings((1, 2))
    Spec.plot_embeddings((2, 3))

    index_to_cluster = {index: label for index, label in enumerate(labels)}
    index_to_contig = {index: df.index[index] for index in range(len(df))}
    indices_to_find = [69, 44, 195, 176, 2, 23, 14, 198, 199, 15, 12, 10, 11, 8, 161, 4, 18, 20]

    contig_names = [f'contig_{i}' for i in indices_to_find]

    # Step 2: Find the actual row indices in df2 for these contig names
    actual_indices = []
    missing_contigs = []
    for contig_name in contig_names:
        try:
            index = df.index.get_loc(contig_name)
            actual_indices.append(index)
        except KeyError:
            missing_contigs.append(contig_name)

    # Check for any missing contigs
    if missing_contigs:
        print("Missing contigs in df2:", missing_contigs)
    else:
        print("All contigs found in df2.")

    contig_to_cluster = {index_to_contig[index]: index_to_cluster[index]
                         for index in actual_indices}

    print(contig_to_cluster)

# {'contig_160': 2, 'contig_194': 0, 'contig_3': 0, 'contig_68': 0, 'contig_197': 2, 'contig_198': 1, 'contig_7': 1, 'contig_9': 0, 'contig_10': 1, 'contig_11': 1, 'contig_43': 0, 'contig_13': 0, 'contig_14': 0, 'contig_175': 2, 'contig_17': 1, 'contig_19': 0, 'contig_22': 0}