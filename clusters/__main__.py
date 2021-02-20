import glob
import pandas as pd
import numpy as np
import os 
import sys
from algs import Ligand, HierarchicalClustering, PartitionClustering,Silhouette_Coefficient, ClusterNode, getEuclidean, Rand_Index
import matplotlib.pyplot as plt
from __future__ import print_function
import time
import numpy as np
import pandas as pd
from sklearn.datasets import fetch_openml
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
%matplotlib inline
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import json

csv = pd.read_csv('../ligand_information.csv', index_col= False)
Data = list(np.random.randint(1,csv.shape[0],size=2000))

Ligands = []

for i in Data:

    a_ligand = Ligand(csv[csv.index==i])
    Ligands.append(a_ligand)
    
IDs = []
scores = []
dataset = []
for ligand in Ligands:
	IDs.append(ligand.get_ID())
	scores.append(ligand.get_scores())
	dataset.append(ligand.get_OnBits())
dataset = np.stack((dataset), axis=0)
onBits = [ 'pixel'+str(i) for i in range(dataset.shape[1]) ]
df = pd.DataFrame(dataset,columns=onBits)
_df = df[onBits].values


time_start = time.time()
tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
tsne_results = tsne.fit_transform(_df)
print('t-SNE done! Time elapsed: {} seconds'.format(time.time()-time_start))

df['scores'] = scores
#score_ir = pd.interval_range(start =-6.2, end = 534.3, freq = 100)

df['bins_scores'] = pd.cut(df['scores'], bins = [-7,37,81,450])
del df['scores']
df['tsne-2d-one'] = tsne_results[:,0]
df['tsne-2d-two'] = tsne_results[:,1]
plt.figure(figsize=(16,10))
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-two",
    hue="bins_scores",
    palette=sns.color_palette("hls", len(set(df['bins_scores']))),
    data=df,
    legend="full",
    alpha=0.3
)

# Trying to find the most suitable cluster numbers
num_Cs = []
Silhouette_Coefficients =å []

for num in range(2,10):
    print(num)
    kmeans = PartitionClustering(Ligands, num, 10)
    kmeans.implement()
    Sil_co = Silhouette_Coefficient(Ligands, kmeans.labels)
    num_Cs.append(num)
    Silhouette_Coefficients.append(Sil_co)


plt.plot(num_Cs, Silhouette_Coefficients,'.r')

kmeans = PartitionClustering(Ligands, 5, 10)
kmeans.implement()

df['bins_scores'] = kmeans.labels
df['tsne-2d-one'] = tsne_results[:,0]
df['tsne-2d-two'] = tsne_results[:,1]
plt.figure(figsize=(16,10))
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-two",
    hue="bins_scores",
    palette=sns.color_palette("hls", len(set(df['bins_scores']))),
    data=df,
    legend="full",
    alpha=0.3
)

Hiera_C = HierarchicalClustering(Ligands, 10)
Hiera_C.implement()

df['bins_scores'] = Hiera_C.labels
df['tsne-2d-one'] = tsne_results[:,0]
df['tsne-2d-two'] = tsne_results[:,1]
plt.figure(figsize=(16,10))
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-two",
    hue="bins_scores",
    palette=sns.color_palette("hls", len(set(df['bins_scores']))),
    data=df,
    legend="full",
    alpha=0.3
)

Silhouette_Coefficient(Ligands, Hiera_C.labels)

Silhouette_Coefficient(Ligands, kmeans.labels)

Rand_Index(kmeans.labels, Hiera_C.labels)

# Kmeans
index = 2000
Cluster_scores = {}
for i in range(index):
    c = kmeans.labels[i]
    if c not in Cluster_scores: 
        Cluster_scores[c] = [kmeans.scores[i]]
    else:
        Cluster_scores[c].append(kmeans.scores[i])
num = 0
keys = []
for index, key in enumerate(Cluster_scores):
    if len(Cluster_scores[key]) == 1:
        pass
    else:
        keys.append(key)
        num = num +1


plt.rcParams["figure.figsize"]=30,20

fig, axs = plt.subplots(nrows = 2, ncols = 3)
p = 0
for ax in axs.flatten():
    if p <5:
        values = Cluster_scores[keys[p]]
        ax.hist(values, 100, density = True)
        ax.set_title('Cluster '+str(keys[p]))
        p = p +1
    else:
        df['bins_scores'] = kmeans.labels
        cols = ['r','b','g','y','m']
        for i  in range(2000):
            
            ax.plot(df["tsne-2d-one"][i], df["tsne-2d-two"][i], cols[kmeans.labels[i]]+'.')
        ax.set_title('TSNE')
        ax.set_xlabel('tsne-2d-one')
        ax.set_ylabel('tsne-2d-two')

    
fig.tight_layout()
plt.show()

num = 0
for index, key in enumerate(Cluster_scores):
    if len(Cluster_scores[key]) == 1:
        continue
    else:
        num = num +1


top_scores = []
for i in range(p):
    values = Cluster_scores[keys[i]]
    score = max(values)
    top_scores.append(score)
    top_scores.sort()


Leads = []
for Ligand in Ligands:
    if Ligand.get_scores() in top_scores:
        Leads.append(Ligand.get_ID())  
