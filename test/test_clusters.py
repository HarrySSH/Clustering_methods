import pytest
import glob
import pandas as pd
import numpy as np
import os 
import sys
sys.path.append("..")

from clusters.algs import Ligand, HierarchicalClustering, PartitionClustering,Silhouette_Coefficient, ClusterNode, getEuclidean, Rand_Index
import matplotlib.pyplot as plt
import time
import numpy as np
import pandas as pd
from sklearn.datasets import fetch_openml
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.datasets.samples_generator import make_blobs
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import json

X, y = make_blobs(n_samples=300, centers=4, cluster_std=0.20, random_state=0)

onBites = []
for i in range(300):
    
    onBites.append(str(list(X[:,0])[i]) +',' + str(list(X[:,0])[i]))
    
# Build a stimulated ligands dataset
df = pd.DataFrame (onBites ,columns = ['OnBits'])
df['LigandID'] = list(range(0,300))
df['Score']='.'
df['SMILES']='.'

class pseudo_Ligand:
	# Initiliza the class
	def __init__(self, df):
		self.ID = df['LigandID'].values[0]
		self.Score = df['Score'].values[0]
		self.SMILES = df['SMILES'].values[0]
		OnBits = np.zeros(1024)
		#print(df['OnBits'].values[0])
		OnBits = [ float(x) for x in df['OnBits'].values[0].split(',')]
		self.OnBits = OnBits
	# return the ID
	def get_ID(self):  # return the ID
		return self.ID   
	def get_scores(self): # return the scores
		return self.Score
	def get_SMILES(self): # return the SMILES
		return self.SMILES
	def get_OnBits(self): # return the onBits information
		return self.OnBits

Data = list(range(df.shape[0]))
Ligands = []

for i in Data:

    a_ligand = pseudo_Ligand(df[df.index==i])
    Ligands.append(a_ligand)


def test_partitioning():
    kmeans_C = PartitionClustering(Ligands, 4,50)
    kmeans_C.implement()

    silhouette_coefficient = Silhouette_Coefficient(Ligands, kmeans_C.labels)
    rand_index = Rand_Index(kmeans_C.labels, y)
    label = [i for i in kmeans_C.labels]
    cols = ['r','b','g','y']
    for i  in range(300):
        plt.plot(X[:,0][i], X[:,1][i], cols[label[i]]+'.')
    assert silhouette_coefficient > 0.75, 'Should give a good cluster'
    assert rand_index > 0.75, 'Should also be similar with the groundtruth'

def test_hierarchical():
    Hirar_C = HierarchicalClustering(Ligands, 4)
    Hirar_C.implement()
    res = Hirar_C.implement_viz()
    Hirar_C.show(res,300)
    silhouette_coefficient = Silhouette_Coefficient(Ligands, Hirar_C.labels)
    rand_index = Rand_Index(Hirar_C.labels, y)
    label = [i for i in Hirar_C.labels]
    cols = ['r','b','g','y']
    for i  in range(300):
        plt.plot(X[:,0][i], X[:,1][i], cols[label[i]]+'.')
    assert silhouette_coefficient > 0.75, 'Should give a good cluster'
    assert rand_index > 0.75, 'Should also be similar with the groundtruth'    
    
    
