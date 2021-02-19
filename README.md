# Clustering_BMI_HW2

Implement cluatering algorithms on ligand data


- [Background](#background)
- [Usage](#usage)
- [API](#api)
- [Example](#example)

## Background
In this assignment, I will implement two classical clustering algorithms and then evaluate each algorithm’s performance with a range of parameters. 
For easy implementation, a ligand class was constructed. For evaluating the clustering quality, Silhouette Coefficient is written.
For comparing the clustering pattern, Rand_index matrix is written.

## Usage 
Common alignment methos (global alignment and local alignment) for comparing the similarity of two sequences and computed the similarity scores. 

## API

### Function def getEuclidean(point1, point2):

    '''
    I use the the Euclidean distance. It can be calculated from the Cartesian coordinates of the points using the Pythagorean theorem, therefore occasionally being called the Pythagorean distance.
    The input are two list/arrays
    '''
### Function Silhouette_Coefficient(Ligands,  labels):
    '''
    Silhouette refers to a method of interpretation and validation of consistency within clusters of data. The technique provides a succinct graphical representation of how well each object has been classified.The silhouette value is a measure of how similar an object is to its own cluster (cohesion) compared to other clusters (separation). The silhouette ranges from −1 to +1, where a high value indicates that the object is well matched to its own cluster and poorly matched to neighboring clusters. If most objects have a high value, then the clustering configuration is appropriate. If many points have a low or negative value, then the clustering configuration may have too many or too few clusters.
    In my analysis, given the one point cluster will siginificant affect the result intepretation, for calculation I only consider the cluster that has more than one point.
    The input are a list of ligand class, and the labels are calculated through the clustering calculation
    '''
### Class and_Index(label_1, label_2)
    '''
    The Rand index or Rand measure (named after William M. Rand) in statistics, and in particular in data clustering, is a measure of the similarity between two data clusterings
    The input are labels calculating from clustering algorithms
    '''
###  class ClusterNode(self, vec, left=None, right=None, distance=-1, id=None, count=1, check = 0)

#### Parameters:

vec: the value for the node, after linking the two nodes, the new node will be generated via the average linkage
left: the smaller numder is put into the left node 
right:  the bigger number is put into the right node
param distance: the distance is sotred for ploting the distance
param id: the ID of the node
param count: the number of nodes that under this node
param check: checking whether this node have been identified

### class Ligand(Ligands)

#### Parameters:

Ligands: the input is a dataframe

#### Attributes:

ID: the ID of the ligand
Score : The score of the ligand
SMILES: SMLES format
OnBits: a array variable representing the Onbits information

#### Method:

get_ID(self):  # return the ID
get_scores(self): # return the scores
get_SMILES(self): # return the SMILES
get_OnBits(self): # return the onBits information

### class HierarchicalClustering(self, ligands,k):

#### paramemters

Ligands : a list of ligands
k : the number of cluster you decided

#### attributes

ID: The ID list
scores: Scores list   
dataset: the score matrix
labels: the cluster label list 
k: the number of cluster you decided
	    
#### methods	    
	    
* Initialization

input are a list of ligands and the decided cluster number

* minDist(self,dataset)

only used when you are using a small sample size to visilize the Hierarchical clustering plot
Go over all the nodes and find the smallest distance

* Implement_viz(self)
* show(self)

Visualize the hirarchical clusering, only for a small sample size.
First use Implement_viz() and then use show()  
Don't use Implement_viz() for the true analysis with a huge sample

* implement(self)

This is the formal implement function, this function use the parameter K, once the clusters merges into the
into K clusters then this function will break and return the labels for all the elements

* calc_label(self)

Return the analysis result, calculating labels

* leaf_traversal(self)

Travel across all nodes

### PartitionClustering(dataset, k, iteration)

K means is much eaiser considering it don't save the linkage between every two nodes

#### paramemters

Ligands : a list of ligands
k : the number of cluster you decided
iteration : the number of itertation

#### attributes

ID: The ID list
scores: Scores list   
dataset: the score matrix
labels: the cluster label list 
k: the number of cluster you decided
C: use to store all the nodes  

#### Methods

* Initilization

input are a list of ligands and the decided cluster number and the number of iteration

* Implement

Implement() and that's it!

## Example

In this turtorial, we will together initialize a global alignment class and a local alignment class, compute the alignment scores and alignment matrix.
Go to the folder test, open a jupyter notebook. ( Or other editors...)

* Initialization

Load packages

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


Take a subset of all the ligands and perform kmeans and hirarchical clustering

    csv = pd.read_csv('../ligand_information.csv', index_col= False)
    Data = list(np.random.randint(1,csv.shape[0],size=100))
    Ligands = []
    for i in Data:

        a_ligand = Ligand(csv[csv.index==i])
        Ligands.append(a_ligand)

		
Perform hirarchical clustering, and visualize the results

    Hiera_C = HierarchicalClustering(Ligands, 5)
    results = Hiera_C.implement_viz()
    Hiera_C.show(results, len(Ligands))
    # run Hiera_C.implement_viz() and Hiera_C.show() to visualize the results. But this can run very slow so please use a small subset.
    
Perform hirarchical clustering, keams clustering

    # Can use a higher 
    Data = list(np.random.randint(1,csv.shape[0],size=1000))
    Ligands = []
    for i in Data:

        a_ligand = Ligand(csv[csv.index==i])
        Ligands.append(a_ligand)
    
    
    Hiera_C = HierarchicalClustering(Ligands, 5)
    Hiera_C.implement()
    
    Kmeans_C = PartitionClustering(Ligands, 5, 10)
    Kmeans_C.implement()
    
    

	  
	




