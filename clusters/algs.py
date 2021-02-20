import pandas as pd
import os
import sys
import glob
import numpy as np
import random
import matplotlib.pyplot as plt
import queue
import math
import copy
import matplotlib.pyplot as plt
import statistics as stat
import json
from scipy.spatial import distance
def getEuclidean(point1, point2):
    '''
    I use the the Euclidean distance. It can be calculated from the Cartesian coordinates of the points using the Pythagorean theorem, therefore occasionally being called the Pythagorean distance.
    '''

    return distance.euclidean(point1,point2)

def Silhouette_Coefficient(Ligands,  labels):
    '''
    Silhouette refers to a method of interpretation and validation of consistency within clusters of data. The technique provides a succinct graphical representation of how well each object has been classified.The silhouette value is a measure of how similar an object is to its own cluster (cohesion) compared to other clusters (separation). The silhouette ranges from −1 to +1, where a high value indicates that the object is well matched to its own cluster and poorly matched to neighboring clusters. If most objects have a high value, then the clustering configuration is appropriate. If many points have a low or negative value, then the clustering configuration may have too many or too few clusters.
    
    
    In my analysis, given the one point cluster will siginificant affect the result intepretation, for calculation I only consider the cluster that has more than one point.
    
    '''
    
    # stacking the onbits information into a matrix
    dataset = []
    for ligand in Ligands:

        dataset.append(ligand.get_OnBits())
    
    dataset = np.stack((dataset), axis=0)
    
    # create the distance matrix, and put it into a distionary for further calculation
    nodes = [ClusterNode(vec=v, id=i) for i,v in enumerate(dataset)]
    Silhouette_coefs = []
    C     = set(labels)
    index = list(range(0, len(dataset)))
    # Just in case you want to build a distance for all ligands, there is a dictionary I made which would save time
    if len(index) != 8524:
        distances = {}
        for i in range(len(index) -1):
            for j in range(i + 1, len(index)):
                d_key = (nodes[i].id,nodes[j].id)
                distances[str(d_key)] = getEuclidean(nodes[i].vec,nodes[j].vec)
    else:
        with open('dic.json') as load_json:
            distances = json.load(load_json)
    
    # interation for each cluster
    for c in C:
        The_Cluster_index = []
        Other_Clusters_index  = []
        # Define the cluster as this cluster and others cluster
        for i in index:
            if labels[i] == c:
                
        #The_Cluster = dataset[labels == c]
                The_Cluster_index.append(i)
            else:
                Other_Clusters_index.append(i)
        The_Cluster = dataset[np.ix_(The_Cluster_index,)]
        Other_Clusters = dataset[np.ix_(Other_Clusters_index,)]
        
        
        
        Silhouette_a = []
        Silhouette_b = []
        # in cluster
        for i in range(The_Cluster.shape[0]-1):
            for j in range(i+1, The_Cluster.shape[0]):
                index_i = The_Cluster_index[i]
                index_j = The_Cluster_index[j]
                d_key = (index_i,index_j)
                
                Silhouette_a.append(distances[str(d_key)])
        # if there is only one ligand in the cluster, skip it  
        if(len(Silhouette_a) != 0):
            
            Silhouette_a = stat.mean(Silhouette_a)
        else:
            continue
        # out cluster
        for i in range(The_Cluster.shape[0]):
            for j in range(Other_Clusters.shape[0]):
                index_i = The_Cluster_index[i]
                index_j = Other_Clusters_index[j]
                if index_i < index_j :
                    d_key = (index_i, index_j)
                    Silhouette_b.append(distances[str(d_key)])
                elif index_i > index_j:
                    d_key = (index_j, index_i)
                    Silhouette_b.append(distances[str(d_key)])                    
        Silhouette_b = stat.mean(Silhouette_b)
        Silhouette_coefs.append( (Silhouette_b - Silhouette_a)/max(Silhouette_b,Silhouette_a))
    return stat.mean(Silhouette_coefs)



def Rand_Index(label_1, label_2):
    '''
    The Rand index or Rand measure (named after William M. Rand) in statistics, and in particular in data clustering, is a measure of the similarity between two data clusterings
    '''
    y_true = label_1
    y_pred = label_2
    a = 0
    b = 0
    assert len(label_1) == len(label_2), 'The nodes are not the same!'  # interupte if the elements are not same
    n = len(label_1)
    # after the iteration, a is counted as the true negative and true positive
    # b is counted as the true positive, false positive, true negative and false negative
    for i in range(len(label_1)):
        for j in range(i+1, len(label_1)):
            if (y_true[i] == y_true[j]) & (y_pred[i] == y_pred[j]):
                a +=1
            elif (y_true[i] != y_true[j]) & (y_pred[i] != y_pred[j]):
                b +=1
            else:
                pass
    RI = (a + b) / (n*(n-1)/2)
    return RI

class ClusterNode(object):
	def __init__(self, vec, left=None, right=None, distance=-1, id=None, count=1, check = 0):
		"""
		param vec: the value for the node, after linking the two nodes, the new node will be generated via the average linkage
		param left: the smaller numder is put into the left node 
		param right:  the bigger number is put into the right node
		param distance: the distance is sotred for ploting the distance
		param id: the ID of the node
		param count: the number of nodes that under this node
		param check: checking whether this node have been identified
		"""
		self.vec = vec
		self.left = left
		self.right = right
		self.distance = distance
		self.id = id
		self.count = count
		self.check = check

'''
param ID: the ID of the ligand
param Score : The score of the ligand
param SMILES: SMLES format
param OnBits: a array variable representing the Onbits information
'''
class Ligand:
	# Initiliza the class
	def __init__(self, df):
		self.ID = df['LigandID'].values[0]
		self.Score = df['Score'].values[0]
		self.SMILES = df['SMILES'].values[0]
		OnBits = np.zeros(1024)
		#print(df['OnBits'].values[0])
		num_list = [ int(x) for x in df['OnBits'].values[0].split(',')]
		OnBits[num_list] = 1
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






class HierarchicalClustering:
	"""
	param ID: The ID list
	param scores: Scores list   
	param dataset: the score matrix
	param labels: the cluster label list 
	param k: the number of cluster you decided    
	"""
	def __init__(self, ligands,k): 
		self.ID = []
		self.scores = []
		dataset = []
		for ligand in ligands:
			self.ID.append(ligand.get_ID())
			self.scores.append(ligand.get_scores())
			dataset.append(ligand.get_OnBits())
		self.dataset = np.stack((dataset), axis=0)
		#self.iteration = iteration
		self.labels = []
		self.k = k
		
		# Initialize the predicted labels
		for i in range(len(self.dataset)):
				self.labels.append(-1)
	"""
	function minDist()
	only used when you are using a small sample size to visilize the Hierarchical clustering plot
	Go over all the nodes and find the smallest distance
	"""
	def minDist(self,dataset):
		
		mindist = math.inf
		for i in range(len(dataset)-1):
			if dataset[i].check == 1:
				#If the check is 1, then skip it
				continue
			for j in range(i+1,len(dataset)):
				if dataset[j].check == 1:
					continue
				dist = getEuclidean(dataset[i].vec,dataset[j].vec)
				if dist < mindist:
					mindist = dist
					x, y = i, j
		return mindist, x, y
	def implement_viz(self):
		# Visualize the hirarchical clusering, only for a small sample size.
		nodes = [ClusterNode(vec=item,id=[(chr(ord('a')+i))],count=1) for i,item in enumerate(self.dataset)]
		#Make all the elements as node cluster
		length = len(nodes)
		while(True):
			mindist, x, y = self.minDist(nodes)
			nodes[x].check = 1
			nodes[y].check = 1
			tmpid = copy.deepcopy(nodes[x].id)
			tmpid.extend(nodes[y].id)
			nodes.append(ClusterNode(vec=(nodes[x].vec+nodes[y].vec)/2,id=tmpid,\
				left=nodes[x],right=nodes[y],distance=mindist,count=nodes[x].count+nodes[y].count))
			#Generate new nodes
			if len(tmpid) == length:
				#When all nodes are reached, quit the loop
				break
		return nodes
	# Use this show function to visualize the plots
	def show(self,dataset,num):
		plt.figure(1)
		showqueue = queue.Queue()
		#build a queue
		showqueue.put(dataset[len(dataset) - 1])
		#first put into the start point
		showqueue.put(num)
		#store the location of the start point
		while not showqueue.empty():
			index = showqueue.get()
			#draw tyhe current point
			i = showqueue.get()
			#draw the x axis value
			left = i - (index.count)/2
			right = i + (index.count)/2
			if index.left != None:
				x = [left,right]
				y = [index.distance,index.distance]
				plt.plot(x,y)
				x = [left,left]
				y = [index.distance,index.left.distance]
				plt.plot(x,y)
				showqueue.put(index.left)
				showqueue.put(left)
			if index.right != None:
				x = [right,right]
				y = [index.distance,index.right.distance]
				plt.plot(x,y)
				showqueue.put(index.right)
				showqueue.put(right)
		plt.show()
	# This is the formal implement function, this function use the parameter K, once the clusters merges into the
	# into K clusters then this function will break and return the labels for all the elements
	def implement(self):       
		#  build all the Cluster nodes
		nodes = [ClusterNode(vec=v, id=i) for i,v in enumerate(self.dataset)]
		point_num, feature_num = self.dataset.shape[0], self.dataset.shape[1]
		
		currentclustid = -1
		#load a existed dictiondary if you want to use all the ligands for clustering, which can save time
		if len(nodes) != 8524:
			print('Building a new one')
			distances = {}
		else:
			with open('dic.json') as load_json:
				distances = json.load(load_json)
		while (len(nodes) > self.k) :
			
			
			
			min_dist = math.inf
			nodes_len = len(nodes)
			closest_part = None
			for i in range(nodes_len -1):
				for j in range(i + 1, nodes_len):
					d_key = (nodes[i].id,nodes[j].id)
					if str(d_key) not in distances:
						distances[str(d_key)] = getEuclidean(nodes[i].vec,nodes[j].vec)
					d = distances[str(d_key)]
					if d<= min_dist:
						min_dist = d
						closest_key = str(d_key)
						closest_part = (i,j)
			            
			part1, part2 = closest_part
			node1, node2 = nodes[part1], nodes[part2]
			# merge two node
			# once you generate a new node, the 
			new_vec = [ (node1.vec[i] * node1.count + node2.vec[i] * node2.count ) / (node1.count + node2.count)
						for i in range(feature_num)]
			new_node = ClusterNode(vec=new_vec,
								left=node1,
								right=node2,
								distance=min_dist,
								id=currentclustid,
								count=node1.count + node2.count)
			currentclustid -= 1
			del nodes[part2], nodes[part1]   # must del the larger index
			nodes.append(new_node)
			
		self.nodes = nodes
		self.calc_label()

	def calc_label(self):
		"""
		Return the analysis result
		"""
		for i, node in enumerate(self.nodes):
				
			self.leaf_traversal(node, i)

	def leaf_traversal(self, node: ClusterNode, label):
		"""
		Travel across all nodes
		"""
		if node.left == None and node.right == None:
			self.labels[node.id] = label
		if node.left:
			self.leaf_traversal(node.left, label)
		if node.right:
			self.leaf_traversal(node.right, label)



class PartitionClustering:
	"""
	
	K means is much eaiser considering it don't save the linkage between every two nodes
	param ID: The ID list
	param scores: Scores list   
	param dataset: the score matrix
	param labels: the cluster label list 
	param k: the number of cluster you decided
	param C: use to store all the nodes  
	    
	    
	"""
	def __init__(self, dataset, k, iteration):
		self.ID = []
		self.scores = []
		self.dataset = []
		for ligand in dataset:
			self.ID.append(ligand.get_ID())
			self.scores.append(ligand.get_scores())
			self.dataset.append(ligand.get_OnBits())
			
		self.dataset = np.stack((self.dataset), axis=0)
		#self.iteration = iteration
		self.iteration = iteration
		self.k = k
		self.labels = []
		self.C = []
		# Initialize the predicted labels
		for i in range(len(dataset)):
			self.labels.append(-1)
			
	def implement(self):
		index = random.sample(list(range(len(self.dataset))), self.k)
		vectors = []
		for i in index:
			vectors.append(self.dataset[i])
		# Perform k means based on the iteration
		while(self.iteration > 0 ):
			#print(self.iteration)
			self.C = []
			for i in range(self.k):
				self.C.append([])
			for labelIndex, item in enumerate(self.dataset):
				# Before deciding which class it came from, set it to -1

				classIndex = -1
				minDist = 1e6
				for i, point in enumerate(vectors):
					dist = getEuclidean(item, point)
					if(dist < minDist):
						classIndex = i
						minDist = dist
				self.C[classIndex].append(item)
				self.labels[labelIndex] = classIndex
			for i, cluster in enumerate(self.C):
				clusterHeart = []
				dimension = len(self.dataset[0])
				for j in range(dimension):
					clusterHeart.append(0)
				for item in cluster:
					for j, coordinate in enumerate(item):
						clusterHeart[j] += coordinate / len(cluster)
				vectors[i] = clusterHeart
			self.iteration -= 1
		
