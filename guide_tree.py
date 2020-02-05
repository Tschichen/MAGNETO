from Bio.Phylo.TreeConstruction import _Matrix, DistanceTreeConstructor, DistanceMatrix
from random import *
from io import StringIO
import networkx as nx
import matplotlib.pyplot as plt
from Bio import Phylo

class Guide_tree():
	"""
	A class used to represent guide trees for the progressive multiple graph alignment (MGA).

		Attributes
		----------
		network : networkX DiGraph representing the guide tree for MGA
		root : node representing the root of the guide tree

		Methods
		-------
		tree_drawing 
			visualizes guide tree to .png file
		newick_output
			string of newick format of the guide tree
		alignment_order
			order for progressive MGA based on the guide tree
	"""

	def __init__(self, network):
		"""
			Parameters
			----------
			network : networkX DiGraph
		"""
		self.network = network
		for v in self.network.nodes:
			if self.network.in_degree(v) == 0:
				self.root = v
				break

	def tree_drawing (self,out_file):
		"""Prints the guide tree.
		"""
		newick = self.newick_output()
		handle = StringIO(newick)
		tree = Phylo.read(handle, "newick")
		tree.ladderize()
		Phylo.draw(tree, ylabel=["graphs"],do_show=False)
		plt.savefig(out_file)
		plt.clf()

	def newick_output(self):
		"""Generates Newick String of the guide tree.
		"""
		newick = ""
		newick += str(self.newick_add_to_string(self.root))
		return newick
			
	def newick_add_to_string(self, node):
		"""Help Function that adds recursively nodes to the newick string.
	
			Parameters
			----------
			node : networkX DiGraph node
		"""
		if self.network.out_degree(node) == 0:
			return str(self.network.nodes[node]['name']) + ":" + str(node.branch_length)
		else:
			subtree = ""
			successorsIter = self.network.successors(node)
			u = next(successorsIter, None)
			if u:
				subtree += "(" + self.newick_add_to_string(u) + ","
			w = next(successorsIter, None)
			if w:
				subtree += self.newick_add_to_string(w) + ")" + str(self.network.nodes[node]['name']) + ":" + str(node.branch_length)
			return subtree

	def add_name_to_node(self, node):
		""" Help function that adds attriute 'name' to nodes of guide tree.
	
			Parameters
			----------
			node : networkX DiGraph node
		"""
		if self.network.out_degree(node) == 0:
			self.network.nodes[node]['name'] = node.name
		else:
			left_name = ""
			right_name = ""
			successorsIter = self.network.neighbors(node)
			u = next(successorsIter, None)
			if u:
				left_name = self.network.nodes[u]['name']
			w = next(successorsIter, None)
			if w:
				right_name = self.network.nodes[w]['name']
			self.network.nodes[node]['name'] = left_name + right_name

	def alignment_order(self):
		"""Returns list for alignment order based on the guide tree.
		"""
		dfs_postorder_iter = nx.dfs_postorder_nodes(self.network, source=self.root)
		for v in dfs_postorder_iter:
			self.add_name_to_node(v)
		dfs_postorder_iter_2 = nx.dfs_postorder_nodes(self.network, source=self.root)
		
		alignment_order = []
		for u in dfs_postorder_iter_2:
			if self.network.out_degree(u) != 0:
				triple = self.make_triple(u)
				alignment_order.append(triple)
		return alignment_order
			
	def make_triple(self, node):
		"""makes value triple for each inner node of the guide tree.
	
			Parameters
			----------
			node : inner node of the guide tree.
		
			Output
			------
			triple (left child name, right child name, inner node name)
		"""
		successorsIter = self.network.successors(node)
		left_name = ""
		right_name = ""
		u = next(successorsIter, None)
		if u:
			left_name = self.network.nodes[u]['name']
		w = next(successorsIter, None)
		if w:
			right_name = self.network.nodes[w]['name']
		node_name = self.network.nodes[node]['name']
		node_triple = (left_name, right_name, node_name)
		return node_triple

	def alignment_order_drawing(self):
		"""help function to represent the alignment order for guide tree drawing.
		"""
		newick = self.newick_output()
		handle = StringIO(newick)
		tree = Phylo.read(handle, "newick")
		leafs = tree.get_terminals(order='level')
		alignment_order_drawing = []
		for entry in leafs:
			alignment_order_drawing.append(entry.name)
		alignment_order_drawing.reverse()
		
		return alignment_order_drawing

class Guide_tree_Generator:
	"""
	A class used to generate guide trees out of results from pairwise graph alignments.

		Methods
		-------
		tree_from_scores
			generates guide tree via a distance matrix and UPGMA algorithm for MGA from pairwise scorings of input graphs
		tree_from_newick
			makes guide tree object from user input of a guide tree in newick format
		tree_from_random
			generates random guide tree
	"""

	def tree_from_scores(list_with_scores):
		"""Generates Guide_tree object from list of pairwise scoring input from graph matching algorithms.
	
			Parameters
			----------
			list with scores : scores from the pairwise alignments of the graphs. Example for three graphs a, b, c: [["a", "b", 2], ["a", "c", 4], ["b", "c", 3]]
		
			Output
			------
			Guide_tree object
		"""
		matrix = Guide_tree_Generator.score_to_matrix(list_with_scores)
		constructor = DistanceTreeConstructor()
		upgmatree = constructor.upgma(matrix)
		tree = Phylo.to_networkx(upgmatree)
		guide_tree = Guide_tree(tree)
		
		return guide_tree

	def tree_from_newick(path):
		"""Generates Guide_tree object from a newick tree entered by the user.
	
			Parameters
			----------
			path : path to newick string representing the desired aligning sequence for MGA
		
			Output
			------
			Guide_tree object
		"""
		tree = next(Phylo.parse(path, 'newick', rooted=True))
		networkx_tree = Phylo.to_networkx(tree)
		guide_tree = Guide_tree(networkx_tree)
		if Guide_tree_Generator.is_binary_tree(guide_tree) == True:
			return guide_tree
		else:
			print("The input is not a binary tree. Please enter a binary tree to get a guide tree")
		
	def is_binary_tree(self):
		"""Help function to check, whether tree input is a binary tree.
		"""
		for v in self.network.nodes:
			if self.network.out_degree(v) > 2:
				return False
		return True
	
	def tree_from_random(list_of_scores):
		"""Generates a random guide tree for MGA.
	
			Parameters
			----------
			list_of_scores : scores from the pairwise alignments of the graphs to get graph names. Example for three graphs a, b, c: [["a", "b", 2], ["a", "c", 4], ["b", "c", 3]]
		
			Output
			------
			Guide_tree object
		"""
		names = Guide_tree_Generator.make_graph_list(list_of_scores)
		matrix = Guide_tree_Generator.random_score_matrix(names)
		constructor = DistanceTreeConstructor()
		upgmatree = constructor.upgma(matrix)
		tree = Phylo.to_networkx(upgmatree)
		guide_tree = Guide_tree(tree)
		
		return guide_tree

	def score_to_matrix(list_with_scores):
		"""Help function that returns a distance matrix for guide tree generation.
		"""
		# lexikographic sort of list and graphs for further proceeding
		for i in range(len(list_with_scores)):
			graphs = [list_with_scores[i][0], list_with_scores[i][1]]
			graphs.sort()
			graphs.reverse()
			list_with_scores[i][0] = graphs[0]
			list_with_scores[i][1] = graphs[1]
		list_with_scores.sort()
		list_with_scores.reverse()
		# create name and score list for generation of distance matrix
		names = []
		scores = []
		i = -1
		for entry in list_with_scores:
			if entry[0] not in names:
				names.append(entry[0])
				scores.append([0, entry[2]])
				i += 1
			else:
				scores[i].append(entry[2])
		last_graph = list_with_scores[-1][1]
		names.append(last_graph)
		names.reverse()
		scores.append([0])
		scores.reverse()
		for line in scores:
			line.reverse()
		# generating distance matrix
		matrix = DistanceMatrix(names, scores)
		return matrix
	
	def random_score_matrix(list_of_graph_names):
		"""Help function that generates random scores for a distance matrix.
		"""
		scores = []
		k = 0
		for i in range (len(list_of_graph_names)):
			scores.append([0])
			if k > 0:
				for j in range (k):
					random_score = randint(1, 100) / 100
					scores[i].append(random_score)
			scores[i].reverse()
			k += 1
		# generating distance matrix
		matrix = DistanceMatrix(list_of_graph_names, scores)
		return matrix

	def make_graph_list(list_with_scores):
		"""Help function that generates lexikographic sorted List of graphs for random_tree function.
		"""
		names = []
		for i in range (len(list_with_scores)):
			graph_1 = list_with_scores[i][0]
			graph_2 = list_with_scores[i][1]
			if graph_1 not in names:
				names.append(graph_1)
			if graph_2 not in names:
				names.append(graph_2)
		names_sorted = sorted(names, key=str.casefold)
		return names_sorted
