import math
import networkx as nx
import matching as m
import GraphIO as io

"""
A script that contains functions to implement scoring of alignments depending on node and / or edge labels and integration of forbidden matches between certain node and / or edge labels.

	Methods
	-------
	parse_keep_input
		takes list with user input for scoring and / or forbidden match combination of labels and returns list of scores for either nodes or edges or both
	modular_product_valid_matchings
		removes forbidden node or edge matches from the modular graph product
	similarity_score_to_distance(graph_1, graph_2, score):
		converts similarity scores depending on the number of matched nodes in the pairwise alignments into distance values for creating the distance matrix
	upgma_score_for_list
		converts similarity scores depending on the scoring of subgraph of the pairwise alignments into distance values for creating the distance matrix
	score_from_BK
		calculates score of the cliques found via bk in the modular graph product and returns the highest scoring clique and the score
"""

def parse_keep_input(path):
	"""Takes a path to an user input file with information for scoring node and / or edge labels and returns list of scores for either nodes or edges or both; additional forbidden label combinations can be realized through this user input.
	
		Parameters
		----------
		path : path to file representing user scoring function for node and / or edge labels as well as forbidden matches if desired
		input must be in the layout described in template_keep_function.txt
		
		Output
		------
		returns array with scores for node- and edge labels [[node_scores], [edge_scores]], if only scores for edge or labels are in input file, array looks like: [[score_list], []]
	"""
	try:
		file = open(path, "r")
	except:
		print("File input needed for keep and scoring function")
		return False

	user_list_node_scores = []
	user_list_edge_scores = []
	node_label_array = []
	node_score_array = []
	edge_label_array = []
	edge_score_array = []

	input = file.readlines()
	hashtag = 0

	for i in range(len(input)):
		line = input[i]
		line_gestrippt = line.strip()
		if line_gestrippt.startswith("#"):
			hashtag += 1
		if hashtag < 2:
			if line_gestrippt.startswith("#"):
				line_gestrippt = line_gestrippt[1:]
				node_label_array = line_gestrippt.split("\t")
				for label in node_label_array:
					user_list_node_scores.append([label])
			else:
				node_scores = line_gestrippt.split("\t")
				node_score_array.append(node_scores)
		else:
			if line_gestrippt.startswith("#"):
				line_gestrippt = line_gestrippt[1:]
				edge_label_array = line_gestrippt.split("\t")
				for label in edge_label_array:
					user_list_edge_scores.append([label])
			else:
				edge_scores = line_gestrippt.split("\t")
				edge_score_array.append(edge_scores)

	# make score array for nodes
	for z in range(len(node_label_array)):
		j = 0
		k = 0
		score_list_per_label = []
		for j in range(len(node_score_array[z])):
			recent_label = node_label_array[k]
			score = node_score_array[z][j]
			recent_node_score_array = [recent_label]
			if score == "N":
				recent_node_score_array.append(-math.inf)
			else:
				recent_node_score_array.append(int(score))
			score_list_per_label.append(recent_node_score_array)
			
			k += 1
			j += 1
		
		if z > 0:
			for g in range(z):
				recent_label = node_label_array[k]
				score = node_score_array[z - g - 1][j - 1]
				recent_node_score_array = [recent_label]
				if score == "N":
					recent_node_score_array.append(-math.inf)
				else:
					recent_node_score_array.append(int(score))
				k += 1
				score_list_per_label.append(recent_node_score_array)
		
		user_list_node_scores[len(node_score_array) - z - 1].append(score_list_per_label)

	# make score array for edges
	z = 0
	for z in range(len(edge_label_array)):
		j = 0
		k = 0
		score_list_per_label = []
		for j in range(len(edge_score_array[z])):
			recent_label = edge_label_array[k]
			score = edge_score_array[z][j]
			recent_edge_score_array = [recent_label]
			if score == "N":
				recent_edge_score_array.append(-math.inf)
			else:
				recent_edge_score_array.append(int(score))
			score_list_per_label.append(recent_edge_score_array)

			k += 1
			j += 1

		if z > 0:
			for g in range(z):
				recent_label = edge_label_array[k]
				score = edge_score_array[z - g - 1][j - 1]
				recent_edge_score_array = [recent_label]
				if score == "N":
					recent_edge_score_array.append(-math.inf)
				else:
					recent_edge_score_array.append(int(score))

				k += 1
				score_list_per_label.append(recent_edge_score_array)

		user_list_edge_scores[len(edge_score_array) - z - 1].append(score_list_per_label)

	file.close()

	return [user_list_node_scores, user_list_edge_scores]

def modular_product_valid_matchings(modular_graphproduct, user_list_node_and_edge_scores, nodes_edges_both):
	"""Deletes nodes or edges or both with forbidden label matchings from modular graph product.

		Parameters
		----------
		modular_graphproduct : modular graph product of two graphs in networkX graph
		user_list_node_and_edge_scores : list generated in function parse_keep_input
		nodes_edges_both : type "nodes" if only forbidden nodes label matches should be deleted, "edges" for only edge deleting, and "both" for deleting forbidden nodes and edges
		
		Output
		------
		modular graph product without the forbidden label matches of nodes and / or edges as networkX graph
	"""
	nodes_delete = []
	edges_delete = []
	if nodes_edges_both == "nodes":
		forbidden_labels_nodes = find_forbidden_matches(user_list_node_and_edge_scores[0])
		for v in modular_graphproduct.nodes(data=True):
			label_1 = v[1]['Label'][0]
			label_2 = v[1]['Label'][1]
			if match_is_permitted(label_1, label_2, forbidden_labels_nodes) == False:
				nodes_delete.append(v[0])
		modular_graphproduct.remove_nodes_from(nodes_delete)
	elif nodes_edges_both == "edges":
		forbidden_labels_edges = find_forbidden_matches(user_list_node_and_edge_scores[1])
		for edge in modular_graphproduct.edges(data=True):
			label_1 = edge[2]['Label'][0]
			label_2 = edge[2]['Label'][1]
			if label_1 != None:
				if match_is_permitted(label_1, label_2, forbidden_labels_edges) == False:
					edges_delete.append(edge)
		modular_graphproduct.remove_edges_from(edges_delete)
	else:
		forbidden_labels_nodes = find_forbidden_matches(user_list_node_and_edge_scores[0])
		forbidden_labels_edges = find_forbidden_matches(user_list_node_and_edge_scores[1])
		for v in modular_graphproduct.nodes(data=True):
			label_1 = v[1]['Label'][0]
			label_2 = v[1]['Label'][1]
			if match_is_permitted(label_1, label_2, forbidden_labels_nodes) == False:
				nodes_delete.append(v[0])
		modular_graphproduct.remove_nodes_from(nodes_delete)
		for edge in modular_graphproduct.edges(data=True):
			label_1 = edge[2]['Label'][0]
			label_2 = edge[2]['Label'][1]
			if label_1 != None:
				if match_is_permitted(label_1, label_2, forbidden_labels_edges) == False:
					edges_delete.append(edge)
		modular_graphproduct.remove_edges_from(edges_delete)

	return modular_graphproduct

def find_forbidden_matches(user_list_scores):
	"""Help function for finding forbidden label matches in the list generated from parse_keep_input.
	
		Output
		------
		returns dict: key = label, value: forbidden dicts for key label
	"""
	forbidden_label_dict = dict()
	for label_list in user_list_scores:
		label_1 = label_list[0]
		forbidden_matches = set()
		for i in range(len(label_list[1])):
			if math.isinf(label_list[1][i][1]) == True:
				forbidden_matches.add(label_list[1][i][0])

		forbidden_label_dict.update({label_1: forbidden_matches})

	return forbidden_label_dict

def match_is_permitted(label_a, label_b, forbidden_label_dict):
	"""Help function for finding forbidden label matches in the list generated from parse_keep_input.
	
		Parameters
		----------
		label_a, label_b : node or edge labels
		forbidden_label_dict : dict generated in find_forbidden_matches
	
		Output
		------
		returns true, is match is allowed, false otherwise
	"""
	try:
		label_set = forbidden_label_dict[label_a]
	except KeyError:
		return True
	if label_b in label_set:
		return False
	else:
		return True

def similarity_score_to_distance(graph_1, graph_2, score):
	"""Converts similarity scores depending on the number of matched nodes in the pairwise alignments into distance values for creating the distance matrix.

		Parameters
		----------
		graph_1, graph_2 : networkX graphs
		score : int value from pairwise alignment of graph_1 and graph_2
	"""
	nodes_1 = len(graph_1.nodes)
	nodes_2 = len(graph_2.nodes)
	min_nodes = min(nodes_1, nodes_2)
	try:
		upgma_score = 1 - (score / min_nodes)
	except:
		print("Error. Can't score graph that is too small.")
		return 1

	return upgma_score

def upgma_score_for_list(list_with_scores):
	"""Takes a list of all pairwise similarity scores depending on the individual scoring of nodes, edges or both and converts them to a distance for generating a distance matrix.
	"""
	maximum = max(list_with_scores)
	upgma_scores = []
	for i in list_with_scores:
		try:
			score = 1 - (i/maximum)
			upgma_scores.append(score)
		except ZeroDivisionError:
			upgma_scores.append(1)

	return upgma_scores

def score_from_BK(list_of_cliques, user_list_or_arrays_scores, nodes_edges_both):
	"""Returns the highest scoring clique found via Bron Kerbosch algorithm (bk) in the modular graph product.

		Parameters
		----------
		list_of_cliques : list of all cliques found in the modular graph product of two graphs via bk algorithm
		user_list_or_arrays_scores : list generated in function parse_keep_input
		nodes_edges_both : type "nodes" if only nodes should be scored, "edges" for only edge scoring, and "both" for scoring nodes and edges
	
		Output
		------
		tuple of max scoring clique and score
	"""
	scores = []
	if nodes_edges_both == "nodes":
		for clique in list_of_cliques:
			score = 0
			for v in clique.nodes(data=True):
				label_1 = v[1]['Label'][0]
				label_2 = v[1]['Label'][1]
				single_score = find_score(label_1, label_2, user_list_or_arrays_scores[0])
				score += single_score
			scores.append((score, clique))
		max_score = max(scores, key=lambda x: x[0])

	elif nodes_edges_both == "edges":
		for clique in list_of_cliques:
			score = 0
			for v in clique.edges(data=True):
				label_1 = v[2]['Label'][0]
				label_2 = v[2]['Label'][1]
				if label_1 != None and label_2 != None:
					single_score = find_score(label_1, label_2, user_list_or_arrays_scores[1])
					score += single_score
			scores.append((score, clique))

		max_score = max(scores, key=lambda x: x[0])


	else:
		for clique in list_of_cliques:
			score = 0
			for v in clique.nodes(data=True):
				label_1 = v[1]['Label'][0]
				label_2 = v[1]['Label'][1]
				single_score = find_score(label_1, label_2, user_list_or_arrays_scores[0])
				score += single_score
			for v in clique.edges(data=True):
				label_1 = v[2]['Label'][0]
				label_2 = v[2]['Label'][1]
				if label_1 != "None" or label_2 != "None":
					single_score = find_score(label_1, label_2, user_list_or_arrays_scores[1])
					score += single_score
			scores.append((score, clique))
		max_score = max(scores, key=lambda x: x[0])

	return max_score



def find_score(label_1, label_2, user_list_scores):
	"""Help function that finds scores for two matched node or edge labels in user_score_list.

		If there is no score in the list, match is scored with 0.
	"""
	if label_1 == None or label_2 == None:
		return 0

	for label_list in user_list_scores:
		if label_list[0] == label_1:
			for i in range(0, len(label_list[1])):
				if label_2 == label_list[1][i][0]:
					return int(label_list[1][i][1])
					break
	else:
		return 0
