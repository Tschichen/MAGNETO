import networkx as nx
import math
import random
import migraine.GraphIO as io
import matplotlib.pyplot as plt
import migraine.keep as keep


class Matching():
    """A class used to execute progressive multiple graph alignment (MGA) via Bron Kerbosch algorithm (bk) and to visualize graph matching objects.
        Attributes
        ----------
        network : networkX Graph or DiGraph
        match_nodes : nodes of all graphs, that match at this node. If no node of graph matches here, match_node = "None"
        match_edges : edges of all graphs, that match at this edge. If no edge of graph matches here, match_edge = "None"
        Methods
        -------
        modular_product
            creates a modular graph product out of two matching objects
        BK
            implementation of bk for pairwise alignments, pivoting strategy: costum pivot element, creates list of all maximal cliques in the modular graph product
        BK_random
            bk algorithm fpr pairwise alignment, handling with random pivot element, creates list of all maximal cliques in the modular graph product
        bk_multiple_alignment
            performing the bk for two matching objects in the MGA, pivoting: either costum or random. Returns Matching object with the two alignes objects.
        bk_multiple_alignment_clique
            performing the bk for two matching objects in the MGA while each matching object contains a specific clique according to user input. Returns Matching object with the two alignes objects.
        MGA_bk
            function that handles the MGA with bk for all input graphs according to the guide tree. Different settings are possible (scoring, input of a pre-clique or random pivoting).
        draw_matching
            visualizes matching object via pyplot
    """

    def __init__(self, network):
        """
            Parameters
            ----------
            network : networkx Graph or DiGraph object
            Attributes
            ----------
            match_nodes : nodes of all graphs, that match at this node are collected here. If no node of graph matches at this node, match_node = "None"
            match_edges : edges of all graphs, that match at this edge are collected here. If no edge of graph matches at this edge, match_edge = "None"
        """
        self.network = network
        for node in self.network.nodes(data=True):
            match_nodes = [(node[0], node[1]['Label'], node[1]['Graph'])]
            self.network.nodes[node[0]]['Matches'] = match_nodes
        for edge in self.network.edges(data=True):
            match_edges = [(edge[0], edge[1], edge[2]['Label'], edge[2]['Graph'])]
            self.network.edges[edge[0], edge[1]]['Matches'] = match_edges

    #Help function, that compares edges in an undirected graph
    edgecompare = lambda x, y: x[0] == y[0] and x[1] == y[1] or x[1] == y[0] and x[0] == y[1]

    def modular_product(self, matching_b):
        """Function that creates modular graph product out of two matching objects and returns it.
        """
        graph_product = nx.Graph()
        if self.network.is_directed() == matching_b.network.is_directed():
            graph_product.add_nodes_from(nx.cartesian_product(self.network, matching_b.network).nodes(True))
            for i in range(len(graph_product.nodes)):
                for j in range(i + 1, len(graph_product.nodes)):
                    # Check whether a node would map to itself
                    if list(graph_product.nodes)[i][0] != list(graph_product.nodes)[j][0] and \
                            list(graph_product.nodes)[i][1] != list(graph_product.nodes)[j][1]:
                        # check whether nodes are connected in the same way in the original graph
                        if self.is_connected_how(list(graph_product.nodes)[i][0], list(graph_product.nodes)[j][0],
                                                 self.network) == \
                                self.is_connected_how(list(graph_product.nodes)[i][1], list(graph_product.nodes)[j][1],
                                                      matching_b.network):
                            label = []
                            # add edge labels to edges in the modular graph product
                            for edge in self.network.edges(data=True):
                                if Matching.edgecompare(edge, [list(graph_product.nodes)[i][0],
                                                               list(graph_product.nodes)[j][0]]):
                                    attributes = self.network.get_edge_data(edge[0], edge[1])
                                    if 'Label' in attributes:
                                        label_1 = attributes['Label']
                                    else:
                                        label_1 = None
                                    label.append(label_1)
                                    break
                            else:
                                label.append(None)
                            for edge in matching_b.network.edges(data=True):
                                if Matching.edgecompare(edge, [list(graph_product.nodes)[i][1],
                                                               list(graph_product.nodes)[j][1]]):
                                    attributes = matching_b.network.get_edge_data(edge[0], edge[1])
                                    if 'Label' in attributes:
                                        label_2 = attributes['Label']
                                    else:
                                        label_2 = None
                                    label.append(label_2)
                                    break
                            else:
                                label.append(None)
                            graph_product.add_edge(list(graph_product.nodes)[i], list(graph_product.nodes)[j],
                                                   Label=label)

            return graph_product

        else:
            return -1

    def is_connected_how(self, node_a, node_b, graph):
        """Help function for modular product that tests, how two nodes are connected in a graph.
        """
        if node_a not in graph.neighbors(node_b) and node_b not in graph.neighbors(node_a):
            return 0
        elif node_a not in graph.neighbors(node_b) and node_b in graph.neighbors(node_a):
            return 1
        elif node_a in graph.neighbors(node_b) and node_b not in graph.neighbors(node_a):
            return 2
        else:
            return 3

    def BK_random(self, graph, clique, candidates, checked, clique_list):
        """Implementation of bk for pairwise alignments, pivoting strategy: random pivot element, calculates recursively all maximal cliques in a graph.

            Parameters
            ----------
            graph : modular graph product
            clique : complete graph that represents the nodes in the clique, that is investigated in the moment
            candidates : set of candidate nodes that are not investigated yet
            checked : already checked nodes, that do not expand the clique
            clique_list : list that collects all found maximal cliques
        """
        if not bool(checked) and not bool(candidates):
            clique_list.append(graph.subgraph(clique.nodes))

        else:
            # select random pivot element
            pivot = self.pivot_random(candidates, checked)
            neighbors_of_p = [x for x in graph.neighbors(pivot)]
            for i in candidates[:]:
                new_clique = nx.Graph()
                new_clique.add_nodes_from(clique.nodes)

                if i not in neighbors_of_p:
                    if i not in new_clique.nodes:
                        new_clique.add_node(i)
                    new_candidates = [x for x in candidates if x in graph.neighbors(i)]
                    new_checked = [x for x in checked if x in graph.neighbors(i)]
                    self.BK(graph, new_clique, new_candidates, new_checked, clique_list)
                    candidates.remove(i)
                    if i not in checked:
                        checked.append(i)

    def pivot_random(self, candidates, checked):
        """Help function that selects a random pivot element out of the candidate node set.
        """
        alle_elem = [x for x in checked]
        for i in candidates:
            if i not in alle_elem:
                alle_elem.append(i)
        pivot = random.choice(alle_elem)
        return pivot

    def BK(self, graph, clique, candidates, checked, clique_list):
        """Implementation of bk for pairwise alignments, pivoting strategy: costum pivot element.

            Creates list of all maximal cliques in the modular graph product.

            Parameters
            ----------
            graph : modular graph product
            clique : complete graph that represents the nodes in the clique, that is investigated in the moment
            candidates : set of candidate nodes that are not investigated yet
            checked : already checked nodes, that do not expand the clique
            clique_list : list that collects all found maximal cliques
        """
        if not bool(checked) and not bool(candidates):
            clique_list.append(graph.subgraph(clique.nodes))

        else:
            pivot = self.pivot_custom(graph, candidates, checked)
            neighbors_of_p = [x for x in graph.neighbors(pivot)]
            for i in candidates[:]:
                new_clique = nx.Graph()
                new_clique.add_nodes_from(clique.nodes)

                if i not in neighbors_of_p:
                    if i not in new_clique.nodes:
                        new_clique.add_node(i)
                    new_candidates = [x for x in candidates if x in graph.neighbors(i)]
                    new_checked = [x for x in checked if x in graph.neighbors(i)]
                    self.BK(graph, new_clique, new_candidates, new_checked, clique_list)
                    candidates.remove(i)
                    if i not in checked:
                        checked.append(i)

    def pivot_custom(self, graph, candidates, checked):
        """Help function that selects a costum pivot element.

            The node with the highest degree in a set union of the candidate and checked node sets is selected as pivoting element.
        """
        alle_elem = [x for x in checked]
        for i in candidates:
            if i not in alle_elem:
                alle_elem.append(i)
        count_neighbours = []
        for i in range(len(alle_elem)):
            neighbors = [i for i in graph.neighbors(alle_elem[i])]
            count_neighbours.append(len(neighbors))

        pivot = alle_elem[count_neighbours.index(max(count_neighbours))]
        return pivot

    def find_largest_clique(self, list_of_cliques):
        """Help function that finds and returns largest maximall clique out of the list of all found maximal cliques.
        """
        max_index = 0
        max_len = 0

        for i in range(len(list_of_cliques)):
            if len(list_of_cliques[i]) > max_len:
                max_index = i
                max_len = len(list_of_cliques[i])

        largest_clique = list_of_cliques[max_index]

        return largest_clique

    def add_nodes_to_alignment(self, matching_b, max_clique):
        """Help function that addes nodes and edges to the matching object according to the found max_clique.

            If node of matching b is in max_clique, then the node is added to the 'Matches' array of the corresponding
            node from self matching object. Otherwise a new node is appended to the matching object. Edges are added in
            the same way. The labels correspond to the label of the first node or edge added to the matching object.
            If a node or edge of the self matching object are not part of the matching clique, None is added to
            the 'Matches' array.
        """
        matched_nodes_matching_b = []
        len_matches_nodes = 0
        len_matches_edges = 0
        matched_edges_b = []
        matched_edges_self = []

        for match in max_clique:
            node_align = match[0]
            node_matching_b = match[1]
            for node in self.network.nodes(data=True):
                match_nodes = node[1]['Matches']
                if node[0] == node_align:
                    for noode in matching_b.network.nodes(data=True):
                        matching_b_name = noode[1]['Graph']
                        match_graph_name = node[1]['Graph']
                        if noode[0] == node_matching_b and matching_b_name != match_graph_name:
                            matches_noode = noode[1]['Matches']
                            for entry in matches_noode:
                                match_nodes.append(entry)
                            matched_nodes_matching_b.append(noode)
                            len_matches_nodes = len(match_nodes)
                            node[1]['Matches'] = match_nodes

        for edge in max_clique.edges(data=True):
            edge_self = (edge[0][0], edge[1][0])
            edge_b = (edge[0][1], edge[1][1])
            if edge_self in self.network.edges:
                match_edges = matching_b.network.edges[edge_b]['Matches']
                edge_self_matches = self.network.edges[edge_self]['Matches']
                for entry in match_edges:
                    edge_self_matches.append(entry)
                self.network.edges[edge_self]['Matches'] = edge_self_matches
                matched_edges_b.append(edge_b)
                matched_edges_self.append(edge_self)
                len_matches_edges = len(edge_self_matches)

        for edge in self.network.edges(data=True):
            match_edges = edge[2]['Matches']
            difference = len_matches_edges - len(match_edges)
            if len(match_edges) < len_matches_edges:
                for i in range(difference):
                    match_edges.append(None)
                    self.network.edges[edge[0], edge[1]]['Matches'] = match_edges

        for node in self.network.nodes(data=True):
            match_nodes = node[1]['Matches']
            difference = len_matches_nodes - len(match_nodes)
            if len(match_nodes) < len_matches_nodes:
                for i in range(difference):
                    match_nodes.append(None)
                    self.network.nodes[node[0]]['Matches'] = match_nodes

        not_matched_nodes_names = []

        for node in matching_b.network.nodes(data=True):
            node_name = node[0]
            graph_name = node[1]['Graph']
            for matched_node in matched_nodes_matching_b:
                node_matched_name = matched_node[0]
                graph_matched_name = matched_node[1]['Graph']
                if node_matched_name == node_name and graph_matched_name == graph_name:
                    break
            else:
                self.network.add_node(node[0], Label=node[1]['Label'], Graph=node[1]['Graph'])
                match_nodes = []
                match_nodes_b = node[1]['Matches']
                for i in range(len_matches_nodes - len(match_nodes_b)):
                    match_nodes.append(None)
                for entry in match_nodes_b:
                    match_nodes.append(entry)
                self.network.nodes[node[0]]['Matches'] = match_nodes
                not_matched_nodes_names.append(node_name)

        for edge in matching_b.network.edges(data=True):
            node_1 = edge[0]
            node_2 = edge[1]
            if 'Label' in edge[2]:
                label = edge[2]['Label']
            else:
                label = None
            graph = edge[2]['Graph']
            # both nodes not in matching subgraph
            if node_1 in not_matched_nodes_names:
                if node_2 in not_matched_nodes_names:
                    match_edges = []
                    match_edges_b = matching_b.network.edges[node_1, node_2]['Matches']
                    for i in range(len_matches_edges - len(match_edges_b)):
                        match_edges.append(None)
                    for entry in match_edges_b:
                        match_edges.append(entry)
                    self.network.add_edge(node_1, node_2, Label=label, Graph=graph, Matches=match_edges)
                else:
                    # node1 not matched, node2 part of matching subgraph
                    for match in max_clique:
                        match_node_align = match[0]
                        match_node_matching_b = match[1]
                        if match_node_matching_b == node_2:
                            match_edges = []
                            match_edges_b = matching_b.network.edges[node_1, node_2]['Matches']
                            for i in range(len_matches_edges - len(match_edges_b)):
                                match_edges.append(None)
                            for entry in match_edges_b:
                                match_edges.append(entry)
                            self.network.add_edge(node_1, match_node_align, Label=label, Graph=graph,
                                                  Matches=match_edges)
            # node2 is not matched, node1 is matched
            if node_2 in not_matched_nodes_names:
                for match in max_clique:
                    match_node_align = match[0]
                    match_node_matching_b = match[1]
                    if match_node_matching_b == node_1:
                        match_edges = []
                        match_edges_b = matching_b.network.edges[node_1, node_2]['Matches']
                        difference = len_matches_edges - len(match_edges_b)
                        for i in range(difference):
                            match_edges.append(None)
                        for entry in match_edges_b:
                            match_edges.append(entry)
                        self.network.add_edge(node_2, match_node_align, Label=label, Graph=graph, Matches=match_edges)

        return self

    def bk_multiple_alignment(self, matching_b, save_all, random_pivot=False, score=False, user_list_score=None,
                              nodes_edges_both=None, forbidden =False):
        """Performing the bk for two matching objects in the MGA, pivoting: either costum or random.

            Parameters
            ----------
            matching_b : second matching object, that is to be aligned with the self matching object
            save_all : boolean, if True, all steps of MGA are saved.
            random_pivot : boolean, if True, pivoting strategy for bk is random
            score : boolean, if True, not largest common subgraph between matching object is used for alignment but highest scoring clique is selected.
            user_list_score : array with label scores that are used for scoring of node label matches, edge label matches or both.
            nodes_edges_both : string. Specifies whether node labels, edge labels or both should be scored.
            forbidden : boolean, if True, forbidden label matches combinations of edges or nodes are removed from the modular graph product.

            Returns
            -------
            Returns Matching object with the two aligned objects.
        """
        mod_prod = self.modular_product(matching_b)
        if user_list_score:
            mod_prod = keep.modular_product_valid_matchings(mod_prod, user_list_score,nodes_edges_both)
        checked = []
        clique_list = []
        clique = nx.Graph()
        if random_pivot:
            self.BK_random(mod_prod, clique, list(mod_prod.nodes), checked, clique_list)
        else:
            self.BK(mod_prod, clique, list(mod_prod.nodes), checked, clique_list)

        if not score or forbidden:
            max_clique = self.find_largest_clique(clique_list)
            self = self.add_nodes_to_alignment(matching_b, max_clique)
        else:
            max_clique = keep.score_from_BK(clique_list, user_list_score, nodes_edges_both)
            self = self.add_nodes_to_alignment(matching_b, max_clique[1])
        if save_all:
            try:
                io.GraphIO.write_graphML_file(self, save_all)
            except:
                print("Could not write to " + save_all + " !")

        return self

    def bk_multiple_alignment_clique(self, matching_b, save_all, all_clique, score=False, user_list_score=None,
                                     nodes_edges_both=None, forbidden =False):
        """Performing the bk for two matching objects in the MGA while each matching object contains a specific clique according to user input.

            Parameters
            ----------
            matching_b : second matching object, that is to be aligned with the self matching object
            save_all : boolean, if True, all steps of MGA are saved as GraphML file
            all_clique : clique, that shpuld be included in every matching object.
            score : boolean, if True, not largest common subgraph between matching object is used for alignment but highest scoring clique is selected.
            user_list_score : array with label scores that are used for scoring of node label matches, edge label matches or both.
            nodes_edges_both : string. Specifies whether node labels, edge labels or both should be scored.

            Returns
            -------
            Returns Matching object with the two aligned objects.
       """
        mod_prod = self.modular_product(matching_b)
        if user_list_score:
            mod_prod = keep.modular_product_valid_matchings(mod_prod, user_list_score, nodes_edges_both)
        checked = []
        clique_list = []
        clique = nx.Graph()
        matching_name_a = list(self.network.nodes(data=True))[0][1]['Graph']
        matching_name_b = list(matching_b.network.nodes(data=True))[0][1]['Graph']

        clique1 = all_clique.clique_dict[matching_name_a]
        clique2 = all_clique.clique_dict[matching_name_b]
        clique = all_clique.get_nodes(clique1, clique2, mod_prod)
        clique = all_clique.get_edges(clique1, clique2, clique)
        candidates = all_clique.filter_candidates_mult(clique1, clique2, clique, mod_prod)

        clique_list.append(clique)
        self.BK(mod_prod, clique, candidates, checked, clique_list)
        if not score or forbidden:
            max_clique = self.find_largest_clique(clique_list)
            self = self.add_nodes_to_alignment(matching_b, max_clique)
        else:
            max_clique = keep.score_from_BK(clique_list, user_list_score, nodes_edges_both)
            self = self.add_nodes_to_alignment(matching_b, max_clique[1])
        if save_all:
            try:
                io.GraphIO.write_graphML_file(self, save_all)
            except:
                print("Could not write to " + save_all + " !")

        return self

    def MGA_bk(self, matching_object_dict, list_alignment_order, save_all=None, pre_clique=False, random_pivot=False,
               score=False, user_list_score=None, nodes_edges_both=None, forbidden =False):
        """Function that handles the MGA with bk for all input graphs according to the guide tree.

            Different settings are possible (scoring, input of a pre-clique or random pivoting).

            Parameters
            ----------
            matching_object_dict :
            list_alignment_order : input list for aligning order according to the guide tree. List must be in the following format: [(matching_a, matching_b, matching_amatching_b), (matching_amatching_b, matching_c, matching_matching_bmatching_c)]
            save_all : boolean, if True, all alginment steps are saved as GraphML file
            pre_clique : Option if a pre-clique should be part of every matching object
            random_pivot : boolean, if True, random pivoting strategy is used.
            score : boolean, if true, scoring strategy is used to determine best common clique between two matching objects.
            user_list_score : array with label scores that are used for scoring of node label matches, edge label matches or both.
            nodes_edges_both : string. Specifies whether node labels, edge labels or both should be scored.
        """

        for i in range(len(list_alignment_order)):
            if save_all:
                save_all_as = save_all + str(i + 1) + ".graphml"
            else:
                save_all_as = None
            matching_name_a = list_alignment_order[i][0]
            matching_name_b = list_alignment_order[i][1]
            print("alignment step: " + str(matching_name_a) + ", " + str(matching_name_b))
            matching_a = matching_object_dict[matching_name_a]
            matching_b = matching_object_dict[matching_name_b]
            if not random_pivot:
                if not pre_clique:
                    if not score:
                        matching_ab = matching_a.bk_multiple_alignment(matching_b, save_all_as, forbidden=forbidden)
                        matching_object_dict[list_alignment_order[i][2]] = matching_ab
                    else:
                        matching_ab = matching_a.bk_multiple_alignment(matching_b, save_all_as, False, score,
                                                                       user_list_score, nodes_edges_both,forbidden=forbidden)
                        matching_object_dict[list_alignment_order[i][2]] = matching_ab
                else:
                    if not score:
                        matching_ab = matching_a.bk_multiple_alignment_clique(matching_b, save_all_as, pre_clique,forbidden=forbidden)
                        matching_object_dict[list_alignment_order[i][2]] = matching_ab
                    else:
                        matching_ab = matching_a.bk_multiple_alignment_clique(matching_b, save_all_as, pre_clique, score, user_list_score, nodes_edges_both,forbidden=forbidden)
                        matching_object_dict[list_alignment_order[i][2]] = matching_ab
            else:
                if not score:
                    matching_ab = matching_a.bk_multiple_alignment(matching_b, save_all_as, random_pivot,forbidden=forbidden)
                    matching_object_dict[list_alignment_order[i][2]] = matching_ab
                else:
                    matching_ab = matching_a.bk_multiple_alignment(matching_b, save_all_as, random_pivot, score,
                                                                   user_list_score, nodes_edges_both, forbidden=forbidden)
                    matching_object_dict[list_alignment_order[i][2]] = matching_ab



    def draw_matching(self, alignment_order_list, label_node_name, file):
        """Function that visualizes the matching object.

            Parameters
            ----------
            alignment_order_list : format of alignment order list must be: [first_graph, second_graph, third_graph, ...]
            label_node_name : string. Represents, what should be printed: "label" is node label, "node_name" is name of the node.
            file : output file
        """

        pos = nx.spring_layout(self.network)
        nodelist_complete_match = []
        nodelist_else = []
        nodelabels = {}
        for node in self.network.nodes(data=True):
            node_name = node[0]
            matchings = []
            matches = node[1]['Matches']
            for match in matches:
                if match == None:
                    matchings.append("-")
                else:
                    try:
                        if label_node_name == "label":
                            matchings.append(match[1])
                        else:
                            matchings.append(match[0])
                    except:
                        pass
                nodelabels[node_name] = matchings
            if "-" not in nodelabels[node_name]:
                nodelist_complete_match.append(node[0])
            else:
                nodelist_else.append(node[0])

        matching_subgraph = self.network.subgraph(nodelist_complete_match)

        edgelist_matching = list(matching_subgraph.edges)
        edgelist = set(self.network.edges)
        for e in edgelist_matching:
            try:
                edgelist.remove(e)
            except:
                edgelist.remove((e[1], e[0]))

        nx.draw_networkx_nodes(self.network, pos, nodelist_complete_match, node_color='r', node_size=500, alpha=0.8)
        nx.draw_networkx_nodes(self.network, pos, nodelist_else, node_color='g', node_size=300, alpha=0.8)
        nx.draw_networkx_edges(self.network, pos, edgelist_matching, width=3, alpha=0.5, edge_color='r')
        nx.draw_networkx_edges(self.network, pos, edgelist, width=3, alpha=0.5, edge_color='b')
        nx.draw_networkx_labels(self.network, pos, nodelabels, font_size=8)
        plt.title(alignment_order_list)
        plt.axis('off')
        plt.savefig(file)
        plt.clf()
