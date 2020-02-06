import networkx as nx
import math
import matplotlib.pyplot as plt
import inspect
import random
import sys
from GraphIO import GraphIO


class DirGraphAlign(object):

    def __init__(self, G1, G2, node_score_list=None, edge_score_list=None,forbid =False):

        self.G1 = G1
        self.G2 = G2
        self.G1_nodes = set(G1.nodes())
        self.G2_nodes = set(G2.nodes())
        self.G2_node_order = {n: i for i, n in enumerate(G2)}
        #self.G2_node_degorder = sorted(list(G2.degree()), key = lambda kgrad: kgrad[1], reverse = True)
        #self.G2_node_order_list = [n[0] for n in self.G2_node_degorder]
        #self.G2_node_order = {n: i for i, n in enumerate(self.G2_node_order_list)}
        self.node_score_list = node_score_list
        if self.node_score_list != None:
            self.evalNodeAttr = True
        else:
            self.evalNodeAttr = False
        self.edge_score_list = edge_score_list
        self.forbidden_labels = self.find_forbidden_matches()
        self.forbidden_edges = self.find_forbidden_edges()
        self.forbidden = forbid

        self.initialize()

    def initialize(self):
        """Create working dicts and attributes.

        core_1[n]: contains the index of the node in G2 paired with n, if n is in the mapping.
        core_2[m]: contains the index of the node in G1 paired with m, if m is in the mapping.
        in_1[n]:  > 0 if n is either in the mapped nodes of G1 (M_1) or in T_1{in}.
        out_1[n]:  > 0 if n is either in the mapped nodes of G1 (M_1) or in T_1{out}.
        in_2[m]:  > 0 if m is either in the mapped nodes of G2 (M_2) or in T_2{in}.
        out_2[m]:  > 0 if m is either in the mapped nodes of G2 (M_2) or in T_2{out}.
        The values stored represent the depth of the search tree when the node was added to the set.

        mapping: a copy of the core_1 dict, returned at the end to the user.
        """
        sys.setrecursionlimit(10000)
        self.core_1 = dict()
        self.core_2 = dict()

        self.in_1 = {}
        self.in_2 = {}
        self.out_1 = {}
        self.out_2 = {}

        self.state = MappingState(self, 0)

        # Provide a convenient way to access the isomorphism mapping.
        self.mapping = {}
        self.current_max_len = 0
        self.score = 0
        self.counting = 0  # variable to store temporary score
        self.runningStatelist = [0]

    def find_forbidden_matches(self):
        """Return dict with forbidden node matches from node label list"""
        if self.node_score_list:
            forbidden_label_dict = dict()
            for label_list in self.node_score_list:
                label_1 = label_list[0]
                forbidden_matches = set()
                for i in range(0, len(label_list[1])):
                    if math.isinf(label_list[1][i][1]) == True:
                        forbidden_matches.add(label_list[1][i][0])
                        forbidden_label_dict.update({label_1: forbidden_matches})
            return forbidden_label_dict
        else:
            return 0

    def find_forbidden_edges(self):
        """Return dict with forbidden node matches from node label list"""
        if self.edge_score_list:
            forbidden_edge_dict = dict()
            for label_list in self.edge_score_list:
                label_1 = label_list[0]
                forbidden_edge_match = set()
                for i in range(1, len(label_list)):
                    if math.isinf(label_list[i][0]) == True:
                        forbidden_edge_match.add(label_list[i][0])
                        forbidden_edge_dict.update({label_1: forbidden_edge_match})
            return forbidden_edge_dict
        else:
            return None

    def candidate_pairs_iter(self):
        """Iterator over candidate pairs of nodes in G1 and G2.

            All computations are done in the current state.
        """

        G1_nodes = self.G1_nodes
        G2_nodes = self.G2_nodes
        min_key = self.G2_node_order.__getitem__
        # max_grad = self.G2_node_order.__getitem__
        #alternative function?

        T1_out = sorted([node for node in self.out_1 if node not in self.core_1])
        T2_out = sorted([node for node in self.out_2 if node not in self.core_2])

        # If T1_out and T2_out contain nodes:
        # P(s) = T1_out x min (T2_out)
        if T1_out and T2_out:
            for node_2 in T2_out:
                #print(node_2)
            #node_2 = min(T2_out, key=min_key)
            # node_2 = alternative?
                for node_1 in T1_out:
                    if self.label_match(node_1, node_2):
                        pair = tuple((node_1, node_2))
                        yield pair
                #else:
                    #yield None
            else:
                yield None


        # In case T1_out and T2_out are both empty:
        elif not (T1_out or T2_out):
            T1_in = sorted([node for node in self.in_1 if node not in self.core_1])
            T2_in = sorted([node for node in self.in_2 if node not in self.core_2])
            #t1os = sorted(T1_in)

            # If T1_in and T2_in contain nodes:
            # P(s) = T1_in x min (T2_in)
            if T1_in and T2_in:
                #print("t1 and t2")
                #node_2 = min(T2_in, key=min_key)
                for node_2i in T2_in:
                    for node_1i in T1_in:
                        if self.label_match(node_1i, node_2i):
                            pair = tuple((node_1i, node_2i))
                            yield pair
                else:
                    yield None


            # In case all terminal sets are empty:
            # P(s) = (N_1 - M_1) x (min (N_2 - M_2))
            elif (len(self.core_2) != len(G2_nodes)):
                #print("all empty")
                #node_2 = min(G2_nodes - set(self.core_2), key=min_key)
                #Ps = list(G2_nodes - set(self.core_2))
                all_g2 = sorted([node for node in G2_nodes if node not in self.core_2])
                all_g1 = sorted([node for node in G1_nodes if node not in self.core_1])
                for node_2 in all_g2:
                    if node_2 not in self.core_2:
                #node_2 = alterantive function?
                        for node_1 in all_g1:
                            #if node_1 not in self.core_1:
                            if self.label_match(node_1, node_2):
                                pair = tuple((node_1, node_2))
                                yield pair
                else:
                    yield None
            else:
                yield None

        else:
            #print("no candidates")
            yield None

    def matcher(self):
        """Generator for alignments.

        Candidate pairs are generated by candidate_pairs_iter() and tested for syntatic feasibility. If feasible, a new state
        (partial mapping) is generated. This represents a new level in the depth-first search. New candidate pairs are calculated
        for that state. When no new candidates can be found, the current mapping is yielded to the calling function, and the
        current state object is deleted. All variables are restored to the values held in the next upper level of the search tree.
        """
        for x in self.candidate_pairs_iter():
            if x == None:
                G1_node = None
                G2_node = None
                #print(len(self.core_1))
                if len(self.core_1) < self.current_max_len:
                    self.state.restore()
                    yield 0
                elif len(self.core_1) != 0:
                    self.mapping = self.core_1.copy()
                    self.current_max_len = len(self.mapping)
                    if self.evalNodeAttr:
                        self.score = self.counting
                        self.state.restore()
                        yield tuple((self.score, self.mapping))
                    else:
                        self.state.restore()
                        yield self.mapping
                else:
                    self.state.restore()
                    yield 0
            else:
                G1_node = x[0]
                G2_node = x[1]
                if self.syntactic_feas(G1_node, G2_node):
                    newstatenum = self.state.statenumber + 1
                    self.runningStatelist.append(newstatenum)
                    newstate = self.state.__class__(self, newstatenum, G1_node, G2_node)
                    self.state = newstate

                    if self.evalNodeAttr and not self.forbidden:
                        self.counting += self.score_pair(G1_node, G2_node)

                    for gollum in self.matcher():
                        yield gollum

    def syntactic_feas(self, G1_node, G2_node):
        """Return True if the addition of G1_node and G2_node lead to a common induced subgraph of G1 and G2.

        Parameters:
        G1_node : node from G1 that is tested for matching with G2_node from G2
        G2_node : node from G2 that is tested for matching with G1_node from G1
        """

        # check self loops
        if self.G1.number_of_edges(G1_node, G1_node) != self.G2.number_of_edges(G2_node, G2_node):
            return False

        # For each predecessor n' of n in the mapping, the
        # corresponding node m' is a predecessor of m, and vice versa. Also,
        # the number of edges must be equal.
        for predecessor in self.G1.pred[G1_node]:
            if predecessor in self.core_1:
                if not (self.core_1[predecessor] in self.G2.pred[G2_node]):
                    return False
                elif self.G1.number_of_edges(predecessor, G1_node) != self.G2.number_of_edges(self.core_1[predecessor],
                                                                                              G2_node):
                    return False
                elif self.match_is_forbidden(self.G1[predecessor][G1_node]['Label'],
                                             self.G2[self.core_1[predecessor]][G2_node]['Label'], self.forbidden_edges):
                    return False

        for predecessor in self.G2.pred[G2_node]:
            if predecessor in self.core_2:
                if not (self.core_2[predecessor] in self.G1.pred[G1_node]):
                    return False
                elif self.match_is_forbidden(self.G2[predecessor][G2_node]['Label'],
                                             self.G1[self.core_2[predecessor]][G1_node]['Label'], self.forbidden_edges):
                    return False
                elif self.G1.number_of_edges(self.core_2[predecessor], G1_node) != self.G2.number_of_edges(predecessor,
                                                                                                           G2_node):
                    return False

        # For each successor n' of n in the mapping, the corresponding
        # node m' is a successor of m, and vice versa. Also, the number of
        # edges must be equal.
        for successor in self.G1[G1_node]:
            if successor in self.core_1:
                if not (self.core_1[successor] in self.G2[G2_node]):
                    return False
                elif self.G1.number_of_edges(G1_node, successor) != self.G2.number_of_edges(G2_node,
                                                                                            self.core_1[successor]):
                    return False
                elif self.match_is_forbidden(self.G1[G1_node][successor]['Label'],
                                             self.G2[G2_node][self.core_1[successor]]['Label'], self.forbidden_edges):
                    return False

        for successor in self.G2[G2_node]:
            if successor in self.core_2:
                if not (self.core_2[successor] in self.G1[G1_node]):
                    return False
                elif self.G1.number_of_edges(self.core_2[successor], G1_node) != self.G2.number_of_edges(successor,
                                                                                                         G2_node):
                    return False
                elif self.match_is_forbidden(self.G2[G2_node][successor]['Label'],
                                             self.G1[G1_node][self.core_2[successor]]['Label'], self.forbidden_edges):
                    return False

        # Otherwise, addition of this node pair constitutes a common induced subgraph!
        return True

    def label_match(self, G1_node, G2_node):
        """Check if labels of node match, in case label comparison is desired

        Parameters:
        G1_node : node from G1
        G2_node : node from G2

        Returns:
        returns True if Labels are not forbidden matches
        """
        if not self.evalNodeAttr:
            return True
        label1 = self.G1.nodes[G1_node]['Label']
        label2 = self.G2.nodes[G2_node]['Label']
        try:
            if label1 and label2:
                if label1 in self.forbidden_labels[label2]:
                    return False
                else:
                    return True
            else:
                print(
                    'you have requested to evaluate node labels, but at least one label was not defined. True is returned as default.')
                return True
        except KeyError:
            return True

    def match_is_forbidden(self, label_a, label_b, forbidden_label_dict):
        """Return True if edge labels are forbidden to match"""
        if forbidden_label_dict != None:
            try:
                label_set = forbidden_label_dict[label_a]
            except:
                return False
            if label_b in label_set:
                return True
            else:
                return False
        else:
            return False

    def score_pair(self, G1_node, G2_node):
        """Calculate score of matched nodes"""
        if not self.evalNodeAttr or self.forbidden:
            return 0
        label1 = self.G1.nodes[G1_node]['Label']
        label2 = self.G2.nodes[G2_node]['Label']
        if label1 and label2:
            for label_list in self.node_score_list:
                if label_list[0] == label1:
                    for i in range(0, len(label_list[1])):
                        if label2 == label_list[1][i][0]:
                            return int(label_list[1][i][1])
        else:
            print(
                'you have requested to evaluate node labels, but at least one label was not defined. Default score 0 is returned.')
            return 0

        print("Some labels on your graph have not been found in the scoring list. Default Score: 0")
        return 0

    def alig_lister(self):
        """Initialize aligner-object and return list with possible alignments"""
        self.initialize()
        allResults = [alig for alig in self.matcher() if alig != 0]
        return (allResults)

    def process_results(self, results):
        """Prune short and low-scored alignments from the result list

        Parameters
        results : a list of possible mappings

        Returns
        Returns one mapping of maximal length
        """
        if self.node_score_list and not self.forbidden:
            max_len_alig = max(results, key=lambda x : len(x[1]))
            cutoffScore = max_len_alig[0]
            cutoffLen = len(max_len_alig[1])
            pruned_Res = [result for result in results if
                          ((len(result[1]) >= cutoffLen) and (result[0] >= cutoffScore))]
            print('There are ' + str(len(pruned_Res)) + ' alignments with a score >= ' + str(
                cutoffScore) + ' and a length >= ' + str(cutoffLen))
        else:
            pruned_Res = sorted(results, key=len, reverse=True)
            number_of_possAligs = len(pruned_Res)
            cutoffLen = len(pruned_Res[0])
            for i in range(0, number_of_possAligs):
                if (len(pruned_Res[i]) < cutoffLen):
                    del pruned_Res[i:(number_of_possAligs + 1)]
                    break
            print('There are ' + str(len(pruned_Res)) + ' possible Alignments.')

        if len(pruned_Res) > 1:
            lucky_number = random.randrange(0, len(pruned_Res))
            selected_result = pruned_Res[lucky_number]
            self.giveOutput(selected_result)
            return (selected_result)
        else:
            self.giveOutput(pruned_Res[0])
            return (pruned_Res[0])

    def getkey(self, element):
        return (element[1])

    def giveOutput(self, table):
        """Print results to terminal. """
        if self.node_score_list and not self.forbidden:

            print('Graph1: ' + self.G1.graph['Name'] + ', ' + 'Graph2: ' + self.G2.graph['Name'] + ' Score: ' + str(
                table[0]) + ", Alignment length = " + str(len(table[1])))

        else:
            print('Graph1: ' + self.G1.graph['Name'] + ', ' + 'Graph2: ' + self.G2.graph[
                'Name'] + ", Alignment length = " + str(len(table)))


class MappingState(object):
    """Internal representation of state for the DirGraphAlign class.

    Attrbutes
    GM: The DirGraphAlign object
    statenumber: in case the user wants to keep track of the states
    G1_node: the node from G1 that was matched to G2_node and last added to the mapping
    G2_node: the node from G2 that was matched to G1_node and last added to the mapping
    depth: the depth of the search tree, indicated by the length of the mapping

    Methods
    restore : when the mapping cannot be extended, the MappingState object is reset to the values before the last pair
    of nodes was added

    This class is used only to store state specific data.
    """

    def __init__(self, GM, statenumber, G1_node=None, G2_node=None):
        """Initializes MappingState object."""
        self.GM = GM
        self.statenumber = statenumber

        # Initialize the last stored node pair.
        self.G1_node = None
        self.G2_node = None
        self.depth = len(GM.core_1)

        if G1_node is None or G2_node is None:
            #reset the class variables
            GM.core_1 = dict()
            GM.core_2 = dict()
            GM.in_1 = {}
            GM.in_2 = {}
            GM.out_1 = {}
            GM.out_2 = {}

        if G1_node is not None and G2_node is not None:
            # Add the node pair to the mapping.
            GM.core_1[G1_node] = G2_node
            GM.core_2[G2_node] = G1_node

            # Store the node that was added last.
            self.G1_node = G1_node
            self.G2_node = G2_node

            # update the other four vectors.
            self.depth = len(GM.core_1)

            for vector in (GM.in_1, GM.out_1):
                if G1_node not in vector:
                    vector[G1_node] = self.depth
            for vector in (GM.in_2, GM.out_2):
                if G2_node not in vector:
                    vector[G2_node] = self.depth

            # Updates for T_1^{in}
            new_nodes = set([])
            for node in GM.core_1:
                new_nodes.update([predecessor for predecessor in GM.G1.predecessors(node)
                                  if predecessor not in GM.core_1])
            for node in new_nodes:
                if node not in GM.in_1:
                    GM.in_1[node] = self.depth

            # Updates for T_2^{in}
            new_nodes = set([])
            for node in GM.core_2:
                new_nodes.update([predecessor for predecessor in GM.G2.predecessors(node)
                                  if predecessor not in GM.core_2])
            for node in new_nodes:
                if node not in GM.in_2:
                    GM.in_2[node] = self.depth

            # Updates for T_1^{out}
            new_nodes = set([])
            for node in GM.core_1:
                new_nodes.update([successor for successor in GM.G1.successors(node) if successor not in GM.core_1])
            for node in new_nodes:
                if node not in GM.out_1:
                    GM.out_1[node] = self.depth

            # Updates for T_2^{out}
            new_nodes = set([])
            for node in GM.core_2:
                new_nodes.update([successor for successor in GM.G2.successors(node) if successor not in GM.core_2])
            for node in new_nodes:
                if node not in GM.out_2:
                    GM.out_2[node] = self.depth

    def restore(self):
        """Delete the MappingState object and restore the class variables."""

        # First we remove the node that was added from the core vectors.
        if self.G1_node is not None and self.G2_node is not None:
            if len(self.GM.core_1) != 0:
                self.GM.core_1.popitem()
            if len(self.GM.core_2) != 0:
                self.GM.core_2.popitem()

        if len(self.GM.runningStatelist) != 0:
            self.GM.runningStatelist.pop()

        # Now we restore the other vectors.
        # Thus, we delete all entries added at this depth level of the search
        for vector in (self.GM.in_1, self.GM.in_2, self.GM.out_1, self.GM.out_2):
            for node in list(vector.keys()):
                if vector[node] == self.depth:
                    del vector[node]

        # counter is reset to value before adding the pair
        if self.GM.evalNodeAttr and not self.GM.forbidden:
            self.GM.counting -= self.GM.score_pair(self.G1_node, self.G2_node)


class VF2_GraphAligner(object):
    """Class to generate paired and multiple alignments of graphs from user input.

    Attributes
    graph_list:
    graph_dict:
    node_score_list:
    edge_score_list:
    mappings:
    scores:
    length:
    edge_scores:

    Methods
    vf2mult:
    vf2_pga:
    score_edges:
    laber_checker:
    find_match:
    disconnected_alig:
    add_matches_to_nodes:
    add_matches_to_edges:
    build_matched_graphs:
    output_generator:
    report_score:
    convert_mga_to_list:
    convert_mga_to_edgelist:


    """

    def __init__(self, graph_list, graph_dict=None, node_score_list=None, edge_score_list=None):
        self.graph_list = graph_list
        self.list_alignment_order = graph_list
        self.graph_dict = graph_dict
        self.node_score_list = node_score_list
        self.edge_score_list = edge_score_list
        self.mappings = [{}]
        self.scores = []
        self.length = 0
        self.edge_scores = []


    def vf2mult(self, save_all='end', save_path=None,forbidden=False):
        """Alignes multiple graphs

        Parameters:
        save_all: option for saving of the results
        save_path: path to which alignment graph is saved
        forbidden: a dictionary of labels which are forbidden to match

        Returns:
            Returns an alignment graph made up of all input graphs
        """
        for i in range(len(self.list_alignment_order)):
            if save_all == 'all':
                save_all_as = save_path + str(i + 1) + ".graphml"
            else:
                save_all_as = False
            graph_name_a = self.list_alignment_order[i][0]
            graph_name_b = self.list_alignment_order[i][1]
            print("alignment step: " + str(graph_name_a) + ", " + str(graph_name_b))
            # a und b aus dem dict suchen
            graph_a = self.graph_dict[graph_name_a]
            graph_b = self.graph_dict[graph_name_b]

            try:
                if nx.is_connected(graph_a) and nx.is_connected(graph_b):
                    graph_ab = self.find_match(graph_a, graph_b, self.node_score_list, self.edge_score_list,
                                               save_all_as, forbidden=forbidden)
                else:
                    graph_ab = self.disconnected_alig(graph_a, graph_b, self.node_score_list, self.edge_score_list,
                                                      save_all_as, forbidden=forbidden)
            except:
                graph_ab = self.find_match(graph_a, graph_b, self.node_score_list, self.edge_score_list,
                                               save_all_as, forbidden=forbidden)

            self.graph_dict[self.list_alignment_order[i][2]] = graph_ab
            if self.edge_score_list != None:
                self.edge_score_list.append(self.score_edges(graph_ab))
            self.length = len(self.mappings[-1])
        return (graph_ab)

    def vf2_pga(self,forbidden=False):
        """Generate a paired alignment of two graphs.

        Parameters
        forbidden: dictionary of labels which are forbidden to match

        Returns
        Returns an alignment graph of two input graphs
        """
        G1 = self.graph_list[0]
        G2 = self.graph_list[1]
        self.label_checker(G1, G2)

        try:
            if nx.is_connected(G1) and nx.is_connected(G2):
                paired_result = self.find_match(G1, G2, self.node_score_list, self.edge_score_list,
                                           forbidden=forbidden)
            else:
                paired_result = self.disconnected_alig(G1, G2, self.node_score_list, self.edge_score_list,
                                                  forbidden=forbidden)
        except:
            paired_result = self.find_match(G1, G2, self.node_score_list, self.edge_score_list,
                                           forbidden=forbidden)

        if self.edge_score_list != None:
            self.edge_scores.append(self.score_edges(paired_result))
            if self.scores:
                self.scores[0] += self.edge_scores[0]
            else:
                self.scores.append(self.edge_scores[0])
        self.length = len(self.mappings[-1])
        return (paired_result)

    def score_edges(self, graph):
        """Calculate score for matched edges.

        Parameters
        graph: the graph which edges are to be scored

        Returns
        Returns the sum of all edge scores
        """
        if self.edge_score_list == None:
            return 0
        sum_edge_score = 0
        for edge in graph.edges(data=True):
            try:
                label1 = graph.edges[edge]['Matchings'][0]  # stimmt das so?
            except:
                print('No matching edges found, default score 0 is returned.')
                return 0
            label2 = graph.edges[edge]['Matchings'][-1]
            for label_list in self.edge_score_list:
                if label_list[0] == label1:
                    for i in range(1, len(label_list)):
                        if label2 == label_list[i][0]:
                            sum_edge_score += int(label_list[i][1])
                            break
                    else:
                        print('No score was provided for this pair of labels')
                    break
            else:
                print('Label 1 was not found in edge scoring list')
        return (sum_edge_score)


    def label_checker(self, G1, G2):
        """Check if edges are correctly labelled, if not add a None label.

        Parameters
        G1 and G2: input graphs which are later to be matched

        """
        for edge in G1.edges():
            try:
                G1.edges[edge]['Label']
            except:
                G1.edges[edge]['Label'] = None
            else:
                continue
        for edge in G2.edges():
            try:
                G2.edges[edge]['Label']
            except:
                G2.edges[edge]['Label'] = None
            else:
                continue
        for node in G1.nodes():
            try:
                G1.nodes[node]['Label']
            except:
                G1.nodes[node]['Label'] = None
            else:
                continue
        for node in G2.nodes():
            try:
                G2.nodes[node]['Label']
            except:
                G2.nodes[node]['Label'] = None
            else:
                continue

    def find_match(self, G1, G2, node_score_list=False, edge_score_list=False, save_all=False, forbidden=False):
        """Create a DirGraphAlign object and generate alignment of input graphs. The alignment is added to the
        self.mappings list.

        Parameters
        G1 and G2: input graphs
        node_score_list: contains scores for node label matches
        edge_score_list: contains scores for edge label matches
        save_all: saving option
        forbidden: dict with labels which are forbidden to match

        Returns
        Returns alignment graph
        """

        if G1.is_directed() and G2.is_directed():
            inGraph1 = G1
            inGraph2 = G2

        elif (not G1.is_directed()) and (not G2.is_directed()):
            inGraph1 = G1.to_directed()
            inGraph2 = G2.to_directed()

        if node_score_list:
            alignment = DirGraphAlign(inGraph1, inGraph2, node_score_list,forbid=forbidden)
            all_alignments = alignment.alig_lister()
            if len(all_alignments) > 1:
                selected_alignment = alignment.process_results(all_alignments)
                self.mappings.append(selected_alignment[1])
                self.scores.append(selected_alignment[0])
            elif len(all_alignments) == 1:
                self.mappings.append(all_alignments[0][1])
                self.scores.append(all_alignments[0][0])
            else:
                print("No Common Induced Subgraph Found")
                self.scores.append(1.0)

        else:
            alignment = DirGraphAlign(inGraph1, inGraph2)
            all_alignments = alignment.alig_lister()
            if len(all_alignments) > 1:
                selected_alignment = alignment.process_results(all_alignments)
                self.mappings.append(selected_alignment)
            elif len(all_alignments) == 1:
                self.mappings.append(all_alignments[0])
            else:
                print("No Common Induced Subgraph Found")

        pair_or_mult = inspect.currentframe().f_back.f_code.co_name
        if pair_or_mult == "vf2mult":
            self.add_matches_to_nodes(self.mappings[-1], inGraph1, inGraph2)
            self.add_matches_to_edges(self.mappings[-1], inGraph1, inGraph2)
        graphAlignment = self.build_matched_graph(self.mappings[-1], inGraph1, inGraph2)
        if save_all:
            GraphIO.write_graphML_file(graphAlignment, save_all)

        if (not G1.is_directed()) and (not G2.is_directed()):
            undir_gA = graphAlignment.to_undirected()
            return (undir_gA)
        else:
            return (graphAlignment)


    def disconnected_alig(self, G1, G2, node_score_list=False, edge_score_list=False, save_all=False,
                          forbidden=False):
        """Create a DirGraphAlign object and generate alignment of disconnected input graphs. The alignment is added
        to the self.mappings list.

        Parameters
        G1 and G2: input graphs
        node_score_list: contains scores for node label matches
        edge_score_list: contains scores for edge label matches
        save_all: saving option
        forbidden: dict with labels which are forbidden to match

        Returns
        Returns alignment graph
        """
        try:
            comp1=nx.number_connected_components(G1)
            comp2=nx.number_connected_components(G2)
        except:
            print("At least one input graph is disconnected and directed, but matching tool for disconnected graphs is"
                  "only implemented for undirected graphs.")
        else:
            G1_attr = G1.graph
            workgraph1 = nx.Graph(G1, **G1_attr)
            G2_attr = G2.graph
            workgraph2 = nx.Graph(G2, **G2_attr)
            orderlist = []
            mapping_list = []
            score_list = []

            for i in range(nx.number_connected_components(G1)):
                if workgraph1.order() == 0 or workgraph2.order() == 0:
                    break
                largest_cc_1 = max(nx.connected_components(workgraph1), key=len)
                cc_graph1 = workgraph1.subgraph(largest_cc_1)
                largest_cc_2 = max(nx.connected_components(workgraph2), key=len)
                cc_graph2 = workgraph2.subgraph(largest_cc_2)

                if cc_graph1.order() >= cc_graph2.order():
                    inGraph1 = cc_graph1.to_directed()
                    inGraph2 = cc_graph2.to_directed()
                    gorder = ["G1", "G2"]
                    orderlist.append(gorder)
                else:
                    inGraph1 = cc_graph2.to_directed()
                    inGraph2 = cc_graph1.to_directed()
                    gorder = ["G2", "G1"]
                    orderlist.append(gorder)

                if node_score_list:
                    disc_alignment = DirGraphAlign(inGraph1, inGraph2, node_score_list, forbid=forbidden)
                    all_alignments = disc_alignment.alig_lister()
                    if len(all_alignments) > 1:
                        selected_alignment = disc_alignment.process_results(all_alignments)
                        mapping_list.append(selected_alignment[1])
                        score_list.append(selected_alignment[0])
                    else:
                        mapping_list.append(all_alignments[0][1])
                        score_list.append(all_alignments[0][0])
                else:
                    disc_alignment = DirGraphAlign(inGraph1, inGraph2)
                    all_alignments = disc_alignment.alig_lister()
                    if len(all_alignments) > 1:
                        selected_alignment = disc_alignment.process_results(all_alignments)
                        mapping_list.append(selected_alignment)
                    else:
                        mapping_list.append(all_alignments[0])

                # remove nodes and neighbors of last matched components
                neighbors_to_delete1 = (set(cc_graph1.nodes))
                for node in cc_graph1:
                    neighbors_to_delete1.update(set(workgraph1.neighbors(node)))
                for delnode in neighbors_to_delete1:
                        workgraph1.remove_node(delnode)

                neighbors_to_delete2 = (set(cc_graph2.nodes))
                for node in cc_graph2:
                    neighbors_to_delete2.update(set(workgraph2.neighbors(node)))
                for delnode in neighbors_to_delete2:
                        workgraph2.remove_node(delnode)

            # mapping zusammenbauen
            joint_map = {}
            for j in range(len(mapping_list)):
                bla = orderlist[j][0]
                if orderlist[j][0] == "G2":
                    switchmap = {}
                    for k, v in mapping_list[j].items():
                        switchmap[v] = k
                    joint_map.update(switchmap)
                else:
                    joint_map.update(mapping_list[j])
            self.mappings.append(joint_map)

            # score addieren:
            if node_score_list:
                scoresum = 0
                for score in score_list:
                    scoresum += score
                self.scores.append(scoresum)
            pair_or_mult = inspect.currentframe().f_back.f_code.co_name
            if pair_or_mult == "vf2mult":
                self.add_matches_to_nodes(self.mappings[-1], G1, G2)
                self.add_matches_to_edges(self.mappings[-1], G1, G2)
            graphAlignment = self.build_matched_graph(self.mappings[-1], G1, G2)
            if save_all:
                GraphIO.write_graphML_file(graphAlignment, save_all)

            undir_gA = graphAlignment.to_undirected()
            print(self.mappings[-1])
            return (undir_gA)

    def add_matches_to_nodes(self, mapping, G1, G2):
        """Add matched nodes as a dict to the nodes of the input graphs.

        Parameters
        mapping: dict with matched nodes
        G1 and G2: matched graphs
        """
        for node in G1:
            if node in mapping:
                if 'Matches' in G1.nodes[node]:
                    alreadyMatched = G1.nodes[node]['Matches']
                    if 'Matches' in G2.nodes[mapping[node]]:
                        matched_matches = G2.nodes[mapping[node]]['Matches']
                        alreadyMatched.extend(matched_matches)
                    else:
                        alreadyMatched.append(mapping[node])
                    G1.nodes[node]['Matches'] = alreadyMatched

                else:
                    matched_nodes = []
                    if 'Matches' in G2.nodes[mapping[node]]:
                        matched_matches = G2.nodes[mapping[node]]['Matches']
                        matched_matches.insert(0, node)
                        matched_nodes.extend(matched_matches)
                    else:
                        matched_nodes.extend([node, mapping[node]])
                    G1.nodes[node]['Matches'] = matched_nodes
                if "None" in G1.nodes[node]['Matches']:
                    print(mapping)
                    print("Stop")
                    print("blabla")


                #for listing of node labels
                if 'Labelmatches' in G1.nodes[node]:
                    already_labelled = G1.nodes[node]['Labelmatches']
                    if 'Labelmatches' in G2.nodes[mapping[node]]:
                        matched_labmatches = G2.nodes[mapping[node]]['Labelmatches']
                        already_labelled.extend(matched_labmatches)
                    else:
                        already_labelled.append(G2.nodes[mapping[node]]['Label'])
                    G1.nodes[node]['Labelmatches'] = already_labelled
                else:
                    matched_labels = []
                    if 'Labelmatches' in G2.nodes[mapping[node]]:
                        matched_labmatches = G2.nodes[mapping[node]]['Labelmatches']
                        matched_labmatches.insert(0, G1.node.node['Label'])
                        matched_labels.extend(matched_labmatches)
                    else:
                        matched_node_label = G2.nodes[mapping[node]]['Label']
                        matched_labels.extend([G1.nodes[node]['Label'], matched_node_label])
                    G1.nodes[node]['Labelmatches'] = matched_labels

            else:
                g2matchtest = list(G2.nodes.data('Matches', default = [1]))
                gapnum = ['-'] * len(g2matchtest[0][1])

                if 'Matches' in G1.nodes[node]:
                    alreadyMatched = G1.nodes[node]['Matches']
                    alreadyMatched.extend(gapnum)
                    G1.nodes[node]['Matches'] = alreadyMatched
                else:
                    matched_nodes = []
                    matched_nodes.append(node)
                    matched_nodes.extend(gapnum)
                    G1.nodes[node]['Matches'] = matched_nodes


                #for listing of labels
                if 'Labelmatches' in G1.nodes[node]:
                    already_labelled = G1.nodes[node]['Labelmatches']
                    already_labelled.extend(gapnum)
                    G1.nodes[node]['Labelmatches'] = already_labelled
                else:
                    matched_labels = []
                    matched_labels.append(G1.nodes[node]['Label'])
                    matched_labels.extend(gapnum)
                    G1.nodes[node]['Labelmatches'] = matched_labels


        maplen = len(G1.nodes[node]['Matches'])

        for node2 in G2:
            if node2 not in mapping.values():
                maplist = ["-"] * (maplen-1)
                maplist.append(node2)
                G2.nodes[node2]['Matches'] = maplist
                labellist = ['-']*(maplen-1)
                labellist.append(G2.nodes[node2]['Label'])
                G2.nodes[node2]['Labelmatches'] = labellist


    def add_matches_to_edges(self, mapping, G1, G2):
        """Add matched edges to edges of input graphs.

        Parameters
        mapping: dict with matched nodes
        G1 and G2: matched graphs
        """
        edge_mapping = {}

        # for Graph1
        for edge in G1.edges(data=True):
            edge_ident = (str(edge[0]) + '-' + str(edge[1]))
            tail_node = edge[0]
            head_node = edge[1]
            edge_label = G1[tail_node][head_node]['Label']
            edge_construct = (tail_node, head_node, edge_label)
            # case: both nodes are in mapping
            if tail_node in mapping and head_node in mapping:
                matched_tail = mapping[tail_node]
                matched_head = mapping[head_node]
                matched_edge_label = G2[matched_tail][matched_head]['Label']
                if 'Matches' in G2[matched_tail][matched_head]:
                    matched_edge_matchings = G2[matched_tail][matched_head]['Matches']
                else:
                    matched_edge_matchings = None
                matched_edge_construct = (matched_tail, matched_head, matched_edge_label)
                edge_mapping[edge_ident] = matched_edge_construct
                if 'Matches' in G1[tail_node][head_node]:
                    already_matched = G1[tail_node][head_node]['Matches']
                    if matched_edge_matchings != None:
                        already_matched.extend(matched_edge_matchings)
                    else:
                        already_matched.append(matched_edge_construct)
                    G1[tail_node][head_node]['Matches'] = already_matched
                else:
                    matched = []
                    if matched_edge_matchings != None:
                        matched.extend([edge_construct, matched_edge_matchings])
                    else:
                        matched.extend([edge_construct, matched_edge_construct])
                    G1[tail_node][head_node]['Matches'] = matched

            # case: tail or head is not in mapping
            else:
                g2matchteste = list(G2.edges.data('Matches', default=[1]))
                gapnum = ['-'] * len(g2matchteste[0][2])

                if 'Matches' in G1[tail_node][head_node]:
                    already_matched = G1[tail_node][head_node]['Matches']
                    already_matched.extend(gapnum)
                    G1[tail_node][head_node]['Matches'] = already_matched
                else:
                    matched = [edge_construct]
                    matched.extend(gapnum)
                    G1[tail_node][head_node]['Matches'] = matched

            maplen = len(G1[tail_node][head_node]['Matches'])

            # for Graph2:
        for edge2 in G2.edges(data=True):
            edge2tuple = (edge2[0], edge2[1], edge2[2]['Label'])
            if edge2tuple not in edge_mapping.values():
                if 'Matches' in G2[edge2[0]][edge2[1]]:
                    g2matched = G2[edge2[0]][edge2[1]]['Matches']
                    maplist = ["-"] * (maplen - len(g2matched))
                    maplist.extend(g2matched)
                else:
                    maplist = ["-"] * (maplen - 1)
                    maplist.append(edge2tuple)
                G2[edge2[0]][edge2[1]]['Matches'] = maplist


    def build_matched_graph(self, mapping, G1, G2):
        """Build a new graph from mapping.

        Parameters
        mapping: a dict of matched nodes
        G1 and G2: graphs that were matched
        """
        ga = nx.DiGraph(Name='Alignment')
        for nodeg1 in G1.nodes():
            nodeattr = G1.nodes[nodeg1]
            ga.add_node(nodeg1, **nodeattr)
        ga.add_edges_from(G1.edges(data=True))

        for node in G2.nodes():
            if node not in mapping.values():
                nodeattr = G2.nodes[node]
                ga.add_node(node, **nodeattr)

        for edge in G2.edges():
            node1 = edge[0]
            node2 = edge[1]
            if (node1 not in mapping) and (node1 not in mapping.values()):
                # both nodes are not in mapping
                if (node2 not in mapping.values()) and (node2 not in mapping):
                    # ga.add_edge(node1, node2)
                    edgeattr = G2.edges[node1, node2]
                    ga.add_edge(edge[0], edge[1], **edgeattr)

                # node 2 is part of mapping, node 1 not
                else:
                    for k,v in mapping.items():
                        if v == node2:
                            edgeattr = G2.edges[edge]
                            ga.add_edge(node1, k, **edgeattr)

            elif (node2 not in mapping.values()) and (node2 not in mapping):
                for ke,val in mapping.items():
                    if val == node1:
                        edgeattr = G2.edges[edge]
                        ga.add_edge(node2, ke, **edgeattr)

        return (ga)

    def output_generator(self, mga, drawing_order, has_label, file):
        """Prepares and saves a graphical representation of the alignment.

        Parameters
        mga: alignment graph
        drawing_order: order in which input graphs were aligned with each other
        has_label: determines if node names or node labels are displayed
        file: where output is saved
        """
        pos = nx.spring_layout(mga)
        nodelist_complete_match = []
        nodelist_else = []
        edgelist_complete_match = []
        edgelist_else = []
        nodenames = {}
        nodelabels = {}

        for node in mga.nodes(data=True):
            node_name = node[0]
            nodenames[node_name] = node[1]['Matches']
            nodelabels[node_name] = node[1]['Labelmatches']
            if '-' in nodenames[node_name]:
                nodelist_else.append(node[0])
            else:
                nodelist_complete_match.append(node[0])

        for edge in mga.edges(data=True):
            if '-' in edge[2]['Matches']:
                edgelist_else.append(edge)
            else:
                edgelist_complete_match.append(edge)

        #here, lists can be generated from the mapping
        Mlist = self.convert_mga_to_list(mga)
        #print(Mlist)
        #EList = self.convert_mga_to_edgelist(mga)

        nx.draw_networkx_nodes(mga, pos, nodelist_complete_match, node_color='r', node_size=300, alpha=0.8)
        nx.draw_networkx_nodes(mga, pos, nodelist_else, node_color='g', node_size=300, alpha=0.8)
        nx.draw_networkx_edges(mga, pos, edgelist_complete_match, width=3, alpha=0.5, edge_color='r')
        nx.draw_networkx_edges(mga, pos, edgelist_else, width=3, alpha=0.5, edge_color='b')
        if has_label == "nolabel":
            nx.draw_networkx_labels(mga, pos, nodenames, font_size=12)
        else:
            nx.draw_networkx_labels(mga, pos, nodelabels, font_size=12)
        plt.title(drawing_order)
        plt.axis('off')
        plt.savefig(file, dpi=300)
        plt.clf()

    def report_score(self):
        if self.scores:
            return (self.scores[-1])
        else:
            pass

    def convert_mga_to_list(self, mga):
        MList = []
        for node in mga.nodes:
            matches = mga.nodes[node]['Matches']
            matchedlabels = mga.nodes[node]['Labelmatches']
            nodes_and_labels = []
            howmany = len(matches)
            for i in range(howmany):
                kombi = [matches[i], matchedlabels[i]]
                nodes_and_labels.append(kombi)
            MList.append(nodes_and_labels)
        return (MList)

    def convert_mga_to_edgelist(self, mga):
        edge_list = []
        for edge in mga.edges:
            matches = mga.edges[edge]['Matches']
            edge_list.append(matches)
        return (edge_list)

