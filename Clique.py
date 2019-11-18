import networkx as nx

class Clique():
    """Handles graph objects that are subgraphs to an input graph and correspond to a user predefined clique.
    Attributes
    ----------
    clique_dict : a dictionary storing graph objects that match the desired pre-defined clique, a subgraph found for each input graph.
    Methods
    ---------
    get_nodes
        gets nodes from modular graph product that correspond to predefined clique and adds them to a clique-graph object
    get_edges
        gets edges from modular graph product that correspond to predefined clique and adds them to a clique-graph object
    filter_candidates_mult
        removes candidates according to rules for BK algorithm, when starting with predefined clique graph object
    """


    def __init__(self):
        """
        Attributes
        --------
        clique_dict : a dictionary which stores graph-name as key and for value where clique was found in this graph in networkx format
        """
        self.clique_dict = {}


    def get_nodes(self,clique1,clique2,mgp):
        """Function gets nodes from modular graph product that correspond to predefined clique and adds them to an graph object to be used for alignment.

        Matches nodes of the previously found subgraphs. Looks for tupels in mod. product nodes that match.
        Parameters
        ---------
        clique1 : subgraph to graph1, corresponding to predefined anchor point
        clique1 : subgraph to graph2, corresponding to predefined anchor point
        mgp : modular graph product of graph1 and graph2
        """
        clique = nx.Graph()

        clique_nodes = []
        for node1 in clique1.nodes(data=True):
            for node2 in clique2.nodes(data=True):
                if node1[0][0] == node2[0][0]:
                    clique_nodes.append((node1[0][1], node2[0][1]))

        for a in mgp.nodes(data=True):
            if a[0] in clique_nodes:
                clique.add_node(a[0], Label=a[1]['Label'], Matches=a[1]['Matches'], Graph=a[1]['Graph'])

        return clique


    def get_edges(self,clique1,clique2,clique):
        """Function gets edges from modular graph product that correspond to to predefined clique and adds them to a graph object to be used for alignment.
        Parameters
        ---------
        clique1 : subgraph to graph1, corresponding to predefined anchor point
        clique1 : subgraph to graph2, corresponding to predefined anchor point
        clique : graph object, nodes were added, is to be used as anchor for BK after edges are added
        """

        for edges1 in clique1.edges(data=True):
            for edges2 in clique2.edges(data=True):
                if edges1[0][0] == edges2[0][0] and edges1[1][0] == edges2[1][0]:
                    label = [edges1[2]['Label'][1], edges2[2]['Label'][1]]
                    clique.add_edge((edges1[0][1], edges2[0][1]), (edges1[1][1], edges2[1][1]), Label=label)
                elif edges1[0][0] == edges2[1][0] and edges1[1][0] == edges2[0][0]:
                    label = [edges1[2]['Label'][1], edges2[2]['Label'][1]]
                    clique.add_edge((edges1[0][1], edges2[1][1]), (edges1[1][1], edges2[0][1]), Label=label)

        return clique


    def filter_candidates_mult(self, clique1,clique2,clique, mgp):
        """Removes candidates for using in BK algorithm when using pre-defined clique.

        Removes nodes already in clique object, nodes that already have a match in clique object and nodes that are not neighbours to clique object
        Parameters
        ----------
        clique1 : subgraph to graph1, corresponding to predefined anchor point
        clique1 : subgraph to graph2, corresponding to predefined anchor point
        clique : graph object, anchor for BK in alignment of graph1 and graph2
        mgp : modular graph product of graph1 and graph2
        """

        candidates = list(mgp.nodes)

        for x in clique.nodes:
            candidates.remove(x)

        clique1_nodes = [x[1] for x in clique1.nodes]
        clique2_nodes = [x[1] for x in clique2.nodes]

        for thing in candidates[:]:
            keep_cand = False
            if thing[1] in clique2_nodes:
                try:
                    candidates.remove(thing)
                except ValueError:
                    pass
            elif thing[0] in clique1_nodes:
                try:
                    candidates.remove(thing)
                except ValueError:
                    pass

            else:
                for clique_node in clique.nodes:
                    if thing in mgp.neighbors(clique_node):
                        keep_cand = True

                if not keep_cand:
                    try:
                        candidates.remove(thing)
                    except ValueError:
                        pass

        return candidates
