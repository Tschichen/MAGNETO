import networkx as nx
import json


class GraphIO():
    """A class for parsing and writing graphs and matches from and in different file formats.

        Methods
        -------
        parse_json_pubchem
            takes a path to a json formatted 2D file of molecules from the pubchem database and returns a networkX Graph object
        parse_graphML
            takes a path to a graph in graphML format and returns a networkX DiGraph object
        parse_file
            takes a path to a graph in graph format and returns a networkX DiGraph object
        write_graph_file
            takes a path to the location of the output file and a networkX Graph or DiGraph object and writes it into a file in graph file format
        write_graphML_file
            takes a path to a graphML output file and a networkX Graph, networkX Digraph or matching object and writes it in graphML format
        write_matches
            takes a path to an output file and generates a MATCHING file with the nodes and edges of the final matching object.
    """

    def parse_json_pubchem(path):
        """A function that takes a path to a json formatted 2D file of molecules from the pubchem database and returns a networkX Graph object.

            The labeling od nodes depends on the atoms of which the molecule consists and corresponds to the numbering in the periodic table.
            The edge labeling depends on the bond between the single atoms or groups and is exactly broken down in the function return_bond.
        """
        try:
            file = open(path, "r")
        except:
            print("File not found")
            exit()

        input = json.load(file)

        input_id = input['PC_Compounds'][0]['id']['id']['cid']
        name = input['PC_Compounds'][0]['props'][9]['value']['sval']
        nodes = input['PC_Compounds'][0]['atoms']['aid']
        node_label = input['PC_Compounds'][0]['atoms']['element']
        node_1 = input['PC_Compounds'][0]['bonds']['aid1']
        node_2 = input['PC_Compounds'][0]['bonds']['aid2']
        bond = input['PC_Compounds'][0]['bonds']['order']

        nodes_with_labels = []
        for i in range(len(nodes)):
            node_name = str(input_id) + "_" + str(nodes[i])
            node = (node_name, {'Label': str(node_label[i]), 'Graph': input_id})
            nodes_with_labels.append(node)

        edges_with_labels = []
        for i in range(len(node_1)):
            name_node_1 = str(input_id) + "_" + str(node_1[i])
            name_node_2 = str(input_id) + "_" + str(node_2[i])
            edge = (name_node_1, name_node_2, {'Label': GraphIO.return_bond(bond[i]), 'Graph': input_id})
            edges_with_labels.append(edge)

        G = nx.Graph(Name=str(input_id))
        G.add_nodes_from(nodes_with_labels)
        G.add_edges_from(edges_with_labels)

        file.close()

        return G

    def return_bond(int):
        """Help function that takes an integer and returns the type of chemical bond.

            Args and Return
            ---------------
                1 : single bond
                2 : double bond
                3 : triple bond
                4 : quadruple bond
                5 : dative bond
                6 : complex bond
                7 : ionic bond
        """
        if int == 1:
            return "single"
        elif int == 2:
            return "double"
        elif int == 3:
            return "triple"
        elif int == 4:
            return "quadruple"
        elif int == 5:
            return "dative"
        elif int == 6:
            return "complex"
        elif int == 7:
            return "ionic"
        else:
            return None

    def parse_graphML_file(path):
        """Function that takes a path to a graph in graphML format and returns a networkX DiGraph object.

            All available informations such as knots and edge labels as well as whether a directed or undirected graph is processed.
        """
        try:
            file = open(path, "r")
        except:
            print("File not found")
            exit()
        G = nx.read_graphml(file)
        file.close()
        for _, d in G.nodes(True):
            for k in d:
                d[k] = json.loads(d[k])
        for _, _, d in G.edges(data=True):
            for k in d:
                d[k] = json.loads(d[k])
        return G

    def parse_file(file_loc):
        """This function takes a path to a graph in graph format and returns a networkX DiGraph object.

            All available informations such as knots and edge labels as well as whether a directed or undirected graph is processed.
        """
        false_var = ["false", "0"]
        true_var = ["true", "1"]
        g = nx.DiGraph()
        i = 0
        anz_attributes = 5
        name = file_loc.split('\\')[-1].split('.')[0]
        try:
            f = open(file_loc)
            fl = f.readlines()
        except:
            print("File not found: " + file_loc)
            exit()

        # find starting point of the data in file, in case of comments at the top of the document
        try:
            while 'start' not in locals():
                x = fl[i]
                x = x.lower()
                # number of Nodes
                if x.startswith("#nodes"):
                    anznodes = int(x.split(";")[-1].strip())
                    start = fl.index(x)
                i += 1

        except:
            print("Error: Couldn't find #nodes.")
            exit()

        try:
            for x in fl[start:start + anz_attributes]:
                x = x.lower()

                # number of edges
                if "#edges" in x:
                    anzedges = int(x.split(";")[-1].strip())

                # nodes labelled or not
                if "nodes labelled" in x:
                    if x.split(";")[-1].lower().strip() in true_var:
                        n_label = True
                    elif x.split(";")[-1].lower().strip() in false_var:
                        n_label = False
                    else:
                        print("Error: Please define if nodes are labelled.")
                        exit()

                # edges labelled or not
                if "edges labelled" in x:
                    if x.split(";")[-1].lower().strip() in true_var:
                        e_label = True
                    elif x.split(";")[-1].lower().strip() in false_var:
                        e_label = False
                    else:
                        print("Error: Please define if edges are labelled.")
                        exit()

                # directed or undirected
                if "directed graph" in x and 'graphdir' not in locals():
                    # check if graph is undirected
                    if x.split(";")[-1].lower().strip() in true_var:
                        graphdir = True
                        has_edges = True
                        end_of_nodes = start + 6 + anznodes
                    # check if graph has no defined edges
                    elif len(fl[start + 5:]) == anznodes + 1 and x.split(";")[-1].lower().strip() in false_var:
                        graphdir = False
                        has_edges = False
                        end_of_nodes = len(fl[:])
                        g = g.to_undirected()
                    # Check if graph has defined edges
                    elif len(fl[start + 5:]) == anznodes + anzedges + 2 and x.split(";")[
                        -1].lower().strip() in false_var:
                        graphdir = False
                        has_edges = True
                        end_of_nodes = start + 6 + anznodes
                        g = g.to_undirected()
                    else:
                        print("Please check var 'Graph directed' and edges " + name)
                        exit()

        except:
            print("Error: Please check File for correct structure: " + name)
            exit()

        # check if all attributes were gathered
        if 'anzedges' not in locals():
            print("Error: Couldn't find '#edges'.")
            exit()
        elif 'graphdir' not in locals():
            print("Error: Couldn't find 'Graph directed'.")
            exit()
        elif 'e_label' not in locals():
            print("Error: Couldn't find 'Edges labelled'.")
            exit()
        elif 'n_label' not in locals():
            print("Error: Couldn't find 'Nodes labelled'.")
            exit()

        # Nodes
        try:
            # check for spacebar
            if fl[start + 5] == "\n" and (has_edges == False or fl[start + 6 + anznodes] == "\n"):
                # add nodes to graph
                for x in fl[start + 6:end_of_nodes]:
                    node = x.strip().split(";")

                    node_name = name + "_" + node[0]
                    # unlabelled
                    if n_label == False:
                        g.add_node(node_name, Label=None, Graph=name)
                    # labelled
                    elif n_label == True:
                        g.add_node(node_name, Label=str(node[1].strip()), Graph=name)
                    else:
                        print("Error: '#nodes' is not correct")
                        exit()

        except:
            print("Error: Nodes couldn't parse. Please check file for correct structure.")
            exit()

        # edges
        try:
            if has_edges == True and fl[start + 6 + anznodes] == "\n" and fl[start + 6 + anznodes + anzedges] == fl[-1]:
                for x in fl[start + 7 + anznodes:]:
                    edge = x.strip().split(";")
                    edge_name1 = name + "_" + edge[0]
                    edge_name2 = name + "_" + edge[1]
                    try:
                        # undirected
                        if graphdir == False:
                            if e_label == True:
                                g.add_edge(edge_name1, edge_name2, Label=edge[2], Graph=name)
                            else:
                                g.add_edge(edge_name1, edge_name2, Label=0, Graph=name)

                        # directed
                        elif graphdir == True:
                            if e_label == True:
                                g.add_edge(edge_name1, edge_name2, Label=edge[2], Graph=name)
                            else:
                                g.add_edge(edge_name1, edge_name2, Label=0, Graph=name)


                    except:
                        print("Error: Please check line: " + str(edge))
                        exit()


        except:
            print("Error: Edges couldn't parse. Please check file for correct structure.")
            exit()

        try:
            if has_edges == True:
                for x in fl[start + 6:start + 6 + anznodes]:
                    node = x.split(";")
                    if len(node) > 1:
                        node[-1] = node[-1].strip()
                        h = g.to_undirected()
                        node_name = name + "_" + node[0]
                        neig = list(h.adj[node_name])

                        if set(neig).issubset(node[2:]) != True and set(node[2:]).issubset(neig) != True:
                            print("Error: Edges are not equal with neighbors")
                            exit()
        except:
            print("Error: Couldn't compare edges with neighbors")
            exit()

        g.graph['Name'] = file_loc.split('\\')[-1].split('.')[0]
        return g

    def write_graph_file(graph, path):
        """Function that takes a path to the location of the output file and a networkX Graph or DiGraph object and writes it into a file in graph file format.
        """
        file = open(path, "w+")
        graph_file = ""
        graph_file = "#" + graph.graph['name']
        graph_file += "#" + str(len(graph.nodes())) + "\n"
        graph_file += "#" + str(len(graph.edges())) + "\n" + "True\n" + "True\n" + "False\n\n"
        for node in graph.nodes():
            if 'Label' in graph.nodes[node]:
                graph_file += str(node) + "; " + graph.nodes[node]['Label'] + "\n"
            else:
                graph_file += str(node) + ";\n"
        graph_file += "\n"
        for edge in graph.edges():
            if 'Label' in graph.edges[edge]:
                graph_file += str(edge[0]) + "; " + str(edge[1]) + "; " + graph.edges[edge]['Label'] + "\n"
            else:
                graph_file += str(edge[0]) + "; " + str(edge[1]) + ";\n"
        file.write(graph_file)
        file.close()

    def write_graphML_file(graph, path):
        """Function that takes a path to the location of the output file and a networkX Graph or DiGraph object and writes it into a file in graphML file format.
        """
        try:
            outgraph = graph.copy()
        except AttributeError:
            outgraph = graph.network.copy()
        for _, d in outgraph.nodes(True):
            for k in d:
                d[k] = json.dumps(d[k])
        for _, _, d in outgraph.edges(data=True):
            for k in d:
                d[k] = json.dumps(d[k])
        nx.write_graphml(outgraph, path)


    def write_matches(matching, path):
        """Function that takes a path to an output file and generates a MATCHING file with the nodes and edges of the
        final matching object.
        """

        matching_obj = "#Node matching\n"

        for node in matching.nodes(data=True):
            matching_obj += 'Node;' + node[0] + '\nLabel;' + str(node[1]['Label']) + '\nGraph;' + str(node[1][
                                                                                                          'Graph']) + '\nMatches;'
            for match in node[1]['Matches']:
                if match != '-':
                    triple = "(" + str(match) + ")"
                else:
                    triple = 'None'
                matching_obj += triple + ";"
            matching_obj += "\n\n"

        matching_obj += "#Edges\n"
        for edge in matching.edges(data=True):
            matching_obj += 'Edge;' + edge[0] + ";" + edge[1] + '\nLabel;' + str(edge[2]['Label']) + '\nGraph;' + \
                            str(edge[2]['Graph']) + '\nMatches;'
            for match in edge[2]['Matches']:
                if match == '-'or match == None:
                    quadruple = 'None'
                else:
                    quadruple = "(" + str(match[0]) + "," + str(match[1]) + "," + str(match[2]) + ")"
                matching_obj += quadruple + ";"
            matching_obj += "\n\n"

        with open(path, 'w') as save_matches:
            save_matches.write(matching_obj)


