#! /usr/bin/env python3

import sys
import getopt
import glob
import networkx as nx
import time
import os
# moduls by our group
from GraphIO import GraphIO
import guide_tree as gt
import matching as mt
import VF2_pair_multi_flash as vf
import keep
import Clique
import graph_gen as gg


def main(argv):
    vf2 = True  # use vf2 algorithm
    pre_clique = False  # using a predefined clique
    windows = False  # catch exceptions if not used under windows
    scoring_for = "both"
    save_graphs = "end"
    forbidden = False

    '''Get all User Options'''

    # Catch some Exceptions for bad user input.
    try:
        opts, args = getopt.getopt(argv, "hi:o:c:n:rts:bqg:m:dvfHp:",
                                   ["help", "bronkerbosch", "score=", "labelling=", "newick=", "randomt", "showt",
                                    "savetn", "drawm", "graphgen", "forbidden"])

        all_opts = [x[:][0] for x in opts]

    except getopt.GetoptError:
        print("Argument Error. Please check --help or -h")
        exit()

    if "--graphgen" in all_opts:
        print("Welcome to the graph generator!")
        gg.main()
        exit()

    # Iterate over user options.
    for opt, arg in opts:
        if opt == '-h' or opt == '--help':
            print("----Multiple Graph Alignment Tool----")
            print("aligns multiple graphs in a progressive way using Bron Kerbosch or VF2 algorithm")
            print()
            print("GRAPH OPTIONS")
            print(
                "-i <input path>\t\t\tfilepath from which all the graph files are read. Files can be .graph, .json or .graphml. MANDATORY OPTION!")
            print("-o <name>\t\t\tdefine name for the current aligment run. MANDATORY OPTION!")
            print(
                "-c <input file>\t\t\tinput a graph file as an anchor for alignment. This graph needs to be common to all graphs that are to be aligned and acts as a basis to the alignment")
            print("-H\t\t\t\tremoves hydrogen molecules from input graphs. Or in general nodes labelled '1'")
            print()
            print("ALIGNMENT OPTIONS")
            print(
                "-b or --bronkerbosch\t\tswitches to multiple alignment via the Bron-Kerbosch algorithm. Default: VF2")
            print("-v\t\t\t\twhen Bron Kerbosch is activated, uses a random pivot element (Note: not usable with clique and scoring options)")
            print(
                "-s or --score <input file>\tinput a list of scores for matching node labels, edge labels or both. Default: Scoring based on the size of the largest common subgraph found between two graphs")
            print(
                "-f or --forbidden\t\twhile inputting a scoring list, score via size of largest common subgraph after excluding forbidden label matches.")
            print(
                "-g or --labelling <string>\tdefine what labels should be used for scoring: node labels, edge labels or both. Default: both")
            print()
            print("GUIDE TREE OPTIONS")
            print(
                "\tDefault: guide tree is generated from scores or size of the largest common subgraph computed via pairwise alignment")
            print("-n or --newick <file>\t\tinput a guide tree file in newick format")
            print("-r or --randomt\t\t\tadditionally creates a random guide tree and saves it in newick format")
            print("-t or --showt\t\t\tsaves the generated guide tree as .png")
            print("-q or --savetn\t\t\tsaves generated guide tree in newick format")
            print()
            print("RESULT OPTIONS")
            print(
                "-m <string>\t\t\tsaves multiple alignment as graphml file. Chose between \"end\" and \"all\". end only saves total matching. all saves all the matchings in between as well.")
            print("-d or --drawm\t\t\tsaves multiple alignment visualization as .png")
            print()
            print("ADDITIONAL FEATURE:")
            print(
                "--graphgen\t\t\ta tool, that makes random graphs with properties that are set through user input")
            print()
            print("Note: Clique option only available for BK. Disconnected directed graphs can only be assessed via BK algorithm. Input graphs must be all directed or undirected.")
            exit()

        if '-i' not in all_opts and '-p' not in all_opts:
            print("Please specify the location where your graph files are stored using -i <location>.")
            exit()

        if '-o' not in all_opts:
            print("Please name your current alignment job using -o <name>.")
            exit()

        if opt == '-i':
            if sys.platform == "win32":
                windows = True
            if arg.endswith("\\") or arg.endswith("/"):
                graph_loc = arg
            else:
                if windows:
                    graph_loc = arg + "\\"
                else:
                    graph_loc = arg + "/"

        if opt == '-p':
            if sys.platform == "win32":
                windows = True

            graphs = arg.split('%')
            for i in range(len(graphs)):
                singlegraph = []
                singlegraph.append(graphs[i])
            if not windows:
                graph_loc_list = arg[0].split('/')
                graph_loc_list.pop()
                graph_loc = '/'.join(graph_loc_list) + '/'
            else:
                graph_loc_list = arg[0].split('\\')
                graph_loc_list.pop()
                graph_loc = '\\'.join(graph_loc_list) + '\\'


        if opt == '-o':
            name = arg
            dir_path = os.getcwd()
            if windows:
                dir_path += "\\"
            else:
                dir_path += "/"
            path = dir_path + name
            try:
                os.mkdir(path)
            except:
                pass

            if windows:
                path += '\\'
            else:
                path += '/'

        if opt == '-b' or opt == '--bronkerbosch':
            vf2 = False

        if opt == '-c':
            if vf2:
                print("There is no clique option when using VF2, please switch to Bron kerbosch using '-b'")
                exit()
            pre_clique = True
            if arg.endswith('graph'):
                clique = GraphIO.parse_file(arg,'clique_temp')
            elif arg.endswith('json'):
                clique = GraphIO.parse_json_pubchem(arg)
            elif arg.endswith('graphml'):
                clique = GraphIO.parse_graphML_file(arg)
            else:
                print("Could not parse Clique file. Please check file format")
                exit()
            if not clique:
                exit()

        if opt == '-s' or opt == '--score':
            try:
                user_list_scores = keep.parse_keep_input(arg)  # parse list of scores
                if '-f' in all_opts or '--forbidden' in all_opts:
                    forbidden = True
                else:
                    forbidden = False
            except:
                print("Error with Score List Input File!")
                exit()

        if opt == '-n' or opt == '--newick':
            newick_path = arg

        if opt == '-g' or opt == '--labelling':
            if arg in ["nodes", "edges", "both"]:
                scoring_for = arg
            else:
                print("Error! input for scoring mode can only be \"nodes\", \"edges\" or \"both\".")
                exit()

        if opt == '-m':
            if arg != "end" and arg != "all":
                print("Please chose between \"end\" or \"all\", when using -m option.")
                exit()
            save_graphs = arg

    if '-s' not in all_opts and '--score' not in all_opts:
        user_list_scores = []

    '''parse all Graph files in Location'''

    # Get all graph file paths.
    if not '-p' in all_opts:
        try:
            graphs = glob.glob(graph_loc + "*.graph")
            js_graphs = glob.glob(graph_loc + "*.json")
            ml_graphs = glob.glob(graph_loc + "*.graphml")
        except:
            print("Input Error, check command line!")
            exit()
    else:
        graphs = []
        js_graphs = []
        ml_graphs = []
        for i in range(len(singlegraph)):
            if singlegraph[i].endswith('.graph'):
                graphs.append(singlegraph[i])
            elif singlegraph[i].endswith('.json'):
                js_graphs.append(singlegraph[i])
            elif singlegraph[i].endswith('.graphml'):
                ml_graphs.append(singlegraph[i])

    graphen = {}
    all_graphs = []

    # Parse all .graph files.
    for i in graphs:
        if windows:
            graph_name = i.split("\\")[-1].split(".")[0]
        else:
            graph_name = i.split("/")[-1].split(".")[0]
        if graph_name not in all_graphs:
            g = GraphIO.parse_file(i, graph_name)
            graphen[graph_name] = g
            all_graphs.append(graph_name)
        else:
            print("Multiple graphs are named: " + graph_name + " Input was skipped.")

    # Parse all .json files.
    for i in js_graphs:
        if windows:
            graph_name = i.split("\\")[-1].split(".")[0]
        else:
            graph_name = i.split("/")[-1].split(".")[0]
        if graph_name not in all_graphs:
            g = GraphIO.parse_json_pubchem(i)
            graphen[graph_name] = g
            all_graphs.append(graph_name)
        else:
            print("Multiple graphs are named: " + graph_name + " Input was skipped.")

    # Parse all .graphml files.
    for i in ml_graphs:
        if windows:
            graph_name = i.split("\\")[-1].split(".")[0]
        else:
            graph_name = i.split("/")[-1].split(".")[0]
        if graph_name not in all_graphs:
            g = GraphIO.parse_graphML_file(i)
            graphen[graph_name] = g
            all_graphs.append(graph_name)
        else:
            print("Multiple graphs are named: " + graph_name + " Input was skipped.")

    if len(all_graphs) < 2:
        print("Too few graphs in directory! " + str(len(all_graphs)) + " graphs in " + graph_loc)
        exit()


    #Check if graphs are all directed or undirected.
    directed_graphs = []
    for i in all_graphs:
        if graphen[i].is_directed:
            directed_graphs.append(i)

    if (len(directed_graphs) != 0) and (len(directed_graphs) != len(all_graphs)):
        print("Graphs are inconsistent in directed and undirected attribute. Please check input and only use graphs that are undirected or directed!")
        print("Those graphs are directed:")
        print(directed_graphs)
        exit()


    # Check if graphs are consistent in labelling.
    node_labels = []
    for i in all_graphs:
        for node in graphen[i].nodes(data=True):
            node_labels.append(node[1]['Label'])
            break
    if None not in node_labels:
        nodes_labelled = True
    else:
        for i in node_labels:
            if i != None:
                print("Warning! Some of your graphs have labeled nodes, some don't.")
                break
        nodes_labelled = False

    edge_labels = []
    for i in all_graphs:
        for edge in graphen[i].edges(data=True):
            edge_labels.append(edge[2]['Label'])
            break
    if int(0) not in edge_labels:
        edge_labels = True
    else:
        for i in edge_labels:
            if i != 0:
                print("Warning! Some of your graphs have labeled edges, some don't.")
                break
        edges_labelled = False
        
    if "-H" in all_opts:
        a = nx.get_node_attributes(graphen[all_graphs[0]],"Label")
        if list(a.items())[0][1] != None:
            for graph in all_graphs:
                h_nodes = []
                for node in graphen[graph].nodes(data=True):
                    if node[1]['Label'] == '1':
                        h_nodes.append(node[0])
                graphen[graph].remove_nodes_from(h_nodes)
            print('Filtered Hydrogen')
        else:
            print("Graphs are not labelled, can't filter for Hydrogen")        

    '''pairwise Alignment and Scoring'''

    print("Multiple Graph Alignment: " + name)
    print("Alignment: " + ", ".join(all_graphs))



    # If a clique is pre-defined all graphs need to be checked if clique is in them. Else return graph that does not have clique.
    if pre_clique == True:
        print("\tChecking if pre-defined clique is present in all graphs.")
        all_clique = Clique.Clique()
        for j in all_graphs:
            m = mt.Matching(clique)
            k = mt.Matching(graphen[j])
            mgp = m.modular_product(k)
            if mgp == -1:
                print("Error!! Your clique does not match in directed/undirected attribut! Check graph " + j + "!")
                exit()
            if user_list_scores:
                mgp = keep.modular_product_valid_matchings(mgp, user_list_scores, scoring_for)
            clique_list = []
            clique_empty = nx.Graph()
            checked = []
            m.BK(mgp, clique_empty, list(mgp.nodes), checked, clique_list)
            largest_clique = m.find_largest_clique(clique_list)

            score = keep.similarity_score_to_distance(clique, graphen[j], len(largest_clique.nodes))
            if score > 0:
                print("Clique Error. Clique not found in graph " + j)
                exit()
            else:
                all_clique.clique_dict[j] = largest_clique
        print("\tdone, success!")

    '''Generating Guide Tree after PWA'''

    if '-n' not in all_opts:
        pre_score_list = [[]]
        for i in range(3):
            pre_score_list.append([])
        pair_align_scores = []
        print("\tEvaluating scores for Guide Tree...")

        if vf2:
            time_vf21 = time.time()
        else:
            time_bk1 = time.time()

        # Pairwise alignment of all graphs.
        if not vf2:
            for i in range(len(all_graphs)):
                matcher = mt.Matching(graphen[all_graphs[i]])
                for j in range(i + 1, len(all_graphs)):
                # Evaluate pairwise scores via Bron Kerbosch.
                    obj = mt.Matching(graphen[all_graphs[j]])
                    h = matcher.modular_product(obj)
                    if h == -1:
                        print("Error!! Your graphs do not match in directed/undirected attribut! Check graphs " +
                              all_graphs[i] + " and " + all_graphs[j] + "!")
                        exit()
                    if user_list_scores:
                        h = keep.modular_product_valid_matchings(h, user_list_scores, scoring_for)
                    clique_list = []
                    checked = []
                    clique = nx.Graph()
                    candidates = list(h.nodes)

                    # Construction of clique for BK algorithm.
                    if pre_clique:
                        clique1 = all_clique.clique_dict[all_graphs[i]]
                        clique2 = all_clique.clique_dict[all_graphs[j]]

                        clique = all_clique.get_nodes(clique1, clique2, h)
                        candidates = all_clique.filter_candidates_mult(clique1, clique2, clique, h)

                        clique = all_clique.get_edges(clique1, clique2, clique)

                        clique_list.append(clique)
                        matcher.BK(h, clique, candidates, checked, clique_list)

                    # Random or custom pivot element choice.
                    elif '-v' not in all_opts:
                        matcher.BK(h, clique, candidates, checked, clique_list)
                    elif '-v' in all_opts:
                        matcher.BK_random(h, clique, candidates, checked, clique_list)

                    largest_clique = matcher.find_largest_clique(clique_list)

                    if not user_list_scores or forbidden:
                        score = keep.similarity_score_to_distance(graphen[all_graphs[i]], graphen[all_graphs[j]],
                                                                  len(largest_clique.nodes))
                        print("Score: " + str(score) + " -  " + all_graphs[i] + " & " + all_graphs[j])
                        print("Best Scoring Clique:")
                        print(largest_clique.nodes)
                        pair_align_scores.append([all_graphs[i], all_graphs[j], score])
                    else:
                        pre_score = keep.score_from_BK(clique_list, user_list_scores, scoring_for)
                        pre_score_list[0].append(pre_score[0])
                        pre_score_list[1].append(all_graphs[i])
                        pre_score_list[2].append(all_graphs[j])
                        pre_score_list[3].append(largest_clique.nodes)
                        print(all_graphs[i] + " & " + all_graphs[j] + " - done!")


        # Pairwise alignment via VF2.
        else:
            for i in range(len(all_graphs)):
                for j in range(i + 1, len(all_graphs)):
                    graph_list = [graphen[all_graphs[i]], graphen[all_graphs[j]]]
                    if user_list_scores and scoring_for == "nodes":
                        aligner = vf.VF2_GraphAligner(graph_list, None, user_list_scores[0])
                    elif user_list_scores and scoring_for == "edges":
                        aligner = vf.VF2_GraphAligner(graph_list, None, False, user_list_scores[1])
                    elif user_list_scores and scoring_for == "both":
                        aligner = vf.VF2_GraphAligner(graph_list, None, user_list_scores[0], user_list_scores[1])
                    else:
                        aligner = vf.VF2_GraphAligner(graph_list)

                    m = aligner.vf2_pga(forbidden)
                    if not user_list_scores or forbidden:
                        largest_clique = aligner.length
                        score = keep.similarity_score_to_distance(graphen[all_graphs[i]], graphen[all_graphs[j]],
                                                                  largest_clique)
                        print("Score: " + str(score) + " - " + all_graphs[i] + " & " + all_graphs[j])
                        print("Matching:")
                        print(aligner.mappings[-1])
                        pair_align_scores.append([all_graphs[i], all_graphs[j], score])

                    else:
                        pre_score_list[0].append(aligner.scores[0])
                        pre_score_list[1].append(all_graphs[i])
                        pre_score_list[2].append(all_graphs[j])
                        pre_score_list[3].append(aligner.mappings[-1])
                        print(all_graphs[i] + " & " + all_graphs[j] + " - done!")

        # Compute upgma scores after all pairwise scores were created when using custom score lists.
        if user_list_scores and not vf2 and not forbidden:
            upgma_scores = keep.upgma_score_for_list(pre_score_list[0])
            for i in range(len(upgma_scores)):
                print("Score: " + str(upgma_scores[i]) + " -  " + pre_score_list[1][i] + " & " + pre_score_list[2][i])
                print("Best Scoring Clique:")
                print(pre_score_list[3][i])
                pair_align_scores.append([pre_score_list[1][i], pre_score_list[2][i], upgma_scores[i]])

        elif user_list_scores and vf2 and not forbidden:
            upgma_scores = keep.upgma_score_for_list(pre_score_list[0])
            for i in range(len(upgma_scores)):
                print("Score: " + str(upgma_scores[i]) + " -  " + pre_score_list[1][i] + " & " + pre_score_list[2][i])
                print("Matching:")
                print(pre_score_list[3][i])
                pair_align_scores.append([pre_score_list[1][i], pre_score_list[2][i], upgma_scores[i]])

        if vf2:
            time_vf2_paar = time.time() - time_vf21
        else:
            time_bk_paar = time.time() - time_bk1

        if True:
            # Generate guide tree from scores.
            guide_tree = gt.Guide_tree_Generator.tree_from_scores(pair_align_scores)

    else:
        # Read guide tree from newick file.
        try:
            time_vf2_paar = 0
            time_bk_paar = 0
            guide_tree = gt.Guide_tree_Generator.tree_from_newick(newick_path)
        except:
            print("Guide Tree Error! A Problem with the newick input file occurred.")
            exit()

    # Write a random guide tree to a newick file.
    if ('-r' in all_opts or '--randomt' in all_opts) and not '-n' in all_opts:
        random_guide_tree = gt.Guide_tree_Generator.tree_from_random(pair_align_scores)
        random_guide_tree.alignment_order()
        file = path + name + "_random.nwtree"

        try:
            with open(file, 'w') as save:
                save.write(random_guide_tree.newick_output())
        except:
            print("Could not write to " + file + ".")

    ao = guide_tree.alignment_order()
    alignment_order = guide_tree.alignment_order_drawing()

    # Brings order for alignment and order for visualization-function to same structure.
    drawing_order = []
    already_seen = set()
    for i in ao:
        if i[0] not in already_seen and i[0] not in drawing_order:
            drawing_order.append(i[0])
            already_seen.add(i[2])
        else:
            already_seen.add(i[0])
        if i[1] not in already_seen and i[1] not in drawing_order:
            drawing_order.append(i[1])
            already_seen.add(i[2])
        else:
            already_seen.add(i[0])
            already_seen.add(i[1])
            already_seen.add(i[2])

    # Switches some problematic leaves in guide tree for optimized alignment order.
    # Note: Does not change the actual order, just switches some neighbouring leaves around.
    new_ao = []
    conflicts = []
    solutions = []
    walkalong = set()
    new_ao.append(ao[0])
    for i in ao[0]:
        walkalong.add(i)
    for i in range(1, len(ao)):
        if ao[i][0] not in conflicts and ao[i][1] not in conflicts:
            if ao[i][0] in walkalong:
                new_ao.append(ao[i])
                walkalong.add(ao[i][1])
                walkalong.add(ao[i][2])
            elif ao[i][1] in walkalong:
                conflicts.append(ao[i][2])
                new_name = ao[i][1] + ao[i][0]
                solutions.append(new_name)
                triple = (ao[i][1], ao[i][0], new_name)
                new_ao.append(triple)
                walkalong.add(ao[i][0])
                walkalong.add(ao[i][2])
            else:
                new_ao.append(ao[i])
                walkalong.add(ao[i][0])
                walkalong.add(ao[i][1])
                walkalong.add(ao[i][2])
        elif ao[i][0] in conflicts and ao[i][1] not in conflicts:
            get_index = conflicts.index(ao[i][0])
            new_name = solutions[get_index] + ao[i][1]
            triple = (solutions[get_index], ao[i][1], new_name)
            new_ao.append(triple)
            conflicts.append(ao[i][2])
            solutions.append(new_name)
            walkalong.add(ao[i][1])
            walkalong.add(ao[i][2])
        elif ao[i][0] not in conflicts and ao[i][1] in conflicts:
            get_index = conflicts.index(ao[i][1])
            new_name = solutions[get_index] + ao[i][0]
            triple = (solutions[get_index], ao[i][0], new_name)
            new_ao.append(triple)
            conflicts.append(ao[i][2])
            solutions.append(new_name)
            walkalong.add(ao[i][0])
            walkalong.add(ao[i][2])
        else:
            get_index1 = conflicts.index(ao[i][0])
            get_index2 = conflicts.index(ao[i][1])
            new_name = solutions[get_index1] + solutions[get_index2]
            triple = (solutions[get_index1], solutions[get_index2], new_name)
            new_ao.append(triple)
            conflicts.append(ao[i][2])
            solutions.append(new_name)
            walkalong.add(ao[i][0])
            walkalong.add(ao[i][1])
            walkalong.add(ao[i][2])

    # Save guide tree as a PNG file.
    if '-t' in all_opts or '--showt' in all_opts:
        print("Saved Guide Tree as PNG")
        png_file = path + name + "_tree.png"
        guide_tree.tree_drawing(png_file)

    # Save guide tree as newick string in file.
    if '-q' in all_opts or '--savetn' in all_opts:
        file = path + name + ".nwtree"
        try:
            with open(file, 'w') as save:
                save.write(guide_tree.newick_output())
            print("Guide tree \"" + name + "\" saved as newick format")
        except:
            print("Could not write to " + name + ".nwtree file!")

    '''Multiple Graph Alignment'''

    print("\tGenerating MGA...")

    ml_path = path + name + "_fin_match.graphml"

    # Multiple Graph Alignment via Bron Kerbosch.
    if not vf2:
        time_bk2 = time.time()
        try:
            matching_dict = {}
            for i in range(0, len(drawing_order)):
                matching_dict[drawing_order[i]] = mt.Matching(graphen[drawing_order[i]])

            m = matching_dict[drawing_order[0]]
            if '-v' not in all_opts:
                if pre_clique == True:
                    if user_list_scores:
                        if save_graphs == "end":
                            m.MGA_bk(matching_dict, new_ao, save_graphs, all_clique, False, True, user_list_scores,
                                     scoring_for, forbidden=True)
                        else:
                            all_matchings = path + name + "_matching"
                            m.MGA_bk(matching_dict, new_ao, all_matchings, all_clique, False, True, user_list_scores,
                                     scoring_for, forbidden=True)

                    else:
                        if save_graphs == "end":
                            m.MGA_bk(matching_dict, new_ao, save_graphs, all_clique, False, False, None, scoring_for,
                                     forbidden=False)
                        else:
                            all_matchings = path + name + "_matching"
                            m.MGA_bk(matching_dict, new_ao, all_matchings, all_clique, False, False, None, scoring_for,
                                     forbidden=False)

                else:
                    if user_list_scores:
                        if save_graphs == "end":
                            # m.MGA_bk(matching_dict, new_ao, save_graphs, all_clique)
                            m.MGA_bk(matching_dict, new_ao, save_graphs, False, False, True, user_list_scores,
                                     scoring_for, forbidden=True)
                        else:
                            all_matchings = path + name + "_matching"
                            m.MGA_bk(matching_dict, new_ao, all_matchings, False, False, True, user_list_scores,
                                     scoring_for, forbidden=True)

                    else:
                        if save_graphs == "end":
                            m.MGA_bk(matching_dict, new_ao, save_graphs, False, False, False, None, scoring_for,
                                     forbidden=False)
                        else:
                            all_matchings = path + name + "_matching"
                            # m.MGA_bk(matching_dict, new_ao, all_matchings, all_clique)
                            m.MGA_bk(matching_dict, new_ao, all_matchings, False, False, False, None, scoring_for,
                                     forbidden=False)


            else:
                m.MGA_bk(matching_dict, new_ao, None, False, True, False, user_list_scores, scoring_for,
                         forbidden=False)


        except KeyError:
            if '-n' in all_opts:
                print("Error! Can't create Multiple Alignment. Check if tree corresponds to input graphs.")
                exit()

        time_bk_mult = time.time() - time_bk2

        # Save visualization of final graph alignment as PNG file.
        if '-d' in all_opts or '--drawm' in all_opts:
            png_matching = path + name + "_matching.png"
            if nodes_labelled:
                m.draw_matching(drawing_order, "label", png_matching)
            else:
                m.draw_matching(drawing_order, "no label", png_matching)
            print("Saved Matching as PNG")

        # Save final graph alignment as GRAPHML file.
        try:
            GraphIO.write_graphML_file(m, ml_path)
            print("Saved Final Matching as Graphml")
        except:
            print("Could not write to " + ml_path + " !")

        matches = path + name + "_matches.txt"
        GraphIO.write_matches(m.network, matches)


    # Multiple Graph Alignment via VF2.
    else:
        time_vf22 = time.time()

        try:
            matching_dict = {}
            for i in range(0, len(drawing_order)):
                matching_dict[drawing_order[i]] = graphen[drawing_order[i]]

            if user_list_scores and scoring_for == "nodes":
                multiple_aligner = vf.VF2_GraphAligner(new_ao, matching_dict, user_list_scores[0])
            elif user_list_scores and scoring_for == "edges":
                multiple_aligner = vf.VF2_GraphAligner(new_ao, matching_dict, None, user_list_scores[1])
            elif user_list_scores and scoring_for == "both":
                multiple_aligner = vf.VF2_GraphAligner(new_ao, matching_dict, user_list_scores[0], user_list_scores[1])
            else:
                multiple_aligner = vf.VF2_GraphAligner(new_ao, matching_dict)

            if save_graphs != "end":
                all_matchings = path + name + "_matching"
                mga = multiple_aligner.vf2mult(save_graphs, all_matchings,forbidden=forbidden)

            else:
                mga = multiple_aligner.vf2mult(save_graphs,forbidden=forbidden)

            # Save final graph alignment as GRAPHML file.
            GraphIO.write_graphML_file(mga, ml_path)
            print("Saved Final Matching as Graphml")
            matches = path + name + "_matches.txt"
            GraphIO.write_matches(mga, matches)

            # Save visualization of final graph alignment as PNG file.
            if '-d' in all_opts or '--drawm' in all_opts:
                png_matching = path + name + "_matching.png"
                if nodes_labelled:
                    multiple_aligner.output_generator(mga, drawing_order, "label", png_matching)
                else:
                    multiple_aligner.output_generator(mga, drawing_order, "nolabel", png_matching)
                print("Saved Matching as PNG")

        except KeyError:
            if '-n' in all_opts:
                print("Error! Can't create Multiple Alignment. Check if tree corresponds to input graphs.")
                exit()
                
        time_vf2_mult = time.time() - time_vf22

    # Benchmarking.
    if vf2:
        total = time_vf2_paar + time_vf2_mult
        print("Graph aligment used up a total time of: " + str(total) + " s.")
    else:
        total = time_bk_paar + time_bk_mult
        print("Graphs aligment used up a total time of: " +str(total) + " s.")

    benchm = path + name

    if vf2:
        benchm += "_vf2.txt"
        with open(benchm, 'a') as bencher:
            bencher.write(graph_loc + "\t" + str(time_vf2_paar) + "\t" + str(time_vf2_mult) + "\n")
    else:
        benchm += "_bk.txt"
        with open(benchm, 'a') as bencher:
            bencher.write(graph_loc + "\t" + str(time_bk_paar) + "\t" + str(time_bk_mult) + "\n")


"""
Authors: Anna Grunwald, Lena Suenke Mortensen, Christiane Gaertner
"""

if __name__ == '__main__':
    main(sys.argv[1:])
