import random
import os
import math
import string

def main():
    directory = os.path.dirname(os.path.abspath(__file__))
    if "\\" in directory:
        windows = True
    else:
        windows = False

    number = int(input("How many graphs do you want to create?"))
    dire = "."
    while dire.lower() != "y" and dire.lower() != "n":
        dire = input("directed graphs? y/n ")
    try:
        nodes = int(input("How many nodes shall the graph have: "))
    except:
        print("Please only use int for the number of nodes!")
        exit()
    nolab = "."
    while nolab.lower() != "y" and nolab.lower() != "n":
        nolab = input("Shall the graph nodes be labelled (labels randomly generated)? y/n ")

    try:
        prob = float(input("What probability shall an edge between two nodes have (between 0 and 1)?: "))
    except:
        print("Input error. Probability set to standard: 0.7")
        prob = 0.7
    edlab = "."
    while edlab.lower() != "y" and edlab.lower() != "n":
        edlab = input("Shall the graph edges be labelled (labels randomly generated)? y/n ")
    mesh = "."
    while mesh.lower() != "y" and mesh.lower() != "n":
        mesh = input("Do you want generate mesh graphs? y/n ")

    proj = input("Lastely name your graphs: ")
    if windows:
        path = directory + "\\" + proj
    else:
        path = directory + "/" + proj
    try:
        os.mkdir(path)
    except:
        path += "_1"
        os.mkdir(path)

    if dire == "y":
        directed = True
    else:
        directed = False

    if nolab == "y":
        nodes_labeled = True
    else:
        nodes_labeled = False

    if edlab == "y":
        edges_labeled = True
    else:
        edges_labeled = False

    if mesh == "y":
        mesh_graph = True
    else:
        mesh_graph = False

    for i in range(number):
        if windows:
            name = path + "\\" + proj + str(i+1) + ".graph"
        else:
            name = path + "/" + proj + str(i + 1) + ".graph"

        if mesh_graph:
            generate_mesh(name,nodes,directed,nodes_labeled,edges_labeled,prob)

        else:
            generate_graph(name,directed,nodes, nodes_labeled,edges_labeled, prob)



def generate_graph(name, dir, nodes, nodlabel, edgelabel, p_for_edge):
    nodelabels = []
    edges = []

    if nodlabel:
        while len(nodelabels) < nodes:
            label = random_labels(5)
            if label not in nodelabels:
                nodelabels.append(label)

    else:
        for i in range(nodes):
            nodelabels.append(str(i + 1))

    if dir:
        for i in range(len(nodelabels)):
            for j in range(len(nodelabels)):
                if i != j:
                    x = float(random.randint(1, 100)) / 100
                    if x <= p_for_edge:
                        edges.append((i + 1, j + 1))

        save_graph(name, dir, nodes, nodlabel, nodelabels, edges, edgelabel)

    else:
        for i in range(len(nodelabels)):
            for j in range(i + 1, len(nodelabels)):
                x = float(random.randint(1, 100) / 100)
                if x <= p_for_edge:
                    edges.append((i + 1, j + 1))

        save_graph(name, dir, nodes, nodlabel, nodelabels, edges, edgelabel)


# schreibt die generierten Daten in verzeichnis (siehe unten)
def save_graph(name, dir, nodes, nodlabel, nodelables, edges, edglabel):
    with open(name, 'a') as graph:
        graph.write("#author;random_generator\n")
        graph.write("#nodes;" + str(nodes) + "\n")
        graph.write("#edges;" + str(len(edges)) + "\n")
        if nodlabel:
            graph.write("Nodes labelled;True\n")
        else:
            graph.write("Nodes labelled;False\n")
        if edglabel:
            graph.write("Edges labelled;True\n")
        else:
            graph.write("Edges labelled;False\n")
        if dir:
            graph.write("Directed Graph;True\n")
        else:
            graph.write("Directed Graph;False\n")
        graph.write("\n")

        if nodlabel:
            for i in range(len(nodelables)):
                graph.write(str(i + 1) + ";" + nodelables[i] + "\n")
        else:
            for i in range(nodes):
                graph.write(str(i + 1) + "\n")

        if edges:
            if not edglabel:
                for i in edges:
                    graph.write("\n" + str(i[0]) + ";" + str(i[1]))
            else:
                for i in edges:
                    label = random_labels(8)
                    graph.write("\n" + str(i[0]) + ";" + str(i[1]) + ";" + label)

    print("saved " + name)


def random_labels(length):

    label = ""

    for i in range(length):
        label = label + random.choice(string.ascii_letters + string.digits)

    return label


def generate_mesh(name, nodes, directed, nodes_labeled, edges_labeled, probability):
    height = random.randint(1,nodes)

    width = math.ceil((nodes/height))

    mesh = []
    count = 1
    node_labels = []

    for j in range(height):
        mesh.append([])
        for n in range(width):
            if count <= nodes:
                mesh[j].append(count)
                node_labels.append(random_labels(5))
                count += 1

            else:
                mesh[j].append("NA")
        if count > nodes:
            break

    edges = []

    for i in range(len(mesh) - 1):
        for j in range(len(mesh[0]) - 1):
            x = float(random.randint(1, 100)) / 100
            if x <= probability:
                edges.append((mesh[i][j], mesh[i + 1][j]))
            x = float(random.randint(1, 100)) / 100
            if x <= probability:
                edges.append((mesh[i][j], mesh[i][j + 1]))

    for i in range(len(mesh) - 1):
        x = float(random.randint(1, 100)) / 100
        if x <= probability:
            edges.append((mesh[i][len(mesh[0]) - 1], mesh[i + 1][len(mesh[0]) - 1]))

    for i in range(len(mesh[0]) - 1):
        x = float(random.randint(1, 100)) / 100
        if x <= probability:
            edges.append((mesh[len(mesh) - 1][i], mesh[len(mesh) - 1][i + 1]))

    if directed:
        for i in range(1, len(mesh)):
            for j in range(1, len(mesh[0])):
                x = float(random.randint(1, 100)) / 100
                if x <= probability:
                    edges.append((mesh[i][j], mesh[i - 1][j]))
                x = float(random.randint(1, 100)) / 100
                if x <= probability:
                    edges.append((mesh[i][j], mesh[i][j - 1]))

        for i in range(1, len(mesh)):
            x = float(random.randint(1, 100)) / 100
            if x <= probability:
                edges.append((mesh[i][0], mesh[i - 1][0]))

        for i in range(1, len(mesh[0])):
            x = float(random.randint(1, 100)) / 100
            if x <= probability:
                edges.append((mesh[0][i], mesh[0][i - 1]))

    for i in edges[:]:
        if "NA" in i:
            edges.remove(i)


    save_graph(name, directed, nodes, nodes_labeled, node_labels, edges, edges_labeled)
