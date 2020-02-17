# MIGRAINE
MIGRAINE = MultIple GRaph AlIgNmEnt tool via Bron Kerbosch and VF2 algorithms

![title](https://github.com/Tschichen/MIGRAINE/blob/master/pictures/CF2708A0-5EB2-48E4-B90E-1A8B449FC132.jpeg "title picture")

## Remarks for practical course
* There is now an option to remove all hydrogen from the chemical graphs. Type -H to use it. Speeds up the algorithms.
* New option --opt will give you an additional output, only showing the perfectly matched subgraph without suboptimals.

## General description
A python3-based tool that takes multiple graph objects and aligns them in a progressive way using either Bron-Kerbosch or VF2 algorithm.
  * Input:
    * Two or more graphs in either JSON, GRAPHML or the tool-specific GRAPH format.
  * Output:
    * A multiple graph matching in GRAPHML format.
    
## Install

```git clone https://github.com/Tschichen/MIGRAINE.git```

Dependencies can be installed using [pipenv](https://github.com/pypa/pipenv).
Simply run:
```
pipenv install
```

and precede all commands with `pipenv run` to execute them in a virtualized environment.

## Input examples
* minimum input  
    ```python3 migraine.py-i users/path/dir_with_graphs/ -o my_name```
	
* reads newick tree from file  
 ```python3 migraine.py -i users/path/dir_with_graphs -o my_name -n users/path/some_dir/tree.nwtree ```
	
* uses BK and predefines a clique, saves all alignment steps  
	```python3 migraine.py -i users/path/dir_with_graphs/ -o my_name -b -c users/path/dir/clique.json -m all```
	
* reads a list with scores and scores labels according to list, saves guide tree as newick  
	```python3 migraine.py -i users/path/dir_with_graphs/ -o my_name -s users/path/some_dir/score_list.txt -q```
	

## Options

Note: Clique option only available for BK. Disconnected directed graphs can only be assessed via BK algorithm. Input graphs must be all directed or undirected.

option| explanation
------------ | -------------
--graphgen | starts built in tool, that builds (multiple) randomly created graphs. Graph properties like node count, edge probability or random labels need to be customized.  
-h or --help |  provides a list of options that can used with this tool and additional information
-i <inputpath> | Use -i and specify a path where all the graph files, you want to be aligned, are stored. All files in JSON, GRAPHML or GRAPH format in this location will be used. **Mandatory Option!**
-o <name> | Name your current run. A directory with that name will be generated in your current location and all outputfiles will be saved there. **Mandatory Option!**
-b or --bronkerbosch | As the tool uses VF2 algorithm by default, use this parameter to switch to using Bron Kerbosch algorithm.
-v | random pivoting strategy for bk algorithm. Combination with pre-clique option is excluded. _Default: costum pivoting._
-c <graph-file> | Clique-option: Uses the specified graph-file as anchor for Bron Kerbosch. Checks if clique is in all graphs that are to be aligned and starts multiple alignment from pre-defined clique. Input is possible in all three formats. _Only works with Bron Kerbosch._
-H | Option to filter hydrogen from chemical graphs
-s or --score <file.txt> | Allows to input a matrix of scores for node and edge labels. Check the bottom of this file for format example. _Default: If no scoring list is given, the score of an alignment is generated via the size of the largest mapping._
-f or --forbidden |	while inputting a scoring list, score via size of largest common subgraph after excluding forbidden label matches.
-g or --labelling <string> | String can be "nodes", "edges" or "both". Defines what labels should be used for scoring.
-n or --newick <newick-file> | Input a guide tree for multiple alignment in newick format. If this option is not enabled the guide tree is built from scores computed via pairwise alignment.
-r or --randomt | saves a randomly generated guide tree for your current run
-t or --showt | saves the generated guide tree as .png
-q or --savetn | saves the generated guide tree as newick file
-m <string> | saves the multiple alignment as GRAPHML file. Options for string are: "end" only end result will be saved, "all" saves every matching after one graph is added to the alignment till the end. _Default: end_
-d or --drawm | saves a visualization of the multiple graph alignment as .png
--opt | will provide an additional output, only showing the perfectly matched subgraph without suboptimals.


	
## Required Python libraries
- networkx
- sys
- os
- getopt
- glob
- math
- json
- matplotlib
- collections
- operator
- random
- pprint
- string
- Bio
- cStringIO
- time

## Format examples

### Graph format
#### in general:

	#author; <string>
	#nodes; <int>
	#edges; <int>
	Nodes labelled; <True/False>
	Edges labelled; <True/False>
	Directed Graph; <True/False>
	SPACE
	node count;Label(optional)
	...
	SPACE
	node1;node2;Label(optional)
	...
	
#### example:		
				    
	#author;someone					
	#nodes;2					      
	#edges;1					      
	Nodes labelled;True			
	Edges labelled;True			
	Directed Graph;False		
							          
	1;Label						      
	2;Label	
	
	1;2;Label					     

### Score list format
***Please type N, if match between two labels is forbidden.***

#### in general:
    
    #nodelabel1 \t nodelabel2 \t nodelabel3...
    score13 \t score23 \t score33
    score12	\t score22	
    score11
    #edgelabel1 \t edgelabel2 \t edgelabel3...
    score13	\t score23 \t score33
    score12	\t score22	
    score11

#### example:

    #C	N	H	B   
    5	3	2	1
    4	N	3
    2	1
    4
    #s	d	i
    3	3	2
    7	N
    2

in the example, a match between a node labeled with "N" and a node labeled with  "H" would be forbidden, score between edge labels "s" and "i" would be 3.
