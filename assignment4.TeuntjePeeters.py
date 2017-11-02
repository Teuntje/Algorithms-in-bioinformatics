#!/usr/bin/env python

"""
Author: Teuntje Peeters
Student nr: 920301645120
Script to: make graphs and genome assembly
"""
import copy


def is_eulerian(graph):
	""" Return if graph is balanced, boolean
	
	Function to determine whether a graph is Eulerian

	graph: dictionary with graph
	"""
	#for each node in keys and values
	for node in graph.items():
		#initialise edges as empty list
		edges = []
		#For each edge in graph.items()
		for edge in graph.items():
			#Check positions in each edge 
			#(because it can have a longer length)
			for pos in edge[1]:
				#add counts to list edges
				edges.append(pos.count(node[0]))
		#If length of the edges corresponding to the node is equal to counts
		#of the occurences in the edges
		if (len(node[1]) == sum(edges)) == False:
			return False
		else:
			return True
			

def has_eulerian_path(graph):
	""" Return if a graph is an Eulerian path, boolean

	Function to determine whether a graph has an Eulerian path
	A connected graph has an Eulerian path if and only if it contains 
	at most two semibalanced vertices and all other vertices are balanced.
	A vertex v in a graph is called semibalanced if 
	indegree(v) - outdegree(v) = 1

	graph: dictionary with graph
	"""
	#initiate variables
	count = 0
	balanced = []
	#semi_balanced = []
	#For each position in the graph
	for position in graph.items():
		vertex = position[0]
		#Outdegree is list going 'out' of the vertex
		outdegree = len(position[1])
		indegree = []
		#For each position, get the edge
		for pos in graph.items():
			edge = pos[1]
			#For each position in the edge
			for e in edge: 
				#Count the number of occurences of vertex in edges
				indegree.append(e.count(vertex))
		#print sum(indegree) == outdegree
		#If the number of coming in edges minus the out-going edges
		#is equal to 1, count till the maximum of 2
		if (sum(indegree) - outdegree) == 1:
		 	count += 1
		 	#semi_balanced.append(position)
		 	#If bigger than 2, return false
		 	if count > 2:
		 		return False #, semi_balanced
		#Extra check if the rest of the vertices is balanced
		elif (sum(indegree) == outdegree) == True:
			balanced.append(1)
		elif (sum(indegree) == outdegree) == False:
			balanced.append(0)
	#If the rest is balanced, return true
	if sum(balanced) >= (len(graph)-2):
		return True #, semi_balanced



def find_eulerian_cycle(graph):
	"""Return Eulerian cycle in a graph as list

	Function to find an Eulerian cycle in a graph

	graph: dictionary with graph
	"""

	#Initiate variables
	#Creation of a copy of the original grap
	copyGraph = copy.deepcopy(graph)
	#Start node
	current = list(copyGraph)[0]
	subcycle = [current]

	correctCycle = []
	#While length of subcycle is bigger than 0
	while len(subcycle) > 0:
		#If the vertex still has edges
		if copyGraph[current] != []:
			#Append the current to the subcycle list
			subcycle.append(copyGraph[current][0])
			#Get the next vertex
			nextV = copyGraph[current].pop(0)
			#Change current vertex to 'next' for the next 'round'
			current = nextV
		else:
			#Add the last item of subcycle to the correct Cycle
			correctCycle.append(subcycle[-1])
			#Delete it from the subcycle
			del subcycle[-1]
			#If subcycle is empty, break
			if subcycle == []:
				break
			else: 
				#Replace current with the last one from subcycle
				current = subcycle[-1]
	#Return the correct Cycle reverse
	return correctCycle[::-1]


def find_eulerian_path(graph, start):
	"""Return an Eulerian path in a graph

	Function to find an Eulerian path in the graf

	graph, dictionary with graph
	"""

	#Initiate variables
	copyGraph = copy.deepcopy(graph)
	if start == "":
		current = list(copyGraph)[0]
	else:
		current = list(copyGraph)[start]
	subpath = [current]

	correctPath = []
	#Check is length subpath is not equal to 0
	while len(subpath) > 0:
		#If current not in graph
		if current not in copyGraph:
			#Append to the correctCycle and delete from subpath
			correctPath.append(subpath[-1])
			del subpath[-1]
			#If subpath is empty, break
			if subpath == []:
				break
			else: 
				#Change current to last element of subpath
				current = subpath[-1]
		#If current in copyGraph is not empty
		if copyGraph[current] != []:
			#Append to subpath
			subpath.append(copyGraph[current][0])
			#Get the next variable
			nextV = copyGraph[current].pop(0)
			current = nextV
		else:
			#Append to correctPath
			correctPath.append(subpath[-1])
			del subpath[-1]
			if subpath == []:
				break
			else: 
				current = subpath[-1]
	
	return correctPath[::-1]	

def create_graph(s):
	"""Return graph as dictionary

	Function to create a graph with a string
	The keys are the first two letters of the kmer
	The value are the last two values of the kmer
	So the middle value (kmer = 3) is in both key and value

	s list with kmers
	"""

	#Initiate variables
	graph = {}
	#For kmer in s
	for kmer in s:
		#If first two letters are not in the graph
		if kmer[0:2] not in graph:
			#Add key to the graph with empty list
			graph[kmer[0:2]] = []
		#Append the rest behind the correct key
		graph[kmer[0:2]].append(kmer[1:])

	return graph

def prettyPrint(s):
	"""Nice printing of the graph

	s, list
	"""
	#First position, create list with strings
	string = [str(s[0])]
	#For element in string
	for el in s[1:]:
		#Append to the list with strings
		string.append(' -> {}'.format(str(el)))

	#Join the whole string and print this
	print "".join(string)

def printSequence(s):
	"""Nice printing of the sequence

	graph, graph of the sequence
	"""
	#Print the sequence
	string = [str(s[0])]
	for el in s[1:]:
		string.append('{}'.format(str(el[1])))
	print "".join(string)


if __name__ == "__main__":

    # GRAPH FROM FIG 8.22
    graph_822 = {'A':['B'],'B':['C'],'I':['H'],'H':['F'],'F':['G','E'],\
        'C':['I','J'],'G':['A'],'E':['J'],'J':['F','D'],'D':['C']}

    # A SLIGHTLY BIGGER GRAPH, NEEDED FOR Q8
    bigger_graph = {5:[6],6:[7],10:[11],11:[4],4:[5,3],\
        7:[10,9],3:[9,1],9:[4,8],8:[7],1:[2], 2:[3]}

    #example_class = {5:[8, 3], 8:[9], 9: [5], 3: [6, 1, 4], 6: [7], 7: [3], 1: [2], 2: [3], 4: [5]}

    #print find_eulerian_path(example_class, "")
    # SPECTRUM FROM FIG 8.20
    s = ['ATG','TGG','TGC','GTG','GGC','GCA','GCG','CGT']

    print "Question 1"
    print "This graph is balanced:", is_eulerian(graph_822)

    print "\nQuestion 2"
    b = has_eulerian_path(graph_822)
    print "Has eulerian path:", b

    cycle = find_eulerian_cycle(graph_822)

    print "\nQuestion 3"
    path = find_eulerian_path(graph_822, "")
    prettyPrint(path)

    print "\nQuestion 4"
    prettyPrint(path)
    prettyPrint(path)
    prettyPrint(path)

    print "\nQuestion 5"
    graph = create_graph(s)

    for k, v in graph.items():
    	print k, v

    print "\nQuestion 6"
    print "This graph is balanced:", is_eulerian(graph)
    b = has_eulerian_path(graph)
    print "Has eulerian path:", b

    print "\nQuestion 7"
    start = list(graph).index('AT')
    path = find_eulerian_path(graph, start)
    prettyPrint(path)
    printSequence(path)

    print "\nQuestion 8"
    path = find_eulerian_path(bigger_graph, "")
    print "path: "
    prettyPrint(path)
    cycle = find_eulerian_cycle(bigger_graph)
    print "cycle: "
    prettyPrint(cycle)

    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file) 
