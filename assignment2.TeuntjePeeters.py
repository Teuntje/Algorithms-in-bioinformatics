#!/usr/bin/env python
"""
Author: Teuntje Peeters
Student number: 920301645120
Implementation of the SSAHA algorithm

Hints:
- write a function to generate a hash table of k-mers from a set of sequences
- write a function to return a sorted list of hits in the hash table for a query sequence
- write a function to find the best hit
- write a fasta parser to read in the Arabidopsis data

"""
# import statements
from sys import argv
import itertools
import pprint
from collections import Counter

def hashTable(seqs, k):
	"""Return a hashtable with non-overlapping kmers

	seqs: list with sequences
	k: integer of length of kmer
	"""

	#Creation of an empty dictionary
	tupHashTable = {}

	nrSeq = 0
	for seq in seqs:
		for n in xrange(0,len(seq)-1, k):
			#Create kmers
			kmer = "".join([seq[n:n+k]])

			#Skip if length is not the same as the k
			if not len(kmer) == k: 
				continue
			#Add kmer to tupHashtable if it is already there
			#Otherwise add this kmer
			if kmer in tupHashTable:
				tupHashTable[kmer].append((nrSeq, n))
			else:
				tupHashTable[kmer] = [(nrSeq, n)]
			
		nrSeq += 1
	return tupHashTable

def printHash(hashTable, indent):
	"""Nice printing of the Hash table, like table 1 from the article

	hashTable: dictionary with keys (kmers) and values (positions)
	indent: integer for formatting the distance
	"""
	for key, value in hashTable.iteritems():
		print '\t' * indent + str(key)
		if isinstance(value, dict):
			pretty(value, indent+1)
		else:
			print '\t' * (indent+1) + str(value)

def dbSearch(seq, db, k):
	"""Return list of hits sorted by index and shift
	
	seq: input sequence
	db: dictionary with kmers as keys and locations as values
	k: integer of length of the kmer
	"""
	match_list = []

	#Creation of kmers of the query sequence
	for q_idx in xrange(0,len(seq), 1):
		kmer = seq[q_idx: q_idx + k]
		#Creation of the list of hits
		if kmer in db:
			for value in db[kmer]:
				match_list.append((value[0], value[1] - q_idx, value[1]))

	#Return a sorted list based on index 0 and index 1
	return sorted(sorted(match_list, key = lambda x: x[1]), key = lambda x: x[0])


def matching(M):
	"""Return the longest match

	M: List with matches of the query sequence (list with tuples inside)
	"""
	allMatches = []
	longestMatch = []

	#Check if list is longer than 1
	if len(M) > 1: 
		for match in range(len(M)-1):
			#If these values are the same, add them to allMatches list
			if (M[match][1] == M[match+1][1]) and (M[match][0] ==  M[match+1][0]):
				allMatches.append(M[match][0:2])
				allMatches.append(M[match+1][0:2])

	#If not longer than one, return this match
	else:
		allMatches = M[0:2]
		return allMatches

	#If allMatches is not empty, calculate the longest match
	#and append to LongestMatch
	if allMatches:
		countMatches = Counter(allMatches)
		longestAl =  max(countMatches.iterkeys(), key=lambda k: countMatches[k])
		for i in M:
			if i[0] == longestAl[0] and i[1] == longestAl[1]:
				longestMatch.append(i)

	return longestMatch

def alignment(seqs, Qs, longestMatch, k):
	"""Return alignment with longest match
	
	seqs: list with sequences in the database
	Qs: query sequence, string
	longestMatch: list with the locations of the longest match
	k: integer of kmer length
	"""
	#Get the sequence and positions in sequence
	seq = seqs[longestMatch[0][0]]
	posinseq = [item[2] for item in longestMatch]

	aligned_seq = seq[posinseq[0]:posinseq[-1]+k]

	#Check of substring
	front =  Qs.index(aligned_seq)
	back = Qs[::-1].index(aligned_seq[::-1])

	aligneed_seq_start = posinseq[0]

	#Print the alignment nicely
	print(seq)
	if len(aligned_seq) != len(Qs):
		spaces = " " * (aligneed_seq_start-front)
		print(spaces + "-" * front + "|" * len(aligned_seq) + "-" * back)
		print(spaces + Qs)

	else: 
		spaces = " " * aligneed_seq_start
		print(spaces + "|" * len(aligned_seq))
		print(spaces + aligned_seq)
	# posinseq =  range(posinseq[0], posinseq[len(posinseq)-1]+k)

def parseFasta(fileName):
    """Return fasta file as dictionary, key is name of the sequence and value is sequence

	fileName: string with filename
    """
    fastaDic = {}
    #Open and read the file
    with open(fileName) as inFile:
    	for line in inFile:
    		line = line.strip()
    		#Creation of key and values
    		if line.startswith(">"):
    			key = line
    			fastaDic[key] = []
    		else:
    			fastaDic[key].append(line)
    #Values are all strings in lists, so they should be joined
    fastaDic = {key: "".join(value) for key, value in fastaDic.items()}

    return fastaDic


if __name__ == "__main__":

    # the code below should produce the results necessary to answer the questions
    # in other words, if we run your code, we should see the data that you used
    # to answer the questions
    
    query = 'TGCAACAT'
    
    s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
    seqs = [s1,s2,s3]

    #print "Seq 1: ", s1
    #print "Seq 2: ", s2
    #print "Seq 3: ", s3

    #Question 1
    print("Question 1: print the Hash table")
    HashTable = hashTable(seqs, k=2)
    printHash(HashTable, indent = 0)

    #Question 2
    print("\nQuestion 2: printing a list of hits sorted by index and shift")
    seq = "TGCAACAT"
    M = dbSearch(seq, HashTable, k=2)
    print "Number of hits: ", len(M)
    print "First hit: ", M[0]
    print "Last hit: ", M[len(M)-1], "\n"

    #Question 3
    print "Question 3: return the longest match between the query and the database"
    longestM = matching(M)
    print "The longest match: ", longestM, "\n"
    print "The alignment without any gaps: "
    alignment(seqs, query, longestM, k=2)
    #Question 4
    print "\nQuestion 4: How many sequences does the genome file contain and what is the total sequence length?"
    genome = parseFasta("TAIR10.fasta")
    print "The genome file contains:", len(genome.keys()), "sequences"
    print "The total sequence length: ", sum(len(v) for v in genome.values())
    print "\nDoes this match your expectation on the Arabidopis genome?"
    print "Not exactly, I'd have expected that it would contain more chromosomes"

    #Question 5
    print "\nQuestion 5: Building a hash table for the aribidopsis chromosomes"
    HashTable = hashTable(genome.values(), k=15)
    #printHash(HashTable, indent = 0)
    print "\nAmount of keys in this hash table: ", len(HashTable.keys())

    #Question 6
    print "\nQuestion 6: find the maximum hit for athal_query.fasta"
    athal = parseFasta("athal_query.fasta")
    #HashTable = hashTable(athal.values(), k=5)
    for q in athal.values():
        M = dbSearch(q, HashTable, k=15)
        print "Number of hits: ", len(M)
        print "First hit: ", M[0]
        print "Last hit: ", M[len(M)-1], "\n"
        longestM = matching(M)
        print "The longest match: ", longestM, "\n"





