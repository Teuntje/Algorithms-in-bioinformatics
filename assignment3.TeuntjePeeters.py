#!/usr/bin/env python

"""
Author: Teuntje Peeters

Description: this is a script to calculate HMM
"""

from sys import argv
from random import random
import collections
from collections import Counter
import numpy as np

# Background amino acid probabilities
pa = { 'A':0.074, 'C':0.025, 'D':0.054, 'E':0.054, 'F':0.047, 'G':0.074, 'H':0.026, 'I':0.068, 'L':0.099, 'K':0.058, 'M':0.025, 'N':0.045, 'P':0.039, 'Q':0.034, 'R':0.052, 'S':0.057, 'T':0.051, 'V':0.073, 'W':0.013, 'Y':0.034 } 


class HMM():
    """HMM object to store an HMM model

    This object is designed to keep track of all HMM states, emissions, and 
    transitions. It may be used in your implementation, but may also be ignored,
    and replaced by a data structure of choice
    """
    # Emission probabilities for the match and insert states
    e_m   = []; e_i   = pa; 
    
    # Transition probabilities from/to matches, inserts and deletions
    t_mm  = []; t_mi  = []; t_md = [];
    t_im  = []; t_ii  = []
    t_dm  = []; t_dd  = []; 
    
    def __init__(self,nmatches):
        """Initialize HMM object with number of match states
        
        nmatches: int, number of match states
        """
    
        self.nmatches = nmatches
        
        self.e_m   = [dict(pa) for i in range(0,nmatches)]
        for i in range(0,nmatches):
            for j in pa.keys():
                self.e_m[i][j] = 0.0
        self.e_i   = pa;

        self.t_mm  = [0.0 for i in range(0,nmatches+1)]
        self.t_mi  = [0.0 for i in range(0,nmatches+1)]
        self.t_im  = [0.0 for i in range(0,nmatches+1)]
        self.t_ii  = [0.0 for i in range(0,nmatches+1)]
        self.t_md  = [0.0 for i in range(0,nmatches+1)]
        self.t_dm  = [0.0 for i in range(0,nmatches+1)]
        self.t_dd  = [0.0 for i in range(0,nmatches+1)]

def sample(events):
    """Return a key from dict based on the probabilities 

    events: dict of {key: probability}, sum of probabilities
    should be 1.0. 
    """
    k = events.keys()
    cum = [0 for i in k]

    cum[0] = events[k[0]]
    for i in range(1,len(events)):
        cum[i] = cum[i-1] + events[k[i]]
    # Should not be necessary, but for safety
    cum[len(cum)-1] = 1.0

    r = random()
    pick = ''
    i = 0
    while (pick == '') and (i < len(cum)):
        if r < cum[i]:
            pick = k[i]
        i = i + 1
    return pick

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
                key = line[1:]
                fastaDic[key] = []
            else:
                fastaDic[key].append(line)
    #Values are all strings in lists, so they should be joined
    fastaDic = {key: "".join(value) for key, value in fastaDic.items()}

    return fastaDic

def matchStates(file, k):
    """Return the number of matchStates, 
    the positions where a match is found and
    the list of indexes with inserts

    file: dictionary with fasta file
    k: integer, percentage for determining if it is a match state or not
    """
    #Split the sequence into amino acids
    l =  [list(i) for i in file.values()]
    #Zip this list
    pos = zip(*l)
    #Count the amout ot "-"
    counts = [i.count("-") for i in pos]
    num = len(l)*k
    #Get list of indexes where it is an insert
    index = [j for j, i in enumerate(counts) if i>num]
    #Get the match states
    matchState = sum(i < num for i in counts)
    pos = [i for j, i in enumerate(pos) if j not in index]

    return pos, index, matchState


def emissions(positions):
    """Return dictionary with emissions probabilities

    positions: tuple in list where index is the position in alignment
    """

    #Create an empty list
    emission = []
    #For each position, count the occurences
    for position in positions:
        c = Counter(position)
        #Delete the gap
        c.pop('-', None)
        #Calculate the probabilities
        total = sum([float(i) for i in c.values()])
        percent = {key:value/total for key, value in c.items()}
        #Add into the list
        emission.append([percent])

    return emission

#These were my attempts of calculation of the transition
# def transitions(file, index):
#     """Return number of transitions

#     """
#     positions = zip(*file.values())
    
#     matches = []
#     insertions = []
#     deletions = []

#     for idx in range(0,len(positions)):
#         print positions[idx]
#         print idx
#         if idx not in index:
#             print "MATCH"
#             c = positions[idx].count("-")
#             deletion = float(c)/float(len(positions[idx]))
#             deletions.append(deletion)
#             matches.append(1.00-deletion)
#             insertions.append(0.0)
#             print "match: ", 1.00-deletion
#             print "deletion: ", deletion
#             print "insertion: ", 0.0
#         else:
#             print "INSERTION"
#             c = positions[idx].count("-")
#             insert = float(c)/float(len(positions[idx]))
#             insertions.append(1.00-insert)
#             insertions.append(1.00)
#             matches.append(insert)
#             deletions.append(0.0)
#             print "match: ", insert
#             print "insertion: ", 1-insert
#             print "deletion: ", 0.0
            
#     print "Matches: ", matches
#     print "Insertions: ", insertions
#     print "Deletions: ", deletions

# def transitions(file, index):
#     """Return the transitions

#     file: parsed fasta file, dictionary
#     index: the positions of the insertions in the sequence
#     """
#     #zip the positions
#     positions = zip(*file.values())

#     #Creation of three dictionaries for matches, deletions and insertions
#     M = {}
#     D = {}
#     I = {}

#     #For each index in the positions
#     for idx in range(len(positions)):
#         #If index is not in 'index', it is not an insertion
#         if idx not in index:
#             #Count the amount of "-"
#             c = positions[idx].count("-")
#             #if the amount is bigger than 0
#             if c > 0:
#                 deletion = float(c)/float(len(positions[idx]))
#                 D["".join(["MD", str(idx)])] = 1-deletion
#                 if idx+1 not in index:
#                     M["".join(["M", str(idx)])] = deletion
#             #If the amount is smaller than 0
#             else:
#                 print positions[idx]
#         #It is an insertion
#         else:
#             c = positions[idx].count("-")
#             insert = float(c)/float(len(positions[idx]))
#             I["".join(["MI", str(idx)])] = 1 - insert
#             M["".join(["M", str(idx)])] = insert
#     print D
#     print M
#     print I

# def transitions(f, index):
#     """

#     """

    
#     print len(f.values()[0])-1

#     HMM = {}

#     state = ""

#     for sequence in f.values():
#         print sequence
#         for pos in range(len(sequence)):
#             if pos not in index:
#                 if sequence[pos] != "-":
#                     match = "".join(["M", str(pos+1)])
#                     state = "match"
#                     if match in HMM:
#                         HMM[match] += 1
#                     else: 
#                         HMM[match] = 1
#                 else:
#                     deletion = "".join(["D", str(pos+1)])
#                     if deletion in HMM:
#                         HMM[deletion] += 1
#                     else:
#                         HMM[deletion] = 1
#                     currentstate = "D"
#             else:
#                 if sequence[pos] != "-":
#                     insertion = "".join(["I", str(pos)])
#                     if insertion in HMM:
#                         HMM[insertion] += 1
#                     else:
#                         HMM[insertion] = 1
#                     currentstate = "I"
#                 else:
#                     print "Go to the next match"
#     print HMM

def determine_state(s_id, i, l, insertstates):
    """Return the state of position

    s_id: sequence identity, which sequence in the list, integer
    i: current position
    l: list with sequences
    insertstates: list with the insert positions
    """

    #Get the amino acid at this position
    amino = l[s_id][i]
    #Is it an insert state? 
    if i in insertstates:
        if amino != "-":
            return("INSERT")
    else:
        #Is it a match or deletion? 
        if amino == "-":
            return("DELETE")
        else:
            return("MATCH")

# def transitions(f, insertstates):
#     """
#     """

#     l = [i for i in f.values()]

#     n_states = len(l[0]) - len(insertstates)
#     matches = [[0,0,0] for i in range(n_states)]
#     insertions = [[0,0,0] for i in range(n_states+1)]
#     deletions = [[0,0,0] for i in range(n_states)]

#     for s_id, s in enumerate(l):
#         match_pos = -1
#         for i, a in enumerate(s[:-1]):
#             if i in insertstates:
#                 pass
#                 #from ins -> x
#             else:
#                 match_pos+=1
#                 for t_i, target in enumerate(s[i+1:]):
#                     if t_i+i+1 not in insertstates:
#                         #Target is match

#                         #Target is delete

#                     else: #insert state
#                         if target == '-':
#                             continue
#                         else: #target is insertion
#                             if determine_state(s_id, i, l, insertstates) == 'MATCH':
#                                 matches[match_pos][1] = matches[match_pos][1] + 1
#                                 break
#                             else:


#     print matches

def transitions(f, insertstates):
    """Return three lists for matches, inserts and deletions
    with the counts.

    In these lists, each position is the index the position of the state
    So M1, M2, M3 etc. Within each position there is another list with 
    the counts of the neighbors. So where does it go? [MATCH, INSERT, DEL]

    f: dictionary of fasta file
    insertstates: list with indexes with insert states
    """

    #initialise all the lists
    l = [i for i in f.values()]

    #Number of states and create the correct lists with values of 0
    n_states = len(l[0]) - len(insertstates)
    matches = [[0,0,0] for i in range(n_states+1)]
    insertions = [[0,0,0] for i in range(n_states+1)]
    deletions = [[0,0,0] for i in range(n_states+1)]

    #Fill the first position in each list
    for idx, seq in enumerate(l):
        #Check the state
        state = determine_state(idx, 0, l, insertstates)
        #Add up all the counts for MATCH, INSERT and DELETION
        if state == "MATCH":
            matches[0][0] += 1
        elif state == "INSERT":
            matches[0][1] += 1
        elif state == "DELETE":
            matches[0][2] += 1

    #Loop over all the sequences and keep track of sequence ID and seq
    for idx_seq, seq in enumerate(l):
        match_pos = 1
        #Go over all the amino acids in all the sequences
        for i in range(len(seq)-1):
            #Check the current state (MATCH, INSERTION or DELETION)
            current_state = determine_state(idx_seq, i, l, insertstates)
            for target in range(i+1,len(seq)):
                #Check the target state, neighbor
                target_state = determine_state(idx_seq, target, l, insertstates)
                #Check for match, insertion and deletions in both the current and
                #target state. Count the number of occurences
                if current_state == "MATCH":
                    if target_state == "MATCH":
                        matches[match_pos][0] += 1
                        match_pos += 1
                        break
                    if target_state == "INSERT":
                        matches[match_pos][1] += 1
                        break
                    if target_state == "DELETE":
                        matches[match_pos][2] += 1
                        match_pos += 1
                        break
                if current_state == "INSERT":
                    if target_state == "MATCH":
                        insertions[match_pos][0] += 1
                        match_pos += 1
                        break
                    if target_state == "INSERT":
                        insertions[match_pos][1] += 1
                        break
                    if target_state == "DELETE":
                        insertions[match_pos][2] += 1
                        match_pos += 1
                        break
                if current_state == "DELETE":
                    if target_state == "MATCH":
                        deletions[match_pos][0] += 1
                        match_pos += 1
                        break
                    if target_state == "INSERT":
                        deletions[match_pos][1] += 1
                        break
                    if target_state == "DELETE":
                        deletions[match_pos][2] += 1
                        match_pos += 1
                        break

    return matches, insertions, deletions

def normalisationTransitions(input_lists):
    """Return a list with normalised counts in the matches, insertions
    and deletions

    input_lists: one list with three lists containing match counts, insertion
    counts and deletion counts
    """
    output_list = []
    for l in input_lists:
        output_list.append([list(np.array(x, dtype=np.float64)/sum(x)) if sum(x) > 0 else [0,0,0] for x in l])

    return output_list

def sampleSequence(emissions, m, ins, d):
    """Return the sample sequence as a string

    emissions: list with at each position a dictionary of the 
    emissions of the amino acids
    m: normalised matches, list
    ins: normalised inserts, list
    d: normalised deletions, list
    """

    #Define the states and beginning of the sequence
    states = ["MATCH", "INSERTION", "DELETE"]
    begin = np.random.choice(states, p = m[0])

    seq = []
    i = 0
    
    #Define the previous list based on one of three conditions
    if begin == "MATCH":
        prev_list = m
    elif begin == "INSERTION":
        prev_list = ins
    elif begin == "DELETE":
        prev_list = d

    #while i is smaller than length of the match list (randomly chosen, equal length)
    while i < len(m)-1:
        #Get a random choice at the ith position in states
        begin = np.random.choice(states, p = prev_list[i])
        #Choose from the emissions an amino acid based on the probabilities
        if begin == 'MATCH':
            #seq.append(sample(pa))
            seq.append(sample(emissions[i][0]))
            prev_list = m
        elif begin == 'INSERTION':
            seq.append(sample(pa))
            prev_list = ins
            continue
        elif begin == "DELETE":
            prev_list = d
        i += 1

    return "".join(seq)

def outputPSSM(emission):
    """Return output of PSSM for the logo

    emission: list with dictionary for each position with emission prob
    """
    for i in emission:
        for j in i:
            for a, v in j.items():
                print a, v

           
if __name__ == "__main__":
    infile = "test_large.fasta"
    f = parseFasta(infile)
    
    #Question 1, 
    #print "Question 1"
    positions, index, matchS = matchStates(f, k = 0.5)
    #print "The number of match states is: ", matchS

    #Question 2
    emission = emissions(positions)
    #outputPSSM(emission)
    print "Emissions: "
    for i, pos in enumerate(emission):
        print "position: ", i+1
        print pos

    matches, insertions, deletions = transitions(f, index)
    print "transitions:"
    nmatches, ninsertions, ndeletions = normalisationTransitions([matches, insertions, deletions])
    print "Normalised matches: ", nmatches
    print "Normalised insertions: ", ninsertions
    print "Normalised deletions: ", ndeletions
    #Normalise the matches, insertions and deletions

    seq = sampleSequence(emission, nmatches, ninsertions, ndeletions)
    print seq



