#!/usr/bin/env python

"""
Author: Teuntje Peeters
Student number: 920301645120

Description: this is a script to make a global alignment of two sequences
"""


# functions between here and __main__
blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""

def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int,parts[1:])))
    return order, blosum_matrix

BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()

def score(res1, res2):
    """Return similarity score from BLOSUM62 matrix for two residues
    
    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]

# write your own functions below here
def gapPentalty(ngap, E):
    """Return total gap penalty
    
    ngap: integer
    E: integer, penalty for gaps
    """
    return (-ngap) * E

def templateMatrix(seq1, seq2, penalty = 4):
    """Return the template matrix 
    
    seq1: sequence 1, string
    seq2: sequence 2, string
    penalty: gap penalty, integer
    """
    matrix = [[None for i in range(len(seq1)+1)] for j in range(len(seq2)+1)]
    matrix[0][0] = 0

    for i in range(len(matrix[0])):
      matrix[0][i] =  gapPentalty(i, penalty)

    for i in range(len(matrix)):
      matrix[i][0] =  gapPentalty(i, penalty)

    #for line in matrix:
    #  print "\t".join(map(str, line))

    return matrix

def NW(seq1, seq2, penalty):
    """Needle Wunsch algorithm
    
    seq1: sequence 1, string
    seq2: sequence 2, string
    penalty: integer for the gap penalty
    """

    #Creation of the alignment scoring matrix and the parent matrix
    #The penalty is overwritten to make it negative
    scoringMatrix = templateMatrix(seq1, seq2, penalty)
    pMatrix = templateMatrix(seq1, seq2, penalty)
    penalty = -(penalty)

    #Looping over matrix anc calculation of the scorings
    try:
        for i in range(1, len(seq2)+1):
            for j in range(1, len(seq1)+1):
                current = scoringMatrix[i-1][j-1]
                right = scoringMatrix[i][j-1] + penalty
                diagonal = current + score(seq2[i-1], seq1[j-1])
                down = scoringMatrix[i-1][j] + penalty
                val = {"right" : right, "diagonal" : diagonal, "down" : down}

                #Get the maximum scoring to put on the right position
                #Fill the parent matrix with the directions for the traceback
                scoringMatrix[i][j] = max(right, diagonal, down)
                pMatrix[i][j] = max(val, key = val.get)
    
                
    except IndexError:
        print "There is an indexerror"

    return pMatrix, scoringMatrix


def alignSeq(parentMatrix, scoringMatrix, seq1, seq2, endGap, penalty):
    """Create the final alignment based on the scoring matrix 
    
    parentMatrix, list in list with the instructions for the right path
    scoringMatrix, list in list with the scorings for each position
    seq1: sequence 1, string
    seq2: sequence 2, string    
    endGap: integer for gap penalty at the beginning and end
    penalty: integer for gap penalty within the sequences
    """

    #Adjust matrix, delete the numbers and replace with directions
    matrix = adjustMatrix(parentMatrix)

    #Creation of variables
    s1, s2 = len(seq1), len(seq2)
    a1, a2 = "", ""
    x, y = 0, 0
    count= 0

    #Loop over the matrix for the traceback
    #The if statements are telling where to go in the matrix
    while matrix[s2-x][s1-y] != matrix[0][0]:
        if matrix[s2-x][s1-y] == "right":
            a1 += seq1[len(seq1)-1-y]
            a2 += "-"
            y = y + 1
            count += 1
        if matrix[s2-x][s1-y] == "diagonal":
            #print scoringMatrix[s2-x][s1-y]
            a1 += seq1[len(seq1)-1-y]
            a2 += seq2[len(seq2)-1-x]
            y = y + 1
            x = x + 1
            count += 1
        if matrix[s2-x][s1-y] == "down":
            #print scoringMatrix[s2-x][s1-y]
            a2 += seq2[len(seq2)-1-x]
            a1 += "-"
            x = x + 1
            count += 1
    scoring = scoringMatrix[s2][s1]

    #Get the alignment
    alignment = [a1[::-1], a2[::-1]]

    #Calculation on the end gaps
    #The end gaps are calculated in this function (endGapCalc)
    #The penalty is added to the alignment again for correcting for the end gaps
    #The penalty of the endgaps is added to the total scoring
    ends = endGapCalc(alignment)
    scoring += ends * penalty
    scoring -= ends*endGap

    #Printing of the output
    print "Alignment score: ", scoring
    return alignment

def endGapCalc(alignment):
    """Calculate the amount of the endgaps
    
    alignment: list in list with the alignments
    """
    counter = 0
    for i in alignment:
        if i.startswith("-"):
            for j in i:
                if j == "-":
                    counter += 1
                else:
                    break
        if i.endswith("-"):
            for j in i[::-1]:
                if j == "-":
                    counter += 1
                else:
                    break

    return counter


def adjustMatrix(matrix):
    """Adjust the align matrix

    Removal of the scoring matrix, replace it with 'instructions', as down and right

    matrix this is a list in list with the instructions for the alignment
    """

    for i in range(len(matrix[0])):
      matrix[0][i] =  "right"

    for i in range(len(matrix)):
      matrix[i][0] =  "down"

    matrix[0][0] = 0

    return matrix

def calcIdentity(alignment):
    """Calculation of the identity between the two alignments
    
    alignment, string in list with the alignment
    """
    listA, listB = list(alignment[0]), list(alignment[1])
    zipLists = zip(listA, listB)

    ident = 0
    lines = ""
    for i in zipLists:
        if i[0] == i[1]:
            ident += 1
            lines += "|"
        else:
            lines += " "

    identity = (float(ident) / float(len(listA))) * 100
    return identity, lines

def overView(seq1, seq2, penalty, endGap):
    """This function creates an overview of the functions that are about to run

    seq1, string with first sequence
    seq2, string with second sequence
    penalty, integer for the gap penalty
    endGap, integer of gap penalty at the end and beginnen of the two sequences
    """

    #Creation of parent matrix (for traceback) and scoring matrix (with alignment scorings)
    parentMatrix, scoringMatrix = NW(seq1, seq2, penalty)

    #Creation of a list with the two alignments and incorporated gaps
    alignment= alignSeq(parentMatrix, scoringMatrix, seq1, seq2, endGap, penalty)

    #Calculation of identity scoring and the lines are for a nice output of the two alignments
    identity, lines = calcIdentity(alignment)

    #Printing the alignment nicely
    print alignment[0]
    print lines
    print alignment[1]
    print "Identity: ", round(identity, 2), "%"

if __name__ == "__main__":
    seq1 = "THISLINE"
    seq2 = "ISALIGNED"
    # seq3: GPA1_ARATH
    seq3 = "MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYMLSSESIAIGEKLSEIGGRLDYPRLTKDIAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIPTKEDVLYARVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHLFEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSFMLFLNKFDIFEKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYRTTALDQKLVKKTFKLVDETLRRRNLLEA"
    # seq4: GPA1_ORYSI
    seq4 = "MGSSCSRSHSLSEAETTKNAKSADIDRRILQETKAEQHIHKLLLLGAGESGKSTIFKQIKLLFQTGFDEAELRSYTSVIHANVYQTIKILYEGAKELSQVESDSSKYVISPDNQEIGEKLSDIDGRLDYPLLNKELVLDVKRLWQDPAIQETYLRGSILQLPDCAQYFMENLDRLAEAGYVPTKEDVLYARVRTNGVVQIQFSPVGENKRGGEVYRLYDVGGQRNERRKWIHLFEGVNAVIFCAAISEYDQMLFEDETKNRMMETKELFDWVLKQRCFEKTSFILFLNKFDIFEKKIQKVPLSVCEWFKDYQPIAPGKQEVEHAYEFVKKKFEELYFQSSKPDRVDRVFKIYRTTALDQKLVKKTFKLIDESMRRSREGT"

    #Autimatically run the functions

    #Question 1
    print "Question 1"
    overView(seq1, seq2, penalty=4, endGap = 1)
    
    #Question 2
    #See report

    #Question 3
    print "\nQuestion 3"
    print "\nPenalty = 1: "
    overView(seq1, seq2, penalty=1, endGap = 0)
    print "\nPenalty = 5: "
    overView(seq1, seq2, penalty=5, endGap = 0)
    print "\nPenalty = 8: "
    overView(seq1, seq2, penalty=8, endGap = 0)


    #Question 4
    print "\nQuestion 4\nSequence 3 and 4"
    print"\nPenalty = 5 and endGap = 1: "
    overView(seq3, seq4, penalty = 5, endGap = 1)
    print"\nPenalty = 5 and endGap = 10: "
    overView(seq3, seq4, penalty = 5, endGap = 10)

