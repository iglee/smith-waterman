# imports
import pandas as pd
import numpy as np
import argparse
from itertools import permutations
from random import sample, randint
import sys
import time

startTime = time.time()
# arguments for the script
parser = argparse.ArgumentParser(description = "Implementation of Smith-Waterman Local Alignment.")
parser.add_argument("--str-input", action="store_true", help="indicates string inputs would be given. i.e. AKA")
parser.add_argument("--file-input", action="store_true", help="indicates fasta file inputs would be given.")
parser.add_argument("-A", metavar="A", action='store', type=str, help="one of the input sequences")
parser.add_argument("-B", metavar="B", action='store', type=str, help="the other one of the input sequences. this is the input that would be used for permutations.")
parser.add_argument("-af", metavar="af", action='store', type=str, help="one of the input sequences *in a fasta file*")
parser.add_argument("-bf", metavar="bf", action='store', type=str, help="one of the other input sequences *in a fasta file*")
parser.add_argument("-p", action="store_true", help = "calculate the p-value. with B string permutations")
parser.add_argument("-o", action="store", type=str, help="output file directory")
args = parser.parse_args()

f_out = open(args.o, 'w')

if args.str_input:
    A = list(args.A.upper())
    B = list(args.B.upper())

# read fasta files
def read_data(filename):
    f = open(filename, "r")
    seq = f.readlines()
    seq_out = []
    f.close()
    for l in seq[1:]:
        seq_out += list(l.strip())
    return seq_out

if args.file_input:
    A = read_data(args.af)
    B = read_data(args.bf)

# BLOSUM scoring matrix as Global variable- read in raw file and convert to dataframe
f = open("BLOSUM62.txt", "r")
blosum_raw = f.readlines()[6:]

temp = []
for x in blosum_raw:
    l = x.strip().split()
    temp.append(l)

blosum_processed = []
for x in temp[1:]:
    blosum_processed.append(x[1:])

BLOSUM = pd.DataFrame(blosum_processed, columns = temp[0], index = temp[0]).apply(pd.to_numeric)


# amino acid codes as global variable
amino_acid_code = temp[0][:20]

# initialize the score matrix
V = np.zeros((len(A)+1,len(B)+1))

# score function with BLOSUM matrix
def score(V, A, B, i, j):
    opt = [0, V[i-1][j]-4, V[i][j-1]-4, V[i-1][j-1]+BLOSUM[A[i-1]][B[j-1]]]
    #print("blosum score:",BLOSUM[A[i-1]][B[j-1]])
    #print(opt)
    return max(opt)

# calculate the score matrix
for i in range(1,len(A)+1):
    for j in range(1,len(B)+1):
        #print(i,j)
        V[i][j] = score(V, A, B, i, j)


if len(A) <= 15 and len(B) <= 15:
    print("printing the score matrix: ")
    print(V)
    print("Score Matrix:\n", file=f_out)
    print(V, file=f_out)

# location of the max score
(max_i,max_j) = np.unravel_index(np.argmax(V), V.shape)
print("location of the maximum score: ")
print(max_i,", " , max_j)

alignment_score = np.max(V)
print("maximum score: ", alignment_score)
print("\nmaximum score: ", alignment_score, file=f_out)

# backtrace and decode sequence
def backtrack(V, i, j):
    subseq = []
    substringA = []
    substringB = []
    
    while V[i][j] != 0:
        subseq.append((i-1,j-1,V[i][j]))
        if V[i-1][j-1] >= V[i-1][j] and V[i-1][j-1] >= V[i][j-1]:
            substringA.append(A[i-1])
            substringB.append(B[j-1])
            i -= 1
            j -=1
        elif V[i-1][j] > V[i-1][j-1] and V[i-1][j] > V[i][j-1]:
            substringA.append(A[i-1])
            substringB.append("-")
            i -= 1
        elif V[i][j-1] > V[i-1][j-1] and V[i][j-1] >= V[i-1][j]:
            substringA.append("-")
            substringB.append(B[j-1])
            j -= 1
    return subseq, ''.join(substringA[::-1]), ''.join(substringB[::-1])

subseq, substringA, substringB = backtrack(V,max_i,max_j)

if args.p:
    # generate permutations of input string B
    N = 999
    # permBs = sample(list(permutations(B)), N)  # use this as a check for fisher yates shuffle -- ONLY FOR SMALL INPUTS < 15 tokens


    # use Fisher Yates shuffle to generate random permutations
    def fs_shuffle(B):
        seq = B.copy()
        n = len(seq)
        for i in range(n-1, 0, -1):
            j = randint(0,i)
            seq[i], seq[j] = seq[j], seq[i]
        return seq

    permBs = []

    for x in range(N):
        permBs.append(fs_shuffle(B))


    def score_perms(A, B_perm):
        B_perm = list(B_perm)
        
        V = np.zeros((len(A)+1,len(B_perm)+1))
        
        for i in range(1,len(A)+1):
            for j in range(1,len(B_perm)+1):
                V[i][j] = score(V, A, B_perm, i, j)
        return np.max(V)

    all_scores = []

    print("Calculating scores for all permutations...")
    m = 0
    lst = time.time()
    for x in permBs:
        all_scores.append(score_perms(A,x))
        m += 1
        if m % 20 == 0:
            #sys.stdout.write(" . " + str(m))

            loopTime = (time.time() - lst)
            print('Calculation time for',m,' score matrices seconds: ' + str(loopTime))
    print("Done!")

    score_opts = set(all_scores)

    k = 0
    for x in score_opts:
        if x >= alignment_score:
            k+=all_scores.count(x)

    p = (k+1)/(N+1)
    print("\nempirical p-value: ", p, file=f_out)
    print("empirical p-value: ", p)

executionTime = (time.time() - startTime)
print('Total execution seconds: ' + str(executionTime))

f_out.close()