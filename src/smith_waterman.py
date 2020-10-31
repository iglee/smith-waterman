# imports
import pandas as pd
import numpy as np
import argparse

# arguments for the script
parser = argparse.ArgumentParser(description = "Implementation of Smith-Waterman Local Alignment.")
parser.add_argument("-A", metavar="A", action='store', type=str, help="one of the input sequences")
parser.add_argument("-B", metavar="B", action='store', type=str, help="the other one of the input sequences")
args = parser.parse_args()

A = args.A
B = args.B

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

# read fasta files
def read_data(filename):
    f = open(filename, "r")
    seq = f.readlines()
    seq_out = []

    for l in seq[1:]:
        seq_out += list(l.strip())
    return seq_out

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

print("printing the score matrix: ")
print(V)
