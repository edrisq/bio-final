#!/usr/bin/env python

from fasta.fasta import parse_fasta, ORFs_in_fasta
from os.path import exists
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--seq1', nargs=1, dest='seq1', type=str, help='Path of input fasta file', required=True)
parser.add_argument('--seq2', nargs=1, dest='seq2', type=str, help='Path of input fasta file', required=True)
parser.add_argument('--type', nargs=1, dest='type', type=str, help='Options: local, global, both', required=True)
parser.add_argument('--matrix', nargs=1, dest='matrix', type=str, help='Options: pam250, blosum62, simple', required=True)
parser.add_argument('--gap', nargs=1, dest='gap', type=int, help='Neg ints are gap penalties', required=True)

# Copied from https://github.com/biopython/biopython/blob/master/Bio/SubsMat/MatrixInfo.py
# http://www.embl-heidelberg.de/~vogt/matrices/blosum62.cmp
blosum62 = {
    ("W", "F"): 1, ("L", "R"): -2, ("S", "P"): -1, ("V", "T"): 0,
    ("Q", "Q"): 5, ("N", "A"): -2, ("Z", "Y"): -2, ("W", "R"): -3,
    ("Q", "A"): -1, ("S", "D"): 0, ("H", "H"): 8, ("S", "H"): -1,
    ("H", "D"): -1, ("L", "N"): -3, ("W", "A"): -3, ("Y", "M"): -1,
    ("G", "R"): -2, ("Y", "I"): -1, ("Y", "E"): -2, ("B", "Y"): -3,
    ("Y", "A"): -2, ("V", "D"): -3, ("B", "S"): 0, ("Y", "Y"): 7,
    ("G", "N"): 0, ("E", "C"): -4, ("Y", "Q"): -1, ("Z", "Z"): 4,
    ("V", "A"): 0, ("C", "C"): 9, ("M", "R"): -1, ("V", "E"): -2,
    ("T", "N"): 0, ("P", "P"): 7, ("V", "I"): 3, ("V", "S"): -2,
    ("Z", "P"): -1, ("V", "M"): 1, ("T", "F"): -2, ("V", "Q"): -2,
    ("K", "K"): 5, ("P", "D"): -1, ("I", "H"): -3, ("I", "D"): -3,
    ("T", "R"): -1, ("P", "L"): -3, ("K", "G"): -2, ("M", "N"): -2,
    ("P", "H"): -2, ("F", "Q"): -3, ("Z", "G"): -2, ("X", "L"): -1,
    ("T", "M"): -1, ("Z", "C"): -3, ("X", "H"): -1, ("D", "R"): -2,
    ("B", "W"): -4, ("X", "D"): -1, ("Z", "K"): 1, ("F", "A"): -2,
    ("Z", "W"): -3, ("F", "E"): -3, ("D", "N"): 1, ("B", "K"): 0,
    ("X", "X"): -1, ("F", "I"): 0, ("B", "G"): -1, ("X", "T"): 0,
    ("F", "M"): 0, ("B", "C"): -3, ("Z", "I"): -3, ("Z", "V"): -2,
    ("S", "S"): 4, ("L", "Q"): -2, ("W", "E"): -3, ("Q", "R"): 1,
    ("N", "N"): 6, ("W", "M"): -1, ("Q", "C"): -3, ("W", "I"): -3,
    ("S", "C"): -1, ("L", "A"): -1, ("S", "G"): 0, ("L", "E"): -3,
    ("W", "Q"): -2, ("H", "G"): -2, ("S", "K"): 0, ("Q", "N"): 0,
    ("N", "R"): 0, ("H", "C"): -3, ("Y", "N"): -2, ("G", "Q"): -2,
    ("Y", "F"): 3, ("C", "A"): 0, ("V", "L"): 1, ("G", "E"): -2,
    ("G", "A"): 0, ("K", "R"): 2, ("E", "D"): 2, ("Y", "R"): -2,
    ("M", "Q"): 0, ("T", "I"): -1, ("C", "D"): -3, ("V", "F"): -1,
    ("T", "A"): 0, ("T", "P"): -1, ("B", "P"): -2, ("T", "E"): -1,
    ("V", "N"): -3, ("P", "G"): -2, ("M", "A"): -1, ("K", "H"): -1,
    ("V", "R"): -3, ("P", "C"): -3, ("M", "E"): -2, ("K", "L"): -2,
    ("V", "V"): 4, ("M", "I"): 1, ("T", "Q"): -1, ("I", "G"): -4,
    ("P", "K"): -1, ("M", "M"): 5, ("K", "D"): -1, ("I", "C"): -1,
    ("Z", "D"): 1, ("F", "R"): -3, ("X", "K"): -1, ("Q", "D"): 0,
    ("X", "G"): -1, ("Z", "L"): -3, ("X", "C"): -2, ("Z", "H"): 0,
    ("B", "L"): -4, ("B", "H"): 0, ("F", "F"): 6, ("X", "W"): -2,
    ("B", "D"): 4, ("D", "A"): -2, ("S", "L"): -2, ("X", "S"): 0,
    ("F", "N"): -3, ("S", "R"): -1, ("W", "D"): -4, ("V", "Y"): -1,
    ("W", "L"): -2, ("H", "R"): 0, ("W", "H"): -2, ("H", "N"): 1,
    ("W", "T"): -2, ("T", "T"): 5, ("S", "F"): -2, ("W", "P"): -4,
    ("L", "D"): -4, ("B", "I"): -3, ("L", "H"): -3, ("S", "N"): 1,
    ("B", "T"): -1, ("L", "L"): 4, ("Y", "K"): -2, ("E", "Q"): 2,
    ("Y", "G"): -3, ("Z", "S"): 0, ("Y", "C"): -2, ("G", "D"): -1,
    ("B", "V"): -3, ("E", "A"): -1, ("Y", "W"): 2, ("E", "E"): 5,
    ("Y", "S"): -2, ("C", "N"): -3, ("V", "C"): -1, ("T", "H"): -2,
    ("P", "R"): -2, ("V", "G"): -3, ("T", "L"): -1, ("V", "K"): -2,
    ("K", "Q"): 1, ("R", "A"): -1, ("I", "R"): -3, ("T", "D"): -1,
    ("P", "F"): -4, ("I", "N"): -3, ("K", "I"): -3, ("M", "D"): -3,
    ("V", "W"): -3, ("W", "W"): 11, ("M", "H"): -2, ("P", "N"): -2,
    ("K", "A"): -1, ("M", "L"): 2, ("K", "E"): 1, ("Z", "E"): 4,
    ("X", "N"): -1, ("Z", "A"): -1, ("Z", "M"): -1, ("X", "F"): -1,
    ("K", "C"): -3, ("B", "Q"): 0, ("X", "B"): -1, ("B", "M"): -3,
    ("F", "C"): -2, ("Z", "Q"): 3, ("X", "Z"): -1, ("F", "G"): -3,
    ("B", "E"): 1, ("X", "V"): -1, ("F", "K"): -3, ("B", "A"): -2,
    ("X", "R"): -1, ("D", "D"): 6, ("W", "G"): -2, ("Z", "F"): -3,
    ("S", "Q"): 0, ("W", "C"): -2, ("W", "K"): -3, ("H", "Q"): 0,
    ("L", "C"): -1, ("W", "N"): -4, ("S", "A"): 1, ("L", "G"): -4,
    ("W", "S"): -3, ("S", "E"): 0, ("H", "E"): 0, ("S", "I"): -2,
    ("H", "A"): -2, ("S", "M"): -1, ("Y", "L"): -1, ("Y", "H"): 2,
    ("Y", "D"): -3, ("E", "R"): 0, ("X", "P"): -2, ("G", "G"): 6,
    ("G", "C"): -3, ("E", "N"): 0, ("Y", "T"): -2, ("Y", "P"): -3,
    ("T", "K"): -1, ("A", "A"): 4, ("P", "Q"): -1, ("T", "C"): -1,
    ("V", "H"): -3, ("T", "G"): -2, ("I", "Q"): -3, ("Z", "T"): -1,
    ("C", "R"): -3, ("V", "P"): -2, ("P", "E"): -1, ("M", "C"): -1,
    ("K", "N"): 0, ("I", "I"): 4, ("P", "A"): -1, ("M", "G"): -3,
    ("T", "S"): 1, ("I", "E"): -3, ("P", "M"): -2, ("M", "K"): -1,
    ("I", "A"): -1, ("P", "I"): -3, ("R", "R"): 5, ("X", "M"): -1,
    ("L", "I"): 2, ("X", "I"): -1, ("Z", "B"): 1, ("X", "E"): -1,
    ("Z", "N"): 0, ("X", "A"): 0, ("B", "R"): -1, ("B", "N"): 3,
    ("F", "D"): -3, ("X", "Y"): -1, ("Z", "R"): 0, ("F", "H"): -1,
    ("B", "F"): -3, ("F", "L"): 0, ("X", "Q"): -1, ("B", "B"): 4
}

# Copied from https://github.com/biopython/biopython/blob/master/Bio/SubsMat/MatrixInfo.py
# http://www.embl-heidelberg.de/~vogt/matrices/pam250.cmp
pam250 = {
    ("W", "F"): 0, ("L", "R"): -3, ("S", "P"): 1, ("V", "T"): 0,
    ("Q", "Q"): 4, ("N", "A"): 0, ("Z", "Y"): -4, ("W", "R"): 2,
    ("Q", "A"): 0, ("S", "D"): 0, ("H", "H"): 6, ("S", "H"): -1,
    ("H", "D"): 1, ("L", "N"): -3, ("W", "A"): -6, ("Y", "M"): -2,
    ("G", "R"): -3, ("Y", "I"): -1, ("Y", "E"): -4, ("B", "Y"): -3,
    ("Y", "A"): -3, ("V", "D"): -2, ("B", "S"): 0, ("Y", "Y"): 10,
    ("G", "N"): 0, ("E", "C"): -5, ("Y", "Q"): -4, ("Z", "Z"): 3,
    ("V", "A"): 0, ("C", "C"): 12, ("M", "R"): 0, ("V", "E"): -2,
    ("T", "N"): 0, ("P", "P"): 6, ("V", "I"): 4, ("V", "S"): -1,
    ("Z", "P"): 0, ("V", "M"): 2, ("T", "F"): -3, ("V", "Q"): -2,
    ("K", "K"): 5, ("P", "D"): -1, ("I", "H"): -2, ("I", "D"): -2,
    ("T", "R"): -1, ("P", "L"): -3, ("K", "G"): -2, ("M", "N"): -2,
    ("P", "H"): 0, ("F", "Q"): -5, ("Z", "G"): 0, ("X", "L"): -1,
    ("T", "M"): -1, ("Z", "C"): -5, ("X", "H"): -1, ("D", "R"): -1,
    ("B", "W"): -5, ("X", "D"): -1, ("Z", "K"): 0, ("F", "A"): -3,
    ("Z", "W"): -6, ("F", "E"): -5, ("D", "N"): 2, ("B", "K"): 1,
    ("X", "X"): -1, ("F", "I"): 1, ("B", "G"): 0, ("X", "T"): 0,
    ("F", "M"): 0, ("B", "C"): -4, ("Z", "I"): -2, ("Z", "V"): -2,
    ("S", "S"): 2, ("L", "Q"): -2, ("W", "E"): -7, ("Q", "R"): 1,
    ("N", "N"): 2, ("W", "M"): -4, ("Q", "C"): -5, ("W", "I"): -5,
    ("S", "C"): 0, ("L", "A"): -2, ("S", "G"): 1, ("L", "E"): -3,
    ("W", "Q"): -5, ("H", "G"): -2, ("S", "K"): 0, ("Q", "N"): 1,
    ("N", "R"): 0, ("H", "C"): -3, ("Y", "N"): -2, ("G", "Q"): -1,
    ("Y", "F"): 7, ("C", "A"): -2, ("V", "L"): 2, ("G", "E"): 0,
    ("G", "A"): 1, ("K", "R"): 3, ("E", "D"): 3, ("Y", "R"): -4,
    ("M", "Q"): -1, ("T", "I"): 0, ("C", "D"): -5, ("V", "F"): -1,
    ("T", "A"): 1, ("T", "P"): 0, ("B", "P"): -1, ("T", "E"): 0,
    ("V", "N"): -2, ("P", "G"): 0, ("M", "A"): -1, ("K", "H"): 0,
    ("V", "R"): -2, ("P", "C"): -3, ("M", "E"): -2, ("K", "L"): -3,
    ("V", "V"): 4, ("M", "I"): 2, ("T", "Q"): -1, ("I", "G"): -3,
    ("P", "K"): -1, ("M", "M"): 6, ("K", "D"): 0, ("I", "C"): -2,
    ("Z", "D"): 3, ("F", "R"): -4, ("X", "K"): -1, ("Q", "D"): 2,
    ("X", "G"): -1, ("Z", "L"): -3, ("X", "C"): -3, ("Z", "H"): 2,
    ("B", "L"): -3, ("B", "H"): 1, ("F", "F"): 9, ("X", "W"): -4,
    ("B", "D"): 3, ("D", "A"): 0, ("S", "L"): -3, ("X", "S"): 0,
    ("F", "N"): -3, ("S", "R"): 0, ("W", "D"): -7, ("V", "Y"): -2,
    ("W", "L"): -2, ("H", "R"): 2, ("W", "H"): -3, ("H", "N"): 2,
    ("W", "T"): -5, ("T", "T"): 3, ("S", "F"): -3, ("W", "P"): -6,
    ("L", "D"): -4, ("B", "I"): -2, ("L", "H"): -2, ("S", "N"): 1,
    ("B", "T"): 0, ("L", "L"): 6, ("Y", "K"): -4, ("E", "Q"): 2,
    ("Y", "G"): -5, ("Z", "S"): 0, ("Y", "C"): 0, ("G", "D"): 1,
    ("B", "V"): -2, ("E", "A"): 0, ("Y", "W"): 0, ("E", "E"): 4,
    ("Y", "S"): -3, ("C", "N"): -4, ("V", "C"): -2, ("T", "H"): -1,
    ("P", "R"): 0, ("V", "G"): -1, ("T", "L"): -2, ("V", "K"): -2,
    ("K", "Q"): 1, ("R", "A"): -2, ("I", "R"): -2, ("T", "D"): 0,
    ("P", "F"): -5, ("I", "N"): -2, ("K", "I"): -2, ("M", "D"): -3,
    ("V", "W"): -6, ("W", "W"): 17, ("M", "H"): -2, ("P", "N"): 0,
    ("K", "A"): -1, ("M", "L"): 4, ("K", "E"): 0, ("Z", "E"): 3,
    ("X", "N"): 0, ("Z", "A"): 0, ("Z", "M"): -2, ("X", "F"): -2,
    ("K", "C"): -5, ("B", "Q"): 1, ("X", "B"): -1, ("B", "M"): -2,
    ("F", "C"): -4, ("Z", "Q"): 3, ("X", "Z"): -1, ("F", "G"): -5,
    ("B", "E"): 3, ("X", "V"): -1, ("F", "K"): -5, ("B", "A"): 0,
    ("X", "R"): -1, ("D", "D"): 4, ("W", "G"): -7, ("Z", "F"): -5,
    ("S", "Q"): -1, ("W", "C"): -8, ("W", "K"): -3, ("H", "Q"): 3,
    ("L", "C"): -6, ("W", "N"): -4, ("S", "A"): 1, ("L", "G"): -4,
    ("W", "S"): -2, ("S", "E"): 0, ("H", "E"): 1, ("S", "I"): -1,
    ("H", "A"): -1, ("S", "M"): -2, ("Y", "L"): -1, ("Y", "H"): 0,
    ("Y", "D"): -4, ("E", "R"): -1, ("X", "P"): -1, ("G", "G"): 5,
    ("G", "C"): -3, ("E", "N"): 1, ("Y", "T"): -3, ("Y", "P"): -5,
    ("T", "K"): 0, ("A", "A"): 2, ("P", "Q"): 0, ("T", "C"): -2,
    ("V", "H"): -2, ("T", "G"): 0, ("I", "Q"): -2, ("Z", "T"): -1,
    ("C", "R"): -4, ("V", "P"): -1, ("P", "E"): -1, ("M", "C"): -5,
    ("K", "N"): 1, ("I", "I"): 5, ("P", "A"): 1, ("M", "G"): -3,
    ("T", "S"): 1, ("I", "E"): -2, ("P", "M"): -2, ("M", "K"): 0,
    ("I", "A"): -1, ("P", "I"): -2, ("R", "R"): 6, ("X", "M"): -1,
    ("L", "I"): 2, ("X", "I"): -1, ("Z", "B"): 2, ("X", "E"): -1,
    ("Z", "N"): 1, ("X", "A"): 0, ("B", "R"): -1, ("B", "N"): 2,
    ("F", "D"): -6, ("X", "Y"): -2, ("Z", "R"): 0, ("F", "H"): -2,
    ("B", "F"): -4, ("F", "L"): 2, ("X", "Q"): -1, ("B", "B"): 3
}

# Add reverse orientation of keys to dictionary
for key in [k for k in blosum62]:
    blosum62[(key[1], key[0])] = blosum62[key]
for key in [k for k in pam250]:
    pam250[(key[1], key[0])] = pam250[key]

alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
            'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
            'W', 'X', 'Y', 'Z']
letter_tuples = [pair for sublist in [[(A, B) for A in alphabet] for B in alphabet] for pair in sublist]
simple_matrix = {lt: (1 if lt[0] == lt[1] else -1) for lt in letter_tuples}

def global_max(f, sources, seq1, seq2):
    # Start at the last character in both sequences
    i = len(seq2)
    j = len(seq1)
    seq1_str = ''
    seq2_str = ''
#     sol_stack = [(i, j)]
    
#     while sol_stack:
#         i, j = sol_stack.pop()
#         current = sources[i][j].pop()
#         seq1_str = seq1_str[j:]
#         seq2_str = seq2_str[i:]

    # Decrement until i=0 or j= 0 (in which case, there are no further sources to trace)
    while (i and j):
#             if (len(sources[i][j][0]) > 1):
#                 sol_stack.append((i, j))

        # Current implementation only looks at sources[i][j][0], which is left-first
        if (sources[i][j][0][0] and sources[i][j][0][1]):
            seq1_str = seq1[j - 1] + seq1_str
            seq2_str = seq2[i - 1] + seq2_str
        elif sources[i][j][0][0]:
            seq1_str = '_' + seq1_str
            seq2_str = seq2[i - 1] + seq2_str
        elif sources[i][j][0][1]:
            seq1_str = seq1[j - 1] + seq1_str
            seq2_str = '_' + seq2_str
        
        incj = sources[i][j][0][1] if i else -1
        i += sources[i][j][0][0] if j else -1
        j += incj
    
    # If seq1 is exhausted, we pair each remaining character in seq2 with a space
    while i:
        seq1_str = '_' + seq1_str
        seq2_str = seq2[i - 1] + seq2_str
        i -= 1
    
    # If seq2 is exhausted, we pair each remaining character in seq1 with a space
    while j:
        seq1_str = seq1[j - 1] + seq1_str
        seq2_str = '_' + seq2_str
        j -= 1
    
    identity = 0
    for i in range(len(seq1_str)):
        if seq1_str[i] == seq2_str[i]:
            identity += 1

#     print(f'Alignment score = {f[len(seq2)][len(seq1)]}')
#     print(f'Percent identity = {identity / len(seq1_str) * 100:.2f}%')
#     print(seq1_str)
#     print(seq2_str)
    
    return (f[len(seq2)][len(seq1)], identity, seq1_str, seq2_str)

def local_max(f, sources, seq1, seq2):
    seq1_str = ''
    seq2_str = ''
    local_max = 0
    
    # Find the largest value across the entire matrix
    for m in range(len(seq2) + 1):
        for n in range(len(seq1) + 1):
            if f[m][n] > local_max:
                i = m
                j = n
                local_max = f[i][j]

    # Trace back until a 0 is encountered
    while (f[i][j]):
        # Current implementation only looks at sources[i][j][0], which is left-first
        if (sources[i][j][0][0] and sources[i][j][0][1]):
            seq1_str = seq1[j - 1] + seq1_str
            seq2_str = seq2[i - 1] + seq2_str
        elif sources[i][j][0][0]:
            seq1_str = '_' + seq1_str
            seq2_str = seq2[i - 1] + seq2_str
        elif sources[i][j][0][1]:
            seq1_str = seq1[j - 1] + seq1_str
            seq2_str = '_' + seq2_str
        
        incj = sources[i][j][0][1] if i else -1
        i += sources[i][j][0][0] if j else -1
        j += incj

    identity = 0
    for i in range(len(seq1_str)):
        if seq1_str[i] == seq2_str[i]:
            identity += 1
  
    return (local_max, identity, seq1_str, seq2_str)
#     print(f'Alignment score = {local_max}')
#     print(f'Percent identity = {identity / len(seq1_str) * 100:.2f}%')
#     print(seq1_str)
#     print(seq2_str)

def mat_print(matrix, seq1, seq2):
    width = 0
    
    # Set width based on widest entry
    for row in matrix:
        for col in row:
            width = max(width, len(str(col)))

    # Adjust width to provide some padding
    width = 5 if width < 5 else width
    width += 2 if width < 20 else 0
    
    # Print the first row (chars in seq1)
    print(" ", end="")
    for letter in seq1:
        print(f"{letter:^{width}}", end="")
    print()
    
    # Print remaining rows (each prepended by char in seq2)
    for i in range(len(seq2)):
        print(f"{seq2[i]}", end="")
        for j in range(len(seq1)):
            print(f"{str(matrix[i+1][j+1]):^{width}}", end="")
        print()

def needleman_wunsch(seq1, seq2, matrix, gap):
    # f matrix stores scores
    f = [[0 for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    
    # Source matrix stores a list of origins of each cell in f matrix 
    sources = [[[] for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    
    # First row and column initialized
    for j in range(len(seq1) + 1):
        f[0][j] = j * gap
    for i in range(len(seq2) + 1):
        f[i][0] = i * gap
    
    # Populate the f matrix, tracking sources
    for j in range(len(seq1)):
        for i in range(len(seq2)):
            up = f[i][j + 1] + gap
            left = f[i + 1][j] + gap
            diag = f[i][j] + matrix[(seq2[i],seq1[j])]
            best = max(up, left, diag)

            # As only the first is read, this is left-first
            if best == left:
                sources[i + 1][j + 1].append((0, -1))
            if best == diag:
                sources[i + 1][j + 1].append((-1, -1))
            if best == up:
                sources[i + 1][j + 1].append((-1, 0))
                
            f[i+1][j+1] = best

    return (f, sources)
#     global_max(f, sources, seq1, seq2)
#     mat_print(f, seq1, seq2)
#     mat_print(sources, seq1, seq2)
    

def smith_waterman(seq1, seq2, matrix, gap):
    f = [[0 for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    sources = [[[] for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    
    for j in range(len(seq1)):
        for i in range(len(seq2)):
            up = f[i][j + 1] + gap
            left = f[i + 1][j] + gap
            diag = f[i][j] + matrix[(seq2[i],seq1[j])]
            best = max(up, left, diag, 0)
            
            if best == left:
                sources[i + 1][j + 1].append((0, -1))
            if best == diag:
                sources[i + 1][j + 1].append((-1, -1))
            if best == up:
                sources[i + 1][j + 1].append((-1, 0))
                
            f[i+1][j+1] = best

    return (f, sources)
#     local_max(f, sources, seq1, seq2)
#     mat_print(f, seq1, seq2)
#     mat_print(sources, seq1, seq2)

def align_sequences(fasta1, fasta2, align_type, matrix, gap):
    # Determine all alignments of pairings of sequences from fasta1 with those from fasta2
    for seq1 in parse_fasta(fasta1):
        for seq2 in parse_fasta(fasta2):
            if (align_type == "both" or align_type == "global"): 
                print(f"GLOBAL left-first alignment for {seq1[0][:-1]} and {seq2[0][:-1]}:")
                
                (f, sources) = needleman_wunsch(seq1[1], seq2[1], matrix, gap)
                (score, identity, seq1_str, seq2_str) = global_max(f, sources, seq1[1], seq2[1])
                
                print(f'Alignment score = {score}')
                print(f'Percent identity = {identity / len(seq1_str) * 100:.2f}%')
                print(seq1_str)
                print(seq2_str)
                
                mat_print(f, seq1[1], seq2[1])
                
            if (align_type == "both"):
                print()
                
            if (align_type == "both" or align_type == "local"):
                
                print(f"LOCAL left-first alignment for {seq1[0][:-1]} and {seq2[0][:-1]}:")
                
                (f, sources) = smith_waterman(seq1[1], seq2[1], matrix, gap)
                (score, identity, seq1_str, seq2_str) = local_max(f, sources, seq1[1], seq2[1])

                print(f'Alignment score = {score}')
                print(f'Percent identity = {identity / len(seq1_str) * 100:.2f}%')
                print(seq1_str)
                print(seq2_str)

                mat_print(f, seq1[1], seq2[1])
            print('---------------------------------------------------------------------------')
            
def global_results(seq1, seq2, matrix, gap):
    (f, sources) = needleman_wunsch(seq1, seq2, matrix, gap)
    return global_max(f, sources, seq1, seq2)

def local_results(seq1, seq2, matrix, gap):
    (f, sources) = smith_waterman(seq1, seq2, matrix, gap)
    return local_max(f, sources, seq1, seq2)

def global_score(seq1, seq2, matrix, gap):
    results = global_results(seq1, seq2, matrix, gap)
    return results[0]

def local_score(seq1, seq2, matrix, gap):
    results = local_results(seq1, seq2, matrix, gap)
    return results[0]

# align.py --seq1 sequence_file.fa --seq2 different_sequence_file.fa --type global --matrix blosum62 --gap -2
def main():
    args = parser.parse_args()
    
    if args.matrix[0] == 'pam250':
        matrix = pam250
    elif args.matrix[0] == 'blosum62':
        matrix = blosum62
    elif args.matrix[0] == 'simple':
        matrix = simple_matrix
    else:
        print("Invalid matrix. Enter pam250, blosum62 or simple")
        return
    
    if args.type[0] not in ['local', 'global', 'both']:
        print("Invalid alignment type. Enter local, global or both")
        return
    
    if not exists(args.seq1[0]):
        print(args.seq1 + ' does not exist.')
        return
    if not exists(args.seq2[0]):
        print(args.seq2 + ' does not exist.')
        return
    
    align_sequences(args.seq1[0], args.seq2[0], args.type[0], matrix, args.gap[0])

if __name__ == "__main__":
    main()