import numpy as np
import json
#import matplotlib.pyplot as plt
#from time import time

from argparse import ArgumentParser
#import typer
parser = ArgumentParser(description="Shittyscorer")
parser.add_argument("-p", action="store", dest="PEPTIDE", type=str, help="Seq of peptide")
parser.add_argument("-mat", action="store", dest="mat_file", type=str, help="File with PSSM")
args = parser.parse_args()
PEPTIDE = args.PEPTIDE
mat_file = args.mat_file

#print(MAT)

#_matrix = np.loadtxt(mat_file, dtype=str).tolist()

#with open(mat_file) as f:
#    mat = json.load(mat_file)
#
#
#print(_matrix)
#
#import json
#x = list(map(str.strip, json.loads( _matrix )))
# ['A', 'B', 'C', 'D']

#print(x)
PEPTIDE = "GCCTGCGCG"
_matrix = [{'A': 0.0, 'T': 0.0, 'G': 1.9770995978899213, 'C': 0.3921370971687649}, {'A': 0.438307278601691, 'T': 0.438307278601691, 'G': 0.3921370971687649, 'C': 0.3921370971687649}, {'A': 0.438307278601691, 'T': 0.438307278601691, 'G': 0.3921370971687649, 'C': 1.9770995978899213}, {'A': 0.438307278601691, 'T': 0.438307278601691, 'G': 0.3921370971687649, 'C': 0.3921370971687649}, {'A': 1.438307278601691, 'T': 0.0, 'G': 0.3921370971687649, 'C': 0.3921370971687649}, {'A': 0.0, 'T': 0.438307278601691, 'G': 0.3921370971687649, 'C': 1.3921370971687648}, {'A': 0.0, 'T': 0.0, 'G': 1.9770995978899213, 'C': 0.3921370971687649}, {'A': 0.438307278601691, 'T': 0.0, 'G': 0.3921370971687649, 'C': 0.3921370971687649}, {'A': 0.438307278601691, 'T': 0.438307278601691, 'G': 0.3921370971687649, 'C': 0.3921370971687649}]

#print( type(matrix1) )

def score_peptide(
    peptide: str,
    matrix: list):
    acum = 0
    hey = [0]
    if( len(peptide) == len(matrix) ):
        for i in range(0, len(matrix) ):
            acum += matrix[i][peptide[i]]
    elif( len(peptide) > len(matrix) ):
        a = len(peptide)-len(matrix)
        for j in range(0, a+1):
            for i in range(0, len(matrix) ):
                acum += matrix[i][peptide[i+j]]
            hey.append(acum)
            acum=0
        acum = max(hey)
    elif( len(matrix) > len(peptide) ):
        a = len(matrix)-len(peptide)
        for j in range(0, a+1,1):
            #print(matrix[j])
            for i in range(0, len(matrix)-a ):
                acum += matrix[i+j][peptide[i]]
            hey.append(acum)
            acum=0
        acum = max(hey)
        #print(hey)
    return acum

print( score_peptide(PEPTIDE,_matrix) )