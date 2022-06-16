

import numpy as np
import json
#import matplotlib.pyplot as plt
#from time import time

from argparse import ArgumentParser
import json
#import typer
#parser = ArgumentParser(description="Shittyscorer")
#parser.add_argument("-p", action="store", dest="PEPTIDE", type=str, help="Seq of peptide") #parser.add_argument("-mat", action="store", dest="mat_file", type=str, help="File with PSSM")
#args = parser.parse_args()
#PEPTIDE = args.PEPTIDE
#mat_file = args.mat_file

# Mat = {"A" : [0, 0,	1,	13,	8,	10,	4,	8,	6,	10,	8,	6,	3,	0,	0,	17,	0,	13,	0,	7,	5],
#         "C"	: [4,	0,	0,	3,	4,	2,	0,	6,	4,	2,	2,	4,	1,	0,	1,	2,	19,	0,	10,	5,	3],
#         "G"	: [0,	19,	3,	0,	0,	1,	7,	2,	3,	2,	1,	7,	1,	1,	2,	0,	0,	4,	5,	2,	3],
#         "T"	: [15,	0,	15,	3,	7,	6,	8,	3,	6,	5,	8,	2,	14,	18,	16,	0,	0,	2,	4,	5,	8]}
#print(Mat["A"] )

def load_matrix(imod_json):
    with open(imod_json) as f:
        return json.load(f)


def matrix_to_dict(matrix: list[list[float]]) -> list[dict[str, float]]:
    # A; T; G; C
    return [{k: v for k, v in zip("ATGC", pos)} for pos in matrix]

Mat = load_matrix("eval/init_eval.json")["MalT"]

#for i in range(0, len(Mat["A"] ))):
#print(Mat)

alphabet = np.array( ["A","T","G","C"])


NTscoring = np.array( [[1,1/10,1/10,1/10],[1/10,1,1/10,1/10],[1/10,1/10,1,1/10],[1/10,1/10,1/10,1]] )

#bg_file = data_dir + "Matrices/bg.freq.fmt"

#_bg = np.loadtxt(bg_file, dtype=float)

bg = {}
#for i in range(0, len(alphabet)):
#    bg[alphabet[i]] = _bg[i]
#
#blosum_file = data_dir + "Matrices//blosum62.freq_rownorm"
#
#_blosum62 = np.loadtxt(blosum_file, dtype=float).reshape((20, 20)).T
#
blosum62 = {}

_blosum62 = NTscoring

for i, letter_1 in enumerate(alphabet):
    
    blosum62[letter_1] = {}

    for j, letter_2 in enumerate(alphabet):
        
        blosum62[letter_1][letter_2] = _blosum62[i, j]


def make_background_freq_vector(GC):
  #  ‘’'A T G C’’'
    A = (1 - (GC))/2
    T = (1 - (GC))/2
    G = GC / 2
    C = GC / 2
    return np.array([A, T, G, C])


_bg = make_background_freq_vector(0.508)


#_bg = make_background_freq_vector(0.5)
#alphabet = np.array([“A”, ‘T’, ‘G’, ‘C’])
bg = dict(zip(alphabet, _bg))

# print(bg)

# ## Calculate log-odd matrix from a given peptide core alignment

# In[40]:

#def initialize_matrix(core_len, alphabet):
#
#    init_matrix = [0]*core_len
#
#    for i in range(0, core_len):
#
#        row = {}
#
#        for letter in alphabet: 
#            row[letter] = 0.0
#
#        init_matrix[i] = row
#        
#    return init_matrix
#
#def put_to_zero(matrix):
#    
#    for i in range(0, len(matrix)):
#
#        for key in matrix[i].keys():
#        
#            matrix[i][key] = 0.0
#    
#    return matrix
#
peptides = Mat
def get_log_odds(peptides = Mat, 
    alphabet = alphabet, 
    bg = bg, 
    scoring_scheme = blosum62, 
    core_len = len(peptides["A"]) ):

    # Amino Acid Count Matrix (c)
  #  c_matrix = initialize_matrix(core_len, alphabet)

    c_matrix = []
    for i in range(0, core_len):
        c_matrix.append( [peptides["A"][i],peptides["T"][i],peptides["G"][i],peptides["C"][i]] )

    #print(c_matrix)
    # Sequence Weighting
    #print(c_matrix)
    weights = {}
    #print(np.array(c_matrix) )
    #for element in peptides:
	#
    #    peptide = element[0]
    #    core_start = element[1]
    #        
    #    # apply sequence weighting
    #    if sequence_weighting:
	#
    #        w = 0.0
    #        neff = 0.0
	#
    #        for position in range(0, core_len):
	#
    #            r = 0
	#
    #            for letter in alphabet:        
	#
    #                if c_matrix[position][letter] != 0:
	#
    #                    r += 1
	#
    #            s = c_matrix[position][peptide[core_start+position]]
	#
    #            w += 1.0/(r * s)
	#
    #            neff += r
	#
    #        neff = neff / core_len
	#
    #    # do not apply sequence weighting
    #    else:
	#
    #        w = 1  
	#
    #        neff = len(peptides)  
	#
	#
    #    weights[peptide] = w

    
    # Observed Frequencies Matrix (f)
    #f_matrix = put_to_zero(f_matrix)
	#
    #for position in range(0, core_len):
	#
    #    n = 0;
	#
    #    for element in peptides:
	#
    #        peptide = element[0]
    #        
    #        core_start = element[1]
    #          
    #        f_matrix[position][peptide[core_start+position]] += weights[peptide]
	#
    #        n += weights[peptide]
	#
    #    for letter in alphabet: 
	#
    #        f_matrix[position][letter] = f_matrix[position][letter]/n
    
	#Calculate f_matrix from c_matrix
    f_matrix = np.array(c_matrix)/sum(c_matrix[0])
    f_matrix = f_matrix.tolist()
   # for i in 
    #print(f_matrix)
    #print(f_matrix)
    # Pseudo Frequencies Matrix (g)
    #print( len(alphabet) )
    #g_matrix = np.array(f_matrix)*0
    #w_matrix = (np.array(f_matrix)*0).tolist()
    #for position in range(0, core_len):
    #
    #    for letter_1 in range(0,len(alphabet),1 ):
    #       
    #        for letter_2 in range(0, len(alphabet),1 ):
    #           #a = f_matrix[position][letter_2] * scoring_scheme[ alphabet[letter_1] ][ alphabet[letter_2] ]
    #           #print(a)
    #           #g_matrix[position][letter_1] =+ a
    #            #print( f_matrix[position][letter_2].item() * scoring_scheme[alphabet[letter_1] ][ alphabet[letter_2] ] )
    #            g_matrix[position][letter_1] += f_matrix[position][letter_2] * scoring_scheme[ alphabet[letter_1] ][ alphabet[letter_2] ]
    #
    #
    ##print(g_matrix)
    ##print(g_matrix)
    ## Combined Frequencies Matrix (p)
    ##print(g_matrix)
    #neff = sum(c_matrix[0])
    #alpha = neff - 1
    #beta = 0
    ##for position in range(0, core_len):
	##
    ##    for letter in range(0, len(alphabet), 1):
	##
    ##        num = alpha*f_matrix[position][letter] + beta*g_matrix[position][letter]
    ##        
    #den = alpha + beta
    ##        #print(alpha)
    ##        p_matrix[position][letter] = num / den
    ##print( np.array(g_matrix)*1 )
    #num = np.array(f_matrix)*alpha + beta*np.array(g_matrix)
    #p_matrix = num/den
    #p_matrix = p_matrix.tolist()
    ##print(  p_matrix )
    ## Log Odds Weight Matrix (w)
    p_matrix = f_matrix
    w_matrix = (np.array(f_matrix)*0).tolist()
    for position in range(0, core_len):

        for letter in range(0, len(alphabet) ):
            if p_matrix[position][letter] != 0:
                #fuck1 = p_matrix[position][letter]/bg[  alphabet[letter] ]
                #print( np.log(fuck1)/np.log(2) )
                w_matrix[position][letter] = np.log( p_matrix[position][letter] / bg[ alphabet[letter] ] ) / np.log(2)
                #w_matrix[position][letter] = np.log(fuck1)/np.log(2)
    
    #print( len( w_matrix[0]) )
    # Calculate the overall score of the peptides to the LO matrix
    _sum = 0
    for position in range(0, core_len):
        for letter in range(0, len(alphabet)):
            #fuck1 = w_matrix[position][letter] 
            _sum += f_matrix[position][letter] * w_matrix[position][letter]
    #print(w_matrix)
    return w_matrix, _sum, p_matrix

w_matrix, _sum, p_matrix = get_log_odds(Mat)
print(matrix_to_dict(w_matrix))
