

import numpy as np
import json

def f2k(Imod):
    
    def load_matrix(imod_json):
        with open(imod_json) as f:
            return json.load(f)
    
    
    #def matrix_to_dict(matrix: list[list[float]]) -> list[dict[str, float]]:
    #    # A; T; G; C
    #    return [{k: v for k, v in zip("ATGC", pos)} for pos in matrix]
    #
    def matrix_to_dict(matrix):
       # A; T; G; C
       return [{k: v for k, v in zip("ATGC", pos)} for pos in matrix]
    Mat = load_matrix("eval/init_eval.json")[Imod]
    
    #print(Mat)
    
    alphabet = np.array( ["A","T","G","C"])
    NTscoring = np.array( [[1,1/10,1/10,1/10],[1/10,1,1/10,1/10],[1/10,1/10,1,1/10],[1/10,1/10,1/10,1]] )
    
    bg = {}

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
    
  
    bg = dict(zip(alphabet, _bg))
 
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
    
        
        #Calculate f_matrix from c_matrix
        f_matrix = np.array(c_matrix)/sum(c_matrix[0])
        f_matrix = f_matrix.tolist()
 
        ## Log Odds Weight Matrix (w)
        p_matrix = f_matrix
        w_matrix = (np.array(f_matrix)*0).tolist()
        for position in range(0, core_len):
    
            for letter in range(0, len(alphabet) ):
                if p_matrix[position][letter] != 0:
                    #print( np.log(fuck1)/np.log(2) )
                    w_matrix[position][letter] = np.log( p_matrix[position][letter] / bg[ alphabet[letter] ] ) / np.log(2)
        
        #print( len( w_matrix[0]) )
        # Calculate the overall score of the peptides to the LO matrix
        _sum = 0
        for position in range(0, core_len):
            for letter in range(0, len(alphabet)):
                #fuck1 = w_matrix[position][letter] 
                _sum += f_matrix[position][letter] * w_matrix[position][letter]
        #print(w_matrix)
        return w_matrix#, _sum, p_matrix
    w_matrix = get_log_odds()
    
    #w_matrix, _sum, p_matrix = get_log_odds(Mat)
    #print(matrix_to_dict(w_matrix))
    return matrix_to_dict(w_matrix)