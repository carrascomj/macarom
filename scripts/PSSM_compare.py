
import numpy as np

# both PSSM_test and PSSM_eval should be bare matrices
#PSSM_test = [{'A': -0.298658315564515, 'T': -0.298658315564515, 'G': 1.9770995978899213, 'C': -0.34482849699744106}, {'A': -0.298658315564515, 'T': -0.298658315564515, 'G': 0.6551715030025589, 'C': 1.240134003723715}, {'A': -0.298658315564515, 'T': -0.298658315564515, 'G': 1.6551715030025589, 'C': -0.34482849699744106}, {'A': -0.298658315564515, 'T': -0.298658315564515, 'G': -0.34482849699744106, 'C': 1.9770995978899213}, {'A': 2.0232697793228476, 'T': -0.298658315564515, 'G': -0.34482849699744106, 'C': -0.34482849699744106}, {'A': 1.286304185156641, 'T': -0.298658315564515, 'G': 0.6551715030025589, 'C': -0.34482849699744106}, {'A': -0.298658315564515, 'T': -0.298658315564515, 'G': -0.34482849699744106, 'C': 1.9770995978899213}, {'A': 0.7013416844354851, 'T': -0.298658315564515, 'G': 1.240134003723715, 'C': -0.34482849699744106}, {'A': -0.298658315564515, 'T': -0.298658315564515, 'G': 1.9770995978899213, 'C': -0.34482849699744106}]
#PSSM_eval = [{'A': -0.298658315564515, 'T': -0.298658315564515, 'G': 1.9770995978899213, 'C': -0.34482849699744106},{'A': -0.298658315564515, 'T': -0.298658315564515, 'G': -0.34482849699744106, 'C': 1.9770995978899213}, {'A': -0.298658315564515, 'T': -0.298658315564515, 'G': 1.240134003723715, 'C': -0.34482849699744106}, {'A': -0.298658315564515, 'T': -0.298658315564515, 'G': 1.240134003723715, 'C': 0.6551715030025589}, {'A': -0.298658315564515, 'T': -0.298658315564515, 'G': -0.34482849699744106, 'C': 1.6551715030025589}, {'A': -0.298658315564515, 'T': -0.298658315564515, 'G': 1.6551715030025589, 'C': -0.34482849699744106}, {'A': -0.298658315564515, 'T': -0.298658315564515, 'G': -0.34482849699744106, 'C': 1.6551715030025589}, {'A': -0.298658315564515, 'T': 0.7013416844354851, 'G': -0.34482849699744106, 'C': -0.34482849699744106}, {'A': -0.298658315564515, 'T': -0.298658315564515, 'G': 1.240134003723715, 'C': -0.34482849699744106}]

#del
#print(  )
#print(PSSM_test)

def PSSM_comp(PSSM_test, PSSM_eval):
    
    def dict2list(mat):
        _hey = []
        for i in mat:
        #print(i.items() )
            _hey.append(list(i.values()) )
        return(_hey)
    
    def make_background_freq_vector(GC):
    #  ‘’'A T G C’’'
        A = (1 - (GC)) / 2
        T = (1 - (GC)) / 2
        G = GC / 2
        C = GC / 2
        return np.array([A, T, G, C])
        
    GC_content = 0.508
    
    #print iwheri
    
    
    _bg = make_background_freq_vector(GC_content)
    alphabet = np.array(["A", "T", "G", "C"])
    bg = dict(zip(alphabet, _bg))
    
        
    #np.asarray(PSSM_test, dtype=float)
    
    def PSSM2Prob(mat):
        mat2 = (2**(np.asarray(mat,dtype=float) ))
        
        for position in range(0, len(mat)):
    
            for letter in range(0,len(alphabet)):
    
                if mat2[position][letter] != 0:
                    mat2[position][letter] = mat2[position][letter] * bg[ alphabet[letter] ]
            if np.sum( mat2[position] ) != 1.0:
            # if sum of probabilities for a given position don't add up to 1, normalize to 1
                mat2[position] = mat2[position] * 1/np.sum(mat2[position] )
                
    
        
        return mat2
    
    
    #x = dict2list(PSSM_test)
    #
    #print( x )
    #
    #x = dict2list(PSSM_test)
    #x = np.asarray(x, dtype=float)
    #
    #print( 2**x )
    ##print( np.linalg.matrix_power(x, 2))
    
    def KL(PSSM_test, PSSM_eval):
        epsilon = 1e-5
        #PSSM_test = np.asarray(PSSM_test, dtype=float)
        #PSSM_eval = np.asarray(PSSM_eval, dtype=float)
        
        #PSSM_test = np.power())
    
        PSSM_test += epsilon
        PSSM_eval += epsilon
    
        KL_value = np.sum(np.where(PSSM_test != 0, PSSM_test * np.log(PSSM_test / PSSM_eval), 0))
    
        return KL_value
    
    #Variable pssm_test and pssm_eval lengths
    #
    #print( len(PSSM_test[1:len(PSSM_test)]) )
    if( len(PSSM_test) == len(PSSM_eval) ):
        scores = ( KL( PSSM2Prob( dict2list(PSSM_test) ) ,  PSSM2Prob( dict2list(PSSM_eval) )  )  )
        #print("Equal size")
        #score = ( (scores) )
        score = ( 1-np.exp(-(scores)) )
    elif( len(PSSM_test) > len(PSSM_eval) ):
        a = len(PSSM_test) - len(PSSM_eval)
        scores = []
        for i in range(0,a+1):
            scores.append( KL( PSSM2Prob( dict2list(PSSM_test[i:(len(PSSM_test)-a+i) ] ) ) ,  PSSM2Prob( dict2list(PSSM_eval) )) )
        #score = ( min(scores) )
        score = ( 1-np.exp(-min(scores)) )
        #print("Test>Eval")
    elif( len(PSSM_eval) > len(PSSM_test) ):
        a = len(PSSM_eval) - len(PSSM_test)
        scores = []
        #print("Eval>Test")
        #print(a)
        for i in range(0,a+1):
            #print(i)
            scores.append( KL( PSSM2Prob( dict2list(PSSM_test ) ),  PSSM2Prob(dict2list(PSSM_eval[i:(len(PSSM_eval)-a+i) ] ) )  ) )
        #score = ( min(scores) )
        score = ( 1-np.exp(-min(scores)) )
    #
    return score

#print(PSSM_comp(PSSM_test, PSSM_eval))
