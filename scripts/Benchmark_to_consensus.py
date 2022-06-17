from PSSM_compare import PSSM_comp
from freq2kullbackleibler import f2k
import pandas as pd
import os
import numpy as np
import csv

#sequence_data_file = 'data/process/gene_seqs.csv'
algo_dir = "data/process/"
eval_pssm_file = "eval/init_eval.json"
#test_pssm_dir = "data/process/gibbssampler1_results/"
#For each imodulon pssm in test_pssm_dir get filenames 

#shameless from stackoverflow
algo_dirs = [f for f in os.listdir(algo_dir) if os.path.isdir(os.path.join(algo_dir, f))]
#print( ( algo_dirs ) )

#for w in algo_dirs:
#        test_pssm_dir =  os.path.join(algo_dir,w,"" )
#        print(test_pssm_dir)
output = []
def b2c():
    for w in algo_dirs:
        test_pssm_dir =  os.path.join(algo_dir,w,"" )
        test_pssm_files = [f for f in os.listdir(test_pssm_dir) if os.path.isfile(os.path.join(test_pssm_dir, f)) ]
        
        real_test_pssm_files = []
        # workaround for many files in pssm dir
        for i in test_pssm_files:
            if not i.find("pssm_"):
                real_test_pssm_files.append(i)
                #print("yay")
        
        #print(real_test_pssm_files)
        test_pssm_files = real_test_pssm_files
        
        test = []
        for i in test_pssm_files:
            test.append(os.path.join(test_pssm_dir, i))
        
        test_pssm = {}
        for i in range(0, len(test)):
            test_pssm.update( { test_pssm_files[i].replace("pssm_",""): pd.read_json(test[i])}  )
            
        a = list(test_pssm.keys())
        
        
        alphabet = ["A","T","G","C"]
                        
        def test2list(mat):
            _hey2 = []
            for i in range(0, len(mat) ):
                _hey = {}
                for j in alphabet:
                    _hey.update( { j : list(pd.DataFrame.to_dict( mat )[j].values())[i] } )
                _hey2.append(_hey)
            return _hey2
        
        
        #print(f2k(a[0]))
        #print("###")
        #print(test2list(test_pssm[a[0]]))
        
        for i in range(0 , len(test_pssm)):
            output.append( [ a[i] , w , PSSM_comp( test2list( test_pssm[a[i]] ) ,  f2k( a[i] ) ) ] )
    return output

#test = b2c() 

#print( str( test[0] ).replace("[","").replace("]","")  )

header = ["Imod","Algorithm","PSSM-PSSM" ]

#with open("testmik.csv","w") as f:
#    writer = csv.writer(f)
#    writer.writerow(header)
#    writer.writerows(test)
#
def lol_to_csv(lol, out_csv):
    with open(out_csv,"w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(lol)

#def lol_to_csv(lol, out_csv):
#    with open(out_csv, "w") as f:
#        f.write("\n".join([",".join(row) for row in lol]))
#
lol_to_csv( b2c() , "testmik2.txt" )
#
#print( b2c() )
