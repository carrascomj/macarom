from gibbssampler_function import run
import pandas as pd
import os

# load the sequence data
# find imodulons in curated list
# loop over imodulons
# save pssm matrixes

sequence_data_file = 'data/process/gene_seqs.csv'
eval_pssm_file = "eval/init_eval.json"

eval_pssm = pd.read_json(eval_pssm_file)
tf_list = eval_pssm.columns.to_list()
output_dir = "data/process/gibbssampler1_results/"

if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

for tf in tf_list[3:]:
    output_file = output_dir + "_" + str(tf)

    print(tf)
    try:
        run(sequences_csv=sequence_data_file, imodulon = tf, pssm_json=output_file)
    except IndexError:
        print("##################### " + str(tf) + " is not found in the iModulons")
