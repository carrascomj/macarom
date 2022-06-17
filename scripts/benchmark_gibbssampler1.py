from ctypes import alignment
from gibbssampler_function import run
import pandas as pd
import os

# load the sequence data
# find imodulons in curated list
# loop over imodulons
# save pssm matrixes


def benchmark_gibbssampler(
    tf_list: str, sequence_csv: str, output_dir: str, *args, **kwargs
):
    """Run through the transcription factors in tf_list with given setting for the gibbssamples. Settings are given in the kwargs"""
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    for tf in tf_list:
        pssm_output = output_dir + "pssm_" + str(tf)
        alignment_output = output_dir + "alignment_" + str(tf)

        print(tf)
        try:
            run(
                sequences_csv=sequence_csv,
                imodulon=tf,
                pssm_json=pssm_output,
                alignment_output=alignment_output,
                **kwargs
            )
        except IndexError:
            print("##################### " + str(tf) + " is not found in the iModulons")


################ INITIATION ###############################
# Get evaluation transcription factors
eval_pssm_file = "eval/init_eval.json"
eval_pssm = pd.read_json(eval_pssm_file)
tf_list = eval_pssm.columns.to_list()

# Get sequences data with all imodulons
sequence_data_file = "data/process/gene_seqs.csv"

################ BENCHMARK ################################
###### CONFIG 1: Fixed core length 9 nucleotides
# parameters
core_min = 9
core_max = 9
t_steps = 200
iters_per_point = 200

# output
output_dir = "data/process/gibbssampler_config1/"

benchmark_gibbssampler(
    tf_list=tf_list,
    sequence_csv=sequence_data_file,
    output_dir=output_dir,
    core_min=core_min,
    core_max=core_max,
    t_steps=t_steps,
    iters_per_point=iters_per_point,
)

###### CONFIG 2: Variable core length
# parameters
core_min = 5
core_max = 30
t_steps = 200
iters_per_point = 200

# output
output_dir = "data/process/gibbssampler_config2/"

benchmark_gibbssampler(
    tf_list=tf_list,
    sequence_csv=sequence_data_file,
    output_dir=output_dir,
    core_min=core_min,
    core_max=core_max,
    t_steps=t_steps,
    iters_per_point=iters_per_point,
)
