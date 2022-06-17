from gibbssampler_function import run
import pandas as pd
import os
from os.path import join
from time import time
from multiprocessing import Pool

# load the sequence data
# find imodulons in curated list
# loop over imodulons
# save pssm matrixes


def run_with_tf_catch(tf_kwargs_output_dir_sequences_csv):
    tf, kwargs, output_dir, sequence_csv = tf_kwargs_output_dir_sequences_csv
    try:
        run(
            sequences_csv=sequence_csv,
            imodulon=tf,
            pssm_json=output_dir + "pssm_" + str(tf),
            alignment_output=output_dir + "alignment_" + str(tf),
            **kwargs,
        )
    except IndexError:
        print("##################### " + str(tf) + " is not found in the iModulons")


def benchmark_gibbssampler(
    tf_list: str, sequence_csv: str, output_dir: str, *args, **kwargs
):
    """Run through the transcription factors in tf_list with given setting for the gibbssamples. Settings are given in the kwargs"""
    output_dir  = join(output_dir, f"output_{time()}")
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    with open(join(output_dir, "config"), "w") as fo:
        fo.write("\n".join([f"{key} = {val}" for key, val in kwargs.items()]))

    t0 = time()

    with Pool(5) as pool:
        pool.map(
            run_with_tf_catch,
            [(tf, kwargs, output_dir, sequence_csv) for tf in tf_list],
        )
        pool.close()
        pool.join()
    t1 = time()
    runtime = (t1 - t0) / 60
    with open(join(output_dir, "config"), "w") as fo:
        fo.write("\n".join([f"{key} = {val}" for key, val in kwargs.items()]))

    with open(output_dir + "runtime.txt", "w") as f:
        f.write("Runetime in min.: " + str(runtime))


################ INITIATION ###############################
# Get evaluation transcription factors
eval_pssm_file = "eval/init_eval.json"
eval_pssm = pd.read_json(eval_pssm_file)
tf_list = eval_pssm.columns.to_list()

# Get sequences data with all imodulons
sequence_data_file = "data/process/gene_seqs.csv"
output_dir = "data/process/"

################ BENCHMARK ################################
###### CONFIG 1: Fixed core length 9 nucleotides
# parameters
core_min = 9
core_max = 9
t_steps = 200
iters_per_point = 50

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
core_min = 9
core_max = 19
t_steps = 200
iters_per_point = 200

# output

benchmark_gibbssampler(
    tf_list=tf_list,
    sequence_csv=sequence_data_file,
    output_dir=output_dir,
    core_min=core_min,
    core_max=core_max,
    t_steps=t_steps,
    iters_per_point=iters_per_point,
)
