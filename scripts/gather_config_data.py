## This script reads the configuration files from the different simulations.
## It outputs a .csv with the algorithm id and the configurations
import pandas as pd
import typer
from os import listdir
from os.path import join

# TEST_FILE = "data/process/output_1655478081.5168605/config"
# DATA_FOLDER = "data/process/"
# OUTPUT_FILE = "data/prcess/benchmark_results/config_summary.csv"

def read_config(config_file: str)-> pd.DataFrame:
    config = pd.read_csv(config_file, sep='=', header=None)
    config.columns = ['parameter', 'value']

    # remove white spaces
    config['parameter'] = config['parameter'].str.strip()
    config['value'] = config['value'].str.strip()
    
    return config

def parse_config_files(DATA_FOLDER: str)->pd.DataFrame:
    '''Finds all config files and parse them into a pandas dataframe.'''
    all_dirs = listdir(DATA_FOLDER)
    output_dirs = [d for d in all_dirs if "output_" in d and "runtime.txt" not in d]


    config_df = pd.DataFrame()
    for d in output_dirs:
        config_df_tmp = read_config(join(DATA_FOLDER, d, "config"))
        config_df_tmp["config_id"] = d
        config_df = config_df.append(config_df_tmp)

    config_df = config_df[["config_id", "parameter", "value"]]
    return config_df

def run(data_folder: str, output_file: str):
    config_df = parse_config_files(data_folder)
    config_df.to_csv(output_file)
    

if __name__ == "__main__":
    typer.run(run)


