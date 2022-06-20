## This script reads the configuration files from the different simulations.
## It outputs a .csv with the algorithm id and the configurations
import pandas as pd

TEST_FILE = "data/process/output_1655478081.5168605/config"

def read_config(config_file: str)-> pd.DataFrame:
    config = pd.read_csv(config_file, sep='=', header=None)
    config.columns = ['parameter', 'value']
    return config

print(read_config(TEST_FILE))

