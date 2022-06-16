import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import typer

data_file = "data/benchmark_plot_test_data.csv"


def benchmark_boxplot(data_file: str, output_file: str):
    dat = pd.read_csv(data_file)
    sns.violinplot(data=dat, x = 'algorithm', y = 'pssm_pssm', alpha = 0.6)
    sns.swarmplot(data=dat, x = 'algorithm', y = 'pssm_pssm', hue='imodulon')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_file,  bbox_inches='tight')

def run(data_file, output_file):
    benchmark_boxplot(data_file, output_file)
    print("Saved boxplot to " + output_file)

if __name__ == "__main__":
    typer.run(run)
