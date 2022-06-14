"""Functions to get the upstream positions for each imodulon in the S matrix."""
import sys
from os.path import join

import pandas as pd

# TODO(jorge): use pymodulon


def get_imodulon(
    name: str,
    S: pd.DataFrame,
    imodulons: pd.DataFrame,
    gene_info: pd.DataFrame,
    show_info: bool = True,
    show_regs: bool = True,
) -> pd.DataFrame:
    """Return pandas dataframe containing i-modulon genes and coefficients.

    Parameters:
        name: name of the imodulon
        S: https://github.com/SBRG/precise-db/blob/master/data/S.csv
        imodulons: https://github.com/SBRG/precise-db/blob/master/data/curated_enrichments.csv
        gene_info: https://github.com/SBRG/precise-db/blob/master/data/gene_info.csv
        show_info: Show extended information about each gene
        show_regs: Show known regulators for each gene (separated by commas)

    Adapted from https://github.com/SBRG/precise-db/
    """

    if name not in imodulons.index:
        raise ValueError(
            "{} is not a valid i-modulon name. See imodulons.index".format(name)
        )

    comp = S[name]
    thresh = imodulons.loc[name, "threshold"]
    genes = comp[abs(comp) > thresh].sort_values(ascending=False)

    if not show_info:
        df = pd.DataFrame(genes)
    else:
        df = gene_info.loc[genes.index]
        df["coefficient"] = genes.values

    return df


def prepare_imodulon(name: str, S, imodulons, gene_info) -> pd.DataFrame:
    df = get_imodulon(name, S, imodulons, gene_info)[
        ["start", "stop", "strand", "coefficient"]
    ]
    df["imodulon"] = name
    return df


def process_matrix(
    S: pd.DataFrame, imodulons: pd.DataFrame, gene_info: pd.DataFrame
) -> pd.DataFrame:
    """Prepare S matrix to be fetched from a genome (id, position)."""
    df = pd.concat(
        [prepare_imodulon(name, S, imodulons, gene_info) for name in S.columns]
    )
    df["position"] = df.apply(
        lambda x: x.start if x.strand == "+" else x.stop + 200, axis=1
    )
    df["reverse"] = df.strand.apply(lambda strand: strand == "-")
    return df[["position", "reverse", "coefficient", "imodulon"]]


if __name__ == "__main__":
    data_dir = sys.argv[1]
    gene_info = pd.read_csv(join(data_dir, "gene_info.csv"), index_col=0)
    S = pd.read_csv(join(data_dir, "S.csv"), index_col=0)
    imodulons = pd.read_csv(join(data_dir, "curated_enrichments.csv"), index_col=0)
    process_matrix(S, imodulons, gene_info).to_csv(
        sys.argv[2] if len(sys.argv) > 2 else "genes.csv"
    )
