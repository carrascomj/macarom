"""Write our csv sequences to a fasta."""
import pandas as pd
import typer


def to_fasta(csv_path: str, fasta_path: str, imodulon: str = "AtoC"):
    """Write a csv to fasta."""
    df: pd.DataFrame = pd.read_csv(csv_path)
    df = df.loc[df.imodulon.str.contains(imodulon), ["geneid", "seq"]]
    with open(fasta_path, "w") as f:
        f.write(
            "\n".join(
                [
                    f">{gene_id}\n{seq}"
                    for _, gene_id, seq in df.itertuples()
                ]
            )
        )

if __name__ == "__main__":
    typer.run(to_fasta)
