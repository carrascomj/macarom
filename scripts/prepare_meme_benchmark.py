import json
import re
from os.path import join

import typer

POS_PAT = re.compile(r"w= (\d+)")


def get_matrix(
    meme_txt: str, matrix_name: str = "letter-probability matrix"
) -> list[dict[str, float]]:
    """Parse the log-odds matrix from meme."""
    matrix: list[dict[str, float]] = []
    with open(meme_txt) as f:
        n_pos = 0
        while n_pos == 0:
            line = f.readline()
            if not line:
                return []
            if line.startswith(matrix_name):
                matched = POS_PAT.search(line)
                if matched is not None:
                    n_pos = int(matched[1])
        for _ in range(n_pos):
            line = f.readline()
            log_odds = [num for num in line.split(" ") if num]
            matrix.append({letter: float(num) for num, letter in zip(log_odds, "ACGT")})
    return matrix


def run_benchmark(eval_json: str, meme_results: str, output_dir: str):
    with open(eval_json) as f:
        eval_tfs = json.load(f)
    for imod in eval_tfs:
        mat = get_matrix(join(meme_results, imod, "meme.txt"), "log-odds matrix")
        if mat:
            with open(join(output_dir, f"pssm_{imod}"), "w") as f:
                json.dump(mat, f)
        

if __name__ == "__main__":
    typer.run(run_benchmark)
