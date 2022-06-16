import json
import sys
from time import time
from typing import Optional

import numpy as np
import pandas as pd
import typer


def initialize_matrix(core_len, alphabet):
    init_matrix = [0] * core_len
    for i in range(0, core_len):
        row = {letter: 0 for letter in alphabet}
        init_matrix[i] = row
    return init_matrix


def put_to_zero(matrix):
    for i in range(0, len(matrix)):
        for key in matrix[i].keys():
            matrix[i][key] = 0.0
    return matrix


def get_log_odds(
    peptides,
    alphabet,
    bg,
    scoring_scheme,
    core_len,
    c_matrix,
    f_matrix,
    g_matrix,
    p_matrix,
    w_matrix,
    sequence_weighting,
    beta,
):
    # Amino Acid Count Matrix (c)
    c_matrix = put_to_zero(c_matrix)

    for position in range(0, core_len):

        # peptides has two elements; element[0] is the peptide sequence, element [1] is the core location
        for element in peptides:

            peptide = element[0]

            core_start = element[1]
            peptide, core_start = adjust_to_complementary(peptide, core_start)

            c_matrix[position][peptide[core_start + position]] += 1

    # Sequence Weighting
    weights = {}

    for element in peptides:

        peptide = element[0]
        core_start = element[1]
        peptide, core_start = adjust_to_complementary(peptide, core_start)

        # apply sequence weighting
        if sequence_weighting:
            w = 0.0
            neff = 0.0

            for position in range(0, core_len):
                r = 0

                for letter in alphabet:
                    if c_matrix[position][letter] != 0:
                        r += 1

                s = c_matrix[position][peptide[core_start + position]]
                w += 1.0 / (r * s)
                neff += r

            neff = neff / core_len

        # do not apply sequence weighting
        else:
            w = 1
            neff = len(peptides)

        weights[peptide] = w

    # Observed Frequencies Matrix (f)
    f_matrix = put_to_zero(f_matrix)

    for position in range(0, core_len):
        n = 0

        for element in peptides:
            peptide = element[0]
            core_start = element[1]
            peptide, core_start = adjust_to_complementary(peptide, core_start)
            f_matrix[position][peptide[core_start + position]] += weights[peptide]
            n += weights[peptide]

        for letter in alphabet:
            f_matrix[position][letter] = f_matrix[position][letter] / n

    # Pseudo Frequencies Matrix (g)
    g_matrix = put_to_zero(g_matrix)

    for position in range(0, core_len):
        for letter_1 in alphabet:
            for letter_2 in alphabet:
                g_matrix[position][letter_1] += (
                    f_matrix[position][letter_2] * scoring_scheme[letter_1][letter_2]
                )

    # Combined Frequencies Matrix (p)

    alpha = neff - 1

    for position in range(0, core_len):

        for letter in alphabet:

            num = alpha * f_matrix[position][letter] + beta * g_matrix[position][letter]
            den = alpha + beta
            p_matrix[position][letter] = num / den

    # Log Odds Weight Matrix (w)
    for position in range(0, core_len):

        for letter in alphabet:

            if p_matrix[position][letter] != 0:

                w_matrix[position][letter] = np.log(
                    p_matrix[position][letter] / bg[letter]
                ) / np.log(2)

    # Calculate the overall score of the peptides to the LO matrix
    _sum = 0
    for position in range(0, core_len):
        for letter in alphabet:
            _sum += f_matrix[position][letter] * w_matrix[position][letter]

    return w_matrix, _sum, p_matrix


def score_peptide(peptide, core_start, core_len, matrix):
    acum = 0
    for i in range(0, core_len):
        acum += matrix[i][peptide[i + core_start]]
    return acum


def adjust_to_complementary(candidate: str, start: int) -> tuple[str, int]:
    if start >= len(candidate) - start:
        return candidate[::-1], start - len(candidate)
    return candidate, start


def make_background_freq_vector(GC):
    #  ‘’'A T G C’’'
    A = (1 - (GC)) / 2
    T = (1 - (GC)) / 2
    G = GC / 2
    C = GC / 2
    return np.array([A, T, G, C])


def load_peptide_data(peptides_list: list[str], core_len: int):
    raw_peptides = [x.upper() for x in peptides_list]
    # only keep peptides with length equal to or longer than core_len
    peptides = []
    for i in range(0, len(raw_peptides)):
        if len(raw_peptides[i]) >= core_len:
            peptides.append(raw_peptides[i])
        else:
            print("Peptide length too short discard", raw_peptides[i])

    peptides = sorted(peptides, key=len)
    min_pep_len = len(peptides[0])
    max_pep_len = len(peptides[-1])

    # random core start
    np.random.shuffle(peptides)
    cores_start = [0] * len(peptides)

    for i in range(0, len(cores_start)):
        if len(peptides[i]) != core_len:
            min_core_start = 0
            # Note bug in code, corrected to +1, 9/6/2021
            max_core_start = len(peptides[i]) - core_len + 1
            cores_start[i] = np.random.randint(min_core_start, max_core_start)

    peptides = list(zip(peptides, cores_start))

    return peptides, min_pep_len, core_len


def gibbs_sampler_dna(
    peptides_list: list,
    alphabet: list,
    NTscoring: np.ndarray,
    GC_content: float,
    beta=0.0,
    iters_per_point=6,
    seed=1,
    T_i=1.0,
    T_f=0.0001,
    T_steps=5,
    sequence_weighting=False,
):

    # create "blosom" matrix
    blosum62 = {}
    for i, letter_1 in enumerate(alphabet):

        blosum62[letter_1] = {}

        for j, letter_2 in enumerate(alphabet):

            blosum62[letter_1][letter_2] = NTscoring[i, j]

    # prepare background frequencies
    _bg = make_background_freq_vector(GC_content)
    bg = dict(zip(alphabet, _bg))

    # print(bg)

    # ## Calculate log-odd matrix from a given peptide core alignment

    np.random.seed(seed)

    core_len = 9

    # get peptides
    peptides, _, core_len = load_peptide_data(peptides_list, core_len)

    print(peptides)

    T = np.linspace(T_i, T_f, T_steps)

    # Define number of iterations per temperature step
    iters = len(peptides) * iters_per_point

    print("# beta:", beta)
    print("# Sequence weighting:", sequence_weighting)
    print("# iters_per_point:", iters_per_point)
    print("# Seed:", seed)
    print("# Temperature:", T)

    # initialize matrices
    c_matrix = initialize_matrix(core_len, alphabet)
    f_matrix = initialize_matrix(core_len, alphabet)
    g_matrix = initialize_matrix(core_len, alphabet)
    p_matrix = initialize_matrix(core_len, alphabet)
    w_matrix = initialize_matrix(core_len, alphabet)

    # get initial log-odds matrix
    log_odds_matrix, peptide_scores, _ = get_log_odds(
        peptides,
        alphabet,
        bg,
        blosum62,
        core_len,
        c_matrix,
        f_matrix,
        g_matrix,
        p_matrix,
        w_matrix,
        sequence_weighting,
        beta,
    )

    # initial kld score
    print("Initial KLD score: " + str(peptide_scores))
    kld = []
    kld.append(peptide_scores)

    # other stuff
    t0 = time()
    debug = False
    # debug = True

    for t in T:

        for i in range(0, iters):

            # extract peptide
            rand_index = np.random.randint(0, len(peptides) - 1)
            peptide = peptides[rand_index][0]
            core_start_original = peptides[rand_index][1]

            if len(peptide) != core_len:

                # Bug in code, missing +1, added 09062021
                max_core_start = len(peptide) * 2 - core_len + 1
                # Maybe add check to sure core_start_shifted != core_start_original ?
                core_start_shifted = np.random.randint(0, max_core_start)

                n = len(peptide)
                while core_start_shifted == core_start_original and (
                    core_start_shifted < n - core_len
                ):
                    core_start_shifted = np.random.randint(0, max_core_start)
                prev_candidate, prev_start = adjust_to_complementary(
                    peptide, core_start_original
                )
                candidate, start = adjust_to_complementary(peptide, core_start_shifted)
                # remove peptide from list
                peptides.remove(peptides[rand_index])

                # get base log_odds
                log_odds_matrix, peptide_scores, p_matrix = get_log_odds(
                    peptides,
                    alphabet,
                    bg,
                    blosum62,
                    core_len,
                    c_matrix,
                    f_matrix,
                    g_matrix,
                    p_matrix,
                    w_matrix,
                    sequence_weighting,
                    beta,
                )

                # score peptide against log_odds
                e_original = score_peptide(
                    prev_candidate, prev_start, core_len, log_odds_matrix
                )
                if debug:
                    print("Energy before shifting: " + str(e_original))

                # score shifted peptide against log_odds
                # e_shift = XX
                e_shift = score_peptide(candidate, start, core_len, log_odds_matrix)
                if debug:
                    print("Energy after shifting: " + str(e_shift))

                # energy differential

                # de_shift = XX
                de = e_shift - e_original

                # probability of accepting move
                if de > 0:
                    p = 1
                else:
                    # e_shift = XX
                    p = np.exp(de / t)

                # throw coin
                coin = np.random.uniform(0.0, 1.0, 1)[0]

                # if coin XX p:
                if coin < p:
                    if debug:
                        print("RNG < P, Move accepted")
                    peptides.append((peptide, core_start_shifted))
                    kld.append(peptide_scores)

                else:
                    if debug:
                        print("RNG >= P, Move rejected")
                    peptides.append((peptide, core_start_original))

            else:
                if debug:
                    print("Can't shift peptide, it is a " + str(core_len) + "mer")

        print("KLD score t: " + str(t) + " KLD: " + str(peptide_scores))

    t1 = time()

    print("Time elapsed (m):", (t1 - t0) / 60)

    alignment = pd.DataFrame(peptides, columns=["sequence", "core_start"])
    alignment["core_length"] = core_len

    return log_odds_matrix, alignment


def show_aligned_cores(peptide, start, core_len):
    peptide, start = adjust_to_complementary(peptide, start)
    return peptide[start : start + core_len]


def run(
    sequences_csv: str,
    imodulon: str = "AtoC",
    sequence_weighting: bool = False,
    pssm_json: Optional[str] = None,
):
    """Fit an alignment with a core length in an iModulon."""
    if pssm_json is None:
        pssm_json = f"{imodulon}.json"
    dat_all = pd.read_csv(sequences_csv)
    pep_df = dat_all.loc[dat_all.imodulon == imodulon, "seq"]
    peptides_list = pep_df.to_list()
    alphabet = np.array(["A", "T", "G", "C"])
    NTscoring = np.array(
        [[1, -1, -1, -1], [-1, 1, -1, -1], [-1, -1, 1, -1], [-1, -1, -1, 1]]
    )
    GC_content = 0.508

    log_odds, df = gibbs_sampler_dna(
        peptides_list,
        alphabet,
        NTscoring,
        GC_content,
        T_steps=50,
        iters_per_point=50,
        sequence_weighting=sequence_weighting,
    )
    print(log_odds)
    with open(pssm_json, "w") as f:
        json.dump(log_odds, f)
    df.to_csv(sys.stdout)
    print(
        df.apply(
            lambda x: show_aligned_cores(
                x["sequence"], x["core_start"], x["core_length"]
            ),
            axis=1,
        )
    )


if __name__ == "__main__":
    typer.run(run)
