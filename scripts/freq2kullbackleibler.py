

import numpy as np
import json
#import matplotlib.pyplot as plt
#from time import time

from argparse import ArgumentParser
#import typer
#parser = ArgumentParser(description="Shittyscorer")
#parser.add_argument("-p", action="store", dest="PEPTIDE", type=str, help="Seq of peptide")
#parser.add_argument("-mat", action="store", dest="mat_file", type=str, help="File with PSSM")
#args = parser.parse_args()
#PEPTIDE = args.PEPTIDE
#mat_file = args.mat_file

Mat = {"A" : [0, 0,	1,	13,	8,	10,	4,	8,	6,	10,	8,	6,	3,	0,	0,	17,	0,	13,	0,	7,	5],
        "C"	: [4,	0,	0,	3,	4,	2,	0,	6,	4,	2,	2,	4,	1,	0,	1,	2,	19,	0,	10,	5,	3],
        "G"	: [0,	19,	3,	0,	0,	1,	7,	2,	3,	2,	1,	7,	1,	1,	2,	0,	0,	4,	5,	2,	3],
        "T"	: [15,	0,	15,	3,	7,	6,	8,	3,	6,	5,	8,	2,	14,	18,	16,	0,	0,	2,	4,	5,	8]}
print(Mat["A"] )

