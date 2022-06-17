
## Include PSSM matrix from Gibss sample, name of the iModulon, and the output directory ##

import os
import logomaker
import pandas as pd
import matplotlib.pyplot as plt

def PSSMtoLOGO(PSSM, name, filepath):
    logo = pd.DataFrame({"A": [d["A"] for d in PSSM], "C": [d["C"] for d in PSSM], "G" : [d["G"] for d in PSSM], "T" : [d["T"] for d in PSSM]})
    logo.index.name = "pos"
    logo_fig = logomaker.Logo(logo,color_scheme="colorblind_safe",baseline_width=2.0,flip_below=False,vsep=0.2)
    logo_fig.style_spines(spines=["bottom","top"],visible=False)
    logo_fig.ax.xaxis.set_ticks_position('none')
    logo_fig.ax.set_ylabel("log odds")
    file = name+".png"
    logo_fig.fig.savefig(os.path.join(filepath+file),dpi=300)

#PSSM_test = [{"A": -2.8836208162856716, "T": -2.8836208162856716, "G": -2.929790997718597, "C": 1.928189997408975}, {"A": -1.298658315564515, "T": -2.8836208162856716, "G": -0.12243607566099331, "C": 1.3921370971687648}, {"A": 1.438307278601691, "T": -0.07626589422806718, "G": -2.929790997718597, "C": -1.9297909977185974}, {"A": 0.7013416844354851, "T": 0.28630418515664097, "G": -2.929790997718597, "C": 0.07020900228140264}, {"A": 1.7602353734890535, "T": -0.8836208162856716, "G": -2.929790997718597, "C": -2.929790997718597}, {"A": -0.8836208162856716, "T": 1.7013416844354852, "G": -1.9297909977185974, "C": -2.929790997718597}, {"A": -0.561692721398309, "T": 1.7013416844354852, "G": -2.929790997718597, "C": -2.929790997718597}, {"A": -2.8836208162856716, "T": 1.6399411397713415, "G": -2.929790997718597, "C": -0.6078629028312352}, {"A": -0.8836208162856716, "T": 1.438307278601691, "G": -0.6078629028312352, "C": -2.929790997718597}]
#PSSMtoLOGO(PSSM_test,"CysB","/home/shtapa/22175/Project/")
