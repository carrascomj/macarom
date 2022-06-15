
import numpy as np

# both PSSM_test and PSSM_eval should be bare matrices

def KL(PSSM_test, PSSM_eval):
    epsilon = 1e-5
    PSSM_test = np.asarray(PSSM_test, dtype=float)
    PSSM_eval = np.asarray(PSSM_eval, dtype=float)

    PSSM_test += epsilon
    PSSM_eval += epsilon

    KL_value = np.sum(np.where(PSSM_test != 0, PSSM_test * np.log(PSSM_test / PSSM_eval), 0))

    return KL_value
