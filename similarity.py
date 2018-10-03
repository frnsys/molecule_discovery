import math
import numpy as np

SQRT2 = math.sqrt(2)

# Using RV coefficient, as described in the WISDM paper
def RV(X, Y):
    # <https://www.utdallas.edu/~herve/Abdi-RV2007-pretty.pdf>
    S_11 = np.matmul(X,X.T)
    S_22 = np.matmul(Y,Y.T)
    S_12 = np.matmul(X,Y.T)

    # S_21 = np.matmul(Y,X.T)
    # n = np.trace(np.dot(S_12, S_21))
    # d = np.trace(S_11**2) * np.trace(S_22**2)

    n = np.trace(S_12)
    d = np.trace(S_11) * np.trace(S_22)
    return n/math.sqrt(d)


def sim(X, Y):
    rv = RV(X,Y)
    # Convert to dissimilarity
    # <https://www.researchgate.net/profile/Francois_Husson/publication/255569780_Testing_the_signiflcance_of_the_RV_coefficient/links/5548e4250cf205bce7abfba0/Testing-the-signiflcance-of-the-RV-coefficient.pdf>
    dis = SQRT2 * math.sqrt(1 - rv)
    return 1 - dis


if __name__ == '__main__':
    a = np.random.random((5,8)).T
    b = np.random.random((5,12)).T
    print(sim(a,b))
    print(sim(a,a))
