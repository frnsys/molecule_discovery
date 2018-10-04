import math
import numpy as np

SQRT2 = math.sqrt(2)

# Using RV coefficient as basis for similarity, as described in the WISDM paper
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


def dist(X, Y):
    rv = RV(X,Y)

    # Sometimes due to floating pt precision,
    # this is slightly above 1
    rv = min(rv, 1)

    # Convert to dissimilarity
    # <https://www.researchgate.net/profile/Francois_Husson/publication/255569780_Testing_the_signiflcance_of_the_RV_coefficient/links/5548e4250cf205bce7abfba0/Testing-the-signiflcance-of-the-RV-coefficient.pdf>
    return SQRT2 * math.sqrt(1 - rv)

def sim(X, Y):
    return 1 - dist(X, Y)


if __name__ == '__main__':
    a = np.random.random((8,5))
    b = np.random.random((12,5))
    c = np.vstack([
        [1,2,3,8,8],
        [4,5,6,8,8],
        [7,8,9,8,8],
        [8,8,8,8,8]
    ])
    print(sim(a,b))
    print(sim(a,a))
    print(sim(a,c))
