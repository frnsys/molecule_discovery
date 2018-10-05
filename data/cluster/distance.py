"""
Defines RV coefficient as basis for similarity, as described in:

Botev, Viktor, Kaloyan Marinov, and Florian Schäfer. "Word importance-based similarity of documents metric (WISDM):
    Fast and scalable document similarity metric for analysis of scientific documents."
    Proceedings of the 6th International Workshop on Mining Scientific Publications. ACM, 2017.
"""

import math
import numpy as np

SQRT2 = math.sqrt(2)

def RV(X, Y):
    """RV coefficient, as described in:

    Josse, Julie, Jérome Pagès, and François Husson. "Testing the significance of the RV coefficient."
        Computational Statistics & Data Analysis 53.1 (2008): 82-91."""

    S_11 = np.matmul(X,X.T)
    S_22 = np.matmul(Y,Y.T)
    S_12 = np.matmul(X,Y.T)

    n = np.trace(S_12)
    d = np.trace(S_11) * np.trace(S_22)
    return n/math.sqrt(d)


def dist(X, Y):
    """Convert RV coefficient to dissimilarity, as described in:

    Josse, Julie, Jérome Pagès, and François Husson. "Testing the significance of the RV coefficient."
        Computational Statistics & Data Analysis 53.1 (2008): 82-91."""

    rv = RV(X,Y)

    # Sometimes due to floating pt precision,
    # this is slightly above 1
    rv = min(rv, 1)

    return SQRT2 * math.sqrt(1 - rv)


if __name__ == '__main__':
    a = np.random.random((8,5))
    b = np.random.random((12,5))
    c = np.vstack([
        [1,2,3,8,8],
        [4,5,6,8,8],
        [7,8,9,8,8],
        [8,8,8,8,8]
    ])
    print(dist(a,b))
    print(dist(a,a))
    print(dist(a,c))
