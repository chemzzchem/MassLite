import numpy as np
from math import exp, sqrt, pi, log10
import time

from scipy.cluster.hierarchy import ward, fcluster
from scipy.spatial.distance import pdist

'''X = [[0, 0], [0, 1], [1, 0],
     [0, 4], [0, 3], [1, 4],
     [4, 0], [3, 0], [4, 1],
     [4, 4], [3, 4], [4, 3]]'''


X = [[0,0],[1,0],[3,0],[10,0]]

print(X)

a = pdist(X)

Z = ward(a)

'''a1 = np.array(range(50,2001))

slope1 = 1/(log10(1.000001))
intercept1 = -(log10(50)/log10(1.000001))

a2 = np.log10(a1)

a3 = a2 * slope1 + intercept1

a4 = np.zeros(len(a3))

a5 = np.stack((a3,a4),axis = 1)

a6 = pdist(a5)

a7 = ward(a6)

a8 = fcluster(a7,t=5,criterion = 'distance')'''

print("done")