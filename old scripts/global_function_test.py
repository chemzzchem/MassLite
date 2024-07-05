from enum import IntEnum
from matplotlib.style import use
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import exp, sqrt, pi
import sys, glob, os
import csv
import time

'''def useless(i):
    global a
    a[i] += 1
    return a[i]

a = np.array(range(5))

print(useless(1))'''

N = 2
a = 0
b = np.array(range(10**N))

st = time.time()

for i in range(len(b)):
    if(b[i]!=0):
        a += 1

en1 = time.time()
print("time1 is", en1 - st)

for smb in b:
    if (smb!=0):
        a+=1

en2 = time.time()
print("time2 is", en2 - en1)



print("done")