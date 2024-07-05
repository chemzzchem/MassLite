import numpy as np
import copy
from math import exp, sqrt, pi,log10
import time

from pyimzml.ImzMLParser import ImzMLParser

global_slope = 1/(log10(1.000001))
global_intercept = -(log10(50)/log10(1.000001))

def read_msi_data(p,thresh = 3000):
    a = max(p.coordinates)
    xmax,ymax,zmax = a
    for idx,(x,y,z) in enumerate(p.coordinates):
        data = np.transpose(p.getspectrum(idx))
        print("ok")
    print("ok")

def combine_peaks(data):

    print("done")

def rela_transform(data):
    if isinstance(data, (float,int)):
        # If the input is an integer, perform some operations for single integer input
        result = log10(data) * global_slope + global_intercept  # Example operation, you can replace this with your desired logic
    elif isinstance(data, (list, np.ndarray)):
        # If the input is a list or a NumPy array, perform some operations for array input
        arr = np.array(data)  # Convert the input to a NumPy array for uniform processing
        if arr.ndim == 1:
            # If the array is 1-dimensional, perform operations for 1D array
            result = np.log10(arr)* global_slope + global_intercept  # Example operation, you can replace this with your desired logic
        elif arr.ndim == 2:
            n_col = arr.shape[1]
            if n_col == 2:
                mz = arr[:,0]
                result = np.log10(mz)* global_slope + global_intercept
            
            elif n_col >= 3:
                raise ValueError("Has more than 2 columns")    

            # If the array is 2-dimensional, perform operations for 2D array
        else:
            raise ValueError("Input array must have either 1 or 2 dimensions.")
    else:
        raise ValueError("Input must be either an integer, a list, or a NumPy array.")

    return result

def rela_transform_rev(data):
    if isinstance(data, (float,int)):
        result = pow(10,((data-global_intercept)/global_slope))
    elif isinstance(data, (list, np.ndarray)):
        # If the input is a list or a NumPy array, perform some operations for array input
        arr = np.array(data)  # Convert the input to a NumPy array for uniform processing
        if arr.ndim == 1:
            # If the array is 1-dimensional, perform operations for 1D array
            result = np.power(10,((data-global_intercept)/global_slope))
        elif arr.ndim == 2:
            n_col = arr.shape[1]
            result = arr.copy()

            rev_mz = arr[:,min(n_col-2,1)]
            mz = np.power(10,((rev_mz-global_intercept)/global_slope))
            result[:,min(n_col-2,1)]=mz
 

            # If the array is 2-dimensional, perform operations for 2D array
        else:
            raise ValueError("Input array must have either 1 or 2 dimensions.")
    else:
        raise ValueError("Input must be either an integer, a list, or a NumPy array.")

    return result

def weighted_avg (a,b):
    return(a[0]*a[1]+)

def three_point_2 (a,b,c):
    if a[1]==0 and c[1]==0:
        d = b
    elif a[1]==0:
# main start from here

st_time = time.time()

file_path = '4_row.imzML'
ms_file = ImzMLParser(filename=file_path)

p0 = read_msi_data(ms_file)

print("end")