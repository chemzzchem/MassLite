from ctypes import alignment
from enum import IntEnum
from matplotlib.style import use
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import exp, sqrt, pi,log10
import sys, glob, os
import csv
import time
import copy

from scipy.cluster.hierarchy import ward, fcluster
from scipy.spatial.distance import pdist
from pyimzml.ImzMLParser import ImzMLParser
from scipy.signal import argrelextrema

import sklearn

def flatten(xss):
    return [x for xs in xss for x in xs]

def timing_test(kw="this"):
    global st
    en = time.time()
    print("till ",kw," step takes",en-st)

def weighted_avg(peaks,mz_col_index=0,int_col_index=1,row_index_st=0,row_index_en=-1):
    mzs = np.take(peaks[row_index_st:row_index_en],mz_col_index,axis = 1)
    intensity = np.take(peaks[row_index_st:row_index_en],int_col_index,axis = 1)
    produc = mzs * intensity
    return(sum(produc)/sum(intensity),sum(intensity))


def mass_difference(mz1,mz2,tol):
    if tol >=1:
        return((2*abs(mz1-mz2)/(mz1+mz2))<(tol*10**-6))
    else:
        return ((2*abs(mz1-mz2)/(mz1+mz2))<tol)

def new_m_int(m_old,int_old,m_add,int_add):
    m_new = (m_old*int_old+m_add*int_add)/(int_old+int_add)
    int_new = int_old+int_add
    return(m_new,int_new)

def three_point2 (x,y):
    if x[0]==0 and x[-1]==0:
        mz = x[1]
        inten = y[1]
    elif x[0]==0 or x[-1]==0:
        a = x[np.where(x>0)]
        b = y[np.where(y>0)]
        center = (a[0]*b[0]+a[1]*b[1])/(b[0]+b[1])
        a=np.append(a,center*2-x[0])
        b=np.append(b,y[0])
        z=np.polyfit(a,b,2)
        mz = -z[1]/2/z[0]
        inten = z[2]-z[1]*z[1]/4/z[0]
        if mz<a[0] or mz>a[1] or inten < np.min(b):
            mz = center
            inten = 1.125*0.5*(b[0]+b[1])
    else:
        z=np.polyfit(x,y,2)
        mz = -z[1]/2/z[0]
        inten = z[2]-z[1]*z[1]/4/z[0]
    return (mz,inten)
        
def read_msi_data(p,thresh = 3000):
    a = max(p.coordinates)
    xmax,ymax,zmax = a
    pixels = []
    for idx, (x,y,z) in enumerate(p.coordinates):
        mzs, intensities = p.getspectrum(idx)
        indexs = (intensities>thresh)
        mzs = mzs[indexs]
        intensities = intensities[indexs]
        peaks = np.stack((mzs,intensities),axis=1)
        tic = sum(intensities)
        pixel = [idx,x,y,tic]
        pixel.append(mzs)
        pixel.append(intensities)
        pixel.append(peaks)
        pixel.append(1)
        pixels.append(pixel)
    return(pixels,xmax,ymax)

def combine_msi_data2(pixel,ppm=10):
    original_mzs = pixel[4]
    intensities = pixel[5]
    l = len(original_mzs)

    replicated_length = []
    cases = np.zeros(6)
    slope1 = 1/(log10(1.000001))
    intercept1 = -(log10(50)/log10(1.000001))

    #logged_mz = np.log10(original_mzs)
    transformed_peaks = np.log10(original_mzs) * slope1 + intercept1
    #transformed_sub_from = np.concatenate((np.zeros(1),transformed_peaks[0:-1]),axis=0)
    #transformed_diff = transformed_peaks - transformed_sub_from
    transformed_diff = np.insert(np.diff(transformed_peaks),0,0)
    #print(np.where(transformed_diff>=5))
    pieces = np.where(transformed_diff>ppm)[0]+1
    sliced_mzs = np.append(np.insert(original_mzs,pieces,0),0)
    sliced_int = np.append(np.insert(intensities,pieces,0),0)

    mask = argrelextrema(sliced_int,np.greater)
    mask = mask[0]
    #real_mzs = sliced_mzs[mask]
    #real_int = sliced_int[mask]

    real_mzs = np.zeros(len(mask))
    real_int = np.zeros(len(mask))
    for i in range(len(mask)):
        a = mask[i]
        real_mzs[i],real_int[i]=three_point2(sliced_mzs[a-1:a+2],sliced_int[a-1:a+2])



    return(real_mzs,real_int,cases)

def ms_align2(pixels,ppm=10):
    pixel_indices = []

    for i in range(len(pixels)):
        pixel_indices.append(len(pixels[i][4]))
    pixel_indices = np.insert(np.cumsum(pixel_indices),0,0)
    mass_extract = np.take(pixels,4,axis = 1)
    inten_extract = np.take(pixels,5,axis = 1)
    sum_mz = np.array(flatten(mass_extract))
    sum_int = np.array(flatten(inten_extract))
    sum_rows = np.zeros(max(pixel_indices))
    
    for i in range(len(pixels)):
        st_i = pixel_indices[i]
        en_i = pixel_indices[i+1]
        sum_rows[st_i:en_i]=i
    sum_of_pixel = np.column_stack((sum_rows,sum_mz,sum_int))

    '''sum_of_pixel = np.zeros(shape=(0,3))

    for i in range(len(pixels)):
        pixel_reading = pixels[i][6]
        pixel_reading = np.concatenate((i*np.ones(shape = (len(pixel_reading),1)),pixel_reading),axis =1)
        sum_of_pixel = np.concatenate((sum_of_pixel,pixel_reading),axis = 0)
    '''
    sorted_sum = sum_of_pixel[sum_of_pixel[:,1].argsort()]
    slope1 = 1/(log10(1.000001))
    intercept1 = -(log10(50)/log10(1.000001))
   

    peaks = np.take(sorted_sum,1,axis = 1)
    logged_peaks = np.log10(peaks)
    transformed_peaks = logged_peaks * slope1 + intercept1

    # Important Efficiency Setting Here!!!!!

    pieces = np.arange(0,len(transformed_peaks),i*2,dtype=int)
    pieces[-1]=len(transformed_peaks)
    groups = []

    for i in range(len(pieces)-1):
        reading = transformed_peaks[pieces[i]:pieces[i+1]]
        logged_combined = np.stack((reading,np.zeros(pieces[i+1]-pieces[i])),axis = 1)
        dendro = ward(pdist(logged_combined))
        grouped_mzs = fcluster(dendro,t=ppm,criterion = 'distance')
        a,b = np.unique(grouped_mzs,return_index = True)
        b.sort()
        b=b+pieces[i]
        groups.append(b)

    '''print("done")
    for i in range(len(groups)-2,-1,-1):
        mz1,int1 = weighted_avg(sorted_sum,mz_col_index=1,int_col_index=2,row_index_st=groups[i][-1],row_index_en=groups[i+1][0]) 
        mz2,int2 = weighted_avg(sorted_sum,mz_col_index=1,int_col_index=2,row_index_st=groups[i+1][0],row_index_en=groups[i+1][1])
        if mass_difference(mz1,mz2,tol=10):
            groups[i+1]=groups[i+1][1:]'''

    all_groups = [x for xs in groups for x in xs]
    all_groups.append(pieces[-1])

    aligned_mzs_group = np.zeros(pieces[-1])
    aligned_mzs = np.zeros(pieces[-1])
    aligned_tot_int = np.zeros(len(all_groups)-1)
    aligned_occ = np.zeros(len(all_groups)-1)
    for i in range(len(all_groups)-1):
        st_i = int(all_groups[i])
        en_i = int(all_groups[i+1])
        aligned_mzs_group[st_i:en_i]=i
        mz,intensity = weighted_avg(sorted_sum,mz_col_index=1,int_col_index=2,row_index_st=st_i,row_index_en=en_i)
        aligned_mzs[st_i:en_i]=mz
        aligned_occ[i] = en_i - st_i
        aligned_tot_int[i]=intensity

    CommonMZ = np.unique(aligned_mzs)
    logged_Common_MZ = np.log10(CommonMZ)*slope1 + intercept1
    review_diff = np.insert(np.diff(logged_Common_MZ),0,1000000)
    bad_indices = np.where(review_diff<5)[0]
    for i in range(len(bad_indices)-1,-1,-1):
        # modify property table first
        j = bad_indices[i]
        CommonMZ[j-1]= (CommonMZ[j-1]*aligned_tot_int[j-1]+CommonMZ[j]*aligned_tot_int[j])/(aligned_tot_int[j-1]+aligned_tot_int[j])
        aligned_tot_int[j-1]+=aligned_tot_int[j]
        aligned_occ[j-1]+=aligned_occ[j]
        CommonMZ[j]=0
        aligned_tot_int[j]=0
        aligned_occ[j]=0

        st_i = int(all_groups[j-1])
        mid_i = int(all_groups[j])
        en_i = int(all_groups[j+1])
        aligned_mzs[st_i:en_i]=CommonMZ[j-1]
        aligned_mzs_group[mid_i:]-=1

    property_mask = (CommonMZ>0)
    
    full_table = np.column_stack((sorted_sum,aligned_mzs,aligned_mzs_group))
    aligned_table = full_table[np.lexsort((full_table[:,2],full_table[:,3],full_table[:,0]))]
    aligned_table = aligned_table[np.invert(np.insert(np.diff(np.take(aligned_table,4,axis = 1))==0,0,0))]
    property_table = np.column_stack((CommonMZ,aligned_tot_int,aligned_occ))
    property_table = property_table[property_mask]
    #np.savetxt("aligned_stat.csv",aligned,delimiter = ',')
    #np.savetxt("cluster_alignment.csv",addon2,delimiter = ',')

    return aligned_table,property_table

def gen_aligned_array(aligned_table,property_table,inten_thresh=0,occ_thresh=3,method='and'):
    n_scans = int(aligned_table[-1][0])+1
    intensities = np.take(property_table,1,axis = 1)
    if (inten_thresh==-1):
        inten_thresh = np.median(intensities)
    inten_mask = intensities > inten_thresh
    occurences = np.take(property_table,2,axis = 1)
    if (occ_thresh == -1):
        occ_thresh = np.median(occurences)
    occ_mask = occurences >= occ_thresh
    if method == 'and':
        group_mask = inten_mask * occ_mask
    elif method == 'or':
        group_mask = inten_mask + occ_mask

    real_peaks = property_table[group_mask]
    real_mz = np.take(real_peaks,0,axis = 1 )
    peak_groups = np.arange(len(property_table))[group_mask]
    aligned_mzs = np.take(aligned_table,3,axis = 1)
    peak_mask = np.isin(aligned_mzs,real_mz)
    filtered_signal = aligned_table[peak_mask]
    

    big_table = np.zeros((n_scans,len(real_mz)))
    idx = np.take(filtered_signal,0,axis = 1)
    idx = idx.astype(int)
    for i in range(n_scans):
        id_mask = (idx==i)
        pixel = filtered_signal[id_mask]
        pixel_mzs = np.take(pixel,3,axis = 1)
        mz_mask = np.isin(real_mz,pixel_mzs)
        big_table[i][mz_mask] = np.take(pixel,2,axis = 1)

    useful_table = np.row_stack((real_peaks.T,big_table))
    useful_table = np.column_stack((np.arange(-3,n_scans),useful_table))

    np.savetxt("alignment_15row_result_v2.csv",useful_table,delimiter = ',')
    return big_table


#main starts from this line

st = time.time()

p = ImzMLParser('4_row.imzML')

if p.coordinates[-1][1] == 1:
    ms_data_type = 'SC'
elif p.coordinates[-1][1] > 1:
    ms_data_type = 'IM'

p0,xmax,ymax= read_msi_data(p,thresh=0)
print(len(p0),"before bad scan removal")
timing_test("reading")

if ms_data_type == 'SC':
    for pixel in p0:
        if max(pixel[5])<15000:
            pixel[-1]=0
            #print("a scan is found to be bad")
            p0.remove(pixel)

print(len(p0),"after bad scan removal")
timing_test("bad scan removal")

c = np.zeros(6)
p1 = copy.copy(p0)

for i in range(len(p1)):
    p1[i][4],p1[i][5],c2 = combine_msi_data2(p1[i])
    p1[i][6] = np.column_stack((p1[i][4],p1[i][5]))

timing_test("combine")

aligned_table,peak_properties = ms_align2(p1)

timing_test("alignment")

aligned_matrix = gen_aligned_array(aligned_table,peak_properties,inten_thresh=-1,occ_thresh=-1)

timing_test("big_table")

print("done!")