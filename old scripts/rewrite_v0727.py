import numpy as np
import copy
from math import exp, sqrt, pi,log10
import time
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
from scipy.spatial.distance import pdist
from pyimzml.ImzMLParser import ImzMLParser
from scipy.signal import argrelextrema
from sklearn.cluster import KMeans

global_slope = 1/(log10(1.000001))
global_intercept = -(log10(50)/log10(1.000001))
global_st = time.time()

def timing_test(kw="this"):
    en = time.time()
    print("till ",kw," step takes",en - global_st)

def read_msi_data(p,thresh = 3000):
    a = max(p.coordinates)
    xmax,ymax,zmax = a
    peaks_all = []
    for idx,(x,y,z) in enumerate(p.coordinates):
        data = np.transpose(p.getspectrum(idx))
        data = data[data[:,1]>1]
        data2 = combine_peaks(mzs=data[:,0],ints=data[:,1])
        n_peaks = len(data2)
        peaks_current = np.column_stack((np.ones(n_peaks)*idx,data2))
        peaks_all.append(peaks_current)
    return(peaks_all)

def combine_peaks(mzs,ints,mzr=[],ppm=5):
    if mzr == []: mzr = rela_transform(mzs)
    data = np.column_stack((mzs,ints,mzr))
    rel_dis = np.diff(mzr)
    data2 = np.insert(data,np.where(np.diff(mzr)>2*ppm)[0]+1,[0,0,1],axis = 0)
    data2 = np.append(data2,[[0,0,1]],axis = 0)
    peak_pos = argrelextrema(data2[:,1],np.greater)[0]
    peaks = np.zeros((len(peak_pos),2))
    for i in range(len(peaks)):
        peaks[i] = three_point_2(data2[peak_pos[i]-1],data2[peak_pos[i]],data2[peak_pos[i]+1])

    return peaks

def rela_transform(mzs):
    if isinstance(mzs, (float,int)):
        # If the input is an integer, perform some operations for single integer input
        result = log10(mzs) * global_slope + global_intercept  # Example operation, you can replace this with your desired logic
    elif isinstance(mzs, (list, np.ndarray)):
        # If the input is a list or a NumPy array, perform some operations for array input
        arr = np.array(mzs)  # Convert the input to a NumPy array for uniform processing
        result = np.log10(arr)* global_slope + global_intercept  # Example operation, you can replace this with your desired logic
        
    return result

def rela_transform_rev(mzr):
    if isinstance(mzr, (float,int)):
        result = pow(10,((mzr-global_intercept)/global_slope))
    elif isinstance(mzr, (list, np.ndarray)):
        # If the input is a list or a NumPy array, perform some operations for array input
        arr = np.array(mzr)  # Convert the input to a NumPy array for uniform processing
        result = np.power(10,((mzr-global_intercept)/global_slope))
    return result

def weighted_avg (mzs,ints):
    produc = mzs*ints
    return(sum(produc)/sum(ints),sum(ints))

def three_point_2 (a,b,c):
    x = np.array((a[0],b[0],c[0]))
    y = np.array((a[1],b[1],c[1]))
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

def classify_scans(p0):
    fixed_bin = np.arange(int(p0[0][0][1]//100),int(p0[0][-1][1]//100)+3)*100
    peak_hist = np.zeros((len(p0),len(fixed_bin)))
    for i in range
    print("done")
    p1 = p0
    return(p1)

def ms_alignment(pixels,ppm=5,missing_value = 2,int_tot_thresh = 0.0, int_max_thresh = 0.0):
    # initiating data array
    total_rows = sum(array.shape[0] for array in pixels)
    data = np.zeros((total_rows,5))
    data_org=np.vstack(pixels)
    lexsort_indices = np.lexsort((data_org[:, 0], data_org[:, 1]))
    data[:,0:3]=data_org[lexsort_indices]
    data[:,3]=rela_transform(data[:,1])

    #basic alignment with dynamic grouping in chunks
    chunk_size = 2* len(pixels)
    chunks = np.array_split(data,data.shape[0] // chunk_size)

    for chunk in chunks:
        linkage_matrix = linkage(chunk[:,3][:,np.newaxis],method = 'average')
        '''plt.figure(figsize=(10, 6))
        dendrogram(linkage_matrix)
        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('Sample Index')
        plt.ylabel('Distance')
        plt.show()'''
        clusters = fcluster(linkage_matrix,t=ppm,criterion = 'distance')
        cl,cl_location = np.unique(clusters,return_index=True)
        cl_location.sort()
        groups_in_chunk = np.array_split(chunk,cl_location[1:])
        for group in groups_in_chunk:
            if len(group)>=2:
                avg_mz,tot_int = weighted_avg(group[:,1],group[:,2])
                group[:,1]=avg_mz

        # print("done")

    # re-evaluate the ppm difference and merge adjacent peaks from neiboring groups
    data[:,3]=rela_transform(data[:,1])
    data[1:,4]=np.diff(data[:,3])
    peak_start = np.where(data[:,4]>0)[0]
    peak_fake = np.where((data[:,4]>0) & (data[:,4]<5))[0]
    n_redos = 0

    while len(peak_fake!=0):
        n_redos += 1
        print("reiteration of combination for ", n_redos, " times")

        mask = np.isin(peak_start,peak_fake)
        fake_pos = np.where(mask)[0]

        for faked in fake_pos:
            peak_index_prev = peak_start[faked-1]
            peak_index_next = peak_start[faked+1]
            avg_mz,tot_int = weighted_avg(data[peak_index_prev:peak_index_next,1],data[peak_index_prev:peak_index_next,2])
            data[peak_index_prev:peak_index_next,1]=avg_mz
        data[:,3]=rela_transform(data[:,1])
        data[1:,4]=np.diff(data[:,3])
        peak_start = np.where(data[:,4]>0)[0]
        peak_fake = np.where((data[:,4]>0) & (data[:,4]<5))[0]
    
    # remove duplicated peaks
    lexsort_indices = np.lexsort((data[:, 2], data[:, 0], data[:, 1]))
    data = data[lexsort_indices]
    data_rev = data[::-1]
    id_mz = data_rev[:,:2]
    unique_mzs, indices_keep = np.unique(id_mz,axis = 0, return_index=True)
    data=data_rev [indices_keep]

    # missing value, total intensity, max intensity filters
    lexsort_indices = np.lexsort((data[:, 0], data[:, 2], data[:, 1]))
    data = data[lexsort_indices]

    # debug with first 100 rows
    data = data[:100,:3]

    data_rev = data[::-1]
    unique_mzs,indices_keep,mz_counts = np.unique(data_rev[:,1],return_counts=True,return_index=True)

    # missing value filter
    if (missing_value>0) and (missing_value<1):
        missing_value=missing_value*len(pixels)
    mask_mv = ( mz_counts >= missing_value)

    # max intensity filter
    max_intensities = data_rev[indices_keep,2]
    mask_maxint = ( max_intensities >= int_max_thresh)

    # sum_intensity filter
    sum_intensities = np.zeros(len(unique_mzs))
    sum_intensities[0]=sum(data_rev[indices_keep[0]:,2])
    for i in range(1,len(unique_mzs)):
        sum_intensities[i]=sum(data_rev[indices_keep[i]:indices_keep[i-1],2])
    mask_totint = (sum_intensities >= int_tot_thresh)

    mask_all = mask_mv * mask_maxint * mask_totint
    mask_by_peak = np.repeat(mask_all,mz_counts)

    data2 = data[mask_by_peak]

    lexsort_indices = np.lexsort((data[:, 1], data[:, 0]))
    data = data[lexsort_indices]
    lexsort_indices = np.lexsort((data2[:, 1], data2[:, 0]))
    data2 = data2[lexsort_indices]

    timing_test("alignment")
    print("done!")
    return data


# main start from here

file_path = '4_row.imzML'
ms_file = ImzMLParser(filename=file_path)

p0 = read_msi_data(ms_file)

p1 = classify_scans(p0)

timing_test("data read-in")
p2 = ms_alignment(p1,missing_value=0.1)

print("end")