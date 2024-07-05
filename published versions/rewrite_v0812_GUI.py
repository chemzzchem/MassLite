import numpy as np
import copy
import csv
import pandas as pd
from math import exp, sqrt, pi,log10
import time
import pymzml
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
from scipy.spatial.distance import pdist
from pyimzml.ImzMLParser import ImzMLParser
from scipy.signal import argrelextrema
from sklearn.cluster import KMeans
import umap

global_slope = 1000000/log10(np.e)
global_intercept = -(log10(50)*global_slope)
global_st = time.time()
# global file_imzML_path
file_imzML_path = ''
mass_threshold = 5
noise_threshold = 3000
mv_threshold = 3
max_int_threshold = 0
tic_threshold = 0
p0 = []
group_label = None
group_int = None
source_scan = []
source_file = []
cell_region = None
bg_region = None
marker_i = None
peak_info = None

aligned_table = None
trimmed_table = None

# create main window

main_window = tk.Tk()
main_window.title("MassLite GUI")

# function part

def timing_test(kw="this"):
    en = time.time()
    print("till ",kw," step takes",en - global_st)

def read_msi_data(p,thresh = 3000,ppm=5):
    a = max(p.coordinates)
    xmax,ymax,zmax = a
    peaks_all = []
    for idx,(x,y,z) in enumerate(p.coordinates):
        data = np.transpose(p.getspectrum(idx))
        data = data[data[:,1]>0]
        data2 = combine_peaks(mzs=data[:,0],ints=data[:,1],ppm=ppm)
        data2 = data2[np.where(data2[:,1]>thresh)[0]]
        n_peaks = len(data2)
        peaks_current = np.column_stack((np.ones(n_peaks)*idx,data2))
        peaks_all.append(peaks_current)
    return(peaks_all)

def read_ms_data(p,thresh = 3000,ppm=5):
    p1 = []
    org_n_scan = len(p0)
    for n, spec in enumerate(p):
        data = spec.centroidedPeaks
        n_peaks = len(data)
        peaks_current = np.column_stack((np.ones(n_peaks)*(n+org_n_scan),data))
        peaks_current = peaks_current[peaks_current[:,2]>thresh]
        p1.append(peaks_current)
    return(p1)

def combine_peaks(mzs,ints,mzr=[],ppm=5,collapse = False):
    if mzr == []: mzr = rela_transform(mzs)
    data = np.column_stack((mzs,ints,mzr))
    rel_dis = np.diff(mzr)
    data2 = np.insert(data,np.where(np.diff(mzr)>2*ppm)[0]+1,[0,0,1],axis = 0)
    data2 = np.insert(data2,0,[0,0,1],axis = 0)
    data2 = np.append(data2,[[0,0,1]],axis = 0)
    peak_pos = argrelextrema(data2[:,1],np.greater)[0]
    peaks = np.zeros((len(peak_pos),2))
    if collapse == False:
        for i in range(len(peaks)):
            peaks[i] = three_point_2(data2[peak_pos[i]-1],data2[peak_pos[i]],data2[peak_pos[i]+1])
    elif collapse == True:
        low_pos = argrelextrema(data2[:,1],np.less_equal)[0]
        for i in range(len(peaks)):
            peak = data2[low_pos[i]:low_pos[i+1]]
            peaks[i][0] = sum(peak[:,0]*peak[:,1])/sum(peak[:,1])
            peaks[i][1] = sum(peak[:,1])

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

def classify_scans(p0,ncl=3):
    fixed_bin = np.arange(int(p0[0][0][1]//100),int(p0[0][-1][1]//100)+2)*100
    fixed_bin[0]=0
    fixed_bin[-1]+=500
    bin_indices = np.zeros(len(fixed_bin),dtype=int)
    tics = np.zeros(len(p0))
    peak_int_hist = np.zeros((len(p0),len(fixed_bin)-1))
    for i in range(len(p0)):
        tics[i] = np.sum(p0[i][:,2])
        hist,bi = np.histogram(p0[i][:,1],bins = fixed_bin)
        bin_indices[1:]=np.cumsum(hist)
        scan_split = np.array_split(p0[i][:,2],bin_indices)
        peak_int_hist[i] = [np.sum (part) for part in scan_split][1:-1]
    peaks1 = np.where(peak_int_hist==0,1,peak_int_hist)
    if(checkbox_log_scaling_var.get()==1):
        peaks1 = np.log10(peaks1)
    if(checkbox_umap_var.get()==1):
        reducer = umap.UMAP()
        peaks1 = reducer.fit_transform(peaks1)
    kmeans = KMeans(n_clusters=ncl, random_state=0).fit(peaks1)
    int_log_label = kmeans.labels_
    tic_avg = np.zeros(ncl)
    for i in range(ncl):
        tic_avg[i] = np.sum(tics[np.where(int_log_label==i)[0]])
        tic_avg[i] = tic_avg[i]/len(np.where(int_log_label==i)[0])

    return(int_log_label,tic_avg)

def locate_peaks(p,marker_range = [(760.40,760.70)]):
    global marker_i, peak_info
    peak_edge_thresh = float(Entry_bg_level.get())/100
    marker_i = np.zeros(len(p))
    for i in range(len(p)):
        mask_all = np.zeros(len(p[i]),dtype=bool)
        for lower,upper in marker_range:
            mask = (p[i][:,1]>=lower) & (p[i][:,1]<=upper)
            mask_all = mask_all + mask
        if mask_all.any():
            marker_i[i]+=np.sum(p[i][mask_all,2])
    peak_max = np.max(marker_i)
    peak_thresh = 0.10
    sigma = 1.0
    peaks, _ = find_peaks(marker_i,height = peak_thresh * peak_max)
    n_peak_old = -1
    while len(peaks)!=n_peak_old and sigma <=10:
        n_peak_old = len(peaks)
        smoothed_signal = gaussian_filter(marker_i,sigma=sigma)
        peaks, _ = find_peaks(smoothed_signal,height = peak_thresh * peak_max)
        print("when sigma =",sigma,", ",len(peaks), " peaks are found")
        sigma += 1
    reversed_signal = (-1)*smoothed_signal
    lows, _ = find_peaks(reversed_signal)
    lows = np.insert(lows,0,0)
    lows = np.append(lows,len(marker_i))

    # store peak region info in 8 columns,peak index, peak intensity, l/r minimum, l/r 20% st, l/r 5% st

    peak_info = np.zeros((len(peaks),8))
    for i in range(len(peaks)):
        l_st = lows[lows<peaks[i]].max()
        r_st = lows[lows>peaks[i]].min()
        peak_info[i,2]=l_st
        peak_info[i,3]=r_st
        peak_info[i,0] = np.argmax(marker_i[l_st:r_st])+l_st
        peak_info[i,1]=marker_i[int(peak_info[i,0])]
        l_index = l_st
        r_index = r_st-1
        while marker_i[l_index]<peak_edge_thresh*peak_info[i,1] and l_index<peak_info[i,0]:
            l_index += 1
        while marker_i[r_index]<peak_edge_thresh*peak_info[i,1] and r_index>peak_info[i,0]:
            r_index -= 1
        peak_info[i,6]=l_index
        peak_info[i,7]=r_index
        while marker_i[l_index]<0.2*peak_info[i,1] and l_index<peak_info[i,0]:
            l_index += 1
        while marker_i[r_index]<0.2*peak_info[i,1] and r_index>peak_info[i,0]:
            r_index -= 1
        peak_info[i,4]=l_index
        peak_info[i,5]=r_index
    global cell_region,bg_region
    cell_region = find_peak_regions(peak_info)
    cell_group_update(org=marker_i,smoothed=smoothed_signal,p_info = peak_info)
    bg_region = find_bg_regions(peak_info)
    radio_in_cell.config(state=tk.NORMAL)


def sort_regions(rg,IsCell = True):    
    if IsCell == True:
        rg[:,1]+=1
        return rg.astype(int)
    elif IsCell == False:
        mask = (rg[:,1]-rg[:,0]==0)
        return rg[~mask].astype(int)

def find_peak_regions(peak_info):
    return (sort_regions(np.take(peak_info,(4,5),axis = 1),IsCell=True))

def find_bg_regions(peak_info):
    peak_region_broad_st = np.take(peak_info,6,axis=1)
    peak_region_broad_en = np.take(peak_info,7,axis=1)
    bg_region_st = np.insert(peak_region_broad_en,0,0)
    bg_region_en = np.append(peak_region_broad_st,len(marker_i))
    return(sort_regions(np.column_stack((bg_region_st,bg_region_en)),IsCell = False))

def norm_peaks(p0,norm_range = [(0,10000)]):
    p = p0
    mz_min = 0
    mz_max = 10000
    norm_i = np.zeros(len(p))
    for i in range(len(p)):
        mask_all = np.zeros(len(p[i]),dtype=bool)
        for lower,upper in norm_range:
            mask = (p[i][:,1]>=lower) & (p[i][:,1]<=upper)
            mask_all = mask_all + mask
        if mask_all.any():
            norm_i[i]+=np.sum(p[i][mask_all,2])
    norm_i = norm_i/1000000
    for i in range (len(p)):
        norm_fa = norm_i[i]
        p[i][:,2]=p[i][:,2]/norm_fa
    return p
    
def ms_alignment(data = p0,ppm=5):
    # initiating data array
    pixels = data
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
        clusters = fcluster(linkage_matrix,t=ppm,criterion = 'distance')
        cl,cl_location = np.unique(clusters,return_index=True)
        cl_location.sort()
        groups_in_chunk = np.array_split(chunk,cl_location[1:])
        for group in groups_in_chunk:
            if len(group)>=2:
                avg_mz,tot_int = weighted_avg(group[:,1],group[:,2])
                group[:,1]=avg_mz

    # re-evaluate the ppm difference and merge adjacent peaks from neiboring groups
    data[:,3]=rela_transform(data[:,1])
    data[1:,4]=np.diff(data[:,3])
    peak_start = np.where(data[:,4]>0)[0]
    peak_fake = np.where((data[:,4]>0) & (data[:,4]<ppm))[0]
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
        peak_fake = np.where((data[:,4]>0) & (data[:,4]<ppm))[0]
    
    # remove duplicated peaks by first sorting m/z, then scan number, then intensity
    lexsort_indices = np.lexsort((data[:, 2], data[:, 0], data[:, 1]))
    data = data[lexsort_indices]
    data_rev = data[::-1]
    id_mz = data_rev[:,:2]
    unique_mzs, indices_keep = np.unique(id_mz,axis = 0, return_index=True)
    data=data_rev [indices_keep]
    data = data[data[:,2]>np.median(data[:,2])]
    global aligned_table
    aligned_table = data[:,:3]
    log_message("Alignment completed.")
    label_al_status.config(text="Alignment completed.")
    checkbox_al.config(state=tk.NORMAL)
    return(data[:,:3])

def bg_removal(data,ratio=1):
    if checkbox_autobg_var.get()==1:
        bgscan_id = [num for start,end in bg_region for num in range(start,end)]
    elif checkbox_autobg_var.get()==0:
        if Entry_bg_manual.get() == '':
            log_message("BG not selected!")
            return data,[]
        else:
            bg_ranges = []
            bg_range_text = Entry_bg_manual.get()
            bg_range_strings = bg_range_text.split(',')
            for bgrange_str in bg_range_strings:
                bgrange_vals = bgrange_str.split('-')
                if len(bgrange_vals) ==2:
                    try:
                        start_value = int(bgrange_vals[0])
                        end_value = int(bgrange_vals[1])
                        bg_ranges.append((start_value,end_value))
                    except ValueError:
                        log_message("BG selection range is invalid.")
                        return data,[]
                else:
                    log_message("BG selection range is invalid.")
                    return data,[]
    
    lexsort_indices = np.lexsort((data[:, 0], data[:, 2], data[:, 1]))
    data = data[lexsort_indices]
    data = data[::-1]
    unique_mz_all,indices_all = np.unique(data[:,1],return_index=True)
    indices_all.sort()
    data_split = np.array_split(data,indices_all[1:])
    for peaks in data_split:
        bg_indices = np.where(np.isin(peaks[:,0],bgscan_id))
        if np.any(bg_indices):
            bg_peaks = peaks[bg_indices]
            if np.max(peaks[:,2])<=(np.max(bg_peaks[:,2])*ratio):
                peaks[:,0]=np.ones(len(peaks))*(-1)
    peaks_in_bg = np.unique(np.take(data[np.where(data[:,0]==-1)],1,axis=1))
    indices = np.where(data[:,0]>=0)
    data = data[indices]
    lexsort_indices = np.lexsort((data[:, 2], data[:, 1], data[:, 0]))
    data = data[lexsort_indices]
    
    print("done")
    return(data,peaks_in_bg)


def trim_data(data,missing_value = 2,int_tot_thresh = 0.0, int_max_thresh = 0.0):

    # missing value, total intensity, max intensity filters
    n_pixels = len(np.unique(data[:,0]))
    lexsort_indices = np.lexsort((data[:, 0], data[:, 2], data[:, 1]))
    data = data[lexsort_indices]

    data_rev = data[::-1]
    unique_mzs,indices_keep,mz_counts = np.unique(data_rev[:,1],return_counts=True,return_index=True)

    # missing value filter
    if (missing_value>0) and (missing_value<1):
        missing_value=missing_value*n_pixels
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
    return data2

def raw_scan2cell(s):
    c = []
    for i in range(len(cell_region)):
        cell = np.vstack(s[cell_region[i][0]:cell_region[i][1]])[:,1:3]
        sort_indices = np.argsort((cell[:,0]))
        cell = cell[sort_indices]
        single_cell = combine_peaks(mzs = cell[:,0],ints = cell[:,1],collapse=True)
        single_cell = np.column_stack((np.ones(len(single_cell))*i,single_cell))
        c.append(single_cell)
    return c

def aligned_scan2cell(s):
    c = np.zeros((0,3))
    for i in range(len(cell_region)):
        mask = np.where(np.isin(s[:,0], range (cell_region[i][0],cell_region[i][1]+1)))[0]
        cell = s[mask,1:3]
        sort_indices = np.argsort((cell[:,0]))
        cell = cell[sort_indices]
        single_cell = combine_peaks(mzs = cell[:,0],ints = cell[:,1],collapse=True)
        single_cell = np.column_stack((np.ones(len(single_cell))*i,single_cell))
        c=np.row_stack((c,single_cell))
    return c

def exp_org_table(p,exp_file_path,datatype = 'scan'):
    max_rows = max(arr.shape[0] for arr in p)
    df_with_nan = [pd.DataFrame(np.vstack((arr[:,1:3], np.full((max_rows - arr.shape[0], 2), np.nan)))) for arr in p]
    # Create a list of headers based on the pattern
    headers = [datatype] * len(p)*2
    max_width = len(str(len(p)))
    for i in range(len(p)):
        headers[2*i]+= "_"+ f"{int(i):0{max_width}}" + "_mz"
        headers[2*i]+= "_"+ f"{int(i):0{max_width}}" + "_intensity"
    # Create a DataFrame by vertically stacking the arrays
    df = pd.concat(df_with_nan,axis = 1)
    df.columns = headers
    output_path = exp_file_path+time.strftime("%Y%m%d-%H%M%S")+".csv"
    df.to_csv(output_path, index=False)
    log_message(output_path + " exported.")

def exp_aligned_table(data,exp_file_path,datatype = 'scan'):
    p = data[:,0:3]
    df = pd.DataFrame(p, columns=[datatype, "mz", "Intensity"])
    pivot_df = df.pivot_table(index=datatype, columns="mz", values="Intensity", aggfunc='sum')
    #pivot_df = df.pivot(index=datatype, columns="mz", values="Intensity")
    output_path = exp_file_path+time.strftime("%Y%m%d-%H%M%S")+".csv"
    pivot_df = pivot_df.T
    pivot_df.to_csv(output_path)
    log_message(output_path + " exported.")

# define actions

def on_browse_button_click():
    global file_imzML_path 
    file_path= filedialog.askopenfilename()
    if file_path:
        file_imzML_path = file_path
        Entry_imzML_file.delete(0, tk.END)
        Entry_imzML_file.insert(0, file_path)

def on_read_button_click():
    global file_imzML_path, p0, mass_threshold, noise_threshold,source_scan,source_file
    log_message("Reading "+ file_imzML_path)
    reply = "Read!"
    if (len(p0)!=0): reply = "New Data Added!"
    old_length = len(p0)
    if file_imzML_path:
        if(file_imzML_path[-5:]=='imzML'):
            ms_file = ImzMLParser(filename=file_imzML_path)
            noise_threshold = float(Entry_noise_thresh.get())
            mass_threshold = float(Entry_ppm_thresh.get())
            p0.extend(read_msi_data(ms_file,thresh=noise_threshold))
        elif (file_imzML_path[-4:]=='mzML'):
            ms_file = pymzml.run.Reader(file_imzML_path)
            # ms_file = ImzMLParser(filename=file_imzML_path)
            noise_threshold = float(Entry_noise_thresh.get())
            mass_threshold = float(Entry_ppm_thresh.get())
            p0.extend(read_ms_data(ms_file,thresh=noise_threshold))
        else:
            log_message("Invalid data type!")
    new_length = len(p0)
    source_scan.append((old_length,new_length))
    source_file.append((file_imzML_path))
    log_message("scan "+str(old_length)+" to "+str(new_length) + " read from "+ file_imzML_path)
    button_read.config(text="Add Data")
    if label_al_status.cget("text") == 'Alignment completed.':
        label_al_status.config(text="New Data Added. New Alignment Needed.")

def on_align_button_click():
    global aligned_table
    log_message(text="Alignment in progress")
    print("align button clicked")
    aligned_table = ms_alignment(p0,ppm=float(Entry_ppm_thresh.get()))
    log_message(text="Alignment done!")
    checkbox_al.config(state=tk.NORMAL)

def on_exp_button_click():
    out_pref = radio_exp_var.get()
    if checkbox_al_var.get() == 0:
        out_pref += '_unaligned_'
        if radio_exp_var.get() == 'scan':
            exp_org_table(p=p0,exp_file_path=out_pref,datatype='scan')
        elif radio_exp_var.get() == 'cell':
            exp_org_table(p=raw_scan2cell(p0),exp_file_path=out_pref,datatype='cell')
    elif checkbox_al_var.get() == 1:
        out_pref += '_aligned_'
        aligned = aligned_table
        if checkbox_bg_var.get() == 1:
            aligned,bg_mz = bg_removal(aligned)
            if len(bg_mz) != 0:
                out_pref += '_bg_removed_'
            elif len(bg_mz) == 0:
                out_pref += '_bg_NOT_removed_'
        if radio_exp_var.get() == 'cell':
            aligned = aligned_scan2cell(aligned)
        if checkbox_fi_var.get() == 1:
            aligned = trim_data(aligned,missing_value=float(Entry_mv_thresh.get()),int_tot_thresh=float(Entry_ticm_thresh.get()),int_max_thresh=float(Entry_mi_thresh.get()))
            out_pref += 'filtered'
        exp_aligned_table(aligned,exp_file_path=out_pref,datatype = radio_exp_var.get())

def on_classify_button_click():
    global group_label,group_int
    ncl=int(Entry_gn.get())
    group_label,group_int = classify_scans(p0,ncl=ncl)
    # combos = [("Group "+ str(i) + " TIC on average: " + str(group_int[i])) for i in range(3)]
    combos = [(str(i) + ": " + str(len(np.where(group_label==i)[0]))+ " @ "+"{:.3e}".format(group_int[i])) for i in range(ncl)]
    combobox_class['values'] = combos
    ax1.clear()
    tics = [np.sum(scan[:,2]) for scan in p0]
    tics = tics/max(tics)
    ax1.bar(x=range(len(tics)),height=tics,width = 0.9,label='all scans',color='g')
    ax1.legend()
    ax1.set_xlabel("scan number")
    ax1.set_ylabel("normalized TIC")
    canvas_class.draw()

def on_update_button_click():
    ax1.clear()
    group = combobox_class.current()
    tics = [np.sum(scan[:,2]) for scan in p0]
    tics = tics/max(tics)
    selected = np.where(group_label==group)[0]
    others = np.where(group_label!=group)[0]
    selected = selected[selected <=600]
    others = others[others <=600]
    ax1.bar(x=selected,height=tics[selected],width = 0.9,label='selected scans',color='r')
    ax1.bar(x=others,height=tics[others],width = 0.9,label='other scans',color='g')
    ax1.legend()
    ax1.set_xlabel("scan number")
    ax1.set_ylabel("normalized TIC")
    ax1.set_title("n=3")
    canvas_class.draw()

def cell_group_update(org,smoothed,p_info):
    
    ax2.clear()
    x = np.arange(len(org))
    mask = np.take(p_info,0,axis = 1)
    mask = mask.astype(int)
    ax2.plot(x,org,label='Raw Data')
    ax2.plot(x,smoothed,label='Smoothed Data')
    ax2.scatter(x[mask],org[mask],color = 'r',marker = 'o')
    for start, end in cell_region:
        ax2.fill_between(x, org, where=(x >= start) & (x <= end), alpha=0.5)
    ax2.legend()
    ax2.set_xlabel("scan number")
    ax2.set_ylabel("normalized TIC")
    ax2.set_title(str(len(cell_region)) + " cells found")
    canvas_group.draw()

def on_drop_button_click():
    group = combobox_class.current()
    global p0
    mask = np.where(group_label==group)[0]
    dropped_p0 = [arr for i,arr in enumerate (p0) if i not in mask]
    p0 = dropped_p0
    log_message(message="Group " + str (group) + " dropped.")

def on_grouping_button_click():
    range_text = Entry_marker_range.get()
    if (range_text == ''):
        locate_peaks(p=p0)
    else:
        ranges = []
        range_strings = range_text.split(',')
        for range_str in range_strings:
            range_values = range_str.split('-')
            if len(range_values) == 2:
                try:
                    start_value = float(range_values[0])
                    end_value = float(range_values[1])
                    ranges.append((start_value,end_value))
                except ValueError:
                    log_message("Marker range is invalid.")
            else:
                log_message("Marker range is invalid.")
        locate_peaks(p=p0,marker_range=ranges)

    log_message(str(len(peak_info))+" cells were found.")

def on_exp_al_toggle():
    if checkbox_al_var.get() == 1:
        checkbox_bg.config(state=tk.NORMAL)
        checkbox_fi.config(state=tk.NORMAL)
    else:
        checkbox_bg.deselect()
        checkbox_bg.config(state=tk.DISABLED)
        checkbox_fi.deselect()
        checkbox_fi.config(state=tk.DISABLED)

def on_autobg_toggle():
    if checkbox_autobg_var.get() == 1:
        Entry_bg_manual.config(state=tk.DISABLED)
    else:
        Entry_bg_manual.config(state=tk.NORMAL)

def on_combobox_select(event):
    data_type = int(combobox_class.get()[0])

def on_combobox_class_select(event):
    data_type = int(combobox_class.get()[0])

def on_combobox_group_select(event):
    selected_sigma = int(combobox_group.get()[6])

def on_radio_exp_select():
    selection = radio_exp_var.get()

def on_debug_button_click():
    log_message("Debug mode started.")
    log_message("Debug mode ended")

def log_message(message):
    listbox_log.insert(tk.END, message)
    listbox_log.see(tk.END)

# define GUI grid

label_imzML_path = tk.Label(main_window, text="Select imzML file")
label_imzML_path.grid(row=0,column=0)
button_browse = tk.Button(main_window, text= "Browse file...", command= on_browse_button_click)
button_browse.grid(row=0,column=1)
Entry_imzML_file = tk.Entry(main_window,width = 60)
Entry_imzML_file.grid(row=0,column=2,columnspan=4)
button_read = tk.Button(main_window, text= "Read Data", command= on_read_button_click)
button_read.grid(row=0,column=6)
separator_l_h_1 = ttk.Separator(main_window,orient='horizontal')
separator_l_h_1.grid(row=2,column=0,columnspan=7,sticky="ew", pady=10)

label_gn_status = tk.Label(main_window,text="# Groups")
label_gn_status.grid(row=3,column=0)
Entry_gn = tk.Entry(main_window)
Entry_gn.insert(0,'3')
Entry_gn.grid(row=3,column=1)
checkbox_log_scaling_var = tk.IntVar()
checkbox_log_scaling = tk.Checkbutton(main_window,text="log scaling",variable = checkbox_log_scaling_var,command=None)
checkbox_log_scaling.grid(row=3,column=2)
checkbox_umap_var = tk.IntVar()
checkbox_umap = tk.Checkbutton(main_window,text="UMAP",variable = checkbox_umap_var,command=None)
checkbox_umap.grid(row=3,column=3)
button_MLC = tk.Button(main_window, text= "Find void scans", command= on_classify_button_click)
button_MLC.grid(row=3,column=4)
combobox_class = ttk.Combobox(main_window,values=[])
combobox_class.grid(row=3,column=5)
combobox_class.bind("<<ComboboxSelected>>",on_combobox_class_select)

button_drop = tk.Button(main_window, text= "Drop", command=on_drop_button_click)
button_drop.grid(row=5,column=6)
button_update = tk.Button(main_window, text= "Update", command=on_update_button_click)
button_update.grid(row=7,column=6)

fig_class = Figure(figsize=(7, 2.5), dpi=100)
ax1 = fig_class.add_subplot(1, 1, 1)
canvas_class = FigureCanvasTkAgg(fig_class, master=main_window)
canvas_class.get_tk_widget().grid(row=4,column=0,columnspan=6,rowspan=9)

separator_l_h_2 = ttk.Separator(main_window,orient='horizontal')
separator_l_h_2.grid(row=13,column=0,columnspan=7,sticky="ew", pady=10)


label_marker_range = tk.Label(main_window,text="Cell marker ranges (e.g. 782.56-782.60,760.54-760.58)")
label_marker_range.grid(row=14,column=0,columnspan=3)
Entry_marker_range = tk.Entry(main_window,width = 50)
Entry_marker_range.grid(row=14,column=3,columnspan=4)

label_bg_level = tk.Label(main_window,text="BG level(%)")
label_bg_level.grid(row=15,column=0)
Entry_bg_level = tk.Entry(main_window)
Entry_bg_level.insert(0,'5')
Entry_bg_level.grid(row=15,column=1)
button_grouping = tk.Button(main_window, text= "Grouping", command=on_grouping_button_click)
button_grouping.grid(row=15,column=2)
combobox_group = ttk.Combobox(main_window,values=[],state=tk.DISABLED)
combobox_group.grid(row=15,column=3)
combobox_group.bind("<<ComboboxSelected>>",on_combobox_group_select)
button_cells = tk.Button(main_window, text= "Confirm Cells", command=None,state=tk.DISABLED)
button_cells.grid(row=15,column=4)

checkbox_autobg_var = tk.IntVar()
checkbox_autobg = tk.Checkbutton(main_window,text="Auto BG",variable = checkbox_autobg_var,command=on_autobg_toggle)
checkbox_autobg.grid(row=15,column=5)

fig_group = Figure(figsize=(7, 2.5), dpi=100)
ax2 = fig_group.add_subplot(1, 1, 1)
canvas_group = FigureCanvasTkAgg(fig_group, master=main_window)
canvas_group.get_tk_widget().grid(row=16,column=0,columnspan=6,rowspan=9)

separator_v_1=ttk.Separator(main_window,orient='vertical')
separator_v_1.grid(row=0,column=7,rowspan=24,sticky="ns", padx=10)

label_al_title = tk.Label(main_window,text="Alignment")
label_al_title.grid(row=0,column=8)
label_al_norm = tk.Label(main_window,text="Normalization range")
label_al_norm.grid(row=1,column=8) 
Entry_al_norm_range = tk.Entry(main_window,borderwidth=5)
Entry_al_norm_range.grid(row=1,column=9)
label_al_norm_range_demo = tk.Label(main_window,text="blank:without normalization, 1: tic normalization")
label_al_norm_range_demo.grid(row=1,column=10) 

label_mass = tk.Label(main_window,text="mass tolerance/ppm")
label_mass.grid(row=2,column=8)
Entry_ppm_thresh = tk.Entry(main_window, borderwidth=5)
Entry_ppm_thresh.insert(0,'5')
Entry_ppm_thresh.grid(row = 2,column=9)
label_intensity = tk.Label(main_window,text="Minimum intensity")
label_intensity.grid(row=3,column=8)
Entry_noise_thresh = tk.Entry(main_window, borderwidth=5)
Entry_noise_thresh.insert(0,'3000')
Entry_noise_thresh.grid(row = 3,column=9)
button_align = tk.Button(main_window, text= "Alignment", command=ms_alignment)
button_align.grid(row=2,column=10)
label_al_status = tk.Label(main_window,text="Alignment not performed")
label_al_status.grid(row=3,column=10)

separator_r_h_1 = ttk.Separator(main_window,orient='horizontal')
separator_r_h_1.grid(row=4,column=8,columnspan=4,sticky="ew", pady=10)

label_bg_manual = tk.Label(main_window,text="Manual BG scan selection")
label_bg_manual.grid(row=5,column=8)
Entry_bg_manual = tk.Entry(main_window, borderwidth=5,state=tk.DISABLED)
Entry_bg_manual.grid(row = 5,column=9)
label_mv = tk.Label(main_window,text="Missing Value")
label_mv.grid(row=6,column=8)
Entry_mv_thresh = tk.Entry(main_window, borderwidth=5)
Entry_mv_thresh.insert(0,'3')
Entry_mv_thresh.grid(row = 6,column=9)
label_mi = tk.Label(main_window,text="Max intensity Threshold")
label_mi.grid(row=7,column=8)
Entry_mi_thresh = tk.Entry(main_window, borderwidth=5)
Entry_mi_thresh.insert(0,'0')
Entry_mi_thresh.grid(row = 7,column=9)
label_ticm = tk.Label(main_window,text="Peak TIC Threshold")
label_ticm.grid(row=8,column=8)
Entry_ticm_thresh = tk.Entry(main_window, borderwidth=5)
Entry_ticm_thresh.insert(0,'0')
Entry_ticm_thresh.grid(row = 8,column=9)

separator_r_h_2 = ttk.Separator(main_window,orient='horizontal')
separator_r_h_2.grid(row=9,column=8,columnspan=4,sticky="ew", pady=10)

label_exp_scan = tk.Label(main_window,text="Export options")
label_exp_scan.grid(row=10,column=8)

radio_exp_var = tk.StringVar(value='scan')
radio_in_scan = tk.Radiobutton(main_window,text='in scans',variable=radio_exp_var,value='scan',command=on_radio_exp_select)
radio_in_scan.grid(row=11,column=8)
radio_in_cell = tk.Radiobutton(main_window,text='in cells',variable=radio_exp_var,value='cell',command=on_radio_exp_select,state=tk.DISABLED)
radio_in_cell.grid(row=11,column=9)
button_exp = tk.Button(main_window,text="Export",command=on_exp_button_click)
button_exp.grid(row = 11,column=10)

checkbox_al_var = tk.IntVar()
checkbox_al = tk.Checkbutton(main_window,text='with alignment',variable=checkbox_al_var,state=tk.DISABLED,command=on_exp_al_toggle)
checkbox_al.grid(row = 12, column=8,sticky="w")
checkbox_fi_var = tk.IntVar()
checkbox_fi = tk.Checkbutton(main_window,text='with peaks filtered',variable=checkbox_fi_var, state=tk.DISABLED)
checkbox_fi.grid(row = 13, column=8,sticky="w")
checkbox_bg_var = tk.IntVar()
checkbox_bg = tk.Checkbutton(main_window,text='with aligned background peaks removed',variable=checkbox_bg_var, state=tk.DISABLED)
checkbox_bg.grid(row = 14, column=8,columnspan=2,sticky="w")

separator_r_h_3 = ttk.Separator(main_window,orient='horizontal')
separator_r_h_3.grid(row=15,column=8,columnspan=4,sticky="ew", pady=10)

label_log = tk.Label(main_window,text='User Log')
label_log.grid(row=16,column=8)
button_debug = tk.Button(main_window,text='Debug',command=on_debug_button_click)
button_debug.grid(row=16, column=10)
listbox_log = tk.Listbox(main_window, width=80, height=12)
listbox_log.grid(row=17, column=8, columnspan=3,rowspan=5,padx=10, pady=10)
scrollbar_log = tk.Scrollbar(main_window, command=listbox_log.yview)
scrollbar_log.grid(row=17, column=11, pady=10, rowspan=5,sticky="ns")
listbox_log.configure(yscrollcommand=scrollbar_log.set)

# main start from here

main_window.mainloop()

main_window.mainloop()

print("end")