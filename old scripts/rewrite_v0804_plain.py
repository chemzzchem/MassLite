import numpy as np
import copy
import csv
import pandas as pd
from math import exp, sqrt, pi,log10
import time
import pymzml
import matplotlib.pyplot as plt
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

global_slope = 1/(log10(1.000001))
global_intercept = -(log10(50)/log10(1.000001))
global_st = time.time()
# global file_imzML_path
file_imzML_path = ''
mass_threshold = 5
noise_threshold = 3000
mv_threshold = 3
max_int_threshold = 0
tic_threshold = 0
p0 = []
dropped_p0 = None
group_label = None
group_int = None
source = []

aligned_table = None
trimmed_table = None

'''fig = Figure(figsize=(5, 4), dpi=100)
ax = fig.add_subplot(1, 1, 1)

# create main window

main_window = tk.Tk()
main_window.title("MassLite GUI")'''

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

'''def read_ms_data(p,thresh = 3000,ppm=5):
    peaks_all = []
    for n, spec in enumerate(p):
        data = spec.centroidedPeaks
        n_peaks = len(data)
        peaks_current = np.column_stack((np.ones(n_peaks)*n,data))
        peaks_current = peaks_current[peaks_current[:,2]>thresh]
        peaks_all.append(peaks_current)
    return(peaks_all)'''

def read_ms_data(p,thresh = 3000,ppm=5):
    global p0
    p1 = []
    org_n_scan = len(p0)
    for n, spec in enumerate(p):
        data = spec.centroidedPeaks
        n_peaks = len(data)
        peaks_current = np.column_stack((np.ones(n_peaks)*(n+org_n_scan),data))
        peaks_current = peaks_current[peaks_current[:,2]>thresh]
        p1.append(peaks_current)
    # p0.extend(p1)
    return(p1)

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

def merge_overlapping_ranges(ranges,peak_gap = 5,peak_width = 20):
    # Sort the ranges by their start values
    ranges.sort(key=lambda x: x[0])

    merged_ranges = []
    current_start, current_end = ranges[0]

    for start, end in ranges[1:]:
        if start <= current_end+peak_gap:  # Check for overlap
            current_end = max(current_end, end)  # Merge overlapping ranges
        else:
            if current_end-current_start>=peak_width:
                merged_ranges.append((current_start, current_end))
            current_start, current_end = start, end

    merged_ranges.append(np.array((current_start, current_end)))
    return merged_ranges

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
    peak_int_hist = np.where(peak_int_hist==0,1,peak_int_hist)
    peak_log_int_hist = np.log10(peak_int_hist)
    kmeans = KMeans(n_clusters=ncl, random_state=0, n_init="auto").fit(peak_log_int_hist)
    int_log_label = kmeans.labels_
    tic_avg = np.zeros(ncl)
    for i in range(ncl):
        tic_avg[i] = np.sum(tics[np.where(int_log_label==i)[0]])
        tic_avg[i] = tic_avg[i]/len(np.where(int_log_label==i)[0])

    return(int_log_label,tic_avg)

def locate_peaks(p,marker_range = [(760.40,760.70)]):
#def locate_peaks(p,marker_range = [(782.54,782.58)]):
#def locate_peaks(p,marker_range = [(782.54,782.58),(760.56,760.60),(734.55,734.59),(788.52,788.56)]):
    marker_i = np.zeros(len(p))
    for i in range(len(p)):
        mask_all = np.zeros(len(p[i]),dtype=bool)
        for lower,upper in marker_range:
            mask = (p[i][:,1]>=lower) & (p[i][:,1]<=upper)
            mask_all = mask_all + mask
        if mask_all.any():
            marker_i[i]+=np.sum(p[i][mask_all,2])
    peak_max = np.max(marker_i)
    sigma = 1.0
    peaks, _ = find_peaks(marker_i,height = 0.01 * peak_max)
    n_peak_old = -1
    print("Without smoothing: ",len(peaks)," peaks")
    while len(peaks)!=n_peak_old:
        n_peak_old = len(peaks)
        smoothed_signal = gaussian_filter(marker_i,sigma=sigma)
        peaks, _ = find_peaks(smoothed_signal,height = 0.01 * peak_max)
        print("when sigma =",sigma,", ",len(peaks), " peaks are found")
        sigma += 1
    reversed_signal = (-1)*smoothed_signal
    lows, _ = find_peaks(reversed_signal)
    lows = np.insert(lows,0,0)
    lows = np.append(lows,len(marker_i))
    # store peak region info in 4 columns, left start, right start, intensity, index
    # peak_info = np.zeros((len(peaks),4))
    peak_info = np.zeros((len(peaks),6))
    for i in range(len(peaks)):
        l_st = lows[lows<peaks[i]].max()
        r_st = lows[lows>peaks[i]].min()
        peak_info[i,4]=l_st
        peak_info[i,5]=r_st
        peak_info[i,0] = np.argmax(marker_i[l_st:r_st])+l_st
        peak_info[i,1]=marker_i[int(peak_info[i,0])]
        l_index = l_st
        r_index = r_st-1
        while marker_i[l_index]<0.2*peak_info[i,1] and l_index<peak_info[i,0]:
            l_index += 1
        while marker_i[r_index]<0.2*peak_info[i,1] and r_index>peak_info[i,0]:
            r_index -= 1
        peak_info[i,2]=l_index
        peak_info[i,3]=r_index
    
    print("done!")
    # debug graphs
    peak_regions = np.take(peak_info,(2,3),axis = 1)
    x = range(len(marker_i))
    plt.plot(x, marker_i, label='Raw Data')
    plt.plot(x, smoothed_signal, label='Smoothed Data')
    for start, end in peak_regions:
        plt.fill_between(x, marker_i, where=(x >= start) & (x <= end), alpha=0.5)
    plt.legend()
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt_path = file_path + 'png'
    plt.title(plt_path)
    plt.show()

    return peak_info

def ms_alignment(pixels,ppm=5):
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

    return(data[:,:3])

def bg_removal(data,bgscan_id,ratio=1):
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

    # debug with first 100 rows
    # data = data[:100,:3]

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

def exp_org_table(p,exp_file_path):
    p0T = [np.hstack((pixel[0,0]*np.ones((2,1)),pixel[:,1:3].T)) for pixel in p]
    with open(exp_file_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(p0T)

def exp_aligned_table(data,exp_file_path,tp = False):
    p = data[:,0:3]
    df = pd.DataFrame(p, columns=["Scan", "mz", "Intensity"])
    pivot_df = df.pivot(index="Scan", columns="mz", values="Intensity")
    if tp == True:
        pivot_df=pivot_df.T
    pivot_df.to_csv(exp_file_path)
'''
# define actions

def on_browse_button_click():
    global file_imzML_path 
    file_path= filedialog.askopenfilename()
    if file_path:
        file_imzML_path = file_path
        Entry_imzML_file.delete(0, tk.END)
        Entry_imzML_file.insert(0, file_path)

def on_read_button_click():
    global file_imzML_path, p0, mass_threshold, noise_threshold
    label_status.config(text="Reading...")
    reply = "Read!"
    if (len(p0)!=0): reply = "New Data Added!"
    if file_imzML_path:
        ms_file = pymzml.run.Reader(file_imzML_path)
        # ms_file = ImzMLParser(filename=file_imzML_path)
        noise_threshold = float(Entry_noise_thresh.get())
        mass_threshold = float(Entry_ppm_thresh.get())
        p0.extend(read_ms_data(ms_file,thresh=noise_threshold))
    
    #p0 = read_msi_data(ms_file, thresh = Entry_noise_thresh.get())
    timing_test("data read-in")
    label_status.config(text=reply)
    button_read.config(text="Add Data")

def on_align_button_click():
    global aligned_table
    label_status.config(text="Alignment in progress")
    print("align button clicked")
    if dropped_p0:
        aligned_table = ms_alignment(dropped_p0,ppm=float(Entry_ppm_thresh.get()))
    else:
        aligned_table = ms_alignment(p0,ppm=float(Entry_ppm_thresh.get()))
    label_status.config(text="Alignment done!")

def on_exp_raw_button_click():
    if dropped_p0:
        exp_org_table(p=dropped_p0,exp_file_path="data_dropped.csv")
        label_status.config(text="cropped data file exported successfully.")
    else:
        exp_org_table(p=p0,exp_file_path="data_raw.csv")
        label_status.config(text="raw data file exported successfully.")
        
def on_exp_aligned_button_click():
    if aligned_table.any():
        exp_aligned_table(aligned_table,"data_aligned.csv")
        label_status.config(text="Aligned data exported successfully.")
    else:
        label_status.config(text="Alignment not done")
    
        
def on_exp_trimmed_button_click():
    if aligned_table.any() :
        global mv_threshold, max_int_threshold, tic_threshold
        global trimmed_table
        mv_threshold = float(Entry_mv_thresh.get())
        max_int_threshold = float(Entry_mi_thresh.get())
        tic_threshold = float(Entry_ticm_thresh.get())
        trimmed_table = trim_data(aligned_table,missing_value=mv_threshold,int_max_thresh=max_int_threshold, int_tot_thresh=tic_threshold)
        exp_aligned_table(trimmed_table,"data_aligned_trimmed.csv")
        label_status.config(text="Trimmed aligned data exported successfully")
    else:
        label_status.config(text="Alignment not done")

def on_classify_button_click():
    global group_label,group_int
    group_label,group_int = classify_scans(p0,ncl=3)
    # combos = [("Group "+ str(i) + " TIC on average: " + str(group_int[i])) for i in range(3)]
    combos = [(str(i) + ": " + str(len(np.where(group_label==i)[0]))+ " @ "+str(group_int[i])) for i in range(3)]
    combobox_class['values'] = combos
    ax.clear()
    tics = [np.sum(scan[:,2]) for scan in p0]
    tics = tics/max(tics)
    ax.bar(x=range(len(tics)),height=tics,width = 0.9,label='all scans',color='g')
    ax.legend()
    ax.set_xlabel("scan number")
    ax.set_ylabel("normalized TIC")
    canvas.draw()

def on_update_button_click():
    ax.clear()
    group = combobox_class.current()
    tics = [np.sum(scan[:,2]) for scan in p0]
    tics = tics/max(tics)
    selected = np.where(group_label==group)[0]
    others = np.where(group_label!=group)[0]
    selected = selected[selected <=600]
    others = others[others <=600]
    ax.bar(x=selected,height=tics[selected],width = 0.9,label='selected scans',color='r')
    ax.bar(x=others,height=tics[others],width = 0.9,label='other scans',color='g')
    ax.legend()
    ax.set_xlabel("scan number")
    ax.set_ylabel("normalized TIC")
    canvas.draw()

def on_drop_button_click():
    group = combobox_class.current()
    global dropped_p0, p0
    mask = np.where(group_label==group)[0]
    dropped_p0 = [arr for i,arr in enumerate (p0) if i not in mask]
    p0 = dropped_p0
    label_status.config(text="Group " + str (group) + " dropped.")

def on_combobox_select(event):
    data_type = int(combobox_class.get()[0])

# define GUI grid
# define labels
label_imzML_path = tk.Label(main_window, text="Select imzML file")
label_imzML_path.grid(row=0,column=0)
label_mass = tk.Label(main_window,text="mass tolerance")
label_mass.grid(row=1,column=0)
label_ppm = tk.Label(main_window,text="ppm")
label_ppm.grid(row=1,column=2)
label_intensity = tk.Label(main_window,text="Intensity")
label_intensity.grid(row=1,column=3)
label_status = tk.Label(main_window,text="")
label_status.grid(row=5,column=1)
label_mv = tk.Label(main_window,text="Missing Value")
label_mv.grid(row=6,column=0)
label_mi = tk.Label(main_window,text="Max intensity Threshold")
label_mi.grid(row=6,column=1)
label_ticm = tk.Label(main_window,text="Peak TIC Threshold")
label_ticm.grid(row=6,column=2)

# define entrys
Entry_imzML_file = tk.Entry(main_window,width=50)
Entry_imzML_file.grid(row=0,column=2)
Entry_ppm_thresh = tk.Entry(main_window, width= 30, borderwidth=5)
Entry_ppm_thresh.insert(0,'5')
Entry_ppm_thresh.grid(row = 1,column=1)
Entry_noise_thresh = tk.Entry(main_window, width= 30, borderwidth=5)
Entry_noise_thresh.insert(0,'3000')
Entry_noise_thresh.grid(row = 1,column=4)
Entry_mv_thresh = tk.Entry(main_window, width= 30, borderwidth=5)
Entry_mv_thresh.insert(0,'0')
Entry_mv_thresh.grid(row = 7,column=0)
Entry_mi_thresh = tk.Entry(main_window, width= 30, borderwidth=5)
Entry_mi_thresh.insert(0,'0')
Entry_mi_thresh.grid(row = 7,column=1)
Entry_ticm_thresh = tk.Entry(main_window, width= 30, borderwidth=5)
Entry_ticm_thresh.insert(0,'0')
Entry_ticm_thresh.grid(row = 7,column=2)

# define buttons
button_browse = tk.Button(main_window, text= "Browse file...", command= on_browse_button_click)
button_browse.grid(row=0,column=1)
button_read = tk.Button(main_window, text= "Read Data", command= on_read_button_click)
button_read.grid(row=1,column=5)
button_classify = tk.Button(main_window, text= "ML Classify", command=on_classify_button_click)
button_classify.grid(row=2,column=0)
button_drop = tk.Button(main_window, text= "Drop", command=on_drop_button_click)
button_drop.grid(row=3,column=1)
button_update = tk.Button(main_window, text= "Update", command=on_update_button_click)
button_update.grid(row=3,column=2)
button_exp_raw = tk.Button(main_window, text= "Export Raw", command=on_exp_raw_button_click)
button_exp_raw.grid(row=4,column=1)
button_align = tk.Button(main_window, text= "Alignment", command=on_align_button_click)
button_align.grid(row=5,column=0)
button_exp_aligned = tk.Button(main_window, text= "Export Aligned", command=on_exp_aligned_button_click)
button_exp_aligned.grid(row=8,column=0)
button_exp_trimmed = tk.Button(main_window, text= "Export Trimmed", command=on_exp_trimmed_button_click)
button_exp_trimmed.grid(row=8,column=1)

# define combobox
combobox_class = ttk.Combobox(main_window,values=[],width = 50)
combobox_class.grid(row=2,column=2)
combobox_class.bind("<<ComboboxSelected>>",on_combobox_select)

#define figure location
canvas = FigureCanvasTkAgg(fig, master=main_window)
canvas.get_tk_widget().grid(row=3,column=0)

# main start from here


main_window.mainloop()
'''

file_path = 'lung-25umID-40mMjiasuanan-1.mzML'
#file_path = 'kai.mzML'
ms_file = pymzml.run.Reader(file_path)

p0.extend(read_ms_data(ms_file))


file_path = 'lung-25umID-40mMjiasuanan-2.mzML'
ms_file = pymzml.run.Reader(file_path)
p0.extend(read_ms_data(ms_file))
file_path = 'lung-25umID-40mMjiasuanan-3.mzML'
ms_file = pymzml.run.Reader(file_path)
p0.extend(read_ms_data(ms_file))

group_labels, avg_tic = classify_scans(p0,ncl=3)

timing_test("classifier")

peak_location = locate_peaks(p0)

p2 = ms_alignment(p0,ppm=5)
p3,bg_peaks = bg_removal(p2,np.arange(0,184),ratio=1)
p4 = trim_data(p3,missing_value=2,int_tot_thresh=0,int_max_thresh=0)

print("end")