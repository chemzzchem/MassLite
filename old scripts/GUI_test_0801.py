import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

# create main window

main_window = tk.Tk()
main_window.title("MassLite GUI")

# define actions

def on_browse_button_click():
    file_imzML_path = filedialog.askopenfilename()
    if file_imzML_path:
        Entry_imzML_file.delete(0, tk.END)
        Entry_imzML_file.insert(0, file_imzML_path)

def on_read_button_click():
    file_path = Entry_imzML_file.get()
    ms_file = ImzMLParser(filename=imzML_path)
    p0 = read_msi_data(ms_file)
    
    #p0 = read_msi_data(ms_file, thresh = Entry_noise_thresh.get())
    timing_test("data read-in")

def on_combobox_select(event):
    data_type = combobox_class.get()

# define GUI grid
# define labels
label_imzML_path = tk.Label(main_window, text="Please browse for the imzML file")
label_imzML_path.grid(row=0,column=0)
label_noise_threshold = tk.Label(main_window,text="Please set noise level")
label_noise_threshold.grid(row=1,column=0)

# define entrys
Entry_imzML_file = tk.Entry(main_window,width=50)
Entry_imzML_file.grid(row=0,column=1)
Entry_noise_thresh = tk.Entry(main_window, width= 50, borderwidth=5)
Entry_noise_thresh.grid(row = 1,column=1)

# define buttons
button_read = tk.Button(main_window, text= "Read Data", command= on_read_button_click)
button_read.grid(row=2,column=0)
button_browse = tk.Button(main_window, text= "Browse file...", command= on_browse_button_click)
button_browse.grid(row=0,column=1)

# define combobox
combobox_class = ttk.Combobox(main_window,values=["Single Cell", "MS Imaging"])
combobox_class.grid(row=4,column=1)
combobox_class.bind("<<ComboboxSelected>>",on_combobox_select)

# main start from here

imzML_path = ''
main_window.mainloop()

print("end")