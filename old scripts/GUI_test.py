import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

# create main window

main_window = tk.Tk()
main_window.title("MassLite GUI")

# define actions

def on_browse_button_click():
    file_imzML_path = filedialog.askopenfilename()
    if (file_imzML_path):
        label_imzML_path.config(text="Selected: "+ file_imzML_path)
    else:
        label_imzML_path.config(text="Please browse for the imzML file. No file selected.")
    
def on_click_close():
    main_window.destroy()

def on_combobox_select(event):
    data_type = combobox_datatype.get()

# define labels
label_imzML_path = tk.Label(main_window, text="Please browse for the imzML file")
label_imzML_path.grid(row=0,column=0)
label_noise_threshold = tk.Label(main_window,text="Please set noise level")
label_noise_threshold.grid(row=2,column=0)
label_datatype = tk.Label(main_window,text="Please select data type")
label_datatype.grid(row = 4,column=0)

# define entrys
Entry_noise_thresh = tk.Entry(main_window, width= 50, borderwidth=5)
Entry_noise_thresh.grid(row = 2,column=1)

# define buttons
browse_button = tk.Button(main_window, text= "Browse file...", command= on_browse_button_click)
browse_button.grid(row=1,column=0)
close_button = tk.Button(main_window, text= 'Close', command=on_click_close )
close_button.grid(row = 1, column=1)

# define combobox
combobox_datatype = ttk.Combobox(main_window,values=["Single Cell", "MS Imaging"])
combobox_datatype.grid(row=4,column=1)
combobox_datatype.bind("<<ComboboxSelected>>",on_combobox_select)

main_window.mainloop()