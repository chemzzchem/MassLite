import tkinter as tk
from tkinter import filedialog
import numpy as np

tx = None

def on_use_button_click():
    global tx
    tx = text_entry.get()
    print(tx)
    print(type(tx))
    result_label.config(text=tx)

# Create the main window
root = tk.Tk()
root.title("Repeat")

# Create buttons to repeat

text_entry = tk.Entry(root,width = 30, borderwidth=5)
text_entry.pack()
use_button = tk.Button(root, text="repeat", command=on_use_button_click)
use_button.pack()

# Create a label to display the result
result_label = tk.Label(root, text="", font=("Helvetica", 12), justify="left")
result_label.pack()

# Start the main event loop
root.mainloop()

print("done")