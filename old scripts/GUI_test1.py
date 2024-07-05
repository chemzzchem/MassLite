from tkinter import *
#import tkinter

main_window = Tk()
main_window.title("First GUI")

# Labels
Label(main_window, text="name").grid(row = 0, column= 0)
Label(main_window, text="age").grid(row = 1, column= 0)

# Text input
my_name = Entry(main_window, width= 50 , borderwidth = 5)
my_name.grid(row = 0, column=1)

my_age = Entry(main_window, width= 50 , borderwidth = 5)
my_age.grid(row = 1, column=1)

def on_click_output():
    print(f"my name is {my_name.get()}, and my age is {my_age.get()}.")
# Buttons
Button(main_window, text = 'click here for output', command= on_click_output).grid(row = 2, column=0)

def on_click_close():
    main_window.destroy()

Button(main_window,text = 'click here to close', command= on_click_close).grid(row = 2, column=1)

main_window.mainloop()