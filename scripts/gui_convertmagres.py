#!python
# -*- coding: utf-8 -*-
import sys
import re
import math
import numpy

from magres.oldmagres import OldMagres

import tkinter.filedialog
from tkinter import *

class ConvertGUI(object):
    """
      The ConvertGUI class encapsulates the methods and window elements required for the program
    """
    
    def __init__(self, master, tot_height):

        self.top_padding = Frame(master, height=tot_height*0.1)
        self.top_padding.grid(row=0, column=0, columnspan=2)
        
        self.file_frame = Frame(master)
        self.file_name = StringVar()
        self.file_name.set("")
        self.file_entry = Entry(self.file_frame, background="white", width="15", textvariable=self.file_name)
        self.file_entry.pack(side=LEFT)
        self.file_label = Label(self.file_frame, text="   ")
        self.file_label.pack(side=LEFT)
        self.file_button = Button(self.file_frame, text="Load file", command=self.file_load)
        self.file_button.pack(side=LEFT)
        self.file_frame.grid(row=1, column=0, columnspan = 2, sticky=W)
        
        self.middle_padding = Frame(master, height=tot_height*0.07)
        self.middle_padding.grid(row=2, column=0, columnspan=2)
        
        self.castep_file_frame = Frame(master, height=tot_height*0.1)
        self.castep_file_name = StringVar()
        self.castep_file_name.set("")
        self.castep_file_entry = Entry(self.castep_file_frame, background="white", width="15", textvariable=self.castep_file_name)
        self.castep_file_entry.pack(side=LEFT)
        self.castep_file_label = Label(self.castep_file_frame, text="   ")
        self.castep_file_label.pack(side=LEFT)
        self.castep_file_button = Button(self.castep_file_frame, text="Load CASTEP file (optional)", command=self.castep_file_load)
        self.castep_file_button.pack(side=LEFT)
        self.castep_file_frame.grid(row=3, column=0, columnspan = 2, sticky=W)
        
        self.middle_padding_2 = Frame(master, height=tot_height*0.07)
        self.middle_padding_2.grid(row=4, column=0, columnspan=2)
        
        self.out_file_frame = Frame(master)
        self.out_file_name = StringVar()
        self.out_file_name.set("")
        self.out_file_entry = Entry(self.out_file_frame, background="white", width="15", textvariable=self.out_file_name)
        self.out_file_entry.pack(side=LEFT)
        self.out_file_label = Label(self.out_file_frame, text="   ")
        self.out_file_label.pack(side=LEFT)
        self.out_file_button = Button(self.out_file_frame, text="Save as", command=self.out_file_pick)
        self.out_file_button.pack(side=LEFT)
        self.out_file_frame.grid(row=5, column=0, columnspan = 2, sticky=W)
        
        self.bottom_padding = Frame(master, height=tot_height*0.07)
        self.bottom_padding.grid(row=6, column=0, columnspan=2)
        
        self.convert_frame = Frame(master)
        self.convert_button = Button(self.convert_frame, text="Convert file", command=self.convert_file)
        self.convert_button.pack(side=BOTTOM)
        self.convert_frame.grid(row=7, column=0, columnspan=2)
        
        self.status_msg = StringVar()
        self.status_msg.set("")
        self.status_label = Label(master, textvariable=self.status_msg)
        self.status_label.grid(row=8, column=0, columnspan=2)
        
    def file_load(self):
        self.file_name.set(tkinter.filedialog.askopenfilename())
        in_name = self.file_name.get()
        dot_i = in_name.rfind('.')
        self.out_file_name.set(in_name[:dot_i] + ".new" + in_name[dot_i:])
        
    def castep_file_load(self):
        self.castep_file_name.set(tkinter.filedialog.askopenfilename())
    
    def out_file_pick(self):
        self.out_file_name.set(tkinter.filedialog.asksaveasfilename())
        
    def convert_file(self):
        
        try:
            magres_file = open(self.file_name.get()).read()
        except IOError:
            self.status_msg.set("ERROR: Couldn't load input file")
            return
        
        try:
            castep_file = open(self.castep_file_name.get()).read()
        except IOError:
            castep_file = None
            
        try:
            out_file = open(self.out_file_name.get(), 'w')
        except IOError:
            self.status_msg.set("ERROR: Couldn't create output file")
            return
        
        oldmagres_file = OldMagres(magres_file, castep_file)
        
        out = oldmagres_file.as_new_format()
        
        out_file.write(str(out))
        
        if castep_file is not None:
            self.status_msg.set("Conversion successfully completed")
        else:
            self.status_msg.set("Conversion successfully completed but .castep file was not found")
        
        
if __name__ == "__main__":
    # Initialize the root window
    root = Tk()

    scr_W = root.winfo_screenwidth()
    scr_H = root.winfo_screenheight()

    scr_prc = 0.22  # Size of window in percentage of the screen

    root.title("Convert old magres")
    root.geometry("{:d}x{:d}+{:d}+{:d}".format(int(scr_W*scr_prc), int(scr_H*scr_prc), int(scr_W*(1.0-scr_prc)/2.0), int(scr_H*(1.0-scr_prc)/2.0)))

    # Create the GUI with the root window as master, then run the main loop
    main_frame = Frame(root)
    gui = ConvertGUI(main_frame, scr_H*scr_prc)
    main_frame.pack()

    root.mainloop()
