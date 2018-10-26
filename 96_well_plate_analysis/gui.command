#!/usr/bin/env python
import os
from tkinter import Tk, Label, Button, filedialog, Entry
import analysis


class GUI:
    def __init__(self, master):
        self.master = master
        self.xtt_file = ""
        self.info_file = ""

        master.title("Analyze plate reader data")

        self.xtt_file_button = Button(master, text="Select XTT plate readout", command=self.select_xtt_file)
        self.xtt_file_button.pack()

        self.info_file_button = Button(master, text="Select text from excel", command=self.select_info_file)
        self.info_file_button.pack()

        self.analyze_button = Button(master, text="Perform analysis", command=self.analyze)
        self.analyze_button.pack()

        self.close_button = Button(master, text="Close", command=master.quit)
        self.close_button.pack()

    def select_xtt_file(self):
        fname = filedialog.askopenfilename()
        if fname:
            self.xtt_file = fname
            return

    def select_info_file(self):
        fname = filedialog.askopenfilename()
        if fname:
            self.info_file = fname
            return

    def analyze(self):
        analysis.main(xtt_file=self.xtt_file,
                      info_file=self.info_file,
                      out_file=os.path.dirname(self.info_file)+"/results.csv")


root = Tk()
my_gui = GUI(root)
root.mainloop()
