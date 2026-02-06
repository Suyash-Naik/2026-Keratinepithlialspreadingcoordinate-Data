import numpy as np
from glob import glob
import os
import matplotlib.pyplot as plt
import pandas as pd
import math
# from scipy import sparse
# from scipy.spatial.distance import cdist
import matplotlib.widgets as mpw
from matplotlib.widgets import SpanSelector
'''

The block gets a GUI (TKinter based) to select a particular folder for importing files in the next step
Folder_path is the variable with the folder name
'''
from tkinter import filedialog
from tkinter import Tk#, StringVar, Label, Button, mainloop


def load(file):
    dat = np.genfromtxt(file, delimiter=',')
    x = dat[:, 0]
    y = dat[:, 1]
    return x, y


def parametereCalc(lasp, lret, Rp=120, P=150, Rcac=367):
    eta = Rp * P / (2 * math.pi * (lasp + np.abs(lret)))
    Pc = P - (3 * math.pi * eta * lasp) / Rp
    gamma = Rp * Pc / 2
    return(eta, Pc, gamma)


class Onselect():

    def __init__(self, ax, c='red'):
        self.coords = {}
        self.ax = ax
        self.lines = []
        self.texts = []
        self.color = c
        self.xmin, self.xmax = 0., 0.
        self.slope = 1.
        self.shift = 0.

    def __call__(self, xmin, xmax):
        indmin, indmax = np.searchsorted(x, (xmin, xmax))
        indmax = min(len(x) - 1, indmax)

        thisx = x[indmin:indmax]
        thisy = y[indmin:indmax]
        m, b = np.polyfit(thisx, thisy, 1)

        if len(self.lines) > 0:
            self.lines.pop(0).remove()
            self.texts.pop(0).remove()
        line, = self.ax.plot(thisx, m * thisx + b, c=self.color)
        text = self.ax.text(thisx[0], (m * thisx + b)[0], 'y=%fx+%f' % (m, b), size=12, color=self.color)
        self.lines.append(line)
        self.texts.append(text)
        fig.canvas.draw()
        self.xmin = x[indmin]
        self.xmax = x[indmax]
        self.slope = m
        self.shift = b

    def get_fit_params(self):
        return self.xmin, self.xmax, self.slope, self.shift

# class  ViscosityMeasure():

    # def __init__():


class ButtonClickProcessor(object):
    def __init__(self, axes, label, filename, onselect1, onselect2, x, y):
        # onselect holds the info about fitted values
        self.filename = filename
        if onselect_asp is not None:
            process = self.process
        else:
            process = self.process_end
        self.button = mpw.Button(axes, label)
        self.button.on_clicked(process)
        self.onselect1 = onselect_asp
        self.onselect2 = onselect2
        self.x = x
        self.y = y

    def process(self, event):
        # when happy with the fitting, close figure
        plt.close()
        lasp = onselect_asp.get_fit_params()[2]
        lret = onselect2.get_fit_params()[2]

    def process_end(self, event):
        # close plot - we skip
        # os.rename(self.filename,self.filename+".skip")
        plt.close()


# def browse_button():
#     global folder_path
#     foldname = filedialog.askdirectory(initialdir=r'H:\PHD_data\Imaging_et_analysis\Sp5imaging\Pipettes\11_nov\26112021\Analysis\ViscoMeasurement', title='Whats up? Witch folder?')
#     folder_path.set(foldname)
#     print(foldname)
#     root.destroy()


root = Tk()
root.withdraw()
folder_path = filedialog.askdirectory()#StringVar()

# button2 = Button(root, text="Browse", command=browse_button)
# button2.grid(row=1, column=3)

# root.mainloop()

mainDir = folder_path
filename = glob(mainDir + "/*Values*.csv")
if len(filename) == 0:
    try:
        print("No files in main dir, looking in subdirs")
        filename = glob(mainDir + "/*/*Values*.csv")
    except:
        print("No files found")     
        exit()

# x, y = load(filename[0])
# fig = plt.figure()
# ax1 = fig.add_subplot(211)
# ax2 = fig.add_subplot(212)
# ax1.plot(x, y)
# ax2.plot(x, y)
# onselect_asp = Onselect(ax1, 'red')

# span_asp = SpanSelector(
#     ax1,
#     onselect=onselect_asp,
#     direction='horizontal',
#     minspan=0,
#     useblit=True,
#     span_stays=True,
#     button=1,
#     rectprops={'facecolor': 'yellow', 'alpha': 0.3}
# )
# onselect2 = Onselect(ax2, 'green')
# span_ret = SpanSelector(
#     ax2,
#     onselect=onselect2,
#     direction='horizontal',
#     minspan=0,
#     useblit=True,
#     span_stays=True,
#     button=1,
#     rectprops={'facecolor': 'yellow', 'alpha': 0.3}
# )

# # make a button to kill the plot once happy
# axdone = plt.axes([0.81, 0.05, 0.1, 0.075])
# # # if curve is really bad - skip
# axskip = plt.axes([0.51, 0.05, 0.1, 0.075])
# bnext = ButtonClickProcessor(axdone, 'Done', filename, onselect_asp, onselect2, x, y)
# bskip = ButtonClickProcessor(axskip, 'Skip', filename, None, None, None, None)
# plt.show()

# lasp = onselect_asp.get_fit_params()[2]
# lret = onselect2.get_fit_params()[2]
# print(parametereCalc(lasp, lret))
dicto = {}

for i in range(len(filename)):
    x, y = load(filename[i])
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.plot(x, y)
    ax2.plot(x, y)
    onselect_asp = Onselect(ax1, 'red')

    span_asp = SpanSelector(
        ax1,
        onselect=onselect_asp,
        direction='horizontal',
        minspan=0,
        useblit=True,
        button=1,
        props={'facecolor': 'yellow', 'alpha': 0.3}
    )
    onselect2 = Onselect(ax2, 'green')
    span_ret = SpanSelector(
        ax2,
        onselect=onselect2,
        direction='horizontal',
        minspan=0,
        useblit=True,
        button=1,
        props={'facecolor': 'green', 'alpha': 0.3}
    )

    # make a button to kill the plot once happy
    axdone = plt.axes([0.81, 0.05, 0.1, 0.075])
    # # if curve is really bad - skip
    # axskip = plt.axes([0.51, 0.05, 0.1, 0.075])
    bnext = ButtonClickProcessor(axdone, 'Done', filename, onselect_asp, onselect2, x, y)

    # bskip = ButtonClickProcessor(axskip, 'Skip', filename, None, None, None, None)
    plt.show()
    lasp = onselect_asp.get_fit_params()[2]
    lret = onselect2.get_fit_params()[2]
    Eta, Pc, Gamma = parametereCalc(lasp, lret)
    dicto[filename[i][-10:-4]] = [lasp, lret, Eta, Pc, Gamma]
    print(dicto)
#save the dictionary in the target folder
dictpd=pd.DataFrame.from_dict(dicto, orient='index', columns=['lasp', 'lret', 'eta', 'Pc', 'gamma'])
dictpd.to_csv(mainDir+'/ViscoResults_120umpip1311.csv')