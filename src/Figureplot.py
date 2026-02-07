#import required libraries
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
from glob import glob
import pandas as pd
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from re import findall as find

from datetime import datetime

today=datetime.today().strftime('%Y%m%d')
print("Figureplot.py run on ", today)

def setup_figure(
    figsize=(7, 5.3),
    dpi=100,
    font_size=24,
    savefig_dpi=300,
    font_family="sans-serif",
    font_sans="Arial",
    hide_spines=True,
    xticks_minor=None,
    yticks_minor=None,
):
    mpl.rcParams.update({
        "figure.dpi": dpi,
        "font.size": font_size,
        "savefig.dpi": savefig_dpi,
        "font.family": font_family,
        "font.sans-serif": font_sans,
    })
    fig, ax = plt.subplots(figsize=figsize)
    if hide_spines:
        ax.spines["right"].set_color("none")
        ax.spines["top"].set_color("none")
    if xticks_minor is not None:
        ax.set_xticks(xticks_minor, minor=True)
    if yticks_minor is not None:
        ax.set_yticks(yticks_minor, minor=True)
    return fig, ax
