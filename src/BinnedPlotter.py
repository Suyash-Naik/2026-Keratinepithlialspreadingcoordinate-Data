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
import re
from datetime import datetime
today=datetime.today().strftime('%Y%m%d')
class BinnedIntensityPlotter:
    def __init__(self, intensity_files, bin_edges, save_folder):
        self.intensity_files = intensity_files
        self.bin_edges = bin_edges
        self.save_folder = save_folder
        self.intensity_data = pd.DataFrame()
        self.binned_intensity_data = pd.DataFrame()
        print("BinnedIntensityPlotter initialized on ", today)
    def process_data(self):
        bin_centers = 0.5 * (self.bin_edges[1:] + self.bin_edges[:-1])
        self.binned_intensity_data["Time (hpf)"] = bin_centers
        bin_averages = np.zeros(len(bin_centers))

        for file in self.intensity_files:
            if file.find("0804") == -1:
                data = pd.read_csv(file)
                data["Label"] = "23022021" + file.split("Pos")[1][0:3]
                data["Time"] = [(x - 1) * 15.5 / 60 + 4 for x in data[" "]]
                self.intensity_data = pd.concat([self.intensity_data, data])
                for i in range(len(bin_centers)):
                    bin_averages[i] = np.mean(data["Mean"][(data["Time"] > self.bin_edges[i]) & (data["Time"] < self.bin_edges[i + 1])])
                self.binned_intensity_data["23022021_" + file.split("Pos")[1][0:3]] = bin_averages

            else:
                data = pd.read_csv(file)
                data["Time"] = [(x - 1) * 5.2 / 60 + 4 for x in data[" "]]
                self.intensity_data = pd.concat([self.intensity_data, data])
                for i in range(len(bin_centers)):
                    bin_averages[i] = np.mean(data["Mean"][(data["Time"] > self.bin_edges[i]) & (data["Time"] < self.bin_edges[i + 1])])
                self.binned_intensity_data["08042021_" + file.split("Pos")[1][0:3]] = bin_averages
        self.binned_intensity_data.to_csv(self.save_folder +today+ "_Fig1D_BinnedIntensityData.csv", index=False)

    def plot_single_data(self,legend=False,save_plots=True):
        fig, ax = plt.subplots(figsize=(7, 5.3))
        plt.rcParams['figure.dpi'] = 100
        plt.rcParams['font.size'] = 24
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = 'Arial'
        plt.gca().spines['right'].set_color('none')
        plt.gca().spines['top'].set_color('none')
        ax.set_xlim(4.5, 8.5)
        ax.set_xticks(np.arange(4.5, 9, 0.25))
        ax.set_xticklabels(["4.5", "", "", "", "5.5", "", "", "", "6.5", "", "", "", "7.5", "", "", "", "8.5", ""])
        ax.set_yticks(np.arange(0, 470, 25), minor=True)
        ax.plot(self.binned_intensity_data["Time (hpf)"], self.binned_intensity_data.iloc[:, 1:5].mean(axis=1), color="#83bb03", linewidth=2)
        ax.scatter(self.binned_intensity_data["Time (hpf)"], self.binned_intensity_data.iloc[:, 1:5].mean(axis=1), color="#83bb03", s=4)
        ax.fill_between(self.binned_intensity_data["Time (hpf)"], self.binned_intensity_data.iloc[:, 1:5].mean(axis=1) - self.binned_intensity_data.iloc[:, 1:5].sem(axis=1), self.binned_intensity_data.iloc[:, 1:5].mean(axis=1) + self.binned_intensity_data.iloc[:, 1:5].sem(axis=1), color="#83bb03", alpha=0.3)
        if legend:
            ax.legend(frameon=False, loc='best', fontsize=20)
        if save_plots:
            plt.savefig(self.save_folder +today+ "_Figure1D_KeratinIntensity_nolab.png", dpi=300, bbox_inches='tight', transparent=True)
            plt.savefig(self.save_folder +today+ "_Figure1D_KeratinIntensity_nolab.pdf", dpi=300, bbox_inches='tight', transparent=True)
            plt.savefig(self.save_folder +today+ "_Figure1D_KeratinIntensity_nolab.svg", dpi=300, bbox_inches='tight', transparent=True)

        plt.show()
        print("Plots saved in ", self.save_folder)
