#Import necessary libraries
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
""" BinnedIntensityPlotter class for processing and plotting binned intensity data made by Suyash v 26/05/2025
This class reads intensity data from CSV files, bins the data based on specified bin edges, and plots the binned data.
It also allows for treatment-specific data handling and visualization.
The class includes methods for:
- Initializing with file paths, bin edges, and treatment mapping
- Processing data from CSV files to extract and bin intensity values
- Importing binned data from a DataFrame
- Plotting the binned intensity data with treatment differentiation
- Getting treatment names based on position numbers
- Averaging different experiments and treatments to create a binned intensity DataFrame
- Changing the parameter for plotting (e.g., Mean, Median, etc.) is also possible as long as it is a Series vs Value function 
"""
class BinnedIntensityPlotter:
    def __init__(self, intensity_files, bin_edges, save_folder, treatment_map,timeframe : float):
        """
        intensity_files: list of csv file paths
        bin_edges: array-like, bin edges for time binning
        save_folder: str, folder to save plots
        treatment_map: dict, e.g. {"Control": [3, 5], "K4K8MO": [6, 7, 9]}
        """
        self.intensity_files = intensity_files
        self.bin_edges = bin_edges
        self.save_folder = save_folder
        self.treatment_map = treatment_map
        self.intensity_data = pd.DataFrame()
        self.binned_intensity_data = pd.DataFrame()
        self.treatment_indices = {k: [] for k in treatment_map}
        self.binned_intensity_data["Time (hpf)"] = []
        self.timeframe = timeframe
        
    def get_treatment(self, pos_number):
        """
        Returns the treatment name for a given position number.
        """
        for treatment, pos_list in self.treatment_map.items():
            if pos_number in pos_list:
                return treatment
        return None
    
    def process_data(self,plotparam:str="Mean"):
        """
        Reads intensity files, assigns treatments, bins data, and stores results.
        """
        bin_centers = 0.5 * (self.bin_edges[1:] + self.bin_edges[:-1])
        self.binned_intensity_data["Time (hpf)"] = bin_centers
        

        for file in self.intensity_files:
            date=file.split("/")[-1].split("_")[0]
           
            posstr=file.split("Pos")[1][0:3]
            posnumber = int(posstr)
            treatment = self.get_treatment(posnumber)
            if treatment is None:
                print(f"Treatment not found for position {posnumber} in file {file}. Skipping.")
                continue
            data= pd.read_csv(file)
            data["Label"] = f"{date}_Pos{posstr}"
            data["Treatment"] = treatment
            # Change the first " " column to time in hours post fertilization (hpf) [timeframe in seconds +4.5 hrs]
            data["Time"] = [(x - 1) * self.timeframe/60/60+4 for x in data[" "]]
            plt.plot(data["Time"], data[plotparam]-data[plotparam][0], label=posstr, alpha=0.5)
            plt.legend()
            self.intensity_data = pd.concat([self.intensity_data, data])
            # Bin the data
            bin_averages = [
                    np.mean(data["Mean"][(data["Time"] > self.bin_edges[i]) & (data["Time"] < self.bin_edges[i + 1])])
                    for i in range(len(bin_centers))
                ]
            col_name = f"{treatment}_{plotparam}_{date}_Pos{posstr}"
            #check if the column is new
            if col_name not in self.binned_intensity_data.columns:
                self.binned_intensity_data[col_name] = bin_averages
                self.treatment_indices[treatment].append(col_name)
            else:
                print(f"Column {col_name} already exists in binned intensity data. Skipping.")
        # Ensure "Time (hpf)" is the first column
        self.binned_intensity_data = self.binned_intensity_data[["Time (hpf)"] + [col for col in self.binned_intensity_data.columns if col != "Time (hpf)"]]
        #export binned intensity data 
        self.binned_intensity_data.to_csv(os.path.join(self.save_folder, "binned_intensity_data.csv"), index=False)
        return self.binned_intensity_data
    
    def import_data(self, binned_intdf):
        """ Imports binned intensity data from a DataFrame.
        binned_intdf: DataFrame containing binned intensity data
        """
        self.binned_intensity_data = binned_intdf
        # Extract treatment indices from the DataFrame columns
        for col in self.binned_intensity_data.columns:
            if col.startswith("Control_"):
                self.treatment_indices["Control"].append(col)
            elif col.startswith("K4K8MO_"):
                self.treatment_indices["K4K8MO"].append(col)
        # Ensure "Time (hpf)" is the first column
        if "Time (hpf)" not in self.binned_intensity_data.columns:
            raise ValueError("DataFrame must contain 'Time (hpf)' column.")
    def plot_data(self):
        "Plotting the data"
        fig, ax = plt.subplots(figsize=(7, 5.3))
        plt.rcParams['figure.dpi'] = 100
        plt.rcParams['font.size'] = 24
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = 'Arial'
        plt.gca().spines['right'].set_color('none')
        plt.gca().spines['top'].set_color('none')
        #split self.binned_intensity_data into control and k4k8
        control_cols = self.treatment_indices.get("Control", [])
        k4k8_cols = self.treatment_indices.get("K4K8MO", [])
        control_data = self.binned_intensity_data[control_cols]
        ax.set_xlim(4.0, 10.0)
        ax.set_xticks(np.arange(4.0, 10.0, 0.5),minor=True)
        #ax.set_xticklabels(["4.5", "", "", "", "5.5", "", "", "", "6.5", "", "", "", "7.5", "", "", "", "8.5", ""])
        ax.set_yticks(np.arange(0, 1500, 25), minor=True)

        ax.plot(self.binned_intensity_data["Time (hpf)"], control_data.mean(axis=1), color="#83bb03", linewidth=2,label="Control")
           
        ax.fill_between(
                self.binned_intensity_data["Time (hpf)"],
                control_data.mean(axis=1) - control_data.sem(axis=1),
                control_data.mean(axis=1) + control_data.sem(axis=1),
                color="#83bb03", alpha=0.3
            )
        if k4k8_cols:
            k4k8_data = self.binned_intensity_data[k4k8_cols]
            ax.plot(self.binned_intensity_data["Time (hpf)"], k4k8_data.mean(axis=1), color="#ff7f00", linewidth=2,label="K4K8MO")
            ax.fill_between(
                self.binned_intensity_data["Time (hpf)"],
                k4k8_data.mean(axis=1) - k4k8_data.sem(axis=1),
                k4k8_data.mean(axis=1) + k4k8_data.sem(axis=1),
                color="#ff7f00", alpha=0.3
            )
        plt.legend(loc="upper left", fontsize=16,frameon=False)
        plt.savefig(os.path.join(self.save_folder, "Binned__intensity_plot.png"), bbox_inches='tight',transparent=True)
        plt.savefig(os.path.join(self.save_folder, "Binned__intensity_plot.pdf"), bbox_inches='tight',transparent=True)
        plt.savefig(os.path.join(self.save_folder, "Binned__intensity_plot.svg"), bbox_inches='tight',transparent=True)
        plt.show()