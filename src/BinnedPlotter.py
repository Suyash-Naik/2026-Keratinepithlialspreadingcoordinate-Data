#Import necessary libraries
import os
from glob import glob
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from src.Figureplot import setup_figure
""" 

BinnedIntensityPlotter class for processing and plotting binned intensity data 
Authors: Dr. Suyash Naik
version: Plot the whole paper version baby 1.0

This class reads intensity data from CSV files, bins the data based on 
specified bin edges, and plots the binned data.

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
    def __init__(self, intensity_files, bin_edges: np.ndarray = np.arange(4, 10, 0.25), save_folder=None, timeframe: float | dict = 1.0):
        """
        intensity_files: list of csv file paths
        bin_edges: array-like, bin edges for time binning
        save_folder: str, folder to save plots
        timeframe: float (seconds per frame) or dict mapping date string -> seconds per frame
        """
        self.intensity_files = intensity_files
        self.bin_edges = bin_edges 
        self.save_folder = save_folder
        self.intensity_data = pd.DataFrame()
        self.binned_intensity_data = pd.DataFrame()
        self.treatment_indices = {}
        self.binned_intensity_data["Time (hpf)"] = []
        self.timeframe = timeframe

    def _resolve_timeframe_seconds(self, date):
        if isinstance(self.timeframe, dict):
            frame_seconds = self.timeframe.get(date, self.timeframe.get("default"))
            if frame_seconds is None:
                raise ValueError(f"No timeframe set for date {date}. Provide a mapping or a default.")
            return frame_seconds
        return self.timeframe
        
    def process_data(self,plotparam:str="Mean"):
        """
        Reads intensity files, assigns treatments, bins data, and stores results.
        """
        bin_centers = 0.5 * (self.bin_edges[1:] + self.bin_edges[:-1])
        self.binned_intensity_data["Time (hpf)"] = bin_centers
        

        for file in self.intensity_files:
            match = re.search(r"(\d{8})", os.path.basename(file))
            date = match.group(1) if match else "default"
           
            posstr=file.split("Pos")[1][0:3]
            posnumber = int(posstr)
            if "Control" in file or "control" in file:
                treatment = "Control"
            elif "MyptYsl" in file or "Myptysl" in file:
                treatment = "caMypt"
            elif "caRhoA" in file:
                treatment = "caRhoA"
            else:
                treatment = None
            if treatment is None:
                print(f"Treatment not found for position {posnumber} in file {file}. Skipping.")
                continue
            data= pd.read_csv(file)
            data["Label"] = f"{date}_Pos{posstr}"
            data["Treatment"] = treatment
            # Change the first " " column to time in hours post fertilization (hpf) [timeframe in seconds +4.5 hrs]
            frame_seconds = self._resolve_timeframe_seconds(date)
            data["Time"] = [(x - 1) * frame_seconds/60/60+4 for x in data[" "]]
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
                self.treatment_indices.setdefault(treatment, []).append(col_name)
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
        self.treatment_indices = {}
        # Extract treatment indices from the DataFrame columns
        for col in self.binned_intensity_data.columns:
            if col == "Time (hpf)":
                continue
            treatment = col.split("_")[0]
            self.treatment_indices.setdefault(treatment, []).append(col)
        # Ensure "Time (hpf)" is the first column
        if "Time (hpf)" not in self.binned_intensity_data.columns:
            raise ValueError("DataFrame must contain 'Time (hpf)' column.")

    def plot_data(self, style: str = "default", xlim=None, xticks=None, xticklabels=None, yticks_minor=None, ylim=None, colors=None, labels=None):
        "Plotting the data"
        fig, ax = setup_figure(figsize=(7, 5.3), dpi=100, font_size=24)
        if colors is None:
            colors = {
                "Control": "#83bb03",
                "K4K8MO": "#0383bb",
                "K4K8mutant": "#bb0383",
            }
        elif colors == "figure2C":
            colors={
                "Control": "#83bb03",
                "caMypt": "#0383bb",
                "caRhoA": "#bb0383",
            }
        if labels is None:
            labels = {}

        if style == "figure1d":
            ax.set_xlim(4.5, 8.5)
            ax.set_xticks(np.arange(4.5, 9.0, 0.25))
            ax.set_xticklabels(["4.5", "", "", "", "5.5", "", "", "", "6.5", "", "", "", "7.5", "", "", "", "8.5", ""])
            ax.set_yticks(np.arange(0, 470, 25), minor=True)
        else:
            ax.set_xlim(4.0, 10.0)
            ax.set_xticks(np.arange(4.0, 10.0, 0.5),minor=True)
            ax.set_yticks(np.arange(0, 1500, 25), minor=True)

        if xlim is not None:
            ax.set_xlim(*xlim)
        if xticks is not None:
            ax.set_xticks(xticks, minor=False)
        if xticklabels is not None:
            ax.set_xticklabels(xticklabels)
        if yticks_minor is not None:
            ax.set_yticks(yticks_minor, minor=True)
        if ylim is not None:
            ax.set_ylim(*ylim)

        treatment_items = list(self.treatment_indices.items())
        if not treatment_items:
            data_cols = [col for col in self.binned_intensity_data.columns if col != "Time (hpf)"]
            self.treatment_indices = {"All": data_cols}
            treatment_items = list(self.treatment_indices.items())

        color_cycle = plt.cm.tab10.colors
        for idx, (treatment, cols) in enumerate(treatment_items):
            if not cols:
                continue
            group_data = self.binned_intensity_data[cols]
            color = colors.get(treatment, color_cycle[idx % len(color_cycle)])
            label = labels.get(treatment, treatment)
            ax.plot(self.binned_intensity_data["Time (hpf)"], group_data.mean(axis=1), color=color, linewidth=2, label=label)
            ax.fill_between(
                self.binned_intensity_data["Time (hpf)"],
                group_data.mean(axis=1) - group_data.sem(axis=1),
                group_data.mean(axis=1) + group_data.sem(axis=1),
                color=color, alpha=0.3
            )
        plt.legend(loc="upper left", fontsize=16, frameon=False)
        plt.savefig(os.path.join(self.save_folder, "Binned__intensity_plot.png"), bbox_inches='tight',transparent=True)
        plt.savefig(os.path.join(self.save_folder, "Binned__intensity_plot.pdf"), bbox_inches='tight',transparent=True)
        plt.savefig(os.path.join(self.save_folder, "Binned__intensity_plot.svg"), bbox_inches='tight',transparent=True)
        plt.show()


class NetworkDensityPlotter:
    def __init__(
        self,
        folders,
        time_minutes,
        time0,
        id_tags,
        timebins,
        save_folder=None,
        skip_patterns=None,
        save_processed=False,
    ):
        self.folders = folders
        self.time_minutes = time_minutes
        self.time0 = time0
        self.id_tags = id_tags
        self.timebins = np.asarray(timebins)
        self.save_folder = save_folder
        self.skip_patterns = skip_patterns or [None] * len(folders)
        self.save_processed = save_processed
        self.nd_bins = pd.DataFrame({
            "Time (hpf)": 0.5 * (self.timebins[1:] + self.timebins[:-1])
        })
        self.nd_data = []

    @staticmethod
    def extract_time_idx(path):
        m = re.findall(r"(\d+)(?=\.csv$)", path)
        if m:
            return int(m[-1])
        m = re.findall(r"TimeSrs(\d+)", path)
        if m:
            return int(m[-1])
        m = re.findall(r"(\d+)", os.path.basename(path))
        return int(m[-1]) if m else 0

    def _build_folder_dataset(self, folder, time_minutes, time0, id_tag, skip_patterns=None):
        files = sorted(glob(os.path.join(folder, "*.csv")), key=self.extract_time_idx)
        if skip_patterns:
            files = [
                f for f in files
                if not any(p in os.path.basename(f) for p in skip_patterns)
            ]
        if not files:
            print(f"No CSV files found in {folder}")
            return None

        timelist, ndlist, ndsemlist = [], [], []
        arealist, areasemlist = [], []
        intensitylist, intensitysemlist = [], []
        for idx, file in enumerate(files):
            data = pd.read_csv(file)
            time = time0 + (idx + 1) * time_minutes / 60
            timelist.append(time)
            ndlist.append(data["%Area"].mean())
            ndsemlist.append(data["%Area"].sem())
            arealist.append(data["Area"].mean())
            areasemlist.append(data["Area"].sem())
            intensitylist.append(data["Mean"].mean())
            intensitysemlist.append(data["Mean"].sem())

        ndpd = pd.DataFrame({
            "Time (hpf)": timelist,
            "NetworkDensity": ndlist,
            "NetworkDensity_SEM": ndsemlist,
            "Area": arealist,
            "Area_SEM": areasemlist,
            "Intensity": intensitylist,
            "Intensity_SEM": intensitysemlist,
        })

        if self.save_processed:
            processed_dir = os.path.join(folder, "Processed")
            os.makedirs(processed_dir, exist_ok=True)
            ndpd.to_csv(
                os.path.join(processed_dir, f"{id_tag}_NetworkPd.csv"),
                index=False,
            )

        return ndpd

    def process_data(self):
        metrics = ["NetworkDensity", "Area", "Intensity"]
        for idx, folder in enumerate(self.folders):
            ndpd = self._build_folder_dataset(
                folder=folder,
                time_minutes=self.time_minutes[idx],
                time0=self.time0[idx],
                id_tag=self.id_tags[idx],
                skip_patterns=self.skip_patterns[idx],
            )
            if ndpd is None:
                continue
            self.nd_data.append(ndpd)

            for metric in metrics:
                col_name = f"{self.id_tags[idx]}_{metric}"
                sem_name = f"{self.id_tags[idx]}_{metric}_SEM"
                bin_means = []
                bin_sems = []
                for i in range(len(self.timebins) - 1):
                    mask = (
                        (ndpd["Time (hpf)"] > self.timebins[i])
                        & (ndpd["Time (hpf)"] < self.timebins[i + 1])
                    )
                    bin_means.append(ndpd[metric][mask].mean())
                    bin_sems.append(ndpd[metric][mask].sem())
                self.nd_bins[col_name] = bin_means
                self.nd_bins[sem_name] = bin_sems

        if self.save_folder:
            os.makedirs(self.save_folder, exist_ok=True)
            self.nd_bins.to_csv(
                os.path.join(self.save_folder, "network_density_binned.csv"),
                index=False,
            )
        return self.nd_bins

    def plot_mean(
        self,
        metric="NetworkDensity",
        xlim=None,
        xticks=None,
        xticklabels=None,
        yticks_minor=None,
        ylim=None,
        color="#83bb03",
        label=None,
    ):
        fig, ax = setup_figure(figsize=(7, 5.3), dpi=100, font_size=24)
        if xlim is not None:
            ax.set_xlim(*xlim)
        if xticks is not None:
            ax.set_xticks(xticks, minor=False)
        if xticklabels is not None:
            ax.set_xticklabels(xticklabels)
        if yticks_minor is not None:
            ax.set_yticks(yticks_minor, minor=True)
        if ylim is not None:
            ax.set_ylim(*ylim)

        metric_cols = [
            f"{tag}_{metric}" for tag in self.id_tags
            if f"{tag}_{metric}" in self.nd_bins.columns
        ]
        if not metric_cols:
            raise ValueError(f"No columns found for metric '{metric}'.")
        mean_vals = self.nd_bins[metric_cols].mean(axis=1)
        sem_vals = self.nd_bins[metric_cols].sem(axis=1)
        ax.plot(self.nd_bins["Time (hpf)"], mean_vals, color=color, linewidth=2, label=label)
        ax.fill_between(
            self.nd_bins["Time (hpf)"],
            mean_vals - sem_vals,
            mean_vals + sem_vals,
            color=color,
            alpha=0.3,
        )
        if label:
            ax.legend(loc="best", frameon=False, fontsize=16)

        if self.save_folder:
            os.makedirs(self.save_folder, exist_ok=True)
            fig.savefig(
                os.path.join(self.save_folder, f"{metric}_mean.png"),
                dpi=300,
                bbox_inches="tight",
                transparent=True,
            )
            fig.savefig(
                os.path.join(self.save_folder, f"{metric}_mean.pdf"),
                dpi=300,
                bbox_inches="tight",
                transparent=True,
            )
            fig.savefig(
                os.path.join(self.save_folder, f"{metric}_mean.svg"),
                dpi=300,
                bbox_inches="tight",
                transparent=True,
            )
        plt.show()
        return fig, ax