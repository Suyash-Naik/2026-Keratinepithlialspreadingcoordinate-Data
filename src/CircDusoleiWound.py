import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
from glob import glob
import pandas as pd
from re import findall as find
from pathlib import Path
from tqdm import tqdm
from typing import Dict, List, Optional, Tuple

def save_curves_to_csv(
    distances_common: np.ndarray,
    per_folder_curves: Dict[int, List[np.ndarray]],
    folders: List[str],
    output_dir: str,
    condition_name: str = "Condition",
):
    """
    Save per-folder curves to a CSV file with organized columns.
    
    Parameters:
      distances_common: array of common bin centers
      per_folder_curves: dict[timepoint_index] -> list of curves (one per folder)
      folders: list of folder paths (used to extract folder names)
      output_dir: directory to save the CSV
      condition_name: condition name for the output filename
    
    Returns:
      None (saves CSV directly)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize dataframe with distances
    df = pd.DataFrame({'distance': distances_common})
    
    # Extract folder identifiers - try to find meaningful names like TimeSrs\d+ or Pos\d+
    # Otherwise fall back to indices
    folder_ids = []
    for folder_path in folders:
        # Try to find TimeSrs\d+ pattern
        match = find(r'TimeSrs\d+', folder_path)
        if match:
            folder_ids.append(match[0])
        else:
            # Try to find Pos\d+ pattern
            match = find(r'Pos\d+', folder_path)
            if match:
                folder_ids.append(match[0])
            else:
                # Fall back to using index
                folder_ids.append(f"Folder_{len(folder_ids)}")
    
    # Populate columns for each folder-timepoint combination
    for t_idx in sorted(per_folder_curves.keys()):
        curves_list = per_folder_curves[t_idx]
        for folder_idx, curve in enumerate(curves_list):
            if folder_idx < len(folder_ids):
                folder_id = folder_ids[folder_idx]
            else:
                folder_id = f"Folder_{folder_idx}"
            
            col_name = f"{folder_id}_TimePoint_{t_idx}"
            df[col_name] = curve
    
    # Save to CSV
    output_path = os.path.join(output_dir, f"{condition_name}_curves.csv")
    df.to_csv(output_path, index=False)
    print(f"Saved curves to: {output_path}")
    return df


def plot_difference_curves(
    distances_common: np.ndarray,
    mean_curves_control: Dict[int, np.ndarray],
    sem_curves_control: Dict[int, np.ndarray],
    mean_curves_treatment: Dict[int, np.ndarray],
    sem_curves_treatment: Dict[int, np.ndarray],
    label_control: str = "Control",
    label_treatment: str = "Treatment",
    cmap_name: str = "spring",
    save_path: Optional[str] = None,
    time_start_sec: int = 20,
    time_step_sec: int = 20,
    time_end_sec: Optional[int] = None,
):
    """
    Plot the difference between control and treatment curves for each timepoint.
    
    Computes: difference = control - treatment
    Error propagation: SEM_diff = sqrt(SEM_control^2 + SEM_treatment^2)
    
    Parameters:
      distances_common: common distance bins
      mean_curves_control: dict[timepoint_index] -> mean control curve
      sem_curves_control: dict[timepoint_index] -> SEM control curve
      mean_curves_treatment: dict[timepoint_index] -> mean treatment curve
      sem_curves_treatment: dict[timepoint_index] -> SEM treatment curve
      label_control: label for control condition
      label_treatment: label for treatment condition
      cmap_name: colormap name for timepoint colors
      save_path: optional base path for saving figures (adds suffixes .png, .pdf, .svg)
      time_start_sec, time_step_sec, time_end_sec: control colorbar range
    """
    fig, ax = plt.subplots(figsize=(7, 5.3))
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['font.size'] = 24
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.ylim(0,0.15)
    # Get common timepoints
    common_timepoints = sorted(set(mean_curves_control.keys()) & set(mean_curves_treatment.keys()))
    
    if len(common_timepoints) == 0:
        print("Warning: No common timepoints between control and treatment")
        return
    cmap = mpl.colormaps[cmap_name]
    cmap = mpl.colors.LinearSegmentedColormap.from_list("Greys_r_trim", cmap(np.linspace(0.0, 0.8, 256)))

    # Colorbar mapping
    max_idx = max(common_timepoints)
    tmax = time_end_sec if time_end_sec is not None else (time_start_sec + max_idx * time_step_sec)
    norm = mpl.colors.Normalize(vmin=time_start_sec, vmax=tmax)
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    linestyles=['-', '--', '-.']
    for idx, t_idx in enumerate(common_timepoints):
        # Compute difference: control - treatment
        y_control = mean_curves_control[t_idx]
        y_treatment = mean_curves_treatment[t_idx]
        y_diff = y_control - y_treatment
        
        # Propagate error: sqrt(sem1^2 + sem2^2)
        sem_control = sem_curves_control[t_idx]
        sem_treatment = sem_curves_treatment[t_idx]
        sem_diff = np.sqrt(sem_control**2 + sem_treatment**2)
        
        # Time-based color
        t_seconds =  ( t_idx -1)* time_step_sec
        color = cmap(norm(t_seconds))

        # Plot
        ax.plot(distances_common, y_diff, color=color, linewidth=2,linestyle=linestyles[idx %len(linestyles) ],label="    "*(idx+1))
        ax.fill_between(distances_common, y_diff - sem_diff, y_diff + sem_diff, 
                        color=color, alpha=0.2)
    
    # Add zero reference line
    #ax.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    plt.legend(frameon=False,fontsize=18, loc="upper left")
    #ax.xlabel = 'Distance from center' if distances_common.max() > 1.0 else 'Normalized distance from center'
    #plt.xlabel(ax.xlabel)
    #plt.ylabel(f'Magnitude difference\n({label_control} - {label_treatment})')
    #plt.title(f'{label_control} vs {label_treatment}', fontsize=20)

    # Inset colorbar inside plot area (smaller legend)
    #cax = inset_axes(ax, width="3%", height="35%", loc="upper right", borderpad=0.3)
    #cbar = plt.colorbar(sm, cax=cax, orientation='vertical')
    #cbar.set_label('Time (s)')
    #cbar.ax.tick_params(labelsize=14, length=0)
    
    # Save if path provided
    if save_path:
        plt.savefig(f"{save_path}.png", dpi=300, bbox_inches='tight', transparent=True, format='png')
        plt.savefig(f"{save_path}.pdf", dpi=300, bbox_inches='tight', transparent=True, format='pdf')
        plt.savefig(f"{save_path}.svg", dpi=300, bbox_inches='tight', transparent=True, format='svg')
    
    plt.show()


def average_timepoint_curves_across_folders(
    folders: List[str],
    timepoint_indices: List[int],
    num_ranges: int = 10,
    folderdepth: int = 20,
    binning: str = "normalized",  # 'normalized' or 'absolute'
    common_max_radius: Optional[float] = None,
) -> Tuple[np.ndarray, Dict[int, np.ndarray], Dict[int, np.ndarray], Dict[int, List[np.ndarray]]]:
    """
    Compute per-timepoint average radial mean-magnitude curves across folders.

    Returns:
      distances_common: array of common bin centers (length num_ranges-1)
      mean_curves: dict[timepoint_index] -> mean curve across folders
      sem_curves: dict[timepoint_index] -> SEM across folders
      per_folder_curves: dict[timepoint_index] -> list of interpolated curves per folder
    """
    if binning not in ("normalized", "absolute"):
        raise ValueError("binning must be 'normalized' or 'absolute'")

    # Establish common bin positions (match original length num_ranges-1)
    if binning == "normalized":
        distances_common = np.linspace(0, 1.0, num_ranges)[:-1]
    else:
        if common_max_radius is None:
            raise ValueError("For absolute binning, provide common_max_radius")
        distances_common = np.linspace(0, common_max_radius, num_ranges)[:-1]

    per_timepoint_curves: Dict[int, List[np.ndarray]] = {t: [] for t in timepoint_indices}

    for folder in folders:
        analyzer = CircularAverageAnalyzer(folder)
        distances_i, curves_by_time, _ = analyzer.analyze_files(folderdepth=folderdepth, show_progress=True)
        if len(distances_i) == 0:
            continue

        if binning == "normalized":
            maxd = distances_i[-1] if distances_i[-1] != 0 else 1.0
            x_interp = distances_i / maxd
        else:
            x_interp = distances_i

        for t_idx in timepoint_indices:
            if t_idx >= len(curves_by_time):
                continue
            curve_i = np.array(curves_by_time[t_idx], dtype=float)
            y_interp = np.interp(distances_common, x_interp, curve_i)
            per_timepoint_curves[t_idx].append(y_interp)

    mean_curves: Dict[int, np.ndarray] = {}
    sem_curves: Dict[int, np.ndarray] = {}
    for t_idx in timepoint_indices:
        curves_list = per_timepoint_curves.get(t_idx, [])
        if len(curves_list) == 0:
            mean_curves[t_idx] = np.full_like(distances_common, np.nan, dtype=float)
            sem_curves[t_idx] = np.full_like(distances_common, np.nan, dtype=float)
            continue
        stack = np.vstack(curves_list)
        mean_curves[t_idx] = np.nanmean(stack, axis=0)
        sem_curves[t_idx] = pd.DataFrame(stack).sem(axis=0, skipna=True).values

    return distances_common, mean_curves, sem_curves, per_timepoint_curves


def plot_averaged_timepoint_curves(
    distances_common: np.ndarray,
    mean_curves: Dict[int, np.ndarray],
    sem_curves: Optional[Dict[int, np.ndarray]] = None,
    label_prefix: str = "TP",
    cmap_name: str = "pink",
    time_start_sec: int = 20,
    time_step_sec: int = 20,
    time_end_sec: Optional[int] = None,
):
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['font.size'] = 24
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Arial'

    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.ylim(0,0.2)
    cmap = mpl.colormaps[cmap_name]
    keys_sorted = sorted(mean_curves.keys())

    # Map timepoints to seconds: t_idx=0 -> 20s, 1 -> 40s, ...
    max_idx = max(keys_sorted) if keys_sorted else 0
    tmax = time_end_sec if time_end_sec is not None else (time_start_sec + max_idx * time_step_sec)
    norm = mpl.colors.Normalize(vmin=time_start_sec-10, vmax=tmax+10)
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    linestyles = ['-', '--', '-.', '-..']
    for idx, t_idx in enumerate(keys_sorted):
        t_seconds = time_start_sec + (t_idx) * time_step_sec
        color = cmap(norm(t_seconds))
        y = mean_curves[t_idx]
        linestyle = linestyles[idx % len(linestyles)]
        plt.plot(distances_common, y, color=color, linestyle=linestyle,label=" "*(idx+1))
        if sem_curves is not None:
            s = sem_curves.get(t_idx, None)
            if s is not None:
                plt.fill_between(distances_common, y - s, y + s, color=color, alpha=0.2)

    #plt.xlabel('Normalized distance from center' if distances_common.max() <= 1.0 else 'Distance from center')
    #plt.ylabel('Mean magnitude')
    plt.legend(frameon=False,fontsize=18,loc="upper left")
    # Inset colorbar inside plot area (smaller legend)
    #cax = inset_axes(ax, width="3%", height="35%", loc="upper right", borderpad=0.3)
    #cbar = plt.colorbar(sm, cax=cax, orientation='vertical')
    #cbar.set_label('Time (s)')
    #cbar.ax.tick_params(labelsize=14, length=0)
    outputfolder="/mnt/b/home/PHD_Data/Imaging_et_analysis/PaperFigures/Figure4/PIVtissueflow/Control_vs_MO_Difference"
    os.makedirs(outputfolder,exist_ok=True)
    plt.savefig(f"{outputfolder}/Averaged_{label_prefix}_closure.pdf",bbox_inches='tight',pad_inches=0,transparent=True,dpi=300)
    plt.savefig(f"{outputfolder}/Averaged_{label_prefix}_closure.png",bbox_inches='tight',pad_inches=0,transparent=True,dpi=300)
    plt.savefig(f"{outputfolder}/Averaged_{label_prefix}_closure.svg",bbox_inches='tight',pad_inches=0,transparent=True,dpi=300)
    plt.show()



class CircularAverageAnalyzer:
    def __init__(self, folder_path):
        self.folder_path = folder_path
        self._verbose = False

    def _load_data(self, file_path):
        """Load data from a text file."""
        return np.loadtxt(file_path, skiprows=3, delimiter=",")

    def _calculate_central_point(self,x, y, u, v):
        """Calculate the central point based on the weighted average of vector endpoints."""
        if getattr(self, "_verbose", False):
            print(os.path.basename(os.path.basename(self.folder_path)))
        if "ResultsCenter.csv" in os.listdir(self.folder_path):
            centerfile = pd.read_csv(self.folder_path+ "ResultsCenter.csv")
            central_x = centerfile['X'][0]
            central_y = centerfile['Y'][0]
        else:
            if getattr(self, "_verbose", False):
                print("\n____Calculating center based on weighted average of vector endpoints._____")
        return central_x, central_y

    def _calculate_circular_mean_magnitude(self, x, y, u, v, central_x, central_y, num_ranges=10):
        """Calculate the circular mean magnitude over distance from the center point."""
        magnitude = np.sqrt(u**2 + v**2)
        max_radius = np.max(np.sqrt((x - central_x)**2 + (y - central_y)**2))
        distances = np.linspace(0, max_radius, num_ranges)
        mean_magnitudes = []
        for i in range(num_ranges - 1):
            indices = np.where((np.sqrt((x - central_x)**2 + (y - central_y)**2) >= distances[i]) &
                               (np.sqrt((x - central_x)**2 + (y - central_y)**2) < distances[i + 1]))
            mean_magnitude = np.nanmean(magnitude[indices])
            mean_magnitudes.append(mean_magnitude)
        return distances[:-1], mean_magnitudes
    def _calculate_mean_magnitude(self, u, v):
        """Calculate the mean magnitude over all vectors."""
        magnitude = np.sqrt(u**2 + v**2)
        return np.nanmean(magnitude)

    def analyze_files(self,folderdepth=20, show_progress=True, verbose=None):
        """Analyze all text files in the folder."""
        if verbose is not None:
            self._verbose = bool(verbose)
        file_paths = sorted([os.path.join(self.folder_path, f) for f in os.listdir(self.folder_path) if f.endswith('.txt')])
        mean_magnitudes_over_distances = []
        timemagnitude=[]
        print(f"Analyzing {min(folderdepth,len(file_paths))} files in {self.folder_path}")
        for _, file_path in enumerate(tqdm(file_paths[:folderdepth], leave=False, disable=(not show_progress), dynamic_ncols=True)):
            data = self._load_data(file_path)
            x, y, u, v = data[:, 0]*0.1, data[:, 1]*0.1, data[:, 2]*1/184.656, data[:, 3]*1/184.656
            central_x, central_y = self._calculate_central_point(x, y, u, v)
            distances, mean_magnitudes = self._calculate_circular_mean_magnitude(x, y, u, v, central_x, central_y)
            mean_magnitudes_over_distances.append(mean_magnitudes)
            timemagnitude.append(self._calculate_mean_magnitude( u, v))
        return distances, mean_magnitudes_over_distances,timemagnitude
    def _average_across_experiment(self, folderdepth=20):
        """Calculate the average mean magnitudes over distances across all files."""
        distances, mean_magnitudes_over_distances, _ = self.analyze_files(folderdepth=folderdepth)
        

    def plot_mean_magnitudes_over_distances(self,show=False, show_progress=False):
        """Plot the mean magnitude over distance from the center point for each time point."""
        distances, mean_magnitudes_over_distances,timemag = self.analyze_files(show_progress=show_progress)
        fig, ax = plt.subplots(figsize=(7, 5.3))
        plt.rcParams['figure.dpi'] = 100 # 200 e.g. is really fine, but slower
        plt.rcParams['font.size'] = 24
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = 'Arial'
        plt.gca().spines['right'].set_color('none')
        plt.gca().spines['top'].set_color('none')
        cmap=mpl.colormaps['viridis']
        for i, mean_magnitudes in enumerate(mean_magnitudes_over_distances):
            plt.plot(distances, mean_magnitudes, label=f'Timepoint {i + 1}', marker='o',color=cmap(1-i/len(mean_magnitudes_over_distances)))
        if find(r'TimeSrs\d+',self.folder_path):   
            positionname=find(r'TimeSrs\d+',self.folder_path)[0]
        elif find(r'Pos\d+',self.folder_path):
            positionname=find(r'Pos\d+',self.folder_path)[0]
        plt.xlabel('Distance from center')
        plt.ylabel('Mean Magnitude')
        #plt.title('Mean Magnitude Over Distance from Center for Different Timepoints')
        #plt.legend()
        outputfolder=Path(f"{self.folder_path}/Corrected/")
        outputfolder.mkdir(parents=True, exist_ok=True)
        plt.savefig(outputfolder/f"{positionname}_20F_Mean_Magnitude_Over_Distance.png",dpi=300,bbox_inches='tight',transparent=True,format='png')
        plt.savefig(outputfolder/f"{positionname}_20F_Mean_Magnitude_Over_Distance.svg",dpi=300,bbox_inches='tight',transparent=True,format='svg')
        plt.savefig(outputfolder/f"{positionname}_20F_Mean_Magnitude_Over_Distance.pdf",dpi=300,bbox_inches='tight',transparent=True,format='pdf')
        if show:    
            plt.show()
        else:
            plt.close()
        
    def plot_mean_magnitudes(self,show=False, show_progress=False):
        """Plot the mean magnitude over consecutive timepoints."""
        _, __,timemag = self.analyze_files(show_progress=show_progress)
        fig, ax = plt.subplots(figsize=(7, 5.3))
        plt.rcParams['figure.dpi'] = 100 # 200 e.g. is really fine, but slower
        plt.rcParams['font.size'] = 24
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = 'Arial'
        plt.gca().spines['right'].set_color('none')
        plt.gca().spines['top'].set_color('none')
        if find(r'TimeSrs\d+',self.folder_path):
            positionname=find(r'TimeSrs\d+',self.folder_path)[0]
        elif find(r'Pos\d+',self.folder_path):
            positionname=find(r'Pos\d+',self.folder_path)[0]
        plt.plot(range(1, len(timemag) + 1), timemag,color="#83bb03", marker='o')
        plt.xlabel('Timepoint')
        plt.ylabel('Mean Magnitude')
        outputfolder=Path(f"{self.folder_path}/Corrected/")
        outputfolder.mkdir(parents=True, exist_ok=True)
        #plt.title('Mean Magnitude Over Consecutive Timepoints')
        plt.savefig(f"{self.folder_path}/Corrected/{positionname}_20F_Mean_Magnitude_Over_Consecutive_Timepoints.png",dpi=300,bbox_inches='tight',transparent=True,format='png')
        plt.savefig(f"{self.folder_path}/Corrected/{positionname}_20F_Mean_Magnitude_Over_Consecutive_Timepoints.svg",dpi=300,bbox_inches='tight',transparent=True,format='svg')
        plt.savefig(f"{self.folder_path}/Corrected/{positionname}_20F_Mean_Magnitude_Over_Consecutive_Timepoints.pdf",dpi=300,bbox_inches='tight',transparent=True,format='pdf')
        if show:
            plt.show()
        else:
            plt.close()


