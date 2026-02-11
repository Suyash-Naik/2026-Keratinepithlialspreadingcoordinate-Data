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
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
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

def setup3d_figure(
    figsize=(7, 5.3),
    dpi=100,
    font_size=16,
    savefig_dpi=300,
    font_family="sans-serif",
    font_sans="Arial",
):
    mpl.rcParams.update({
        "figure.dpi": dpi,
        "font.size": font_size,
        "savefig.dpi": savefig_dpi,
        "font.family": font_family,
        "font.sans-serif": font_sans,
    })
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    return fig, ax

def extract_time_idx(path):
    matches = find(r"(\d+)(?=\.csv$)", path)
    if matches:
        return int(matches[-1])
    matches = find(r"TimeSrs(\d+)", path)
    if matches:
        return int(matches[-1])
    matches = find(r"(\d+)", os.path.basename(path))
    return int(matches[-1]) if matches else 0


def intensity_data_sorting(
    folders,
    earlyfiles,
    file_glob="*.csv",
    delimiter=r"[;,]",
    mean_col="Mean",
    skip_rows=10,
):
    data = []
    earlyindex = []
    current_index = 0
    for folder in folders:
        files = sorted(glob(os.path.join(folder, file_glob)))
        for file in files:
            pos_match = find(r"Pos\d{2,3}", file)
            pos = pos_match[0] if pos_match else "Pos000"

            date_match = find(r"\d{8}", file)
            date = date_match[0] if date_match else "unknown"

            fid = f"{date}_{pos}"
            if fid in earlyfiles:
                earlyindex.append(current_index)
            current_index += 1

            df = pd.read_csv(file, delimiter=delimiter, engine="python")
            df = df.iloc[1::2].reset_index(drop=True)
            if skip_rows:
                df = df.iloc[skip_rows:]

            pipmax = np.max(df[mean_col])
            pipmin = np.min(df[mean_col])

            data.append([fid, pipmin, pipmax])

    data = np.array(data, dtype=object).T
    return data, earlyindex


def build_nd_raw_dataset(
    ndfolders,
    time_minutes,
    time0,
    id_tags=None,
    file_glob="*.csv",
    skip_patterns=None,
    nd_col="%Area",
    area_col="Area",
    intensity_col="Mean",
    aggregate=None,
):


    if id_tags is None:
        id_tags = [os.path.basename(os.path.normpath(folder)) for folder in ndfolders]
    if isinstance(time_minutes, (int, float)):
        time_minutes = [time_minutes] * len(ndfolders)
    if isinstance(time0, (int, float)):
        time0 = [time0] * len(ndfolders)
    if skip_patterns is None:
        skip_patterns = []

    records = []
    for folder, tag, dt_min, t0 in zip(ndfolders, id_tags, time_minutes, time0):
        files = sorted(glob(os.path.join(folder, file_glob)), key=extract_time_idx)
        if skip_patterns:
            files = [
                f for f in files
                if not any(pat in os.path.basename(f) for pat in skip_patterns)
            ]
        for idx, file in enumerate(files):
            data = pd.read_csv(file)
            time_hpf = t0 + (idx + 1) * dt_min / 60
            if aggregate == "mean":
                record = {
                    "Sample": tag,
                    "File": os.path.basename(file),
                    "Time (hpf)": time_hpf,
                    "NetworkDensity": data[nd_col].mean(),
                    "Area": data[area_col].mean(),
                    "Intensity": data[intensity_col].mean(),
                }
                records.append(record)
            else:
                data = data.copy()
                data["Sample"] = tag
                data["File"] = os.path.basename(file)
                data["Time (hpf)"] = time_hpf
                data["NetworkDensity"] = data[nd_col]
                data["Intensity"] = data[intensity_col]
                data["Area"] = data[area_col]
                records.append(data[[
                    "Sample",
                    "File",
                    "Time (hpf)",
                    "NetworkDensity",
                    "Area",
                    "Intensity",
                ]])

    if not records:
        return pd.DataFrame()
    if aggregate == "mean":
        return pd.DataFrame.from_records(records)
    return pd.concat(records, ignore_index=True)


def plot_nd_surface(
    data,
    x_col,
    y_col,
    z_col,
    grid_res=60,
    method="linear",
    cmap="viridis",
    figsize=(7, 5.3),
    dpi=100,
    font_size=16,
    sigma=10
):
    """Create a 3D surface from long-form data using grid interpolation."""
    x = data[x_col].to_numpy()
    y = data[y_col].to_numpy()
    z = data[z_col].to_numpy()

    valid_mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x = x[valid_mask]
    y = y[valid_mask]
    z = z[valid_mask]
    if x.size == 0:
        raise ValueError("No finite points available after filtering NaNs.")

    xi = np.linspace(np.nanmin(x), np.nanmax(x), grid_res)
    yi = np.linspace(np.nanmin(y), np.nanmax(y), grid_res)
    grid_x, grid_y = np.meshgrid(xi, yi)
    grid_z = griddata((x, y), z, (grid_x, grid_y), method=method)
    grid_z = gaussian_filter(grid_z, sigma=sigma) 
    fig, ax = setup3d_figure(figsize=figsize, dpi=dpi, font_size=font_size)
    surf = ax.plot_surface(
        grid_x,
        grid_y,
        grid_z,
        cmap=cmap,
        linewidth=0,
        antialiased=True,
    )
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_zlabel(z_col)
    fig.colorbar(surf, ax=ax, shrink=0.6, aspect=12, pad=0.1)
    return fig, ax, surf

