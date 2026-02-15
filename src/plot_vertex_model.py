import os, sys
import numpy as np
from glob import glob

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap

# export plots to csv

def export(ax, oname):
    # get lines
    out = []
    for l in ax.lines:
        out += [
            ["x - %s" % l.get_label(), *l.get_xdata()],
            ["y - %s" % l.get_label(), *l.get_ydata()]]
    # transpose
    tout = []
    for i in range(max([len(_) for _ in out])):
        tout += [[]]
        for j in range(len(out)):
            if len(out[j]) > i: tout[-1] += [out[j][i]]
            else: tout[-1] += [""]
    # export
    np.savetxt(oname, tout, delimiter=",", fmt="%s")

# DATA

script_dname = os.path.dirname(os.path.realpath(__file__))
data_dname = os.path.join(script_dname, "..", "data_vertex_model")
os.chdir(data_dname)

dt = 2. # shift in time scale in hours

alpha = 4000

# stretching simulations
k, z, dzdt, pk = {}, {}, {}, {}
for beta in (0, 0.005):
    k[beta], z[beta], dzdt[beta], pk[beta] = {}, {}, {}, {}
    #for fpull in (0.1425, 0.285, 0.57, 1.14, 2.28):
    for fpull in (0.1425, 0.285, 0.57, 1.14):
        k[beta][fpull] = np.genfromtxt(
            "keratin_beta%.3f_fpull%s.csv" % (beta, fpull),
            delimiter=",", skip_header=1)
        z[beta][fpull] = np.genfromtxt(
            "z_beta%.3f_fpull%s.csv" % (beta, fpull),
            delimiter=",", skip_header=1)
        dzdt[beta][fpull] = np.genfromtxt(
            "dzdt_beta%.3f_fpull%s.csv" % (beta, fpull),
            delimiter=",", skip_header=1)
        try:
            pk[beta][fpull] = np.genfromtxt(
                "pressure_keratin.beta%.3f.fpull%s.csv" % (beta, fpull),
                delimiter=",", skip_header=1)
        except FileNotFoundError: pass
        # shift time
        k[beta][fpull][:, 0] -= dt
        z[beta][fpull][:, 0] -= dt
        dzdt[beta][fpull][:, 0] -= dt

# cutting simulation
tcut, acut, kcut = {}, {}, {}
for beta, label in zip((0, 0.005), ("without", "with")):
    tcut[beta], acut[beta], kcut[beta] = np.transpose(np.genfromtxt(
        "closure_%s.csv" % label,
        delimiter=",", skip_header=1))
vel_row, disp_row, ang_row = {}, {}, {}
for title in (
    "closure_without", "closure_with", "recoil_without", "recoil_with"):
    vel_row[title], disp_row[title], ang_row[title] = {}, {}, {}
    for fname in glob("vel_row.%s.t*.csv" % title):
        time = float(fname.split(".t")[1].split(".csv")[0])
        vel_row[title][time] = np.genfromtxt(fname,
            skip_header=1, delimiter=",")
        disp_row[title][time] = np.genfromtxt(fname.replace("vel", "disp"),
            skip_header=1, delimiter=",")
        ang_row[title][time] = np.genfromtxt(fname.replace("vel", "ang"),
            skip_header=1, delimiter=",")

# experimental keratin distribution
dist_exp, sdist_exp = {}, {}
for fname in glob("Results_Pos001_19012023_F*.dist.csv"):
    index = int(fname.split("F")[1].split(".")[0])
    time = (615.45*index)/3600 + 4.5                            # shift time
    dist_exp[time] = np.genfromtxt(fname,
        delimiter=",", skip_header=1)
    sdist_exp[time] = np.genfromtxt(fname.replace("dist", "sdist"),
        delimiter=",", skip_header=1)

# numerical keratin distribution
dist_sim, heter_dist_sim, sdist_sim, heter_sdist_sim = {}, {}, {}, {}
for fname in glob("0.57.t*.dist.csv"):
    time = float(fname.split(".t")[1].split(".dist")[0]) - dt   # shift time
    if time < 2: continue
    dist_sim[time] = np.genfromtxt(
        fname,
        delimiter=",", skip_header=1)
    sdist_sim[time] = np.genfromtxt(
        fname.replace("dist", "sdist"),
        delimiter=",", skip_header=1)
    heter_dist_sim[time] = np.genfromtxt(
        "hetera.%s" % fname,
        delimiter=",", skip_header=1)
    heter_sdist_sim[time] = np.genfromtxt(
        ("hetera.%s" % fname).replace("dist", "sdist"),
        delimiter=",", skip_header=1)

# numerical gradients
hist_tension_radius_sim, hist_keratin_radius_sim, hist_pressure_radius_sim, \
    hist_keratin_area_sim = {}, {}, {}, {}
for fname in glob("0.57.hist_tension_radius.t*.csv"):
    time = float(fname.split(".t")[1].split(".csv")[0]) - dt    # shift time
    if time < 2: continue
    hist_tension_radius_sim[time] = np.genfromtxt(
        fname,
        delimiter=",", skip_header=1)
    hist_keratin_radius_sim[time] = np.genfromtxt(
        fname.replace("tension_radius", "keratin_radius"),
        delimiter=",", skip_header=1)
    hist_pressure_radius_sim[time] = np.genfromtxt(
        fname.replace("tension_radius", "pressure_radius"),
        delimiter=",", skip_header=1)
    hist_pressure_radius_sim[time][:, 1:] = \
        -hist_pressure_radius_sim[time][:, 1:]/alpha
    hist_keratin_area_sim[time] = np.genfromtxt(
        fname.replace("tension_radius", "keratin_area"),
        delimiter=",", skip_header=1)
assert sorted(dist_sim.keys()) == sorted(hist_tension_radius_sim.keys())

# keratin ratio of standard deviation and mean
stat_exp = np.genfromtxt("stat_ker_exp.csv", delimiter=",", skip_header=1)
stat_sim = {
    float(fname.split(".stat")[0]):
        np.genfromtxt(fname, delimiter=",", skip_header=1)
    for fname in glob("*.stat_ker_sim.csv")}
for fpull in stat_sim: stat_sim[fpull][:, 0] -= dt              # shift time

# pipette experiment
keratin_pipette, height_pipette = {}, {}
for fname in glob("*_keratin_pipette.csv"):
    time = float(fname.split("_keratin")[0])
    keratin_pipette[time] = np.genfromtxt(fname,
        delimiter=",", skip_header=True)
    height_pipette[time] = np.genfromtxt(fname.replace("keratin", "height"),
        delimiter=",", skip_header=True)
t_pip, h0, v0, vm, vp, taur, taurp, tau, taup, tauK, xi, springk, T = \
    np.transpose(np.genfromtxt(
        "pipette_values.csv", delimiter=",", skip_header=1,
        usecols=(0, 1, 2, 3, 4, 5, 6, 9, 10, 7, 11, 12, 13)))

# PLOT PARAMETERS

# style
plt.rcParams["figure.figsize"] = 7, 5.3
plt.rcParams["figure.dpi"] = 100
plt.rcParams["font.size"] = 24
plt.rcParams["font.weight"] = "normal"
plt.rcParams["axes.labelsize"] = 24
plt.rcParams["xtick.labelsize"] = 24
plt.rcParams["ytick.labelsize"] = 24
plt.rcParams["legend.fontsize"] = 16
plt.rcParams["savefig.dpi"] = 300
# plt.rcParams["font.family"] = "sans-serif"
# plt.rcParams["font.sans-serif"] = "Arial"
# plt.rcParams["mathtext.default"] = "regular"
# plt.rcParams["mathtext.fontset"] = "stixsans"
plt.rcParams["axes.prop_cycle"] = 'cycler(color='\
    '["#83BB03", "#0383bb", "#BB9703", "#D51B66", "#5C2352", "#0173B2"])'
plt.rcParams["legend.fancybox"] = False
plt.rcParams["legend.frameon"] = False
plt.rcParams["legend.borderaxespad"] = 0
plt.rcParams["legend.handletextpad"] = 0.4
plt.rcParams["legend.labelspacing"] = 0.25
plt.rcParams["legend.handlelength"] = 1.125
plt.rcParams["legend.numpoints"] = 1
plt.rcParams["lines.linewidth"] = 4
plt.rcParams["xtick.top"] = False
plt.rcParams["ytick.right"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["xtick.direction"] = "out"
plt.rcParams["ytick.direction"] = "out"
plt.rcParams["figure.subplot.left"] = 0.195
plt.rcParams["figure.subplot.bottom"] = 0.15
plt.rcParams["figure.subplot.top"] = 0.95
plt.rcParams["text.usetex"] = False

# Arial font
import matplotlib.font_manager as font_manager
fpath = os.path.join(os.path.dirname(__file__), "Arial.ttf")
font_manager.fontManager.addfont(fpath)
prop = font_manager.FontProperties(fname=fpath)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = prop.get_name()

# rainbow colourmap
# rainbow_cmap = LinearSegmentedColormap.from_list("rainbow", (
#     (0/5, (0.467, 0.000, 0.533)),
#     (1/5, (0.000, 0.298, 1.000)),
#     (2/5, (0.008, 0.506, 0.129)),
#     (3/5, (1.000, 0.933, 0.000)),
#     (4/5, (1.000, 0.553, 0.000)),
#     (5/5, (0.898, 0.000, 0.000))))
# rainbow_cmap = plt.cm.viridis
sys.path.append(script_dname)
from batlow import batlow_map as rainbow_cmap

# OUTPUT DIRECTORY

from shutil import rmtree
output_dname = os.path.join(data_dname, "output")
try: rmtree(output_dname)
except FileNotFoundError: pass
os.mkdir(output_dname)
os.chdir(output_dname)

# CUSTOM

# velocities beta = 0

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (h)")
ax.set_ylabel(r"Velocity (µm/min)")
# ax.set_ylim([0, 18])

beta = 0
for fpull in sorted(dzdt[beta]):
    time, meandzdt, stddzdt = np.transpose(dzdt[beta][fpull])
    l, = ax.plot(time, meandzdt,
        label=r"$F_{\mathrm{YSL}}$=%sµN" % fpull)
    ax.fill_between(time, meandzdt - stddzdt/2, meandzdt + stddzdt/2,
        color=l.get_color(), alpha=0.3)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper left",
    handles=[
        Line2D([0], [0], lw=0, label=r"$\beta=%s$" % beta),
        *ax.get_legend_handles_labels()[0]]))

fig.savefig("velocities_beta0.pdf")
export(ax, "velocities_beta0.csv")

# velocities beta = 0.005

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (h)")
ax.set_ylabel(r"Velocity (µm/min)")
# ax.set_ylim([0, 18])

beta = 0.005
for fpull in sorted(dzdt[beta]):
    time, meandzdt, stddzdt = np.transpose(dzdt[beta][fpull])
    l, = ax.plot(time, meandzdt,
        label=r"$F_{\mathrm{YSL}}$=%sµN" % fpull)
    ax.fill_between(time, meandzdt - stddzdt/2, meandzdt + stddzdt/2,
        color=l.get_color(), alpha=0.3)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper right",
    handles=[
        Line2D([0], [0], lw=0, label=r"$\beta=%s$" % beta),
        *ax.get_legend_handles_labels()[0]]))

fig.savefig("velocities_beta0.005.pdf")
export(ax, "velocities_beta0.005.csv")

# velocities sm

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (h)")
ax.set_ylabel(r"Velocity (µm/min)")

l = {}
for beta in (0, 0.005):
    l[beta] = []
    for fpull in (0.1425, 0.285, 1.14):
        #if beta == 0 and fpull == 1.14: continue
        time, meandzdt, stddzdt = np.transpose(dzdt[beta][fpull])
        l[beta] += [ax.plot(time, meandzdt,
            label=r"$F_{\mathrm{YSL}}$=%sµN" % fpull,
            linestyle="--" if beta == 0 else "-")[0]]
        ax.fill_between(time, meandzdt - stddzdt/2, meandzdt + stddzdt/2,
            color=l[beta][-1].get_color(), alpha=0.3)

plt.sca(ax)
handles = []
for beta in sorted(l):
    handles += [Line2D([0], [0], lw=0, label=r"$\beta$=%s" % beta), *l[beta]]
ax.add_artist(plt.legend(loc="upper left",
    handles=handles))

fig.savefig("velocities_sm.pdf")
export(ax, "velocities_sm.csv")

# distribution exp

fig, ax = plt.subplots()
ax.set_xlabel(r"Keratin intensity")
ax.set_ylabel(r"Probability")

for i, time in enumerate(sorted(dist_exp)):
    ax.plot(dist_exp[time][:, 0], dist_exp[time][:, 1],
        color=rainbow_cmap(i/len(dist_exp)), label=r"$T=%.1f$h" % time)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper right"))

fig.savefig("dist_exp.pdf")
export(ax, "dist_exp.csv")

# distribution sim

fig, ax = plt.subplots()
ax.set_xlabel(r"Keratin intensity")
ax.set_xlim([0, 600])
ax.set_ylabel(r"Probability")

for i, time in enumerate(sorted(dist_sim)):
    ax.plot(dist_sim[time][:, 0], dist_sim[time][:, 1],
        color=rainbow_cmap(i/len(dist_sim)), label=r"$T=%.1f$h" % time)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper right"))

fig.savefig("dist_sim.pdf")
export(ax, "dist_sim.csv")
ylim = ax.get_ylim()

# distribution heterogeneous sim

fig, ax = plt.subplots()
ax.set_xlabel(r"Keratin intensity")
ax.set_xlim([0, 600])
ax.set_ylabel(r"Probability")
ax.set_ylim(ylim)

for i, time in enumerate(sorted(dist_sim)):
    ax.plot(heter_dist_sim[time][:, 0], heter_dist_sim[time][:, 1],
        color=rainbow_cmap(i/len(heter_dist_sim)), label=r"$T=%.1f$h" % time)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper right"))

fig.savefig("heter_dist_sim.pdf")
export(ax, "heter_dist_sim.csv")

# PLOTS STRETCH

for beta in (0, 0.005):

    # keratin

    fig, ax = plt.subplots()
    ax.set_xlabel(r"Time (h)")
    ax.set_ylabel(r"Keratin intensity")
    if beta != 0: ax.set_ylim([0, 800])

    for fpull in sorted(k[beta]):
        time, meank, stdk = np.transpose(k[beta][fpull])
        l, = ax.plot(time, meank,
            label=r"$F_{\mathrm{YSL}}$=%sµN" % fpull)
        ax.fill_between(time, meank - stdk/2, meank + stdk/2,
            color=l.get_color(), alpha=0.3)
    ax.axhline(y=150, color="black", linestyle="--")

    plt.sca(ax)
    ax.add_artist(plt.legend(loc="upper left"))

    fig.savefig("keratin_beta%.3f.pdf" % beta)
    export(ax, "keratin_beta%.3f.csv" % beta)

    # height

    fig, ax = plt.subplots()
    ax.set_xlabel(r"Time (h)")
    ax.set_ylabel(r"Height (mm)")

    for fpull in sorted(z[beta]):
        time, meanz, stdz = np.transpose(z[beta][fpull])
        l, = ax.plot(time, meanz,
            label=r"$F_{\mathrm{YSL}}$=%sµN" % fpull)
        ax.fill_between(time, meanz - stdz/2, meanz + stdz/2,
            color=l.get_color(), alpha=0.3)

    plt.sca(ax)
    ax.add_artist(plt.legend(loc="upper left"))

    fig.savefig("z_beta%.3f.pdf" % beta)
    export(ax, "z_beta%.3f.csv" % beta)

    # velocity

    fig, ax = plt.subplots()
    ax.set_xlabel(r"Time (h)")
    ax.set_ylabel(r"Velocity (µm/min)")
    ax.set_ylim([0, 8])

    for fpull in sorted(dzdt[beta]):
        time, meandzdt, stddzdt = np.transpose(dzdt[beta][fpull])
        l, = ax.plot(time, meandzdt,
            label=r"$F_{\mathrm{YSL}}$=%sµN" % fpull)
        ax.fill_between(time, meandzdt - stddzdt/2, meandzdt + stddzdt/2,
            color=l.get_color(), alpha=0.3)

    plt.sca(ax)
    ax.add_artist(plt.legend(loc="upper right",
        handles=[
            Line2D([0], [0], lw=0),
            Line2D([0], [0], lw=0, label=r"$\beta=%s$" % beta),
            *ax.get_legend_handles_labels()[0]]))

    fig.savefig("dzdt_beta%.3f.pdf" % beta)
    export(ax, "dzdt_beta%.3f.csv" % beta)

# pressure vs. keratin

fig, ax = plt.subplots()
ax.set_xlabel(r"$\alpha p_i$")
ax.set_ylabel(r"$K_i$")

for beta, marker, color in zip((0, 0.005), ("o", "s"), ("#83BB03", "#ff7f0e")):
    fpull = 0.57
    ax.scatter(-pk[beta][fpull][:, 0], pk[beta][fpull][:, 1],
        marker=marker, color=color, label=r"$\beta=%s$" % beta, s=50,
        rasterized=True)

ax.set_xlim(ax.get_xlim())
ax.set_ylim(ax.get_ylim())
ax.plot(ax.get_xlim(), ax.get_xlim(),
    color="black", linestyle="--", label=r"$\alpha p_i = K_i$")

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper left"))

fig.savefig("pressure_vs_keratin.pdf")
export(ax, "pressure_vs_keratin.csv")

# tension gradient

fig, ax = plt.subplots()
ax.set_xlabel(r"Scaled radius")
ax.set_ylabel(r"Tension")

for i, time in enumerate(sorted(hist_tension_radius_sim)):
    colour = rainbow_cmap(i/len(hist_tension_radius_sim))
    ax.plot(
        hist_tension_radius_sim[time][:, 0],
        hist_tension_radius_sim[time][:, 1],
        color=colour,
        label=r"$T=%.1f$h" % time)
    ax.fill_between(
        hist_tension_radius_sim[time][:, 0],
        hist_tension_radius_sim[time][:, 1]
            - hist_tension_radius_sim[time][:, 2]/2,
        hist_tension_radius_sim[time][:, 1]
            + hist_tension_radius_sim[time][:, 2]/2,
        color=colour, alpha=0.2)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper right", ncols=3))

fig.savefig("hist_tension_radius_sim.pdf")
export(ax, "hist_tension_radius_sim.csv")

# keratin gradient

fig, ax = plt.subplots()
ax.set_xlabel(r"Scaled radius")
ax.set_ylabel(r"Keratin intensity")

for i, time in enumerate(sorted(hist_keratin_radius_sim)):
    colour = rainbow_cmap(i/len(hist_keratin_radius_sim))
    ax.plot(
        hist_keratin_radius_sim[time][:, 0],
        hist_keratin_radius_sim[time][:, 1],
        color=colour,
        label=r"$T=%.1f$h" % time)
    ax.fill_between(
        hist_keratin_radius_sim[time][:, 0],
        hist_keratin_radius_sim[time][:, 1]
            - hist_keratin_radius_sim[time][:, 2]/2,
        hist_keratin_radius_sim[time][:, 1]
            + hist_keratin_radius_sim[time][:, 2]/2,
        color=colour, alpha=0.2)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper left", ncols=3))

fig.savefig("hist_keratin_radius_sim.pdf")
export(ax, "hist_keratin_radius_sim.csv")

# pressure gradient

fig, ax = plt.subplots()
ax.set_xlabel(r"Scaled radius")
ax.set_ylabel(r"Stress (x$10^6$ Pa)")

scale = 30
for i, time in enumerate(sorted(hist_pressure_radius_sim)):
    colour = rainbow_cmap(i/len(hist_pressure_radius_sim))
    ax.plot(
        hist_pressure_radius_sim[time][:, 0],
        scale*hist_pressure_radius_sim[time][:, 1],
        color=colour,
        label=r"$T=%.1f$h" % time)
    ax.fill_between(
        hist_pressure_radius_sim[time][:, 0],
        scale*(hist_pressure_radius_sim[time][:, 1]
            - hist_pressure_radius_sim[time][:, 2]/2),
        scale*(hist_pressure_radius_sim[time][:, 1]
            + hist_pressure_radius_sim[time][:, 2]/2),
        color=colour, alpha=0.2)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper left", ncols=3))

fig.savefig("hist_pressure_radius_sim.pdf")
export(ax, "hist_pressure_radius_sim.csv")

# keratin on area gradient

fig, ax = plt.subplots()
ax.set_xlabel(r"Scaled area")
ax.set_ylabel(r"Keratin intensity")

for i, time in enumerate(sorted(hist_keratin_area_sim)):
    colour = rainbow_cmap(i/len(hist_keratin_area_sim))
    ax.plot(
        hist_keratin_area_sim[time][:, 0],
        hist_keratin_area_sim[time][:, 1],
        color=colour,
        label=r"$T=%.1f$h" % time)
    ax.fill_between(
        hist_keratin_area_sim[time][:, 0],
        hist_keratin_area_sim[time][:, 1]
            - hist_keratin_area_sim[time][:, 2]/2,
        hist_keratin_area_sim[time][:, 1]
            + hist_keratin_area_sim[time][:, 2]/2,
        color=colour, alpha=0.2)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper left", ncols=3))

fig.savefig("hist_keratin_area_sim.pdf")
export(ax, "hist_keratin_area_sim.csv")

# PLOTS CUT

# height

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (s)")
ax.set_xlim([0, 100])
ax.set_ylabel(r"Hole area (µm$^{\mathdefault{2}}$)", labelpad=-8.75)

for beta in sorted(acut)[::-1]:
    ax.plot(tcut[beta], acut[beta],
        label=r"$\beta=%s$" % beta)

plt.sca(ax)
ax.add_artist(plt.legend())

fig.savefig("area_cut.pdf")
export(ax, "area_cut.csv")

# velocity closure

for title in ("closure_with", "closure_without"):

    beta = 0 if "without" in title else 0.005
    gamma = 0 if "recoil" in title else 2

    fig, ax = plt.subplots()
    ax.set_xlabel(r"Row index")
    ax.set_xlim([0, 3])
    ax.set_xticks(range(4))
    ax.set_xticklabels(range(1, 5))
    ax.set_ylabel(r"Average velocity (µm/s)")
    ax.set_ylim([0, 0.5])

    axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47],
        xlim=ax.get_xlim(), ylim=[0, 4])
    axins.set_xticks(ax.get_xticks())
    axins.set_xticklabels(ax.get_xticklabels())
    axins.set_yticks(range(5))

    for t in (
        (sorted(vel_row[title])[0],) if "recoil" in title else
        sorted(vel_row[title])[::2]):
        c = rainbow_cmap(
            list(sorted(vel_row[title])).index(t)/(len(vel_row[title]) - 1))
        ax.plot(vel_row[title][t],
            color=c, marker="o", label=r"t=%is" % round(t))
        axins.plot(vel_row[title][t],
            color=c, marker="o")

    plt.sca(ax)
    ax.add_artist(plt.legend(loc="lower left"))
    ax.add_artist(plt.legend(loc="upper right", fontsize="small",
        handles=[
            Line2D([0], [0], lw=0, label=r"β=%s, γ=%s" % (beta, gamma))]))

    fig.savefig("vel_row.%s.pdf" % title)
    export(ax, "vel_row.%s.csv" % title)

# velocity recoil

fig, ax = plt.subplots()
ax.set_xlabel(r"Row index")
ax.set_xlim([0, 3])
ax.set_xticks(range(4))
ax.set_xticklabels(range(1, 5))
ax.set_ylabel(r"Average velocity (µm/s)")
# ax.set_ylim([0, 0.5])

# axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47],
#     xlim=ax.get_xlim(), ylim=[0, 4])
# axins.set_xticks(ax.get_xticks())
# axins.set_yticks(range(5))

for title in ("recoil_with", "recoil_without"):
    beta = 0 if "without" in title else 0.005
    for t in (sorted(vel_row[title])[0],):
        try: assert t0 == t
        except NameError: t0 = t
        l, = ax.plot(vel_row[title][t],
            marker="o", label=r"β=%s" % beta)
#         axins.plot(vel_row[title][t],
#             color=l.get_color(), marker="o")

plt.sca(ax)
ax.add_artist(plt.legend(loc="lower left"))
ax.add_artist(plt.legend(loc="upper right", fontsize="small",
    handles=[
        Line2D([0], [0], lw=0, label=r"γ=0, t=%is" % round(t0))]))
del t0

fig.savefig("vel_row.recoil.pdf")
export(ax, "vel_row.recoil.csv")

# displacement closure

for title in ("closure_with", "closure_without"):

    beta = 0 if "without" in title else 0.005
    gamma = 0 if "recoil" in title else 2

    fig, ax = plt.subplots()
    ax.set_xlabel(r"Row index")
    ax.set_xlim([0, 3])
    ax.set_xticks(range(4))
    ax.set_xticklabels(range(1, 5))
    ax.set_ylabel(r"displacement (µm)")
    ax.set_ylim([0, 35])

    for t in (
        (sorted(disp_row[title])[0],) if "recoil" in title else
        sorted(disp_row[title])[::2]):
        c = rainbow_cmap(
            list(sorted(disp_row[title])).index(t)/(len(disp_row[title]) - 1))
        ax.plot(disp_row[title][t],
            color=c, marker="o", label=r"t=%is" % round(t))

    plt.sca(ax)
    ax.add_artist(plt.legend(loc="center right"))
    ax.add_artist(plt.legend(loc="upper right", fontsize="small",
        handles=[
            Line2D([0], [0], lw=0, label=r"β=%s, γ=%s" % (beta, gamma))]))

    fig.savefig("disp_row.%s.pdf" % title)
    export(ax, "disp_row.%s.csv" % title)

# difference displacement closure

fig, ax = plt.subplots()
ax.set_xlabel(r"Row index")
ax.set_xlim([0, 3])
ax.set_xticks(range(4))
ax.set_xticklabels(range(1, 5))
ax.set_ylabel(r"displacement diff.""\n"r"[WT - MO] (µm)")

ax.axhline(y=0, color="black", linestyle="--")
for t in sorted(disp_row[title])[::2]:
    c = rainbow_cmap(
        list(sorted(disp_row[title])).index(t)/(len(disp_row[title]) - 1))
    ax.plot(disp_row["closure_with"][t] - disp_row["closure_without"][t],
        color=c, marker="o", label=r"t=%is" % round(t))

plt.sca(ax)
ax.add_artist(plt.legend(loc="lower right"))
ax.add_artist(plt.legend(loc="upper center", fontsize="small",
    handles=[
        Line2D([0], [0], lw=0, label=r"γ=2")]))

fig.savefig("diff_disp_row.closure.pdf")
export(ax, "diff_disp_row.closure.csv")

# displacement recoil

fig, ax = plt.subplots()
ax.set_xlabel(r"Row index")
ax.set_xlim([0, 3])
ax.set_xticks(range(4))
ax.set_xticklabels(range(1, 5))
ax.set_ylabel(r"displacement (µm)")

for title in ("recoil_with", "recoil_without"):
    beta = 0 if "without" in title else 0.005
    for t in (sorted(disp_row[title])[0],):
        try: assert t0 == t
        except NameError: t0 = t
        l, = ax.plot(disp_row[title][t],
            marker="o", label=r"β=%s" % beta)

plt.sca(ax)
ax.add_artist(plt.legend(loc="lower left"))
ax.add_artist(plt.legend(loc="upper right", fontsize="small",
    handles=[
        Line2D([0], [0], lw=0, label=r"γ=0, t=%is" % round(t0))]))
del t0

fig.savefig("disp_row.recoil.pdf")
export(ax, "disp_row.recoil.csv")

# difference displacement recoil

fig, ax = plt.subplots()
ax.set_xlabel(r"Row index")
ax.set_xlim([0, 3])
ax.set_xticks(range(4))
ax.set_xticklabels(range(1, 5))
ax.set_ylabel(r"displacement diff.""\n"r"[WT - MO] (µm)")

ax.axhline(y=0, color="black", linestyle="--")
for t in sorted(disp_row[title])[::2]:
    c = rainbow_cmap(
        list(sorted(disp_row[title])).index(t)/(len(disp_row[title]) - 1))
    ax.plot(disp_row["recoil_with"][t] - disp_row["recoil_without"][t],
        color=c, marker="o", label=r"t=%is" % round(t))

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper right"))
ax.add_artist(plt.legend(loc="upper center", fontsize="small",
    handles=[
        Line2D([0], [0], lw=0, label=r"γ=0")]))

fig.savefig("diff_disp_row.recoil.pdf")
export(ax, "diff_disp_row.recoil.csv")

# angle closure

fig, ax = plt.subplots()
ax.set_xlabel(r"Row index")
ax.set_xlim([0, 3])
ax.set_xticks(range(4))
ax.set_xticklabels(range(1, 5))
ax.set_ylabel(r"Alignment angle (°)")
ax.set_ylim([0, 90])

ax.axhline(y=45, color="black", linestyle="--")

for title in ("closure_with", "closure_without"):
    beta = 0 if "without" in title else 0.005
    for t in sorted(ang_row[title]):
        if round(t) != 60: continue
        try: assert t0 == t
        except NameError: t0 = t
        l, = ax.plot(ang_row[title][t],
            marker="o", label=r"β=%s" % beta)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper left"))
ax.add_artist(plt.legend(loc="upper right", fontsize="small",
    handles=[
        Line2D([0], [0], lw=0, label=r"γ=2, t=%is" % round(t0))]))
del t0

fig.savefig("ang_row.closure.pdf")
export(ax, "ang_row.closure.csv")

# angle recoil

fig, ax = plt.subplots()
ax.set_xlabel(r"Row index")
ax.set_xlim([0, 3])
ax.set_xticks(range(4))
ax.set_xticklabels(range(1, 5))
ax.set_ylabel(r"Alignment angle (°)")
ax.set_ylim([0, 90])

ax.axhline(y=45, color="black", linestyle="--")

for title in ("recoil_with", "recoil_without"):
    beta = 0 if "without" in title else 0.005
    for t in (sorted(ang_row[title])[0],):
        try: assert t0 == t
        except NameError: t0 = t
        l, = ax.plot(ang_row[title][t],
            marker="o", label=r"β=%s" % beta)

plt.sca(ax)
ax.add_artist(plt.legend(loc="lower left"))
ax.add_artist(plt.legend(loc="upper right", fontsize="small",
    handles=[
        Line2D([0], [0], lw=0, label=r"γ=0, t=%is" % round(t0))]))
del t0

fig.savefig("ang_row.recoil.pdf")
export(ax, "ang_row.recoil.csv")

# PLOTS PIPETTE

plt.rcParams["axes.prop_cycle"] = \
    'cycler(color=["0173b2","de8f05","029e73","d55e00"])'

cmap_pipette = rainbow_cmap
norm_pipette = Normalize(4, 9)
scalarMap_pipette = ScalarMappable(cmap=cmap_pipette, norm=norm_pipette)

# keratin

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (s)")
ax.set_ylabel(r"Keratin $K$ (a.u.)")

for i, time in enumerate(sorted(keratin_pipette)):
    if i == 1: continue
    ax.plot(keratin_pipette[time][:, 0], keratin_pipette[time][:, 1],
        color=scalarMap_pipette.to_rgba(time))

cax = make_axes_locatable(ax).append_axes(
    "top", size="5%", pad=0.5)
cbar = ColorbarBase(
    cax, cmap=cmap_pipette, norm=norm_pipette, orientation="horizontal")
cax.xaxis.set_label_position("top")
cax.xaxis.set_ticks_position("bottom")
cbar.set_label(r"Hours post fertilisation", labelpad=5)
ax.figure.subplots_adjust(top=0.925)

fig.savefig("keratin_pipette.pdf")
export(ax, "keratin_pipette.csv")

# height

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (s)")
ax.set_ylabel(r"Height $\ell$ (µm)")

for i, time in enumerate(sorted(keratin_pipette)):
    if i == 1: continue
    ax.plot(height_pipette[time][:, 0], height_pipette[time][:, 1],
        color=scalarMap_pipette.to_rgba(time))

cax = make_axes_locatable(ax).append_axes(
    "top", size="5%", pad=0.5)
cbar = ColorbarBase(
    cax, cmap=cmap_pipette, norm=norm_pipette, orientation="horizontal")
cax.xaxis.set_label_position("top")
cax.xaxis.set_ticks_position("bottom")
cbar.set_label(r"Hours post fertilisation", labelpad=5)
ax.figure.subplots_adjust(top=0.925)

fig.savefig("height_pipette.pdf")
export(ax, "height_pipette.csv")

# velocities

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (hpf)")
ax.set_xlim([4, 8.5])
ax.set_ylabel(r"Velocities (µm/s)")

for v, marker, label in zip(
    (v0, vm, vp),
    ("s", "o", "d"),
    (
        r"$\dot{\ell}(0^+)$",
        r"$\dot{\ell}(t_{\mathrm{release}}^-)$",
        r"$\dot{\ell}(t_{\mathrm{release}}^+)$")):
    l, = ax.plot(t_pip, v, lw=0, marker=marker)
    ax.axhline(y=v.mean(), color=l.get_color(),
        label=r"%s=%.1f±%.2fµm/s" % (label, v.mean(), v.std()/2))
    ax.fill_between(ax.get_xlim(), v.mean() - v.std()/2, v.mean() + v.std()/2,
        color=l.get_color(), alpha=0.2)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper right"))

fig.savefig("velocities_pipette.pdf")
export(ax, "velocities_pipette.csv")

# times measured

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (hpf)")
ax.set_xlim([4, 8.5])
ax.set_ylabel(r"Times (s)")
ax.set_yscale("log")

for t, marker, label in zip(
    (taur, taurp, tauK),
    ("s", "o", "d"),
    (r"$\tau_r$", r"$\tau_r^{\prime}$", r"$\tau_{\mathrm{K}}$")):
    l, = ax.plot(t_pip, t, lw=0, marker=marker)
    ax.axhline(y=t.mean(), color=l.get_color(),
        label=r"%s=%i±%.1fs" % (label, t.mean(), t.std()/2))
    ax.fill_between(ax.get_xlim(), t.mean() - t.std()/2, t.mean() + t.std()/2,
        color=l.get_color(), alpha=0.2)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper right"))

fig.savefig("times_measured_pipette.pdf")
export(ax, "times_measured_pipette.csv")

# times computed

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (hpf)")
ax.set_xlim([4, 8.5])
ax.set_ylabel(r"Times (s)")
ax.set_yscale("log")

for t, marker, label in zip(
    (tau, taup),
    ("s", "o", "d"),
    (r"$\tau$", r"$\tau^{\prime}$")):
    l, = ax.plot(t_pip, t, lw=0, marker=marker)
    ax.axhline(y=t.mean(), color=l.get_color(),
        label=r"%s=%i±%.1fs" % (label, t.mean(), t.std()/2))
    ax.fill_between(ax.get_xlim(), t.mean() - t.std()/2, t.mean() + t.std()/2,
        color=l.get_color(), alpha=0.2)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper right"))

fig.savefig("times_computed_pipette.pdf")
export(ax, "times_computed_pipette.csv")

# elastic constant

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (hpf)")
ax.set_xlim([4, 8.5])
ax.set_ylabel(
    r"Elastic constant $k$ (10$^{\mathdefault{-3}}$kg/s$^{\mathdefault{2}}$)")
ax.set_ylim([3.5, 8.0])

for sk in (1e3*springk,):
    l, = ax.plot(t_pip, sk, lw=0, marker="s")
    ax.axhline(y=sk.mean(), color=l.get_color(), label=
        r"$k$=%.1f±%.2f$\times$10$^{\mathdefault{-3}}$kg/s$^{\mathdefault{2}}$"
        % (sk.mean(), sk.std()/2))
    ax.fill_between(ax.get_xlim(),
        sk.mean() - sk.std()/2, sk.mean() + sk.std()/2,
        color=l.get_color(), alpha=0.2)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper left"))

fig.savefig("springk_pipette.pdf")
export(ax, "springk_pipette.csv")

# residual tension

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (hpf)")
ax.set_xlim([4, 8.5])
ax.set_ylabel(r"Tension $T$ (µN)")

l, = ax.plot(t_pip, T, lw=0, marker="s")
ax.axhline(y=T.mean(), color=l.get_color(),
    label=r"$T$=%.1f±%.2fµN" % (T.mean(), T.std()/2))
ax.fill_between(ax.get_xlim(),
    T.mean() - T.std()/2, T.mean() + T.std()/2,
    color=l.get_color(), alpha=0.2)
ax.axhline(y=0, color="black", linestyle="--")

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper right"))

fig.savefig("tension_pipette.pdf")
export(ax, "tension_pipette.csv")

# substrate friction

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (hpf)")
ax.set_xlim([4, 8.5])
ax.set_ylabel(r"Friction coefficient $\zeta$ (kg/s)")

l, = ax.plot(t_pip, xi, lw=0, marker="s")
ax.axhline(y=xi.mean(), color=l.get_color(),
    label=r"$\zeta$=%.2f±%.3fkg/s" % (xi.mean(), xi.std()/2))
ax.fill_between(ax.get_xlim(),
    xi.mean() - xi.std()/2, xi.mean() + xi.std()/2,
    color=l.get_color(), alpha=0.2)

plt.sca(ax)
ax.add_artist(plt.legend(loc="upper right"))

fig.savefig("friction_pipette.pdf")
export(ax, "friction_pipette.csv")

# fit

fig, ax = plt.subplots()
ax.set_xlabel(r"Time (s)")
ax.set_ylabel(r"Height $\ell$ (µm)")

index = 8
time = sorted(height_pipette)[index]
assert time == t_pip[index]
tmax = height_pipette[time][np.argmax(height_pipette[time][:, 1]), 0]
hmax = np.max(height_pipette[time][:, 1])

def fit(t):
    if t < tmax:
        return (h0[index]
            + (v0[index] - vm[index])*taur[index]*(1 - np.exp(-t/taur[index]))
            + vm[index]*t)
    else:
        return (hmax
            + (height_pipette[time][-1, 1] - hmax)/(
                np.exp(-(height_pipette[time][-1, 0] - tmax)/taurp[index]) - 1)
                *(-1 + np.exp(-(t - tmax)/taurp[index])))

ax.plot(height_pipette[time][:, 0], height_pipette[time][:, 1],
    marker="s", markeredgecolor=None, label=r"Pipette experiment")
(lambda tt:
    ax.plot(tt, list(map(fit, tt)),
        label=r"Fit with $\tau_r$"))(
    (lambda t: t[t <= tmax])(height_pipette[time][:, 0]))
(lambda tt:
    ax.plot(tt, list(map(fit, tt)),
        label=r"Fit with $\tau_r^{\prime} > \tau_r$"))(
    (lambda t: t[t >= tmax])(height_pipette[time][:, 0]))

plt.sca(ax)
ax.add_artist(plt.legend(loc="lower right"))
fig.savefig("fit_pipette.pdf")
export(ax, "fit_pipette.csv")

# SHOW

exit()
plt.show()

