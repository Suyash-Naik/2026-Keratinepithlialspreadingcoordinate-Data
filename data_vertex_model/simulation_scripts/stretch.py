from cells.bind import VertexModel, getMaxLengthBoundaries
from cells.read import Read
from cells.init import movie_sh_fname
from cells.keratin import plot_keratin as plot
from cells.plot import _update_canvas

import numpy as np
import argparse
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
if mpl.get_backend() == "agg": mpl.use("ps")
from tempfile import mkdtemp
import os, sys, subprocess, time

# --- pipette
xi = 0.32       # kg s^{-1}
tauK = 93       # s
sigma = 0
tau = 500       # s
A0 = 600        # (micro m)^2
alpha = 8e3     # [ker] kg^{-1} (micro m) s^2
k = 5.4e-3      # kg s^{-2}
Gamma = k       # kg s^{-2}
# K = k/A0        # kg (micro m)^{-2} s^{-2}
beta = 5e-3
kth = 150
ron = 0
fpull = 0.57    # kg (micro m) s^{-2}
p0 = 3.72

kmax = 400

# show tensions and keratin even below threshold

# parameters

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-xi", type=float, default=xi)
parser.add_argument("-tau", type=float, default=tau)
parser.add_argument("-tauK", type=float, default=tauK)
parser.add_argument("-p0", type=float, default=p0)
parser.add_argument("-alpha", type=float, default=alpha)
parser.add_argument("-Gamma", type=float, default=Gamma)
parser.add_argument("-beta", type=float, default=beta)
parser.add_argument("-kth", type=float, default=kth)
parser.add_argument("-sigma", type=float, default=sigma)
parser.add_argument("-ron", type=float, default=ron)

parser.add_argument("-fpull", type=float, default=fpull)
parser.add_argument("-gamma", type=float, default=0)
parser.add_argument("-ramp", action=argparse.BooleanOptionalAction)

parser.add_argument("-seed", type=int, default=np.random.randint(int(1e7)))

parser.add_argument("input", nargs="?", default="out.p")

parser.add_argument("-dname", type=str, nargs="?")

args = parser.parse_args()

# initialisation

def get_perimeters(vm):
    return np.array(list(map(
        lambda i: vm.getVertexToNeighboursPerimeter(i),
        vm.getVertexIndicesByType("centre"))))
def get_areas(vm):
    return np.array(list(map(
        lambda i: vm.getVertexToNeighboursArea(i),
        vm.getVertexIndicesByType("centre"))))

# --- rescale distances such that mean force is low (dichotomic search)
scalemin, scalemax = 0, 2
try: r = Read(args.input)
except AssertionError: pass
while scalemax - scalemin > 1e-6:

    try:
        vm = r[r.frames[-1]]
    except NameError:
        with open(args.input, "rb") as dump:
            vm = pickle.load(dump)

    time0 = vm.time

    # --- set forces
    for _ in vm.vertexForces: vm.removeVertexForce(_)
    for _ in vm.halfEdgeForces: vm.removeHalfEdgeForce(_)
    vm.setOverdampedIntegrator(
        args.xi)
    vm.addKeratinModel("keratin",
        args.Gamma/A0, A0, args.tau,
        args.Gamma, args.p0,
        args.alpha, args.beta, args.kth,
        args.tauK, args.sigma, args.ron)

    # --- compute forces
    scale = (scalemin + scalemax)/2
    vm.scale(scale*np.sqrt(A0/get_areas(vm).mean()))    # scale distances
    vm.nintegrate(1, 0)                                 # integrate with time step 0 to get forces
    meanForce = 0
    vertices = vm.getVertexIndicesByType("vertex")
    positions = (
        lambda p: np.array(list(map(
            lambda i: p[i],
            vertices))))(
        vm.getPositions())
    posCM = positions.mean(axis=0)
    forces = vm.forces
    for i, index in enumerate(vertices):
        meanForce += np.dot(positions[i] - posCM, forces[index])/len(vertices)
    if args.fpull > 0:
        if np.abs(meanForce) < 1e-2*args.fpull: break
    elif args.gamma > 0:
        if np.abs(meanForce) < 1e-2*args.gamma: break
    else:
        break
    if meanForce > 0: scalemin = scale
    if meanForce < 0: scalemax = scale
    
# --- set seed
vm.setSeed(args.seed)

# --- set system size
l = getMaxLengthBoundaries(vm)
assert len(l) == 1
l = list(l.values())[0]
vm.setSystemSize([10*l, 10*l])

if "DISPLAY" in os.environ:
    fig, ax = plot(vm, kmax=kmax)
    fig.axes[1].set_ylabel(r"$k_i$", labelpad=30)
    #fig.axes[2].set_ylabel(r"scaled $t_i$", labelpad=30)
    plt.ion()
    plt.show()
else:
    fig, ax = None, None

dt = 2e-3*min(args.tau, args.tauK)
iterations = 2000
# frames_per_phase = int(np.ceil(18000/(dt*iterations)))
frames_per_phase = int(np.ceil(40000/(dt*iterations)))

# dt = 1e-3
# iterations = 2000
# frames_per_phase = int(np.ceil(100/(dt*iterations)))

# --- set pulling force
if args.fpull != 0:
    vm.addPressureForce("pull", 0, False)   # FPULL = PRESSURE x NUMBER OF BOUNDARY VERTICES x PERIMETER

# --- set boundary tension
if args.gamma != 0:
    vm.addBoundaryTension("boundary_tension", args.gamma)

class Run:

    def __init__(self, time0, dt, iterations, fpull, ramp):

        self.time0 = time0
        self.dt = dt
        self.iterations = iterations
        self.fpull = fpull
        self.ramp = ramp

        if not(args.dname is None): dname = args.dname
        else: dname = time.strftime("%Y-%m-%d_%H-%M-%S", time.gmtime())
        count = 0
        while True:
            try:
                self.dname = "%s.%i" % (dname, count)
                os.mkdir(self.dname)
                break
            except FileExistsError:
                count += 1
        self.fname = os.path.join(self.dname, "out.p")
        with open(self.fname, "wb") as dump: pass

        if "DISPLAY" in os.environ: self.frames_dir = mkdtemp()
        self.count = 0

    def integrate(self, vm):
        vm.nintegrate(self.iterations, dt=self.dt)
        param = vm.vertexForces["keratin"].parameters
#         #####
#         A0 = param["A0"]
#         P0 = np.sqrt(param["A0"])*param["p0"] - param["T"]/(2*param["Gamma"])
#         print("areas [A0=%s]" % A0)
#         print(get_areas(vm)/A0)
#         print("perimeters [P0=%s]" % P0)
#         print(get_perimeters(vm)/P0)
#         #####

    def run(self, n, vm, fig, ax):

        for i in range(n):
            if "DISPLAY" in os.environ:
                alpha, beta = (lambda param: (param["alpha"], param["beta"]))(
                    vm.vertexForces["keratin"].parameters)
                ax.set_title(
                    r"t$=%i$s, lengths in $\mu$m, $\alpha=%.2e$, $\beta=%.2e$"
                        % (vm.time - self.time0, alpha, beta),
                    pad=25, size=30)
#                 ax.set_xlim([100, 400])
#                 ax.set_xticks([100, 200, 300, 400])
#                 ax.set_xticklabels([r"$0$", r"$100$", r"$200$", r"$300$"])
#                 ax.set_ylim([100, 400])
#                 ax.set_yticks([100, 200, 300, 400])
#                 ax.set_yticklabels([r"$300$", r"$200$", r"$100$", r"$0$"])
                _update_canvas(fig)
                fig.savefig(
                    os.path.join(self.frames_dir, "%05d.png" % self.count))
            with open(self.fname, "ab") as dump: pickle.dump(vm, dump)
            self.count += 1
            if self.fpull != 0:
                vm.removeVertexForce("pull")
                if self.ramp:
                    fpull = min(self.fpull, ((i + 1.)/(n/1.))*self.fpull)
                else:
                    fpull = self.fpull
                vm.addPressureForce("pull", fpull, False)   # FPULL = PRESSURE x NUMBER OF BOUNDARY VERTICES x PERIMETER
            self.integrate(vm)
            if "DISPLAY" in os.environ:
                plot(vm, fig=fig, ax=ax, time0=self.time0, update=False)
#             pressure = list(vm.vertexForces["keratin"].pressure.values())
#             tension = list(vm.vertexForces["keratin"].tension.values())
#             keratin = list(vm.vertexForces["keratin"].keratin.values())
#             print("min(pressure) = %s, max(pressure) = %s"
#                 % (min(pressure), max(pressure)))
#             print("min(tension) = %s, max(tension) = %s"
#                 % (min(tension), max(tension)))
#             print("min(keratin) = %s, max(keratin) = %s"
#                 % (min(keratin), max(keratin)))
#             print()
#             velocities = np.array(list(vm.velocities.values()))
#             print("max(v*dt) = %s"
#                 % (np.abs(velocities).max()*self.dt))
#             print()

# save

obj = Run(time0, dt, iterations, args.fpull, args.ramp)

# # initialise
# 
# obj.integrate(vm)

# pull

count = obj.run(frames_per_phase, vm, fig, ax)
# try: count = obj.run(frames_per_phase, vm, fig, ax)
# except: pass

# release

# vm.removeVertexForce("pull")
# count = obj.run(frames_per_phase, vm, fig, ax)

# movie

if "DISPLAY" in os.environ:
    subprocess.call(
        [movie_sh_fname, "-d", obj.frames_dir, "-p", sys.executable, "-y"],
        cwd=obj.dname)

