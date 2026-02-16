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
import os, sys, subprocess, time, signal, atexit, __main__

# DEFAULT KERATIN MODEL PARAMETERS

xi = 0.32       # kg s^{-1}
tauK = 93       # s
sigma = 0
tau = 500       # s
A0 = 600        # (micro m)^2
alpha = 8e3     # [ker] kg^{-1} (micro m) s^2
k = 5.4e-3      # kg s^{-2}
Gamma = k       # kg s^{-2}
beta = 5e-3
kth = 150
ron = 0
fpull = 0.57    # kg (micro m) s^{-2}
p0 = 3.72

kmax = 400      # maximum keratin value on colourmap

# COMMAND-LINE ARGUMENTS

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# keratin model parameters
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

# external forces parameters
parser.add_argument("-fpull", type=float, default=fpull)            # pulling force
parser.add_argument("-gamma", type=float, default=0)                # boundary tension
parser.add_argument("-ramp", action=argparse.BooleanOptionalAction) # activate force ramp

# random seed
parser.add_argument("-seed", type=int, default=np.random.randint(int(1e7)))

# input/output
parser.add_argument("input", nargs="?", default="disc.p")   # input keratin configuration
parser.add_argument("-dname", type=str, nargs="?")          # output directory

args = parser.parse_args()

# INITIALISATION

# helper functions
def get_perimeters(vm):
    """
    Returns perimeters of all cells.
    """
    return np.array(list(map(
        lambda i: vm.getVertexToNeighboursPerimeter(i),
        vm.getVertexIndicesByType("centre"))))
def get_areas(vm):
    """
    Returns areas of all cells.
    """
    return np.array(list(map(
        lambda i: vm.getVertexToNeighboursArea(i),
        vm.getVertexIndicesByType("centre"))))

# --- rescale distances in the system such that the mean force is low
# --- (dichotomic search) at the initial step
scalemin, scalemax = 0, 2
try: r = Read(args.input)                       # try to load input file as usual simulation output file
except AssertionError: pass
while scalemax - scalemin > 1e-6:

    try:
        vm = r[r.frames[-1]]                    # load last frame from simulation output file
    except NameError:
        with open(args.input, "rb") as dump:    # otherwise load input file as pickle file with unique configuration
            vm = pickle.load(dump)

    time0 = vm.time

    # --- set forces
    for _ in vm.vertexForces: vm.removeVertexForce(_)       # remove all prior forces
    for _ in vm.halfEdgeForces: vm.removeHalfEdgeForce(_)   # remove all prior forces
    vm.setOverdampedIntegrator(                             # overdamped dynamics
        args.xi)
    vm.addKeratinModel("keratin",                           # keratin model
        args.Gamma/A0, A0, args.tau,
        args.Gamma, args.p0,
        args.alpha, args.beta, args.kth,
        args.tauK, args.sigma, args.ron)

    # --- compute average force projected on the vector from the tissue centre to the cell centre
    scale = (scalemin + scalemax)/2
    vm.scale(scale*np.sqrt(A0/get_areas(vm).mean()))    # scale distances
    vm.nintegrate(1, 0)                                 # integrate with time step 0 to get forces
    meanForce = 0
    vertices = vm.getVertexIndicesByType("vertex")
    positions = (                                       # position of all cell centres
        lambda p: np.array(list(map(
            lambda i: p[i],
            vertices))))(
        vm.getPositions())
    posCM = positions.mean(axis=0)                      # position of tissue cnetre
    forces = vm.forces
    for i, index in enumerate(vertices):                # compute average projected force
        meanForce += np.dot(positions[i] - posCM, forces[index])/len(vertices)
    # breaking and update conditions for dichotomic search
    if args.fpull > 0:
        if np.abs(meanForce) < 1e-2*args.fpull: break
    elif args.gamma > 0:
        if np.abs(meanForce) < 1e-2*args.gamma: break
    else:
        break
    if meanForce > 0: scalemin = scale
    if meanForce < 0: scalemax = scale
    
# --- set random number generator seed
vm.setSeed(args.seed)

# --- set system size (since the system will be stretched we have to ensure the simulation box is large enough)
l = getMaxLengthBoundaries(vm)
assert len(l) == 1
l = list(l.values())[0]
vm.setSystemSize([10*l, 10*l])

# --- plotting options
if "DISPLAY" in os.environ:
    fig, ax = plot(vm, kmax=kmax)
    fig.axes[1].set_ylabel(r"$k_i$", labelpad=30)
    plt.ion()
    plt.show()
else:
    fig, ax = None, None

# --- integration parameters
dt = 2e-3*min(args.tau, args.tauK)                      # simulation time step
iterations = 2000                                       # total number of iterations between frames
frames_per_phase = int(np.ceil(40000/(dt*iterations)))  # total number of frames

# --- set pulling force
if args.fpull != 0:
    vm.addPressureForce("pull", 0, False)   # FPULL = PRESSURE x NUMBER OF BOUNDARY VERTICES x PERIMETER

# --- set boundary tension
if args.gamma != 0:
    vm.addBoundaryTension("boundary_tension", args.gamma)

# PERFORM SIMULATION

class Run:
    """
    Object to perform and save simulation.
    """

    def __init__(self, time0, dt, iterations, fpull, ramp):

        self.time0 = time0              # initial time of simulation
        self.dt = dt                    # integration time step
        self.iterations = iterations    # number of iterations between frames
        self.fpull = fpull              # reference pulling force
        self.ramp = ramp                # ramp up pulling force or keep it constant

        # output directory name with date and time
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

        # make movie of simulation if it is being plotted in real time
        if "DISPLAY" in os.environ: self.frames_dir = mkdtemp()
        self.count = 0  # frame count

    def integrate(self, vm):
        """
        Integrate vertex model one frame forward.
        """

        vm.nintegrate(self.iterations, dt=self.dt)

    def run(self, n, vm, fig, ax):
        """
        Run routine which integrates model, ramps up the pulling force, and
        plots the model.
        """

        for i in range(n):  # loop over the desired number of frames
            # save plotted frame
            if "DISPLAY" in os.environ:
                alpha, beta = (lambda param: (param["alpha"], param["beta"]))(
                    vm.vertexForces["keratin"].parameters)
                ax.set_title(
                    r"t$=%i$s, lengths in $\mu$m, $\alpha=%.2e$, $\beta=%.2e$"
                        % (vm.time - self.time0, alpha, beta),
                    pad=25, size=30)
                _update_canvas(fig)
                fig.savefig(
                    os.path.join(self.frames_dir, "%05d.png" % self.count))
            # save configuration to output file
            with open(self.fname, "ab") as dump: pickle.dump(vm, dump)
            self.count += 1
            # ramp up pulling force
            if self.fpull != 0:
                vm.removeVertexForce("pull")
                if self.ramp:
                    fpull = min(self.fpull, ((i + 1.)/(n/1.))*self.fpull)
                else:
                    fpull = self.fpull
                vm.addPressureForce("pull", fpull, False)   # FPULL = PRESSURE x NUMBER OF BOUNDARY VERTICES x PERIMETER
            # integrate vertex model
            self.integrate(vm)
            # update plotted frame
            if "DISPLAY" in os.environ:
                plot(vm, fig=fig, ax=ax, time0=self.time0, update=False)

# make movie from frames at the end of the script
def _exit_handler(*_args, **_kwargs):
    if "DISPLAY" in os.environ:
        subprocess.call(                            # compile frames into movie
            [movie_sh_fname, "-d", obj.frames_dir, "-p", sys.executable, "-y"],
            cwd=obj.dname)
    if hasattr(__main__, "__file__"): os._exit(0)   # not a python console: exit
# use exit handlers in case the program is terminated before completion
signal.signal(signal.SIGINT, _exit_handler)
signal.signal(signal.SIGTERM, _exit_handler)
atexit.register(_exit_handler)

# run simulation
obj = Run(time0, dt, iterations, args.fpull, args.ramp)
count = obj.run(frames_per_phase, vm, fig, ax)

