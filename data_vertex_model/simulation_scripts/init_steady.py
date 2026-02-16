from cells.plot import WindowClosedException
from cells.keratin import plot_keratin
from cells.run import run
from cells.init import init_vm

import pickle

args, _ = init_vm()                                     # load command-line arguments
with open("init_keratin.p", "rb") as dump:              # load configuration with keratin model and desired areas
    vm = pickle.load(dump)
parameters = vm.vertexForces["keratin"].parameters
dt = 2e-3*min(parameters["tau"], parameters["taur"])    # integration time step

args.dt = dt
args.save = True                                                            # activate the saving of configurations during real-time simulation
def plot(*args, **kwargs): return plot_keratin(*args, **kwargs, kmax=400)   # custom plotting function
run(args, vm, plot_function=plot)                                           # launch real-time simulation

