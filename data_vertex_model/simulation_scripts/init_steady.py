from cells.plot import WindowClosedException
from cells.keratin import plot_keratin
from cells.run import run
from cells.init import init_vm

import pickle

args, _ = init_vm()
with open("init_keratin.p", "rb") as dump:
    vm = pickle.load(dump)
parameters = vm.vertexForces["keratin"].parameters
dt = 2e-3*min(parameters["tau"], parameters["taur"])

args.dt = dt
args.save = True
def plot(*args, **kwargs): return plot_keratin(*args, **kwargs, kmax=400)
run(args, vm, plot_function=plot)

