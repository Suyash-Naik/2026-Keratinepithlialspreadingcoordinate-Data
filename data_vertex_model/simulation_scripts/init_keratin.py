from cells.plot import _update_canvas
from cells.keratin import plot_keratin as plot

import pickle
import numpy as np
import matplotlib.pyplot as plt

# READ

# load initial disordered configuration
with open("init_disordered.p", "rb") as dump:
    vm = pickle.load(dump)

# remove all forces
for _ in vm.vertexForces: vm.removeVertexForce(_)
for _ in vm.halfEdgeForces: vm.removeHalfEdgeForce(_)

# compute areas of cells in input configuration
areas = np.array(list(map(
    lambda i: vm.getVertexToNeighboursArea(i),
    vm.getVertexIndicesByType("centre"))))

# MODIFY

# keratin default parameters
xi = 0.32
tauK = 93
sigma = 0
tau = 500
A0 = 600
alpha = 4e3
k = 5.4e-3
Gamma = k
beta = 5e-3
kth = 150
ron = 0
p0 = 3.72

# initialisation parameters
meanK = 100     # target mean keratin value
meanA = 1750    # target mean area
meanRatio = (-1 + np.sqrt((4*meanK)/(alpha*Gamma)))/8

# scale system to change mean area
vm.scale(np.sqrt(meanA/areas.mean()))

# add keratin
beta, tau = 0, np.inf           # do not close loop and do not relax target areas
vm.setOverdampedIntegrator(     # add overdamped integrator with drag coefficient xi
    xi)
vm.addKeratinModel("keratin",   # load keratin model
    Gamma/A0, A0, tau, Gamma, p0, alpha, beta, kth, tauK, sigma, ron)

# change target areas
vm.nintegrate(1, 0)
vm.vertexForces["keratin"].targetArea = dict(map(
    lambda i: (i, vm.vertexForces["keratin"].area[i]/(1 + meanRatio)),
    vm.getVertexIndicesByType("centre")))

# change keratin
vm.vertexForces["keratin"].keratin = dict(map(
    lambda i: (i, meanK),
    vm.getVertexIndicesByType("centre")))

# SAVE

# output configuration with keratin model and desired areas
with open("init_keratin.p", "wb") as dump:
    pickle.dump(vm, dump)

# PLOT

fig, ax = plot(vm, kmax=400)
plt.show()

