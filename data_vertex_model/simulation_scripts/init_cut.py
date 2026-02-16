from cells.plot import _update_canvas
from cells.keratin import plot_keratin as plot

import pickle
import matplotlib.pyplot as plt

# READ

# load disordered configuration with desired keratin levels
with open("init_steady.p", "rb") as dump:
    vm = pickle.load(dump)

# MODIFY

# keratin model parameters
xi = 0.32       # kg s^{-1}
tauK = 93       # s
sigma = 0
tau = 500       # s
A0 = 600        # (micro m)^2
alpha = 4e3     # [ker] kg^{-1} (micro m) s^2
k = 5.4e-3      # kg s^{-2}
Gamma = k       # kg s^{-2}
beta = 5e-3
kth = 150
ron = 0
p0 = 3.72
keratin = vm.vertexForces["keratin"].keratin        # conserve initial keratin levels
targetArea = vm.vertexForces["keratin"].targetArea  # conserve initial target areas
vm.removeVertexForce("keratin")
vm.addKeratinModel("keratin",                       # set keratin model
    Gamma/A0, A0, tau,
    Gamma, p0,
    alpha, beta, kth,
    tau, sigma, ron)
vm.vertexForces["keratin"].keratin = keratin        # load initial keratin levels
vm.vertexForces["keratin"].targetArea = targetArea  # load initial target areas

# create hole by...
# ... first merging cells...
vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(639, 682))[1])
vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(682, 808))[1])
vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(808, 894))[1])
vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(894, 934))[1])
vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(934, 1020))[1])
# vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(1020, 1023))[1])
# vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(1023, 980))[1])
# vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(980, 1109))[1])
# ... then changing the merged cell to a boundary vertex so as to delete it and create an exterior
vm.changeToBoundary(1020)

# add boundary tension
gamma = 2
vm.addBoundaryTension("boundary_tension",
    gamma)

# SAVE

vm.nintegrate(1, 0)
with open("init_hole.p", "wb") as dump: # output configuration with keratin model, boundary tension, and hole
    pickle.dump(vm, dump)

# PLOT

fig, ax = plot(vm, kmax=400)
plt.show()

