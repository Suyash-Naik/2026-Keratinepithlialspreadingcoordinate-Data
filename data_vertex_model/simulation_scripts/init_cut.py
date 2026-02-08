from cells.plot import _update_canvas
from cells.keratin import plot_keratin as plot

import pickle
import matplotlib.pyplot as plt

# READ

with open("init_steady.p", "rb") as dump:
    vm = pickle.load(dump)

# MODIFY

# add keratin
# --- pipette
xi = 0.32       # kg s^{-1}
tauK = 93       # s
sigma = 0
tau = 500       # s
A0 = 600        # (micro m)^2
alpha = 4e3     # [ker] kg^{-1} (micro m) s^2
k = 5.4e-3      # kg s^{-2}
Gamma = k       # kg s^{-2}
# K = k/A0        # kg (micro m)^{-2} s^{-2}
beta = 5e-3
# beta = 0
kth = 150
ron = 0
p0 = 3.72
keratin = vm.vertexForces["keratin"].keratin
targetArea = vm.vertexForces["keratin"].targetArea
vm.removeVertexForce("keratin")
vm.addKeratinModel("keratin",
    Gamma/A0, A0, tau,
    Gamma, p0,
    alpha, beta, kth,
    tau, sigma, ron)
vm.vertexForces["keratin"].keratin = keratin
vm.vertexForces["keratin"].targetArea = targetArea

# create hole
vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(639, 682))[1])
vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(682, 808))[1])
vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(808, 894))[1])
vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(894, 934))[1])
vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(934, 1020))[1])
# vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(1020, 1023))[1])
# vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(1023, 980))[1])
# vm.mergeVertices(vm.swapEdge(vm.getHalfEdgeBetweenIndex(980, 1109))[1])
vm.changeToBoundary(1020)

# add boundary tension
gamma = 2
vm.addBoundaryTension("boundary_tension",
    gamma)

# SAVE

vm.nintegrate(1, 0)
with open("init_hole.p", "wb") as dump:
    pickle.dump(vm, dump)

# PLOT

fig, ax = plot(vm, kmax=400)

# (
#     lambda vertices: list(map(
#         lambda i: ax.text(*vertices[i].position, i),
#         vm.getVertexIndicesByType("centre"))))(
#     vm.vertices)
# _update_canvas(fig)

plt.show()

