import os, sys
import numpy as np
import matplotlib.pyplot as plt

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

def save(name, ax):
    ax.get_figure().savefig("%s.pdf" % name)
    export(ax, "%s.csv" % name)

# OUTPUT DIRECTORY

from shutil import rmtree
output_dname = os.path.join(data_dname, "output")
try: rmtree(output_dname)
except FileNotFoundError: pass
os.mkdir(output_dname)
os.chdir(output_dname)

# PLOTS

plt.rcParams["figure.figsize"] = 7, 5.3
plt.rcParams["figure.dpi"] = 100
plt.rcParams["font.size"] = 24
plt.rcParams["font.weight"] = "normal"
plt.rcParams["axes.labelsize"] = 24
plt.rcParams["xtick.labelsize"] = 24
plt.rcParams["ytick.labelsize"] = 24
plt.rcParams["legend.fontsize"] = 16
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["axes.prop_cycle"] = 'cycler(color='\
    '["#83BB03", "#ff7f0e", "#BB9703", "#D51B66", "#5C2352", "#0173B2"])'
plt.rcParams["legend.fancybox"] = False
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

# include Keratin feedback
feedback=True
# contrast knockout model
# feedback=False

# Epiboly fixed constants from pipette matching
zeta = 0.32
tauK = 93
k0 = 0.054

# The relaxation time is in the end a fit factor
tau0 = 500 

Kth = 150 # in epiboly units already
conv=3600  # seconds to hours
conv2 = 60 # seconds to minutes
# From Yann
tmax = 40176
h0 = 10 # initial tissue height

# inital values
l0ini = np.sqrt(600) # square root of initial cell area
V0ini = l0ini**2*h0 # initial cell volume - constant
lini = l0ini 
Kini = Kth

# coupling constants
alpha = 120000
if feedback:
	beta = 0.005
else:
	beta = 0.0

def kval(K):
    return k0*(1+beta*max(0,K-Kth))

def tauval(K):
    return tau0*(1+beta*max(0,K-Kth))

# proper 3d definition of a pressure
def pressure3(l,l0,K):
    return -kval(K)*(l-l0)*l/V0ini

def ldot(l,l0,K,Fpull):
    return -kval(K)*(l-l0)+Fpull

def l0dot(l,l0,K):
    return (l-l0)/tauval(K)

def Kdot(l,l0,K):
    return alpha*max(0,-pressure3(l,l0,K))/tauK-K/tauK


# nullcline equations
def ucline(DeltaK,Fpull):
    us=np.zeros((len(DeltaK),))
    isneg=np.where(DeltaK<0)
    ispos=np.where(DeltaK>=0)
    # Derived the second piece, equivalent to setting DeltaK to 0
    us[ispos]=Fpull*tau0/zeta*(1+beta*DeltaK[ispos])/(1+k0*tau0/zeta*(1+beta*DeltaK[ispos])**2)
    us[isneg]=Fpull*tau0/zeta/(1+k0*tau0/zeta)
    return us

def Kcline(DeltaK,a):
    us=np.zeros((len(DeltaK),))
    isneg=np.where(DeltaK<0)
    ispos=np.where(DeltaK>=0)
    us[ispos]=(DeltaK[ispos]+Kth)/(alpha*a*k0/V0ini*(1+beta*DeltaK[ispos]))
    us[isneg]=(DeltaK[isneg]+Kth)/(alpha*a*k0/V0ini)
    return us

Fevals = [0.1425, 0.285,0.57,1.14]
dt=1
 # Total runtime: 9 hours
N = int(1.0/dt*conv*9)
timeval=np.zeros((N,))
Kout = np.zeros((N,len(Fevals)))
ldotout = np.zeros((N,len(Fevals)))
uout= np.zeros((N,len(Fevals)))
lout= np.zeros((N,len(Fevals)))
n=0
for Fe in Fevals:
    lval=np.zeros((N,))
    ldotval=np.zeros((N,)) # for plotting edge speed
    l0val=np.zeros((N,))
    Kval=np.zeros((N,))
    pressval=np.zeros((N,))
    Fcurr=np.zeros((N,))

    lval[0] = lini
    l0val[0] = l0ini
    Kval[0] = 0 # actual 0

    # Estupido Euler to start with
    for k in range(1,N):
        if k%1000==0:
            print(k)
        timeval[k] = timeval[k-1]+dt
        if timeval[k]<tmax:
            Fcurr[k] = Fe*timeval[k]/tmax
        else:
            Fcurr[k] = Fe
        ldotval[k] = ldot(lval[k-1],l0val[k-1],Kval[k-1],Fcurr[k])
        ldotout[k,n] = ldotval[k]
        lval[k] = lval[k-1]+ldotval[k]*dt
        l0val[k] = l0val[k-1]+l0dot(lval[k-1],l0val[k-1],Kval[k-1])*dt
        Kval[k] = Kval[k-1]+Kdot(lval[k-1],l0val[k-1],Kval[k-1])*dt
        Kout[k,n] = Kval[k]
        uout[k,n] = lval[k]-l0val[k]
        lout[k,n] = lval[k]

    n+=1

cols = ['g','r','m','b']

plt.figure()
for n in range(len(Fevals)):
    plt.plot(timeval/conv,Kout[:,n],color=cols[n],label=r'$F_{YSL} = $' + str(Fevals[n])+'μN')
plt.xlabel('Time (h)')
plt.ylabel('Keratin intensity')
plt.legend()
if feedback: save("figSM9B", plt.gca())

plt.figure()
for n in range(len(Fevals)):
    plt.plot(timeval/conv,ldotout[:,n]*conv2,color=cols[n],label=r'$F_{YSL} = $' + str(Fevals[n])+'μN')
plt.xlabel('Time (h)')
plt.ylabel(r"$ \dot{\ell}$ (μm/min)")
plt.legend(loc="upper right" if feedback else "upper left")
if feedback: save("figSM9C", plt.gca())
else: save("figSM9D", plt.gca())

# Nullclines
plt.figure()
DKpts = np.linspace(-150,1200,500)
for n in range(len(Fevals)):
    Fe=Fevals[n]
    uval = uout[:,n]
    lval = lout[:,n]
    DeltaKval = Kout[:,n]-Kth
    plt.plot(DeltaKval[0:N:400],uval[0:N:400],marker='.',markersize=10,linestyle='none',color=cols[n],label=r'$F_{YSL} = $' + str(Fevals[n])+'μN')
    Fend = Fe*timeval[N-1]/tmax
    plt.plot(DKpts,ucline(DKpts,Fend),linestyle='--',lw=2,color=cols[n])
    plt.plot(DKpts,Kcline(DKpts,lval[N-1]),linestyle='-',lw=2,color=cols[n])
plt.xlabel(r'$\Delta K$')
plt.ylabel('u (μm)')
plt.ylim(-0.3,4.5)
plt.legend()
#plt.title('Force ' + str(Fe) + ' alpha ' + str(alpha) + ' beta ' + str(beta))
if feedback: save("figSM9A", plt.gca())

plt.show()

 
