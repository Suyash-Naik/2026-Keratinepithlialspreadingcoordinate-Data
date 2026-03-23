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
#plt.rcParams["figure.figsize"] = 9, 6.5
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

# Epiboly fixed constants from pipette matching
zetaAV = 0.32 # \pm 0.066
tauKAV = 93 # \pm 19.5

# Conversion factor back to pipette keratin units (~ number of slices contributing to max projection)
kconv=15

k0AV = 0.054 # \pm 0.0054
tau0AV = 24 # \pm 3.9 (result of inversion)
KthAV = 150 # \pm 50 - mostly experimental variability
tpull=200
h0 = 10 # initial tissue height

# inital values
l0ini = np.sqrt(600) # square root of initial cell area
V0ini = l0ini**2*h0 # initial cell volume - constant
lini = l0ini 

# coupling constants
alphaAV = 120000 # the conversion to 3d volume # \pm  at least 2000
betaAV = 0.005 # and open: \pm 0.002 or something



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

Fe = 0.57 # fixed by pipette machine
TAV=0.1 # + 0.3 - 0.1
dt=0.01
 # Total runtime: 500s
N = int(1.0/dt*500)

n=0

Nsample = 40
zetaval=[]
tauKval=[]
k0val=[]
tau0val=[]
Kthval=[]
Tval=[]
alphaval=[]
betaval=[]
for s in range(Nsample):
    zetaval.append(zetaAV + np.random.uniform(-0.066,0.066))
    tauKval.append(tauKAV+ np.random.uniform(-19.5,19.5))
    k0val.append(k0AV+ np.random.uniform(-0.0054,0.0054))
    tau0val.append(tau0AV+ np.random.uniform(-3.9,3.9))
    Kthval.append(KthAV+ np.random.uniform(-50,50))
    Tval.append(TAV+ np.random.uniform(-0.1,0.1))
    alphaval.append(alphaAV+ np.random.uniform(-20000,20000))
    betaval.append(betaAV+ np.random.uniform(-0.0002,0.0002))
# Last one is the mean
zetaval.append(zetaAV)
tauKval.append(tauKAV)
k0val.append(k0AV)
tau0val.append(tau0AV)
Kthval.append(KthAV)
Tval.append(TAV)
alphaval.append(alphaAV)
betaval.append(betaAV)

fig_len = plt.figure()
fig_K = plt.figure()
fig_tau = plt.figure()
fig_kstiff = plt.figure()
# Combined spectrum of extreme values of parameters
for s in range(Nsample+1):
    print(s)
    zeta=zetaval[s]
    tauK=tauKval[s]
    k0=k0val[s]
    tau0=tau0val[s]
    Kth=Kthval[s]
    T=Tval[s]
    alpha=alphaval[s]
    beta=betaval[s]

    def kval(K):
        if hasattr(K, "__len__"):
            outkval=k0*np.ones((len(K),))
            islargeK=np.where(K>Kth)[0]
            outkval[islargeK] = k0*(1+beta*(K[islargeK]-Kth))
            return outkval
        else:
            return k0*(1+beta*max(0,K-Kth))

    def tauval(K):
        if hasattr(K, "__len__"):
            outtauval=tau0*np.ones((len(K),))
            islargeK=np.where(K>Kth)[0]
            outtauval[islargeK] = tau0*(1+beta*(K[islargeK]-Kth))
            return outtauval
        else:
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


    lval=np.zeros((N,))
    ldotval=np.zeros((N,)) # for plotting edge speed
    l0val=np.zeros((N,))
    Kval=np.zeros((N,))
    pressval=np.zeros((N,))
    Fcurr=np.zeros((N,))
    timeval=np.zeros((N,))

    lval[0] = lini+0.1*lini*np.random.uniform(-1,1)
    l0val[0] = l0ini
    Kval[0] = kconv*5 + kconv*np.random.uniform(-1,1) # eyeballing graphs

    # Estupido Euler to start with
    for k in range(1,N):
        #if k%1000==0:
        #    print(k)
        timeval[k] = timeval[k-1]+dt
        if timeval[k]<tpull:
            Fcurr[k] = Fe-T
        else:
            Fcurr[k] = -T
        ldotval[k] = ldot(lval[k-1],l0val[k-1],Kval[k-1],Fcurr[k])
        lval[k] = lval[k-1]+ldotval[k]*dt
        l0val[k] = l0val[k-1]+l0dot(lval[k-1],l0val[k-1],Kval[k-1])*dt
        Kval[k] = Kval[k-1]+Kdot(lval[k-1],l0val[k-1],Kval[k-1])*dt

    if s<Nsample:
        fig_K.gca().plot(timeval,Kval/kconv,'-g',alpha=0.2,lw=1)
        fig_len.gca().plot(timeval,lval,color='r',alpha=0.2,lw=1)
        fig_tau.gca().plot(timeval,tauval(Kval),color='b',alpha=0.2,lw=1)
        fig_kstiff.gca().plot(timeval,kval(Kval),color='m',alpha=0.2,lw=1)
    else:
        print('fitted mean value')
        fig_K.gca().plot(timeval,Kval/kconv,'-g')
        fig_len.gca().plot(timeval,lval,color='r')
        fig_tau.gca().plot(timeval,tauval(Kval),color='b')
        fig_kstiff.gca().plot(timeval,kval(Kval),color='m')



fig_K.gca().set_xlabel('time (s)')
fig_K.gca().set_ylabel('Keratin intensity')
fig_K.gca().set_ylim(0,35)
save("figSM9F", fig_K.gca())

fig_len.gca().set_xlabel('time (s)')
fig_len.gca().set_ylabel('junction length')
fig_len.gca().set_ylim(0,80)
save("figSM9E", fig_len.gca())

fig_tau.gca().set_xlabel('time (s)')
fig_tau.gca().set_ylabel(r"$ \tau(K)$ ")
fig_tau.gca().set_ylim(0,80)
save("figSM9H", fig_tau.gca())

fig_kstiff.gca().set_xlabel('time (s)')
fig_kstiff.gca().set_ylabel(r"$ k(K)$ ")
fig_kstiff.gca().set_ylim(0,0.2)
save("figSM9G", fig_kstiff.gca())


plt.show()

 
