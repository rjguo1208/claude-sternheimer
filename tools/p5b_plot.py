#!/usr/bin/env python3
"""Plot real-space decay of H_W(R), M^W, Sigma^W, V~^W, and run the gauge test:
reconstruct H^W(R) from the SAME U(k) pipeline used for V~ and compare to hr.dat.
If my H^W(R) matches hr.dat (decays), the Wannier rotation/FT is correct."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994
HR="/anvil/projects/x-che190065/rjguo/qe-7.5/T-matrix/data/mos2_hr.dat"

fb=FortranFile("vtilde_block.dat","r")
N_A,nkstot,nk_use,nbndskip=fb.read_ints(np.int32)
omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2k=idx[:N_A]; a2band=idx[N_A:]
eact=fb.read_reals(np.float64)
M =fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
Sg=fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
V =fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
fb.close()
fw=FortranFile("wann_data.dat","r")
nk2,nbndep,nw=fw.read_ints(np.int32); ibndkept=fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nkstot),order="F")
U=fw.read_record(np.complex128).reshape((nbndep,nw,nkstot),order="F")   # [band,wann,k]
fw.close()
nb=nw; Nk=nkstot; Ub=np.transpose(U,(2,0,1))   # [k,band,wann]

# hr.dat
L=open(HR).read().split("\n"); nr=int(L[2]); li=3; deg=[]
while len(deg)<nr: deg+=[int(x) for x in L[li].split()]; li+=1
deg=np.array(deg[:nr],float); Rhr=np.zeros((nr,3),int); Hr=np.zeros((nr,nb,nb),complex)
d=[x for x in L[li:] if x.strip()]; k=0
for ir in range(nr):
    for n in range(nb):
        for m in range(nb):
            p=d[k].split(); k+=1; Rhr[ir]=(int(p[0]),int(p[1]),int(p[2]))
            Hr[ir,int(p[3])-1,int(p[4])-1]=float(p[5])+1j*float(p[6])
Hr_norm=np.array([np.linalg.norm(Hr[i]/deg[i]) for i in range(nr)])   # eV

# --- band Hamiltonian in slot order, per k ---
Hband=np.zeros((Nk,nb,nb),complex)
for ik in range(Nk):
    e=eact[a2k==(ik+1)]*RY        # slot order, eV
    Hband[ik]=np.diag(e)
HW_k=np.einsum("kbm,kbc,kcn->kmn", Ub.conj(), Hband, Ub)   # U^dag diag U -> H^W(k)
# FT to the hr R-vectors (try sign that matches hr.dat)
def HW_R(Rvec, sign):
    ph=np.exp(sign*2j*np.pi*(xk.T@Rvec))      # [k]
    return np.einsum("k,kmn->mn", ph, HW_k)/Nk
HWmine=np.array([np.linalg.norm(HW_R(Rhr[i],-1.0)) for i in range(nr)])
# pick sign by best match
err_m=max(np.abs(HW_R(Rhr[i],-1.0)-Hr[i]/deg[i]).max() for i in range(nr))
err_p=max(np.abs(HW_R(Rhr[i],+1.0)-Hr[i]/deg[i]).max() for i in range(nr))
sign=-1.0 if err_m<err_p else 1.0
HWmine=np.array([np.linalg.norm(HW_R(Rhr[i],sign)) for i in range(nr)])
gerr=min(err_m,err_p)*RY
print(f"GAUGE TEST: max|H^W_mine(R) - hr.dat H_W(R)| = {gerr:.3e} eV  "
      f"({'MATCH -> rotation correct' if gerr<1e-2 else 'MISMATCH -> rotation bug'})")

# --- Wannierize M, Sigma, V~ ; decay of ||O^W(R',R=0)|| vs |R'| ---
ms=np.arange(-6,6); Rg=np.array([(m,n,0) for m in ms for n in ms]); nR=len(Rg)
i0=[i for i in range(nR) if Rg[i,0]==0 and Rg[i,1]==0][0]
F=np.exp(2j*np.pi*(Rg@xk))    # [R,k]
def wann_decay(O):
    O4=O.reshape(Nk,nb,Nk,nb)
    OW=np.einsum("fsp,fsgt,gtw->fpgw",Ub.conj(),O4,Ub)
    OWR=np.einsum("Af,fpgw,Bg->ApBw",F,OW,F.conj())/Nk**2     # [R',wf,R,wi] Ry
    # coupling from origin cell: ||O^W(R', R=0)||
    col=np.linalg.norm(OWR[:,:,i0,:],axis=(1,2))             # vs R'
    diagR=np.linalg.norm(np.einsum("ApAw->Apw",OWR),axis=(1,2))  # on-site ||O^W(R,R)|| vs R
    return col, diagR
Mc,Md=wann_decay(M); Sc,Sd=wann_decay(Sg); Vc,Vd=wann_decay(V)
dist=np.sqrt(Rg[:,0]**2+Rg[:,1]**2)
dhr=np.sqrt(Rhr[:,0]**2+Rhr[:,1]**2)

# --- plots ---
fig,ax=plt.subplots(1,2,figsize=(12,4.6))
o=np.argsort(dhr)
ax[0].semilogy(dhr[o], Hr_norm[o]+1e-12,"o-",ms=3,label="hr.dat  H_W(R)")
ax[0].semilogy(dhr[o], HWmine[o]+1e-12,"x--",label="my pipeline H^W(R)")
ax[0].set_xlabel("|R| (prim. cells)"); ax[0].set_ylabel("||H_W(R)|| (eV)")
ax[0].set_title(f"Gauge test: H_W(R) decay  (max diff {gerr:.1e} eV)"); ax[0].legend(); ax[0].grid(alpha=.3)
o2=np.argsort(dist)
ax[1].semilogy(dist[o2], Mc[o2]+1e-12,".",label="Born  M^W(R',0)")
ax[1].semilogy(dist[o2], Sc[o2]+1e-12,".",label="rest  Σ^W(R',0)")
ax[1].semilogy(dist[o2], Vc[o2]+1e-12,".",label="V~^W(R',0)")
ax[1].set_xlabel("|R'| (prim. cells)"); ax[1].set_ylabel("||O^W(R',R=0)|| (Ry)")
ax[1].set_title("Downfolded potential: coupling from origin cell"); ax[1].legend(); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig("p5b_decay.png",dpi=110)
print("wrote p5b_decay.png")
print(f"  H_W(R):  |R|=0 {Hr_norm[dhr==0][0]:.3e} -> |R|>=5 {Hr_norm[dhr>=5].max():.3e} eV")
print(f"  M^W(R',0): R'=0 {Mc[i0]:.3e} -> |R'|>=5 {Mc[dist>=5].max():.3e} Ry")
print(f"  V~^W(R',0): R'=0 {Vc[i0]:.3e} -> |R'|>=5 {Vc[dist>=5].max():.3e} Ry")
