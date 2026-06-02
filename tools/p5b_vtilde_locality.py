#!/usr/bin/env python3
"""Wannierization of the downfolded potential Vtilde = M + Sigma, and its locality.
Left : both-index ||Vtilde^W(R',R)||_F by shell max(|R'|,|R|), OLD gauge vs NEW (150-band) gauge.
Right: electron-index ||Vtilde^W(R_e;q)|| decay (new gauge) for a few q -> Vtilde^W is localized."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994
fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
fb.read_ints(np.int32); fb.read_reals(np.float64)
M =fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
Sg=fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
V =fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
nb=N_A//nk; n1=int(round(nk**0.5)); V4=V.reshape(nk,nb,nk,nb)

def load(fn):
    fw=FortranFile(fn,"r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
    xk=fw.read_reals(np.float64).reshape((3,nk),order="F")
    U=fw.read_record(np.complex128).reshape((nb,nb,nk),order="F"); fw.close()
    return np.transpose(U,(2,0,1)), xk

ms=np.arange(-n1//2,n1//2); Rg=np.array([(m,n,0) for m in ms for n in ms]); nR=len(Rg)
sh=np.maximum(np.abs(Rg[:,0]),np.abs(Rg[:,1]))

def vtilde_W(Ub,xk):                              # both-index Vtilde^W(R',R)
    VW=np.einsum("fsp,fsgt,gtw->fpgw",Ub.conj(),V4,Ub)
    F=np.exp(2j*np.pi*(Rg@xk))
    return np.einsum("Af,fpgw,Bg->ApBw",F,VW,F.conj())/nk**2, VW

def shells(VWR):
    bl=np.linalg.norm(VWR,axis=(1,3))
    return np.array([bl[(sh[:,None]==d)|(sh[None,:]==d)].max() for d in range(7)])

Ub_o,xk_o=load("wann_data_17gauge.dat"); VWR_o,_=vtilde_W(Ub_o,xk_o)
Ub_n,xk_n=load("wann_data.dat");        VWR_n,VW_n=vtilde_W(Ub_n,xk_n)
sh_o, sh_n = shells(VWR_o), shells(VWR_n)
print("shell  oldVtilde^W   newVtilde^W (Ry)")
for d in range(7): print(f"  {d}    {sh_o[d]:.4e}   {sh_n[d]:.4e}")

# electron-index decay of Vtilde (new gauge): FT VW over electron k at fixed q
ij=np.rint(np.array([xk_n[0]*n1,xk_n[1]*n1]).T).astype(int)%n1
idx={(ij[k,0],ij[k,1]):k for k in range(nk)}
kadd=lambda k,q: idx[((ij[k,0]+ij[q,0])%n1,(ij[k,1]+ij[q,1])%n1)]
def edecay(qi):
    edms=np.array([VW_n[kadd(k,qi),:,k,:] for k in range(nk)])
    F=np.exp(-2j*np.pi*(Rg@xk_n)); MR=np.einsum("Rk,kpw->Rpw",F,edms)/nk
    return np.sqrt(Rg[:,0]**2+Rg[:,1]**2), np.linalg.norm(MR,axis=(1,2))

fig,ax=plt.subplots(1,2,figsize=(12,4.6))
ax[0].semilogy(range(7), sh_o+1e-12,"s--",c="#d62728",label="old filukk (gauge mismatch)")
ax[0].semilogy(range(7), sh_n+1e-12,"o-", c="#1f77b4",label="filukk_150 (consistent gauge)")
ax[0].set_xlabel("shell  max(|R'|,|R|)  (prim. cells)"); ax[0].set_ylabel(r"max $\|\tilde V^W(R',R)\|_F$ (Ry)")
ax[0].set_title(r"Both-index locality of $\tilde V^W$"); ax[0].legend(fontsize=9); ax[0].grid(alpha=.3)
for q,lab in [(0,"q=(0,0)"),(idx[(2,0)],"q=(2,0)"),(idx[(3,3)],"q=(3,3)"),(idx[(6,6)],"q=(6,6)")]:
    d,nv=edecay(q); o=np.argsort(d); ax[1].semilogy(d[o],nv[o]+1e-14,".-",ms=4,label=lab)
ax[1].set_xlabel("|R_e| (prim. cells)"); ax[1].set_ylabel(r"$\|\tilde V^W(R_e;q)\|$ (Ry)")
ax[1].set_title(r"Electron-index decay of $\tilde V^W$ (new gauge)"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig("p5b_vtilde_locality.png",dpi=120); print("wrote p5b_vtilde_locality.png")
