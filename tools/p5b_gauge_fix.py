#!/usr/bin/env python3
"""Before/after of the gauge fix: electron-index decay ||M^W(R_e;q)|| with the
OLD filukk (from the separate 17-band Wannierization) vs the NEW filukk_150
(re-Wannierized on the 150-band NSCF, same Wannier space)."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994
fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
fb.read_ints(np.int32); fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
nb=N_A//nk; n1=int(round(nk**0.5)); M4=M.reshape(nk,nb,nk,nb)

def load_U(fn):
    fw=FortranFile(fn,"r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
    xk=fw.read_reals(np.float64).reshape((3,nk),order="F")
    U=fw.read_record(np.complex128).reshape((nb,nb,nk),order="F"); fw.close()
    return np.transpose(U,(2,0,1)), xk

def decay(Ub, xk, qi, kadd):
    VW=np.einsum("fsp,fsgt,gtw->fpgw", Ub.conj(), M4, Ub)
    edms=np.array([VW[kadd(k,qi),:,k,:] for k in range(nk)])
    ms=np.arange(-n1//2,n1//2); Rg=np.array([(m,n,0) for m in ms for n in ms])
    F=np.exp(-2j*np.pi*(Rg@xk)); MR=np.einsum("Rk,kpw->Rpw",F,edms)/nk
    return np.sqrt(Rg[:,0]**2+Rg[:,1]**2), np.linalg.norm(MR,axis=(1,2))

fig,ax=plt.subplots(1,2,figsize=(12,4.6),sharey=True)
for p,(fn,ttl) in enumerate([("wann_data_17gauge.dat","OLD: filukk from separate 17-band run"),
                              ("wann_data.dat","NEW: re-Wannierized on 150-band NSCF (same space)")]):
    Ub,xk=load_U(fn)
    ij=np.rint(np.array([xk[0]*n1,xk[1]*n1]).T).astype(int)%n1
    idx={(ij[k,0],ij[k,1]):k for k in range(nk)}
    kadd=lambda k,q: idx[((ij[k,0]+ij[q,0])%n1,(ij[k,1]+ij[q,1])%n1)]
    for q,lab in [(0,"q=(0,0)"),(idx[(2,0)],"q=(2,0)"),(idx[(3,3)],"q=(3,3)"),(idx[(6,6)],"q=(6,6)")]:
        d,n=decay(Ub,xk,q,kadd); o=np.argsort(d)
        ax[p].semilogy(d[o],n[o]+1e-14,".-",ms=4,label=lab)
    ax[p].set_xlabel("|R_e| (prim. cells)"); ax[p].set_title(ttl); ax[p].grid(alpha=.3); ax[p].legend(fontsize=8)
ax[0].set_ylabel(r"$\|M^W(R_e;q)\|$  (Ry)")
plt.tight_layout(); plt.savefig("p5b_gauge_fix.png",dpi=115); print("wrote p5b_gauge_fix.png")
