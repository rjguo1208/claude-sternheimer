#!/usr/bin/env python3
"""Reproduce EDI's edbloch2wane decay (electron-index Wannierization at fixed q)
from our M block, and contrast with the both-index (Koster-Slater) transform.

EDI: edms(k) = cu(k)^dag M(k,k+q) cuq(k+q);  M^W(R_e;q) = (1/Nk) sum_k e^{-ik.R} edms(k).
decay.M = max_{m'm}|M^W(R_e;q)| vs |R_e|.  This is the *electron* Wannier decay (fast).
"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994
fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); eact=fb.read_reals(np.float64)
M =fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
Sg=fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
V =fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
nb=N_A//nk
fw=FortranFile("wann_data.dat","r"); _=fw.read_ints(np.int32); _=fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F")
U=fw.read_record(np.complex128).reshape((nb,nb,nk),order="F"); fw.close()
Ub=np.transpose(U,(2,0,1))                      # [k,band,wann]
n1=int(round(nk**0.5))                          # 12

# k-index -> integer (i1,i2) on the 12x12 grid, and back
ij=np.rint(np.array([xk[0]*n1, xk[1]*n1]).T).astype(int)%n1     # [k,2]
idx_of={(ij[k,0],ij[k,1]):k for k in range(nk)}
def kadd(ki,qi):                                # index of k_i + q
    return idx_of[((ij[ki,0]+ij[qi,0])%n1,(ij[ki,1]+ij[qi,1])%n1)]

# verify k-addition table: xk[kadd(ki,q)] - xk[ki] - xk[q] = 0 (mod 1)
kadd_err=0.0
for qi in [0, idx_of[(2,0)], idx_of[(3,3)]]:
    for ki in range(nk):
        d=xk[:,kadd(ki,qi)]-xk[:,ki]-xk[:,qi]
        kadd_err=max(kadd_err, np.max(np.abs(d-np.rint(d))))
print(f"kadd verification: max|xk[k+q]-xk[k]-xk[q] mod 1| = {kadd_err:.2e}  ({'OK' if kadd_err<1e-6 else 'BAD'})")

# rotate the WHOLE block to Wannier gauge once: VW[kf,wf,ki,wi]=U^dag(kf) O(kf,ki) U(ki)
def rotate(O):
    O4=O.reshape(nk,nb,nk,nb)
    return np.einsum("fsp,fsgt,gtw->fpgw", Ub.conj(), O4, Ub)
VWm=rotate(M)

# EDI-style electron-index decay for a fixed q: FT over electron k=ki of VW[ki+q,:,ki,:]
def edi_decay(VW, qi):
    edms=np.array([VW[kadd(ki,qi),:,ki,:] for ki in range(nk)])   # [ki,wf,wi]
    # FT over ki -> R_e on the 12x12 R-grid
    ms=np.arange(-n1//2,n1//2); Rg=np.array([(m,n,0) for m in ms for n in ms])
    F=np.exp(-2j*np.pi*(Rg@xk))                                   # [R,k]
    MR=np.einsum("Rk,kpw->Rpw", F, edms)/nk                       # [R,wf,wi]
    dist=np.sqrt(Rg[:,0]**2+Rg[:,1]**2)
    nrm=np.linalg.norm(MR,axis=(1,2))
    return dist, nrm

print("EDI-style electron-index decay  ||M^W(R_e; q)||  (max over R_e per shell, Ry):")
qlist=[0, idx_of[(2,0)], idx_of[(3,3)], idx_of[(6,6)]]
fig,ax=plt.subplots(figsize=(6.4,4.6))
for qi in qlist:
    d,n=edi_decay(VWm,qi)
    o=np.argsort(d);
    lab="q=(%d,%d)/12"%(ij[qi,0],ij[qi,1])
    ax.semilogy(d[o],n[o]+1e-14,".-",ms=4,label=lab)
    # shell summary
    s0=n[d==0][0] if (d==0).any() else 0.0
    far=n[d>=5].max() if (d>=5).any() else 0.0
    print(f"  {lab:14s}: R_e=0 {s0:.4e}  ->  |R_e|>=5 {far:.4e}  Ry   (ratio {s0/max(far,1e-30):.1f})")
ax.set_xlabel("|R_e|  (prim. cells)"); ax.set_ylabel(r"$\|M^W(R_e;q)\|$  (Ry)")
ax.set_title("EDI-style electron-index decay of M (our block)"); ax.legend(); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig("p5b_edi_decay.png",dpi=120); print("wrote p5b_edi_decay.png")
