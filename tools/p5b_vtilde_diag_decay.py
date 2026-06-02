#!/usr/bin/env python3
"""On-site (R=R') downfolded potential: |Vtilde^W_ij(R,R)| for a fixed Wannier pair (i,j),
vs distance of cell R from the defect center, with the gauge-consistent filukk_150."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994; BOHR2ANG=0.529177
fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); fb.read_reals(np.float64)
fb.read_ints(np.int32); fb.read_reals(np.float64)
fb.read_record(np.complex128); fb.read_record(np.complex128)        # skip M, Sigma
V=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
nb=N_A//nk; n1=int(round(nk**0.5)); V4=V.reshape(nk,nb,nk,nb)
fw=FortranFile("wann_data.dat","r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F")
U=fw.read_record(np.complex128).reshape((nb,nb,nk),order="F"); fw.close()
Ub=np.transpose(U,(2,0,1))

# rotate to Wannier, FT both indices -> Vtilde^W(R',i, R,j)
ms=np.arange(-n1//2,n1//2); Rg=np.array([(m,n) for m in ms for n in ms]); nR=len(Rg)
VW=np.einsum("fsp,fsgt,gtw->fpgw",Ub.conj(),V4,Ub,optimize=True)
F=np.exp(2j*np.pi*(Rg@xk[:2]))                                      # [R,k]
VWR=np.einsum("Af,fpgw,Bg->ApBw",F,VW,F.conj(),optimize=True)/nk**2 # [R',i,R,j]
diag=np.einsum("AiAj->Aij",VWR)                                     # on-site block [R,i,j] (Ry)

# defect center = cell where the on-site block is largest
nrmR=np.linalg.norm(diag,axis=(1,2)); Rd=int(np.argmax(nrmR))
# primitive lattice from the V_d.cube (supercell / 6), Bohr -> Ang
a1=np.array([240*0.150478, 0.0])/6*BOHR2ANG
a2=np.array([240*-0.075239, 240*0.130318])/6*BOHR2ANG
# minimum-image distance from the defect (the 12x12 cell repeats it periodically)
def mindist(R):
    d=np.inf
    for p in (-1,0,1):
        for q in (-1,0,1):
            dm=Rg[R,0]-Rg[Rd,0]+n1*p; dn=Rg[R,1]-Rg[Rd,1]+n1*q
            d=min(d, np.linalg.norm(dm*a1+dn*a2))
    return d
dist=np.array([mindist(R) for R in range(nR)])                     # Ang, minimum image
print(f"defect center cell R_d = ({Rg[Rd,0]},{Rg[Rd,1]}); |a1|={np.linalg.norm(a1):.3f} Ang; "
      f"max min-image dist = {dist.max():.1f} Ang")

# dominant on-site element at the defect, and a diagonal i=i
ij_abs=np.abs(diag[Rd]); i0,j0=np.unravel_index(np.argmax(ij_abs),ij_abs.shape)
print(f"dominant on-site pair at R_d: (i,j)=({i0},{j0}), |Vtilde^W_ij(Rd,Rd)|={ij_abs[i0,j0]:.4f} Ry")

o=np.argsort(dist)
fig,ax=plt.subplots(figsize=(6.8,4.8))
for (i,j,c,lab) in [(i0,j0,"#1f77b4",f"fixed (i,j)=({i0+1},{j0+1})"),
                    (0,0,"#2ca02c","(i,j)=(1,1)"),
                    (0,1,"#ff7f0e","(i,j)=(1,2)")]:
    y=np.abs(diag[:,i,j])
    ax.semilogy(dist[o], y[o]+1e-12, "o-" if (i,j)==(i0,j0) else ".--", ms=5 if (i,j)==(i0,j0) else 3,
                c=c, lw=1.8 if (i,j)==(i0,j0) else 1, label=lab)
# log-linear decay length on the dominant element (over the decaying tail)
y0=np.abs(diag[:,i0,j0]); m=(dist>0.1)&(y0>1e-6)
if m.sum()>3:
    p=np.polyfit(dist[m], np.log(y0[m]),1); lam=-1/p[0]
    ax.plot(np.sort(dist[m]), np.exp(p[1]+p[0]*np.sort(dist[m])),"k:",lw=1,label=f"exp fit, $\\lambda$≈{lam:.1f} Å")
ax.set_xlabel("distance of cell $R$ from defect center  (Å)")
ax.set_ylabel(r"$|\tilde V^W_{ij}(R,R)|$  (Ry)")
ax.set_title(r"On-site downfolded potential $|\tilde V^W_{ij}(R,R)|$ vs distance from the defect")
ax.legend(fontsize=9); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig("p5b_vtilde_diag_decay.png",dpi=120); print("wrote p5b_vtilde_diag_decay.png")
