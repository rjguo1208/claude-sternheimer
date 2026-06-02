#!/usr/bin/env python3
"""Diagnose why V~^W is not localized: check (i) H_W(R) decay (Wannier-function
locality / gauge), (ii) Born M, rest Sigma, and V~ Wannierized shell decay."""
import numpy as np
from scipy.io import FortranFile
RY=13.605693122994
HR="/anvil/projects/x-che190065/rjguo/qe-7.5/T-matrix/data/mos2_hr.dat"

fb=FortranFile("vtilde_block.dat","r")
N_A,nkstot,nk_use,nbndskip=fb.read_ints(np.int32)
omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32)
eact=fb.read_reals(np.float64)
M =fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
Sg=fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
V =fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
fb.close()
fw=FortranFile("wann_data.dat","r")
nk2,nbndep,nw=fw.read_ints(np.int32); ibndkept=fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nkstot),order="F")
U=fw.read_record(np.complex128).reshape((nbndep,nw,nkstot),order="F")
fw.close()
nb=nw; Nk=nkstot; Ub=np.transpose(U,(2,0,1))

# hr -> H_W(R) decay
L=open(HR).read().split("\n"); nwh=int(L[1]); nr=int(L[2]); li=3; deg=[]
while len(deg)<nr: deg+=[int(x) for x in L[li].split()]; li+=1
deg=np.array(deg[:nr],float); Rhr=np.zeros((nr,3),int); Hr=np.zeros((nr,nb,nb),complex)
d=[x for x in L[li:] if x.strip()]; k=0
for ir in range(nr):
    for n in range(nb):
        for m in range(nb):
            p=d[k].split(); k+=1; Rhr[ir]=(int(p[0]),int(p[1]),int(p[2]))
            Hr[ir,int(p[3])-1,int(p[4])-1]=float(p[5])+1j*float(p[6])
print("=== H_W(R) decay (Wannier-function locality), eV ===")
shr=np.maximum(np.abs(Rhr[:,0]),np.abs(Rhr[:,1]))
for dsh in range(7):
    sel=shr==dsh
    if sel.any(): print(f"   |R|={dsh}: max||H_W(R)/ndeg|| = {np.max([np.linalg.norm(Hr[i]/deg[i]) for i in np.where(sel)[0]])*RY:.4e}")

ms=np.arange(-6,6); Rg=np.array([(m,n,0) for m in ms for n in ms]); nR=len(Rg)
F=np.exp(2j*np.pi*(Rg@xk)); sh=np.maximum(np.abs(Rg[:,0]),np.abs(Rg[:,1]))
def shells(O,name):
    O4=O.reshape(nkstot,nb,nkstot,nb)
    OW=np.einsum("fsp,fsgt,gtw->fpgw",Ub.conj(),O4,Ub)
    OWR=np.einsum("Af,fpgw,Bg->ApBw",F,OW,F.conj())/Nk**2
    bl=np.linalg.norm(OWR,axis=(1,3))
    print(f"\n=== {name}: ||O^W(R',R)||_F by shell max(|R'|,|R|), Ry  (total {np.linalg.norm(OWR):.3f}) ===")
    for dd in range(7):
        mask=(sh[:,None]==dd)|(sh[None,:]==dd)
        if mask.any(): print(f"   d={dd}: max {bl[mask].max():.4e}   diag(R,R) {bl[np.diag_indices(nR)][sh==dd].max() if (sh==dd).any() else 0:.4e}")
shells(M,"Born M")
shells(Sg,"rest Sigma")
shells(V,"Vtilde = M+Sigma")
