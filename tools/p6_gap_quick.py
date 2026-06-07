#!/usr/bin/env python3
"""Quick un-projected check: ||T_sub(omega)|| at a handful of energies across the gap, V~ vs M.
If a mid-gap resonance exists in the full Wannier T (but was lost in the band-projection of A),
||T_sub|| will spike there.  Only ~9 energies x 2 inputs -> a few solves."""
import numpy as np
from scipy.io import FortranFile
RY=13.605693122994; ETA=0.05/RY; RCUT=4; NF=48
fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2band=idx[N_A:]; eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.read_record(np.complex128)
V=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
nb=N_A//nk; n1=int(round(nk**0.5)); M4=M.reshape(nk,nb,nk,nb); V4=V.reshape(nk,nb,nk,nb); Ek=eact.reshape(nk,nb)
VBM=Ek[:,6].max(); CBM=Ek[:,7].min()
fw=FortranFile("wann_data.dat","r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F")
U=fw.read_record(np.complex128).reshape((nb,nb,nk),order="F"); fw.close()
Ub=np.transpose(U,(2,0,1))
print(f"VBM={VBM*RY:.3f} CBM={CBM*RY:.3f} gap={(CBM-VBM)*RY:.3f} eV (omega0={omega0*RY:.3f})",flush=True)
ms=np.arange(-n1//2,n1//2); Rg=np.array([(m,n) for m in ms for n in ms]); nR=len(Rg)
F=np.exp(2j*np.pi*(Rg@xk[:2]))
def toR(O4):
    OW=np.einsum("fsp,fsgt,gtw->fpgw",Ub.conj(),O4,Ub,optimize=True)
    return np.einsum("Af,fpgw,Bg->ApBw",F,OW,F.conj(),optimize=True)/nk**2
VWR=toR(V4); MWR=toR(M4)
HWk=np.einsum("ksp,ks,ksw->kpw",Ub.conj(),Ek,Ub,optimize=True)
HWR=np.einsum("Rk,kpw->Rpw",np.exp(-2j*np.pi*(Rg@xk[:2])),HWk,optimize=True)/nk
Rd=int(np.argmax([np.linalg.norm(VWR[r,:,r,:]) for r in range(nR)]))
def msh(R): return ((Rg[R,0]-Rg[Rd,0]+n1//2)%n1)-n1//2, ((Rg[R,1]-Rg[Rd,1]+n1//2)%n1)-n1//2
sel=[r for r in range(nR) if max(abs(msh(r)[0]),abs(msh(r)[1]))<=RCUT]; ns=len(sel); dim=nb*ns
def sub(WR):
    S=np.zeros((dim,dim),complex)
    for a in range(ns):
        for b in range(ns): S[a*nb:(a+1)*nb,b*nb:(b+1)*nb]=WR[sel[a],:,sel[b],:]
    return S
Vsub=sub(VWR); Msub=sub(MWR)
kfg=np.array([[i/NF,j/NF] for i in range(NF) for j in range(NF)])
Hf=np.einsum("kR,Rpw->kpw",np.exp(2j*np.pi*(kfg@Rg.T)),HWR); efine,Wfine=np.linalg.eigh(Hf)
dRs=sorted({(Rg[sel[a],0]-Rg[sel[b],0],Rg[sel[a],1]-Rg[sel[b],1]) for a in range(ns) for b in range(ns)})
dRidx={dR:i for i,dR in enumerate(dRs)}
phdR=np.array([np.exp(2j*np.pi*(kfg[:,0]*dm+kfg[:,1]*dn)) for (dm,dn) in dRs])
dRmat=np.array([[dRidx[(Rg[sel[a],0]-Rg[sel[b],0],Rg[sel[a],1]-Rg[sel[b],1])] for b in range(ns)] for a in range(ns)])
def gt_of(omega):
    D=1.0/((omega+1j*ETA)-efine); G=np.einsum("kpn,kn,kqn->kpq",Wfine,D,Wfine.conj())
    return (np.einsum("dk,kpq->dpq",phdR,G)/NF**2)[dRmat].transpose(0,2,1,3).reshape(dim,dim)
print("setup done; scanning",flush=True)
oms=np.linspace(VBM*RY-0.1, CBM*RY+0.1, 9)/RY
for X,tag in [(Vsub,"V~"),(Msub,"M ")]:
    print(f" {tag}:",flush=True)
    for om in oms:
        T=np.linalg.solve(np.eye(dim)-X@gt_of(om),X)
        sv=np.linalg.svd(T,compute_uv=False)[:3]
        print(f"   omega={om*RY:+.3f}eV ({(om-VBM)*RY:+.2f}>VBM)  ||T||={np.linalg.norm(T):8.1f}  top3SVD={np.array2string(sv,precision=1)}",flush=True)
