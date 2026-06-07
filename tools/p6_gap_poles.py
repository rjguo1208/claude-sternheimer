#!/usr/bin/env python3
"""Fast Koster-Slater gap scan: how many S-vacancy defect resonances sit in the gap, and does the
rest dressing change the count?  For with-rest (V~) and no-rest (M):
  - ||T_sub(omega)||_F and max|Im T_ii(omega)|  (cheap resonance indicators, no eigendecomposition)
  - at each detected peak, the leading singular values of T_sub give the degeneracy
    (a1 singlet -> 1 large SV; e doublet -> 2 comparable large SVs).
This is exactly the machinery behind the spectral figure (Rcut=4, fine NF=48 host G)."""
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
print(f"VBM={VBM*RY:.3f}  CBM={CBM*RY:.3f}  gap={(CBM-VBM)*RY:.3f} eV  (rest-ref omega0={omega0*RY:.3f} eV)",flush=True)
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
def Tof(om,X): return np.linalg.solve(np.eye(dim)-X@gt_of(om),X)

oms=np.linspace(VBM*RY-0.3, CBM*RY+0.2, 121)/RY; omeV=oms*RY
print(f"\nKoster-Slater gap scan [{omeV[0]:.2f},{omeV[-1]:.2f}] eV, {len(oms)} pts (dim={dim}):",flush=True)
for X,tag in [(Vsub,"V~ with-rest"),(Msub,"M  no-rest ")]:
    nrm=np.empty(len(oms)); imt=np.empty(len(oms))
    for j,om in enumerate(oms):
        T=Tof(om,X); nrm[j]=np.linalg.norm(T); imt[j]=np.abs(np.imag(np.diag(T))).max()
    med=np.median(nrm); medi=max(np.median(imt),1e-9)
    pk=[i for i in range(2,len(oms)-2) if imt[i]==max(imt[i-2:i+3]) and imt[i]>3*medi]
    print(f"\n {tag}:  median||T||={med:.0f}",flush=True)
    print(f"   resonances (peaks in max|Im T_ii|):",flush=True)
    if pk:
        for i in pk:
            sv=np.linalg.svd(Tof(oms[i],X),compute_uv=False)[:4]
            deg=int(np.sum(sv>0.5*sv[0]))
            print(f"     omega={omeV[i]:+.3f} eV ({(omeV[i]-VBM*RY):+.2f} above VBM)  |ImT|max={imt[i]:7.1f}  topSVD={np.array2string(sv,precision=1)} ~deg{deg}",flush=True)
    else:
        print("     none above threshold",flush=True)
    print("   profile  (omega:  ||T|| / |ImT|max):",flush=True)
    for i in range(0,len(oms),4):
        b="#"*int(min(imt[i]/medi*1.0,60))
        print(f"     {omeV[i]:+.3f}  {nrm[i]:8.1f} / {imt[i]:8.2f}  {b}",flush=True)
