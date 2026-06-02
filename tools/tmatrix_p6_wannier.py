#!/usr/bin/env python3
"""Active-space T-matrix in the Wannier basis (Koster-Slater) with a FINE-grid host G^A.

  T(omega) = [1 - V~^W G^A(omega)]^{-1} V~^W

V~^W : localized downfolded potential (coarse 12x12 block rotated with filukk_150, FT'd,
       truncated to the defect's compact support) -- the SMALL object that is inverted.
G^A  : host Green's function, G^A_{w'w}(dR;w) = (1/Nf^2) sum_{k fine} e^{2pi i k.dR}
       [(w+i eta) - H_W(k)]^{-1}, with H_W(k) Wannier-interpolated from hr.dat on an Nf x Nf grid.
The defect (V~^W) stays coarse/local; only the host (G^A) is taken to the fine grid -- the point
of the Wannier interpolation. Convergence of ||T|| with Nf is the deliverable.
"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994
HR="/anvil/scratch/x-rg47749/claude/sternheimerED/rewann_150/mos2_hr.dat"
ETA=(float(sys.argv[1]) if len(sys.argv)>1 else 0.05)/RY            # Ry
RCUT=int(sys.argv[2]) if len(sys.argv)>2 else 4                     # defect-centered cutoff (cells)

fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
fb.read_ints(np.int32); fb.read_reals(np.float64)
fb.read_record(np.complex128); fb.read_record(np.complex128)
V=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()
nb=N_A//nk; n1=int(round(nk**0.5)); V4=V.reshape(nk,nb,nk,nb)
fw=FortranFile("wann_data.dat","r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F")
U=fw.read_record(np.complex128).reshape((nb,nb,nk),order="F"); fw.close()
Ub=np.transpose(U,(2,0,1))

# --- V~^W(R',R) on the coarse R-grid, then defect-centered truncation ---
ms=np.arange(-n1//2,n1//2); Rg=np.array([(m,n) for m in ms for n in ms]); nR=len(Rg)
VW=np.einsum("fsp,fsgt,gtw->fpgw",Ub.conj(),V4,Ub,optimize=True)
F=np.exp(2j*np.pi*(Rg@xk[:2]))
VWR=np.einsum("Af,fpgw,Bg->ApBw",F,VW,F.conj(),optimize=True)/nk**2          # [R',i,R,j] Ry
onsite=np.array([np.linalg.norm(np.einsum("AiAj->ij",VWR[[r]][:, :, [r]])) for r in range(nR)])
Rd=int(np.argmax([np.linalg.norm(VWR[r,:,r,:]) for r in range(nR)]))         # defect cell
def mshift(R):                                # min-image displacement of R from Rd (cells)
    dm=((Rg[R,0]-Rg[Rd,0]+n1//2)%n1)-n1//2; dn=((Rg[R,1]-Rg[Rd,1]+n1//2)%n1)-n1//2
    return dm,dn
sel=[r for r in range(nR) if max(abs(mshift(r)[0]),abs(mshift(r)[1]))<=RCUT]
ns=len(sel); dim=nb*ns
print(f"defect cell R_d=({Rg[Rd,0]},{Rg[Rd,1]}); Rcut={RCUT} -> {ns} cells, subspace dim {dim}")
Vsub=np.zeros((dim,dim),complex)
for a,Ra in enumerate(sel):
    for b,Rb in enumerate(sel):
        Vsub[a*nb:(a+1)*nb,b*nb:(b+1)*nb]=VWR[Ra,:,Rb,:]

# --- hr.dat -> H_W(R) (eV) ---
L=open(HR).read().split("\n"); nw=int(L[1]); nr=int(L[2]); li=3; deg=[]
while len(deg)<nr: deg+=[int(x) for x in L[li].split()]; li+=1
deg=np.array(deg[:nr],float); Rhr=np.zeros((nr,3),int); Hr=np.zeros((nr,nw,nw),complex)
d=[x for x in L[li:] if x.strip()]; k=0
for ir in range(nr):
    for jn in range(nw):
        for im in range(nw):
            p=d[k].split(); k+=1; Rhr[ir]=(int(p[0]),int(p[1]),int(p[2]))
            Hr[ir,int(p[3])-1,int(p[4])-1]=float(p[5])+1j*float(p[6])
Rhr2=Rhr[:,:2].astype(float); Hr_ry=Hr/deg[:,None,None]/RY                     # Ry

def gA_fine(Nf, omega):
    kf=np.array([[i/Nf,j/Nf] for i in range(Nf) for j in range(Nf)])           # [Nf^2,2]
    ph=np.exp(2j*np.pi*(kf@Rhr2.T))                                            # [k, nr]
    Hk=np.einsum("kr,rmn->kmn", ph, Hr_ry)                                     # H_W(k) Ry
    Ginv=np.linalg.inv((omega+1j*ETA)*np.eye(nb)[None]-Hk)                     # [k,nb,nb]
    # G^A(dR) for the unique dR = R'-R within the subspace
    dRset={}
    for Ra in sel:
        for Rb in sel:
            dRset[(Rg[Ra,0]-Rg[Rb,0],Rg[Ra,1]-Rg[Rb,1])]=None
    out={}
    for (dm,dn) in dRset:
        phr=np.exp(2j*np.pi*(kf[:,0]*dm+kf[:,1]*dn))
        out[(dm,dn)]=np.einsum("k,kmn->mn",phr,Ginv)/Nf**2
    return out

def Tnorm(Nf, omega):
    g=gA_fine(Nf, omega); Gsub=np.zeros((dim,dim),complex)
    for a,Ra in enumerate(sel):
        for b,Rb in enumerate(sel):
            Gsub[a*nb:(a+1)*nb,b*nb:(b+1)*nb]=g[(Rg[Ra,0]-Rg[Rb,0],Rg[Ra,1]-Rg[Rb,1])]
    T=np.linalg.solve(np.eye(dim)-Vsub@Gsub, Vsub)
    return np.linalg.norm(T)

print(f"\nWannier-basis active T-matrix, omega=omega0(VBM)={omega0*RY:.4f} eV, eta={ETA*RY:.3g} eV")
print("  Nf x Nf     k-points    ||T||_F (Ry, Wannier conv.)")
Nfs=[12,24,36,48,72,96]; Tn=[]
for Nf in Nfs:
    t=Tnorm(Nf, omega0); Tn.append(t)
    tag=" (= coarse P5-a)" if Nf==12 else ""
    print(f"   {Nf:3d}x{Nf:<3d}    {Nf*Nf:7d}     {t:10.4f}{tag}")

fig,ax=plt.subplots(figsize=(6.4,4.4))
ax.plot(Nfs, Tn, "o-", c="#1f77b4")
ax.axvline(12, ls=":", c="gray"); ax.annotate("coarse 12x12\n(= P5-a)", (12, Tn[0]), xytext=(20,Tn[0]*0.9), fontsize=9)
ax.set_xlabel("fine k-grid $N_f$ (host $G^A$ on $N_f\\times N_f$)")
ax.set_ylabel(r"$\|T_{PP}(\omega_0)\|_F$  (Ry, Wannier conv.)")
ax.set_title("Wannier-basis active T-matrix: convergence with the host $G^A$ k-grid")
ax.grid(alpha=.3); plt.tight_layout(); plt.savefig("p6_wannier_converge.png",dpi=120)
print("wrote p6_wannier_converge.png")
