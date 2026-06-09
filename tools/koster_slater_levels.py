#!/usr/bin/env python3
"""Koster-Slater defect levels (a1 + e) from the BARE Wannier potential M^W.

Levels = in-gap roots of det[1 - g(E) dV_D] = 0:
  dV_D = bare M (1st record of vtilde_block.dat) rotated to Wannier (filukk/wann_data) + FT'd to R
         + truncated to the defect-centered block (R_cut).   [NOT the over-screening dressed Vtilde]
  g(E) = host Green's function in that block, fine-grid Wannier-interpolated from hr.dat:
         g(dR;E) = (1/Nf^2) sum_k e^{2pi i k.dR} [(E+i eta) - H_W(k)]^{-1}.
Scan E across the gap; track the smallest singular values of [1 - g(E)dV_D] (dip at a level;
# of small sv = degeneracy) and the Krein-Friedel count dN(E) = -(1/pi) Im ln det (staircase:
+1 step at a1, +2 step at e).  Args: eta_eV[=0.02] Rcut[=3] Nf[=48].
"""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994
HR="/anvil/scratch/x-rg47749/claude/sternheimerED/rewann_150/mos2_hr.dat"
ETA=(float(sys.argv[1]) if len(sys.argv)>1 else 0.02)/RY          # Ry
RCUT=int(sys.argv[2]) if len(sys.argv)>2 else 3
Nf=int(sys.argv[3]) if len(sys.argv)>3 else 48

# --- vtilde_block.dat: header, idx, eact, BARE M (1st complex record) ---
fb=FortranFile("vtilde_block.dat","r")
N_A,nk,nk_use,nbndskip=fb.read_ints(np.int32); omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2band=idx[N_A:2*N_A]
eact=fb.read_reals(np.float64)
M=fb.read_record(np.complex128).reshape((N_A,N_A),order="F"); fb.close()   # BARE M (Ry)
eps=eact*RY; VBM=eps[a2band<=13].max(); CBM=eps[a2band>=14].min()
nb=N_A//nk; n1=int(round(nk**0.5)); M4=M.reshape(nk,nb,nk,nb)
print(f"N_A={N_A} nk={nk} nb={nb}  VBM={VBM:.3f} CBM={CBM:.3f} gap={CBM-VBM:.3f} eV  (DFT: a1 -5.83, e -4.77)")

# --- bare M^W(R',R), defect-block truncation (recipe from tmatrix_p6_wannier.py) ---
fw=FortranFile("wann_data.dat","r"); fw.read_ints(np.int32); fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F")
U=fw.read_record(np.complex128).reshape((nb,nb,nk),order="F"); fw.close()
Ub=np.transpose(U,(2,0,1))
ms=np.arange(-n1//2,n1//2); Rg=np.array([(m,n) for m in ms for n in ms]); nR=len(Rg)
MW=np.einsum("fsp,fsgt,gtw->fpgw",Ub.conj(),M4,Ub,optimize=True)
F=np.exp(2j*np.pi*(Rg@xk[:2]))
MWR=np.einsum("Af,fpgw,Bg->ApBw",F,MW,F.conj(),optimize=True)/nk**2          # [R',i,R,j] Ry (bare)
Rd=int(np.argmax([np.linalg.norm(MWR[r,:,r,:]) for r in range(nR)]))
def mshift(R):
    dm=((Rg[R,0]-Rg[Rd,0]+n1//2)%n1)-n1//2; dn=((Rg[R,1]-Rg[Rd,1]+n1//2)%n1)-n1//2
    return dm,dn
sel=[r for r in range(nR) if max(abs(mshift(r)[0]),abs(mshift(r)[1]))<=RCUT]
ns=len(sel); dim=nb*ns
print(f"defect cell R_d=({Rg[Rd,0]},{Rg[Rd,1]})  Rcut={RCUT} -> {ns} cells, block dim {dim}")
Vsub=np.zeros((dim,dim),complex)
for a,Ra in enumerate(sel):
    for b,Rb in enumerate(sel):
        Vsub[a*nb:(a+1)*nb,b*nb:(b+1)*nb]=MWR[Ra,:,Rb,:]                     # bare M^W block (Ry)

# --- hr.dat -> H_W(R) (Ry) ---
L=open(HR).read().split("\n"); nw=int(L[1]); nr=int(L[2]); li=3; deg=[]
while len(deg)<nr: deg+=[int(x) for x in L[li].split()]; li+=1
deg=np.array(deg[:nr],float); Rhr=np.zeros((nr,3),int); Hr=np.zeros((nr,nw,nw),complex)
dd=[x for x in L[li:] if x.strip()]; kk=0
for ir in range(nr):
    for jn in range(nw):
        for im in range(nw):
            p=dd[kk].split(); kk+=1; Rhr[ir]=(int(p[0]),int(p[1]),int(p[2]))
            Hr[ir,int(p[3])-1,int(p[4])-1]=float(p[5])+1j*float(p[6])
Rhr2=Rhr[:,:2].astype(float); Hr_ry=Hr/deg[:,None,None]/RY
kf=np.array([[i/Nf,j/Nf] for i in range(Nf) for j in range(Nf)])
phk=np.exp(2j*np.pi*(kf@Rhr2.T)); Hk=np.einsum("kr,rmn->kmn",phk,Hr_ry)     # H_W(k) Ry
dRset=sorted({(Rg[Ra,0]-Rg[Rb,0],Rg[Ra,1]-Rg[Rb,1]) for Ra in sel for Rb in sel})
phdR={dR: np.exp(2j*np.pi*(kf[:,0]*dR[0]+kf[:,1]*dR[1])) for dR in dRset}

def Gsub(E_ry):                                # host GF in the defect block at energy E
    Ginv=np.linalg.inv((E_ry+1j*ETA)*np.eye(nb)[None]-Hk)
    g={dR: np.einsum("k,kmn->mn",phdR[dR],Ginv)/Nf**2 for dR in dRset}
    G=np.zeros((dim,dim),complex)
    for a,Ra in enumerate(sel):
        for b,Rb in enumerate(sel):
            G[a*nb:(a+1)*nb,b*nb:(b+1)*nb]=g[(Rg[Ra,0]-Rg[Rb,0],Rg[Ra,1]-Rg[Rb,1])]
    return G

# --- scan E over the gap: singular values of [1 - g(E)Vsub] + Krein-Friedel dN ---
NE=400; Es=np.linspace(VBM-0.05, CBM+0.05, NE)
sv1=np.empty(NE); sv2=np.empty(NE); ph=np.empty(NE)
for j,E in enumerate(Es):
    A=np.eye(dim)-Gsub(E/RY)@Vsub
    s=np.linalg.svd(A,compute_uv=False); ss=np.sort(s)
    sv1[j],sv2[j]=ss[0],ss[1]
    sgn,_=np.linalg.slogdet(A); ph[j]=np.angle(sgn)
dN=-np.unwrap(ph)/np.pi; dN-=dN[0]                                          # Krein-Friedel count

# levels = local minima of sv1 inside the gap; degeneracy from sv2/sv1
gap=(Es>VBM+2e-3)&(Es<CBM-2e-3)
mins=[j for j in range(1,NE-1) if gap[j] and sv1[j]<sv1[j-1] and sv1[j]<=sv1[j+1] and sv1[j]<0.5]
print("\nKoster-Slater in-gap roots of det[1 - g(E) M^W_D]=0  (eta=%.3g eV, Rcut=%d, Nf=%d):"%(ETA*RY,RCUT,Nf))
for j in mins:
    deg=2 if sv2[j]<2.5*sv1[j] else 1
    lab={1:"a1 (singlet)",2:"e (doublet)"}[deg]
    print(f"   E = {Es[j]:.3f} eV  ({Es[j]-VBM:+.2f} above VBM)   sv1={sv1[j]:.3f} sv2={sv2[j]:.3f} -> {lab}")
print(f"   total Krein-Friedel in-gap count dN(CBM)-dN(VBM) = {dN[gap][-1]-dN[gap][0]:.2f}  (a1+e doublet = 3)")

# --- plot ---
fig,(ax1,ax2)=plt.subplots(2,1,figsize=(8,6.4),sharex=True)
ax1.axvspan(Es[0],VBM,color="#cfe0f5",alpha=.5); ax1.axvspan(CBM,Es[-1],color="#f6d3cf",alpha=.5)
ax1.plot(Es,sv1,c="#1f3b73",lw=1.7,label=r"$\sigma_{\min}[\,1-g(E)M^W_D\,]$")
ax1.plot(Es,sv2,c="#8aa6d6",lw=1.0,ls="--",label=r"$\sigma_{2}$")
for j in mins:
    ax1.axvline(Es[j],c="#e36414",lw=1,ls=":")
for e,nm in [(-5.83,"DFT $a_1$"),(-4.77,"DFT $e$")]:
    ax1.axvline(e,c="#999",lw=.8,ls=(0,(1,2)))
ax1.set_ylabel("singular value"); ax1.legend(fontsize=8,loc="upper center")
ax1.set_title(f"Koster–Slater from bare $M^W$ (Rcut={RCUT}, $N_f$={Nf}, $\\eta$={ETA*RY:.2g} eV): in-gap roots = $a_1$+$e$")
ax2.axvspan(Es[0],VBM,color="#cfe0f5",alpha=.5); ax2.axvspan(CBM,Es[-1],color="#f6d3cf",alpha=.5)
ax2.plot(Es,dN,c="#1f3b73",lw=1.7)
ax2.set_xlabel("$E$ (eV)"); ax2.set_ylabel(r"Krein–Friedel count $\Delta N(E)$")
ax2.set_xlim(Es[0],Es[-1])
plt.tight_layout(); plt.savefig("koster_slater_levels.png",dpi=140,bbox_inches="tight")
print("wrote koster_slater_levels.png")
