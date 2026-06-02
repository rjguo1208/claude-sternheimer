#!/usr/bin/env python3
"""P5-b: Wannierize V~, check locality, Koster-Slater T-matrix in the defect subspace.

  1) rotate V~(k_f,k_i) band->Wannier with U(k); FT to V~^W(R',R)  -> check decay (locality)
  2) build host G^A from H_W(k) (hr.dat):  G^A_{w'w}(dR;w) = (1/N_k) sum_k e^{2pi i k.dR}[(w+ieta)-H_W(k)]^{-1}
  3) truncated Koster-Slater:  T = [1 - V~^W G^A]^{-1} V~^W  on |R|<=Rcut  -> ||T|| converges to P5-a
Inputs: vtilde_block.dat, wann_data.dat (xk_cryst+U(k)), mos2_hr.dat .
"""
import sys, numpy as np
RY = 13.605693122994
from scipy.io import FortranFile
HR = sys.argv[1] if len(sys.argv) > 1 else \
     "/anvil/projects/x-che190065/rjguo/qe-7.5/T-matrix/data/mos2_hr.dat"
ETA = (float(sys.argv[2]) if len(sys.argv) > 2 else 0.05)/RY    # Ry

# ---- block ----
fb = FortranFile("vtilde_block.dat","r")
N_A,nkstot,nk_use,nbndskip = fb.read_ints(np.int32)
omega0,wmin,wmax = fb.read_reals(np.float64)
idx = fb.read_ints(np.int32); a2k=idx[:N_A]; a2band=idx[N_A:]
eact = fb.read_reals(np.float64)
M  = fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
Sg = fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
V  = fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
fb.close()
Nk=nkstot

# ---- wann_data: U(k), xk_cryst ----
fw = FortranFile("wann_data.dat","r")
nk2,nbndep,nw = fw.read_ints(np.int32)
ibndkept = fw.read_ints(np.int32)
xk = fw.read_reals(np.float64).reshape((3,nkstot),order="F")     # crystal coords
U  = fw.read_record(np.complex128).reshape((nbndep,nw,nkstot),order="F")  # [ib,iw,ik]
fw.close()
nb=nw
assert N_A==nb*nkstot, (N_A,nb,nkstot)
Ub = np.transpose(U,(2,0,1))      # [k, band, wann]

# ---- hr.dat -> H_W(R) ----
def read_hr(path):
    L=open(path).read().split("\n")
    nw_=int(L[1]); nr=int(L[2]); li=3; deg=[]
    while len(deg)<nr: deg+=[int(x) for x in L[li].split()]; li+=1
    deg=np.array(deg[:nr],float); Rv=np.zeros((nr,3),int); H=np.zeros((nr,nw_,nw_),complex)
    d=[x for x in L[li:] if x.strip()]; k=0
    for ir in range(nr):
        for n in range(nw_):
            for m in range(nw_):
                p=d[k].split(); k+=1
                Rv[ir]=(int(p[0]),int(p[1]),int(p[2]))
                H[ir,int(p[3])-1,int(p[4])-1]=float(p[5])+1j*float(p[6])
    return nw_,nr,deg,Rv,H
nwh,nr,ndeg,Rhr,Hr = read_hr(HR)
assert nwh==nb
def Hk(kc):
    ph=np.exp(2j*np.pi*(Rhr@kc))/ndeg
    return np.einsum("r,rmn->mn",ph,Hr)        # eV

# ---- (A) sanity: H_W(k) eigenvalues vs stored active energies ----
print("=== (A) H_W(k) eig vs eact  (validates k-order + hr<->block), eV ===")
maxd=0.0
for ik in [0, 37, 100, 143]:
    e_hk=np.sort(np.linalg.eigvalsh(Hk(xk[:,ik])).real)
    e_ed=np.sort(eact[a2k==(ik+1)]*RY)
    maxd=max(maxd, np.max(np.abs(e_hk-e_ed)))
    print(f"  k={ik+1:3d}  max|eig(H_W)-eact| = {np.max(np.abs(e_hk-e_ed)):.3e}")
print(f"  overall max = {maxd:.3e} eV  ({'OK' if maxd<1e-2 else 'MISMATCH -> k-order/hr problem'})")

# ---- rotate V~ to Wannier gauge, FT to real space ----
V4 = V.reshape(nkstot,nb,nkstot,nb)                         # [kf,sf,ki,si]
VW = np.einsum("fsp,fsgt,gtw->fpgw", Ub.conj(), V4, Ub)     # [kf,wf,ki,wi]
# real-space grid R=(m,n,0), m,n in -6..5  (144 vectors dual to the 12x12 k-grid)
ms=np.arange(-6,6); Rg=np.array([(m,n,0) for m in ms for n in ms])   # (144,3)
nR=len(Rg)
F=np.exp(2j*np.pi*(Rg@xk))            # [R,k]
VWR=np.einsum("Af,fpgw,Bg->ApBw", F, VW, F.conj())/Nk**2    # [R',wf,R,wi]  (Ry)

# ---- (B) locality: ||V~^W(R',R)|| by shell d=max(|R'|,|R|) ----
def shell(R): return int(max(abs(Rg[R,0]),abs(Rg[R,1])))
print("\n=== (B) Wannierized V~ locality:  ||V~^W(R',R)||_F by shell max(|R'|,|R|) ===")
bl=np.linalg.norm(VWR,axis=(1,3))     # [R',R]
sh=np.array([shell(i) for i in range(nR)])
print("   shell d   max||V~^W(R',R)||   (Ry)")
for d in range(7):
    mask=(sh[:,None]==d)|(sh[None,:]==d)
    if mask.any(): print(f"     {d:2d}        {bl[mask].max():.4e}")
i0 = [i for i in range(nR) if Rg[i,0]==0 and Rg[i,1]==0][0]
print(f"  on-site ||V~^W(0,0)|| = {bl[i0,i0]:.4e} Ry   (cf. ||V~||_F total over all R = {np.linalg.norm(VWR):.3f})")

# ---- (C/D) Koster-Slater: G^A(dR) then truncated inversion vs P5-a ----
# G^A_{w'w}(dR) = (1/Nk) sum_k e^{2pi i k.dR} [(omega0+i eta) - H_W(k)/?]^{-1}  ; work in Ry
Hk_all=np.array([Hk(xk[:,k])/RY for k in range(Nk)])          # [k,nb,nb] Ry
z=omega0+1j*ETA
Ginv=np.linalg.inv(z*np.eye(nb)[None]-Hk_all)                 # [k,nb,nb]
# G^A(dR): index by integer (dm,dn)
def gA(dm,dn):
    ph=np.exp(2j*np.pi*(xk[0]*dm+xk[1]*dn))                   # [k]
    return np.einsum("k,kmn->mn",ph,Ginv)/Nk
print("\n=== (C/D) truncated Koster-Slater  T=[1-V~^W G^A]^-1 V~^W  (omega0, eta=%.3g eV) ==="%(ETA*RY))
print("   Rcut   dim    ||T||_F (Ry)     (P5-a full ||T_PP||=297.1)")
for Rcut in range(7):
    sel=[i for i in range(nR) if shell(i)<=Rcut]
    ns=len(sel); dim=nb*ns
    Vs=np.zeros((dim,dim),complex); Gs=np.zeros((dim,dim),complex)
    for a,Ra in enumerate(sel):
        for b,Rb in enumerate(sel):
            Vs[a*nb:(a+1)*nb, b*nb:(b+1)*nb]=VWR[Ra,:,Rb,:]
            dm=Rg[Ra,0]-Rg[Rb,0]; dn=Rg[Ra,1]-Rg[Rb,1]
            Gs[a*nb:(a+1)*nb, b*nb:(b+1)*nb]=gA(dm,dn)
    T=np.linalg.solve(np.eye(dim)-Vs@Gs, Vs)
    print(f"    {Rcut:2d}   {dim:5d}   {np.linalg.norm(T):12.3f}")
