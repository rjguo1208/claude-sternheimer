#!/usr/bin/env python3
"""Diagnose the q-structure of the Born matrix element M(k_f,k_i) from the block.
Isolated localized defect -> M smooth in q=k_f-k_i (peaked at q=0, decaying).
Periodic-array / non-smooth -> comb or flat -> delocalized M^W."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.io import FortranFile
RY=13.605693122994
fb=FortranFile("vtilde_block.dat","r")
N_A,nkstot,nk_use,nbndskip=fb.read_ints(np.int32)
omega0,wmin,wmax=fb.read_reals(np.float64)
idx=fb.read_ints(np.int32); a2k=idx[:N_A]; a2band=idx[N_A:]
eact=fb.read_reals(np.float64)
M =fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
Sg=fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
V =fb.read_record(np.complex128).reshape((N_A,N_A),order="F")
fb.close()
nk=nkstot; nb=N_A//nk
M4=M.reshape(nk,nb,nk,nb)            # [kf,bf,ki,bi]
fw=FortranFile("wann_data.dat","r"); _=fw.read_ints(np.int32); _=fw.read_ints(np.int32)
xk=fw.read_reals(np.float64).reshape((3,nk),order="F"); fw.close()

# trace over bands: |M|(kf,ki) = ||M4[kf,:,ki,:]||_F
Mkk=np.linalg.norm(M4,axis=(1,3))    # [kf,ki]  Ry
# k-index -> (i1,i2) on 12x12 (QE order: i1=(n)//12, i2=n%12)
n1=int(round(nk**0.5))               # 12
def img(vec):                        # 144 -> 12x12
    return vec.reshape(n1,n1)

iG=0                                  # Gamma = k-index 1
# does M depend only on q=kf-ki ?  build M as fn of q by averaging over ki
# k crystal -> integer (i1,i2); q index = (i1f-i1i mod 12, i2f-i2i mod 12)
ij=np.array([( (n//n1), (n%n1) ) for n in range(nk)])
Mq=np.zeros((n1,n1)); cnt=np.zeros((n1,n1))
var_num=0.0; var_den=0.0
Mq_acc=np.zeros((n1,n1),complex)
# use band-diagonal average element for the q-collapse test
b0=nb-1
for kf in range(nk):
    for ki in range(nk):
        q=((ij[kf,0]-ij[ki,0])%n1, (ij[kf,1]-ij[ki,1])%n1)
        Mq[q]+=abs(M4[kf,b0,ki,b0]); cnt[q]+=1
Mq/=np.maximum(cnt,1)
# how well does |M| collapse onto a function of q? (spread within fixed q)
resid=0.0; tot=0.0
for kf in range(nk):
    for ki in range(nk):
        q=((ij[kf,0]-ij[ki,0])%n1, (ij[kf,1]-ij[ki,1])%n1)
        resid+=(abs(M4[kf,b0,ki,b0])-Mq[q])**2; tot+=abs(M4[kf,b0,ki,b0])**2
print(f"M depends-only-on-q test:  residual/total = {resid/tot:.3f}  (0 => M=f(q) exactly, like a periodic potential)")
print(f"|M(kf,ki)| diag(kf=ki) mean = {np.mean([Mkk[k,k] for k in range(nk)]):.4f} Ry ; off-diag mean = {(Mkk.sum()-sum(Mkk[k,k] for k in range(nk)))/(nk*nk-nk):.4f} Ry")
print(f"|M|(q=0) = {Mq[0,0]:.4f} ; max|M|(q!=0) = {Mq.copy().ravel()[1:].max() if Mq.size>1 else 0:.4f} ; min = {Mq.min():.4f} Ry")

fig,ax=plt.subplots(1,3,figsize=(15,4.4))
im0=ax[0].imshow(Mkk,origin="lower",aspect="auto"); plt.colorbar(im0,ax=ax[0])
ax[0].set_title("||M(k_f,k_i)||_F  (Ry)"); ax[0].set_xlabel("k_i index"); ax[0].set_ylabel("k_f index")
im1=ax[1].imshow(img(Mkk[:,iG]),origin="lower"); plt.colorbar(im1,ax=ax[1])
ax[1].set_title("||M(k_f, k_i=Gamma)|| on 12x12 BZ"); ax[1].set_xlabel("k_f,2"); ax[1].set_ylabel("k_f,1")
im2=ax[2].imshow(Mq,origin="lower"); plt.colorbar(im2,ax=ax[2])
ax[2].set_title("|M| collapsed onto q=k_f-k_i (band %d)"%b0); ax[2].set_xlabel("q2"); ax[2].set_ylabel("q1")
plt.tight_layout(); plt.savefig("p5b_qstruct.png",dpi=110); print("wrote p5b_qstruct.png")
