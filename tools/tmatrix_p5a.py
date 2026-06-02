#!/usr/bin/env python3
"""P5-a: coarse-grid active-space T-matrix from the downfolded potential.

  T_PP(omega) = [1 - V~ G^A(omega)]^{-1} V~        (beyond-Born: V~ = M + Sigma)
  T_M (omega) = [1 - M  G^A(omega)]^{-1} M          (Born input, for comparison)

  G^A is the HOST active Green's function, diagonal in the active Bloch basis,
  carrying the BZ measure 1/N_k (the same 1/N_k that lives in Sigma):
      G^A_a(omega) = (1/N_k) / (omega - eps_a + i*eta)

Reads vtilde_block.dat (eps_a = eact is stored there) so P5-a needs nothing else.

Validations:
  (1) weak-coupling  : V~ -> lambda V~ , lambda->0   =>  T/lambda -> V~   (code check)
  (2) large-eta      : eta -> inf  => G^A -> 0  => T -> V~                 (code check)
  (3) Born hierarchy : ||M|| , ||V~|| , ||T_M|| , ||T_PP||  +  ||T||/||V~|| (resummation)
  (4) resonance scan : ||T_PP(omega)|| over the active window
"""
import sys, numpy as np
from scipy.io import FortranFile

RY = 13.605693122994
ETA_EV = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05      # eV
path = sys.argv[1] if len(sys.argv) > 1 else "vtilde_block.dat"

f = FortranFile(path, "r")
N_A, nkstot, nk_use, nbndskip = f.read_ints(np.int32)
omega0, wmin, wmax = f.read_reals(np.float64)
idx = f.read_ints(np.int32); a2k = idx[:N_A]; a2band = idx[N_A:]
eact = f.read_reals(np.float64)                                  # active energies (Ry)
M  = f.read_record(np.complex128).reshape((N_A, N_A), order="F")
Sg = f.read_record(np.complex128).reshape((N_A, N_A), order="F")
V  = f.read_record(np.complex128).reshape((N_A, N_A), order="F")
f.close()

Nk = nkstot
eta = ETA_EV / RY
nrm = lambda A: np.linalg.norm(A)
print(f"N_A={N_A}  N_k={Nk}  omega0(VBM)={omega0*RY:.4f} eV  eta={ETA_EV} eV")
print(f"||M||={nrm(M):.3f}  ||Sigma||={nrm(Sg):.3f}  ||V~||={nrm(V):.3f} Ry")

def gA(omega):
    return (1.0/Nk) / (omega - eact + 1j*eta)        # length-N_A diagonal vector

def tmat(Vmat, g):
    A = np.eye(N_A) - Vmat * g[None, :]              # 1 - V~ diag(g)
    return np.linalg.solve(A, Vmat)

# ---- (3) Born hierarchy at omega0 (VBM, where carriers live) ----
g = gA(omega0)
Tpp = tmat(V, g)
Tm  = tmat(M, g)
print("\n--- at omega = omega0 (VBM) ---")
print(f"  ||T_M (Born input)|| = {nrm(Tm):.3f} Ry   ||T_M||/||M||  = {nrm(Tm)/nrm(M):.3f}")
print(f"  ||T_PP(beyond-Born)|| = {nrm(Tpp):.3f} Ry   ||T_PP||/||V~|| = {nrm(Tpp)/nrm(V):.3f}")
print(f"  ||T_PP - V~||/||V~|| = {nrm(Tpp-V)/nrm(V):.3f}   (active multiple-scattering effect)")
print(f"  ||T_PP - T_M||/||T_M|| = {nrm(Tpp-Tm)/nrm(Tm):.3f}   (beyond-Born vs Born, fully resummed)")

# ---- (1) weak-coupling: T/lambda -> V~ ----
lam = 1e-4
Tw = tmat(lam*V, g)
print("\n--- validation (1) weak coupling lambda=1e-4:  T/lambda -> V~ ---")
print(f"  ||T/lambda - V~|| / ||V~|| = {nrm(Tw/lam - V)/nrm(V):.2e}   (should be ~lambda)")

# ---- (2) large-eta: G^A -> 0  => T -> V~ ----
g_big = (1.0/Nk) / (omega0 - eact + 1j*(1.0e4/RY))
Tle = tmat(V, g_big)
print("\n--- validation (2) eta=1e4 eV (G^A->0):  T -> V~ ---")
print(f"  ||T - V~|| / ||V~|| = {nrm(Tle - V)/nrm(V):.2e}   (should be ~0)")

# ---- (4) resonance scan over the active window ----
print("\n--- resonance scan: ||T_PP(omega)|| / ||V~||  (peaks = active resonances) ---")
print("    omega[eV]   ||T_PP||/||V~||   ||T_M||/||M||")
for om in np.linspace(wmin, wmax, 9):
    g = gA(om)
    rpp = nrm(tmat(V, g))/nrm(V)
    rm  = nrm(tmat(M, g))/nrm(M)
    mark = "  <-- omega0" if abs(om-omega0) < (wmax-wmin)/18 else ""
    print(f"   {om*RY:9.3f}   {rpp:12.3f}   {rm:12.3f}{mark}")
