#!/bin/bash
# ===========================================================================
#  EDT P3 — full off-diagonal downfolded-potential block  Ṽ_{(m k_f),(n k_i)}
#  over the whole active manifold.  Pool-parallel over the k'-sum:
#  128 MPI ranks = 128 pools (1 rank/pool) -> full G-grid per rank, native
#  h_psi at each pool-local channel k'; blocks mp_sum'd over inter_pool_comm.
#
#  Memory-bound: each rank replicates the supercell KB projectors
#  vkb_d/vkb_p (nkb=1926/1944 x npwx=40779 ~2.5 GB), evc_act_all (~1 GB),
#  the per-channel source kets S (~1 GB), and V_d/V_p/V_colin (~0.4 GB)
#  => ~5.8 GB/rank.  A wholenode node is 257 GB, so ~40 ranks fit per node.
#  We run 36 ranks on ONE wholenode node (36 pools, 144 k / 36 = 4 channels
#  per pool, balanced; ~209 GB, safe).  A 1-node job backfills immediately
#  into idle nodes, whereas a 4-node/8h job stalls behind multi-node
#  reservations for higher-priority (throttled) jobs.  Runtime ~5 h.
#  (For the fastest run on a free queue use highmem -N 1 -n 128 -nk 128.)
# ===========================================================================
#SBATCH -J edt_block
#SBATCH -A che190065
#SBATCH -p wholenode
#SBATCH -N 1
#SBATCH --ntasks-per-node=36
#SBATCH -n 36
#SBATCH -t 07:00:00
#SBATCH -o edt_block.%j.out
#SBATCH -e edt_block.%j.err

module reset
module load aocc/3.1.0 openmpi/4.1.6
module load amdblis/3.0 amdlibflame/3.0 amdlibm/3.0 fftw

export OMP_NUM_THREADS=1
ulimit -s unlimited

RUNDIR=/anvil/scratch/x-rg47749/claude/sternheimerED/edt/run
EDTX=/anvil/scratch/x-rg47749/claude/sternheimerED/edt/src/edt.x
cd "$RUNDIR"

echo "================================================================"
echo " host   = $(hostname)"
echo " nodes  = $SLURM_JOB_NODELIST"
echo " ntasks = $SLURM_NTASKS"
echo " start  = $(date)"
echo "================================================================"

# 36 ranks, 36 pools (-nk 36): one k-point pool per rank, no PW split.
srun --mpi=pmi2 -n 36 "$EDTX" -nk 36 -inp edt_block.in

echo "================================================================"
echo " end    = $(date)  rc=$?"
echo "================================================================"
