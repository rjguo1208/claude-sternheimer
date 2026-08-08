MODULE edt_twolevel
  !---------------------------------------------------------------------------
  !  Two-level Sternheimer (Schur-partitioned / "deflated") rest dressing.
  !
  !    Sigma(w0) = V_A1 D1 V_1A  +  <s~_b| G22 |s~_a>
  !    D1  = (w0 - H1 - V11)^-1          exact (zheevd on the R1 block)
  !    s~a = P_T V (1 + D1 P1 V) |psi_a> dressed entrance channel
  !    G22 = [w0 - H2 - W22]^-1,  W22 = V22 + V21 D1 V12
  !    ladder: x^(n+1) = D2 [ s~ + W22 x^(n) ]
  !
  !  Discrete-grid units: raw matrix elements M (supercell-summed, Ry); every
  !  operator V insertion carries 1/N_k (H = eps + M/N_k).  Output Sg block is
  !  in the vtilde_block_mpi file convention (Sg_file = N_k * Sigma_phys), so
  !  the downstream python keeps using (M+Sg)*RY/N_k unchanged.
  !
  !  THIS VERSION: MODE A only (tail_split_band = NS in (0,nbnd)): the "model
  !  tail" T = bands NS+1..nbnd is explicit, so the entire ladder is dense
  !  algebra and the exact Schur answer (zgesv) is computed alongside — the
  !  end-to-end validation gate for D1 / sources / W22 / rung structure /
  !  unit bookkeeping, cross-checkable against the python reference and
  !  against the internally computed bare limit.  The physical-tail mode
  !  (tail_split_band=0, CG-based D2 with the 150-band lift) plugs into the
  !  same skeleton next.
  !---------------------------------------------------------------------------
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: vtilde_block_twolevel, vtilde_block_lanczos

CONTAINS

  SUBROUTINE vtilde_block_twolevel(xkcr, omega0_ry, win_min_ry, win_max_ry, &
                                   nbndskip_in, nsplit, nrung, edmat_file, outfile)
    USE io_global,        ONLY : ionode, stdout
    USE wvfct,            ONLY : nbnd, et
    USE klist,            ONLY : nkstot, nks
    USE constants,        ONLY : rytoev
    USE edt_sternheimer,  ONLY : edmat_fill_or_check
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xkcr(3,nkstot), omega0_ry, win_min_ry, win_max_ry
    INTEGER,  INTENT(IN) :: nbndskip_in, nsplit, nrung
    CHARACTER(LEN=*), INTENT(IN) :: edmat_file, outfile

    INTEGER :: kg, ib, n, NU, NA_, NR1, NT, iu, info, it, lwork, lrwork, liwork, i, b
    REAL(dp) :: w0, eR, xnk, dstep, dprev
    INTEGER,  ALLOCATABLE :: uk(:), ub(:), selA(:), sel1(:), selT(:), ipiv(:), iwork(:)
    REAL(dp), ALLOCATABLE :: et_all(:,:), epsU(:), mu(:), rwork(:), g2(:)
    COMPLEX(dp), ALLOCATABLE :: MU_(:,:), H1(:,:), work(:)
    COMPLEX(dp), ALLOCATABLE :: VAA(:,:), MA1(:,:), M1A(:,:), M11(:,:), M1T(:,:)
    COMPLEX(dp), ALLOCATABLE :: MTA(:,:), MT1(:,:), MTT(:,:)
    COMPLEX(dp), ALLOCATABLE :: cco(:,:), st(:,:), W22(:,:), Aex(:,:)
    COMPLEX(dp), ALLOCATABLE :: x(:,:), xnew(:,:), xex(:,:), tmpR(:,:), tmpT(:,:)
    COMPLEX(dp), ALLOCATABLE :: SgR1(:,:), Sgt(:,:), Sgex(:,:), Sgbare(:,:), Vout(:,:)
    COMPLEX(dp) :: cone, czero, vij

    cone = (1.0_dp,0.0_dp); czero = (0.0_dp,0.0_dp)
    w0 = omega0_ry; xnk = DBLE(nkstot)
    IF (nsplit == 0) THEN
       CALL twolevel_physical(xkcr, omega0_ry, win_min_ry, win_max_ry, &
                              nbndskip_in, nrung, edmat_file, outfile)
       RETURN
    ENDIF
    IF (nsplit < 0 .OR. nsplit >= nbnd) CALL errore('vtilde_block_twolevel', &
         'need tail_split_band in [0, nbnd): 0 = physical tail, >0 = model tail', 1)

    ! ---- eigenvalues (pool-collected) + union manifold over ALL bands ----
    ALLOCATE(et_all(nbnd,nkstot))
    CALL poolcollect(nbnd, nks, et, nkstot, et_all)
    NU = nbnd*nkstot
    ALLOCATE(uk(NU), ub(NU), epsU(NU))
    n = 0
    DO kg = 1, nkstot
       DO ib = 1, nbnd
          n = n+1; uk(n) = kg; ub(n) = ib; epsU(n) = et_all(ib,kg)
       ENDDO
    ENDDO

    ! ---- A / R1 / T split ----
    ALLOCATE(selA(NU), sel1(NU), selT(NU))
    NA_ = 0; NR1 = 0; NT = 0
    DO n = 1, NU
       eR = epsU(n)
       IF (ub(n) > nbndskip_in .AND. eR >= win_min_ry .AND. eR <= win_max_ry) THEN
          NA_ = NA_+1; selA(NA_) = n
       ELSEIF (ub(n) <= nsplit) THEN
          NR1 = NR1+1; sel1(NR1) = n
       ELSE
          NT = NT+1; selT(NT) = n
       ENDIF
    ENDDO
    IF (ionode) THEN
       WRITE(stdout,'(/,5X,A)') REPEAT('=',64)
       WRITE(stdout,'(5X,A)') 'TWO-LEVEL rest dressing — MODE A (model tail = validation gate)'
       WRITE(stdout,'(5X,A,3I8)')    'N_A, N_R1, N_T       = ', NA_, NR1, NT
       WRITE(stdout,'(5X,A,I4,A,I3)') 'split band = ', nsplit, '    ladder rungs = ', nrung
       WRITE(stdout,'(5X,A,F10.4,A)') 'omega0 = ', w0*rytoev, ' eV'
       FLUSH(stdout)
    ENDIF

    ! ---- union block from the validated edmat reader (ionode) ----
    ALLOCATE(MU_(NU,NU)); MU_ = czero
    IF (ionode) THEN
       CALL edmat_fill_or_check(edmat_file, xkcr, uk, ub, NU, MU_, .FALSE.)

       ALLOCATE(VAA(NA_,NA_), MA1(NA_,NR1), M1A(NR1,NA_), M11(NR1,NR1), M1T(NR1,NT))
       ALLOCATE(MTA(NT,NA_), MT1(NT,NR1), MTT(NT,NT))
       VAA = MU_(selA(1:NA_), selA(1:NA_))
       MA1 = MU_(selA(1:NA_), sel1(1:NR1));   M1A = MU_(sel1(1:NR1), selA(1:NA_))
       M11 = MU_(sel1(1:NR1), sel1(1:NR1));   M1T = MU_(sel1(1:NR1), selT(1:NT))
       MTA = MU_(selT(1:NT),  selA(1:NA_));   MT1 = MU_(selT(1:NT),  sel1(1:NR1))
       MTT = MU_(selT(1:NT),  selT(1:NT))
       DEALLOCATE(MU_)

       ! ---- bare reference (2nd order, both R1 and T through bare denominators) ----
       ALLOCATE(Sgbare(NA_,NA_), tmpR(NR1,NA_), tmpT(NT,NA_))
       DO n = 1, NR1
          tmpR(n,:) = M1A(n,:) / (w0 - epsU(sel1(n)))
       ENDDO
       CALL ZGEMM('N','N', NA_, NA_, NR1, cone/xnk, MA1, NA_, tmpR, NR1, czero, Sgbare, NA_)
       DO n = 1, NT
          tmpT(n,:) = MTA(n,:) / (w0 - epsU(selT(n)))
       ENDDO
       CALL ZGEMM('C','N', NA_, NA_, NT, cone/xnk, MTA, NT, tmpT, NT, cone, Sgbare, NA_)

       ! ---- D1: zheevd of H1 = diag(eps_R1) + M11/N_k ----
       ALLOCATE(H1(NR1,NR1), mu(NR1))
       H1 = M11 / xnk
       DO n = 1, NR1
          H1(n,n) = H1(n,n) + epsU(sel1(n))
       ENDDO
       H1 = (H1 + CONJG(TRANSPOSE(H1))) / 2.0_dp
       lwork = 2*NR1 + NR1*NR1 + 32
       lrwork = 1 + 5*NR1 + 2*NR1*NR1
       liwork = 3 + 5*NR1
       ALLOCATE(work(lwork), rwork(lrwork), iwork(liwork))
       CALL ZHEEVD('V', 'U', NR1, H1, NR1, mu, work, lwork, rwork, lrwork, iwork, liwork, info)
       IF (info /= 0) CALL errore('vtilde_block_twolevel', 'zheevd failed on H1', ABS(info))
       DEALLOCATE(work, rwork, iwork, M11)
       WRITE(stdout,'(5X,A,2F10.3,A,F8.3,A)') 'R1 dressed spectrum: ', &
            MINVAL(mu)*rytoev, MAXVAL(mu)*rytoev, ' eV;  min|w0-mu| = ', &
            MINVAL(ABS(w0-mu))*rytoev, ' eV'

       ! ---- c = D1 (M1A / N_k)  (dressed-source coefficients) ----
       ALLOCATE(cco(NR1,NA_))
       CALL ZGEMM('C','N', NR1, NA_, NR1, cone, H1, NR1, M1A, NR1, czero, tmpR, NR1)
       DO n = 1, NR1
          tmpR(n,:) = tmpR(n,:) / ((w0 - mu(n)) * xnk)
       ENDDO
       CALL ZGEMM('N','N', NR1, NA_, NR1, cone, H1, NR1, tmpR, NR1, czero, cco, NR1)

       ! ---- dressed entrance sources: st = MTA + MT1 * c   (raw) ----
       ALLOCATE(st(NT,NA_))
       st = MTA
       CALL ZGEMM('N','N', NT, NA_, NR1, cone, MT1, NT, cco, NR1, cone, st, NT)

       ! ---- W22 = MTT/N_k + (1/N_k^2) MT1 D1 M1T ----
       ALLOCATE(W22(NT,NT))
       BLOCK
         COMPLEX(dp), ALLOCATABLE :: B1(:,:)
         ALLOCATE(B1(NR1,NT))
         CALL ZGEMM('C','N', NR1, NT, NR1, cone, H1, NR1, M1T, NR1, czero, B1, NR1)
         DO n = 1, NR1
            B1(n,:) = B1(n,:) / ((w0 - mu(n)) * xnk * xnk)
         ENDDO
         CALL ZGEMM('N','N', NR1, NT, NR1, cone, H1, NR1, B1, NR1, czero, M1T, NR1)
         W22 = MTT / xnk
         CALL ZGEMM('N','N', NT, NT, NR1, cone, MT1, NT, M1T, NR1, cone, W22, NT)
         DEALLOCATE(B1)
       END BLOCK

       ! ---- D2 diag + the ladder ----
       ALLOCATE(g2(NT), x(NT,NA_), xnew(NT,NA_))
       DO n = 1, NT
          g2(n) = 1.0_dp / (w0 - epsU(selT(n)))
       ENDDO
       DO n = 1, NT
          x(n,:) = g2(n) * st(n,:)
       ENDDO
       dprev = -1.0_dp
       DO it = 1, nrung
          CALL ZGEMM('N','N', NT, NA_, NT, cone, W22, NT, x, NT, czero, xnew, NT)
          DO n = 1, NT
             xnew(n,:) = g2(n) * (st(n,:) + xnew(n,:))
          ENDDO
          dstep = SQRT(SUM(ABS(xnew - x)**2)) / MAX(SQRT(SUM(ABS(x)**2)), 1.d-30)
          WRITE(stdout,'(5X,A,I3,A,ES12.4,A,F8.4)') 'rung ', it, ':  |dx|/|x| = ', dstep, &
               '   ratio vs prev = ', MERGE(dstep/dprev, 0.0_dp, dprev > 0.0_dp)
          x = xnew
          dprev = dstep
       ENDDO

       ! ---- exact Schur reference: (w0 - eps_T - W22) xex = st ----
       ALLOCATE(Aex(NT,NT), xex(NT,NA_), ipiv(NT))
       Aex = -W22
       DO n = 1, NT
          Aex(n,n) = Aex(n,n) + (w0 - epsU(selT(n)))
       ENDDO
       xex = st
       CALL ZGESV(NT, NA_, Aex, NT, ipiv, xex, NT, info)
       IF (info /= 0) CALL errore('vtilde_block_twolevel', 'zgesv failed on the T block', ABS(info))

       ! ---- Sigma assembly (file convention: Sg_file = N_k * Sigma_phys) ----
       ALLOCATE(SgR1(NA_,NA_), Sgt(NA_,NA_), Sgex(NA_,NA_))
       CALL ZGEMM('N','N', NA_, NA_, NR1, cone, MA1, NA_, cco, NR1, czero, SgR1, NA_)
       CALL ZGEMM('C','N', NA_, NA_, NT, cone/xnk, st, NT, x,   NT, czero, Sgt,  NA_)
       CALL ZGEMM('C','N', NA_, NA_, NT, cone/xnk, st, NT, xex, NT, czero, Sgex, NA_)

       WRITE(stdout,'(5X,A)') REPEAT('-',64)
       WRITE(stdout,'(5X,A,ES12.4,A)') 'max|Sg_ladder - Sg_exact|   = ', &
            MAXVAL(ABS(Sgt - Sgex)), ' Ry  (ladder truncation)'
       WRITE(stdout,'(5X,A,3ES12.4)') 'max|SgR1|, |Sg_tail|, |Sg_bare| = ', &
            MAXVAL(ABS(SgR1)), MAXVAL(ABS(Sgex)), MAXVAL(ABS(Sgbare))
       WRITE(stdout,'(5X,A)') REPEAT('-',64)

       ! ---- write .dat (standard 7-record layout; Sg = dressed R1 + EXACT tail) ----
       ALLOCATE(Vout(NA_,NA_))
       Vout = VAA + SgR1 + Sgex
       DO n = 1, NA_
          DO b = n, NA_
             vij = 0.5_dp*(Vout(b,n) + CONJG(Vout(n,b)))
             Vout(b,n) = vij; Vout(n,b) = CONJG(vij)
          ENDDO
       ENDDO
       iu = 78
       OPEN(unit=iu, file=TRIM(outfile), form='unformatted', status='replace')
       WRITE(iu) NA_, nkstot, nkstot, nbndskip_in
       WRITE(iu) w0, win_min_ry, win_max_ry
       WRITE(iu) uk(selA(1:NA_)), ub(selA(1:NA_))
       WRITE(iu) epsU(selA(1:NA_))
       WRITE(iu) VAA
       WRITE(iu) SgR1 + Sgex
       WRITE(iu) Vout
       CLOSE(iu)
       WRITE(stdout,'(5X,3A,I6,A)') 'wrote ', TRIM(outfile), '  (A-block, N_A = ', NA_, ')'
       WRITE(stdout,'(5X,A)') REPEAT('=',64)
       FLUSH(stdout)
    ENDIF
    IF (ALLOCATED(MU_)) DEALLOCATE(MU_)
    DEALLOCATE(et_all, uk, ub, epsU, selA, sel1, selT)
  END SUBROUTINE vtilde_block_twolevel


  SUBROUTINE twolevel_physical(xkcr, omega0_ry, win_min_ry, win_max_ry, &
                               nbndskip_in, nrung, edmat_file, outfile)
    !-------------------------------------------------------------------------
    ! MODE B v2 (memory-lean): physical tail with local-channel slicing.
    !   st/rhs/x live only on the pool's channels; x is streamed per source-k
    !   (owner-pool bcast) during the fold sweeps; phi/xi are synthesized from
    !   coefficients on the fly (never materialized); sources st are kept for
    !   ALL sources locally (cheap when channel-sliced) so cross-batch Sigma
    !   blocks are exact; the ladder runs in source batches of NBB.
    ! Per rung two fold sweeps (beta -> gco -> xi dependency): sweep 1 gives
    ! (V x) on local channels (term1 + beta rows), sweep 2 gives V xi.
    !-------------------------------------------------------------------------
    USE io_global,        ONLY : ionode, ionode_id, stdout
    USE wvfct,            ONLY : nbnd, et, npwx
    USE klist,            ONLY : nkstot, nks, ngk, igk_k, xk
    USE fft_base,         ONLY : dffts
    USE fft_interfaces,   ONLY : invfft, fwfft
    USE noncollin_module, ONLY : npol
    USE pw_restart_new,   ONLY : read_collected_wfc
    USE io_files,         ONLY : restart_dir
    USE constants,        ONLY : rytoev
    USE gvect,            ONLY : g, ngm
    USE gvecw,            ONLY : gcutw
    USE edic_mod,         ONLY : V_d, V_p
    USE edt_source,       ONLY : build_V_folded, count_nkb, make_coeff
    USE edt_sternheimer,  ONLY : edmat_fill_or_check, hpsi_setup_k, solve_rest_cg
    USE mp,               ONLY : mp_sum, mp_bcast
    USE mp_pools,         ONLY : inter_pool_comm, my_pool_id, npool
    USE mp_world,         ONLY : world_comm
    USE becmod,           ONLY : becp, allocate_bec_type, deallocate_bec_type
    USE uspp,             ONLY : nkb
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xkcr(3,nkstot), omega0_ry, win_min_ry, win_max_ry
    INTEGER,  INTENT(IN) :: nbndskip_in, nrung
    CHARACTER(LEN=*), INTENT(IN) :: edmat_file, outfile

    INTEGER, PARAMETER :: NBB = 64
    INTEGER :: kg, ib, n, NU, NA_, NR1, iu, info, it, lwork, lrwork, liwork
    INTEGER :: ik, iks, a, b, i, a0, na_b, ibatch, nbatch, npw_kp, iters, kgs
    INTEGER :: nkb_d, nkb_p, ikb, owner, sweep
    REAL(dp) :: w0, eR, xnk, alpha_ry, resid
    INTEGER,  ALLOCATABLE :: uk(:), ub(:), selA(:), sel1(:), pos1(:), iwork(:)
    INTEGER,  ALLOCATABLE :: igk_all(:,:), ngk_all(:), kowner(:)
    REAL(dp), ALLOCATABLE :: et_all(:,:), epsU(:), mu(:), rwork(:), gk_tmp(:)
    COMPLEX(dp), ALLOCATABLE :: MU_(:,:), H1(:,:), work(:), VAA(:,:), MA1(:,:), M1A(:,:)
    COMPLEX(dp), ALLOCATABLE :: cco(:,:), SgR1(:,:), Sgt(:,:), Vout(:,:), tR(:,:)
    COMPLEX(dp), ALLOCATABLE :: evc_all(:,:,:), stt(:,:,:), x(:,:,:), rhs3(:,:,:)
    COMPLEX(dp), ALLOCATABLE :: acc(:,:,:), xbuf(:,:)
    COMPLEX(dp), ALLOCATABLE :: beta1(:,:), gco(:,:), Vfc(:,:), psic(:), rhs(:), chi(:)
    COMPLEX(dp), ALLOCATABLE :: vkb_d(:,:), vkb_p(:,:), bd(:), bp(:), tmpd(:), tmpp(:)
    COMPLEX(dp), ALLOCATABLE :: cdsum(:,:), cpsum(:,:), evtmp(:,:), proj(:), gam(:), sloc(:,:)
    COMPLEX(dp), ALLOCATABLE :: kcd(:,:), kcp(:,:)
    REAL(dp) :: xk3(3)
    COMPLEX(dp) :: cone, czero, vij
    INTEGER, EXTERNAL :: global_kpoint_index

    cone = (1.0_dp,0.0_dp); czero = (0.0_dp,0.0_dp)
    w0 = omega0_ry; xnk = DBLE(nkstot)
    IF (npol /= 1) CALL errore('twolevel_physical','npol=1 only',1)
    iks = global_kpoint_index(nkstot, 1)

    ! ---- eigenvalues, manifold, k ownership ----
    ALLOCATE(et_all(nbnd,nkstot))
    CALL poolcollect(nbnd, nks, et, nkstot, et_all)
    NU = nbnd*nkstot
    ALLOCATE(uk(NU), ub(NU), epsU(NU))
    n = 0
    DO kg = 1, nkstot
       DO ib = 1, nbnd
          n = n+1; uk(n) = kg; ub(n) = ib; epsU(n) = et_all(ib,kg)
       ENDDO
    ENDDO
    ALLOCATE(selA(NU), sel1(NU), pos1(NU)); pos1 = 0
    NA_ = 0; NR1 = 0
    DO n = 1, NU
       eR = epsU(n)
       IF (ub(n) > nbndskip_in .AND. eR >= win_min_ry .AND. eR <= win_max_ry) THEN
          NA_ = NA_+1; selA(NA_) = n
       ELSE
          NR1 = NR1+1; sel1(NR1) = n; pos1(n) = NR1
       ENDIF
    ENDDO
    ALLOCATE(kowner(nkstot)); kowner = 0
    DO ik = 1, nks
       kowner(ik + iks - 1) = my_pool_id
    ENDDO
    CALL mp_sum(kowner, inter_pool_comm)
    alpha_ry = 2.0_dp*ABS(w0 - MINVAL(et_all)) + 1.0_dp
    IF (ionode) THEN
       WRITE(stdout,'(/,5X,A)') REPEAT('=',64)
       WRITE(stdout,'(5X,A)') 'TWO-LEVEL rest dressing — MODE B v2 (physical tail, memory-lean)'
       WRITE(stdout,'(5X,A,2I8,A,I3,A,F8.3,A)') 'N_A, N_R1 = ', NA_, NR1, &
            '   rungs = ', nrung, '   alpha = ', alpha_ry*rytoev, ' eV'
       FLUSH(stdout)
    ENDIF

    ! ---- explicit blocks + D1 (ionode), bcast the source coefficients ----
    ALLOCATE(VAA(NA_,NA_), MA1(NA_,NR1), M1A(NR1,NA_), cco(NR1,NA_), SgR1(NA_,NA_))
    IF (ionode) THEN
       ALLOCATE(MU_(NU,NU)); MU_ = czero
       CALL edmat_fill_or_check(edmat_file, xkcr, uk, ub, NU, MU_, .FALSE.)
       VAA = MU_(selA(1:NA_), selA(1:NA_))
       MA1 = MU_(selA(1:NA_), sel1(1:NR1))
       M1A = MU_(sel1(1:NR1), selA(1:NA_))
       ALLOCATE(H1(NR1,NR1), mu(NR1))
       H1 = MU_(sel1(1:NR1), sel1(1:NR1)) / xnk
       DEALLOCATE(MU_)
       DO n = 1, NR1
          H1(n,n) = H1(n,n) + epsU(sel1(n))
       ENDDO
       H1 = (H1 + CONJG(TRANSPOSE(H1))) / 2.0_dp
       lwork = 2*NR1 + NR1*NR1 + 32; lrwork = 1 + 5*NR1 + 2*NR1*NR1; liwork = 3 + 5*NR1
       ALLOCATE(work(lwork), rwork(lrwork), iwork(liwork))
       CALL ZHEEVD('V','U', NR1, H1, NR1, mu, work, lwork, rwork, lrwork, iwork, liwork, info)
       IF (info /= 0) CALL errore('twolevel_physical','zheevd failed',ABS(info))
       DEALLOCATE(work, rwork, iwork)
       WRITE(stdout,'(5X,A,F8.3,A)') 'R1 dressed: min|w0-mu| = ', MINVAL(ABS(w0-mu))*rytoev, ' eV'
       ALLOCATE(tR(NR1,NA_))
       CALL ZGEMM('C','N', NR1, NA_, NR1, cone, H1, NR1, M1A, NR1, czero, tR, NR1)
       DO n = 1, NR1
          tR(n,:) = tR(n,:) / ((w0 - mu(n)) * xnk)
       ENDDO
       CALL ZGEMM('N','N', NR1, NA_, NR1, cone, H1, NR1, tR, NR1, czero, cco, NR1)
       DEALLOCATE(tR)
       CALL ZGEMM('N','N', NA_, NA_, NR1, cone, MA1, NA_, cco, NR1, czero, SgR1, NA_)
    ENDIF
    CALL mp_bcast(cco, ionode_id, world_comm)

    ! ---- replicate wavefunctions (v2 keeps this; v3 streams it too) ----
    ALLOCATE(evc_all(npwx,nbnd,nkstot)); evc_all = czero
    ALLOCATE(evtmp(npwx,nbnd))
    DO ik = 1, nks
       kg = ik + iks - 1
       CALL read_collected_wfc(restart_dir(), ik, evtmp)
       evc_all(:,:,kg) = evtmp
    ENDDO
    DO kg = 1, nkstot
       CALL mp_sum(evc_all(:,:,kg), inter_pool_comm)
    ENDDO
    DEALLOCATE(evtmp)
    ALLOCATE(igk_all(npwx,nkstot), ngk_all(nkstot), gk_tmp(npwx))
    DO kg = 1, nkstot
       xk3 = xk_dummy(kg, xkcr)
       CALL gk_sort(xk3, ngm, g, gcutw, ngk_all(kg), igk_all(1,kg), gk_tmp)
    ENDDO

    CALL count_nkb(V_d%nat, V_d%ityp, V_d%ntyp, nkb_d)
    CALL count_nkb(V_p%nat, V_p%ityp, V_p%ntyp, nkb_p)
    ALLOCATE(vkb_d(npwx,nkb_d), vkb_p(npwx,nkb_p), bd(nkb_d), bp(nkb_p))
    ALLOCATE(tmpd(nkb_d), tmpp(nkb_p))
    ALLOCATE(psic(dffts%nnr), rhs(npwx), chi(npwx), proj(nbnd), gam(nbnd))
    ALLOCATE(Vfc(dffts%nnr, nks))
    ALLOCATE(Sgt(NA_,NA_)); Sgt = czero

    ! ---- channel-local fold cache: Vfc(:,ik) for q = k_gs - k_ch is rebuilt per
    !      kgs inside the sweeps; here we only allocate the per-channel slot ----

    ! ================= sources: st(:, local ch, ALL sources) =================
    ALLOCATE(stt(npwx, nks, NA_)); stt = czero
    ALLOCATE(cdsum(nkb_d,NA_), cpsum(nkb_p,NA_)); cdsum = czero; cpsum = czero
    ALLOCATE(acc(dffts%nnr, nks, NBB))
    nbatch = (NA_ + NBB - 1) / NBB
    DO ibatch = 1, nbatch
       a0 = (ibatch-1)*NBB; na_b = MIN(NBB, NA_ - a0)
       acc(:,:,1:na_b) = czero
       DO kgs = 1, nkstot
          ! folds for this source-k against my channels
          DO ik = 1, nks
             CALL build_V_folded(xkcr(:,kgs) - xkcr(:,ik+iks-1), Vfc(:,ik))
          ENDDO
          ! synthesize phi^{kgs} per source; KB bec on the owner pool only
          IF (kowner(kgs) == my_pool_id) THEN
             CALL get_betavkb(ngk_all(kgs), igk_all(1,kgs), xk_dummy(kgs,xkcr), vkb_d, &
                              V_d%nat,V_d%ityp,V_d%tau, nkb_d)
             CALL get_betavkb(ngk_all(kgs), igk_all(1,kgs), xk_dummy(kgs,xkcr), vkb_p, &
                              V_p%nat,V_p%ityp,V_p%tau, nkb_p)
          ENDIF
          DO a = 1, na_b
             gam = czero
             DO ib = 1, nbnd
                n = (kgs-1)*nbnd + ib
                IF (pos1(n) > 0) gam(ib) = cco(pos1(n), a0+a)
             ENDDO
             IF (uk(selA(a0+a)) == kgs) gam(ub(selA(a0+a))) = gam(ub(selA(a0+a))) + cone
             CALL ZGEMV('N', npwx, nbnd, cone, evc_all(1,1,kgs), npwx, gam, 1, czero, rhs, 1)
             IF (kowner(kgs) == my_pool_id) THEN
                DO ikb = 1, nkb_d
                   bd(ikb) = SUM(CONJG(vkb_d(1:ngk_all(kgs),ikb))*rhs(1:ngk_all(kgs)))
                ENDDO
                DO ikb = 1, nkb_p
                   bp(ikb) = SUM(CONJG(vkb_p(1:ngk_all(kgs),ikb))*rhs(1:ngk_all(kgs)))
                ENDDO
                CALL make_coeff(V_d%nat,V_d%ityp,V_d%ntyp,nkb_d,bd,tmpd)
                CALL make_coeff(V_p%nat,V_p%ityp,V_p%ntyp,nkb_p,bp,tmpp)
                cdsum(:,a0+a) = cdsum(:,a0+a) + tmpd
                cpsum(:,a0+a) = cpsum(:,a0+a) + tmpp
             ENDIF
             psic = czero
             DO n = 1, ngk_all(kgs)
                psic(dffts%nl(igk_all(n,kgs))) = rhs(n)
             ENDDO
             CALL invfft('Wave', psic, dffts)
             DO ik = 1, nks
                acc(:,ik,a) = acc(:,ik,a) + Vfc(:,ik) * psic
             ENDDO
          ENDDO
       ENDDO
       ! close: fwfft per (channel, source), add KB, project P_T
       DO ik = 1, nks
          kg = ik + iks - 1
          DO a = 1, na_b
             psic = acc(:,ik,a)
             CALL fwfft('Wave', psic, dffts)
             DO n = 1, ngk_all(kg)
                stt(n,ik,a0+a) = psic(dffts%nl(igk_all(n,kg)))
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL mp_sum(cdsum, inter_pool_comm)
    CALL mp_sum(cpsum, inter_pool_comm)
    DO ik = 1, nks
       kg = ik + iks - 1
       CALL get_betavkb(ngk_all(kg), igk_all(1,kg), xk(1,ik), vkb_d, V_d%nat,V_d%ityp,V_d%tau, nkb_d)
       CALL get_betavkb(ngk_all(kg), igk_all(1,kg), xk(1,ik), vkb_p, V_p%nat,V_p%ityp,V_p%tau, nkb_p)
       CALL ZGEMM('N','N', ngk_all(kg), NA_, nkb_d,  cone, vkb_d, npwx, cdsum, nkb_d, cone, stt(1,ik,1), npwx*nks)
       CALL ZGEMM('N','N', ngk_all(kg), NA_, nkb_p, -cone, vkb_p, npwx, cpsum, nkb_p, cone, stt(1,ik,1), npwx*nks)
       DO a = 1, NA_
          CALL ZGEMV('C', ngk_all(kg), nbnd, cone, evc_all(1,1,kg), npwx, stt(1,ik,a), 1, czero, proj, 1)
          CALL ZGEMV('N', ngk_all(kg), nbnd, -cone, evc_all(1,1,kg), npwx, proj, 1, cone, stt(1,ik,a), 1)
       ENDDO
    ENDDO
    DEALLOCATE(cdsum, cpsum)
    IF (ionode) THEN
       WRITE(stdout,'(5X,A)') 'sources built (local-channel slices).'
       FLUSH(stdout)
    ENDIF

    ! ================= ladder (batched sources) =================
    ALLOCATE(x(npwx,nks,NBB), rhs3(npwx,nks,NBB), xbuf(npwx,NBB))
    ALLOCATE(beta1(NR1,NBB), gco(NR1,NBB), sloc(NA_,NBB))
    CALL allocate_bec_type(nkb, 1, becp)
    DO ibatch = 1, nbatch
       a0 = (ibatch-1)*NBB; na_b = MIN(NBB, NA_ - a0)
       x(:,:,1:na_b) = czero
       DO it = 0, nrung
          IF (it == 0) THEN
             rhs3(:,:,1:na_b) = stt(:,:,a0+1:a0+na_b)
          ELSE
             ! ---- sweep 1: (V x) on local channels -> term1 + beta rows ----
             acc(:,:,1:na_b) = czero
             DO kgs = 1, nkstot
                owner = kowner(kgs)
                IF (owner == my_pool_id) xbuf(:,1:na_b) = x(:, kgs-iks+1, 1:na_b)
                CALL mp_bcast(xbuf(:,1:na_b), owner, inter_pool_comm)
                DO ik = 1, nks
                   CALL build_V_folded(xkcr(:,kgs) - xkcr(:,ik+iks-1), Vfc(:,ik))
                ENDDO
                DO a = 1, na_b
                   psic = czero
                   DO n = 1, ngk_all(kgs)
                      psic(dffts%nl(igk_all(n,kgs))) = xbuf(n,a)
                   ENDDO
                   CALL invfft('Wave', psic, dffts)
                   DO ik = 1, nks
                      acc(:,ik,a) = acc(:,ik,a) + Vfc(:,ik) * psic
                   ENDDO
                ENDDO
             ENDDO
             ! KB for V x: accumulate coefficients like the source pass
             ! (recomputed per rung; small)
             CALL vx_kb_coeffs(x, na_b)
             beta1 = czero
             rhs3(:,:,1:na_b) = czero
             DO ik = 1, nks
                kg = ik + iks - 1
                DO a = 1, na_b
                   psic = acc(:,ik,a)
                   CALL fwfft('Wave', psic, dffts)
                   DO n = 1, ngk_all(kg)
                      rhs3(n,ik,a) = psic(dffts%nl(igk_all(n,kg)))
                   ENDDO
                ENDDO
             ENDDO
             CALL add_kb_local(rhs3, na_b)
             DO ik = 1, nks
                kg = ik + iks - 1
                DO a = 1, na_b
                   CALL ZGEMV('C', ngk_all(kg), nbnd, cone, evc_all(1,1,kg), npwx, rhs3(1,ik,a), 1, czero, proj, 1)
                   DO ib = 1, nbnd
                      n = (kg-1)*nbnd + ib
                      IF (pos1(n) > 0) beta1(pos1(n),a) = proj(ib)
                   ENDDO
                ENDDO
             ENDDO
             CALL mp_sum(beta1, inter_pool_comm)
             IF (ionode) THEN
                ALLOCATE(tR(NR1,na_b))
                CALL ZGEMM('C','N', NR1, na_b, NR1, cone, H1, NR1, beta1, NR1, czero, tR, NR1)
                DO n = 1, NR1
                   tR(n,:) = tR(n,:) / ((w0 - mu(n)) * xnk * xnk)
                ENDDO
                CALL ZGEMM('N','N', NR1, na_b, NR1, cone, H1, NR1, tR, NR1, czero, gco(:,1:na_b), NR1)
                DEALLOCATE(tR)
             ENDIF
             CALL mp_bcast(gco, ionode_id, world_comm)
             ! ---- sweep 2: V xi (xi synthesized from gco) ----
             rhs3(:,:,1:na_b) = rhs3(:,:,1:na_b)/xnk + stt(:,:,a0+1:a0+na_b)
             CALL sweep_synth(gco, na_b, rhs3)
             ! P_T projection of the rhs
             DO ik = 1, nks
                kg = ik + iks - 1
                DO a = 1, na_b
                   CALL ZGEMV('C', ngk_all(kg), nbnd, cone, evc_all(1,1,kg), npwx, rhs3(1,ik,a), 1, czero, proj, 1)
                   CALL ZGEMV('N', ngk_all(kg), nbnd, -cone, evc_all(1,1,kg), npwx, proj, 1, cone, rhs3(1,ik,a), 1)
                ENDDO
             ENDDO
          ENDIF
          ! ---- D2: pool-local CG;  x = -chi ----
          DO ik = 1, nks
             kg = ik + iks - 1
             CALL hpsi_setup_k(ik)
             DO a = 1, na_b
                rhs = czero; rhs(1:ngk_all(kg)) = rhs3(1:ngk_all(kg),ik,a)
                CALL solve_rest_cg(ik, rhs, w0, alpha_ry, nbnd, evc_all(:,:,kg), 1.d-8, 500, chi, iters, resid)
                x(1:ngk_all(kg),ik,a) = -chi(1:ngk_all(kg))
             ENDDO
          ENDDO
          IF (ionode) THEN
             WRITE(stdout,'(5X,A,I3,A,I3,A)') 'batch ', ibatch, '  rung ', it, ' done'
             FLUSH(stdout)
          ENDIF
       ENDDO
       ! ---- Sigma columns: <st_b | x_a>, all b, batch a ----
       sloc(:,1:na_b) = czero
       DO ik = 1, nks
          kg = ik + iks - 1
          CALL ZGEMM('C','N', NA_, na_b, ngk_all(kg), cone, stt(1,ik,1), npwx*nks, &
                     x(1,ik,1), npwx*nks, cone, sloc, NA_)
       ENDDO
       CALL mp_sum(sloc(:,1:na_b), inter_pool_comm)
       IF (ionode) Sgt(:, a0+1:a0+na_b) = sloc(:,1:na_b)/xnk
    ENDDO
    CALL deallocate_bec_type(becp)

    ! ---- assemble + write (ionode) ----
    IF (ionode) THEN
       ALLOCATE(Vout(NA_,NA_))
       Vout = VAA + SgR1 + Sgt
       DO n = 1, NA_
          DO b = n, NA_
             vij = 0.5_dp*(Vout(b,n) + CONJG(Vout(n,b)))
             Vout(b,n) = vij; Vout(n,b) = CONJG(vij)
          ENDDO
       ENDDO
       iu = 78
       OPEN(unit=iu, file=TRIM(outfile), form='unformatted', status='replace')
       WRITE(iu) NA_, nkstot, nkstot, nbndskip_in
       WRITE(iu) w0, win_min_ry, win_max_ry
       WRITE(iu) uk(selA(1:NA_)), ub(selA(1:NA_))
       WRITE(iu) epsU(selA(1:NA_))
       WRITE(iu) VAA
       WRITE(iu) SgR1 + Sgt
       WRITE(iu) Vout
       CLOSE(iu)
       WRITE(stdout,'(5X,3A,I6,A)') 'wrote ', TRIM(outfile), '  (MODE B v2, N_A = ', NA_, ')'
       WRITE(stdout,'(5X,A)') REPEAT('=',64)
       FLUSH(stdout)
    ENDIF

  CONTAINS

    FUNCTION xk_dummy(kg_, xkcr_) RESULT(xkc3)
      USE cell_base, ONLY : bg
      INTEGER,  INTENT(IN) :: kg_
      REAL(dp), INTENT(IN) :: xkcr_(3,nkstot)
      REAL(dp) :: xkc3(3)
      xkc3 = MATMUL(bg, xkcr_(:,kg_))
    END FUNCTION xk_dummy

    SUBROUTINE vx_kb_coeffs(xin, nsrc)
      !! KB coefficients of (V_NL x): bec at each local k of x, make_coeff, sum;
      !! stored in module-batch arrays kcd/kcp for add_kb_local.
      COMPLEX(dp), INTENT(IN) :: xin(:,:,:)
      INTEGER, INTENT(IN) :: nsrc
      INTEGER :: ikl, kgl, aa, ikb2
      IF (.NOT. ALLOCATED(kcd)) ALLOCATE(kcd(nkb_d,NBB), kcp(nkb_p,NBB))
      kcd = czero; kcp = czero
      DO ikl = 1, nks
         kgl = ikl + iks - 1
         CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xk(1,ikl), vkb_d, V_d%nat,V_d%ityp,V_d%tau, nkb_d)
         CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xk(1,ikl), vkb_p, V_p%nat,V_p%ityp,V_p%tau, nkb_p)
         DO aa = 1, nsrc
            DO ikb2 = 1, nkb_d
               bd(ikb2) = SUM(CONJG(vkb_d(1:ngk_all(kgl),ikb2))*xin(1:ngk_all(kgl),ikl,aa))
            ENDDO
            DO ikb2 = 1, nkb_p
               bp(ikb2) = SUM(CONJG(vkb_p(1:ngk_all(kgl),ikb2))*xin(1:ngk_all(kgl),ikl,aa))
            ENDDO
            CALL make_coeff(V_d%nat,V_d%ityp,V_d%ntyp,nkb_d,bd,tmpd)
            CALL make_coeff(V_p%nat,V_p%ityp,V_p%ntyp,nkb_p,bp,tmpp)
            kcd(:,aa) = kcd(:,aa) + tmpd
            kcp(:,aa) = kcp(:,aa) + tmpp
         ENDDO
      ENDDO
      CALL mp_sum(kcd, inter_pool_comm)
      CALL mp_sum(kcp, inter_pool_comm)
    END SUBROUTINE vx_kb_coeffs

    SUBROUTINE add_kb_local(v3, nsrc)
      !! add the KB nonlocal part (coefficients kcd/kcp) at local channels
      COMPLEX(dp), INTENT(INOUT) :: v3(:,:,:)
      INTEGER, INTENT(IN) :: nsrc
      INTEGER :: ikl, kgl
      DO ikl = 1, nks
         kgl = ikl + iks - 1
         CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xk(1,ikl), vkb_d, V_d%nat,V_d%ityp,V_d%tau, nkb_d)
         CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xk(1,ikl), vkb_p, V_p%nat,V_p%ityp,V_p%tau, nkb_p)
         CALL ZGEMM('N','N', ngk_all(kgl), nsrc, nkb_d,  cone, vkb_d, npwx, kcd, nkb_d, cone, v3(1,ikl,1), SIZE(v3,1)*SIZE(v3,2))
         CALL ZGEMM('N','N', ngk_all(kgl), nsrc, nkb_p, -cone, vkb_p, npwx, kcp, nkb_p, cone, v3(1,ikl,1), SIZE(v3,1)*SIZE(v3,2))
      ENDDO
    END SUBROUTINE add_kb_local

    SUBROUTINE sweep_synth(gc, nsrc, v3out)
      !! v3out += V xi at local channels, xi synthesized per (kgs,a) from gc
      !! (KB part included via owner-pool becs, added locally afterwards).
      COMPLEX(dp), INTENT(IN)    :: gc(:,:)
      INTEGER,     INTENT(IN)    :: nsrc
      COMPLEX(dp), INTENT(INOUT) :: v3out(:,:,:)
      INTEGER :: kgs2, ikl, aa, ib2, n2, ikb2, kgl
      acc(:,:,1:nsrc) = czero
      IF (.NOT. ALLOCATED(kcd)) ALLOCATE(kcd(nkb_d,NBB), kcp(nkb_p,NBB))
      kcd = czero; kcp = czero
      DO kgs2 = 1, nkstot
         DO ikl = 1, nks
            CALL build_V_folded(xkcr(:,kgs2) - xkcr(:,ikl+iks-1), Vfc(:,ikl))
         ENDDO
         IF (kowner(kgs2) == my_pool_id) THEN
            CALL get_betavkb(ngk_all(kgs2), igk_all(1,kgs2), xk_dummy(kgs2,xkcr), vkb_d, &
                             V_d%nat,V_d%ityp,V_d%tau, nkb_d)
            CALL get_betavkb(ngk_all(kgs2), igk_all(1,kgs2), xk_dummy(kgs2,xkcr), vkb_p, &
                             V_p%nat,V_p%ityp,V_p%tau, nkb_p)
         ENDIF
         DO aa = 1, nsrc
            gam = czero
            DO ib2 = 1, nbnd
               n2 = (kgs2-1)*nbnd + ib2
               IF (pos1(n2) > 0) gam(ib2) = gc(pos1(n2), aa)
            ENDDO
            CALL ZGEMV('N', npwx, nbnd, cone, evc_all(1,1,kgs2), npwx, gam, 1, czero, rhs, 1)
            IF (kowner(kgs2) == my_pool_id) THEN
               DO ikb2 = 1, nkb_d
                  bd(ikb2) = SUM(CONJG(vkb_d(1:ngk_all(kgs2),ikb2))*rhs(1:ngk_all(kgs2)))
               ENDDO
               DO ikb2 = 1, nkb_p
                  bp(ikb2) = SUM(CONJG(vkb_p(1:ngk_all(kgs2),ikb2))*rhs(1:ngk_all(kgs2)))
               ENDDO
               CALL make_coeff(V_d%nat,V_d%ityp,V_d%ntyp,nkb_d,bd,tmpd)
               CALL make_coeff(V_p%nat,V_p%ityp,V_p%ntyp,nkb_p,bp,tmpp)
               kcd(:,aa) = kcd(:,aa) + tmpd
               kcp(:,aa) = kcp(:,aa) + tmpp
            ENDIF
            psic = czero
            DO n2 = 1, ngk_all(kgs2)
               psic(dffts%nl(igk_all(n2,kgs2))) = rhs(n2)
            ENDDO
            CALL invfft('Wave', psic, dffts)
            DO ikl = 1, nks
               acc(:,ikl,aa) = acc(:,ikl,aa) + Vfc(:,ikl) * psic
            ENDDO
         ENDDO
      ENDDO
      CALL mp_sum(kcd, inter_pool_comm)
      CALL mp_sum(kcp, inter_pool_comm)
      DO ikl = 1, nks
         kgl = ikl + iks - 1
         DO aa = 1, nsrc
            psic = acc(:,ikl,aa)
            CALL fwfft('Wave', psic, dffts)
            DO n2 = 1, ngk_all(kgl)
               v3out(n2,ikl,aa) = v3out(n2,ikl,aa) + psic(dffts%nl(igk_all(n2,kgl)))
            ENDDO
         ENDDO
         CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xk(1,ikl), vkb_d, V_d%nat,V_d%ityp,V_d%tau, nkb_d)
         CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xk(1,ikl), vkb_p, V_p%nat,V_p%ityp,V_p%tau, nkb_p)
         CALL ZGEMM('N','N', ngk_all(kgl), nsrc, nkb_d,  cone, vkb_d, npwx, kcd, nkb_d, cone, v3out(1,ikl,1), SIZE(v3out,1)*SIZE(v3out,2))
         CALL ZGEMM('N','N', ngk_all(kgl), nsrc, nkb_p, -cone, vkb_p, npwx, kcp, nkb_p, cone, v3out(1,ikl,1), SIZE(v3out,1)*SIZE(v3out,2))
      ENDDO
    END SUBROUTINE sweep_synth

  END SUBROUTINE twolevel_physical


  SUBROUTINE vtilde_block_lanczos(xkcr, omega0_ry, win_min_ry, win_max_ry, &
                                  nbndskip_in, nstep, edmat_file, outfile)
    !-------------------------------------------------------------------------
    ! MODE C: omega-RESOLVED rest self-energy via GLOBAL BLOCK LANCZOS on the
    ! full rest space,
    !     Sigma(w)_ab = <chi_b| [w - P_R H P_R]^-1 |chi_a>,
    !     chi_a = P_R V_raw |psi_a>,   P_R = 1 - P_A  (A = active manifold),
    !     H = H0 (h_psi) + V_raw/N_k   (discrete-grid convention, Ry).
    ! One chain serves every omega: after nstep block steps the file carries
    ! (R0, A_j, B_j); python evaluates the block continued fraction
    !     S_N = (w-A_N)^-1,  S_j = [w - A_j - B_j^dag S_{j+1} B_j]^-1,
    !     Sigma_file(w) = R0^dag S_1(w) R0 / N_k
    ! (same Sg_file scale as MODE A/B: H_eff = eps + (VAA + Sg)*RY/N_k).
    ! Recurrence convention:  H Q_j = Q_{j-1} B_{j-1}^dag + Q_j A_j + Q_{j+1} B_j
    ! with B_j upper triangular from the block QR of the residual.
    ! Efficiency: folds cached once per (local channel, source-k); KB bec sums
    ! as ZGEMM; full re-orth (CGS x2) as ONE strided ZGEMM against all stored
    ! blocks; h_psi batched over the whole 396-wide block per channel.
    !-------------------------------------------------------------------------
    USE io_global,        ONLY : ionode, ionode_id, stdout
    USE wvfct,            ONLY : nbnd, et, npwx
    USE klist,            ONLY : nkstot, nks, ngk, igk_k, xk
    USE fft_base,         ONLY : dffts
    USE fft_interfaces,   ONLY : invfft, fwfft
    USE noncollin_module, ONLY : npol
    USE pw_restart_new,   ONLY : read_collected_wfc
    USE io_files,         ONLY : restart_dir
    USE constants,        ONLY : rytoev
    USE gvect,            ONLY : g, ngm
    USE gvecw,            ONLY : gcutw
    USE edic_mod,         ONLY : V_d, V_p
    USE edt_source,       ONLY : build_V_folded, count_nkb, make_coeff
    USE edt_sternheimer,  ONLY : edmat_fill_or_check, hpsi_setup_k
    USE mp,               ONLY : mp_sum, mp_bcast, mp_alltoall
    USE mp_pools,         ONLY : inter_pool_comm, my_pool_id, npool
    USE mp_world,         ONLY : world_comm
    USE becmod,           ONLY : becp, allocate_bec_type, deallocate_bec_type
    USE uspp,             ONLY : nkb
    USE constants,        ONLY : tpi
    USE edt_input,        ONLY : fold_col, col_chunk, zslab_tol
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xkcr(3,nkstot), omega0_ry, win_min_ry, win_max_ry
    INTEGER,  INTENT(IN) :: nbndskip_in, nstep
    CHARACTER(LEN=*), INTENT(IN) :: edmat_file, outfile

    INTEGER :: kg, ib, n, NU, NA_, iu, info, ik, iks, a, i, j, jj, npw
    INTEGER :: nkb_d, nkb_p, owner, nact_max, nstep_done, istat
    REAL(dp) :: w0, eR, xnk, t0, t1, tfold, thpsi, tortho, herm
    INTEGER,  ALLOCATABLE :: uk(:), ub(:), selA(:), nact_k(:), actb(:,:)
    INTEGER,  ALLOCATABLE :: igk_all(:,:), ngk_all(:), kowner(:)
    REAL(dp), ALLOCATABLE :: et_all(:,:), epsU(:), gk_tmp(:)
    COMPLEX(dp), ALLOCATABLE :: MU_(:,:), VAA(:,:), R0(:,:), Aj(:,:), Bj(:,:), prevB(:,:)
    COMPLEX(dp), ALLOCATABLE :: Ast(:,:,:), Bst(:,:,:), Gm(:,:), Cm(:,:), projm(:,:)
    COMPLEX(dp), ALLOCATABLE :: Qs(:,:,:,:), W3(:,:,:), xbuf(:,:), acc3(:,:,:)
    COMPLEX(dp), ALLOCATABLE :: Vfa(:,:,:), psic(:), evtmp(:,:), actA(:,:,:)
    COMPLEX(dp), ALLOCATABLE :: vkb_d(:,:), vkb_p(:,:), kcd(:,:), kcp(:,:)
    COMPLEX(dp), ALLOCATABLE :: BDl(:,:), BPl(:,:), tmpd(:), tmpp(:), psi_c(:,:), hps_c(:,:)
    COMPLEX(dp), ALLOCATABLE :: evc_all_l(:,:,:)
    ! ---- column-distributed fold (fold_col): each rank owns nloc columns and
    !      ALL k-channels, so every FFT is done once instead of npool times ----
    COMPLEX(dp), ALLOCATABLE :: Vq(:,:), phw(:,:), xcol(:,:,:), acol(:,:,:)
    COMPLEX(dp), ALLOCATABLE :: sbuf(:,:,:), rbuf(:,:,:), psia(:,:,:), xcol_out(:,:,:)
    INTEGER,     ALLOCATABLE :: map_iq(:,:), map_ip(:,:)
    INTEGER :: nloc, col0, iq_, ip_, ia, ga, ncol, ib0c, nbc
    LOGICAL :: use_col
    ! ---- 2D slab: the vacuum carries no wavefunction, so the r-space contraction
    !      can skip those z-slices entirely (the FFTs still need the full grid) ----
    INTEGER, PARAMETER :: CBLKZ = 1024
    INTEGER, ALLOCATABLE :: blk_lo(:), blk_ce(:)
    INTEGER :: nblk, nxy_, nzkeep, nzseg, zseg(2,2)
    REAL(dp) :: zdrop
    REAL(dp) :: xk3(3)
    COMPLEX(dp) :: cone, czero
    INTEGER, EXTERNAL :: global_kpoint_index
    EXTERNAL :: h_psi

    cone = (1.0_dp,0.0_dp); czero = (0.0_dp,0.0_dp)
    w0 = omega0_ry; xnk = DBLE(nkstot)
    IF (npol /= 1) CALL errore('vtilde_block_lanczos','npol=1 only',1)
    iks = global_kpoint_index(nkstot, 1)

    ! ---- eigenvalues, active manifold, k ownership ----
    ALLOCATE(et_all(nbnd,nkstot))
    CALL poolcollect(nbnd, nks, et, nkstot, et_all)
    NU = nbnd*nkstot
    ALLOCATE(uk(NU), ub(NU), epsU(NU))
    n = 0
    DO kg = 1, nkstot
       DO ib = 1, nbnd
          n = n+1; uk(n) = kg; ub(n) = ib; epsU(n) = et_all(ib,kg)
       ENDDO
    ENDDO
    ALLOCATE(selA(NU))
    NA_ = 0
    DO n = 1, NU
       eR = epsU(n)
       IF (ub(n) > nbndskip_in .AND. eR >= win_min_ry .AND. eR <= win_max_ry) THEN
          NA_ = NA_+1; selA(NA_) = n
       ENDIF
    ENDDO
    ALLOCATE(nact_k(nkstot)); nact_k = 0
    DO a = 1, NA_
       nact_k(uk(selA(a))) = nact_k(uk(selA(a))) + 1
    ENDDO
    nact_max = MAXVAL(nact_k)
    ALLOCATE(actb(nact_max, nkstot)); actb = 0; nact_k = 0
    DO a = 1, NA_
       kg = uk(selA(a)); nact_k(kg) = nact_k(kg) + 1
       actb(nact_k(kg), kg) = ub(selA(a))
    ENDDO
    ALLOCATE(kowner(nkstot)); kowner = 0
    DO ik = 1, nks
       kowner(ik + iks - 1) = my_pool_id
    ENDDO
    CALL mp_sum(kowner, inter_pool_comm)
    IF (ionode) THEN
       WRITE(stdout,'(/,5X,A)') REPEAT('=',64)
       WRITE(stdout,'(5X,A)') 'OMEGA-RESOLVED rest — MODE C (global block Lanczos, full rest)'
       WRITE(stdout,'(5X,A,I8,A,I4)') 'N_A (block width) = ', NA_, '   block steps = ', nstep
       FLUSH(stdout)
    ENDIF

    ! ---- VAA from the all-band block file (ionode); MU_ kept for the
    !      operator unit test ----
    ALLOCATE(VAA(NA_,NA_))
    IF (ionode) THEN
       ALLOCATE(MU_(NU,NU)); MU_ = czero
       CALL edmat_fill_or_check(edmat_file, xkcr, uk, ub, NU, MU_, .FALSE.)
       VAA = MU_(selA(1:NA_), selA(1:NA_))
    ENDIF

    ! ---- replicate wavefunctions; igk tables; active-column buffers ----
    ALLOCATE(evc_all_l(npwx,nbnd,nkstot), STAT=istat)
    IF (istat /= 0) CALL errore('vtilde_block_lanczos','evc_all alloc failed',1)
    evc_all_l = czero
    ALLOCATE(evtmp(npwx,nbnd))
    DO ik = 1, nks
       kg = ik + iks - 1
       CALL read_collected_wfc(restart_dir(), ik, evtmp)
       evc_all_l(:,:,kg) = evtmp
    ENDDO
    DO kg = 1, nkstot
       CALL mp_sum(evc_all_l(:,:,kg), inter_pool_comm)
    ENDDO
    DEALLOCATE(evtmp)
    ALLOCATE(igk_all(npwx,nkstot), ngk_all(nkstot), gk_tmp(npwx))
    DO kg = 1, nkstot
       xk3 = xkc3_of(kg)
       CALL gk_sort(xk3, ngm, g, gcutw, ngk_all(kg), igk_all(1,kg), gk_tmp)
    ENDDO
    ALLOCATE(actA(npwx, nact_max, nks))
    DO ik = 1, nks
       kg = ik + iks - 1
       actA(:,:,ik) = czero
       DO i = 1, nact_k(kg)
          actA(:, i, ik) = evc_all_l(:, actb(i,kg), kg)
       ENDDO
    ENDDO

    CALL count_nkb(V_d%nat, V_d%ityp, V_d%ntyp, nkb_d)
    CALL count_nkb(V_p%nat, V_p%ityp, V_p%ntyp, nkb_p)
    ALLOCATE(vkb_d(npwx,nkb_d), vkb_p(npwx,nkb_p), tmpd(nkb_d), tmpp(nkb_p))
    ALLOCATE(kcd(nkb_d,NA_), kcp(nkb_p,NA_), BDl(nkb_d,NA_), BPl(nkb_p,NA_))
    ALLOCATE(psic(dffts%nnr))

    ! ---- fold cache: ALL (source-k, local channel) folds built once ----
    use_col = fold_col .AND. (nks*npool == nkstot) .AND. (MOD(NA_,npool) == 0)
    IF (fold_col .AND. .NOT.use_col .AND. ionode) WRITE(stdout,'(5X,A)') &
         'fold_col requested but layout unsuitable (need nks*npool = nkstot and NA_ divisible by npool) — using the k-layout'
    IF (.NOT.use_col) THEN
       ALLOCATE(Vfa(dffts%nnr, nks, nkstot), STAT=istat)
       IF (istat /= 0) CALL errore('vtilde_block_lanczos','fold cache alloc failed (reduce npool?)',1)
       DO kg = 1, nkstot
          DO ik = 1, nks
             CALL build_V_folded(xkcr(:,kg) - xkcr(:,ik+iks-1), Vfa(:,ik,kg))
          ENDDO
       ENDDO
    ELSE
       CALL build_qcanon()
       nloc = NA_/npool; col0 = my_pool_id*nloc
       ALLOCATE(xcol(npwx,nloc,nkstot), STAT=istat)
       IF (istat /= 0) CALL errore('vtilde_block_lanczos','column-layout alloc failed',1)
       ALLOCATE(sbuf(npwx*nks,nloc,npool), rbuf(npwx*nks,nloc,npool), xcol_out(npwx,nloc,nkstot))
       ncol = nloc
       IF (col_chunk > 0) ncol = MIN(col_chunk, nloc)
       IF (ncol > 64) CALL errore('vtilde_block_lanczos','fold_col: col_chunk>64 (raise the accb bound)',1)
       ALLOCATE(psia(dffts%nnr,nkstot,ncol), acol(dffts%nnr,nkstot,ncol), STAT=istat)
       IF (istat /= 0) CALL errore('vtilde_block_lanczos','psia/acol alloc failed',1)
       IF (ionode) WRITE(stdout,'(5X,A,I4,A,I4,A,I5,A,F7.2,A)') &
            'fold_col ON: ', nloc, ' columns/rank (chunk ', ncol, '), all ', nkstot, &
            ' channels local; fold buffers ', &
            2.d0*DBLE(dffts%nnr)*DBLE(nkstot)*DBLE(ncol)*16.d0/1.d9, ' GB/rank'
       CALL pick_zslab()
    ENDIF

    ! ---- Krylov storage ----
    ALLOCATE(Qs(npwx,nks,NA_,0:nstep), STAT=istat)
    IF (istat /= 0) CALL errore('vtilde_block_lanczos','Qs alloc failed (lower n_lancz)',1)
    ALLOCATE(W3(npwx,nks,NA_), xbuf(npwx,NA_))
    IF (.NOT.use_col) THEN
       ALLOCATE(acc3(dffts%nnr,nks,NA_))          ! k-layout accumulator (3 GB here)
    ELSE
       ALLOCATE(acc3(1,1,1))                      ! unused in the column layout
    ENDIF

    ! ---- OPERATOR UNIT TEST: apply_dV on known band states must reproduce
    !      the edmat columns (assumption-free discriminator vs MODE A/B input) ----
    CALL operator_unit_test()
    IF (ionode) DEALLOCATE(MU_)
    ALLOCATE(Gm(NA_,NA_), Aj(NA_,NA_), Bj(NA_,NA_), prevB(NA_,NA_), R0(NA_,NA_))
    ALLOCATE(projm(nact_max,NA_))
    IF (ionode) THEN
       ALLOCATE(Ast(NA_,NA_,nstep), Bst(NA_,NA_,nstep))
    ENDIF
    CALL allocate_bec_type(nkb, NA_, becp)
    ALLOCATE(psi_c(npwx,NA_), hps_c(npwx,NA_))

    ! ================= sources: chi_a = P_A-projected V_raw |psi_a> ==========
    W3 = czero
    DO a = 1, NA_
       kg = uk(selA(a))
       IF (kowner(kg) == my_pool_id) W3(:, kg-iks+1, a) = evc_all_l(:, ub(selA(a)), kg)
    ENDDO
    CALL apply_dV(W3, .TRUE.)     ! W3 <- V_raw * (planted bundle); raw scale
    CALL project_PA(W3)
    CALL block_gram(W3, W3, Gm)
    Gm = 0.5_dp*(Gm + CONJG(TRANSPOSE(Gm)))
    CALL ZPOTRF('U', NA_, Gm, NA_, info)
    IF (info /= 0) CALL errore('vtilde_block_lanczos','source QR breakdown',ABS(info))
    R0 = czero
    DO i = 1, NA_
       R0(1:i,i) = Gm(1:i,i)
    ENDDO
    CALL apply_trsm(W3, Gm)
    Qs(:,:,:,0) = W3
    IF (ionode) THEN
       WRITE(stdout,'(5X,A,ES12.4)') 'sources built; max diag R0 (Ry) = ', MAXVAL(ABS(R0))
       FLUSH(stdout)
    ENDIF

    ! ================= block Lanczos =================
    prevB = czero; nstep_done = nstep
    DO j = 0, nstep-1
       t0 = wtime()
       ! W = H_phys Q_j
       W3 = Qs(:,:,:,j)
       CALL apply_dV(W3, .FALSE.)          ! W3 <- V_raw Q_j (bundle sweep)
       W3 = W3 / xnk
       t1 = wtime(); tfold = t1-t0; t0 = t1
       DO ik = 1, nks
          kg = ik + iks - 1
          npw = ngk_all(kg)
          CALL hpsi_setup_k(ik)
          psi_c = czero
          psi_c(:,1:NA_) = Qs(:,ik,1:NA_,j)
          hps_c = czero
          CALL h_psi(npwx, npw, NA_, psi_c, hps_c)
          W3(1:npw,ik,1:NA_) = W3(1:npw,ik,1:NA_) + hps_c(1:npw,1:NA_)
       ENDDO
       CALL project_PA(W3)
       t1 = wtime(); thpsi = t1-t0; t0 = t1
       ! A_j
       CALL block_gram(Qs(:,:,:,j), W3, Aj)
       herm = MAXVAL(ABS(Aj - CONJG(TRANSPOSE(Aj))))
       Aj = 0.5_dp*(Aj + CONJG(TRANSPOSE(Aj)))
       ! residual: W - Q_j A_j - Q_{j-1} B_{j-1}^dag
       CALL gemm_sub(W3, Qs(:,:,:,j), Aj, 'N')
       IF (j > 0) CALL gemm_sub(W3, Qs(:,:,:,j-1), prevB, 'C')
       ! full re-orthogonalization, CGS x2, one strided GEMM per pass
       CALL reorth_all(W3, j)
       CALL reorth_all(W3, j)
       t1 = wtime(); tortho = t1-t0
       ! QR -> B_j, Q_{j+1}
       CALL block_gram(W3, W3, Gm)
       Gm = 0.5_dp*(Gm + CONJG(TRANSPOSE(Gm)))
       CALL ZPOTRF('U', NA_, Gm, NA_, info)
       IF (info /= 0) THEN
          nstep_done = j+1
          IF (ionode) WRITE(stdout,'(5X,A,I4)') 'invariant subspace at step ', j+1
          IF (ionode) Ast(:,:,j+1) = Aj
          IF (ionode) Bst(:,:,j+1) = czero
          EXIT
       ENDIF
       Bj = czero
       DO i = 1, NA_
          Bj(1:i,i) = Gm(1:i,i)
       ENDDO
       CALL apply_trsm(W3, Gm)
       Qs(:,:,:,j+1) = W3
       prevB = Bj
       IF (ionode) THEN
          Ast(:,:,j+1) = Aj; Bst(:,:,j+1) = Bj
          WRITE(stdout,'(5X,A,I4,A,ES9.2,A,3F8.1,A)') 'step ', j+1, &
               '  herm(A) = ', herm, '   t(fold,hpsi,orth) = ', tfold, thpsi, tortho, ' s'
          FLUSH(stdout)
       ENDIF
    ENDDO
    CALL deallocate_bec_type(becp)

    ! ---- write chains (ionode) ----
    IF (ionode) THEN
       iu = 79
       OPEN(unit=iu, file=TRIM(outfile), form='unformatted', status='replace')
       WRITE(iu) NA_, nkstot, nbndskip_in, nstep_done
       WRITE(iu) w0, win_min_ry, win_max_ry
       WRITE(iu) uk(selA(1:NA_)), ub(selA(1:NA_))
       WRITE(iu) epsU(selA(1:NA_))
       WRITE(iu) VAA
       WRITE(iu) R0
       WRITE(iu) Ast(:,:,1:nstep_done)
       WRITE(iu) Bst(:,:,1:nstep_done)
       CLOSE(iu)
       WRITE(stdout,'(5X,3A,I4,A)') 'wrote ', TRIM(outfile), '  (MODE C chains, nstep = ', nstep_done, ')'
       WRITE(stdout,'(5X,A)') REPEAT('=',64)
       FLUSH(stdout)
    ENDIF

  CONTAINS

    SUBROUTINE operator_unit_test()
      !! plant psi_(n,kt) for test bands, apply_dV (raw V, no projector),
      !! project on ALL (m,kg): must equal the edmat columns M(m kg, n kt).
      COMPLEX(dp), ALLOCATABLE :: tb(:,:,:), PRJ(:,:,:)
      INTEGER, PARAMETER :: NTB = 4
      INTEGER :: tbands(NTB), tt, kt, ikl, kgl, mrow
      REAL(dp) :: dmax, vmax
      ALLOCATE(tb(npwx,nks,NA_), PRJ(nbnd,nkstot,NTB))
      tbands(1) = MIN(5,nbnd); tbands(2) = MIN(nbndskip_in+12, nbnd)
      tbands(3) = MIN(60, nbnd); tbands(4) = nbnd
      kt = MIN(29, nkstot)
      tb = czero
      DO tt = 1, NTB
         IF (kowner(kt) == my_pool_id) tb(:, kt-iks+1, tt) = evc_all_l(:, tbands(tt), kt)
      ENDDO
      CALL apply_dV(tb, .FALSE.)
      PRJ = czero
      DO ikl = 1, nks
         kgl = ikl + iks - 1
         DO tt = 1, NTB
            CALL ZGEMV('C', ngk_all(kgl), nbnd, cone, evc_all_l(1,1,kgl), npwx, &
                       tb(1,ikl,tt), 1, czero, PRJ(1,kgl,tt), 1)
         ENDDO
      ENDDO
      CALL mp_sum(PRJ, inter_pool_comm)
      IF (ionode) THEN
         dmax = 0.d0; vmax = 0.d0
         DO tt = 1, NTB
            DO kgl = 1, nkstot
               DO mrow = 1, nbnd
                  dmax = MAX(dmax, ABS(PRJ(mrow,kgl,tt) - MU_((kgl-1)*nbnd+mrow, (kt-1)*nbnd+tbands(tt))))
                  vmax = MAX(vmax, ABS(MU_((kgl-1)*nbnd+mrow, (kt-1)*nbnd+tbands(tt))))
               ENDDO
            ENDDO
         ENDDO
         WRITE(stdout,'(5X,A,ES11.3,A,ES11.3,A,4I5)') 'OPERATOR UNIT TEST: max|proj - M| = ', dmax, &
              '   max|M| = ', vmax, '   bands:', tbands
         FLUSH(stdout)
      ENDIF
      DEALLOCATE(tb, PRJ)
    END SUBROUTINE operator_unit_test

    FUNCTION xkc3_of(kg_) RESULT(xkc3)
      USE cell_base, ONLY : bg
      INTEGER, INTENT(IN) :: kg_
      REAL(dp) :: xkc3(3)
      xkc3 = MATMUL(bg, xkcr(:,kg_))
    END FUNCTION xkc3_of

    FUNCTION wtime() RESULT(t_)
      REAL(dp) :: t_
      INTEGER(KIND=8) :: cnt, rate
      CALL SYSTEM_CLOCK(cnt, rate)
      t_ = DBLE(cnt)/DBLE(MAX(rate,1_8))
    END FUNCTION wtime

    SUBROUTINE apply_dV(v3, src_pass)
      !! v3 <- (V_local-fold + V_KB) v3   (raw scale), bundle over all k.
      !! src_pass: skip exactly-zero columns per source-k (planted bundles).
      COMPLEX(dp), INTENT(INOUT) :: v3(:,:,:)
      LOGICAL, INTENT(IN) :: src_pass
      INTEGER :: kgs, ikl, kgl, aa, n2, ig
      REAL(dp) :: tsc, tft, tml, tcl, tz, tzr, tbc, tkb
      INTEGER, SAVE :: ndv = 0
      tsc = 0.d0; tft = 0.d0; tml = 0.d0; tcl = 0.d0; tbc = 0.d0; tkb = 0.d0
      tz = wtime()
      acc3 = czero
      kcd = czero; kcp = czero
      tzr = wtime()-tz
      ! ---- KB coefficients: every rank does its OWN sources, in parallel, BEFORE
      !      the broadcast loop.  (Leaving this inside the loop serialised the
      !      whole sweep: 35 ranks idled at each bcast while one owner worked.)
      tz = wtime()
      DO ikl = 1, nks
         kgl = ikl + iks - 1
         CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xkc3_of(kgl), vkb_d, &
                          V_d%nat,V_d%ityp,V_d%tau, nkb_d)
         CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xkc3_of(kgl), vkb_p, &
                          V_p%nat,V_p%ityp,V_p%tau, nkb_p)
         CALL ZGEMM('C','N', nkb_d, NA_, ngk_all(kgl), cone, vkb_d, npwx, &
                    v3(1,ikl,1), SIZE(v3,1)*SIZE(v3,2), czero, BDl, nkb_d)
         CALL ZGEMM('C','N', nkb_p, NA_, ngk_all(kgl), cone, vkb_p, npwx, &
                    v3(1,ikl,1), SIZE(v3,1)*SIZE(v3,2), czero, BPl, nkb_p)
         DO aa = 1, NA_
            IF (src_pass .AND. uk(selA(aa)) /= kgl) CYCLE
            CALL make_coeff(V_d%nat,V_d%ityp,V_d%ntyp,nkb_d,BDl(:,aa),tmpd)
            CALL make_coeff(V_p%nat,V_p%ityp,V_p%ntyp,nkb_p,BPl(:,aa),tmpp)
            kcd(:,aa) = kcd(:,aa) + tmpd
            kcp(:,aa) = kcp(:,aa) + tmpp
         ENDDO
      ENDDO
      tkb = wtime()-tz
      CALL mp_sum(kcd, inter_pool_comm)
      CALL mp_sum(kcp, inter_pool_comm)

      IF (use_col) THEN
         tz = wtime()
         CALL k2col(v3, xcol)
         tbc = wtime()-tz
         CALL do_fold_col(src_pass, tsc, tft, tml)
         tz = wtime()
         CALL col2k(xcol, v3)
         tbc = tbc + (wtime()-tz); tz = wtime()
         DO ikl = 1, nks
            kgl = ikl + iks - 1
            CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xk(1,ikl), vkb_d, V_d%nat,V_d%ityp,V_d%tau, nkb_d)
            CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xk(1,ikl), vkb_p, V_p%nat,V_p%ityp,V_p%tau, nkb_p)
            CALL ZGEMM('N','N', ngk_all(kgl), NA_, nkb_d,  cone, vkb_d, npwx, kcd, nkb_d, cone, v3(1,ikl,1), SIZE(v3,1)*SIZE(v3,2))
            CALL ZGEMM('N','N', ngk_all(kgl), NA_, nkb_p, -cone, vkb_p, npwx, kcp, nkb_p, cone, v3(1,ikl,1), SIZE(v3,1)*SIZE(v3,2))
         ENDDO
         tcl = wtime()-tz
         ndv = ndv + 1
         IF (ionode .AND. ndv <= 3) THEN
            WRITE(stdout,'(5X,A,I2,A,5(A,F6.1))') 'apply_dV[col] #',ndv,' [s]:', &
                 ' acc0 ', tzr, ' transp ', tbc, ' KB ', tkb, &
                 ' fft ', tft, ' V*psi ', tml
            FLUSH(stdout)
         ENDIF
         RETURN
      ENDIF

      DO kgs = 1, nkstot
         owner = kowner(kgs)
         tz = wtime()
         IF (owner == my_pool_id) xbuf(:,1:NA_) = v3(:, kgs-iks+1, 1:NA_)
         CALL mp_bcast(xbuf(:,1:NA_), owner, inter_pool_comm)
         tbc = tbc + (wtime()-tz)
         DO aa = 1, NA_
            IF (src_pass .AND. uk(selA(aa)) /= kgs) CYCLE
            tz = wtime()
!$omp parallel do simd schedule(static)
            DO ig = 1, dffts%nnr
               psic(ig) = czero
            ENDDO
!$omp end parallel do simd
            DO n2 = 1, ngk_all(kgs)
               psic(dffts%nl(igk_all(n2,kgs))) = xbuf(n2,aa)
            ENDDO
            tsc = tsc + (wtime()-tz); tz = wtime()
            CALL invfft('Wave', psic, dffts)
            tft = tft + (wtime()-tz); tz = wtime()
            DO ikl = 1, nks
!$omp parallel do simd schedule(static)
               DO ig = 1, dffts%nnr
                  acc3(ig,ikl,aa) = acc3(ig,ikl,aa) + Vfa(ig,ikl,kgs) * psic(ig)
               ENDDO
!$omp end parallel do simd
            ENDDO
            tml = tml + (wtime()-tz)
         ENDDO
      ENDDO
      tz = wtime()
      v3 = czero
      DO ikl = 1, nks
         kgl = ikl + iks - 1
         DO aa = 1, NA_
            psic = acc3(:,ikl,aa)
            CALL fwfft('Wave', psic, dffts)
            DO n2 = 1, ngk_all(kgl)
               v3(n2,ikl,aa) = psic(dffts%nl(igk_all(n2,kgl)))
            ENDDO
         ENDDO
         CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xk(1,ikl), vkb_d, V_d%nat,V_d%ityp,V_d%tau, nkb_d)
         CALL get_betavkb(ngk_all(kgl), igk_all(1,kgl), xk(1,ikl), vkb_p, V_p%nat,V_p%ityp,V_p%tau, nkb_p)
         CALL ZGEMM('N','N', ngk_all(kgl), NA_, nkb_d,  cone, vkb_d, npwx, kcd, nkb_d, cone, v3(1,ikl,1), SIZE(v3,1)*SIZE(v3,2))
         CALL ZGEMM('N','N', ngk_all(kgl), NA_, nkb_p, -cone, vkb_p, npwx, kcp, nkb_p, cone, v3(1,ikl,1), SIZE(v3,1)*SIZE(v3,2))
      ENDDO
      tcl = wtime()-tz
      ndv = ndv + 1
      IF (ionode .AND. ndv <= 3) THEN
         WRITE(stdout,'(5X,A,I2,A,7(A,F6.1))') 'apply_dV #',ndv,' [s]:', &
              ' acc3zero ', tzr, ' bcast ', tbc, ' KBown ', tkb, ' scatter ', tsc, &
              ' invfft ', tft, ' V*psi ', tml, ' close ', tcl
         FLUSH(stdout)
      ENDIF
    END SUBROUTINE apply_dV

    SUBROUTINE build_qcanon()
      !! 36 canonical folds V_q (q on the coarse grid, wrapped to [0,1)) plus the
      !! four wrap phases exp(2*pi*i*G.r), G in {0,-1}^2.  Exact identity:
      !!   V_{q_w+G}(r) = exp(2*pi*i*G.r) V_{q_w}(r)   (G.R integer kills the cell sum)
      INTEGER :: kq, ik2, kg2, g1, g2, ig1, ig2, i1, i2, i3, ii, ip
      REAL(dp) :: qw(3), qe(3), a1, a2
      COMPLEX(dp) :: p1, p2
      ALLOCATE(Vq(dffts%nnr, nkstot), phw(dffts%nnr, 4))
      ALLOCATE(map_iq(nkstot,nkstot), map_ip(nkstot,nkstot))
      DO kq = 1, nkstot
         qw = xkcr(:,kq) - xkcr(:,1)
         qw = qw - FLOOR(qw + 1.d-8)
         CALL build_V_folded(qw, Vq(:,kq))
      ENDDO
      DO kg2 = 1, nkstot
         DO ik2 = 1, nkstot
            qe = xkcr(:,kg2) - xkcr(:,ik2)
            qw = qe - FLOOR(qe + 1.d-8)
            map_iq(ik2,kg2) = 0
            DO kq = 1, nkstot
               a1 = xkcr(1,kq)-xkcr(1,1) - FLOOR(xkcr(1,kq)-xkcr(1,1)+1.d-8)
               a2 = xkcr(2,kq)-xkcr(2,1) - FLOOR(xkcr(2,kq)-xkcr(2,1)+1.d-8)
               IF (ABS(a1-qw(1)) < 1.d-6 .AND. ABS(a2-qw(2)) < 1.d-6) map_iq(ik2,kg2) = kq
            ENDDO
            IF (map_iq(ik2,kg2) == 0) CALL errore('build_qcanon','q not found on the grid',1)
            g1 = NINT(qe(1)-qw(1)); g2 = NINT(qe(2)-qw(2))
            IF (g1 < -1 .OR. g1 > 0 .OR. g2 < -1 .OR. g2 > 0) &
                 CALL errore('build_qcanon','unexpected wrap vector',1)
            map_ip(ik2,kg2) = 1 + (-g1) + 2*(-g2)
         ENDDO
      ENDDO
      DO ip = 1, 4
         ig1 = -MOD(ip-1,2); ig2 = -((ip-1)/2)
         ii = 0
         DO i3 = 0, dffts%nr3-1
            DO i2 = 0, dffts%nr2-1
               p2 = CMPLX(COS(tpi*ig2*DBLE(i2)/dffts%nr2), SIN(tpi*ig2*DBLE(i2)/dffts%nr2), dp)
               DO i1 = 0, dffts%nr1-1
                  p1 = CMPLX(COS(tpi*ig1*DBLE(i1)/dffts%nr1), SIN(tpi*ig1*DBLE(i1)/dffts%nr1), dp)
                  ii = ii + 1
                  IF (ii <= dffts%nnr) phw(ii,ip) = p1*p2
               ENDDO
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE build_qcanon

    SUBROUTINE k2col(xk_, xc_)
      !! (npwx, nks, NA_) [my channels, all columns] -> (npwx, nloc, nkstot)
      !! [my columns, all channels].  Pool p owns the contiguous channel block
      !! (p-1)*nks+1 .. p*nks, so the alltoall slot index carries the channel.
      COMPLEX(dp), INTENT(IN)  :: xk_(:,:,:)
      COMPLEX(dp), INTENT(OUT) :: xc_(:,:,:)
      INTEGER :: p, il, ic, o
      DO p = 1, npool
         DO il = 1, nks
            o = (il-1)*npwx
            DO ic = 1, nloc
               sbuf(o+1:o+npwx, ic, p) = xk_(:, il, (p-1)*nloc+ic)
            ENDDO
         ENDDO
      ENDDO
      CALL mp_alltoall(sbuf, rbuf, inter_pool_comm)
      DO p = 1, npool
         DO il = 1, nks
            o = (il-1)*npwx
            DO ic = 1, nloc
               xc_(:, ic, (p-1)*nks+il) = rbuf(o+1:o+npwx, ic, p)
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE k2col

    SUBROUTINE col2k(xc_, xk_)
      COMPLEX(dp), INTENT(IN)  :: xc_(:,:,:)
      COMPLEX(dp), INTENT(OUT) :: xk_(:,:,:)
      INTEGER :: p, il, ic, o
      DO p = 1, npool
         DO il = 1, nks
            o = (il-1)*npwx
            DO ic = 1, nloc
               sbuf(o+1:o+npwx, ic, p) = xc_(:, ic, (p-1)*nks+il)
            ENDDO
         ENDDO
      ENDDO
      CALL mp_alltoall(sbuf, rbuf, inter_pool_comm)
      xk_ = czero
      DO p = 1, npool
         DO il = 1, nks
            o = (il-1)*npwx
            DO ic = 1, nloc
               xk_(:, il, (p-1)*nloc+ic) = rbuf(o+1:o+npwx, ic, p)
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE col2k

    SUBROUTINE pick_zslab()
      !! 2D systems: most of the cell is vacuum, where the wavefunctions have no
      !! weight, so the r-space contraction V(r)psi(r) there is multiplying zeros.
      !! Measure the z-profile of the active density, keep the shortest cyclic
      !! window holding 1-zslab_tol of it, and hand the contraction a block list
      !! covering only that window.  The FFTs still run on the full grid; only the
      !! multiply is restricted, and the multiply is 5/6 of the fold.
      INTEGER  :: iz, ib, kg, n2, lo, ln, i0, i1, best_lo, best_len, cnt, ibk, c0
      REAL(dp) :: tot, acc_, want, got
      REAL(dp), ALLOCATABLE :: rz(:)
      nxy_ = dffts%nnr / MAX(dffts%nr3,1)
      nzkeep = dffts%nr3
      zdrop  = 0.0_dp
      best_lo = 1; best_len = dffts%nr3
      IF (zslab_tol > 0.0_dp .AND. nxy_*dffts%nr3 == dffts%nnr) THEN
         ALLOCATE(rz(dffts%nr3)); rz = 0.0_dp
         DO kg = iks, iks+nks-1                      ! this pool's own channels
            DO ib = 1, nbnd
               psic = czero
               DO n2 = 1, ngk_all(kg)
                  psic(dffts%nl(igk_all(n2,kg))) = evc_all_l(n2,ib,kg)
               ENDDO
               CALL invfft('Wave', psic, dffts)
               DO iz = 1, dffts%nr3
                  rz(iz) = rz(iz) + SUM(ABS(psic((iz-1)*nxy_+1:iz*nxy_))**2)
               ENDDO
            ENDDO
         ENDDO
         CALL mp_sum(rz, inter_pool_comm)
         tot = SUM(rz); want = (1.0_dp - zslab_tol)*tot
         best_len = dffts%nr3 + 1
         DO lo = 1, dffts%nr3                        ! shortest cyclic window >= want
            acc_ = 0.0_dp; cnt = 0
            DO ln = 1, dffts%nr3
               iz = MOD(lo+ln-2, dffts%nr3) + 1
               acc_ = acc_ + rz(iz); cnt = ln
               IF (acc_ >= want) EXIT
            ENDDO
            IF (acc_ >= want .AND. cnt < best_len) THEN
               best_len = cnt; best_lo = lo; got = acc_
            ENDIF
         ENDDO
         IF (best_len > dffts%nr3) THEN
            best_lo = 1; best_len = dffts%nr3; got = tot
         ENDIF
         nzkeep = best_len
         zdrop  = 1.0_dp - got/MAX(tot,1.d-300)
         DEALLOCATE(rz)
      ENDIF
      ! flat index segments (z is the slowest index, so each is contiguous)
      i0 = best_lo; i1 = best_lo + best_len - 1
      IF (i1 <= dffts%nr3) THEN
         nzseg = 1
         zseg(1,1) = (i0-1)*nxy_ + 1;  zseg(2,1) = i1*nxy_
      ELSE
         nzseg = 2
         zseg(1,1) = (i0-1)*nxy_ + 1;  zseg(2,1) = dffts%nr3*nxy_
         zseg(1,2) = 1;                zseg(2,2) = (i1 - dffts%nr3)*nxy_
      ENDIF
      nblk = 0
      DO ibk = 1, nzseg
         DO c0 = zseg(1,ibk), zseg(2,ibk), CBLKZ
            nblk = nblk + 1
         ENDDO
      ENDDO
      ALLOCATE(blk_lo(nblk), blk_ce(nblk))
      nblk = 0
      DO ibk = 1, nzseg
         DO c0 = zseg(1,ibk), zseg(2,ibk), CBLKZ
            nblk = nblk + 1
            blk_lo(nblk) = c0
            blk_ce(nblk) = MIN(CBLKZ, zseg(2,ibk)-c0+1)
         ENDDO
      ENDDO
      IF (ionode .AND. nzkeep < dffts%nr3) THEN
         WRITE(stdout,'(5X,A,I4,A,I4,A,F5.1,A,ES9.2,A,I3,A)') &
              'zslab: keeping ', nzkeep, ' of ', dffts%nr3, ' z-slices (', &
              100.d0*DBLE(nzkeep)/DBLE(dffts%nr3), '% of the grid), dropped weight ', &
              zdrop, ' in ', nzseg, ' segment(s)'
         FLUSH(stdout)
      ENDIF
    END SUBROUTINE pick_zslab

    SUBROUTINE do_fold_col(src_pass, tsc_, tft_, tml_)
      !! local-potential fold, entirely local: 396/npool*nkstot FFTs per rank
      LOGICAL, INTENT(IN) :: src_pass
      REAL(dp), INTENT(INOUT) :: tsc_, tft_, tml_
      INTEGER :: kgs2, ik2, n2, iq2, ip2, ig, c0, ce, ibk
      COMPLEX(dp) :: accb(CBLKZ, 64), vlo(CBLKZ)
      REAL(dp) :: tz2
      ! Column batches: psia/acol hold only `ncol` columns at a time, which is what
      ! makes large k-grids fit (the buffers are 2*nnr*nkstot*ncol complex).
      xcol_out = czero
      DO ib0c = 1, nloc, ncol
         nbc = MIN(ncol, nloc-ib0c+1)
         ! (a) FFT every (source, column-in-batch) ONCE into psia
         psia = czero
         DO kgs2 = 1, nkstot
            DO ia = 1, nbc
               ga = col0 + ib0c + ia - 1
               IF (src_pass .AND. uk(selA(ga)) /= kgs2) CYCLE
               tz2 = wtime()
!$omp parallel do simd schedule(static)
               DO ig = 1, dffts%nnr
                  psic(ig) = czero
               ENDDO
!$omp end parallel do simd
               DO n2 = 1, ngk_all(kgs2)
                  psic(dffts%nl(igk_all(n2,kgs2))) = xcol(n2,ib0c+ia-1,kgs2)
               ENDDO
               tsc_ = tsc_ + (wtime()-tz2); tz2 = wtime()
               CALL invfft('Wave', psic, dffts)
               psia(:,kgs2,ia) = psic
               tft_ = tft_ + (wtime()-tz2)
            ENDDO
         ENDDO
         ! (b) grid-blocked contraction (everything for this batch stays in cache)
         tz2 = wtime()
         IF (nzkeep < dffts%nr3) acol = czero    ! skipped slices must stay zero
!$omp parallel do private(ibk,c0,ce,ik2,kgs2,ia,iq2,ip2,ig,accb,vlo) schedule(static)
         DO ibk = 1, nblk
            c0 = blk_lo(ibk); ce = blk_ce(ibk)
            DO ik2 = 1, nkstot
               accb(1:ce,1:nbc) = czero
               DO kgs2 = 1, nkstot
                  iq2 = map_iq(ik2,kgs2); ip2 = map_ip(ik2,kgs2)
                  DO ig = 1, ce
                     vlo(ig) = phw(c0+ig-1,ip2)*Vq(c0+ig-1,iq2)
                  ENDDO
                  DO ia = 1, nbc
                     DO ig = 1, ce
                        accb(ig,ia) = accb(ig,ia) + vlo(ig)*psia(c0+ig-1,kgs2,ia)
                     ENDDO
                  ENDDO
               ENDDO
               DO ia = 1, nbc
                  DO ig = 1, ce
                     acol(c0+ig-1,ik2,ia) = accb(ig,ia)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
!$omp end parallel do
         tml_ = tml_ + (wtime()-tz2)
         ! (c) back to G space for this batch
         DO ik2 = 1, nkstot
            DO ia = 1, nbc
               psic = acol(:,ik2,ia)
               CALL fwfft('Wave', psic, dffts)
               DO n2 = 1, ngk_all(ik2)
                  xcol_out(n2,ib0c+ia-1,ik2) = psic(dffts%nl(igk_all(n2,ik2)))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      xcol = xcol_out
    END SUBROUTINE do_fold_col

    SUBROUTINE project_PA(v3)
      !! v3 <- (1 - P_A) v3 per local channel (subtract active bands only)
      COMPLEX(dp), INTENT(INOUT) :: v3(:,:,:)
      INTEGER :: ikl, kgl, na_k
      DO ikl = 1, nks
         kgl = ikl + iks - 1
         na_k = nact_k(kgl)
         IF (na_k == 0) CYCLE
         CALL ZGEMM('C','N', na_k, NA_, ngk_all(kgl),  cone, actA(1,1,ikl), npwx, &
                    v3(1,ikl,1), SIZE(v3,1)*SIZE(v3,2), czero, projm, nact_max)
         CALL ZGEMM('N','N', ngk_all(kgl), NA_, na_k, -cone, actA(1,1,ikl), npwx, &
                    projm, nact_max, cone, v3(1,ikl,1), SIZE(v3,1)*SIZE(v3,2))
      ENDDO
    END SUBROUTINE project_PA

    SUBROUTINE block_gram(u3, v3, G_)
      !! G_ = u3^dag v3 summed over local channels + pools (npwx rows, zero-padded)
      COMPLEX(dp), INTENT(IN)  :: u3(:,:,:), v3(:,:,:)
      COMPLEX(dp), INTENT(OUT) :: G_(:,:)
      INTEGER :: ikl
      G_ = czero
      DO ikl = 1, nks
         CALL ZGEMM('C','N', NA_, NA_, npwx, cone, u3(1,ikl,1), SIZE(u3,1)*SIZE(u3,2), &
                    v3(1,ikl,1), SIZE(v3,1)*SIZE(v3,2), cone, G_, NA_)
      ENDDO
      CALL mp_sum(G_, inter_pool_comm)
    END SUBROUTINE block_gram

    SUBROUTINE gemm_sub(w_, q_, m_, tr)
      !! w_ <- w_ - q_ * op(m_),  op = 'N' or 'C' (conjg-transpose)
      COMPLEX(dp), INTENT(INOUT) :: w_(:,:,:)
      COMPLEX(dp), INTENT(IN)    :: q_(:,:,:), m_(:,:)
      CHARACTER(LEN=1), INTENT(IN) :: tr
      INTEGER :: ikl
      DO ikl = 1, nks
         CALL ZGEMM('N',tr, npwx, NA_, NA_, -cone, q_(1,ikl,1), SIZE(q_,1)*SIZE(q_,2), &
                    m_, NA_, cone, w_(1,ikl,1), SIZE(w_,1)*SIZE(w_,2))
      ENDDO
    END SUBROUTINE gemm_sub

    SUBROUTINE reorth_all(w_, jmax)
      !! classical Gram-Schmidt of w_ against ALL stored blocks 0..jmax
      !! (single strided ZGEMM: Qs columns (a,j) share a uniform stride)
      COMPLEX(dp), INTENT(INOUT) :: w_(:,:,:)
      INTEGER, INTENT(IN) :: jmax
      INTEGER :: ikl, ncols
      ncols = NA_*(jmax+1)
      IF (.NOT. ALLOCATED(Cm)) ALLOCATE(Cm(NA_*(nstep+1), NA_))
      Cm(1:ncols,:) = czero
      DO ikl = 1, nks
         CALL ZGEMM('C','N', ncols, NA_, npwx, cone, Qs(1,ikl,1,0), SIZE(Qs,1)*SIZE(Qs,2), &
                    w_(1,ikl,1), SIZE(w_,1)*SIZE(w_,2), cone, Cm, SIZE(Cm,1))
      ENDDO
      CALL mp_sum(Cm(1:ncols,:), inter_pool_comm)
      DO ikl = 1, nks
         CALL ZGEMM('N','N', npwx, NA_, ncols, -cone, Qs(1,ikl,1,0), SIZE(Qs,1)*SIZE(Qs,2), &
                    Cm, SIZE(Cm,1), cone, w_(1,ikl,1), SIZE(w_,1)*SIZE(w_,2))
      ENDDO
    END SUBROUTINE reorth_all

    SUBROUTINE apply_trsm(w_, R_)
      !! w_ <- w_ * R_^-1  (R_ upper triangular, from ZPOTRF 'U')
      COMPLEX(dp), INTENT(INOUT) :: w_(:,:,:)
      COMPLEX(dp), INTENT(IN)    :: R_(:,:)
      INTEGER :: ikl
      DO ikl = 1, nks
         CALL ZTRSM('R','U','N','N', npwx, NA_, cone, R_, NA_, w_(1,ikl,1), SIZE(w_,1)*SIZE(w_,2))
      ENDDO
    END SUBROUTINE apply_trsm

  END SUBROUTINE vtilde_block_lanczos

END MODULE edt_twolevel
