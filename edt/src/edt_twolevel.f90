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
  PUBLIC :: vtilde_block_twolevel

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
    ! MODE B: physical tail.  T = plane-wave complement of ALL nbnd explicit
    ! bands.  Same two-level algebra as MODE A, but T-space objects live as
    ! PW vectors: D2 = Sternheimer CG with the alpha-lift on all nbnd bands,
    ! V applied via fold-multiply FFTs (+ KB nonlocal), P_T via wavefunction
    ! projection.  v1 (6x6-sized): full evc replicated on every pool
    ! (nbnd*nk*npwx complex), folds computed per (k,channel) on the fly,
    ! D1 dense algebra centralized on ionode (only small coefficient vectors
    ! cross ranks).  Multi-k vectors are replicated via mp_sum after each
    ! pool-local update.
    !-------------------------------------------------------------------------
    USE io_global,        ONLY : ionode, stdout
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
    USE mp_pools,         ONLY : inter_pool_comm
    USE mp_world,         ONLY : world_comm
    USE io_global,        ONLY : ionode_id
    USE becmod,           ONLY : becp, allocate_bec_type, deallocate_bec_type
    USE uspp,             ONLY : nkb
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xkcr(3,nkstot), omega0_ry, win_min_ry, win_max_ry
    INTEGER,  INTENT(IN) :: nbndskip_in, nrung
    CHARACTER(LEN=*), INTENT(IN) :: edmat_file, outfile

    INTEGER, PARAMETER :: NB = 400         ! v1: single batch (cross-batch Sigma
                                           ! blocks require all sources resident)
    INTEGER :: kg, ib, n, NU, NA_, NR1, iu, info, it, lwork, lrwork, liwork
    INTEGER :: ik, kgl, iks, a, b, i, a0, na_b, ibatch, nbatch, npw_kp, iters
    INTEGER :: nkb_d, nkb_p, ikb
    REAL(dp) :: w0, eR, xnk, alpha_ry, resid, rnorm, dnrm
    INTEGER,  ALLOCATABLE :: uk(:), ub(:), selA(:), sel1(:), pos1(:), iwork(:)
    INTEGER,  ALLOCATABLE :: igk_all(:,:), ngk_all(:)
    REAL(dp), ALLOCATABLE :: et_all(:,:), epsU(:), mu(:), rwork(:), gk_tmp(:)
    COMPLEX(dp), ALLOCATABLE :: MU_(:,:), H1(:,:), work(:), VAA(:,:), MA1(:,:), M1A(:,:)
    COMPLEX(dp), ALLOCATABLE :: cco(:,:), SgR1(:,:), Sgt(:,:), Vout(:,:)
    COMPLEX(dp), ALLOCATABLE :: evc_all(:,:,:), phi(:,:,:), stt(:,:,:), x(:,:,:), tv(:,:,:)
    COMPLEX(dp), ALLOCATABLE :: beta1(:,:), gco(:,:), Vf(:), psic(:), gbuf(:), rhs(:), chi(:)
    COMPLEX(dp), ALLOCATABLE :: vkb_d(:,:), vkb_p(:,:), becd(:,:), becp_(:,:)
    COMPLEX(dp), ALLOCATABLE :: cofd(:,:), cofp(:,:), evtmp(:,:), proj(:)
    COMPLEX(dp), ALLOCATABLE :: gam(:), tR(:,:), xi(:,:,:), tv2(:,:,:), sloc(:,:)
    REAL(dp) :: xk3(3)
    COMPLEX(dp) :: cone, czero, vij
    INTEGER, EXTERNAL :: global_kpoint_index

    cone = (1.0_dp,0.0_dp); czero = (0.0_dp,0.0_dp)
    w0 = omega0_ry; xnk = DBLE(nkstot)
    IF (npol /= 1) CALL errore('twolevel_physical','npol=1 only',1)
    iks = global_kpoint_index(nkstot, 1)

    ! ---- eigenvalues, manifold ----
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
    alpha_ry = 2.0_dp*ABS(w0 - MINVAL(et_all)) + 1.0_dp
    IF (ionode) THEN
       WRITE(stdout,'(/,5X,A)') REPEAT('=',64)
       WRITE(stdout,'(5X,A)') 'TWO-LEVEL rest dressing — MODE B (PHYSICAL tail, CG D2)'
       WRITE(stdout,'(5X,A,2I8,A,I3,A,F8.3,A)') 'N_A, N_R1 = ', NA_, NR1, &
            '   rungs = ', nrung, '   alpha = ', alpha_ry*rytoev, ' eV'
       FLUSH(stdout)
    ENDIF

    ! ---- explicit blocks + D1 (ionode) ----
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
       ! c = D1 M1A / N_k ;  SgR1_file = MA1 * c
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

    ! ---- replicate ALL wavefunctions + canonical igk (v1, 6x6-sized) ----
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
    ALLOCATE(vkb_d(npwx,nkb_d), vkb_p(npwx,nkb_p))
    ALLOCATE(Vf(dffts%nnr), psic(dffts%nnr), gbuf(dffts%nnr))
    ALLOCATE(rhs(npwx), chi(npwx), proj(nbnd), gam(nbnd))
    ALLOCATE(Sgt(NA_,NA_)); Sgt = czero

    nbatch = (NA_ + NB - 1) / NB
    DO ibatch = 1, nbatch
       a0 = (ibatch-1)*NB
       na_b = MIN(NB, NA_ - a0)

       ! ---- phi coefficients on the explicit basis: gamma(b,k;a) ----
       !  phi_a = psi_a + sum_{i in R1} c_i^a |i>
       ALLOCATE(phi(npwx,nkstot,na_b)); phi = czero
       DO a = 1, na_b
          DO kg = 1, nkstot
             gam = czero
             DO ib = 1, nbnd
                n = (kg-1)*nbnd + ib
                IF (pos1(n) > 0) gam(ib) = cco(pos1(n), a0+a)
             ENDDO
             IF (uk(selA(a0+a)) == kg) gam(ub(selA(a0+a))) = gam(ub(selA(a0+a))) + cone
             CALL ZGEMV('N', npwx, nbnd, cone, evc_all(1,1,kg), npwx, gam, 1, czero, phi(1,kg,a), 1)
          ENDDO
       ENDDO

       ! ---- dressed sources st = P_T V phi  (pool-local channels) ----
       ALLOCATE(stt(npwx,nkstot,na_b)); stt = czero
       CALL apply_dV_multik(phi, na_b, stt, .TRUE.)
       DEALLOCATE(phi)

       ! ---- ladder ----
       ALLOCATE(x(npwx,nkstot,na_b), tv(npwx,nkstot,na_b))
       ALLOCATE(beta1(NR1,na_b), gco(NR1,na_b))
       x = czero
       CALL allocate_bec_type(nkb, 1, becp)
       DO it = 0, nrung
          IF (it == 0) THEN
             tv = stt
          ELSE
             ! t = V x ; term2 via R1 coefficients
             CALL apply_dV_multik(x, na_b, tv, .FALSE.)
             beta1 = czero
             DO ik = 1, nks
                kg = ik + iks - 1
                DO a = 1, na_b
                   CALL ZGEMV('C', ngk_all(kg), nbnd, cone, evc_all(1,1,kg), npwx, &
                              tv(1,kg,a), 1, czero, proj, 1)
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
                CALL ZGEMM('N','N', NR1, na_b, NR1, cone, H1, NR1, tR, NR1, czero, gco, NR1)
                DEALLOCATE(tR)
             ENDIF
             CALL mp_bcast(gco, ionode_id, world_comm)
             ! xi-states from gco, apply V, add (1/Nk)*t : rhs per channel below
             ALLOCATE(xi(npwx,nkstot,na_b), tv2(npwx,nkstot,na_b))
             xi = czero
               DO a = 1, na_b
                  DO kg = 1, nkstot
                     gam = czero
                     DO ib = 1, nbnd
                        n = (kg-1)*nbnd + ib
                        IF (pos1(n) > 0) gam(ib) = gco(pos1(n), a)
                     ENDDO
                     CALL ZGEMV('N', npwx, nbnd, cone, evc_all(1,1,kg), npwx, gam, 1, czero, xi(1,kg,a), 1)
                  ENDDO
               ENDDO
             CALL apply_dV_multik(xi, na_b, tv2, .FALSE.)
             tv = stt + tv/xnk + tv2            ! s~ + W22 x   (tv2 already carries 1/Nk^2 via gco)
             DEALLOCATE(xi, tv2)
             ! P_T projection of the rhs
             DO ik = 1, nks
                kg = ik + iks - 1
                DO a = 1, na_b
                   CALL ZGEMV('C', ngk_all(kg), nbnd, cone, evc_all(1,1,kg), npwx, tv(1,kg,a), 1, czero, proj, 1)
                   CALL ZGEMV('N', ngk_all(kg), nbnd, -cone, evc_all(1,1,kg), npwx, proj, 1, cone, tv(1,kg,a), 1)
                ENDDO
             ENDDO
          ENDIF
          ! D2: pool-local CG solves;  x = -chi
          x = czero
          DO ik = 1, nks
             kg = ik + iks - 1
             CALL hpsi_setup_k(ik)
             DO a = 1, na_b
                rhs = czero; rhs(1:ngk_all(kg)) = tv(1:ngk_all(kg),kg,a)
                CALL solve_rest_cg(ik, rhs, w0, alpha_ry, nbnd, evc_all(:,:,kg), 1.d-8, 500, chi, iters, resid)
                x(1:ngk_all(kg),kg,a) = -chi(1:ngk_all(kg))
             ENDDO
          ENDDO
          DO kg = 1, nkstot
             CALL mp_sum(x(:,kg,:), inter_pool_comm)
          ENDDO
          IF (ionode) THEN
             WRITE(stdout,'(5X,A,I3,A,I3,A)') 'batch ', ibatch, '  rung ', it, ' done'
             FLUSH(stdout)
          ENDIF
       ENDDO
       CALL deallocate_bec_type(becp)

       ! ---- Sigma_tail contribution: (1/Nk) sum_k' <st_b|x_a> ----
       ALLOCATE(sloc(NA_,na_b)); sloc = czero
         DO ik = 1, nks
            kg = ik + iks - 1
            ! need st for ALL b at this channel: recompute? -> stored stt covers only batch
            ! v1: accumulate batch-diagonal style requires all-b sources; instead compute
            ! <st_b|x_a> with b restricted to the SAME batch and fill off-batch via
            ! hermiticity after all batches (valid since Sigma is Hermitian at real w0).
            CALL ZGEMM('C','N', na_b, na_b, ngk_all(kg), cone, stt(1,kg,1), npwx*nkstot, &
                       x(1,kg,1), npwx*nkstot, cone, sloc(a0+1,1), NA_)
         ENDDO
       CALL mp_sum(sloc, inter_pool_comm)
       IF (ionode) Sgt(:, a0+1:a0+na_b) = Sgt(:, a0+1:a0+na_b) + sloc/xnk
       DEALLOCATE(sloc)
       DEALLOCATE(stt, x, tv, beta1, gco)
    ENDDO

    ! ---- assemble + write (ionode) ----
    IF (ionode) THEN
       ! fill cross-batch blocks by hermiticity
       DO a = 1, NA_
          DO b = 1, NA_
             IF (ABS(Sgt(b,a)) == 0.0_dp .AND. ABS(Sgt(a,b)) > 0.0_dp) Sgt(b,a) = CONJG(Sgt(a,b))
          ENDDO
       ENDDO
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
       WRITE(stdout,'(5X,3A,I6,A)') 'wrote ', TRIM(outfile), '  (MODE B physical, N_A = ', NA_, ')'
       WRITE(stdout,'(5X,A)') REPEAT('=',64)
    ENDIF

  CONTAINS

    FUNCTION xk_dummy(kg_, xkcr_) RESULT(xkc3)
      !! cartesian 2pi/alat k-vector reconstructed from crystal coords
      USE cell_base, ONLY : bg
      INTEGER,  INTENT(IN) :: kg_
      REAL(dp), INTENT(IN) :: xkcr_(3,nkstot)
      REAL(dp) :: xkc3(3)
      xkc3 = MATMUL(bg, xkcr_(:,kg_))
    END FUNCTION xk_dummy

    SUBROUTINE apply_dV_multik(vin, nsrc, vout, project_T)
      !! vout^{k'} = [Delta-V vin]^{k'} for pool-local channels k' (zero elsewhere),
      !! local part via fold-multiply FFTs, KB nonlocal via summed home-k becs.
      !! project_T: subtract the nbnd explicit-band components at each channel.
      COMPLEX(dp), INTENT(IN)  :: vin(:,:,:)
      INTEGER,     INTENT(IN)  :: nsrc
      COMPLEX(dp), INTENT(OUT) :: vout(:,:,:)
      LOGICAL,     INTENT(IN)  :: project_T
      INTEGER :: ikl, kgc, kgs, aa, igl, ikb2
      COMPLEX(dp), ALLOCATABLE :: cd(:,:), cp(:,:), bd(:), bp(:), tmpd(:), tmpp(:)

      ! summed KB coefficients per source (home-k loop over local k, then global sum)
      ALLOCATE(cd(nkb_d,nsrc), cp(nkb_p,nsrc), bd(nkb_d), bp(nkb_p))
      ALLOCATE(tmpd(nkb_d), tmpp(nkb_p))
      cd = czero; cp = czero
      DO ikl = 1, nks
         kgs = ikl + iks - 1
         CALL get_betavkb(ngk_all(kgs), igk_all(1,kgs), xk(1,ikl), vkb_d, V_d%nat,V_d%ityp,V_d%tau, nkb_d)
         CALL get_betavkb(ngk_all(kgs), igk_all(1,kgs), xk(1,ikl), vkb_p, V_p%nat,V_p%ityp,V_p%tau, nkb_p)
         DO aa = 1, nsrc
            DO ikb2 = 1, nkb_d
               bd(ikb2) = SUM(CONJG(vkb_d(1:ngk_all(kgs),ikb2))*vin(1:ngk_all(kgs),kgs,aa))
            ENDDO
            DO ikb2 = 1, nkb_p
               bp(ikb2) = SUM(CONJG(vkb_p(1:ngk_all(kgs),ikb2))*vin(1:ngk_all(kgs),kgs,aa))
            ENDDO
            CALL make_coeff(V_d%nat,V_d%ityp,V_d%ntyp,nkb_d,bd,tmpd)
            CALL make_coeff(V_p%nat,V_p%ityp,V_p%ntyp,nkb_p,bp,tmpp)
            cd(:,aa) = cd(:,aa) + tmpd
            cp(:,aa) = cp(:,aa) + tmpp
         ENDDO
      ENDDO
      CALL mp_sum(cd, inter_pool_comm)
      CALL mp_sum(cp, inter_pool_comm)

      vout = czero
      DO ikl = 1, nks                       ! local channels
         kgc = ikl + iks - 1
         npw_kp = ngk_all(kgc)
         CALL get_betavkb(npw_kp, igk_all(1,kgc), xk(1,ikl), vkb_d, V_d%nat,V_d%ityp,V_d%tau, nkb_d)
         CALL get_betavkb(npw_kp, igk_all(1,kgc), xk(1,ikl), vkb_p, V_p%nat,V_p%ityp,V_p%tau, nkb_p)
         DO aa = 1, nsrc
            gbuf = czero
            DO kgs = 1, nkstot              ! source-k sum (fold)
               IF (MAXVAL(ABS(vin(1:ngk_all(kgs),kgs,aa))) < 1.d-14) CYCLE
               psic = czero
               DO igl = 1, ngk_all(kgs)
                  psic(dffts%nl(igk_all(igl,kgs))) = vin(igl,kgs,aa)
               ENDDO
               CALL invfft('Wave', psic, dffts)
               CALL build_V_folded(xkcr(:,kgs) - xkcr(:,kgc), Vf)
               gbuf = gbuf + Vf * psic
            ENDDO
            CALL fwfft('Wave', gbuf, dffts)
            DO igl = 1, npw_kp
               vout(igl,kgc,aa) = gbuf(dffts%nl(igk_all(igl,kgc)))
            ENDDO
            DO ikb2 = 1, nkb_d
               vout(1:npw_kp,kgc,aa) = vout(1:npw_kp,kgc,aa) + vkb_d(1:npw_kp,ikb2)*cd(ikb2,aa)
            ENDDO
            DO ikb2 = 1, nkb_p
               vout(1:npw_kp,kgc,aa) = vout(1:npw_kp,kgc,aa) - vkb_p(1:npw_kp,ikb2)*cp(ikb2,aa)
            ENDDO
            IF (project_T) THEN
               CALL ZGEMV('C', npw_kp, nbnd, cone, evc_all(1,1,kgc), npwx, vout(1,kgc,aa), 1, czero, proj, 1)
               CALL ZGEMV('N', npw_kp, nbnd, -cone, evc_all(1,1,kgc), npwx, proj, 1, cone, vout(1,kgc,aa), 1)
            ENDIF
         ENDDO
      ENDDO
      DO kgs = 1, nkstot
         CALL mp_sum(vout(:,kgs,1:nsrc), inter_pool_comm)
      ENDDO
      DEALLOCATE(cd, cp, bd, bp, tmpd, tmpp)
    END SUBROUTINE apply_dV_multik

  END SUBROUTINE twolevel_physical

END MODULE edt_twolevel
