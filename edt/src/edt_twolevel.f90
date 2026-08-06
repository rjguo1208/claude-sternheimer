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
    IF (nsplit <= 0 .OR. nsplit >= nbnd) CALL errore('vtilde_block_twolevel', &
         'this version implements MODE A only: need 0 < tail_split_band < nbnd', 1)

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

END MODULE edt_twolevel
