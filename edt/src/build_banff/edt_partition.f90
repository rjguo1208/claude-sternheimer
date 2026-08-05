MODULE edt_partition
  !---------------------------------------------------------------------------
  !  Active/rest partition of the band structure (plan.md §3).
  !
  !  Active A(k) = kept bands whose energy lies in [win_min, win_max] (near omega0).
  !  Rest   R(k) = kept bands outside that window (treated by the Sternheimer solve).
  !  Q(k) = 1 - P(k) is applied as "project out the active manifold at k".
  !---------------------------------------------------------------------------
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: build_active_set, apply_Qproj, print_partition_audit

CONTAINS

  SUBROUTINE build_active_set(nbnd, nks, et_ev, kept, win_min, win_max, &
                               is_active, n_act, n_rest_k, &
                               n_active_tot, n_rest_tot, vbm, cbm, omega0_io, gap)
    !-------------------------------------------------------------------------
    !  Classify each (band, k) as active / rest / excluded.
    !  et_ev    : eigenvalues in eV (nbnd, nks)
    !  kept(ib) : .true. if band ib participates (not a semicore/excluded band)
    !  omega0_io: if > 9000 (sentinel) on entry, set to the VBM on exit
    !-------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,  INTENT(IN)  :: nbnd, nks
    REAL(dp), INTENT(IN)  :: et_ev(nbnd, nks)
    LOGICAL,  INTENT(IN)  :: kept(nbnd)
    REAL(dp), INTENT(IN)  :: win_min, win_max
    LOGICAL,  INTENT(OUT) :: is_active(nbnd, nks)
    INTEGER,  INTENT(OUT) :: n_act(nks), n_rest_k(nks)
    INTEGER,  INTENT(OUT) :: n_active_tot, n_rest_tot
    REAL(dp), INTENT(OUT) :: vbm, cbm, gap
    REAL(dp), INTENT(INOUT) :: omega0_io

    INTEGER  :: ik, ib
    REAL(dp) :: e
    REAL(dp), PARAMETER :: BIG = 1.0d30

    is_active = .FALSE.
    n_act = 0; n_rest_k = 0
    vbm = -BIG; cbm = BIG

    DO ik = 1, nks
       DO ib = 1, nbnd
          IF (.NOT. kept(ib)) CYCLE
          e = et_ev(ib, ik)
          IF (e >= win_min .AND. e <= win_max) THEN
             is_active(ib, ik) = .TRUE.
             n_act(ik) = n_act(ik) + 1
             IF (e > vbm) vbm = e
          ELSE
             n_rest_k(ik) = n_rest_k(ik) + 1
             IF (e > win_max .AND. e < cbm) cbm = e
          ENDIF
       ENDDO
    ENDDO
    n_active_tot = SUM(n_act)
    n_rest_tot   = SUM(n_rest_k)

    ! default omega0 = top of the active manifold (VBM)
    IF (omega0_io > 9000.0_dp) omega0_io = vbm

    ! gap = distance from omega0 to the nearest REST state
    gap = BIG
    DO ik = 1, nks
       DO ib = 1, nbnd
          IF (.NOT. kept(ib)) CYCLE
          IF (is_active(ib, ik)) CYCLE
          gap = MIN(gap, ABS(et_ev(ib, ik) - omega0_io))
       ENDDO
    ENDDO
  END SUBROUTINE build_active_set


  SUBROUTINE apply_Qproj(npw, npol, nact, evc_act, psi)
    !-------------------------------------------------------------------------
    !  In place: psi <- Q psi = psi - sum_{a in active} |a><a|psi>
    !  evc_act(ld, nact) : active Bloch states at this k (periodic part, G-space)
    !  Norm-conserving: overlap S = 1.  (plan.md §3)
    !-------------------------------------------------------------------------
    USE mp,        ONLY : mp_sum
    USE mp_bands,  ONLY : intra_bgrp_comm
    IMPLICIT NONE
    INTEGER,     INTENT(IN)    :: npw, npol, nact
    COMPLEX(dp), INTENT(IN)    :: evc_act(:, :)
    COMPLEX(dp), INTENT(INOUT) :: psi(:)
    COMPLEX(dp), ALLOCATABLE :: ovlp(:)
    INTEGER :: a, n

    IF (nact <= 0) RETURN
    n = npw * npol
    ALLOCATE(ovlp(nact))
    ! ovlp(a) = <evc_act(:,a) | psi>
    CALL ZGEMV('C', n, nact, (1.0_dp,0.0_dp), evc_act, SIZE(evc_act,1), &
               psi, 1, (0.0_dp,0.0_dp), ovlp, 1)
    CALL mp_sum(ovlp, intra_bgrp_comm)
    DO a = 1, nact
       psi(1:n) = psi(1:n) - ovlp(a) * evc_act(1:n, a)
    ENDDO
    DEALLOCATE(ovlp)
  END SUBROUTINE apply_Qproj


  SUBROUTINE print_partition_audit(nbnd, nks, nbndep, nbndsub, nexcl, &
                                    win_min, win_max, omega0, eta, &
                                    n_act, n_rest_k, n_active_tot, n_rest_tot, &
                                    vbm, cbm, gap)
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: nbnd, nks, nbndep, nbndsub, nexcl
    REAL(dp), INTENT(IN) :: win_min, win_max, omega0, eta, vbm, cbm, gap
    INTEGER,  INTENT(IN) :: n_act(nks), n_rest_k(nks), n_active_tot, n_rest_tot
    INTEGER :: nmin_a, nmax_a, nmin_r, nmax_r
    REAL(dp):: avg_a

    nmin_a = MINVAL(n_act);    nmax_a = MAXVAL(n_act)
    nmin_r = MINVAL(n_rest_k); nmax_r = MAXVAL(n_rest_k)
    avg_a  = DBLE(n_active_tot) / DBLE(MAX(1,nks))

    WRITE(stdout,'(/,5X,A)') REPEAT('=',64)
    WRITE(stdout,'(5X,A)')   'EDT active / rest partition audit'
    WRITE(stdout,'(5X,A)')   REPEAT('=',64)
    WRITE(stdout,'(5X,A,I6)')      'total bands  nbnd        = ', nbnd
    WRITE(stdout,'(5X,A,I6)')      'semicore excluded        = ', nexcl
    WRITE(stdout,'(5X,A,I6)')      'bands in play (kept)     = ', nbndep
    WRITE(stdout,'(5X,A,I6)')      'k-points     nks         = ', nks
    WRITE(stdout,'(5X,A,I6)')      'Wannier manifold N_W     = ', nbndsub
    WRITE(stdout,'(5X,A)')   REPEAT('-',64)
    WRITE(stdout,'(5X,A,2F10.4,A)')'active window            = [', win_min, win_max, ' ] eV'
    WRITE(stdout,'(5X,A,F10.4,A)') 'omega0 (reference)       = ', omega0, ' eV'
    WRITE(stdout,'(5X,A,F10.5,A)') 'eta (complex shift)      = ', eta, ' eV'
    WRITE(stdout,'(5X,A,F10.4,A)') 'VBM (top of active)      = ', vbm, ' eV'
    WRITE(stdout,'(5X,A,F10.4,A)') 'CBM (bottom of rest>win) = ', cbm, ' eV'
    WRITE(stdout,'(5X,A,F10.4,A)') 'gap  omega0 -> nearest R = ', gap, ' eV'
    WRITE(stdout,'(5X,A)')   REPEAT('-',64)
    WRITE(stdout,'(5X,A,I8)')      'active states (total)    = ', n_active_tot
    WRITE(stdout,'(5X,A,I4,A,I4,A,F6.2)') &
                                   'active per k  min/max/avg= ', nmin_a, ' /', nmax_a, ' /', avg_a
    WRITE(stdout,'(5X,A,I8)')      'rest states  (total)     = ', n_rest_tot
    WRITE(stdout,'(5X,A,I4,A,I4)') 'rest per k    min/max    = ', nmin_r, ' /', nmax_r
    IF (nmin_a /= nmax_a) THEN
       WRITE(stdout,'(5X,A)') 'NOTE: active count varies across k (window crosses a band).'
    ENDIF
    IF (nbndsub > 0 .AND. nmax_a /= nbndsub) THEN
       WRITE(stdout,'(5X,A,I4,A,I4,A)') &
            'NOTE: max active per k (', nmax_a, ') /= N_W from filukk (', nbndsub, ').'
    ENDIF
    WRITE(stdout,'(5X,A)')   REPEAT('=',64)
    FLUSH(stdout)
  END SUBROUTINE print_partition_audit

END MODULE edt_partition
