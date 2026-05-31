MODULE edt_source
  !---------------------------------------------------------------------------
  !  Stage B (plan.md §5): the Sternheimer source  |s> = Q ΔV |n k_i>.
  !
  !  ΔV |ψ_{n,k_src}> is built as a KET, resolved on an output channel k_out,
  !  by reusing EDI's exact-q real-space fold of V_colin:
  !       ζ(r) = V_folded^{q=k_src-k_out}(r) · u_{n,k_src}(r),     fwfft -> k_out basis.
  !  Overlapping with <ψ_{m,k_out}| reproduces EDI's local matrix element
  !       <m k_out | ΔV_loc | n k_src> = (1/nnr) Σ_r u*_{m,k_out} V_folded u_{n,k_src},
  !  which is the Born-limit check (T1, local part).  (Nonlocal KB and the
  !  Q-projection follow.)
  !---------------------------------------------------------------------------
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: build_V_folded, test_source_local

  REAL(dp), PARAMETER :: tpi = 6.283185307179586_dp

CONTAINS

  SUBROUTINE build_V_folded(q_cryst, Vf)
    !-------------------------------------------------------------------------
    ! Vf(dffts%nnr): V_colin folded onto the primitive grid with the exact
    ! phase exp(i q·r), q in crystal coords.  Mirrors ed_coarse_full_q's
    ! V_folded_kf (exact q, no BZ wrap — the double_ft_bug fix).
    !-------------------------------------------------------------------------
    USE edic_mod, ONLY : V_colin, V_d
    USE fft_base, ONLY : dffts
    IMPLICIT NONE
    REAL(dp),    INTENT(IN)  :: q_cryst(3)
    COMPLEX(dp), INTENT(OUT) :: Vf(:)
    INTEGER  :: irx, iry, irz, inr, irnmod, ir1mod, ir2mod, ir3mod
    REAL(dp) :: d1, d2, d3, arg
    COMPLEX(dp) :: phase
    LOGICAL  :: isq0

    Vf = (0.0_dp, 0.0_dp)
    isq0 = (ABS(q_cryst(1)) < 1.0d-8 .AND. ABS(q_cryst(2)) < 1.0d-8 .AND. ABS(q_cryst(3)) < 1.0d-8)
    d1 = tpi * q_cryst(1) / dffts%nr1
    d2 = tpi * q_cryst(2) / dffts%nr2
    d3 = tpi * q_cryst(3) / dffts%nr3
    inr = 0
    DO irz = 0, V_d%nr3 - 1
       ir3mod = irz - (irz / dffts%nr3) * dffts%nr3
       DO iry = 0, V_d%nr2 - 1
          ir2mod = iry - (iry / dffts%nr2) * dffts%nr2
          DO irx = 0, V_d%nr1 - 1
             ir1mod = irx - (irx / dffts%nr1) * dffts%nr1
             inr = inr + 1
             irnmod = ir3mod * dffts%nr1 * dffts%nr2 + ir2mod * dffts%nr1 + ir1mod + 1
             IF (isq0) THEN
                Vf(irnmod) = Vf(irnmod) + CMPLX(V_colin(inr), 0.0_dp, dp)
             ELSE
                arg = irx * d1 + iry * d2 + irz * d3
                phase = CMPLX(COS(arg), SIN(arg), dp)
                Vf(irnmod) = Vf(irnmod) + V_colin(inr) * phase
             ENDIF
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE build_V_folded


  SUBROUTINE test_source_local(ik_a, ik_b, q_cryst_ba, blist, nb)
    !-------------------------------------------------------------------------
    ! T1 (local).  For bands `blist` at (k_a,k_b), compare
    !   M_ref(i,j) = (1/nnr) Σ_r u*_{i,k_a} V_folded^{k_b-k_a} u_{j,k_b}   [EDI kernel]
    !   M_ket(i,j) = <ψ_{i,k_a}| ζ_j>,  ζ_j = fwfft(V_folded·u_{j,k_b}) on k_a basis
    ! q_cryst_ba = xk_cryst(k_b) - xk_cryst(k_a).  Scalar (npol=1).
    !-------------------------------------------------------------------------
    USE io_global,        ONLY : ionode, stdout
    USE wvfct,            ONLY : npwx, nbnd
    USE klist,            ONLY : ngk, igk_k
    USE fft_base,         ONLY : dffts
    USE fft_interfaces,   ONLY : invfft, fwfft
    USE noncollin_module, ONLY : npol
    USE pw_restart_new,   ONLY : read_collected_wfc
    USE io_files,         ONLY : restart_dir
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: ik_a, ik_b, nb, blist(nb)
    REAL(dp), INTENT(IN) :: q_cryst_ba(3)

    COMPLEX(dp), ALLOCATABLE :: evc_a(:,:), evc_b(:,:)
    COMPLEX(dp), ALLOCATABLE :: ua(:,:), ub(:,:), Vf(:), g(:), psic(:), zeta(:)
    COMPLEX(dp), ALLOCATABLE :: Mref(:,:), Mket(:,:)
    INTEGER  :: i, j, ig, ir, nnr, npw_a, npw_b
    REAL(dp) :: maxdiff
    COMPLEX(dp) :: acc

    IF (.NOT. ionode) RETURN          ! single-rank validation
    nnr = dffts%nnr
    npw_a = ngk(ik_a); npw_b = ngk(ik_b)

    ALLOCATE(evc_a(npwx*npol, nbnd), evc_b(npwx*npol, nbnd))
    CALL read_collected_wfc(restart_dir(), ik_a, evc_a)
    CALL read_collected_wfc(restart_dir(), ik_b, evc_b)

    ALLOCATE(ua(nnr,nb), ub(nnr,nb), psic(nnr))
    DO i = 1, nb
       psic = (0.0_dp,0.0_dp)
       DO ig = 1, npw_a
          psic(dffts%nl(igk_k(ig,ik_a))) = evc_a(ig, blist(i))
       ENDDO
       CALL invfft('Wave', psic, dffts); ua(:,i) = psic
       psic = (0.0_dp,0.0_dp)
       DO ig = 1, npw_b
          psic(dffts%nl(igk_k(ig,ik_b))) = evc_b(ig, blist(i))
       ENDDO
       CALL invfft('Wave', psic, dffts); ub(:,i) = psic
    ENDDO

    ALLOCATE(Vf(nnr)); CALL build_V_folded(q_cryst_ba, Vf)

    ALLOCATE(Mref(nb,nb), Mket(nb,nb), g(nnr), zeta(npwx))
    DO j = 1, nb
       DO i = 1, nb
          acc = (0.0_dp,0.0_dp)
          DO ir = 1, nnr
             acc = acc + CONJG(ua(ir,i)) * Vf(ir) * ub(ir,j)
          ENDDO
          Mref(i,j) = acc / DBLE(nnr)
       ENDDO
    ENDDO
    DO j = 1, nb
       g = Vf(:) * ub(:,j)
       CALL fwfft('Wave', g, dffts)
       zeta = (0.0_dp,0.0_dp)
       DO ig = 1, npw_a
          zeta(ig) = g(dffts%nl(igk_k(ig,ik_a)))
       ENDDO
       DO i = 1, nb
          acc = (0.0_dp,0.0_dp)
          DO ig = 1, npw_a
             acc = acc + CONJG(evc_a(ig,blist(i))) * zeta(ig)
          ENDDO
          Mket(i,j) = acc
       ENDDO
    ENDDO

    maxdiff = 0.0_dp
    DO j = 1, nb
       DO i = 1, nb
          maxdiff = MAX(maxdiff, ABS(Mref(i,j) - Mket(i,j)))
       ENDDO
    ENDDO

    WRITE(stdout,'(/,5X,A)') REPEAT('=',64)
    WRITE(stdout,'(5X,A)') 'T1 (local) — source-ket vs EDI local kernel'
    WRITE(stdout,'(5X,A,I4,A,I4,A,3F8.4)') 'k_a=', ik_a, '  k_b=', ik_b, '  q(cryst)=', q_cryst_ba
    WRITE(stdout,'(5X,A,I3)')      'bands tested        = ', nb
    WRITE(stdout,'(5X,A,ES12.4,A)')'|M_ref| max         = ', MAXVAL(ABS(Mref)), ' Ry'
    WRITE(stdout,'(5X,A,ES12.4)')  'max|M_ref - M_ket|  = ', maxdiff
    IF (maxdiff < 1.0d-10) THEN
       WRITE(stdout,'(5X,A)') 'PASS: local source ket reproduces EDI local M to < 1e-10.'
    ELSE
       WRITE(stdout,'(5X,A)') 'CHECK: difference exceeds 1e-10 (normalization/convention).'
    ENDIF
    WRITE(stdout,'(5X,A)') REPEAT('=',64)
    FLUSH(stdout)

    DEALLOCATE(evc_a, evc_b, ua, ub, psic, Vf, g, zeta, Mref, Mket)
  END SUBROUTINE test_source_local

END MODULE edt_source
