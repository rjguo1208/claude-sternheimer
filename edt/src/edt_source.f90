MODULE edt_source
  !---------------------------------------------------------------------------
  !  Stage B (plan.md §5): the Sternheimer source  |s> = Q ΔV |n k_src>.
  !
  !  ΔV |ψ_{n,k_src}> is built as a KET resolved on an output channel k_out:
  !    local:    ζ_loc(r) = V_folded^{q=k_src-k_out}(r) · u_{n,k_src}(r)   (fwfft -> k_out)
  !    nonlocal: ζ_nl(G)  = Σ_a Σ_{ih,jh} β^{d}_{a,ih}(k_out,G) D_{ih,jh} <β^{d}_{a,jh}(k_src)|ψ_{n,k_src}>
  !                       - (defect -> pristine)
  !  Overlapping with <ψ_{m,k_out}| reproduces EDI's matrix element
  !    <m k_out|ΔV|n k_src> = M_loc + (M_nl^d - M_nl^p),
  !  i.e. exactly the local kernel of ed_coarse_full_q and the KB contraction with
  !  dvan.  This is the Born-limit check (T1).  (Q-projection + grid assembly next.)
  !---------------------------------------------------------------------------
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: build_V_folded, test_source

  REAL(dp), PARAMETER :: tpi = 6.283185307179586_dp

CONTAINS

  SUBROUTINE build_V_folded(q_cryst, Vf)
    !-------------------------------------------------------------------------
    ! Vf(dffts%nnr): V_colin folded onto the primitive grid with phase
    ! exp(i q·r), q crystal coords.  Mirrors ed_coarse_full_q V_folded_kf
    ! (exact q, no BZ wrap — the double_ft_bug fix).
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


  SUBROUTINE count_nkb(vf_nat, vf_ityp, vf_ntyp, nkb)
    USE uspp_param, ONLY : nh
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: vf_nat, vf_ntyp, vf_ityp(vf_nat)
    INTEGER, INTENT(OUT) :: nkb
    INTEGER :: na, nt
    nkb = 0
    DO nt = 1, vf_ntyp
       DO na = 1, vf_nat
          IF (vf_ityp(na) == nt) nkb = nkb + nh(nt)
       ENDDO
    ENDDO
  END SUBROUTINE count_nkb


  SUBROUTINE test_source(ik_a, ik_b, q_cryst_ba, blist, nb)
    !-------------------------------------------------------------------------
    ! T1.  For bands `blist` at (k_a=out, k_b=src) compare the full matrix
    ! element from EDI's kernels (M_ref) with the source-ket route (M_ket):
    !   local:    M_loc_ref = (1/nnr) Σ_r u*_{i,a} V_folded^{b-a} u_{j,b}
    !   nonlocal: M_nl_ref  = Σ_ikb conj(bec_a(ikb,i)) Σ_jkb dvan bec_b(jkb,j)  (d - p)
    ! M_ket builds ζ = ζ_loc + ζ_nl on the a-basis and overlaps with ψ_{i,a}.
    !-------------------------------------------------------------------------
    USE io_global,        ONLY : ionode, stdout
    USE wvfct,            ONLY : npwx, nbnd
    USE klist,            ONLY : ngk, igk_k, xk
    USE fft_base,         ONLY : dffts
    USE fft_interfaces,   ONLY : invfft, fwfft
    USE noncollin_module, ONLY : npol
    USE pw_restart_new,   ONLY : read_collected_wfc
    USE io_files,         ONLY : restart_dir
    USE edic_mod,         ONLY : V_d, V_p
    USE uspp,             ONLY : dvan
    USE uspp_param,       ONLY : nh
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: ik_a, ik_b, nb, blist(nb)
    REAL(dp), INTENT(IN) :: q_cryst_ba(3)

    COMPLEX(dp), ALLOCATABLE :: evc_a(:,:), evc_b(:,:)
    COMPLEX(dp), ALLOCATABLE :: ua(:,:), ub(:,:), Vf(:), g(:), psic(:), zeta(:)
    COMPLEX(dp), ALLOCATABLE :: Mloc_ref(:,:), Mloc_ket(:,:), Mnl_ref(:,:), Mnl_ket(:,:)
    COMPLEX(dp), ALLOCATABLE :: vkb_d_a(:,:), vkb_d_b(:,:), vkb_p_a(:,:), vkb_p_b(:,:)
    COMPLEX(dp), ALLOCATABLE :: becd_a(:,:), becd_b(:,:), becp_a(:,:), becp_b(:,:)
    COMPLEX(dp), ALLOCATABLE :: coeff(:), znl(:)
    INTEGER  :: i, j, ig, ir, nnr, npw_a, npw_b, nkb_d, nkb_p
    REAL(dp) :: dloc, dnl, dfull
    COMPLEX(dp) :: acc

    IF (.NOT. ionode) RETURN
    nnr = dffts%nnr
    npw_a = ngk(ik_a); npw_b = ngk(ik_b)

    ALLOCATE(evc_a(npwx*npol, nbnd), evc_b(npwx*npol, nbnd))
    CALL read_collected_wfc(restart_dir(), ik_a, evc_a)
    CALL read_collected_wfc(restart_dir(), ik_b, evc_b)

    ! ---------- local part ----------
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
    ALLOCATE(Mloc_ref(nb,nb), Mloc_ket(nb,nb), g(nnr), zeta(npwx))
    DO j = 1, nb
       DO i = 1, nb
          acc = (0.0_dp,0.0_dp)
          DO ir = 1, nnr
             acc = acc + CONJG(ua(ir,i)) * Vf(ir) * ub(ir,j)
          ENDDO
          Mloc_ref(i,j) = acc / DBLE(nnr)
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
          Mloc_ket(i,j) = acc
       ENDDO
    ENDDO

    ! ---------- nonlocal part (KB, defect - pristine) ----------
    CALL count_nkb(V_d%nat, V_d%ityp, V_d%ntyp, nkb_d)
    CALL count_nkb(V_p%nat, V_p%ityp, V_p%ntyp, nkb_p)
    ALLOCATE(vkb_d_a(npwx,nkb_d), vkb_d_b(npwx,nkb_d), vkb_p_a(npwx,nkb_p), vkb_p_b(npwx,nkb_p))
    CALL get_betavkb(ngk(ik_a), igk_k(1,ik_a), xk(1,ik_a), vkb_d_a, V_d%nat, V_d%ityp, V_d%tau, nkb_d)
    CALL get_betavkb(ngk(ik_b), igk_k(1,ik_b), xk(1,ik_b), vkb_d_b, V_d%nat, V_d%ityp, V_d%tau, nkb_d)
    CALL get_betavkb(ngk(ik_a), igk_k(1,ik_a), xk(1,ik_a), vkb_p_a, V_p%nat, V_p%ityp, V_p%tau, nkb_p)
    CALL get_betavkb(ngk(ik_b), igk_k(1,ik_b), xk(1,ik_b), vkb_p_b, V_p%nat, V_p%ityp, V_p%tau, nkb_p)

    ALLOCATE(becd_a(nkb_d,nb), becd_b(nkb_d,nb), becp_a(nkb_p,nb), becp_b(nkb_p,nb))
    CALL proj_bec(npw_a, npwx, nkb_d, nb, nbnd, vkb_d_a, evc_a, blist, becd_a)
    CALL proj_bec(npw_b, npwx, nkb_d, nb, nbnd, vkb_d_b, evc_b, blist, becd_b)
    CALL proj_bec(npw_a, npwx, nkb_p, nb, nbnd, vkb_p_a, evc_a, blist, becp_a)
    CALL proj_bec(npw_b, npwx, nkb_p, nb, nbnd, vkb_p_b, evc_b, blist, becp_b)

    ALLOCATE(Mnl_ref(nb,nb), Mnl_ket(nb,nb), coeff(MAX(nkb_d,nkb_p)), znl(npwx))
    Mnl_ref = (0.0_dp,0.0_dp); Mnl_ket = (0.0_dp,0.0_dp)
    ! defect (+) then pristine (-)
    CALL nl_contrib(+1.0_dp, V_d%nat, V_d%ityp, V_d%ntyp, npw_a, npwx, nkb_d, nb, nbnd, &
                    becd_a, becd_b, vkb_d_a, evc_a, blist, Mnl_ref, Mnl_ket)
    CALL nl_contrib(-1.0_dp, V_p%nat, V_p%ityp, V_p%ntyp, npw_a, npwx, nkb_p, nb, nbnd, &
                    becp_a, becp_b, vkb_p_a, evc_a, blist, Mnl_ref, Mnl_ket)

    ! ---------- compare ----------
    dloc = 0.0_dp; dnl = 0.0_dp; dfull = 0.0_dp
    DO j = 1, nb
       DO i = 1, nb
          dloc  = MAX(dloc,  ABS(Mloc_ref(i,j) - Mloc_ket(i,j)))
          dnl   = MAX(dnl,   ABS(Mnl_ref(i,j)  - Mnl_ket(i,j)))
          dfull = MAX(dfull, ABS((Mloc_ref(i,j)+Mnl_ref(i,j)) - (Mloc_ket(i,j)+Mnl_ket(i,j))))
       ENDDO
    ENDDO

    WRITE(stdout,'(/,5X,A)') REPEAT('=',64)
    WRITE(stdout,'(5X,A)') 'T1 (Born limit) — full source-ket vs EDI kernels (local + nonlocal KB)'
    WRITE(stdout,'(5X,A,I4,A,I4,A,3F8.4)') 'k_a=', ik_a, '  k_b=', ik_b, '  q(cryst)=', q_cryst_ba
    WRITE(stdout,'(5X,A,I4,A,I5,A,I5)') 'bands=', nb, '   nkb_d=', nkb_d, '   nkb_p=', nkb_p
    WRITE(stdout,'(5X,A,ES12.4,A)') '|M_loc| max         = ', MAXVAL(ABS(Mloc_ref)), ' Ry'
    WRITE(stdout,'(5X,A,ES12.4,A)') '|M_nl|  max         = ', MAXVAL(ABS(Mnl_ref)),  ' Ry'
    WRITE(stdout,'(5X,A,ES12.4)')   'max|dM| local       = ', dloc
    WRITE(stdout,'(5X,A,ES12.4)')   'max|dM| nonlocal    = ', dnl
    WRITE(stdout,'(5X,A,ES12.4)')   'max|dM| full        = ', dfull
    IF (dfull < 1.0d-10) THEN
       WRITE(stdout,'(5X,A)') 'PASS: full source ket reproduces EDI M (local+nonlocal) to < 1e-10.'
    ELSE
       WRITE(stdout,'(5X,A)') 'CHECK: difference exceeds 1e-10.'
    ENDIF
    WRITE(stdout,'(5X,A)') REPEAT('=',64)
    FLUSH(stdout)

    DEALLOCATE(evc_a, evc_b, ua, ub, psic, Vf, g, zeta, Mloc_ref, Mloc_ket)
    DEALLOCATE(vkb_d_a, vkb_d_b, vkb_p_a, vkb_p_b, becd_a, becd_b, becp_a, becp_b)
    DEALLOCATE(Mnl_ref, Mnl_ket, coeff, znl)
  END SUBROUTINE test_source


  SUBROUTINE proj_bec(npw, ld, nkb, nb, nbnd, vkb, evc, blist, bec)
    ! bec(ikb,i) = <beta_ikb | psi_{blist(i)}> = Σ_G conj(vkb(G,ikb)) evc(G,blist(i))
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: npw, ld, nkb, nb, nbnd, blist(nb)
    COMPLEX(dp), INTENT(IN)  :: vkb(ld,nkb), evc(ld,nbnd)
    COMPLEX(dp), INTENT(OUT) :: bec(nkb,nb)
    INTEGER :: ikb, i, ig
    COMPLEX(dp) :: acc
    DO i = 1, nb
       DO ikb = 1, nkb
          acc = (0.0_dp,0.0_dp)
          DO ig = 1, npw
             acc = acc + CONJG(vkb(ig,ikb)) * evc(ig,blist(i))
          ENDDO
          bec(ikb,i) = acc
       ENDDO
    ENDDO
  END SUBROUTINE proj_bec


  SUBROUTINE nl_contrib(sgn, vf_nat, vf_ityp, vf_ntyp, npw_a, ld, nkb, nb, nbnd, &
                         bec_a, bec_b, vkb_a, evc_a, blist, Mref, Mket)
    ! Accumulate sgn * [ KB nonlocal block ] into Mref (direct) and Mket (ket route).
    USE uspp,       ONLY : dvan
    USE uspp_param, ONLY : nh
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: sgn
    INTEGER,  INTENT(IN) :: vf_nat, vf_ntyp, vf_ityp(vf_nat), npw_a, ld, nkb, nb, nbnd, blist(nb)
    COMPLEX(dp), INTENT(IN) :: bec_a(nkb,nb), bec_b(nkb,nb), vkb_a(ld,nkb), evc_a(ld,nbnd)
    COMPLEX(dp), INTENT(INOUT) :: Mref(nb,nb), Mket(nb,nb)
    COMPLEX(dp), ALLOCATABLE :: coeff(:), znl(:)
    INTEGER :: i, j, ig, na, nt, ih, jh, ijkb0
    COMPLEX(dp) :: acc

    ALLOCATE(coeff(nkb), znl(npw_a))
    DO j = 1, nb
       ! coeff(ikb) = Σ_jkb dvan(ih,jh,nt) bec_b(jkb,j)   for source band j
       coeff = (0.0_dp,0.0_dp)
       ijkb0 = 0
       DO nt = 1, vf_ntyp
          DO na = 1, vf_nat
             IF (vf_ityp(na) == nt) THEN
                DO ih = 1, nh(nt)
                   DO jh = 1, nh(nt)
                      coeff(ijkb0+ih) = coeff(ijkb0+ih) + dvan(ih,jh,nt) * bec_b(ijkb0+jh, j)
                   ENDDO
                ENDDO
                ijkb0 = ijkb0 + nh(nt)
             ENDIF
          ENDDO
       ENDDO
       ! direct: Mref(i,j) += sgn Σ_ikb conj(bec_a(ikb,i)) coeff(ikb)
       DO i = 1, nb
          acc = (0.0_dp,0.0_dp)
          DO ig = 1, nkb
             acc = acc + CONJG(bec_a(ig,i)) * coeff(ig)
          ENDDO
          Mref(i,j) = Mref(i,j) + sgn * acc
       ENDDO
       ! ket: znl(G) = Σ_ikb vkb_a(G,ikb) coeff(ikb);  Mket(i,j) += sgn <ψ_{i,a}|znl>
       znl = (0.0_dp,0.0_dp)
       DO ig = 1, nkb
          znl(1:npw_a) = znl(1:npw_a) + vkb_a(1:npw_a,ig) * coeff(ig)
       ENDDO
       DO i = 1, nb
          acc = (0.0_dp,0.0_dp)
          DO ig = 1, npw_a
             acc = acc + CONJG(evc_a(ig,blist(i))) * znl(ig)
          ENDDO
          Mket(i,j) = Mket(i,j) + sgn * acc
       ENDDO
    ENDDO
    DEALLOCATE(coeff, znl)
  END SUBROUTINE nl_contrib

END MODULE edt_source
