module hydro_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_acc_module, only: acc_stream

  implicit none

contains

  subroutine apply_av(lo, hi, idir, dx, &
                      div, div_lo, div_hi, &
                      uin, uin_lo, uin_hi, &
                      flux, f_lo, f_hi) bind(c, name="apply_av")

    use amrex_constants_module, only: ZERO, FOURTH
    use castro_module, only: NVAR, UTEMP

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: idir

    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    integer :: i, j, k, n

    real(rt) :: div1

    real(rt), parameter :: difmag = 0.1d0

    do n = 1, NVAR

       if ( n == UTEMP ) cycle

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                if (idir .eq. 1) then

                   div1 = FOURTH * (div(i,j,k  ) + div(i,j+1,k) + &
                                    div(i,j,k+1) + div(i,j+1,k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i-1,j,k,n))

                else if (idir .eq. 2) then

                   div1 = FOURTH * (div(i,j,k  ) + div(i+1,j,k) + &
                                    div(i,j,k+1) + div(i+1,j,k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i,j-1,k,n))

                else

                   div1 = FOURTH * (div(i,j  ,k) + div(i+1,j  ,k) + &
                                    div(i,j+1,k) + div(i+1,j+1,k))
                   div1 = difmag * min(ZERO, div1)
                   div1 = div1 * (uin(i,j,k,n) - uin(i,j,k-1,n))

                end if

                flux(i,j,k,n) = flux(i,j,k,n) + dx(idir) * div1

             end do
          end do
       end do

    end do

  end subroutine apply_av

  
  subroutine normalize_species_fluxes(lo, hi, flux, f_lo, f_hi) bind(c, name="normalize_species_fluxes")
    ! Normalize the fluxes of the mass fractions so that
    ! they sum to 0.  This is essentially the CMA procedure that is
    ! defined in Plewa & Muller, 1999, A&A, 342, 179.

    use network, only: nspec
    use amrex_constants_module, only: ZERO, ONE
    use amrex_fort_module, only: rt => amrex_real
    use castro_module, only: NVAR, URHO, UFS

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    ! Local variables
    integer  :: i, j, k, n
    real(rt) :: sum, fac

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             sum = ZERO

             do n = UFS, UFS+nspec-1
                sum = sum + flux(i,j,k,n)
             end do

             if (sum .ne. ZERO) then
                fac = flux(i,j,k,URHO) / sum
             else
                fac = ONE
             end if

             do n = UFS, UFS+nspec-1
                flux(i,j,k,n) = flux(i,j,k,n) * fac
             end do

          end do
       end do
    end do

  end subroutine normalize_species_fluxes



  subroutine calc_pdivu(lo, hi, &
                        q1, q1_lo, q1_hi, &
                        area1, a1_lo, a1_hi, &
                        q2, q2_lo, q2_hi, &
                        area2, a2_lo, a2_hi, &
                        q3, q3_lo, q3_hi, &
                        area3, a3_lo, a3_hi, &
                        vol, v_lo, v_hi, &
                        dx, pdivu, div_lo, div_hi)
    ! this computes the cell-centered p div(U) term from the
    ! edge-centered Godunov state.  This is used in the internal energy
    ! update

    use castro_module, only: NGDNV, GDPRES, GDU, GDV, GDW
    use amrex_constants_module, only: HALF
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: div_lo(3), div_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: pdivu(div_lo(1):div_hi(1),div_lo(2):div_hi(2),div_lo(3):div_hi(3))
    integer, intent(in) :: q1_lo(3), q1_hi(3)
    integer, intent(in) :: a1_lo(3), a1_hi(3)
    real(rt), intent(in) :: q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
    real(rt), intent(in) :: area1(a1_lo(1):a1_hi(1),a1_lo(2):a1_hi(2),a1_lo(3):a1_hi(3))
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: a2_lo(3), a2_hi(3)
    real(rt), intent(in) :: q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
    real(rt), intent(in) :: area2(a2_lo(1):a2_hi(1),a1_lo(2):a1_hi(2),a1_lo(3):a1_hi(3))
    integer, intent(in) :: q3_lo(3), q3_hi(3)
    integer, intent(in) :: a3_lo(3), a3_hi(3)
    real(rt), intent(in) :: q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)
    real(rt), intent(in) :: area3(a3_lo(1):a3_hi(1),a1_lo(2):a1_hi(2),a1_lo(3):a1_hi(3))
    integer, intent(in) :: v_lo(3), v_hi(3)
    real(rt), intent(in) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))

    integer  :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             pdivu(i,j,k) = &
                  HALF*(q1(i+1,j,k,GDPRES) + q1(i,j,k,GDPRES)) * &
                  (q1(i+1,j,k,GDU) - q1(i,j,k,GDU))/dx(1) + &
                  HALF*(q2(i,j+1,k,GDPRES) + q2(i,j,k,GDPRES)) * &
                  (q2(i,j+1,k,GDV) - q2(i,j,k,GDV))/dx(2) + &
                  HALF*(q3(i,j,k+1,GDPRES) + q3(i,j,k,GDPRES)) * &
                  (q3(i,j,k+1,GDW) - q3(i,j,k,GDW))/dx(3)

          enddo
       enddo
    enddo

  end subroutine calc_pdivu

  

  subroutine ctu_consup(lo, hi, &
                        uin, uin_lo, uin_hi, &
                        q, q_lo, q_hi, &
                        update, updt_lo, updt_hi, &
                        flux1, flux1_lo, flux1_hi, &
                        flux2, flux2_lo, flux2_hi, &
                        flux3, flux3_lo, flux3_hi, &
                        qx, qx_lo, qx_hi, &
                        qy, qy_lo, qy_hi, &
                        qz, qz_lo, qz_hi, &
                        area1, area1_lo, area1_hi, &
                        area2, area2_lo, area2_hi, &
                        area3, area3_lo, area3_hi, &
                        vol, vol_lo, vol_hi, &
                        pdivu, pdivu_lo, pdivu_hi, &
                        dx, dt) bind(C, name="ctu_consup")

    use castro_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, &
                             UEINT, UTEMP, NGDNV, QVAR, GDPRES
    use amrex_constants_module, only: ZERO, ONE, TWO, FOURTH, HALF

    integer, intent(in) ::       lo(3),       hi(3)
    integer, intent(in) ::   uin_lo(3),   uin_hi(3)
    integer, intent(in) ::     q_lo(3),     q_hi(3)
    integer, intent(in) ::  updt_lo(3),  updt_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) ::    qy_lo(3),    qy_hi(3)
    integer, intent(in) :: flux3_lo(3), flux3_hi(3)
    integer, intent(in) :: area3_lo(3), area3_hi(3)
    integer, intent(in) ::    qz_lo(3),    qz_hi(3)
    integer, intent(in) ::    qx_lo(3),    qx_hi(3)
    integer, intent(in) ::   vol_lo(3),   vol_hi(3)
    integer, intent(in) ::   pdivu_lo(3),   pdivu_hi(3)

    real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),updt_lo(3):updt_hi(3),NVAR)
    real(rt), intent(in) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt), intent(in) ::    qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    real(rt), intent(in) ::    qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt), intent(in) :: flux3(flux3_lo(1):flux3_hi(1),flux3_lo(2):flux3_hi(2),flux3_lo(3):flux3_hi(3),NVAR)
    real(rt), intent(in) :: area3(area3_lo(1):area3_hi(1),area3_lo(2):area3_hi(2),area3_lo(3):area3_hi(3))
    real(rt), intent(in) ::    qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: pdivu(pdivu_lo(1):pdivu_hi(1),pdivu_lo(2):pdivu_hi(2),pdivu_lo(3):pdivu_hi(3))
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in), value :: dt

    integer :: i, j, g, k, n
    integer :: domlo(3), domhi(3)
    real(rt) :: volInv

    call calc_pdivu(lo, hi, &
                    qx, qx_lo, qx_hi, &
                    area1, area1_lo, area1_hi, &
                    qy, qy_lo, qy_hi, &
                    area2, area2_lo, area2_hi, &
                    qz, qz_lo, qz_hi, &
                    area3, area3_lo, area3_hi, &
                    vol, vol_lo, vol_hi, &
                    dx, pdivu, pdivu_lo, pdivu_hi)

    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update.

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                volinv = ONE / vol(i,j,k)

                update(i,j,k,n) = update(i,j,k,n) + &
                     ( flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) &
                     + flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) &
                     + flux3(i,j,k,n) * area3(i,j,k) - flux3(i,j,k+1,n) * area3(i,j,k+1) &
                     ) * volinv

                ! Add the p div(u) source term to (rho e).
                if (n .eq. UEINT) then
                   update(i,j,k,n) = update(i,j,k,n) - pdivu(i,j,k)
                endif

             enddo
          enddo
       enddo
    enddo

  end subroutine ctu_consup


  subroutine scale_flux(lo, hi, &
                        flux, f_lo, f_hi, &
                        area, a_lo, a_hi, dt) bind(c, name="scale_flux")

    use castro_module, only: NVAR, GDPRES, UMX, NGDNV

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    real(rt), intent(in   ), value :: dt

    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3))

    integer :: i, j, k, n

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                flux(i,j,k,n) = dt * flux(i,j,k,n) * area(i,j,k)
             enddo
          enddo
       enddo
    enddo

  end subroutine scale_flux

  CASTRO_FORT_DEVICE subroutine ca_construct_flux(lo, hi, domlo, domhi, dx, dt, &
                                                  idir, stage_weight, &
                                                  uin, uin_lo, uin_hi, &
                                                  div, div_lo, div_hi, &
                                                  qaux, qa_lo, qa_hi, &
                                                  qm, qm_lo, qm_hi, &
                                                  qp, qp_lo, qp_hi, &
                                                  qint, qe_lo, qe_hi, &
                                                  flux, f_lo, f_hi, &
                                                  fluxes, fs_lo, fs_hi, &
                                                  area, a_lo, a_hi) &
                                                  bind(C, name='ca_construct_flux')

    use amrex_fort_module, only: rt => amrex_real
    use castro_module, only: NVAR, QVAR, NGDNV, NQAUX
    use riemann_module, only: cmpflx

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ), value :: idir
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)
    integer,  intent(in   ) :: qe_lo(3), qe_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: fs_lo(3), fs_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)
    real(rt), intent(in   ), value :: stage_weight

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1), div_lo(2):div_hi(2), div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),QVAR,3)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),QVAR,3)
    real(rt), intent(inout) :: qint(qe_lo(1):qe_hi(1), qe_lo(2):qe_hi(2), qe_lo(3):qe_hi(3), NGDNV)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NVAR)
    real(rt), intent(inout) :: fluxes(fs_lo(1):fs_hi(1), fs_lo(2):fs_hi(2), fs_lo(3):fs_hi(3), NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt

    call cmpflx(lo, hi, domlo, domhi, idir, qm, qm_lo, qm_hi, qp, qp_lo, qp_hi, &
                qint, qe_lo, qe_hi, flux, f_lo, f_hi, qaux, qa_lo, qa_hi)
!    call apply_av(lo, hi, idir, dx, div, div_lo, div_hi, uin, uin_lo, uin_hi, flux, f_lo, f_hi)
!    call normalize_species_fluxes(lo, hi, flux, f_lo, f_hi)
!    call scale_flux(lo, hi, flux, f_lo, f_hi, area, a_lo, a_hi, dt)
    call update_fluxes(lo, hi, fluxes, fs_lo, fs_hi, flux, f_lo, f_hi, stage_weight)

  end subroutine ca_construct_flux



  CASTRO_FORT_DEVICE subroutine ca_ctoprim(lo, hi, &
                                           uin, uin_lo, uin_hi, &
                                           q,     q_lo,   q_hi, &
                                           qaux, qa_lo,  qa_hi) &
                                           bind(C, name='ca_ctoprim')

    use network, only: nspec
    use eos_module, only: eos_t, eos_input_re, eos
    use amrex_constants_module, only: ZERO, HALF, ONE
    use castro_module, only: NVAR, URHO, UMX, UMZ, &
                             UEDEN, UEINT, UTEMP, UFS, &
                             QRHO, QU, QV, QW, &
                             QREINT, QPRES, QTEMP, QGAME, QFS, &
                             QVAR, QC, QGAMC, QDPDR, QDPDE, NQAUX, &
                             small_dens

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), parameter :: small = 1.e-8_rt
    real(rt), parameter :: dual_energy_eta1 = 1.e0_rt

    integer  :: i, j, k, g
    integer  :: n, iq, ispec
    real(rt) :: kineng, rhoinv
    real(rt) :: vel(3)

    type (eos_t) :: eos_state

    !$acc parallel loop gang vector collapse(3) deviceptr(uin, q, qaux) private(vel, eos_state) async(acc_stream)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             q(i,j,k,QRHO) = uin(i,j,k,URHO)
             rhoinv = ONE/q(i,j,k,QRHO)

             vel = uin(i,j,k,UMX:UMZ) * rhoinv

             q(i,j,k,QU:QW) = vel

             ! Get the internal energy, which we'll use for
             ! determining the pressure.  We use a dual energy
             ! formalism. If (E - K) < eta1 and eta1 is suitably
             ! small, then we risk serious numerical truncation error
             ! in the internal energy.  Therefore we'll use the result
             ! of the separately updated internal energy equation.
             ! Otherwise, we'll set e = E - K.

             kineng = HALF * q(i,j,k,QRHO) * (q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)

             if ( (uin(i,j,k,UEDEN) - kineng) / uin(i,j,k,UEDEN) .gt. dual_energy_eta1) then
                q(i,j,k,QREINT) = (uin(i,j,k,UEDEN) - kineng) * rhoinv
             else
                q(i,j,k,QREINT) = uin(i,j,k,UEINT) * rhoinv
             endif

             ! If we're advecting in the rotating reference frame,
             ! then subtract off the rotation component here.

             q(i,j,k,QTEMP) = uin(i,j,k,UTEMP)

             ! Load passively advected quantities into q

             do ispec = 1, nspec
                n  = UFS + ispec - 1
                iq = QFS + ispec - 1
                q(i,j,k,iq) = uin(i,j,k,n)/q(i,j,k,QRHO)
             enddo

             ! get gamc, p, T, c, csml using q state

             eos_state % T   = q(i,j,k,QTEMP )
             eos_state % rho = q(i,j,k,QRHO  )
             eos_state % e   = q(i,j,k,QREINT)
             eos_state % xn  = q(i,j,k,QFS:QFS+nspec-1)

             call eos(eos_input_re, eos_state)

             q(i,j,k,QTEMP)  = eos_state % T
             q(i,j,k,QREINT) = eos_state % e * q(i,j,k,QRHO)
             q(i,j,k,QPRES)  = eos_state % p
             q(i,j,k,QGAME)  = q(i,j,k,QPRES) / q(i,j,k,QREINT) + ONE

             qaux(i,j,k,QDPDR)  = eos_state % dpdr_e
             qaux(i,j,k,QDPDE)  = eos_state % dpde

             qaux(i,j,k,QGAMC)  = eos_state % gam1
             qaux(i,j,k,QC   )  = eos_state % cs

          enddo
       enddo
    enddo

  end subroutine ca_ctoprim



  CASTRO_FORT_DEVICE subroutine ca_divu(lo, hi, dx, q, q_lo, q_hi, div, d_lo, d_hi) bind(C, name='ca_divu')

    use amrex_constants_module, only: FOURTH, ONE
    use castro_module, only: QU, QV, QW, QVAR

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: d_lo(3), d_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),QVAR)
    real(rt), intent(inout) :: div(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3))

    integer  :: i, j, k
    real(rt) :: ux, vy, wz, dxinv, dyinv, dzinv

    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    !$acc parallel loop gang vector collapse(3) deviceptr(q, div) async(acc_stream)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ux = FOURTH*( &
                    + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                    + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                    + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                    + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) ) * dxinv

             vy = FOURTH*( &
                    + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                    + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                    + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                    + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) ) * dyinv

             wz = FOURTH*( &
                    + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                    + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                    + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                    + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) ) * dzinv

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo

  end subroutine ca_divu



  CASTRO_FORT_DEVICE subroutine ca_construct_hydro_update(lo, hi, dx, dt, stage_weight, &
                                                          q1, q1_lo, q1_hi, &
                                                          q2, q2_lo, q2_hi, &
                                                          q3, q3_lo, q3_hi, &
                                                          f1, f1_lo, f1_hi, &
                                                          f2, f2_lo, f2_hi, &
                                                          f3, f3_lo, f3_hi, &
                                                          a1, a1_lo, a1_hi, &
                                                          a2, a2_lo, a2_hi, &
                                                          a3, a3_lo, a3_hi, &
                                                          vol, vol_lo, vol_hi, &
                                                          update, u_lo, u_hi) &
                                                          bind(C, name='ca_construct_hydro_update')

    use amrex_constants_module, only: HALF, ONE
    use castro_module, only: NVAR, UEINT, NGDNV, GDPRES, GDU, GDV, GDW

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q1_lo(3), q1_hi(3)
    integer,  intent(in   ) :: q2_lo(3), q2_hi(3)
    integer,  intent(in   ) :: q3_lo(3), q3_hi(3)
    integer,  intent(in   ) :: f1_lo(3), f1_hi(3)
    integer,  intent(in   ) :: f2_lo(3), f2_hi(3)
    integer,  intent(in   ) :: f3_lo(3), f3_hi(3)
    integer,  intent(in   ) :: a1_lo(3), a1_hi(3)
    integer,  intent(in   ) :: a2_lo(3), a2_hi(3)
    integer,  intent(in   ) :: a3_lo(3), a3_hi(3)
    integer,  intent(in   ) :: vol_lo(3), vol_hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, stage_weight

    real(rt), intent(in   ) :: q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
    real(rt), intent(in   ) :: q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
    real(rt), intent(in   ) :: q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)
    real(rt), intent(in   ) :: f1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3),NVAR)
    real(rt), intent(in   ) :: f2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3),NVAR)
    real(rt), intent(in   ) :: f3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3),NVAR)
    real(rt), intent(in   ) :: a1(a1_lo(1):a1_hi(1),a1_lo(2):a1_hi(2),a1_lo(3):a1_hi(3))
    real(rt), intent(in   ) :: a2(a2_lo(1):a2_hi(1),a2_lo(2):a2_hi(2),a2_lo(3):a2_hi(3))
    real(rt), intent(in   ) :: a3(a3_lo(1):a3_hi(1),a3_lo(2):a3_hi(2),a3_lo(3):a3_hi(3))
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: update(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

    integer  :: i, j, k, n
    real(rt) :: pdivu, dxinv(3), dtinv

    dtinv = ONE / dt
    dxinv = ONE / dx

    !$acc parallel loop gang vector collapse(4) deviceptr(q1, q2, q3, f1, f2, f3, a1, a2, a3, vol, update) async(acc_stream)
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Note that the fluxes have already been scaled by dt * dA.

                update(i,j,k,n) = update(i,j,k,n) + stage_weight * dtinv * (f1(i,j,k,n) - f1(i+1,j,k,n) + &
                                                                            f2(i,j,k,n) - f2(i,j+1,k,n) + &
                                                                            f3(i,j,k,n) - f3(i,j,k+1,n) ) / vol(i,j,k)

                ! Add the p div(u) source term to (rho e).
                if (n .eq. UEINT) then

                   pdivu = HALF * (q1(i+1,j,k,GDPRES) + q1(i,j,k,GDPRES)) * &
                                  (q1(i+1,j,k,GDU) - q1(i,j,k,GDU)) * dxinv(1) + &
                           HALF * (q2(i,j+1,k,GDPRES) + q2(i,j,k,GDPRES)) * &
                                  (q2(i,j+1,k,GDV) - q2(i,j,k,GDV)) * dxinv(2) + &
                           HALF * (q3(i,j,k+1,GDPRES) + q3(i,j,k,GDPRES)) * &
                                  (q3(i,j,k+1,GDW) - q3(i,j,k,GDW)) * dxinv(3)

                   update(i,j,k,n) = update(i,j,k,n) - stage_weight * pdivu

                endif

             enddo
          enddo
       enddo
    enddo

  end subroutine ca_construct_hydro_update



  CASTRO_FORT_DEVICE subroutine update_fluxes(lo, hi, fluxes, fs_lo, fs_hi, flux, f_lo, f_hi, stage_weight)

    use castro_module, only: NVAR

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: fs_lo(3), fs_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    real(rt), intent(in   ), value :: stage_weight

    real(rt), intent(inout) :: fluxes(fs_lo(1):fs_hi(1),fs_lo(2):fs_hi(2),fs_lo(3):fs_hi(3),NVAR)
    real(rt), intent(in   ) :: flux(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3),NVAR)

    integer :: i, j, k, n

    ! Update the total fluxes for the timestep with the contribution from this stage.

    !$acc parallel loop gang vector collapse(4) deviceptr(fluxes, flux) async(acc_stream)
    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                fluxes(i,j,k,n) = fluxes(i,j,k,n) + stage_weight * flux(i,j,k,n)
             enddo
          enddo
       enddo
    enddo

  end subroutine update_fluxes

end module hydro_module
