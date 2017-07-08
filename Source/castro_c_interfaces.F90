module c_interface_modules

  use meth_params_module, only: NVAR, NQAUX, NQ, QVAR, NGDNV
  use amrex_fort_module, only: rt => amrex_real
#ifdef CUDA
  use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

contains

  subroutine ca_initdata(level, lo, hi, &
                         state, s_lo, s_hi, dx, &
                         xlo, xhi) bind(C, name="ca_initdata")

    use initdata_module, only: initdata

    implicit none

    integer,  intent(in   ), value :: level
    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(in   ) :: xlo(3), xhi(3), dx(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

#ifdef CUDA
    attributes(device) :: lo, hi, s_lo, s_hi, xlo, xhi, dx, state
#endif

    call initdata &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (level, lo, hi, state, s_lo, s_hi, dx, xlo, xhi)

  end subroutine ca_initdata


  subroutine ca_enforce_consistent_e(lo,hi,state,s_lo,s_hi) &
                                     bind(c, name='ca_enforce_consistent_e')

    use castro_util_module, only: enforce_consistent_e

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: state, lo, hi, s_lo, s_hi
#endif

    call enforce_consistent_e &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, state, s_lo, s_hi)

  end subroutine ca_enforce_consistent_e


  subroutine ca_compute_temp(lo, hi, state, s_lo, s_hi) &
                             bind(C, name="ca_compute_temp")

    use castro_util_module, only: compute_temp

    implicit none

    integer,  intent(in   ) :: lo(3),hi(3)
    integer,  intent(in   ) :: s_lo(3),s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: state, lo, hi, s_lo, s_hi
#endif

    call compute_temp &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, state, s_lo, s_hi)

  end subroutine ca_compute_temp


  subroutine ca_reset_internal_e(lo, hi, u, u_lo, u_hi, verbose) &
                                 bind(C, name="ca_reset_internal_e")

    use castro_util_module, only: reset_internal_e

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ), value :: verbose
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: u, lo, hi, u_lo, u_hi
#endif

    call reset_internal_e &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, u, u_lo, u_hi, verbose)

  end subroutine ca_reset_internal_e


  subroutine ca_normalize_species(u, u_lo, u_hi, lo, hi) &
                                  bind(C, name="ca_normalize_species")

    use castro_util_module, only: normalize_species

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(inout) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: u, lo, hi, u_lo, u_hi
#endif

    call normalize_species &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (u, u_lo, u_hi, lo, hi)

  end subroutine ca_normalize_species


  subroutine ca_enforce_minimum_density(uin, uin_lo, uin_hi, &
                                        uout, uout_lo, uout_hi, &
                                        vol, vol_lo, vol_hi, &
                                        lo, hi, frac_change, verbose) &
                                        bind(C, name="ca_enforce_minimum_density")

    use advection_util_module, only: enforce_minimum_density

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ), value :: verbose
    integer,  intent(in   ) ::  uin_lo(3),  uin_hi(3)
    integer,  intent(in   ) :: uout_lo(3), uout_hi(3)
    integer,  intent(in   ) ::  vol_lo(3),  vol_hi(3)

    real(rt), intent(in   ) ::  uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2), uin_lo(3): uin_hi(3),NVAR)
    real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt), intent(in   ) ::  vol( vol_lo(1): vol_hi(1), vol_lo(2): vol_hi(2), vol_lo(3): vol_hi(3))
    real(rt), intent(inout) :: frac_change

#ifdef CUDA
    attributes(managed) :: uin, uout, vol, lo, hi, uin_lo, uin_hi, uout_lo, uout_hi, vol_lo, vol_hi, frac_change
#endif

    call enforce_minimum_density &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (uin, uin_lo, uin_hi, &
          uout, uout_lo, uout_hi, &
          vol, vol_lo, vol_hi, &
          lo, hi, frac_change, verbose)

  end subroutine ca_enforce_minimum_density


  subroutine ca_check_initial_species(lo, hi, state, state_lo, state_hi) &
                                      bind(C, name="ca_check_initial_species")

    use castro_util_module, only: check_initial_species

    implicit none

    integer,  intent(in) :: lo(3), hi(3)
    integer,  intent(in) :: state_lo(3), state_hi(3)
    real(rt), intent(in) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

#ifdef CUDA
    attributes(managed) :: lo, hi, state, state_lo, state_hi
#endif

    call check_initial_species &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, state, state_lo, state_hi)

  end subroutine ca_check_initial_species


  subroutine ca_ctoprim(lo, hi, &
                        uin, uin_lo, uin_hi, &
                        q,     q_lo,   q_hi, &
                        qaux, qa_lo,  qa_hi) bind(C, name = "ca_ctoprim")

    use advection_util_module, only: ctoprim

    implicit none

    integer,  intent(in) :: lo(3), hi(3)
    integer,  intent(in) :: uin_lo(3), uin_hi(3)
    integer,  intent(in) :: q_lo(3), q_hi(3)
    integer,  intent(in) :: qa_lo(3), qa_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

#ifdef CUDA
    attributes(managed) :: uin, q, qaux, lo, hi, uin_lo, uin_hi, q_lo, q_hi, qa_lo, qa_hi
#endif

    call ctoprim &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, uin, uin_lo, uin_hi, q, q_lo, q_hi, qaux, qa_lo, qa_hi)

  end subroutine ca_ctoprim


  subroutine ca_dervel(vel,v_lo,v_hi,nv, &
                       dat,d_lo,d_hi,nc,lo,hi,domlo, &
                       domhi,dx,xlo,time,dt,bc,level,grid_no) &
                       bind(c, name="ca_dervel")

    use castro_util_module, only: dervel

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3), nv
    integer,  intent(in   ) :: d_lo(3), d_hi(3), nc
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(3,2,nc)
    real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
    real(rt), intent(inout) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3),nv)
    real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer,  intent(in   ) :: level, grid_no

#ifdef CUDA
    attributes(device) :: vel, dat, lo, hi, v_lo, v_hi, d_lo, d_hi, dx, xlo, domlo, domhi, nv, nc, bc, time, dt, level, grid_no
#endif

    call dervel &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (vel,v_lo,v_hi,nv, &
          dat,d_lo,d_hi,nc, &
          lo,hi,domlo,domhi, &
          dx,xlo,time,dt, &
          bc,level,grid_no)

  end subroutine ca_dervel



  subroutine ca_derpres(p,p_lo,p_hi,np, &
                        u,u_lo,u_hi,nc,lo,hi,domlo, &
                        domhi,dx,xlo,time,dt,bc,level,grid_no) &
                        bind(c, name="ca_derpres")

    use castro_util_module, only: derpres

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: p_lo(3), p_hi(3), np
    integer,  intent(in   ) :: u_lo(3), u_hi(3), nc
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),np)
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),nc)
    real(rt), intent(in   ) :: dx(3), xlo(3), time, dt
    integer,  intent(in   ) :: bc(3,2,nc), level, grid_no

#ifdef CUDA
    attributes(managed) :: p, u, lo, hi, p_lo, p_hi, u_lo, u_hi, dx, xlo, domlo, domhi, np, nc, bc, time, dt, level, grid_no
#endif

    call derpres &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (p,p_lo,p_hi,np, &
          u,u_lo,u_hi,nc, &
          lo,hi,domlo,domhi, &
          dx,xlo,time,dt, &
          bc,level,grid_no)

  end subroutine ca_derpres

  subroutine ca_estdt(lo,hi,u,u_lo,u_hi,dx,dt) bind(C, name="ca_estdt")

    use timestep_module, only: estdt

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(inout) :: dt

#ifdef CUDA
    attributes(device) :: u, lo, hi, u_lo, u_hi, dx, dt
#endif

    call estdt &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, u, u_lo, u_hi, dx, dt)

  end subroutine ca_estdt




  subroutine ca_prepare_for_fluxes(lo, hi, &
                                   dx, dt, &
                                   q, q_lo, q_hi, &
                                   qaux, qa_lo, qa_hi, &
                                   flatn, flatn_lo, flatn_hi, &
                                   div, div_lo, div_hi, &
                                   sxm, sxm_lo, sxm_hi, &
                                   sxp, sxp_lo, sxp_hi, &
                                   sym, sym_lo, sym_hi, &
                                   syp, syp_lo, syp_hi, &
                                   szm, szm_lo, szm_hi, &
                                   szp, szp_lo, szp_hi) bind(C, name="ca_prepare_for_fluxes")

    use mol_module, only: prepare_for_fluxes

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: flatn_lo(3), flatn_hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: sxm_lo(3), sxm_hi(3)
    integer,  intent(in   ) :: sxp_lo(3), sxp_hi(3)
    integer,  intent(in   ) :: sym_lo(3), sym_hi(3)
    integer,  intent(in   ) :: syp_lo(3), syp_hi(3)
    integer,  intent(in   ) :: szm_lo(3), szm_hi(3)
    integer,  intent(in   ) :: szp_lo(3), szp_hi(3)

    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: flatn(flatn_lo(1):flatn_hi(1), flatn_lo(2):flatn_hi(2), flatn_lo(3):flatn_hi(3), NQ)
    real(rt), intent(inout) :: div(div_lo(1):div_hi(1), div_lo(2):div_hi(2), div_lo(3):div_hi(3))
    real(rt), intent(inout) :: sxm(sxm_lo(1):sxm_hi(1), sxm_lo(2):sxm_hi(2), sxm_lo(3):sxm_hi(3), NQ)
    real(rt), intent(inout) :: sxp(sxp_lo(1):sxp_hi(1), sxp_lo(2):sxp_hi(2), sxp_lo(3):sxp_hi(3), NQ)
    real(rt), intent(inout) :: sym(sym_lo(1):sym_hi(1), sym_lo(2):sym_hi(2), sym_lo(3):sym_hi(3), NQ)
    real(rt), intent(inout) :: syp(syp_lo(1):syp_hi(1), syp_lo(2):syp_hi(2), syp_lo(3):syp_hi(3), NQ)
    real(rt), intent(inout) :: szm(szm_lo(1):szm_hi(1), szm_lo(2):szm_hi(2), szm_lo(3):szm_hi(3), NQ)
    real(rt), intent(inout) :: szp(szp_lo(1):szp_hi(1), szp_lo(2):szp_hi(2), szp_lo(3):szp_hi(3), NQ)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt

#ifdef CUDA
    attributes(managed) :: lo, hi, dx, q, q_lo, q_hi, qaux, qa_lo, qa_hi, &
                           flatn, flatn_lo, flatn_hi, div, div_lo, div_hi, &
                           sxm, sxm_lo, sxm_hi, sxp, sxp_lo, sxp_hi, &
                           sym, sym_lo, sym_hi, syp, syp_lo, syp_hi, &
                           szm, szm_lo, szm_hi, szp, szp_lo, szp_hi
#endif

    ! Construct edge states as inputs to flux construction

    call prepare_for_fluxes &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dt, dx, &
          q, q_lo, q_hi, &
          qaux, qa_lo, qa_hi, &
          flatn, flatn_lo, flatn_hi, &
          div, div_lo, div_hi, &
          sxm, sxp, sym, syp, szm, szp, sxm_lo, sxm_hi)

  end subroutine ca_prepare_for_fluxes



  subroutine ca_prepare_profile(lo, hi, &
                                q, q_lo, q_hi, &
                                qm, qm_lo, qp_hi, &
                                qp, qp_lo, qm_hi, &
                                sxm, sxm_lo, sxm_hi, &
                                sxp, sxp_lo, sxp_hi, &
                                sym, sym_lo, sym_hi, &
                                syp, syp_lo, syp_hi, &
                                szm, szm_lo, szm_hi, &
                                szp, szp_lo, szp_hi) bind(C, name="ca_prepare_profile")

    use mol_module, only: prepare_profile

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)
    integer,  intent(in   ) :: sxm_lo(3), sxm_hi(3)
    integer,  intent(in   ) :: sxp_lo(3), sxp_hi(3)
    integer,  intent(in   ) :: sym_lo(3), sym_hi(3)
    integer,  intent(in   ) :: syp_lo(3), syp_hi(3)
    integer,  intent(in   ) :: szm_lo(3), szm_hi(3)
    integer,  intent(in   ) :: szp_lo(3), szp_hi(3)

    real(rt), intent(in   ) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1), qm_lo(2):qm_hi(2), qm_lo(3):qm_hi(3), NQ, 3)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1), qp_lo(2):qp_hi(2), qp_lo(3):qp_hi(3), NQ, 3)
    real(rt), intent(in   ) :: sxm(sxm_lo(1):sxm_hi(1), sxm_lo(2):sxm_hi(2), sxm_lo(3):sxm_hi(3), NQ)
    real(rt), intent(in   ) :: sxp(sxp_lo(1):sxp_hi(1), sxp_lo(2):sxp_hi(2), sxp_lo(3):sxp_hi(3), NQ)
    real(rt), intent(in   ) :: sym(sym_lo(1):sym_hi(1), sym_lo(2):sym_hi(2), sym_lo(3):sym_hi(3), NQ)
    real(rt), intent(in   ) :: syp(syp_lo(1):syp_hi(1), syp_lo(2):syp_hi(2), syp_lo(3):syp_hi(3), NQ)
    real(rt), intent(in   ) :: szm(szm_lo(1):szm_hi(1), szm_lo(2):szm_hi(2), szm_lo(3):szm_hi(3), NQ)
    real(rt), intent(in   ) :: szp(szp_lo(1):szp_hi(1), szp_lo(2):szp_hi(2), szp_lo(3):szp_hi(3), NQ)

#ifdef CUDA
    attributes(managed) :: lo, hi, q, q_lo, q_hi, &
                           qm, qm_lo, qm_hi, qp, qp_lo, qp_hi, &
                           sxm, sxm_lo, sxm_hi, sxp, sxp_lo, sxp_hi, &
                           sym, sym_lo, sym_hi, syp, syp_lo, syp_hi, &
                           szm, szm_lo, szm_hi, szp, szp_lo, szp_hi
#endif

    call prepare_profile &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, &
          q, q_lo, q_hi, &
          sxm, sxp, sym, syp, szm, szp, sxm_lo, sxm_hi, &
          qm, qp, qm_lo, qm_hi)

  end subroutine ca_prepare_profile



  subroutine ca_construct_flux(lo, hi, &
                               domlo, domhi, &
                               dx, dt, &
                               idir, &
                               uin, uin_lo, uin_hi, &
                               div, div_lo, div_hi, &
                               qaux, qa_lo, qa_hi, &
                               qm, qm_lo, qp_hi, &
                               qp, qp_lo, qm_hi, &
                               qe, qe_lo, qe_hi, &
                               flux, f_lo, f_hi, &
                               area, a_lo, a_hi) bind(C, name="ca_construct_flux")

    use mol_module, only: construct_flux

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ), value :: idir
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)
    integer,  intent(in   ) :: qe_lo(3), qe_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1), div_lo(2):div_hi(2), div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(in   ) :: qm(qm_lo(1):qm_hi(1), qm_lo(2):qm_hi(2), qm_lo(3):qm_hi(3), NQ, 3)
    real(rt), intent(in   ) :: qp(qp_lo(1):qp_hi(1), qp_lo(2):qp_hi(2), qp_lo(3):qp_hi(3), NQ, 3)
    real(rt), intent(inout) :: qe(qe_lo(1):qe_hi(1), qe_lo(2):qe_hi(2), qe_lo(3):qe_hi(3), NGDNV)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt

#ifdef CUDA
    attributes(managed) :: lo, hi, domlo, domhi, dx, &
                           uin, uin_lo, uin_hi, &
                           div, div_lo, div_hi, &
                           qaux, qa_lo, qa_hi, &
                           qm, qm_lo, qm_hi, &
                           qp, qp_lo, qp_hi, &
                           qe, qe_lo, qe_hi, &
                           flux, f_lo, f_hi, &
                           area, a_lo, a_hi
#endif

    call construct_flux &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, domlo, domhi, dx, dt, idir, &
          div, div_lo, div_hi, &
          uin, uin_lo, uin_hi, &
          qm, qp, qm_lo, qm_hi, &
          flux, qe, f_lo, f_hi, &
          area, a_lo, a_hi, &
          qaux, qa_lo, qa_hi)

  end subroutine ca_construct_flux



  subroutine ca_construct_hydro_update(lo, hi, &
                                       dx, dt, &
                                       stage_weight, &
                                       q1, q1_lo, q1_hi, &
                                       q2, q2_lo, q2_hi, &
                                       q3, q3_lo, q3_hi, &
                                       flux1, f1_lo, f1_hi, &
                                       flux2, f2_lo, f2_hi, &
                                       flux3, f3_lo, f3_hi, &
                                       area1, a1_lo, a1_hi, &
                                       area2, a2_lo, a2_hi, &
                                       area3, a3_lo, a3_hi, &
                                       vol, vol_lo, vol_hi, &
                                       update, updt_lo, updt_hi) bind(C, name="ca_construct_hydro_update")

    use advection_util_module, only: construct_hydro_update

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
    integer,  intent(in   ) :: updt_lo(3), updt_hi(3)

    real(rt), intent(in   ) :: q1(q1_lo(1):q1_hi(1), q1_lo(2):q1_hi(2), q1_lo(3):q1_hi(3), NGDNV)
    real(rt), intent(in   ) :: q2(q2_lo(1):q2_hi(1), q2_lo(2):q2_hi(2), q2_lo(3):q2_hi(3), NGDNV)
    real(rt), intent(in   ) :: q3(q3_lo(1):q3_hi(1), q3_lo(2):q3_hi(2), q3_lo(3):q3_hi(3), NGDNV)
    real(rt), intent(in   ) :: flux1(f1_lo(1):f1_hi(1), f1_lo(2):f1_hi(2), f1_lo(3):f1_hi(3), NVAR)
    real(rt), intent(in   ) :: flux2(f2_lo(1):f2_hi(1), f2_lo(2):f2_hi(2), f2_lo(3):f2_hi(3), NVAR)
    real(rt), intent(in   ) :: flux3(f3_lo(1):f3_hi(1), f3_lo(2):f3_hi(2), f3_lo(3):f3_hi(3), NVAR)
    real(rt), intent(in   ) :: area1(a1_lo(1):a1_hi(1), a1_lo(2):a1_hi(2), a1_lo(3):a1_hi(3))
    real(rt), intent(in   ) :: area2(a2_lo(1):a2_hi(1), a2_lo(2):a2_hi(2), a2_lo(3):a2_hi(3))
    real(rt), intent(in   ) :: area3(a3_lo(1):a3_hi(1), a3_lo(2):a3_hi(2), a3_lo(3):a3_hi(3))
    real(rt), intent(in   ) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
    real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt, stage_weight

#ifdef CUDA
    attributes(managed) :: lo, hi, dx, &
                           q1, q1_lo, q1_hi, &
                           q2, q2_lo, q2_hi, &
                           q3, q3_lo, q3_hi, &
                           flux1, f1_lo, f1_hi, &
                           flux2, f2_lo, f2_hi, &
                           flux3, f3_lo, f3_hi, &
                           area1, a1_lo, a1_hi, &
                           area2, a2_lo, a2_hi, &
                           area3, a3_lo, a3_hi, &
                           vol, vol_lo, vol_hi, &
                           update, updt_lo, updt_hi
#endif

    ! Create an update source term based on the flux divergence.

    call construct_hydro_update &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dx, dt, stage_weight, &
          flux1, q1, f1_lo, f1_hi, &
          flux2, q2, f2_lo, f2_hi, &
          flux3, q3, f3_lo, f3_hi, &
          area1, a1_lo, a1_hi, &
          area2, a2_lo, a2_hi, &
          area3, a3_lo, a3_hi, &
          vol, vol_lo, vol_hi, &
          update, updt_lo, updt_hi)

  end subroutine ca_construct_hydro_update



  subroutine ca_summass(lo,hi,rho,r_lo,r_hi,dx, &
                        vol,v_lo,v_hi,mass) bind(c, name='ca_summass')

    use castro_util_module, only: summass

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: v_lo(3), v_hi(3)
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ) :: rho(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: vol(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    real(rt), intent(inout) :: mass

#ifdef CUDA
    attributes(device) :: rho, vol, lo, hi, r_lo, r_hi, v_lo, v_hi, dx, mass
#endif

    call summass &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, rho, r_lo, r_hi, dx, vol, v_lo, v_hi, mass)

  end subroutine ca_summass



  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,dx,xlo,time,bc) bind(C, name="ca_hypfill")

    use bc_fill_module, only: hypfill

    implicit none

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

#ifdef CUDA
    attributes(device) :: adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, bc, dx, xlo, time, domlo, domhi
#endif

    call hypfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, domlo, domhi, dx, xlo, time, bc)

  end subroutine ca_hypfill



  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,dx,xlo,time,bc) bind(C, name="ca_denfill")

    use amrex_fort_module, only: rt => amrex_real
    use bc_fill_module, only: denfill

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,1)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)

#ifdef CUDA
    attributes(device) :: adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, bc, dx, xlo, time, domlo, domhi
#endif

    call denfill &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, domlo, domhi, dx, xlo, time, bc)

  end subroutine ca_denfill

end module c_interface_modules
