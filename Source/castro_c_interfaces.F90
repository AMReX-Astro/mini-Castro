module c_interface_modules

  use meth_params_module, only: NVAR, NQAUX, NQ, QVAR, NGDNV
  use amrex_fort_module, only: rt => amrex_real
#ifdef AMREX_USE_CUDA
  use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

contains

  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                        adv_h3,domlo,domhi,dx,xlo,time,bc) bind(C, name="ca_hypfill")

    use bc_fill_module, only: hypfill

    implicit none

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3
    integer,  intent(in   ) :: bc(3,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: dx(3), xlo(3), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

#ifdef AMREX_USE_CUDA
    attributes(device) :: adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, bc, dx, xlo, time, domlo, domhi
#endif

    call hypfill &
#ifdef AMREX_USE_CUDA
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

#ifdef AMREX_USE_CUDA
    attributes(device) :: adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, bc, dx, xlo, time, domlo, domhi
#endif

    call denfill &
#ifdef AMREX_USE_CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (adv, adv_l1, adv_l2, adv_l3, adv_h1, adv_h2, adv_h3, domlo, domhi, dx, xlo, time, bc)

  end subroutine ca_denfill

end module c_interface_modules
