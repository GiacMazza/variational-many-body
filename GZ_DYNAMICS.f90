MODULE GZ_DYNAMICS  
  USE GZ_VARS_GLOBAL
  USE GZ_MATRIX_BASIS
  USE GZ_EFFECTIVE_HOPPINGS
  USE GZ_neqAUX_FUNX
  implicit none
  private
  !
  public :: setup_neq_hamiltonian
  public :: setup_neq_dynamics
  public :: gz_neq_measure
  !  
  real(8),dimension(:,:),allocatable    :: gz_neq_local_density_matrix !
  real(8),dimension(:,:),allocatable    :: gz_neq_local_dens_dens
  real(8),dimension(4)                  :: gz_neq_local_angular_momenta
  real(8),dimension(3)                  :: gz_neq_energies
  complex(8),dimension(:,:),allocatable :: gz_neq_dens_constr_slater
  complex(8),dimension(:,:),allocatable :: gz_neq_dens_constr_gzproj
  real(8)                               :: gz_neq_unitary_constr
  complex(8),dimension(:,:),allocatable :: gz_neq_Rhop         
  !
CONTAINS
  !
  subroutine setup_neq_dynamics
    !
    ! myabe put here also other stufss about time grids etc etc..
    !
    allocate(gz_neq_local_density_matrix(Ns,Ns)); gz_neq_local_density_matrix = 0.d0
    allocate(gz_neq_local_dens_dens(Ns,Ns)); gz_neq_local_dens_dens = 0.d0
    gz_neq_local_angular_momenta = 0.d0
    gz_neq_energies = 0.d0
    allocate(gz_neq_dens_constr_slater(Ns,Ns)); gz_neq_dens_constr_slater = 0.d0
    allocate(gz_neq_dens_constr_gzproj(Ns,Ns)); gz_neq_dens_constr_gzproj = 0.d0
    gz_neq_unitary_constr = 0.d0
    allocate(gz_neq_Rhop(Ns,Ns)); gz_neq_Rhop = 0.d0
    !
  end subroutine setup_neq_dynamics

  subroutine get_neq_local_dens(is,js,x)
    integer :: is,js
    real(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_neq_local_dens wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_neq_local_dens wrong indeces'
    x = gz_neq_local_density_matrix(is,js)    
  end subroutine get_neq_local_dens
  !
  subroutine get_neq_local_dens_dens(is,js,x)
    integer :: is,js
    real(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_neq_local_dens_dens wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_neq_local_dens_dens wrong indeces'
    x = gz_neq_local_dens_dens(is,js)    
  end subroutine get_neq_local_dens_dens
  !
  subroutine get_neq_local_angular_momenta(x)
    real(8),dimension(4) :: x
    x = gz_neq_local_angular_momenta
  end subroutine get_neq_local_angular_momenta
  !
  subroutine get_neq_energies(x)
    real(8),dimension(3) :: x
    x = gz_neq_energies
  end subroutine get_neq_energies
  !
  subroutine get_neq_dens_constr_slater(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_neq_dens_constr_slater wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_neq_dens_constr_slater wrong indeces'
    x = gz_neq_dens_constr_slater(is,js)
  end subroutine get_neq_dens_constr_slater
  !
  subroutine get_neq_dens_constr_gzproj(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_neq_dens_constr_slater wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_neq_dens_constr_slater wrong indeces'
    x = gz_neq_dens_constr_gzproj(is,js)    
  end subroutine get_neq_dens_constr_gzproj
  !
  subroutine get_neq_unitary_constr(x)
    real(8) :: x
    x =gz_neq_unitary_constr
  end subroutine get_neq_unitary_constr
  !
  subroutine get_neq_Rhop(is,js,x)
    integer :: is,js
    complex(8) :: x
    x=gz_neq_Rhop(is,js)
  end subroutine get_neq_Rhop
  !
  subroutine setup_neq_hamiltonian(Uloc_t_,Ust_t_,Jh_t_,Jph_t_,Jsf_t_,eLevels_t_,Hk_tb_t_)
    real(8),dimension(3,Nt_aux),optional           :: Uloc_t_
    real(8),dimension(Nt_aux),optional             :: Ust_t_
    real(8),dimension(Nt_aux),optional             :: Jh_t_
    real(8),dimension(Nt_aux),optional             :: Jph_t_
    real(8),dimension(Nt_aux),optional             :: Jsf_t_
    real(8),dimension(Ns,Nt_aux),optional          :: eLevels_t_
    complex(8),dimension(Ns,Ns,Lk,Nt_aux),optional :: Hk_tb_t_
    !
    integer                                    :: it
    !
    allocate(Uloc_t(3,Nt_aux)); forall(it=1:Nt_aux) Uloc_t(:,it) = Uloc
    if(present(Uloc_t_)) Uloc_t = Uloc_t_
    !
    allocate(Ust_t(Nt_aux)); Ust_t = Ust
    if(present(Ust_t_)) Ust_t = Ust_t_
    !
    allocate(Jh_t(Nt_aux)); Jh_t = Jh
    if(present(Jh_t_)) Jh_t = Jh_t_
    !
    allocate(Jph_t(Nt_aux)); Jph_t = Jph
    if(present(Jph_t_)) Jph_t = Jph_t_
    !
    allocate(Jsf_t(Nt_aux)); Jsf_t = Jsf
    if(present(Jsf_t_)) Jsf_t = Jsf_t_
    !
    allocate(eLevels_t(Ns,Nt_aux)); forall(it=1:Nt_aux) eLevels_t(:,it) = eLevels
    if(present(eLevels_t_)) eLevels_t = eLevels_t_
    !
    allocate(Hk_tb_t(Ns,Ns,Lk,Nt_aux)); forall(it=1:Nt_aux) Hk_tb_t(:,:,:,it) = Hk_tb
    if(present(Hk_tb_t_)) Hk_tb_t = Hk_tb_t_
    !
  end subroutine setup_neq_hamiltonian
  

  
  subroutine gz_neq_measure(psi_t,time)
    complex(8),dimension(nDynamics) :: psi_t
    real(8)                         :: time
    complex(8),dimension(Ns,Ns,Lk)  :: slater
    complex(8),dimension(Nphi)      :: gzproj
    real(8)                         :: Estar,Eloc,Egz
    real(8),dimension(Ns)           :: vdm_diag
    
    integer                         :: is,js,ik,it,iphi

    it=0!qlcs
    !
    call dynamicalVector_2_wfMatrix(psi_t,slater,gzproj)  
    !
    Estar=0.d0
    do is=1,Ns
       do js=1,Ns
          gz_neq_local_density_matrix(is,js) = &
               trace_phi_basis(gzproj,phi_traces_basis_local_dens(is,js,:,:))          
          gz_neq_local_dens_dens(is,js) = &
               trace_phi_basis(gzproj,phi_traces_basis_dens_dens(is,js,:,:))
          !
          gz_neq_dens_constr_slater(is,js)=0.d0
          do ik=1,Lk
             gz_neq_dens_constr_slater(is,js) = gz_neq_dens_constr_slater(is,js) + &
                  slater(is,js,ik)*wtk(ik)
             Estar = Estar + Hk_tb_t(is,js,ik,it)*slater(is,js,ik)*wtk(ik)
          end do
          !
          gz_neq_dens_constr_gzproj(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
       end do
       vdm_diag(is) = gz_neq_dens_constr_gzproj(is,is)
    end do
    !
    gz_neq_local_angular_momenta(1) = trace_phi_basis(gzproj,phi_traces_basis_spin2(:,:))
    gz_neq_local_angular_momenta(2) = trace_phi_basis(gzproj,phi_traces_basis_spinZ(:,:))
    gz_neq_local_angular_momenta(3) = trace_phi_basis(gzproj,phi_traces_basis_isoSpin2(:,:))
    gz_neq_local_angular_momenta(4) = trace_phi_basis(gzproj,phi_traces_basis_isoSpinZ(:,:))
    !
    gz_neq_Rhop = hopping_renormalization_normal(gzproj,vdm_diag)
    !
    gz_neq_unitary_constr = 0.d0
    do iphi=1,Nphi
       gz_neq_unitary_constr = gz_neq_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
    end do
    !+-> LOCAL DENSITY MATRIX
    
    !+-> LOCAL DENSITY-DENSITY

    !+-> HOPPING RENORMALIZATION

    !+-> LOCAL ANGULAR MOMENTA

    !+-> ENERGY

    !+-> DENSITY-CONSTRAINTS

    !
  end subroutine gz_neq_measure



  

  

END MODULE GZ_DYNAMICS
