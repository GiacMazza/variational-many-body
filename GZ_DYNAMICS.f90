MODULE GZ_DYNAMICS  
  ! scifor
  USE SF_LINALG
  USE SF_OPTIMIZE
  USE RK_IDE
  ! GZ rooutines
  USE GZ_VARS_GLOBAL
  USE GZ_MATRIX_BASIS
  USE GZ_EFFECTIVE_HOPPINGS
  USE GZ_neqAUX_FUNX
  USE GZ_LOCAL_HAMILTONIAN
  implicit none
  private
  !
  public :: gz_equations_of_motion
  public :: gz_equations_of_motion_superc
  public :: gz_equations_of_motion_superc_lgr
  public :: step_dynamics_td_lagrange_superc
  public :: bcs_equations_of_motion
  !
  public :: setup_neq_hamiltonian
  !
  public :: setup_neq_dynamics
  public :: setup_neq_dynamics_superc
  !
  public :: gz_neq_measure
  public :: gz_neq_measure_superc
  !
  public :: get_neq_local_dens
  public :: get_neq_local_dens_dens
  public :: get_neq_energies  
  public :: get_neq_dens_constr_slater
  public :: get_neq_dens_constr_gzproj
  public :: get_neq_dens_constrA_slater
  public :: get_neq_dens_constrA_gzproj
  public :: get_neq_unitary_constr
  public :: get_neq_Rhop
  public :: get_neq_Qhop
  public :: get_neq_local_angular_momenta
  public :: get_neq_local_sc
  public :: get_neq_nqp
  !  
  complex(8),dimension(:,:),allocatable :: gz_neq_local_density_matrix !
  real(8),dimension(:,:),allocatable    :: gz_neq_local_dens_dens
  real(8),dimension(4)                  :: gz_neq_local_angular_momenta
  real(8),dimension(3)                  :: gz_neq_energies
  complex(8),dimension(:,:),allocatable :: gz_neq_dens_constr_slater
  complex(8),dimension(:,:),allocatable :: gz_neq_dens_constr_gzproj
  complex(8),dimension(:,:),allocatable :: gz_neq_dens_constrA_slater
  complex(8),dimension(:,:),allocatable :: gz_neq_dens_constrA_gzproj
  real(8)                               :: gz_neq_unitary_constr
  complex(8),dimension(:,:),allocatable :: gz_neq_Rhop         
  complex(8),dimension(:,:),allocatable :: gz_neq_Qhop         
  complex(8),dimension(:,:),allocatable :: gz_neq_local_sc_order 
  real(8),dimension(:,:),allocatable    :: gz_neq_nqp        
  !
  complex(8),dimension(:,:,:),allocatable :: neq_lgr
CONTAINS
  !
  include 'gz_eom.f90'
  include 'bcs_eom.f90'
  !
  subroutine setup_neq_dynamics
    allocate(gz_neq_local_density_matrix(Ns,Ns)); gz_neq_local_density_matrix = 0.d0
    allocate(gz_neq_local_dens_dens(Ns,Ns)); gz_neq_local_dens_dens = 0.d0
    gz_neq_local_angular_momenta = 0.d0
    gz_neq_energies = 0.d0
    allocate(gz_neq_dens_constr_slater(Ns,Ns)); gz_neq_dens_constr_slater = 0.d0
    allocate(gz_neq_dens_constr_gzproj(Ns,Ns)); gz_neq_dens_constr_gzproj = 0.d0
    gz_neq_unitary_constr = 0.d0
    allocate(gz_neq_Rhop(Ns,Ns)); gz_neq_Rhop = 0.d0
  end subroutine setup_neq_dynamics
  !
  subroutine setup_neq_dynamics_superc    
    allocate(gz_neq_local_density_matrix(Ns,Ns)); gz_neq_local_density_matrix = 0.d0
    allocate(gz_neq_local_dens_dens(Ns,Ns)); gz_neq_local_dens_dens = 0.d0
    gz_neq_local_angular_momenta = 0.d0
    gz_neq_energies = 0.d0
    allocate(gz_neq_dens_constr_slater(Ns,Ns)); gz_neq_dens_constr_slater = 0.d0
    allocate(gz_neq_dens_constr_gzproj(Ns,Ns)); gz_neq_dens_constr_gzproj = 0.d0
    !
    allocate(gz_neq_dens_constrA_slater(Ns,Ns)); gz_neq_dens_constrA_slater = 0.d0
    allocate(gz_neq_dens_constrA_gzproj(Ns,Ns)); gz_neq_dens_constrA_gzproj = 0.d0
    !
    gz_neq_unitary_constr = 0.d0
    allocate(gz_neq_Rhop(Ns,Ns)); gz_neq_Rhop = 0.d0
    allocate(gz_neq_Qhop(Ns,Ns)); gz_neq_Qhop = 0.d0
    allocate(gz_neq_local_sc_order(Ns,Ns)); gz_neq_local_sc_order = 0.d0
    allocate(gz_neq_nqp(Ns,Lk)); gz_neq_nqp = 0.d0
  end subroutine setup_neq_dynamics_superc


  subroutine get_neq_local_dens(is,js,x)
    integer :: is,js
    complex(8) :: x
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
  subroutine get_neq_local_sc(is,js,x)
    integer :: is,js
    complex(8) :: x
    x=gz_neq_local_sc_order(is,js)
  end subroutine get_neq_local_sc
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
  subroutine get_neq_dens_constrA_slater(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_neq_dens_constr_slater wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_neq_dens_constr_slater wrong indeces'
    x = gz_neq_dens_constrA_slater(is,js)
  end subroutine get_neq_dens_constrA_slater
  !
  subroutine get_neq_dens_constrA_gzproj(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_neq_dens_constr_slater wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_neq_dens_constr_slater wrong indeces'
    x = gz_neq_dens_constrA_gzproj(is,js)    
  end subroutine get_neq_dens_constrA_gzproj
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
  subroutine get_neq_Qhop(is,js,x)
    integer :: is,js
    complex(8) :: x
    x=gz_neq_Qhop(is,js)
  end subroutine get_neq_Qhop
  !
  subroutine get_neq_nqp(is,ik,x)
    integer :: is,ik
    real(8) :: x
    x = gz_neq_nqp(is,ik)
  end subroutine get_neq_nqp
  !
  subroutine setup_neq_hamiltonian(Uloc_t_,Ust_t_,Jh_t_,Jph_t_,Jsf_t_,eLevels_t_)
    real(8),dimension(3,Nt_aux),optional           :: Uloc_t_
    real(8),dimension(Nt_aux),optional             :: Ust_t_
    real(8),dimension(Nt_aux),optional             :: Jh_t_
    real(8),dimension(Nt_aux),optional             :: Jph_t_
    real(8),dimension(Nt_aux),optional             :: Jsf_t_
    real(8),dimension(Ns,Nt_aux),optional          :: eLevels_t_
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
    ! allocate(Hk_tb_t(Ns,Ns,Lk,Nt_aux)); forall(it=1:Nt_aux) Hk_tb_t(:,:,:,it) = Hk_tb
    ! if(present(Hk_tb_t_)) Hk_tb_t = Hk_tb_t_
    !
  end subroutine setup_neq_hamiltonian


  ! subroutine setup_neq_hamiltonian(Uloc_t_,Ust_t_,Jh_t_,Jph_t_,Jsf_t_,eLevels_t_,Hk_tb_t_)
  !   real(8),dimension(3,Nt_aux),optional           :: Uloc_t_
  !   real(8),dimension(Nt_aux),optional             :: Ust_t_
  !   real(8),dimension(Nt_aux),optional             :: Jh_t_
  !   real(8),dimension(Nt_aux),optional             :: Jph_t_
  !   real(8),dimension(Nt_aux),optional             :: Jsf_t_
  !   real(8),dimension(Ns,Nt_aux),optional          :: eLevels_t_
  !   complex(8),dimension(Ns,Ns,Lk,Nt_aux),optional :: Hk_tb_t_
  !   !
  !   integer                                    :: it
  !   !
  !   allocate(Uloc_t(3,Nt_aux)); forall(it=1:Nt_aux) Uloc_t(:,it) = Uloc
  !   if(present(Uloc_t_)) Uloc_t = Uloc_t_
  !   !
  !   allocate(Ust_t(Nt_aux)); Ust_t = Ust
  !   if(present(Ust_t_)) Ust_t = Ust_t_
  !   !
  !   allocate(Jh_t(Nt_aux)); Jh_t = Jh
  !   if(present(Jh_t_)) Jh_t = Jh_t_
  !   !
  !   allocate(Jph_t(Nt_aux)); Jph_t = Jph
  !   if(present(Jph_t_)) Jph_t = Jph_t_
  !   !
  !   allocate(Jsf_t(Nt_aux)); Jsf_t = Jsf
  !   if(present(Jsf_t_)) Jsf_t = Jsf_t_
  !   !
  !   allocate(eLevels_t(Ns,Nt_aux)); forall(it=1:Nt_aux) eLevels_t(:,it) = eLevels
  !   if(present(eLevels_t_)) eLevels_t = eLevels_t_
  !   !
  !   allocate(Hk_tb_t(Ns,Ns,Lk,Nt_aux)); forall(it=1:Nt_aux) Hk_tb_t(:,:,:,it) = Hk_tb
  !   if(present(Hk_tb_t_)) Hk_tb_t = Hk_tb_t_
  !   !
  ! end subroutine setup_neq_hamiltonian



  subroutine gz_neq_measure(psi_t,time)
    complex(8),dimension(nDynamics) :: psi_t
    real(8)                         :: time
    complex(8),dimension(Ns,Ns,Lk)  :: slater
    complex(8),dimension(Ns,Ns)     :: Hk,Hk_tmp
    complex(8),dimension(Nphi)      :: gzproj
    real(8)                         :: Estar,Eloc,Egz
    real(8),dimension(Ns)           :: vdm_diag

    integer                         :: is,js,ik,it,iphi,iis,jjs

    it=t2it(time,tstep)
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
          gz_neq_dens_constr_gzproj(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
          !
       end do
       vdm_diag(is) = gz_neq_dens_constr_gzproj(is,is)
    end do
    !
    Eloc = trace_phi_basis(gzproj,phi_traces_basis_Hloc)
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

    !+- SLATER -+!
    Estar=0.d0
    gz_neq_dens_constr_slater=0.d0
    do ik=1,Lk
       call get_Hk_t(Hk,ik,time)       
       do is=1,Ns
          do js=1,Ns
             gz_neq_dens_constr_slater(is,js) = gz_neq_dens_constr_slater(is,js) + slater(is,js,ik)*wtk(ik)
             !
             do iis=1,Ns
                do jjs=1,Ns
                   !Estar = Estar + conjg(gz_neq_Rhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*gz_neq_Rhop(jjs,js)*slater(is,js,ik)*wtk(ik)
                   Estar = Estar + conjg(gz_neq_Rhop(iis,is))*Hk(iis,jjs)*gz_neq_Rhop(jjs,js)*slater(is,js,ik)*wtk(ik)
                end do
             end do
             !
          end do
       end do
    end do
    !
    gz_neq_energies(1) = Estar+Eloc
    gz_neq_energies(2) = Estar
    gz_neq_energies(3) = Eloc
    !  
  end subroutine gz_neq_measure



  subroutine gz_neq_measure_superc(psi_t,time)
    complex(8),dimension(nDynamics) :: psi_t
    real(8)                         :: time
    complex(8),dimension(2,Ns,Ns,Lk)  :: slater
    complex(8),dimension(3,Ns,Ns,Lk)  :: slater_
    complex(8),dimension(Ns,Ns)     :: Hk,Hk_tmp
    complex(8),dimension(Ns,Ns)     :: Rhop,Qhop,Rhop_dag,Qhop_dag
    complex(8),dimension(2*Ns,2*Ns)     :: Hks
    real(8),dimension(2*Ns)         :: eks
    complex(8),dimension(Nphi)      :: gzproj
    real(8)                         :: Estar,Eloc,Egz
    real(8),dimension(Ns)           :: vdm_diag
    real(8)                         :: nqp

    integer                         :: is,js,ik,it,iphi,iis,jjs

    it=t2it(time,tstep)
    !
    call dynamicalVector_2_wfMatrix_superc(psi_t,slater,gzproj)  
    slater_(1:2,:,:,:) = slater
    slater_(3,:,:,:) = zero
    do is=1,Ns
       slater_(3,is,is,:) = one 
       do js=1,Ns
          slater_(3,is,js,:) = slater_(3,is,js,:) - slater(1,js,is,:)
       end do
    end do
    !
    Estar=0.d0
    do is=1,Ns
       do js=1,Ns
          gz_neq_local_density_matrix(is,js) = &
               trace_phi_basis(gzproj,phi_traces_basis_local_dens(is,js,:,:))          
          gz_neq_local_dens_dens(is,js) = &
               trace_phi_basis(gzproj,phi_traces_basis_dens_dens(is,js,:,:))
          !
          !
          gz_neq_dens_constr_gzproj(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
          gz_neq_dens_constrA_gzproj(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens_anomalous(is,js,:,:))
          !
          gz_neq_local_sc_order(is,js) = trace_phi_basis(gzproj,phi_traces_basis_sc_order(is,js,:,:))
       end do
       vdm_diag(is) = gz_neq_dens_constr_gzproj(is,is)
    end do
    !
    Eloc = trace_phi_basis(gzproj,phi_traces_basis_Hloc)
    !
    gz_neq_local_angular_momenta(1) = trace_phi_basis(gzproj,phi_traces_basis_spin2(:,:))
    gz_neq_local_angular_momenta(2) = trace_phi_basis(gzproj,phi_traces_basis_spinZ(:,:))
    gz_neq_local_angular_momenta(3) = trace_phi_basis(gzproj,phi_traces_basis_isoSpin2(:,:))
    gz_neq_local_angular_momenta(4) = trace_phi_basis(gzproj,phi_traces_basis_isoSpinZ(:,:))
    !
    gz_neq_Rhop = hopping_renormalization_normal(gzproj,vdm_diag)
    gz_neq_Qhop = hopping_renormalization_anomalous(gzproj,vdm_diag)
    !
    Rhop=gz_neq_Rhop
    Qhop=gz_neq_Qhop
    !
    do is=1,Ns
       do js=1,Ns
          Rhop_dag(is,js) = conjg(Rhop(js,is))
          Qhop_dag(is,js) = conjg(Qhop(js,is))
       end do
    end do

    !
    gz_neq_unitary_constr = 0.d0
    do iphi=1,Nphi
       gz_neq_unitary_constr = gz_neq_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
    end do

    !+- SLATER
    Estar=0.d0
    gz_neq_dens_constr_slater=0.d0
    gz_neq_dens_constrA_slater=0.d0
    do ik=1,Lk
       call get_Hk_t(Hk,ik,time)
       !+- define Hk_renormalized -+!
       Hk_tmp=matmul(Hk,Rhop)
       Hk_tmp=matmul(Rhop_dag,Hk_tmp)
       Hks(1:Ns,1:Ns) = Hk_tmp 
       !
       Hk_tmp=matmul(Hk,Qhop)
       Hk_tmp=matmul(Rhop_dag,Hk_tmp)
       Hks(1:Ns,Ns+1:2*Ns) = Hk_tmp 
       !
       Hk_tmp=matmul(Hk,Rhop)
       Hk_tmp=matmul(Qhop_dag,Hk_tmp)
       Hks(Ns+1:2*Ns,1:Ns) = Hk_tmp
       !
       Hk_tmp=matmul(Hk,Qhop)
       Hk_tmp=matmul(Qhop_dag,Hk_tmp)
       Hks(Ns+1:2*Ns,Ns+1:2*Ns) = Hk_tmp
       !
       !+- get eigenstates
       call matrix_diagonalize(Hks,eks)
       !
       !+- project onto the natural basis
       ! do is=1,Ns
       !    gz_neq_nqp(is,ik) = 0.d0
       !    do iis=1,Ns
       !       do jjs=1,Ns
       !          gz_neq_nqp(is,ik) = gz_neq_nqp(is,ik) + &
       !               Hks(iis,is)*conjg(Hks(jjs,is))*slater_(1,iis,jjs,ik)
       !          gz_neq_nqp(is,ik) = gz_neq_nqp(is,ik) + &
       !               Hks(iis,is)*conjg(Hks(jjs+Ns,is))*slater_(2,iis,jjs,ik)
       !          gz_neq_nqp(is,ik) = gz_neq_nqp(is,ik) + &
       !               Hks(iis+Ns,is)*conjg(Hks(jjs,is))*conjg(slater_(2,jjs,iis,ik))
       !          gz_neq_nqp(is,ik) = gz_neq_nqp(is,ik) + &
       !               Hks(iis+Ns,is)*conjg(Hks(jjs+Ns,is))*slater_(3,iis,jjs,ik)
                
       !       enddo
       !    enddo
       ! enddo



       do is=1,Ns
          do js=1,Ns
             !
             gz_neq_dens_constr_slater(is,js) = gz_neq_dens_constr_slater(is,js) + slater(1,is,js,ik)*wtk(ik)
             gz_neq_dens_constrA_slater(is,js) = gz_neq_dens_constrA_slater(is,js) + slater(2,is,js,ik)*wtk(ik)
             !
             do iis=1,Ns
                do jjs=1,Ns
                   !
                   Estar = Estar + conjg(gz_neq_Rhop(iis,is))*Hk(iis,jjs)*gz_neq_Rhop(jjs,js)*slater(1,is,js,ik)*wtk(ik)
                   Estar = Estar + conjg(gz_neq_Rhop(iis,is))*Hk(iis,jjs)*gz_neq_Qhop(jjs,js)*slater(2,is,js,ik)*wtk(ik)
                   Estar = Estar + conjg(gz_neq_Qhop(iis,is))*Hk(iis,jjs)*gz_neq_Rhop(jjs,js)*conjg(slater(2,js,is,ik))*wtk(ik)
                   Estar = Estar - conjg(gz_neq_Qhop(iis,is))*Hk(iis,jjs)*gz_neq_Qhop(jjs,js)*slater(1,js,is,ik)*wtk(ik)                   
                   if(is.eq.js) then
                      Estar = Estar + conjg(gz_neq_Qhop(iis,is))*Hk(iis,jjs)*gz_neq_Qhop(jjs,js)*wtk(ik)
                   end if
                   ! Estar = Estar + conjg(gz_neq_Rhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*gz_neq_Rhop(jjs,js)*slater(1,is,js,ik)*wtk(ik)
                   ! Estar = Estar + conjg(gz_neq_Rhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*gz_neq_Qhop(jjs,js)*slater(2,is,js,ik)*wtk(ik)
                   ! Estar = Estar + conjg(gz_neq_Qhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*gz_neq_Rhop(jjs,js)*conjg(slater(2,js,is,ik))*wtk(ik)
                   ! Estar = Estar - conjg(gz_neq_Qhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*gz_neq_Qhop(jjs,js)*slater(1,js,is,ik)*wtk(ik)                   
                   ! if(is.eq.js) then
                   !    Estar = Estar + conjg(gz_neq_Qhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*gz_neq_Qhop(jjs,js)*wtk(ik)
                   ! end if
                   !
                end do
             end do
             !
          end do
       end do
    end do
    !
    gz_neq_energies(1) = Estar+Eloc
    gz_neq_energies(2) = Estar
    gz_neq_energies(3) = Eloc
    !
  end subroutine gz_neq_measure_superc




  subroutine step_dynamics_td_lagrange_superc(nsys,tstep,t,yt,td_lgr,eom_funct) 
    integer :: nsys
    real(8) :: tstep,t
    complex(8),dimension(nsys),intent(inout) :: yt
    complex(8),dimension(nsys) :: yt_old,yt_new
    complex(8),dimension(2,Ns,Ns),intent(inout) :: td_lgr
    real(8),dimension(:),allocatable            :: lgr,delta_out
    complex(8),dimension(:),allocatable :: lgr_cmplx
    integer :: iter,Nopt
    integer :: i,i0
    real(8) :: delta
    interface
       function eom_funct(t,y,Nsys)
         implicit none
         integer                    :: Nsys
         real(8)                    :: t   
         complex(8),dimension(Nsys) :: eom_funct
         complex(8),dimension(Nsys) :: y
       end function eom_funct
    end interface

    Nopt=2*Nvdm_AC_opt
    allocate(lgr_cmplx(Nvdm_AC_opt))
    allocate(lgr(2*Nopt));allocate(delta_out(2*Nopt))
    !
    call vdm_AC_stride_m2v(td_lgr(1,:,:),lgr_cmplx)    
    do i=1,Nvdm_AC_opt
       lgr(i) = dreal(lgr_cmplx(i))
       lgr(i+Nvdm_AC_opt) = dimag(lgr_cmplx(i))
    end do
    call vdm_AC_stride_m2v(td_lgr(2,:,:),lgr_cmplx)    
    i0 = 2*Nvdm_AC_opt
    do i=1,Nvdm_AC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
       lgr(i0+i+Nvdm_AC_opt) = dimag(lgr_cmplx(i))
    end do
    !
    yt_old = yt
    !
    delta=fix_anomalous_vdm_(lgr)
    !
    if(delta.gt.1.d-10) then
       write(*,*) 'td-fixing of lagrange parameters   ',delta
       !call fsolve(fix_anomalous_vdm,lgr,tol=1.d-06,info=iter)    
       call fmin_cg(lgr,fix_anomalous_vdm_,iter,delta,itmax=20)
       delta_out = fix_anomalous_vdm(lgr)
       delta=0.d0
       do i=1,2*Nopt
          delta = delta + delta_out(i)**2.d0
       end do
       write(*,*) 'Time Dependent lagrange parameters'
       write(*,*) lgr
       write(*,*) 'Time Dependent lagrange parameters: error'
       write(*,*) delta
       write(*,*)
    end if
    !
    yt=yt_new
    !
    lgr_cmplx=zero
    do i=1,Nvdm_AC_opt
       lgr_cmplx(i) = lgr(i)+xi*lgr(i+Nvdm_AC_opt)
    end do
    call vdm_AC_stride_v2m(lgr_cmplx,td_lgr(1,:,:))
    i0=2*Nvdm_AC_opt
    do i=1,Nvdm_AC_opt
       lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_AC_opt)
    end do
    call vdm_AC_stride_v2m(lgr_cmplx,td_lgr(2,:,:))  

    !    
  contains

    function fix_anomalous_vdm(lgr) result(delta)
      real(8),dimension(:) :: lgr
      real(8),dimension(size(lgr)) :: delta
      complex(8),dimension(:),allocatable  :: lgr_cmplx,delta_cmplx
      complex(8),dimension(Ns,Ns) :: anomalous_constrGZ,anomalous_constrSL
      integer :: i0,i,is,js
      real(8) :: tmp_test
      !
      if(allocated(neq_lgr)) deallocate(neq_lgr)
      allocate(neq_lgr(2,Ns,Ns)); neq_lgr=zero
      !
      !+- dump slater_lgr_multipliers
      allocate(lgr_cmplx(Nvdm_AC_opt))
      do i=1,Nvdm_AC_opt
         lgr_cmplx(i) = lgr(i) + xi*lgr(i+Nvdm_AC_opt)
      end do
      call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(1,:,:))
      i0 = 2*Nvdm_AC_opt
      !+- dump gzproj_lgr_multipliers
      do i=1,Nvdm_AC_opt
         lgr_cmplx(i) = lgr(i+i0) + xi*lgr(i+i0+Nvdm_AC_opt)
      end do
      call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(2,:,:))
      !
      yt_new = RK_step(nDynamics,4,tstep,t,yt_old,gz_equations_of_motion_superc_lgr)
      !
      call gz_neq_measure_superc(yt_new,t)
      !
      do is=1,Ns
         do js=1,Ns
            call get_neq_dens_constrA_slater(is,js,anomalous_constrSL(is,js))
            call get_neq_dens_constrA_gzproj(is,js,anomalous_constrGZ(is,js))
         end do
      end do
      !
      delta=0.d0
      allocate(delta_cmplx(Nvdm_AC_opt))
      call vdm_AC_stride_m2v(anomalous_constrSL,delta_cmplx)
      do i=1,Nvdm_AC_opt
         delta(i) = dreal(delta_cmplx(i))
         delta(i+Nvdm_AC_opt) = dimag(delta_cmplx(i))
      end do
      deallocate(delta_cmplx)
      i0 = 2*Nvdm_AC_opt
      allocate(delta_cmplx(Nvdm_AC_opt))
      call vdm_AC_stride_m2v(anomalous_constrGZ,delta_cmplx)
      do i=1,Nvdm_AC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
         delta(i0+i+Nvdm_AC_opt) = dimag(delta_cmplx(i))       
      end do
      deallocate(delta_cmplx)
      tmp_test=0.d0
      do i=1,size(lgr)
         tmp_test=tmp_test+delta(i)**2.d0
      end do
      write(*,*) delta
      !
    end function fix_anomalous_vdm













    function fix_anomalous_vdm_(lgr) result(delta)
      real(8),dimension(:) :: lgr
      real(8) :: delta
      complex(8),dimension(:),allocatable  :: lgr_cmplx,delta_cmplx
      complex(8),dimension(Ns,Ns) :: anomalous_constrGZ,anomalous_constrSL
      integer :: i0,i,is,js
      !
      if(allocated(neq_lgr)) deallocate(neq_lgr)
      allocate(neq_lgr(2,Ns,Ns)); neq_lgr=zero

      !+- dump slater_lgr_multipliers
      allocate(lgr_cmplx(Nvdm_AC_opt))
      do i=1,Nvdm_AC_opt
         lgr_cmplx(i) = lgr(i) + xi*lgr(i+Nvdm_AC_opt)
      end do
      call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(1,:,:))
      i0 = 2*Nvdm_AC_opt
      !+- dump gzproj_lgr_multipliers
      do i=1,Nvdm_AC_opt
         lgr_cmplx(i) = lgr(i+i0) + xi*lgr(i+i0+Nvdm_AC_opt)
      end do
      call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(2,:,:))
      !
      yt_new = RK_step(nDynamics,4,tstep,t,yt_old,gz_equations_of_motion_superc_lgr)
      !
      call gz_neq_measure_superc(yt_new,t)
      !
      do is=1,Ns
         do js=1,Ns
            call get_neq_dens_constrA_slater(is,js,anomalous_constrSL(is,js))
            call get_neq_dens_constrA_gzproj(is,js,anomalous_constrGZ(is,js))
         end do
      end do
      !
      delta=0.d0
      do is=1,Ns
         do js=1,Ns
            delta = delta + anomalous_constrSL(is,js)*conjg(anomalous_constrSL(is,js))
            delta = delta + anomalous_constrGZ(is,js)*conjg(anomalous_constrGZ(is,js))
         end do
      end do

      write(*,*) 'deviation from constraint conservation',delta
      !
      ! allocate(delta_cmplx(Nvdm_AC_opt))
      ! call vdm_AC_stride_m2v(anomalous_constrSL,delta_cmplx)
      ! do i=1,Nvdm_AC_opt
      !    delta(i) = dreal(delta_cmplx(i))
      !    delta(i+Nvdm_AC_opt) = dimag(delta_cmplx(i))
      ! end do
      ! deallocate(delta_cmplx)
      ! i0 = 2*Nvdm_AC_opt
      ! allocate(delta_cmplx(Nvdm_AC_opt))
      ! call vdm_AC_stride_m2v(anomalous_constrGZ,delta_cmplx)
      ! do i=1,Nvdm_AC_opt
      !    delta(i0+i) = dreal(delta_cmplx(i))
      !    delta(i0+i+Nvdm_AC_opt) = dimag(delta_cmplx(i))       
      ! end do
      ! deallocate(delta_cmplx)
      !
    end function fix_anomalous_vdm_





  end subroutine step_dynamics_td_lagrange_superc
  

  

END MODULE GZ_DYNAMICS

