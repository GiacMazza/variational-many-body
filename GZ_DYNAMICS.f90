MODULE GZ_DYNAMICS  
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
  !
CONTAINS
  !
  include 'gz_eom.f90'
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
          !
          gz_neq_dens_constr_gzproj(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
          !
          !          gz_neq_local_sc_order(is,js) = trace_phi_basis(gzproj,phi_traces_basis_sc_order(is,js,:,:))
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
    !    gz_neq_Qhop = hopping_renormalization_anomalous(gzproj,vdm_diag)
    !
    gz_neq_unitary_constr = 0.d0
    do iphi=1,Nphi
       gz_neq_unitary_constr = gz_neq_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
    end do

    !+- SLATER
    Estar=0.d0
    gz_neq_dens_constr_slater=0.d0
    do ik=1,Lk
       do is=1,Ns
          do js=1,Ns
             gz_neq_dens_constr_slater(is,js) = gz_neq_dens_constr_slater(is,js) + slater(is,js,ik)*wtk(ik)
             !
             do iis=1,Ns
                do jjs=1,Ns
                   Estar = Estar + conjg(gz_neq_Rhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*gz_neq_Rhop(jjs,js)*slater(is,js,ik)*wtk(ik)
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
    complex(8),dimension(Ns,Ns)     :: Hk
    complex(8),dimension(Nphi)      :: gzproj
    real(8)                         :: Estar,Eloc,Egz
    real(8),dimension(Ns)           :: vdm_diag

    integer                         :: is,js,ik,it,iphi,iis,jjs

    it=t2it(time,tstep)
    !
    call dynamicalVector_2_wfMatrix_superc(psi_t,slater,gzproj)  
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
       do is=1,Ns
          do js=1,Ns
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



  

  

END MODULE GZ_DYNAMICS
