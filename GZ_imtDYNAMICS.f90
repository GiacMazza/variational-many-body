MODULE GZ_imtDYNAMICS  
  ! scifor
  USE SF_LINALG
  USE SF_OPTIMIZE
  USE SF_SPECIAL
  USE RK_IDE
  ! GZ rooutines
  USE GZ_VARS_GLOBAL
  USE GZ_MATRIX_BASIS
  USE GZ_EFFECTIVE_HOPPINGS
  USE GZ_neqAUX_FUNX
  USE GZ_LOCAL_HAMILTONIAN
  USE MATRIX_SPARSE
  implicit none
  private
  !

  !
  public :: step_imt_dynamics
  public :: gz_imt_eom
  !
  public :: setup_imt_hamiltonian
  public :: init_imt_qpH
  !
  public :: setup_imt_dynamics
  public :: setup_imt_dynamics_superc
  !
  public :: gz_imt_measure!,gz_imt_measure_sp
  !public :: gz_imt_measure_superc,gz_imt_measure_superc_sp
  !
  public :: get_imt_local_dens
  public :: get_imt_local_dens_dens
  public :: get_imt_energies  
  public :: get_imt_dens_constr_slater
  public :: get_imt_dens_constr_gzproj
  public :: get_imt_dens_constrA_slater
  public :: get_imt_dens_constrA_gzproj
  public :: get_imt_lgrA_slater
  public :: get_imt_lgrA_gzproj
  public :: get_imt_unitary_constr
  public :: get_imt_Rhop
  public :: get_imt_Qhop
  public :: get_imt_local_angular_momenta
  public :: get_imt_local_sc
  public :: get_imt_nqp
  !  
  complex(8),dimension(:,:),allocatable :: gz_imt_local_density_matrix !
  real(8),dimension(:,:),allocatable    :: gz_imt_local_dens_dens
  real(8),dimension(4)                  :: gz_imt_local_angular_momenta
  real(8),dimension(3)                  :: gz_imt_energies
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_constr_slater
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_constr_gzproj
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_constrA_slater
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_constrA_gzproj
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_lgrA_slater
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_lgrA_gzproj
  real(8)                               :: gz_imt_unitary_constr
  complex(8),dimension(:,:),allocatable :: gz_imt_Rhop         
  complex(8),dimension(:,:),allocatable :: gz_imt_Qhop         
  complex(8),dimension(:,:),allocatable :: gz_imt_local_sc_order 
  real(8),dimension(:,:),allocatable    :: gz_imt_nqp        
  !

  complex(8),dimension(:,:,:),allocatable :: gz_imt_qpH,Hksave

  complex(8),dimension(:,:),allocatable,public :: imt_lgrNC
  real(8) :: imt_lgrU
  !

  real(8) :: imt_tstep
CONTAINS
  !
  include 'gz_IMTeom.f90'
  !
  subroutine setup_imt_dynamics
    allocate(gz_imt_local_density_matrix(Ns,Ns)); gz_imt_local_density_matrix = 0.d0
    allocate(gz_imt_local_dens_dens(Ns,Ns)); gz_imt_local_dens_dens = 0.d0
    gz_imt_local_angular_momenta = 0.d0
    gz_imt_energies = 0.d0
    allocate(gz_imt_dens_constr_slater(Ns,Ns)); gz_imt_dens_constr_slater = 0.d0
    allocate(gz_imt_dens_constr_gzproj(Ns,Ns)); gz_imt_dens_constr_gzproj = 0.d0
    gz_imt_unitary_constr = 0.d0
    allocate(gz_imt_Rhop(Ns,Ns)); gz_imt_Rhop = 0.d0
  end subroutine setup_imt_dynamics

  subroutine init_imt_qpH(beta0,gzproj,lgrNC_)
    real(8) :: beta0
    complex(8),dimension(Nphi)  :: gzproj
    complex(8),dimension(Ns,Ns),optional :: lgrNC_
    complex(8),dimension(Ns,Ns) :: lgrNC,Hk,Rhop,Rhop_hc
    real(8),dimension(Ns) :: vdm_diag
    integer :: is,js,ik
    !
    lgrNC = zero
    if(present(lgrNC_)) lgrNC = lgrNC_
    allocate(gz_imt_qpH(Lk,Ns,Ns)); gz_imt_qpH=zero
    allocate(Hksave(Lk,Ns,Ns)); HKsave=zero
    !
    do is=1,Ns
       vdm_diag(is) = dreal(trace_phi_basis(gzproj,phi_traces_basis_dens(is,is,:,:)))
    end do
    !
    Rhop = hopping_renormalization_normal_sp(gzproj,vdm_diag)
    Rhop_hc = conjg(transpose(Rhop))
    !
    do ik=1,Lk
       Hk=matmul(Hk_tb(:,:,ik),Rhop)
       Hk=matmul(Rhop_hc,Hk)
       Hk=Hk+lgrNC       
       !
       gz_imt_qpH(ik,:,:) = 2.d0*beta0*Hk  !+- il due, capra!!!
       !
    end do
  end subroutine init_imt_qpH
  !
  subroutine setup_imt_dynamics_superc    
    allocate(gz_imt_local_density_matrix(Ns,Ns)); gz_imt_local_density_matrix = 0.d0
    allocate(gz_imt_local_dens_dens(Ns,Ns)); gz_imt_local_dens_dens = 0.d0
    gz_imt_local_angular_momenta = 0.d0
    gz_imt_energies = 0.d0
    allocate(gz_imt_dens_constr_slater(Ns,Ns)); gz_imt_dens_constr_slater = 0.d0
    allocate(gz_imt_dens_constr_gzproj(Ns,Ns)); gz_imt_dens_constr_gzproj = 0.d0
    !
    allocate(gz_imt_dens_constrA_slater(Ns,Ns)); gz_imt_dens_constrA_slater = 0.d0
    allocate(gz_imt_dens_constrA_gzproj(Ns,Ns)); gz_imt_dens_constrA_gzproj = 0.d0
    allocate(gz_imt_dens_lgrA_slater(Ns,Ns)); gz_imt_dens_lgrA_slater = 0.d0
    allocate(gz_imt_dens_lgrA_gzproj(Ns,Ns)); gz_imt_dens_lgrA_gzproj = 0.d0
    !
    gz_imt_unitary_constr = 0.d0
    allocate(gz_imt_Rhop(Ns,Ns)); gz_imt_Rhop = 0.d0
    allocate(gz_imt_Qhop(Ns,Ns)); gz_imt_Qhop = 0.d0
    allocate(gz_imt_local_sc_order(Ns,Ns)); gz_imt_local_sc_order = 0.d0
    allocate(gz_imt_nqp(Ns,Lk)); gz_imt_nqp = 0.d0
  end subroutine setup_imt_dynamics_superc





  subroutine get_imt_local_dens(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_local_dens wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_local_dens wrong indeces'
    x = gz_imt_local_density_matrix(is,js)    
  end subroutine get_imt_local_dens
  !
  subroutine get_imt_local_dens_dens(is,js,x)
    integer :: is,js
    real(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_local_dens_dens wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_local_dens_dens wrong indeces'
    x = gz_imt_local_dens_dens(is,js)    
  end subroutine get_imt_local_dens_dens
  !
  subroutine get_imt_local_angular_momenta(x)
    real(8),dimension(4) :: x
    x = gz_imt_local_angular_momenta
  end subroutine get_imt_local_angular_momenta
  !
  subroutine get_imt_local_sc(is,js,x)
    integer :: is,js
    complex(8) :: x
    x=gz_imt_local_sc_order(is,js)
  end subroutine get_imt_local_sc
  !
  subroutine get_imt_energies(x)
    real(8),dimension(3) :: x
    x = gz_imt_energies
  end subroutine get_imt_energies
  !
  subroutine get_imt_dens_constr_slater(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_constr_slater wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_constr_slater wrong indeces'
    x = gz_imt_dens_constr_slater(is,js)
  end subroutine get_imt_dens_constr_slater
  !
  subroutine get_imt_dens_constr_gzproj(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_constr_slater wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_constr_slater wrong indeces'
    x = gz_imt_dens_constr_gzproj(is,js)    
  end subroutine get_imt_dens_constr_gzproj
  !
  subroutine get_imt_dens_constrA_slater(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_constr_slater wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_constr_slater wrong indeces'
    x = gz_imt_dens_constrA_slater(is,js)
  end subroutine get_imt_dens_constrA_slater
  !
  subroutine get_imt_dens_constrA_gzproj(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_constr_gzproj wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_constr_gzproj wrong indeces'
    x = gz_imt_dens_constrA_gzproj(is,js)    
  end subroutine get_imt_dens_constrA_gzproj
  !

  subroutine get_imt_lgrA_slater(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_lgrA_sl wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_lgrA_sl wrong indeces'
    x = gz_imt_dens_lgrA_slater(is,js)
  end subroutine get_imt_lgrA_slater
  !
  subroutine get_imt_lgrA_gzproj(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_lgrA_GZ wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_lgrA_GZ wrong indeces'
    x = gz_imt_dens_lgrA_gzproj(is,js)    
  end subroutine get_imt_lgrA_gzproj
  !

  subroutine get_imt_unitary_constr(x)
    real(8) :: x
    x =gz_imt_unitary_constr
  end subroutine get_imt_unitary_constr
  !
  subroutine get_imt_Rhop(is,js,x)
    implicit none
    integer :: is,js
    complex(8),intent(inout) :: x
    x=gz_imt_Rhop(is,js)
  end subroutine get_imt_Rhop
  !
  subroutine get_imt_Qhop(is,js,x)
    integer :: is,js
    complex(8) :: x
    x=gz_imt_Qhop(is,js)
  end subroutine get_imt_Qhop
  !
  subroutine get_imt_nqp(is,ik,x)
    integer :: is,ik
    real(8) :: x
    x = gz_imt_nqp(is,ik)
  end subroutine get_imt_nqp
  !
  subroutine setup_imt_hamiltonian(Uloc_t_,Ust_t_,Jh_t_,Jph_t_,Jsf_t_,eLevels_t_)
    real(8),dimension(3,Nit_aux),optional           :: Uloc_t_
    real(8),dimension(Nit_aux),optional             :: Ust_t_
    real(8),dimension(Nit_aux),optional             :: Jh_t_
    real(8),dimension(Nit_aux),optional             :: Jph_t_
    real(8),dimension(Nit_aux),optional             :: Jsf_t_
    real(8),dimension(Ns,Nit_aux),optional          :: eLevels_t_
    !
    integer                                    :: it
    !
    allocate(Uloc_t(3,Nit_aux)); forall(it=1:Nit_aux) Uloc_t(:,it) = Uloc
    if(present(Uloc_t_)) Uloc_t = Uloc_t_
    !
    allocate(Ust_t(Nit_aux)); Ust_t = Ust
    if(present(Ust_t_)) Ust_t = Ust_t_
    !
    allocate(Jh_t(Nit_aux)); Jh_t = Jh
    if(present(Jh_t_)) Jh_t = Jh_t_
    !
    allocate(Jph_t(Nit_aux)); Jph_t = Jph
    if(present(Jph_t_)) Jph_t = Jph_t_
    !
    allocate(Jsf_t(Nit_aux)); Jsf_t = Jsf
    if(present(Jsf_t_)) Jsf_t = Jsf_t_
    !
    allocate(eLevels_t(Ns,Nit_aux)); forall(it=1:Nit_aux) eLevels_t(:,it) = eLevels
    if(present(eLevels_t_)) eLevels_t = eLevels_t_
    !
  end subroutine setup_imt_hamiltonian


  subroutine gz_imt_measure(psi_t,itime)
    complex(8),dimension(nDynamics) :: psi_t
    real(8)                         :: itime
    complex(8),dimension(Ns,Ns,Lk)  :: slater
    complex(8),dimension(Ns,Ns)     :: tmpHk,Hk
    complex(8),dimension(Nphi)      :: gzproj,gztmp
    real(8)                         :: Estar,Eloc,Egz
    real(8),dimension(Ns)           :: vdm_diag,tmp_eHk

    integer                         :: is,js,ik,iphi,iis,jjs,itt,ks
    !
    if(Nphi.ne.nDynamics) stop "gz_imt_measure/ wrong dimensions"
    !call dynamicalVector_2_wfMatrix(psi_t,slater,gzproj)  
    gzproj=psi_t
    !
    do is=1,Ns
       do js=1,Ns
          gz_imt_local_density_matrix(is,js) = &
               trace_phi_basis_sp(gzproj,phi_spTraces_basis_local_dens(is,js))          
          gz_imt_local_dens_dens(is,js) = &
               trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens_dens(is,js))
          !
          gz_imt_dens_constr_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens(is,js))
          !
       end do
       vdm_diag(is) = gz_imt_dens_constr_gzproj(is,is)
    end do
    !
    itt=imt_t2it(itime,0.5d0*itstep,beta_init)
    ! Uloc=Uloc_t(:,itt)
    ! Ust =Ust_t(itt)
    ! Jh=Jh_t(itt)
    ! Jsf=Jsf_t(itt)
    ! Jph=Jph_t(itt)
    ! eLevels = eLevels_t(:,itt)
    eLevels=0.d0
    !
    gztmp = get_local_hamiltonian_HLOCphi(gzproj,eLevels)
    Eloc = 0.d0
    do iphi=1,Nphi
       Eloc = Eloc + conjg(gzproj(iphi))*gztmp(iphi)
    end do
    !
    gz_imt_local_angular_momenta(1) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_spin2)
    gz_imt_local_angular_momenta(2) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_spinZ)
    gz_imt_local_angular_momenta(3) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_isoSpin2)
    gz_imt_local_angular_momenta(4) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_isoSpinZ)
    !
    gz_imt_Rhop = hopping_renormalization_normal_sp(gzproj,vdm_diag)
    !
    gz_imt_unitary_constr = 0.d0
    do iphi=1,Nphi
       gz_imt_unitary_constr = gz_imt_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
    end do
    !+- SLATER -+!
    Estar=0.d0
    !+- first thing first I have to exctract slater(is,js,ik) from the immaginary time integral of the QP hamiltonian
    do ik=1,Lk
       !

       Hksave(ik,:,:) = matmul(Hk_tb(:,:,ik),gz_imt_rhop)
       Hksave(ik,:,:) = matmul(conjg(transpose(gz_imt_rhop)),Hksave(ik,:,:))

       tmpHk= gz_imt_qpH(ik,:,:) 
       call matrix_diagonalize(tmpHk,tmp_eHk)       
       !
       do is=1,Ns
          do js=1,Ns
             slater(is,js,ik) = zero
             do ks=1,Ns
                slater(is,js,ik) = slater(is,js,ik) + conjg(tmpHk(is,ks))*tmpHk(js,ks)*fermi(tmp_eHk(ks),1.d0)
             end do
          end do
       end do
       !

       if(ik.eq.10) write(343,'(18F18.10)') itime,slater(1,1,ik)
    end do
    !
    gz_imt_dens_constr_slater=0.d0
    do ik=1,Lk
       Hk = Hk_tb(:,:,ik)
       do is=1,Ns
          do js=1,Ns
             gz_imt_dens_constr_slater(is,js) = gz_imt_dens_constr_slater(is,js) + slater(is,js,ik)*wtk(ik)
             do iis=1,Ns
                do jjs=1,Ns
                   Estar = Estar + conjg(gz_imt_Rhop(iis,is))*Hk(iis,jjs)*gz_imt_Rhop(jjs,js)*slater(is,js,ik)*wtk(ik)
                end do
             end do
          end do
       end do
    end do
    !
    gz_imt_energies(1) = Estar+Eloc
    gz_imt_energies(2) = Estar
    gz_imt_energies(3) = Eloc
    !  
  end subroutine gz_imt_measure
  !
  subroutine gz_imt_measure_constr(psi_t,itime)
    complex(8),dimension(nDynamics) :: psi_t
    real(8)                         :: itime
    complex(8),dimension(Ns,Ns,Lk)  :: slater
    complex(8),dimension(Ns,Ns)     :: tmpHk,Hk
    complex(8),dimension(Nphi)      :: gzproj
    real(8)                         :: Estar,Eloc,Egz,b0
    real(8),dimension(Ns)           :: vdm_diag,tmp_eHk
    integer                         :: is,js,ik,it,iphi,iis,jjs,itt,ks
    !
    it=imt_t2it(itime,itstep)
    !
    if(Nphi.ne.nDynamics) stop "gz_imt_measure_constr/ wrong dimensions"
    gzproj=psi_t
    !
    do is=1,Ns
       do js=1,Ns
          gz_imt_dens_constr_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens(is,js))
       end do
    end do
    !
    gz_imt_unitary_constr = 0.d0
    do iphi=1,Nphi
       gz_imt_unitary_constr = gz_imt_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
    end do
    !
    Estar=0.d0
    !+- first thing first I have to exctract slater(is,js,ik) from the immaginary time integral of the QP hamiltonian
    do ik=1,Lk
       !
       tmpHk= gz_imt_qpH(ik,:,:) 
       call matrix_diagonalize(tmpHk,tmp_eHk)       
       !
       do is=1,Ns
          do js=1,Ns
             slater(is,js,ik) = zero
             do ks=1,Ns
                slater(is,js,ik) = slater(is,js,ik) + conjg(tmpHk(is,ks))*tmpHk(js,ks)*fermi(tmp_eHk(ks),1.d0)
             end do
          end do
       end do
       !       
    end do
    !
    gz_imt_dens_constr_slater=0.d0
    do ik=1,Lk
       Hk = Hk_tb(:,:,ik)
       do is=1,Ns
          do js=1,Ns
             gz_imt_dens_constr_slater(is,js) = gz_imt_dens_constr_slater(is,js) + slater(is,js,ik)*wtk(ik)
          end do
       end do
    end do
    ! !
  end subroutine gz_imt_measure_constr












  ! subroutine gz_imt_measure_superc(psi_t,time)
  !   complex(8),dimension(nDynamics) :: psi_t
  !   real(8)                         :: time
  !   complex(8),dimension(2,Ns,Ns,Lk)  :: slater
  !   complex(8),dimension(3,Ns,Ns,Lk)  :: slater_
  !   complex(8),dimension(Ns,Ns)     :: Hk,Hk_tmp
  !   complex(8),dimension(Ns,Ns)     :: Rhop,Qhop,Rhop_dag,Qhop_dag
  !   complex(8),dimension(2*Ns,2*Ns)     :: Hks
  !   real(8),dimension(2*Ns)         :: eks
  !   complex(8),dimension(Nphi)      :: gzproj  !+--> probably there some conflicts with this name variable here when also lgr multipliers are there +-!
  !   real(8)                         :: Estar,Eloc,Egz
  !   real(8),dimension(Ns)           :: vdm_diag
  !   real(8)                         :: nqp

  !   integer                         :: is,js,ik,it,iphi,iis,jjs,itt

  !   it=imt_t2it(time,tstep)
  !   !
  !   call dynamicalVector_2_wfMatrix_superc(psi_t,slater,gzproj)  
  !   slater_(1:2,:,:,:) = slater
  !   slater_(3,:,:,:) = zero
  !   do is=1,Ns
  !      slater_(3,is,is,:) = one 
  !      do js=1,Ns
  !         slater_(3,is,js,:) = slater_(3,is,js,:) - slater(1,js,is,:)
  !      end do
  !   end do
  !   !
  !   Estar=0.d0
  !   do is=1,Ns
  !      do js=1,Ns
  !         gz_imt_local_density_matrix(is,js) = &
  !              trace_phi_basis(gzproj,phi_traces_basis_local_dens(is,js,:,:))          
  !         gz_imt_local_dens_dens(is,js) = &
  !              trace_phi_basis(gzproj,phi_traces_basis_dens_dens(is,js,:,:))
  !         !
  !         !
  !         gz_imt_dens_constr_gzproj(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
  !         gz_imt_dens_constrA_gzproj(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens_anomalous(is,js,:,:))
  !         !
  !         gz_imt_local_sc_order(is,js) = trace_phi_basis(gzproj,phi_traces_basis_sc_order(is,js,:,:))
  !      end do
  !      vdm_diag(is) = gz_imt_dens_constr_gzproj(is,is)
  !   end do
  !   !

  !   itt=imt_t2it(time,0.5d0*tstep)
  !   Uloc=Uloc_t(:,itt)
  !   Ust =Ust_t(itt)
  !   Jh=Jh_t(itt)
  !   Jsf=Jsf_t(itt)
  !   Jph=Jph_t(itt)
  !   eLevels = eLevels_t(:,itt)
  !   call get_local_hamiltonian_trace(eLevels)  

  !   Eloc = trace_phi_basis(gzproj,phi_traces_basis_Hloc)
  !   !
  !   gz_imt_local_angular_momenta(1) = trace_phi_basis(gzproj,phi_traces_basis_spin2(:,:))
  !   gz_imt_local_angular_momenta(2) = trace_phi_basis(gzproj,phi_traces_basis_spinZ(:,:))
  !   gz_imt_local_angular_momenta(3) = trace_phi_basis(gzproj,phi_traces_basis_isoSpin2(:,:))
  !   gz_imt_local_angular_momenta(4) = trace_phi_basis(gzproj,phi_traces_basis_isoSpinZ(:,:))
  !   !
  !   gz_imt_Rhop = hopping_renormalization_normal(gzproj,vdm_diag)
  !   gz_imt_Qhop = hopping_renormalization_anomalous(gzproj,vdm_diag)
  !   !
  !   Rhop=gz_imt_Rhop
  !   Qhop=gz_imt_Qhop
  !   !
  !   do is=1,Ns
  !      do js=1,Ns
  !         Rhop_dag(is,js) = conjg(Rhop(js,is))
  !         Qhop_dag(is,js) = conjg(Qhop(js,is))
  !      end do
  !   end do
  !   !
  !   gz_imt_unitary_constr = 0.d0
  !   do iphi=1,Nphi
  !      gz_imt_unitary_constr = gz_imt_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
  !   end do
  !   !+- SLATER
  !   Estar=0.d0
  !   gz_imt_dens_constr_slater=0.d0
  !   gz_imt_dens_constrA_slater=0.d0
  !   do ik=1,Lk
  !      call get_Hk_t(Hk,ik,time)
  !      !+- define Hk_renormalized -+!
  !      Hk_tmp=matmul(Hk,Rhop)
  !      Hk_tmp=matmul(Rhop_dag,Hk_tmp)
  !      !
  !      Hk_tmp=matmul(Hk,Qhop)
  !      Hk_tmp=matmul(Rhop_dag,Hk_tmp)
  !      !
  !      Hk_tmp=matmul(Hk,Rhop)
  !      Hk_tmp=matmul(Qhop_dag,Hk_tmp)
  !      !
  !      Hk_tmp=matmul(Hk,Qhop)
  !      Hk_tmp=matmul(Qhop_dag,Hk_tmp)
  !      !       
  !      do is=1,Ns
  !         do js=1,Ns
  !            !
  !            gz_imt_dens_constr_slater(is,js) = gz_imt_dens_constr_slater(is,js) + slater(1,is,js,ik)*wtk(ik)
  !            gz_imt_dens_constrA_slater(is,js) = gz_imt_dens_constrA_slater(is,js) + slater(2,is,js,ik)*wtk(ik)
  !            !
  !            do iis=1,Ns
  !               do jjs=1,Ns
  !                  !
  !                  Estar = Estar + conjg(gz_imt_Rhop(iis,is))*Hk(iis,jjs)*gz_imt_Rhop(jjs,js)*slater(1,is,js,ik)*wtk(ik)
  !                  Estar = Estar + conjg(gz_imt_Rhop(iis,is))*Hk(iis,jjs)*gz_imt_Qhop(jjs,js)*slater(2,is,js,ik)*wtk(ik)
  !                  Estar = Estar + conjg(gz_imt_Qhop(iis,is))*Hk(iis,jjs)*gz_imt_Rhop(jjs,js)*conjg(slater(2,js,is,ik))*wtk(ik)
  !                  Estar = Estar - conjg(gz_imt_Qhop(iis,is))*Hk(iis,jjs)*gz_imt_Qhop(jjs,js)*slater(1,js,is,ik)*wtk(ik)                   
  !                  if(is.eq.js) then
  !                     Estar = Estar + conjg(gz_imt_Qhop(iis,is))*Hk(iis,jjs)*gz_imt_Qhop(jjs,js)*wtk(ik)
  !                  end if
  !                  !
  !               end do
  !            end do
  !            !
  !         end do
  !      end do
  !   end do
  !   !
  !   gz_imt_energies(1) = Estar+Eloc
  !   gz_imt_energies(2) = Estar
  !   gz_imt_energies(3) = Eloc
  !   !
  ! end subroutine gz_imt_measure_superc




  ! subroutine gz_imt_measure_superc_sp(psi_t,time,read_slater,read_gzproj)
  !   implicit none
  !   complex(8),dimension(nDynamics)   :: psi_t
  !   real(8)                           :: time
  !   complex(8),dimension(Nphi),optional   :: read_gzproj
  !   complex(8),dimension(2,Ns,Ns,Lk),optional  :: read_slater
  !   complex(8),dimension(2,Ns,Ns,Lk)  :: slater
  !   complex(8),dimension(3,Ns,Ns,Lk)  :: slater_
  !   complex(8),dimension(Ns,Ns)     :: Hk,Hk_tmp
  !   complex(8),dimension(Ns,Ns)     :: Rhop,Qhop,Rhop_dag,Qhop_dag
  !   complex(8),dimension(2*Ns,2*Ns)     :: Hks
  !   real(8),dimension(2*Ns)         :: eks
  !   complex(8),dimension(Nphi)      :: gzproj,gztmp
  !   real(8)                         :: Estar,Eloc,Egz
  !   real(8),dimension(Ns)           :: vdm_diag
  !   real(8)                         :: nqp
  !   integer                         :: is,js,ik,it,iphi,iis,jjs
  !   !
  !   it=imt_t2it(time,tstep)
  !   !
  !   call dynamicalVector_2_wfMatrix_superc(psi_t,slater,gzproj)
  !   if(present(read_slater)) read_slater=slater
  !   if(present(read_gzproj)) read_gzproj=gzproj
  !   slater_(1:2,:,:,:) = slater
  !   slater_(3,:,:,:) = zero
  !   do is=1,Ns
  !      slater_(3,is,is,:) = one 
  !      do js=1,Ns
  !         slater_(3,is,js,:) = slater_(3,is,js,:) - slater(1,js,is,:)
  !      end do
  !   end do
  !   !
  !   gz_imt_local_density_matrix=0.d0
  !   gz_imt_local_dens_dens=0.d0
  !   gz_imt_dens_constr_gzproj=0.d0
  !   gz_imt_dens_constrA_gzproj=0.d0
  !   gz_imt_local_sc_order=0.d0
  !   vdm_diag=0.d0
  !   !
  !   Estar=0.d0
  !   do is=1,Ns
  !      do js=1,Ns
  !         gz_imt_local_density_matrix(is,js) = &
  !              trace_phi_basis_sp(gzproj,phi_spTraces_basis_local_dens(is,js))          
  !         gz_imt_local_dens_dens(is,js) = &
  !              trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens_dens(is,js))
  !         !
  !         !
  !         gz_imt_dens_constr_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens(is,js))
  !         gz_imt_dens_constrA_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens_anomalous(is,js))
  !         !
  !         gz_imt_local_sc_order(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_sc_order(is,js))
  !      end do
  !      vdm_diag(is) = gz_imt_dens_constr_gzproj(is,is)
  !   end do
  !   !
  !   gztmp=zero
  !   gztmp = get_local_hamiltonian_HLOCphi(gzproj,eLevels)
  !   Eloc = 0.d0
  !   do iphi=1,Nphi
  !      Eloc = Eloc + conjg(gzproj(iphi))*gztmp(iphi)
  !   end do
  !   !
  !   gz_imt_local_angular_momenta=0.d0
  !   gz_imt_local_angular_momenta(1) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_spin2)
  !   gz_imt_local_angular_momenta(2) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_spinZ)
  !   gz_imt_local_angular_momenta(3) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_isoSpin2)
  !   gz_imt_local_angular_momenta(4) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_isoSpinZ)
  !   !
  !   write(*,*) 'measuring using spMv representation'
  !   gz_imt_Rhop=zero
  !   gz_imt_Rhop = hopping_renormalization_normal_sp(gzproj,vdm_diag)   
  !   gz_imt_Qhop=zero
  !   gz_imt_Qhop = hopping_renormalization_anomalous_sp(gzproj,vdm_diag)
  !   !
  !   !
  !   Rhop=gz_imt_Rhop
  !   Qhop=gz_imt_Qhop
  !   !
  !   !
  !   do is=1,Ns
  !      do js=1,Ns
  !         Rhop_dag(is,js) = conjg(Rhop(js,is))
  !         Qhop_dag(is,js) = conjg(Qhop(js,is))
  !      end do
  !   end do
  !   !
  !   gz_imt_unitary_constr = 0.d0
  !   do iphi=1,Nphi
  !      gz_imt_unitary_constr = gz_imt_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
  !   end do
  !   !
  !   !+- SLATER
  !   Estar=0.d0
  !   gz_imt_dens_constr_slater=0.d0
  !   gz_imt_dens_constrA_slater=0.d0
  !   do ik=1,Lk
  !      call get_Hk_t(Hk,ik,time)
  !      !+- define Hk_renormalized -+!
  !      Hk_tmp=matmul(Hk,Rhop)
  !      Hk_tmp=matmul(Rhop_dag,Hk_tmp)
  !      !Hks(1:Ns,1:Ns) = Hk_tmp 
  !      !
  !      Hk_tmp=matmul(Hk,Qhop)
  !      Hk_tmp=matmul(Rhop_dag,Hk_tmp)
  !      !Hks(1:Ns,Ns+1:2*Ns) = Hk_tmp 
  !      !
  !      Hk_tmp=matmul(Hk,Rhop)
  !      Hk_tmp=matmul(Qhop_dag,Hk_tmp)
  !      !Hks(Ns+1:2*Ns,1:Ns) = Hk_tmp
  !      !
  !      Hk_tmp=matmul(Hk,Qhop)
  !      Hk_tmp=matmul(Qhop_dag,Hk_tmp)
  !      !Hks(Ns+1:2*Ns,Ns+1:2*Ns) = Hk_tmp
  !      !       
  !      do is=1,Ns
  !         do js=1,Ns
  !            !
  !            gz_imt_dens_constr_slater(is,js) = gz_imt_dens_constr_slater(is,js) + slater(1,is,js,ik)*wtk(ik)
  !            gz_imt_dens_constrA_slater(is,js) = gz_imt_dens_constrA_slater(is,js) + slater(2,is,js,ik)*wtk(ik)
  !            !
  !            do iis=1,Ns
  !               do jjs=1,Ns
  !                  !
  !                  Estar = Estar + conjg(gz_imt_Rhop(iis,is))*Hk(iis,jjs)*gz_imt_Rhop(jjs,js)*slater(1,is,js,ik)*wtk(ik)
  !                  Estar = Estar + conjg(gz_imt_Rhop(iis,is))*Hk(iis,jjs)*gz_imt_Qhop(jjs,js)*slater(2,is,js,ik)*wtk(ik)
  !                  Estar = Estar + conjg(gz_imt_Qhop(iis,is))*Hk(iis,jjs)*gz_imt_Rhop(jjs,js)*conjg(slater(2,js,is,ik))*wtk(ik)
  !                  Estar = Estar - conjg(gz_imt_Qhop(iis,is))*Hk(iis,jjs)*gz_imt_Qhop(jjs,js)*slater(1,js,is,ik)*wtk(ik)                   
  !                  if(is.eq.js) then
  !                     Estar = Estar + conjg(gz_imt_Qhop(iis,is))*Hk(iis,jjs)*gz_imt_Qhop(jjs,js)*wtk(ik)
  !                  end if
  !                  !
  !               end do
  !            end do
  !            !
  !         end do
  !      end do
  !   end do
  !   !
  !   gz_imt_energies(1) = Estar+Eloc
  !   gz_imt_energies(2) = Estar
  !   gz_imt_energies(3) = Eloc
  !   !
  ! end subroutine gz_imt_measure_superc_sp








  ! subroutine gz_imt_measure_constr_superc_sp(psi_t,time)
  !   complex(8),dimension(nDynamics) :: psi_t
  !   real(8) :: time
  !   complex(8),dimension(2,Ns,Ns,Lk) :: slater
  !   complex(8),dimension(3,Ns,Ns,Lk) :: slater_
  !   complex(8),dimension(Ns,Ns) :: Hk,Hk_tmp
  !   complex(8),dimension(Ns,Ns) :: Rhop,Qhop,Rhop_dag,Qhop_dag
  !   complex(8),dimension(2*Ns,2*Ns) :: Hks
  !   real(8),dimension(2*Ns) ::eks
  !   complex(8),dimension(Nphi) :: gzproj
  !   real(8) :: Estar,Eloc,Egz
  !   real(8),dimension(Ns) :: vdm_diag
  !   real(8) :: nqp

  !   integer                         :: is,js,ik,it,iphi,iis,jjs

  !   it=imt_t2it(time,tstep)
  !   !
  !   call dynamicalVector_2_wfMatrix_superc(psi_t,slater,gzproj)  
  !   slater_(1:2,:,:,:) = slater
  !   slater_(3,:,:,:) = zero
  !   do is=1,Ns
  !      slater_(3,is,is,:) = one 
  !      do js=1,Ns
  !         slater_(3,is,js,:) = slater_(3,is,js,:) - slater(1,js,is,:)
  !      end do
  !   end do
  !   !

  !   !
  !   do is=1,Ns
  !      do js=1,Ns
  !         !
  !         gz_imt_dens_constr_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens(is,js))
  !         gz_imt_dens_constrA_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens_anomalous(is,js))          
  !         !
  !      end do
  !      vdm_diag(is) = gz_imt_dens_constr_gzproj(is,is)
  !   end do
  !   !
  !   gz_imt_Rhop = hopping_renormalization_normal_sp(gzproj,vdm_diag)
  !   gz_imt_Qhop = hopping_renormalization_anomalous_sp(gzproj,vdm_diag)
  !   !
  !   Rhop=gz_imt_Rhop
  !   Qhop=gz_imt_Qhop
  !   !
  !   do is=1,Ns
  !      do js=1,Ns
  !         Rhop_dag(is,js) = conjg(Rhop(js,is))
  !         Qhop_dag(is,js) = conjg(Qhop(js,is))
  !      end do
  !   end do
  !   !
  !   gz_imt_unitary_constr = 0.d0
  !   do iphi=1,Nphi
  !      gz_imt_unitary_constr = gz_imt_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
  !   end do

  !   !+- SLATER
  !   Estar=0.d0
  !   gz_imt_dens_constr_slater=0.d0
  !   gz_imt_dens_constrA_slater=0.d0
  !   do ik=1,Lk
  !      call get_Hk_t(Hk,ik,time)
  !      !+- define Hk_renormalized -+!
  !      Hk_tmp=matmul(Hk,Rhop)
  !      Hk_tmp=matmul(Rhop_dag,Hk_tmp)
  !      !
  !      Hk_tmp=matmul(Hk,Qhop)
  !      Hk_tmp=matmul(Rhop_dag,Hk_tmp)
  !      !
  !      Hk_tmp=matmul(Hk,Rhop)
  !      Hk_tmp=matmul(Qhop_dag,Hk_tmp)
  !      !
  !      Hk_tmp=matmul(Hk,Qhop)
  !      Hk_tmp=matmul(Qhop_dag,Hk_tmp)
  !      !
  !      do is=1,Ns
  !         do js=1,Ns
  !            !
  !            gz_imt_dens_constr_slater(is,js) = gz_imt_dens_constr_slater(is,js) + slater(1,is,js,ik)*wtk(ik)
  !            gz_imt_dens_constrA_slater(is,js) = gz_imt_dens_constrA_slater(is,js) + slater(2,is,js,ik)*wtk(ik)
  !            !
  !         end do
  !      end do
  !   end do
  !   !
  ! end subroutine gz_imt_measure_constr_superc_sp




  subroutine step_imt_dynamics(nsys,tstep,t,yt,itd_lgrNC,itd_lgrU,eom_funct) 
    integer                                   :: nsys
    real(8)                                   :: tstep,t
    complex(8),dimension(nsys),intent(inout)  :: yt
    complex(8),dimension(nsys)                :: yt_old,yt_new
    complex(8),dimension(Ns,Ns)                ::lgrNC_old,lgrNC_new
    complex(8),dimension(Ns,Ns),intent(inout) :: itd_lgrNC
    real(8),intent(inout)                     :: itd_lgrU

    real(8),dimension(:),allocatable          :: lgr,delta_out
    complex(8),dimension(:),allocatable       :: lgr_cmplx
    integer                                   :: iter,Nopt
    integer                                   :: i,i0
    real(8)                                   :: delta
    interface
       function eom_funct(t,y,Nsys)
         implicit none
         integer                              :: Nsys
         real(8)                              :: t   
         complex(8),dimension(Nsys)           :: eom_funct
         complex(8),dimension(Nsys)           :: y
       end function eom_funct
    end interface
    !
    !+- defines the initial time of the step (used for computing the increment integral)
    imt_tstep = t
    !+- number of lagrange params to be optimized -+!
    Nopt = Nvdm_NC_opt
    allocate(lgr_cmplx(Nvdm_NC_opt))
    allocate(lgr(2*Nopt));allocate(delta_out(2*Nopt))
    !
    call vdm_NC_stride_m2v(itd_lgrNC,lgr_cmplx)    
    do i=1,Nvdm_NC_opt
       lgr(i) = dreal(lgr_cmplx(i))
       lgr(i+Nvdm_NC_opt) = dimag(lgr_cmplx(i))
    end do
    !
    yt_old = yt

    lgrNC_old = itd_lgrNC
    if(.not.allocated(imt_lgrNC)) allocate(imt_lgrNC(Ns,Ns))
    imt_lgrNC=itd_lgrNC
    !
    ! delta_out = fix_anomalous_vdm(lgr)
    ! delta=0.d0
    ! do i=1,2*Nopt
    !    delta = delta + delta_out(i)**2.d0
    ! end do   
    !    
    yt_new = RK4_step(nDynamics,4,tstep,t,yt_old,eom_funct)
    yt = yt_new

    ! select case(lgr_method)
    ! case('CG_min')
    !    call fmin_cgminimize(lgr,imt_fix_constr_,iter,delta,ftol=1.d-04,itmax=20)    
    ! case('f_zero')
    !    call fsolve(imt_fix_constr,lgr,tol=1.d-04,info=iter)    
    ! case default
    !    call fsolve(imt_fix_constr,lgr,tol=1.d-04,info=iter)    
    ! end select
    ! !
    ! delta_out = imt_fix_constr(lgr)
    ! delta=0.d0
    ! do i=1,2*Nopt
    !    delta = delta + delta_out(i)**2.d0
    ! end do
    ! !
    ! write(*,*) 'Immaginary Time Dynamics -> lagrange parameters'
    ! write(*,*) lgr
    ! write(*,*) delta
    ! write(*,*)
    ! !
    ! yt=yt_new
    ! !
    ! lgr_cmplx=zero
    ! do i=1,Nvdm_AC_opt
    !    lgr_cmplx(i) = lgr(i)+xi*lgr(i+Nvdm_NC_opt)
    ! end do
    ! call vdm_NC_stride_v2m(lgr_cmplx,itd_lgrNC)
    ! itd_lgrU = imt_lgrU !(computed in the derivative)
    ! !
    ! !
    ! lgrNC_new = itd_lgrNC    
    lgrNC_new = lgrNC_old
    write(600,'(10F18.10)') t,gz_imt_qpH(10,1,1),gz_imt_qpH(30,1,1)
    call update_qpH_integral(yt_old,yt_new,lgrNC_old,lgrNC_new)
    write(601,'(10F18.10)') t,gz_imt_qpH(10,1,1),gz_imt_qpH(30,1,1)
    !    
    !
  contains

    function imt_fix_constr(lgr) result(delta)
      real(8),dimension(:) :: lgr
      real(8),dimension(size(lgr)) :: delta
      complex(8),dimension(:),allocatable  :: lgr_cmplx,delta_cmplx
      complex(8),dimension(Ns,Ns) :: normal_constrGZ,normal_constrSL,delta_constr
      real(8) :: unitary_constrGZ
      integer :: i0,i,is,js
      real(8) :: tmp_test
      !
      if(allocated(imt_lgrNC)) deallocate(imt_lgrNC)
      allocate(imt_lgrNC(Ns,Ns)); 
      imt_lgrNC=zero
      !
      !+- dump slater_lgr_multipliers
      allocate(lgr_cmplx(Nvdm_NC_opt))
      do i=1,Nvdm_NC_opt
         lgr_cmplx(i) = lgr(i) + xi*lgr(i+Nvdm_NC_opt)
      end do
      call vdm_NC_stride_v2m(lgr_cmplx,imt_lgrNC)
      !
      yt_new = RK4_step(nDynamics,4,itstep,t,yt_old,eom_funct)
      !
      call gz_imt_measure_constr(yt_new,t) 
      !
      do is=1,Ns
         do js=1,Ns
            call get_imt_dens_constr_slater(is,js,normal_constrSL(is,js))
            call get_imt_dens_constr_gzproj(is,js,normal_constrGZ(is,js))
            delta_constr(is,js) = normal_constrGZ(is,js)-normal_constrSL(is,js)
         end do
      end do
      call get_imt_unitary_constr(unitary_constrGZ)
      !
      !
      !Re-construct solution
      !      
      !
      delta=0.d0
      allocate(delta_cmplx(Nvdm_NC_opt))
      call vdm_NC_stride_m2v(delta_constr,delta_cmplx)
      do i=1,Nvdm_NC_opt
         delta(i) = dreal(delta_cmplx(i))
         delta(i+Nvdm_NC_opt) = dimag(delta_cmplx(i))
      end do
      deallocate(delta_cmplx)
      ! tmp_test=0.d0
      ! do i=1,size(lgr)
      !    tmp_test=tmp_test+delta(i)**2.d0
      ! end do
      if(GZneq_verbose) then
         write(*,*) 'imt-dep LAGRANGE'
         do is=1,Ns
            write(*,'(20F18.10)') imt_lgrNC(is,:)         
         end do
         write(*,'(20F18.10)') imt_lgrU         
         write(*,*)
         write(*,*) 'deviation from constraint conservation'
         write(*,'(20F18.10)') delta
      end if
      !
    end function imt_fix_constr
    !


    function imt_fix_constr_(lgr) result(delta)
      real(8),dimension(:) :: lgr
      real(8) :: delta
      complex(8),dimension(:),allocatable  :: lgr_cmplx,delta_cmplx
      complex(8),dimension(Ns,Ns) :: normal_constrGZ,normal_constrSL,delta_constr
      real(8) :: unitary_constrGZ
      integer :: i0,i,is,js
      real(8) :: tmp_test
      !
      if(allocated(imt_lgrNC)) deallocate(imt_lgrNC)
      allocate(imt_lgrNC(Ns,Ns)); 
      imt_lgrNC=zero;imt_lgrU = 0.d0
      !
      !+- dump slater_lgr_multipliers
      allocate(lgr_cmplx(Nvdm_NC_opt))
      do i=1,Nvdm_NC_opt
         lgr_cmplx(i) = lgr(i) + xi*lgr(i+Nvdm_NC_opt)
      end do
      call vdm_NC_stride_v2m(lgr_cmplx,imt_lgrNC)
      !
      yt_new = RK_step(nDynamics,4,itstep,t,yt_old,eom_funct)
      !
      call gz_imt_measure_constr(yt_new,t) 
      !
      do is=1,Ns
         do js=1,Ns
            call get_imt_dens_constr_slater(is,js,normal_constrSL(is,js))
            call get_imt_dens_constr_gzproj(is,js,normal_constrGZ(is,js))
            delta_constr(is,js) = normal_constrGZ(is,js)-normal_constrSL(is,js)
         end do
      end do
      call get_imt_unitary_constr(unitary_constrGZ)
      !
      !
      !Re-construct solution
      !      
      !
      delta=0.d0
      do is=1,Ns
         do js=1,Ns
            delta = delta + delta_constr(is,js)*conjg(delta_constr(is,js))
         end do
      end do
      if(GZneq_verbose) then
         write(*,*) 'imt-dep LAGRANGE CG_min'
         do is=1,Ns
            write(*,'(20F18.10)') imt_lgrNC(is,:)         
         end do
         write(*,'(20F18.10)') imt_lgrU         
         write(*,*)
         write(*,*) 'deviation from constraint conservation'
         write(*,'(20F18.10)') delta
      end if
      !
    end function imt_fix_constr_

  end subroutine step_imt_dynamics







  function imt_t2it(time,delta_t,t0_) result(it)
    real(8) :: time,delta_t
    real(8),optional :: t0_
    real(8) :: t0
    integer :: it
    real(8) :: eps=1.d-8
    t0=beta_init
    if(present(t0_)) t0=t0_
    it = (time-t0+eps*delta_t)/delta_t + 1
  end function imt_t2it





END MODULE GZ_IMTDYNAMICS





  ! subroutine step_imt_dynamics(nsys,tstep,t,yt,itd_lgr,eom_funct) 
  !   integer :: nsys
  !   real(8) :: tstep,t
  !   complex(8),dimension(nsys),intent(inout) :: yt
  !   complex(8),dimension(nsys) :: yt_old,yt_new
  !   complex(8),dimension(Ns,Ns),intent(inout) :: itd_lgr
  !   real(8),dimension(:),allocatable            :: lgr,delta_out
  !   complex(8),dimension(:),allocatable :: lgr_cmplx
  !   integer :: iter,Nopt
  !   integer :: i,i0
  !   real(8) :: delta
  !   interface
  !      function eom_funct(t,y,Nsys)
  !        implicit none
  !        integer                    :: Nsys
  !        real(8)                    :: t   
  !        complex(8),dimension(Nsys) :: eom_funct
  !        complex(8),dimension(Nsys) :: y
  !      end function eom_funct
  !   end interface
  !   !+- number of lagrange params to be optimized -+!
  !   Nopt = 2*Nvdm_NC_opt
  !   allocate(lgr_cmplx(Nvdm_NC_opt))
  !   allocate(lgr(2*Nopt+1));allocate(delta_out(2*Nopt+1))
  !   !
  !   call vdm_NC_stride_m2v(itd_lgr(:,:),lgr_cmplx)    
  !   do i=1,Nvdm_NC_opt
  !      lgr(i) = dreal(lgr_cmplx(i))
  !      lgr(i+Nvdm_NC_opt) = dimag(lgr_cmplx(i))
  !   end do    
  !   !
  !   yt_old = yt
  !   !
  !   ! delta_out = fix_anomalous_vdm(lgr)
  !   ! delta=0.d0
  !   ! do i=1,2*Nopt
  !   !    delta = delta + delta_out(i)**2.d0
  !   ! end do   
  !   !
  !   select case(lgr_method)
  !   case('CG_min')
  !      call fmin_cgminimize(lgr,fix_anomalous_vdm_,iter,delta,ftol=1.d-04,itmax=20)    
  !   case('f_zero')
  !      call fsolve(fix_anomalous_vdm,lgr,tol=1.d-04,info=iter)    
  !   case default
  !      call fsolve(fix_anomalous_vdm,lgr,tol=1.d-04,info=iter)    
  !   end select
  !   !
  !   delta_out = fix_anomalous_vdm(lgr)
  !   delta=0.d0
  !   do i=1,2*Nopt
  !      delta = delta + delta_out(i)**2.d0
  !   end do
  !   !
  !   write(*,*) 'Time Dependent lagrange parameters: error'
  !   write(*,*) delta
  !   write(*,*)
  !   !
  !   yt=yt_new
  !   !
  !   lgr_cmplx=zero
  !   do i=1,Nvdm_NC_opt
  !      lgr_cmplx(i) = lgr(i)+xi*lgr(i+Nvdm_NC_opt)
  !   end do
  !   call vdm_NC_stride_v2m(lgr_cmplx,td_lgr(1,:,:))
  !   lgr_cmplx(+2*Nvdm_NC_opt)

  !   i0=2*Nvdm_AC_opt
  !   do i=1,Nvdm_AC_opt
  !      lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_AC_opt)
  !   end do
  !   call vdm_AC_stride_v2m(lgr_cmplx,td_lgr(2,:,:))  
  !   !    
  ! contains

  !   function fix_anomalous_vdm(lgr) result(delta)
  !     real(8),dimension(:) :: lgr
  !     real(8),dimension(size(lgr)) :: delta
  !     complex(8),dimension(:),allocatable  :: lgr_cmplx,delta_cmplx
  !     complex(8),dimension(Ns,Ns) :: anomalous_constrGZ,anomalous_constrSL
  !     integer :: i0,i,is,js
  !     real(8) :: tmp_test
  !     !
  !     if(allocated(neq_lgr)) deallocate(neq_lgr)
  !     allocate(neq_lgr(2,Ns,Ns)); neq_lgr=zero
  !     !
  !     !+- dump slater_lgr_multipliers
  !     allocate(lgr_cmplx(Nvdm_AC_opt))
  !     do i=1,Nvdm_AC_opt
  !        lgr_cmplx(i) = lgr(i) + xi*lgr(i+Nvdm_AC_opt)
  !     end do
  !     call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(1,:,:))
  !     i0 = 2*Nvdm_AC_opt
  !     !+- dump gzproj_lgr_multipliers
  !     do i=1,Nvdm_AC_opt
  !        lgr_cmplx(i) = lgr(i+i0) + xi*lgr(i+i0+Nvdm_AC_opt)
  !     end do
  !     call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(2,:,:))
  !     !

  !     !+- major things to change are here -+!

      
  !     yt_new = RK_step(nDynamics,4,tstep,t,yt_old,eom_funct)
  !     !
  !     call gz_imt_measure_constr_superc_sp(yt_new,t) 
  !     !
  !     do is=1,Ns
  !        do js=1,Ns
  !           call get_imt_dens_constrA_slater(is,js,anomalous_constrSL(is,js))
  !           call get_imt_dens_constrA_gzproj(is,js,anomalous_constrGZ(is,js))
  !        end do
  !     end do
  !     !
  !     delta=0.d0
  !     allocate(delta_cmplx(Nvdm_AC_opt))
  !     call vdm_AC_stride_m2v(anomalous_constrSL,delta_cmplx)
  !     do i=1,Nvdm_AC_opt
  !        delta(i) = dreal(delta_cmplx(i))
  !        delta(i+Nvdm_AC_opt) = dimag(delta_cmplx(i))
  !     end do
  !     deallocate(delta_cmplx)
  !     i0 = 2*Nvdm_AC_opt
  !     allocate(delta_cmplx(Nvdm_AC_opt))
  !     call vdm_AC_stride_m2v(anomalous_constrGZ,delta_cmplx)
  !     do i=1,Nvdm_AC_opt
  !        delta(i0+i) = dreal(delta_cmplx(i))
  !        delta(i0+i+Nvdm_AC_opt) = dimag(delta_cmplx(i))       
  !     end do
  !     deallocate(delta_cmplx)
  !     tmp_test=0.d0
  !     do i=1,size(lgr)
  !        tmp_test=tmp_test+delta(i)**2.d0
  !     end do
  !     write(*,*) 't-dep LAGRANGE fsolve'
  !     do is=1,Ns
  !        write(*,'(20F18.10)') neq_lgr(1,is,:)         
  !     end do
  !     write(*,*) 
  !     do is=1,Ns
  !        write(*,'(20F18.10)') neq_lgr(2,is,:)         
  !     end do
  !     write(*,*) 'deviation from constraint conservation'
  !     write(*,'(20F18.10)') delta
  !     !
  !   end function fix_anomalous_vdm
  !   !
  !   function fix_anomalous_vdm_(lgr) result(delta)
  !     real(8),dimension(:) :: lgr
  !     real(8) :: delta
  !     complex(8),dimension(:),allocatable  :: lgr_cmplx,delta_cmplx
  !     complex(8),dimension(Ns,Ns) :: anomalous_constrGZ,anomalous_constrSL
  !     integer :: i0,i,is,js
  !     !
  !     if(allocated(neq_lgr)) deallocate(neq_lgr)
  !     allocate(neq_lgr(2,Ns,Ns)); neq_lgr=zero

  !     !+- dump slater_lgr_multipliers
  !     allocate(lgr_cmplx(Nvdm_AC_opt))
  !     do i=1,Nvdm_AC_opt
  !        lgr_cmplx(i) = lgr(i) + xi*lgr(i+Nvdm_AC_opt)
  !     end do
  !     call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(1,:,:))
  !     i0 = 2*Nvdm_AC_opt
  !     !+- dump gzproj_lgr_multipliers
  !     do i=1,Nvdm_AC_opt
  !        lgr_cmplx(i) = lgr(i+i0) + xi*lgr(i+i0+Nvdm_AC_opt)
  !     end do
  !     call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(2,:,:))
  !     !
  !     !yt_new = RK_step(nDynamics,4,tstep,t,yt_old,gz_equations_of_motion_superc_lgr_sp)
  !     yt_new = RK_step(nDynamics,4,tstep,t,yt_old,eom_funct)
  !     !
  !     call gz_imt_measure_constr_superc_sp(yt_new,t)
  !     !
  !     do is=1,Ns
  !        do js=1,Ns
  !           call get_imt_dens_constrA_slater(is,js,anomalous_constrSL(is,js))
  !           call get_imt_dens_constrA_gzproj(is,js,anomalous_constrGZ(is,js))
  !        end do
  !     end do
  !     !
  !     delta=0.d0
  !     do is=1,Ns
  !        do js=1,Ns
  !           delta = delta + anomalous_constrSL(is,js)*conjg(anomalous_constrSL(is,js))
  !           delta = delta + anomalous_constrGZ(is,js)*conjg(anomalous_constrGZ(is,js))
  !        end do
  !     end do

  !     write(*,*) 'deviation from constraint conservation',delta
  !     !
  !   end function fix_anomalous_vdm_

  ! end subroutine step_imt_dynamics
