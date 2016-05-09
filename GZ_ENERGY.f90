MODULE GZ_ENERGY_MINIMIZATION 
  ! scifor routines
  USE SF_LINALG
  USE SF_IOTOOLS
  USE SF_SPECIAL
  USE SF_OPTIMIZE
  USE DMFT_MISC
  ! lancelot routines
  USE LANCELOT_simple_double
  ! GZ routines
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_AUX_FUNX
  USE GZ_MATRIX_BASIS
  USE GZ_EFFECTIVE_HOPPINGS
  USE MIN_AMOEBA
  !
  USE MATRIX_SPARSE
  !
  implicit none
  private

  public :: gz_get_energy_vdm
  !
  public :: gz_energy_vdm
  !
  !
  interface gz_proj_minimization_lgr
     module procedure gz_projectors_minimization_nlep
     module procedure gz_projectors_minimization_cmin
  end interface gz_proj_minimization_lgr
  public :: gz_proj_minimization_lgr
  public :: gz_proj_minimization_fixed_lgr         
  public :: gz_proj_minimization_fixed_lgr_hop     
  !
  public :: gz_proj_minimization_lgr_superc 
  public :: gz_proj_minimization_fixed_lgr_superc 
  public :: gz_proj_minimization_fixed_lgr_hop_superc 
  !
  

  public :: slater_minimization_lgr
  public :: slater_minimization_fixed_lgr
  !
  public :: slater_minimization_lgr_superc
  public :: slater_minimization_fixed_lgr_superc
  !
  
  !
  public :: get_gz_ground_state
  public :: get_gz_ground_state_superc
  !
  
contains
  !
  include "slater_min_routines.f90"
  include "gzproj_min_routines.f90"
  !
  subroutine get_gz_ground_state(phi_vec)
    complex(8),dimension(Nphi) :: phi_vec
    real(8)                                    :: E_Hstar
    real(8),dimension(Ns,Ns)            :: slater_lgr_multip
    real(8),dimension(Ns) :: n0
    integer                                    :: iorb,jorb,istate,jstate
    integer :: ifock,iphi,jphi,ifock_,is,js
    !
    real(8),dimension(Nphi) :: phi_vector_test    
    !
    if(.not.allocated(GZ_opt_Rhop)) allocate(GZ_opt_Rhop(Ns,Ns))
    if(.not.allocated(GZ_opt_VDM)) allocate(GZ_opt_VDM(Ns,Ns))
    if(.not.allocated(GZ_opt_slater)) allocate(GZ_opt_slater(Ns,Ns,Lk))
    !
    do is=1,Ns
       do js=1,Ns
          GZ_opt_VDM(is,js) = trace_phi_basis(phi_vec,phi_traces_basis_dens(is,js,:,:))
       end do
       n0(is) = GZ_opt_VDM(is,is)
    end do
    GZ_opt_Rhop=hopping_renormalization_normal(GZ_vector,n0)            
    !
    !
    call slater_minimization_fixed_lgr(GZ_opt_Rhop,GZ_opt_slater_lgr,E_Hstar,GZ_opt_VDM,store=slater_store)
    !
    !    
    !+- GET OBSERVABLES -+!
    ! physical density !
    if(.not.allocated(gz_dens)) allocate(gz_dens(Ns))
    if(.not.allocated(gz_dens_matrix)) allocate(gz_dens_matrix(Ns,Ns))
    do istate=1,Ns
       gz_dens(istate) = trace_phi_basis(phi_vec,phi_traces_basis_local_dens(istate,istate,:,:))
       do jstate=1,Ns
          gz_dens_matrix(istate,jstate) = trace_phi_basis(phi_vec,phi_traces_basis_local_dens(istate,jstate,:,:))
       end do
    end do
    ! density-density same orbital -aka orbital doubly occupancy-!
    if(.not.allocated(gz_docc)) allocate(gz_docc(Norb))
    do iorb=1,Norb
       gz_docc(iorb) = trace_phi_basis(phi_vec,phi_traces_basis_docc_orb(iorb,:,:))
    end do
    ! density-density different orbitals !
    if(.not.allocated(gz_dens_dens_orb)) allocate(gz_dens_dens_orb(Norb,Norb))
    do iorb=1,Norb
       do jorb=1,Norb
          gz_dens_dens_orb(iorb,jorb)=trace_phi_basis(phi_vec,phi_traces_basis_dens_dens_orb(iorb,jorb,:,:))
       end do
    end do
    !+-
    gz_spin2  = trace_phi_basis(phi_vec,phi_traces_basis_spin2)
    gz_spinZ  = trace_phi_basis(phi_vec,phi_traces_basis_spinZ)
    gz_isospin2 = trace_phi_basis(phi_vec,phi_traces_basis_isospin2)
    gz_isospinZ = trace_phi_basis(phi_vec,phi_traces_basis_isospinZ)
    !+-
  end subroutine get_gz_ground_state
  


  subroutine get_gz_ground_state_superc(phi_vec)
    complex(8),dimension(Nphi) :: phi_vec
    real(8)                                    :: E_Hstar
    real(8),dimension(Ns) :: n0
    integer                                    :: iorb,jorb,istate,jstate
    integer :: ifock,iphi,jphi,ifock_,is,js
    !
    real(8),dimension(Nphi) :: phi_vector_test    
    !
    if(.not.allocated(GZ_opt_Rhop)) allocate(GZ_opt_Rhop(Ns,Ns))
    if(.not.allocated(GZ_opt_Qhop)) allocate(GZ_opt_Qhop(Ns,Ns))
    !
    if(.not.allocated(GZ_opt_VDM_superc)) allocate(GZ_opt_VDM_superc(2,Ns,Ns))
    if(.not.allocated(GZ_opt_slater_superc)) allocate(GZ_opt_slater_superc(2,Ns,Ns,Lk))
    !
    do is=1,Ns
       do js=1,Ns
          GZ_opt_VDM_superc(1,is,js) = trace_phi_basis(phi_vec,phi_traces_basis_dens(is,js,:,:))
          GZ_opt_VDM_superc(2,is,js) = trace_phi_basis(phi_vec,phi_traces_basis_dens_anomalous(is,js,:,:))
       end do
       n0(is) = GZ_opt_VDM_superc(1,is,is)
    end do
    !
    GZ_opt_Rhop=hopping_renormalization_normal(GZ_vector,n0)            
    GZ_opt_Qhop=hopping_renormalization_anomalous(GZ_vector,n0)            
    !
    !
    call slater_minimization_fixed_lgr_superc(GZ_opt_Rhop,GZ_opt_Qhop,GZ_opt_slater_lgr_superc,E_Hstar,GZ_opt_VDM_superc,store=slater_store)
    !
    !    
    !+- GET OBSERVABLES -+!
    ! physical density !
    if(.not.allocated(gz_dens)) allocate(gz_dens(Ns))
    if(.not.allocated(gz_dens_matrix)) allocate(gz_dens_matrix(Ns,Ns))
    do istate=1,Ns
       gz_dens(istate) = trace_phi_basis(phi_vec,phi_traces_basis_local_dens(istate,istate,:,:))
       do jstate=1,Ns
          gz_dens_matrix(istate,jstate) = trace_phi_basis(phi_vec,phi_traces_basis_local_dens(istate,jstate,:,:))
       end do
    end do
    ! density-density same orbital -aka orbital doubly occupancy-!
    if(.not.allocated(gz_docc)) allocate(gz_docc(Norb))
    do iorb=1,Norb
       gz_docc(iorb) = trace_phi_basis(phi_vec,phi_traces_basis_docc_orb(iorb,:,:))
    end do
    ! density-density different orbitals !
    if(.not.allocated(gz_dens_dens_orb)) allocate(gz_dens_dens_orb(Norb,Norb))
    do iorb=1,Norb
       do jorb=1,Norb
          gz_dens_dens_orb(iorb,jorb)=trace_phi_basis(phi_vec,phi_traces_basis_dens_dens_orb(iorb,jorb,:,:))
       end do
    end do
    !+-
    if(.not.allocated(gz_sc_order)) allocate(gz_sc_order(Ns,Ns))
    do is=1,Ns
       do js=1,Ns
          gz_sc_order(is,js) = trace_phi_basis(phi_vec,phi_traces_basis_sc_order(is,js,:,:))
       end do
    end do
    !+-    
    gz_spin2  = trace_phi_basis(phi_vec,phi_traces_basis_spin2)
    gz_spinZ  = trace_phi_basis(phi_vec,phi_traces_basis_spinZ)
    gz_isospin2 = trace_phi_basis(phi_vec,phi_traces_basis_isospin2)
    gz_isospinZ = trace_phi_basis(phi_vec,phi_traces_basis_isospinZ)
    !+-
  end subroutine get_gz_ground_state_superc
  !
  function gz_energy_vdm(vdm) result(GZenergy)
    real(8),dimension(:),intent(in) :: vdm
    real(8)               :: GZenergy
    select case(min_method)
    case('nlep')
       if(.not.gz_superc) then
          GZenergy=gz_energy_recursive_nlep(vdm)
       else
          GZenergy=gz_energy_recursive_nlep_superc(vdm)
       end if
    case('cmin')
       GZenergy=gz_energy_recursive_cmin(vdm)
    end select
  end function gz_energy_vdm
  !
  subroutine gz_get_energy_vdm(x,GZ_energy,i) 
    real(8),intent(in) :: x(:)
    real(8),intent(out)              :: GZ_energy
    integer,intent(in),optional :: i
    real(8),dimension(Ns) :: vdm
    complex(8),dimension(Ns,Ns) :: Rhop,vdm_matrix
    complex(8),dimension(Ns,Ns) :: slater_derivatives    
    complex(8),dimension(size(x)) :: x_cmplx
    real(8),dimension(Ns)           :: slater_lgr_multip,R_diag
    real(8),dimension(Ns,Ns,2) :: GZproj_lgr_multip  ! 
    real(8),dimension(Ns,Ns) :: GZproj_lgr_multip_  ! 
    real(8)                                :: E_Hstar,E_Hloc
    complex(8),dimension(nPhi)               :: GZvect_iter  ! GZ vector (during iterations)      
    integer :: is,imap
    !
    do is=1,Ns
    end do
    !
    x_cmplx = x
    call vdm_NC_stride_v2m(x_cmplx,vdm_matrix)
    do is=1,Ns
       vdm(is) = dreal(vdm_matrix(is,is))
    end do
    !
    GZ_energy = gz_energy_vdm(vdm)
    !      
  end subroutine gz_get_energy_vdm
  !
  function gz_energy_recursive_nlep(n0)   result(GZ_energy)
    real(8),dimension(:),intent(in) :: n0 !INPUT: Variational Density Matrix (VDM) (diagonal in istate)       
    !
    ! complex(8),dimension(Ns,Ns),optional     :: R_init_        ! initial guess for the NORMAL hopping renormalization matrices
    ! complex(8),dimension(Ns,Ns),optional     :: Q_init_        ! initial guess for the ANOMALOUS hopping renormalization matri
    !
    real(8)                         :: GZ_energy !INPUT: Optimized GZ energy at fixed 
    real(8)                         :: GZ_energy_old,energy_err     ! Value of the GZ energy functional
    !
    complex(8),dimension(Ns,Ns)     :: R_init        ! initial guess for the NORMAL hopping renormalization matrices
    complex(8),dimension(Ns,Ns)     :: R_iter,R_old ! hopping matrix renormalization (during iterations)
    !
    complex(8),dimension(Ns,Ns)     :: slater_derivatives    
    !
    real(8),dimension(Ns)           :: R_diag
    complex(8),dimension(Ns,Ns)        :: slater_lgr_multip  ! 
    complex(8),dimension(Ns,Ns) :: GZproj_lgr_multip
    real(8)                         :: E_Hstar,E_Hloc
    complex(8),dimension(Nphi)      :: GZvect_iter  ! GZ vector (during iterations)
    !
    integer                         :: istate,iter,jstate,ifock,jfock,iphi,jphi,is,js,imap,jmap    
    integer                         :: unit,iorb,jorb,ispin,jspin
    logical                         :: bound
    !
    write(*,*) '********************'
    write(*,*) 'INPUT DENSITY',n0(:)
    bound=.false.
    do istate=1,Ns
       if(n0(istate).lt.0.d0.or.n0(istate).gt.1.d0) bound=.true.
    end do
    !
    if(.not.bound) then
       !+- initialize Rhop and Qhop from seed file (if no seed file initialization to the not interacting values)       
       call init_Rhop_seed(R_init)
       !
       GZ_energy=0.d0 
       R_iter = R_init
       GZvect_iter=0.d0
       !
       do iter=1,Niter_self
          !+- update hopping matrices -+!
          GZ_energy_old=GZ_energy
          R_old = R_iter          
          !+----------------------------+!
          !+- SLATER STEP MINIMIZATION -+!
          !+----------------------------+!    
          call slater_minimization_lgr(R_iter,n0,E_Hstar,slater_lgr_multip, &
               slater_derivatives=slater_derivatives,iverbose=GZmin_verbose)       
          !
          !+----------------------------+!
          !+- GZproj STEP MINIMIZATION -+!
          !+----------------------------+!    
          !
          call gz_proj_minimization_lgr(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,iverbose=GZmin_verbose)
          !
          R_iter=hopping_renormalization_normal(GZvect_iter,n0)
          !
          R_iter=Rmix*R_iter+(1.d0-Rmix)*R_old
          !
          do istate=1,Ns
             R_diag(istate)=R_iter(istate,istate)
          end do
          ! update GZ-energy
          GZ_energy=E_Hstar+E_Hloc
          if(iter.lt.2) then 
             energy_err=1.d0
          else
             energy_err=abs(GZ_energy-GZ_energy_old)
          end if
          if(GZmin_verbose) then
             write(GZmin_unit,'(20F18.10)') dble(iter),energy_err,GZ_energy,E_Hstar,E_Hloc,R_diag(1:Ns)
          end if
          if(energy_err.lt.err_self) exit
       end do
       if(GZmin_verbose) write(GZmin_unit,*) 
       if(GZmin_verbose) write(GZmin_unit,*) 
       if(iter-1.eq.Niter_self) then
          write(*,*) 'Self consistent Gutzwiller minimization'
          write(*,*) 'Input VDM',n0
          write(*,*) 'final error',energy_err
          write(*,*) "Not converged after",Niter_self,'iterations: exiting'
          stop 
       end if
       write(opt_GZ_unit,*) n0
       write(opt_GZ_unit,*)
       !
       do iphi=1,Nphi
          write(opt_GZ_unit,*) GZvect_iter(iphi)
       end do
       !
       write(opt_GZ_unit,*)
       write(opt_GZ_unit,*)
       write(opt_energy_unit,*) n0,GZ_energy,E_Hloc,E_Hstar
       !
       write(opt_rhop_unit,*) n0
       do istate=1,Ns
          write(opt_rhop_unit,'(20F8.4)') dreal(R_iter(istate,:)),dimag(R_iter(istate,:))
       end do
       write(opt_rhop_unit,*)
       !
    else
       GZ_energy=100.d0
    end if
    if(optimization_flag) then
       !+- store final informations to global variables -+!              
       !
       GZ_vector      = GZvect_iter
       GZ_opt_energy  = GZ_energy
       GZ_opt_kinetic = E_Hstar
       GZ_opt_Eloc    = E_Hloc
       GZ_opt_slater_lgr = slater_lgr_multip
       !
    end if
    !
  end function gz_energy_recursive_nlep
  !
  function gz_energy_recursive_nlep_superc(n0)   result(GZ_energy)
    real(8),dimension(:),intent(in) :: n0 !INPUT: Variational Density Matrix (VDM) (diagonal in istate)       
    real(8)                         :: GZ_energy !INPUT: Optimized GZ energy at fixed 
    real(8)                         :: GZ_energy_old,energy_err     ! Value of the GZ energy functional
    !
    complex(8),dimension(Ns,Ns)     :: R_init        ! initial guess for the NORMAL hopping renormalization matrices
    complex(8),dimension(Ns,Ns)     :: Q_init        ! initial guess for the ANOMALOUS hopping renormalization matrices
    complex(8),dimension(Ns,Ns)     :: R_iter,R_old,Q_iter,Q_old ! hopping matrix renormalization (during iterations)
    !
    complex(8),dimension(2,Ns,Ns)     :: slater_derivatives    

    real(8),dimension(Ns)           :: R_diag
    complex(8),dimension(2,Ns,Ns)        :: slater_lgr_multip  ! 
    real(8),dimension(Ns,Ns)      :: tmp_lgr_multip,tmp_gz_lgr
    complex(8),dimension(2,Ns,Ns) :: GZproj_lgr_multip
    real(8)                         :: E_Hstar,E_Hloc
    complex(8),dimension(Nphi)      :: GZvect_iter  ! GZ vector (during iterations)
    !
    integer                         :: istate,iter,jstate,ifock,jfock,iphi,jphi,is,js,imap,jmap    
    integer                         :: unit,iorb,jorb,ispin,jspin
    logical                         :: bound
    !
    write(*,*) '********************'
    write(*,*) 'INPUT DENSITY',n0(:)
    bound=.false.
    do istate=1,Ns
       if(n0(istate).lt.0.d0.or.n0(istate).gt.1.d0) bound=.true.
    end do
    !
    if(.not.bound) then
       !+- initialize Rhop and Qhop from seed file (if no seed file initialization to the not interacting values)
       call init_Rhop_seed(R_init)
       call init_Qhop_seed(Q_init)
       !
       GZ_energy=0.d0    
       R_iter = R_init
       Q_iter = Q_init
       GZvect_iter=0.d0
       do iter=1,Niter_self
          !+- update hopping matrices -+!
          GZ_energy_old=GZ_energy
          R_old = R_iter          
          Q_old = Q_iter
          !+----------------------------+!
          !+- SLATER STEP MINIMIZATION -+!
          !+----------------------------+!    
          !
          call slater_minimization_lgr_superc(R_iter,Q_iter,n0,E_Hstar,slater_lgr_multip, &
               slater_derivatives=slater_derivatives,iverbose=GZmin_verbose)       
          !
          !+----------------------------+!
          !+- GZproj STEP MINIMIZATION -+!
          !+----------------------------+!    
          !
          call gz_proj_minimization_lgr_superc(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,iverbose=GZmin_verbose)
          !
          R_iter=hopping_renormalization_normal(GZvect_iter,n0)
          Q_iter=hopping_renormalization_anomalous(GZvect_iter,n0)
          !
          R_iter=Rmix*R_iter+(1.d0-Rmix)*R_old
          Q_iter=Rmix*Q_iter+(1.d0-Rmix)*Q_old
          !
          do istate=1,Ns
             R_diag(istate)=R_iter(istate,istate)
          end do
          ! update GZ-energy
          GZ_energy=E_Hstar+E_Hloc
          if(iter.lt.2) then 
             energy_err=1.d0
          else
             energy_err=abs(GZ_energy-GZ_energy_old)
          end if
          if(GZmin_verbose) then
             write(GZmin_unit,'(20F18.10)') dble(iter),energy_err,GZ_energy,E_Hstar,E_Hloc,R_diag(1:Ns)
          end if
          if(energy_err.lt.err_self) exit
       end do
       if(GZmin_verbose) write(GZmin_unit,*) 
       if(GZmin_verbose) write(GZmin_unit,*) 
       if(iter-1.eq.Niter_self) then
          write(*,*) 'Self consistent Gutzwiller minimization'
          write(*,*) 'Input VDM',n0
          write(*,*) 'final error',energy_err
          write(*,*) "Not converged after",Niter_self,'iterations: exiting'
          stop 
       end if
       write(opt_GZ_unit,*) n0
       write(opt_GZ_unit,*)
       !
       do iphi=1,Nphi
          write(opt_GZ_unit,*) GZvect_iter(iphi)
       end do
       !
       write(opt_GZ_unit,*)
       write(opt_GZ_unit,*)
       write(opt_energy_unit,*) n0,GZ_energy,E_Hloc,E_Hstar
       !
       write(opt_rhop_unit,*) n0
       do istate=1,Ns
          write(opt_rhop_unit,'(20F8.4)') dreal(R_iter(istate,:)),dimag(R_iter(istate,:))
       end do
       write(opt_rhop_unit,*)
       !
       if(gz_superc) then
          write(opt_qhop_unit,*) n0
          do istate=1,Ns
             write(opt_qhop_unit,'(20F8.4)') dreal(Q_iter(istate,:)),dimag(Q_iter(istate,:))
          end do
          write(opt_qhop_unit,*)
       end if
    else
       GZ_energy=100.d0
    end if
    if(optimization_flag) then
       !+- store final informations to global variables -+!              
       !
       GZ_vector      = GZvect_iter
       GZ_opt_energy  = GZ_energy
       GZ_opt_kinetic = E_Hstar
       GZ_opt_Eloc    = E_Hloc
       GZ_opt_slater_lgr_superc = slater_lgr_multip
       !
    end if
    !
  end function gz_energy_recursive_nlep_superc
  !
  function gz_energy_recursive_cmin(n0)  result(GZ_energy)
    real(8),dimension(:),intent(in)           :: n0 !INPUT: Variational Density Matrix (VDM) (diagonal in istate)    
    real(8)                                   :: GZ_energy !INPUT: Optimized GZ energy at fixed 
    real(8)                                   :: GZ_energy_old,energy_err     ! Value of the GZ energy functional
    complex(8),dimension(Nphi)                  :: GZvect_iter,GZvect_iter_  ! GZ vector (during iterations)
    complex(8),dimension(Ns,Ns)    :: R_iter,R_old ! hopping matrix renormalization (during iterations)

    complex(8),dimension(Ns,Ns)    :: R_init        ! initial guess for the hopping renormalization matrix    
    real(8),dimension(Ns)              :: R_diag
    complex(8),dimension(Ns,Ns,Lk) :: slater_matrix_el    

    complex(8),dimension(Ns,Ns)              :: slater_lgr_multip
    real(8),dimension(Ns,Ns)    :: GZproj_lgr_multip  
    real(8)                                   :: E_Hstar,E_Hloc,tmp_norm
    !
    integer                                   :: istate,iter,jstate,ifock,jfock,i_ind
    integer                                   :: iphi,jphi
    integer                                   :: unit
    logical                                   :: bound
    !
    write(*,*) '*************************'
    write(*,*) 'INPUT DENSITY',n0(:)
    bound=.false.
    do istate=1,Ns
       if(n0(istate).lt.0.d0.or.n0(istate).gt.1.d0) bound=.true.
    end do
    !
    if(.not.bound) then
       R_iter=0.d0
       do istate=1,Ns
          R_iter(istate,istate) = Rseed
       end do
       !
       !call initialize_GZprojectors(GZvect_iter,n0)
       GZ_energy=0.d0
       do iter=1,Niter_self
          !
          !+- update phi_vectors -+!
          !
          GZ_energy_old=GZ_energy
          R_old = R_iter

          !GZvect_iter_ = GZvect_iter

          !
          !+----------------------------+!
          !+- SLATER STEP MINIMIZATION -+!
          !+----------------------------+!    
          !         
          call slater_minimization_lgr(R_iter,n0,E_Hstar,slater_lgr_multip,&
               slater_matrix_el=slater_matrix_el,iverbose=.true.)       
          !
          !+----------------------------+!
          !+- GZproj STEP MINIMIZATION -+!
          !+----------------------------+!    
          !                  
          call gz_proj_minimization_lgr(slater_matrix_el,n0,GZ_energy,GZvect_iter,GZproj_lgr_multip,iverbose=GZmin_verbose)
          !
          R_iter=hopping_renormalization_normal(GZvect_iter,n0)
          R_iter=Rmix*R_iter+(1.d0-Rmix)*R_old
          !GZvect_iter = 0.5d0*GZvect_iter + 0.5d0*GZvect_iter_

          !
          do istate=1,Ns
             R_diag(istate)=R_iter(istate,istate)
          end do
          !
          if(iter.lt.2) then 
             energy_err=1.d0
          else
             energy_err=abs(GZ_energy-GZ_energy_old)
          end if
          E_HLoc=GZ_energy-E_Hstar
          write(*,*) GZ_energy
          if(GZmin_verbose) then
             write(GZmin_unit,*) dble(iter),energy_err,GZ_energy,E_Hstar,E_Hloc,R_diag(1:Ns)
          end if
          if(energy_err.lt.err_self) exit          
       end do
       if(GZmin_verbose) write(GZmin_unit,*) 
       if(iter-1.eq.Niter_self) then
          write(*,*) 'Recursive Gutzwiller minimization using Constrained minimization of the GZ projectors'
          write(*,*) 'Input VDM',n0
          write(*,*) 'final error',energy_err
          write(*,*) "Not converged after",Niter_self,'iterations: exiting'
          stop 
       end if
       write(opt_GZ_unit,*) n0
       write(opt_GZ_unit,*)
       !
       do iphi=1,Nphi
          write(opt_GZ_unit,*) GZvect_iter(iphi)
       end do
       !
       write(opt_GZ_unit,*)
       write(opt_GZ_unit,*)
       write(opt_energy_unit,*) n0,GZ_energy,E_Hstar,E_Hloc
       !
       R_iter=hopping_renormalization_normal(GZvect_iter,n0)
       do istate=1,Ns
          R_diag(istate)=R_iter(istate,istate)
       end do
       write(opt_rhop_unit,*) n0,R_diag(1:Ns)

    else
       !
       GZ_energy=100.d0
       !
    end if
    if(optimization_flag) then
       !+- store final informations to global variables -+!              
       GZ_vector = GZvect_iter
       GZ_opt_energy         = GZ_energy
       GZ_opt_kinetic        = E_Hstar
       GZ_opt_Eloc           = E_Hloc
       GZ_opt_slater_lgr     = Slater_lgr_multip
    end if
    !
  end function gz_energy_recursive_cmin
  !
  subroutine initialize_GZprojectors(GZvect_iter,n0)
    complex(8),dimension(Nphi) :: GZvect_iter
    real(8),dimension(Ns) :: n0
    complex(8),dimension(Ns,Ns) :: R_init        
    complex(8),dimension(Ns,Ns) :: slater_derivatives
    real(8),dimension(Ns,Ns)           :: slater_lgr_multip
    real(8),dimension(Ns,Ns) :: GZproj_lgr_multip  
    real(8) :: E_Hstar,E_HLoc
    logical :: iverbose
    integer :: istate

    iverbose=.false.
    R_init=0.d0
    do istate=1,Ns
       R_init(istate,istate)=1.d0
    end do
    ! call slater_minimization_lgr(R_init,n0,E_Hstar,slater_lgr_multip, &
    !      slater_derivatives=slater_derivatives,iverbose=GZmin_verbose)       
    ! call gz_proj_minimization_lgr(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,ifree=.true.,iverbose=GZmin_verbose)   
  end subroutine initialize_GZprojectors
  !
END MODULE GZ_ENERGY_MINIMIZATION



!+- GALAHAD WANTS THIS FUNCTION HERE...VERY WEIRED FUNCTIONS ALLOCATIONS -+!
subroutine fun ( x, f, i )
  !.............................................................................                                                                               
  real(8), intent( in )   :: x( : )
  real(8), intent( out )  :: f
  integer, intent( in ), optional   :: i  
  if ( .not. present( i ) ) then
     !       the objective function value (user defined)                                                                                                       
     !==============================================================================                                                                           
     f = 100.0d0*(x(2)-x(1)**2)**2 +(1.0d0-x(1))**2                      !                                                                                     
     !==============================================================================                                                                           
  else
     select case ( i )
     case ( 1 )
        !               the equality constraint value (user defined)                                                                                           
	!==============================================================================                                                                        
        f = x(1)+3.0d0*x(2)-3.0d0                                   !                                                                                          
        !==============================================================================                                                                        
     case ( 2 )
	!               the inequality constraint value (user defined)                                                                                         
        !==============================================================================                                                                        
        f = x(1)**2+x(2)**2-4.0d0                                    !                                                                                         
        !==============================================================================       
     end select
  end if
  return
end subroutine fun













