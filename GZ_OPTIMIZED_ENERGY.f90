MODULE GZ_OPTIMIZED_ENERGY
  USE SCIFOR
  !USE SF_OPTIMIZE
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_PROJECTORS
  USE GZ_ENERGY_MINIMIZATION
  USE MIN_AMOEBA
  implicit none
  private

  public :: gz_optimization_simplex
  public :: get_gz_ground_state_estimation
  logical,public :: optimization_flag

CONTAINS

  !+-----------------------------------------------------------------------------------------------------+!
  !+- PURPOSE: Minimize the GUTZWILLER ENERGY FUNCTIONAL WITH RESPECT TO THE VARIATIONAL DENSITY MATRIX -+!
  !+-----------------------------------------------------------------------------------------------------+!
  subroutine get_gz_ground_state_estimation(optimized_vdm)
    real(8),dimension(state_dim),intent(in)    :: optimized_vdm
    real(8)                                    :: energy
    real(8),dimension(nFock,nFock)             :: local_operator
    integer                                    :: iorb,jorb,istate,jstate
    !
    optimization_flag=.true.
    allocate(GZ_opt_projector_diag(nFock))
    allocate(GZ_opt_Rhop(state_dim,state_dim))
    !allocate()
    select case(min_method)
    case('nlep')
       energy=gz_energy_recursive_nlep(optimized_vdm)
    case('cmin')
       energy=gz_energy_recursive_cmin(optimized_vdm)
    end select
    !

    !+- get renormalization matrices -+!

    GZ_opt_Rhop=Rhop_matrix(GZ_opt_projector_diag,optimized_vdm)

    !+- GET OBSERVABLES -+!
    ! physical density !
    allocate(gz_dens(state_dim))
    do istate=1,state_dim
       gz_dens(istate)=gz_local_diag(GZ_opt_projector_diag,dens(istate,:,:))
    end do

    ! density-density same orbital -aka orbital doubly occupancy-!
    allocate(gz_docc(Norb))
    do iorb=1,Norb
       gz_docc(iorb)=gz_local_diag(GZ_opt_projector_diag,docc(iorb,:,:))
    end do

    ! density-density different orbitals !
    allocate(gz_dens_dens_orb(Norb,Norb))
    do iorb=1,Norb
       do jorb=1,Norb
          gz_dens_dens_orb(iorb,jorb)=&
               gz_local_diag(GZ_opt_projector_diag,dens_dens_orb(iorb,jorb,:,:))
       end do
    end do
    !+-
    ! place for other observables... SPINS,ISO-SPINS,...bla bla bla
    !+-
  end subroutine get_gz_ground_state_estimation

  subroutine gz_optimization_simplex(simplex_init,optimized_vdm) 
    real(8),dimension(state_dim+1,state_dim),intent(inout) :: simplex_init
    real(8),dimension(state_dim),intent(out)               :: optimized_vdm
    !real(8),intent(out)                                    :: optimized_energy    
    !+- amoeba_variables-+!
    real(8),allocatable,dimension(:,:)                     :: p
    real(8),allocatable,dimension(:)                       :: y
    real(8)                                                :: ftol
    integer                                                :: np,mp,i_vertex,j_vertex,i_dim,iter
    integer,allocatable,dimension(:)                       :: idum
    integer                                                :: tmp_dum
    real(8)                                                :: rnd,tot_dens,tmp_pol
    integer                                                :: amoeba_unit    
    !+-------------------+!
    optimization_flag=.false.
    amoeba_unit=free_unit()
    open(amoeba_unit,file='Amoeba_minimization.out')
    opt_energy_unit=free_unit()
    open(opt_energy_unit,file='GZ_OptEnergy_VS_vdm.out')
    opt_rhop_unit=free_unit()
    open(opt_rhop_unit,file='GZ_OptRhop_VS_vdm.out')
    opt_GZ_unit=free_unit()
    open(opt_GZ_unit,file='GZ_OptProj_VS_vdm.out')
    if(GZmin_verbose) then
       GZmin_unit=free_unit()
       open(GZmin_unit,file='GZ_SelfCons_min_verbose.out')
       GZmin_unit_=free_unit()
       open(GZmin_unit_,file='GZ_proj_min.out')
    end if
    NP=state_dim
    MP=NP+1  
    allocate(y(MP),p(MP,NP))
    !+- initialize simplex -+!
    p=simplex_init
    write(amoeba_unit,*) 'Initializing simplex verteces'
    !+----------------------+!  
    do i_vertex=1,MP
       select case(min_method)
       case('nlep')
          y(i_vertex)=gz_energy_recursive_nlep(p(i_vertex,:))
       case('cmin')
          y(i_vertex)=gz_energy_recursive_cmin(p(i_vertex,:))
       end select
       write(amoeba_unit,*) p(i_vertex,:),y(i_vertex)
    end do
    !stop
    ftol=amoeba_min_tol
    select case(min_method)
    case('nlep')
       call amoeba(p(1:MP,1:NP),y(1:MP),ftol,gz_energy_recursive_nlep,iter,amoeba_verbose)
    case('cmin')
       call amoeba(p(1:MP,1:NP),y(1:MP),ftol,gz_energy_recursive_cmin,iter,amoeba_verbose)
    end select
    !+- Optimized Variational Density Matrix -+!
    optimized_vdm=p(1,:) 
    !+- Optimized Variational Energy -+!
    simplex_init=p
    write(amoeba_unit,*) 
    write(amoeba_unit,*) 'Last loop simplex verteces'
    write(amoeba_unit,*) 
    do i_vertex=1,MP
       write(amoeba_unit,*) p(i_vertex,:),y(i_vertex)
    end do
    write(amoeba_unit,*) 
    write(amoeba_unit,*) 
    deallocate(y,p)
    close(amoeba_unit)
    close(opt_energy_unit)
  end subroutine gz_optimization_simplex



  function gz_energy_recursive_nlep(n0)   result(GZ_energy)
    real(8),dimension(:),intent(in)        :: n0 !INPUT: Variational Density Matrix (VDM) (diagonal in istate)    
    real(8)                                :: GZ_energy !INPUT: Optimized GZ energy at fixed 
    real(8)                                :: GZ_energy_old,energy_err     ! Value of the GZ energy functional
    real(8),dimension(state_dim,state_dim) :: R_init        ! initial guess for the hopping renormalization matrix    
    real(8),dimension(state_dim,state_dim) :: slater_derivatives    
    real(8),dimension(state_dim,state_dim) :: R_iter,R_old ! hopping matrix renormalization (during iterations)
    real(8),dimension(state_dim)           :: slater_lgr_multip,R_diag
    real(8),dimension(state_dim,state_dim) :: GZproj_lgr_multip  ! GZ vector (during iterations)
    real(8)                                :: E_Hstar,E_Hloc
    real(8),dimension(nFock)               :: GZvect_iter  ! GZ vector (during iterations)
    !
    integer                                :: istate,iter,jstate,ifock,jfock
    integer                                :: unit
    logical                                :: bound
    !
    write(*,*) '********************'
    write(*,*) 'INPUT DENSITY',n0(:)
    bound=.false.
    do istate=1,state_dim
       if(n0(istate).le.1.d-10.or.n0(istate).ge.1.d0-1.d-10) bound=.true.
    end do
    !
    if(.not.bound) then
       R_init=0.d0 
       do istate=1,state_dim
          R_init(istate,istate)=Rseed
       end do
       !
       GZ_energy=0.d0    
       R_iter=R_init
       GZvect_iter=0.d0
       do iter=1,Niter_self
          !+- update hopping matrices -+!
          GZ_energy_old=GZ_energy
          R_old=R_iter          
          !+----------------------------+!
          !+- SLATER STEP MINIMIZATION -+!
          !+----------------------------+!    
          call slater_determinant_minimization_nlep(R_iter,n0,E_Hstar,slater_lgr_multip,slater_derivatives,GZmin_verbose)       
          !+----------------------------+!
          !+- GZproj STEP MINIMIZATION -+!
          !+----------------------------+!    
          call gz_projectors_minimization_nlep(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,GZmin_verbose)          
          ! update R-matrix
          R_iter=Rmix*Rhop_matrix(GZvect_iter,n0)+(1.d0-Rmix)*R_old
          do istate=1,state_dim
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
             write(GZmin_unit,'(20F18.10)') dble(iter),energy_err,GZ_energy,E_Hstar,E_Hloc,R_diag(1:state_dim)
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
       do ifock=1,nFock
          write(opt_GZ_unit,*) GZvect_iter(ifock)
       end do
       write(opt_GZ_unit,*)
       write(opt_GZ_unit,*)
       write(opt_energy_unit,*) n0,GZ_energy,E_Hloc,E_Hstar
       write(opt_rhop_unit,*) n0,R_diag(1:state_dim)
    else
       GZ_energy=100.d0
    end if
    if(optimization_flag) then
       !+- store final informations to global variables -+!              
       GZ_opt_projector_diag = GZvect_iter
       GZ_opt_energy         = GZ_energy
       GZ_opt_kinetic        = E_Hstar
       GZ_opt_Eloc           = E_Hloc
    end if
    !
  end function gz_energy_recursive_nlep







  function gz_energy_recursive_cmin(n0)  result(GZ_energy)
    real(8),dimension(:),intent(in)           :: n0 !INPUT: Variational Density Matrix (VDM) (diagonal in istate)    
    real(8)                                   :: GZ_energy !INPUT: Optimized GZ energy at fixed 
    real(8)                                   :: GZ_energy_old,energy_err     ! Value of the GZ energy functional
    real(8),dimension(nFock)                  :: GZvect_iter  ! GZ vector (during iterations)
    real(8),dimension(nFock_indep)            :: GZvect_iter_,GZvect_iter_tmp  ! GZ vector (during iterations)
    real(8),dimension(state_dim,state_dim)    :: R_iter ! hopping matrix renormalization (during iterations)

    real(8),dimension(state_dim,state_dim)    :: R_init        ! initial guess for the hopping renormalization matrix    
    real(8),dimension(state_dim)              :: R_diag
    real(8),dimension(state_dim,state_dim,Lk) :: slater_matrix_el    

    real(8),dimension(state_dim)              :: slater_lgr_multip
    real(8),dimension(state_dim,state_dim)    :: GZproj_lgr_multip  
    real(8)                                   :: E_Hstar,E_Hloc
    !
    integer                                   :: istate,iter,jstate,ifock,jfock,i_ind
    integer                                   :: unit
    logical                                   :: bound
    !
    write(*,*) '*************************'
    write(*,*) 'INPUT DENSITY',n0(:)
    bound=.false.
    do istate=1,state_dim
       if(n0(istate).le.1.d-10.or.n0(istate).ge.1.d0-1.d-10) bound=.true.
    end do
    !

    if(.not.bound) then
       !+- get not-interacting GZprojectors corresponding to this density matrix -+!
       call initialize_GZprojectors(GZvect_iter,n0)
       !
       do i_ind=1,nFock_indep
          GZvect_iter_(i_ind) = GZvect_iter(fock_indep(i_ind))
       end do
       R_iter=Rhop_matrix(GZvect_iter,n0)
       !
       GZ_energy=0.d0
       do iter=1,Niter_self
          !+- update phi_vectors -+!
          GZ_energy_old=GZ_energy
          GZvect_iter_tmp=GZvect_iter_
          !
          !+----------------------------+!
          !+- SLATER STEP MINIMIZATION -+!
          !+----------------------------+!    
          !
          call slater_determinant_minimization_cmin(R_iter,n0,E_Hstar,slater_lgr_multip,slater_matrix_el,GZmin_verbose)       
          !
          !+----------------------------+!
          !+- GZproj STEP MINIMIZATION -+!
          !+----------------------------+!    
          !
          call gz_projectors_minimization_cmin(slater_matrix_el,n0,GZvect_iter_,GZ_energy,GZproj_lgr_multip,GZmin_verbose)
          !
          R_iter=Rhop_matrix(GZvect_iter,n0)          
          do istate=1,state_dim
             R_diag(istate)=R_iter(istate,istate)
          end do
          !
          do ifock=1,nFock
             GZvect_iter(ifock)=GZvect_iter_(full2indep_fock(ifock))     
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
             write(GZmin_unit,*) dble(iter),energy_err,GZ_energy,E_Hstar,E_Hloc,R_diag(1:state_dim)
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
       do ifock=1,nFock
          write(opt_GZ_unit,*) GZvect_iter(ifock)
       end do
       write(opt_GZ_unit,*)
       write(opt_GZ_unit,*)
       write(opt_energy_unit,*) n0,GZ_energy,E_Hstar,E_Hloc
       R_iter=Rhop_matrix(GZvect_iter,n0)          
       do istate=1,state_dim
          R_diag(istate)=R_iter(istate,istate)
       end do
       write(opt_rhop_unit,*) n0,R_diag(1:state_dim)
    else
       GZ_energy=100.d0
    end if
    if(optimization_flag) then
       !+- store final informations to global variables -+!              
       GZ_opt_projector_diag = GZvect_iter
       GZ_opt_energy         = GZ_energy
       GZ_opt_kinetic        = E_Hstar
       GZ_opt_Eloc           = E_Hloc
    end if
    !
  end function gz_energy_recursive_cmin
  !



  subroutine initialize_GZprojectors(GZvect_iter,n0)
    real(8),dimension(nFock) :: GZvect_iter
    real(8),dimension(state_dim) :: n0
    real(8),dimension(state_dim,state_dim) :: R_init        
    real(8),dimension(state_dim,state_dim) :: slater_derivatives
    real(8),dimension(state_dim)           :: slater_lgr_multip
    real(8),dimension(state_dim,state_dim) :: GZproj_lgr_multip  
    real(8) :: E_Hstar,E_HLoc
    logical :: iverbose
    integer :: istate

    iverbose=.false.
    R_init=0.d0
    do istate=1,state_dim
       R_init(istate,istate)=1.d0
    end do
    call slater_determinant_minimization_nlep(R_init,n0,E_Hstar,slater_lgr_multip,slater_derivatives,iverbose)
    call free_gz_projectors_init(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,iverbose)
    
  end subroutine initialize_GZprojectors

END MODULE GZ_OPTIMIZED_ENERGY
  
