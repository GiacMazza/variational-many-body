MODULE GZ_ENERGY_MINIMIZATION 
  ! scifor routines
  USE SF_LINALG
  USE SF_IOTOOLS
  USE SF_SPECIAL
  USE SF_OPTIMIZE
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
  public :: get_gz_ground_state_estimation !+- get the GZground state estimation at fixed n0
  !public :: get_gz_ground_state_estimation_vdm_R !+- get the GZground state estimation at fixed n0 and Rhop (TO BE CODED) ...and overload the two..
  !
  public :: gz_energy_vdm
  !
contains
  !
  include "slater_min_routines.f90"
  include "GZproj_min_routines.f90"
  !
  subroutine get_gz_ground_state_estimation(optimized_vdm)
    real(8),dimension(Ns),intent(in)    :: optimized_vdm
    real(8)                                    :: energy
    integer                                    :: iorb,jorb,istate,jstate
    integer :: ifock,iphi,jphi,ifock_
    !
    real(8),dimension(Nphi) :: phi_vector_test
    complex(8),dimension(Nphi) :: phi_vec
    !
    optimization_flag=.true.
    allocate(GZ_opt_projector_diag(Nphi))
    allocate(GZ_opt_Rhop(Ns,Ns))
    !

    !+- IN THE PLACE OF THIS TWO LINES I'LL PUT THE FUNCTIONS 
    !   get_optimized_energy_as_a_functions_of_free_parameters 
    !   (to be used also in the simplex minimization)
    !   
    select case(min_method)
    case('nlep')
       energy=gz_energy_recursive_nlep(optimized_vdm)
    case('cmin')
       energy=gz_energy_recursive_cmin(optimized_vdm)
    end select
    !
    GZ_opt_Rhop=hopping_renormalization_normal(GZ_opt_projector_diag,optimized_vdm)
    !    
    phi_vec=GZ_opt_projector_diag
    !
    !+- GET OBSERVABLES -+!
    ! physical density !
    allocate(gz_dens(Ns))
    do istate=1,Ns
       gz_dens(istate) = trace_phi_basis(phi_vec,phi_traces_basis_local_dens(istate,istate,:,:))
    end do
    ! density-density same orbital -aka orbital doubly occupancy-!
    allocate(gz_docc(Norb))
    do iorb=1,Norb
       gz_docc(iorb) = trace_phi_basis(phi_vec,phi_traces_basis_docc_orb(iorb,:,:))
    end do
    ! density-density different orbitals !
    allocate(gz_dens_dens_orb(Norb,Norb))
    do iorb=1,Norb
       do jorb=1,Norb
          gz_dens_dens_orb(iorb,jorb)=trace_phi_basis(phi_vec,phi_traces_basis_dens_dens_orb(iorb,jorb,:,:))
       end do
    end do
    !+-
    ! place for other observables... SPINS,ISO-SPINS,...bla bla bla
    !+-
  end subroutine get_gz_ground_state_estimation



  !
  function gz_energy_vdm(vdm) result(GZenergy)
    real(8),dimension(:),intent(in) :: vdm
    real(8)               :: GZenergy
    select case(min_method)
    case('nlep')
       GZenergy=gz_energy_recursive_nlep(vdm)
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
    complex(8),dimension(Ns,Ns) :: Rhop
    complex(8),dimension(Ns,Ns) :: slater_derivatives    
    real(8),dimension(Ns)           :: slater_lgr_multip,R_diag
    real(8),dimension(Ns,Ns,2) :: GZproj_lgr_multip  ! 
    real(8),dimension(Ns,Ns) :: GZproj_lgr_multip_  ! 
    real(8)                                :: E_Hstar,E_Hloc
    complex(8),dimension(nPhi)               :: GZvect_iter  ! GZ vector (during iterations)      
    integer :: is,imap
    !
    do is=1,Ns
       imap = vdm_map(is)
       vdm(is) = x(imap)
    end do
    !
    GZ_energy = gz_energy_vdm(vdm)
    !      
  end subroutine gz_get_energy_vdm
  
  function gz_energy_recursive_nlep(n0)   result(GZ_energy)
    real(8),dimension(:),intent(in) :: n0 !INPUT: Variational Density Matrix (VDM) (diagonal in istate)    
    real(8)                         :: GZ_energy !INPUT: Optimized GZ energy at fixed 
    real(8)                         :: GZ_energy_old,energy_err     ! Value of the GZ energy functional
    complex(8),dimension(Ns,Ns)     :: R_init        ! initial guess for the hopping renormalization matrix    
    complex(8),dimension(Ns,Ns)     :: slater_derivatives    
    complex(8),dimension(Ns,Ns)     :: R_iter,R_old ! hopping matrix renormalization (during iterations)
    real(8),dimension(Ns)           :: R_diag
    real(8),dimension(Ns,Ns)        :: GZproj_lgr_multip,slater_lgr_multip  ! 
    real(8)                         :: E_Hstar,E_Hloc
    complex(8),dimension(nPhi)      :: GZvect_iter  ! GZ vector (during iterations)
    !
    integer                         :: istate,iter,jstate,ifock,jfock,iphi,jphi,is
    integer                         :: unit
    logical                         :: bound
    !
    write(*,*) '********************'
    write(*,*) 'INPUT DENSITY',n0(:)
    bound=.false.
    do istate=1,Ns
       if(n0(istate).lt.1.d-11.or.n0(istate).gt.1.d0-1.d-11) bound=.true.
    end do
    !
    if(.not.bound) then
       !+- initialize Rhop according to a given wanted symmetry
       R_init=0.d0 
       do istate=1,Ns
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
          select case(lgr_method)
          case('amoeba')
             call gz_projectors_minimization_nlep(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,GZmin_verbose)   
          case('fsolve')
             call gz_projectors_minimization_nlep_fsolve(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,GZmin_verbose)   
          end select
          !
          R_iter=hopping_renormalization_normal(GZvect_iter,n0)
          R_iter=Rmix*R_iter+(1.d0-Rmix)*R_old
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
       write(opt_rhop_unit,*) n0,R_diag(1:Ns)
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
    complex(8),dimension(Nphi)                  :: GZvect_iter  ! GZ vector (during iterations)

    complex(8),dimension(Ns,Ns)    :: R_iter ! hopping matrix renormalization (during iterations)

    complex(8),dimension(Ns,Ns)    :: R_init        ! initial guess for the hopping renormalization matrix    
    real(8),dimension(Ns)              :: R_diag
    complex(8),dimension(Ns,Ns,Lk) :: slater_matrix_el    

    real(8),dimension(Ns,Ns)              :: slater_lgr_multip
    real(8),dimension(Ns,Ns)    :: GZproj_lgr_multip  
    real(8)                                   :: E_Hstar,E_Hloc
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
       if(n0(istate).le.1.d-10.or.n0(istate).ge.1.d0-1.d-10) bound=.true.
    end do
    !

    if(.not.bound) then
       !+- get not-interacting GZprojectors corresponding to this density matrix -+!
       call initialize_GZprojectors(GZvect_iter,n0)
       !
       R_iter=hopping_renormalization_normal(GZvect_iter,n0)
       !
       GZ_energy=0.d0
       do iter=1,Niter_self
          !+- update phi_vectors -+!
          GZ_energy_old=GZ_energy
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
          call gz_projectors_minimization_cmin(slater_matrix_el,n0,GZvect_iter,GZ_energy,GZproj_lgr_multip,GZmin_verbose)
          !
          R_iter=hopping_renormalization_normal(GZvect_iter,n0)
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
       !
    else
       !
       GZ_energy=100.d0
       !
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
    call slater_determinant_minimization_nlep(R_init,n0,E_Hstar,slater_lgr_multip,slater_derivatives,iverbose)
    call free_gz_projectors_init(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,iverbose)    
  end subroutine initialize_GZprojectors
  !
  !+- AS IT IS THIS ROUTINE IS WRONG, INDEED IN ORDER TO IMPLEMENT THE FIX-R MINIMIZATION NO SLATER DERIVATIVES HAVE TO BE CONSIDERED....I'LL COME BACK LATER ON THIS... -+!
  !
  subroutine gz_projectors_minimization_fixR(slater_derivatives,n0_target,R_target,E_Hloc,GZvect,lgr_multip,iverbose)
    complex(8),dimension(Ns,Ns),intent(in)  :: slater_derivatives !input:  Slater Deter GZ energy derivatives
    real(8),dimension(Ns),intent(in)     :: n0_target          !input:  Variational density matrix
    complex(8),dimension(Ns,Ns),intent(in)  :: R_target          !input:  Renormalization matrices
    real(8),dimension(2*Ns)                :: lgr
    real(8),dimension(Ns) :: lgr_tmp
    real(8),dimension(1) :: lgr_tmp_
    real(8),dimension(Ns,Ns,2)                :: lgr_multip
    complex(8),dimension(nPhi)              :: GZvect   !output: GZvector
    real(8)                              :: E_Hloc,Emin   !output: optimized local energy
    real(8)                              :: lgr_symm(1)
    logical,optional                     :: iverbose
    real(8),dimension(2*Ns)                :: err_dens
    real(8),dimension(Ns)                :: err_dens_tmp
    real(8),dimension(Ns*Ns)             :: err_dens_full
    real(8)                              :: off_constraints
    real(8)                              :: delta_out
    logical                              :: iverbose_
    integer                              :: info,istate,i,j,iter
    !+- amoeba_variables-+!
    real(8),allocatable,dimension(:,:)   :: p
    real(8),allocatable,dimension(:)     :: y
    real(8)                              :: ftol
    integer                              :: np,mp,i_vertex,j_vertex,i_dim,is

    iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose    
    !
    lgr=0.d0; lgr_tmp=0.d0
    do istate=1,Ns
       lgr_tmp(istate)=(0.5d0-n0_target(istate))*2.d0
       lgr(istate)=(0.5d0-n0_target(istate))*2.d0
       !lgr(istate+Ns)=0.2d0
    end do

    !
    lgr_tmp=0.d0
    lgr_tmp_=0.d0
    ! call fsolve(constraints_deviation_tmp,lgr_tmp_,tol=1.d-15,info=info)    
    ! stop
    !call fsolve(get_delta_proj_variational_density_diag,lgr_tmp,tol=1.d-15,info=info)    


    !Simplex minimization minimizing of the deviation from the variational target density no_target
    NP=2*NS
    MP=NP+1
    allocate(y(MP),p(MP,NP))
    !+- initialize simplex -+!
    p(1,:) = 0.d0
    do i=1,2*NS
       do j=1,Ns
          p(i+1,j) = (0.5d0-n0_target(j))!*dble(i)*0.1
          p(i+1,j+Ns) = dble(i)*dble(j)*0.5
       end do
    end do
    !
    do i_vertex=1,MP
       y(i_vertex)=constraints_deviation_amb(p(i_vertex,:))
    end do
    call amoeba_er(p(1:MP,1:NP),y(1:MP),1.d-8,constraints_deviation_amb,iter,verbose=.false.)
    lgr=p(1,:); deallocate(y,p)

    !call fsolve(constraints_deviation,lgr,tol=1.d-15,info=info)    

    !    
    lgr_multip=0.d0
    do istate=1,Ns
       lgr_multip(istate,istate,1)=lgr(istate)
       lgr_multip(istate,istate,2)=0.d0!lgr(istate+Ns)
    end do
    write(*,*) lgr_multip(:,:,1)
    call get_GZproj_ground_state_fixR(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect)
    !
    if(iverbose_) then
       write(*,*)
       write(*,*) "GZ projectors: Lagrange Parameters -diagonal case-"
       write(*,'(10F18.10)') lgr
       !
       err_dens=constraints_deviation_amb(lgr)
       !err_dens_tmp=get_delta_proj_variational_density_diag(lgr_tmp)
       write(*,*) "GZ projectors: Variational density matrix error"
       write(*,'(10F18.10)') err_dens
       !
       write(*,*) "GZ projectors: Optimized Local Energy"
       write(*,'(10F18.10)') E_Hloc
       write(*,*)
    end if
    !
  contains
    !
    !include 'self_minimization_GZproj_routines.f90'
    include 'fixR_min_routines.f90'
    !
  end subroutine gz_projectors_minimization_fixR









  ! gz_projectors (lancelot subroutine)
  subroutine gz_projectors_minimization_cmin(slater_matrix_el,n0,GZvect_indep,GZenergy,GZproj_lgr_multip,GZmin_verbose)
    complex(8),dimension(Ns,Ns,Lk),intent(in) :: slater_matrix_el
    real(8),dimension(Ns),intent(in)       :: n0
    complex(8),dimension(Nphi),intent(inout)  :: GZvect_indep    
    real(8),dimension(Nphi)   :: GZvect_indep_
    real(8),intent(inout)                  :: GZenergy
    real(8),dimension(Ns,Ns),intent(inout) :: GZproj_lgr_multip
    logical                                :: GZmin_verbose
    real(8),dimension(Ns)                  :: vdm
    integer                                :: iter,iter_max,istate,i,iter_
    integer                                :: iorb,ispin,i_ind,ifock
    !
    !LANCELOT VARIABLES
    !
    integer                                :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                                :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_
    real(8),allocatable                    :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                                :: iunit,err_unit,ene_unit,Nsuccess    
    !
    if(allocated(slater_matrix_elements)) deallocate(slater_matrix_elements)
    allocate(slater_matrix_elements(Ns,Ns,Lk))
    slater_matrix_elements=slater_matrix_el
    vdm=n0
    !
    ! LANCELOT configuration parameters 
    n_min       = Nphi ! number of minimization parameters
    neq         = Norb + 1    ! number of equality constraints                   
    nin         = 0           ! number of in-equality constraints                   
    maxit       = 1000        ! maximum iteration number 
    gradtol     = 1.d-7       ! maximum norm of the gradient at convergence 
    feastol     = 1.d-7       ! maximum violation of parameters at convergence  
    print_level = lancelot_verbose           ! verbosity
    allocate(bl(n_min),bu(n_min),cx(neq+nin),y(neq+nin))
    bL = 0.d0               ! lower bounds for minimization parameters
    bU = 1.d0-bL                 ! upper bounds for minimization parameters          
    if(GZmin_verbose) then
       write(GZmin_unit_,*) '!+---------------------------------------------------+!'
       write(GZmin_unit_,*) '!+---------------------------------------------------+!'
       write(GZmin_unit_,*) 'Number of minimization parameters', n_min
       write(GZmin_unit_,*) 'Number of equality constrains', neq
       write(GZmin_unit_,*) 'Number of inquality constrains', nin
       write(GZmin_unit_,*) 'Maximum number of iterations', maxit
       write(GZmin_unit_,*) 'Gradient tolerance', gradtol
       write(GZmin_unit_,*) 'Constraints tolerance', feastol
       write(GZmin_unit_,*) 
       write(GZmin_unit_,*) 
       write(GZmin_unit_,*) 'INPUT PARAMETERS'
       do i=1,n_min
          write(GZmin_unit_,*) GZvect_indep(i)
       end do
       write(GZmin_unit_,*) 
       write(GZmin_unit_,*)           
    end if
    GZvect_indep_=GZvect_indep
    call lancelot_simple(n_min,GZvect_indep_,GZenergy,exit_code,my_fun=energy_GZproj_functional, &
         bl = bl, bu = bu,                                                                      &
         neq = neq, nin = nin,                                                                  &
         cx = cx, y = y, iters  = iter, maxit = maxit,                                          &
         gradtol = gradtol, feastol = feastol,                                                  &
         print_level = print_level )
    !+--------------------------------------------------------------------------------------+!    
    GZproj_lgr_multip=0.d0
    do iorb=1,Norb
       do ispin=1,2
          istate=index(ispin,iorb)
          GZproj_lgr_multip(istate,istate)=cx(iorb)
       end do
    end do
    !
    if(GZmin_verbose) then
       write(GZmin_unit_,*) 'LANCELOT EXIT STATUS'
       write(GZmin_unit_,*)
       write(GZmin_unit_,*) exit_code
       write(GZmin_unit_,*) 'OPTIMIZED PARAMETERS'
       do i=1,n_min
          write(GZmin_unit_,*) GZvect_indep_(i)
       end do
       write(GZmin_unit_,*) 
       write(GZmin_unit_,*) 
       write(GZmin_unit_,*) 'FINAL VALUE'
       write(GZmin_unit_,*) GZenergy
       write(GZmin_unit_,*) 
       write(GZmin_unit_,*) 
       write(GZmin_unit_,*)  'CONSTRAINTS'
       do i=1,neq+nin
          write(GZmin_unit_,*) cx(i)
       end do
       write(GZmin_unit_,*) 
       write(GZmin_unit_,*) 
       write(GZmin_unit_,*)  'LAGRANGE MULTIPLIERS'
       do i=1,neq+nin
          write(GZmin_unit_,*) y(i)
       end do
       write(GZmin_unit_,*) 
       write(GZmin_unit_,*) 
       write(GZmin_unit_,*)  'NUMBER OF ITERATION'
       write(GZmin_unit_,*) iter,'/',maxit       
       write(GZmin_unit_,*) 'Eiteration',GZenergy
    end if
    deallocate(bl,bu,cx,y)
    deallocate(slater_matrix_elements)
    GZvect_indep=GZvect_indep_

  contains

    subroutine energy_GZproj_functional(x,f,i)
      implicit none
      !+- routine variables -+!
      real(8), intent(in)           :: x(:)
      real(8), intent(out)          :: f
      integer, intent(in), optional :: i
      complex(8),allocatable           :: phi_(:)
      real(8),allocatable           :: niorb(:)
      real(8)                       :: Estar
      real(8),dimension(Ns)  :: local_dens_
      real(8),dimension(Ns,Ns)  :: Rhop    
      real(8),dimension(Ns,Ns)  :: Hk
      !
      integer                       :: iorb,ispin,istate,jstate,ik,ifock,jfock,jorb,jspin,iphi,jphi
      !
      allocate(phi_(Nphi))
      phi_=x
      allocate(niorb(Norb))
      !
      do iorb=1,Norb
         niorb(iorb) =  0.d0
         do ispin=1,2
            istate=index(ispin,iorb)
            niorb(iorb) = niorb(iorb) + vdm(istate)
         end do
      end do
      !
      do istate=1,Ns
         do jstate=1,Ns
            Rhop(istate,jstate) = 0.d0
            do iphi=1,Nphi
               do jphi=1,Nphi
                  Rhop(istate,jstate) = &
                       Rhop(istate,jstate) + 1.d0/sqrt(vdm(jstate)*(1-vdm(jstate)))*phi_(iphi)*phi_(jphi)*phi_traces_basis_Rhop(istate,jstate,iphi,jphi)
               end do
            end do
         end do
      end do
      !
      if (.not.present(i)) then
         !+- FREE ENERGY ESTIMATION WITHIN GZ APPROX -+!
         Estar=0.d0
         do ik=1,Lk
            Hk=0.d0
            ! hopping renormalization !
            Hk=matmul(Hk_tb(:,:,ik),Rhop)
            Hk=matmul(Rhop,Hk)
            do istate=1,Ns
               do jstate=1,Ns
                  Estar = Estar + slater_matrix_elements(istate,jstate,ik)*wtk(ik)*Hk(istate,jstate)
               end do
            end do
         end do
         !
         f=Estar
         !
         f=f+trace_phi_basis(phi_,phi_traces_basis_Hloc)
         !
      else
         !+- CONSTRAINTS ON GUTZWILLER PARAMETERS -+!
         select case(Norb)
         case(1)
            select case(i)
            case(1)
               f=0.d0
               iorb=1
               do ispin=1,2
                  istate=index(ispin,iorb)
                  !
                  f = f + trace_phi_basis(phi_,phi_traces_basis_dens(istate,istate,:,:))
                  !
               end do
               f=f-niorb(iorb)
            case(2)
               f=0.d0
               !
               do iphi=1,Nphi
                  f = f + phi_(iphi)*phi_(iphi)
               end do
               !
               f=f-1.d0
            end select
         case(2)
            select case(i)
            case(1)
               f=0.d0
               iorb=1
               do ispin=1,2
                  istate=index(ispin,iorb)
                  !
                  f = f + trace_phi_basis(phi_,phi_traces_basis_dens(istate,istate,:,:))
                  !
               end do
               f=f-niorb(iorb)
            case(2)
               f=0.d0
               iorb=2
               do ispin=1,2
                  istate=index(ispin,iorb)
                  !
                  f = f + trace_phi_basis(phi_,phi_traces_basis_dens(istate,istate,:,:))
                  !
               end do
               f=f-niorb(iorb)
            case(3)
               f=0.d0
               !
               do iphi=1,Nphi
                  f = f + phi_(iphi)*phi_(iphi)
               end do
               !
               f=f-1.d0
            end select
         end select
      end if
    end subroutine energy_GZproj_functional


  end subroutine gz_projectors_minimization_cmin
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













