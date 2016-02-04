MODULE GZ_OPTIMIZED_ENERGY
  USE SCIFOR
  !
  !USE SF_OPTIMIZE
  !
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_EFFECTIVE_HOPPINGS
  USE GZ_LOCAL_FOCK  
  !
  USE LANCELOT_simple_double
  !
  !  USE GZ_PROJECTORS
  !
  USE GZ_MATRIX_BASIS
  USE GZ_ENERGY_MINIMIZATION
  USE MIN_AMOEBA
  implicit none
  private


  public :: get_gz_ground_state_estimation
  logical,public :: optimization_flag
  !
  public :: gz_optimization_simplex
  public :: gz_optimization_vdm_Rhop
  public :: gz_optimization_vdm


  public :: gz_optimization_vdm_nlsq
  !
  public :: gz_energy_broyden
  public :: gz_energy_recursive_nlep
  public :: gz_energy_recursive_cmin
  !
CONTAINS

  !+-----------------------------------------------------------------------------------------------------+!
  !+- PURPOSE: Minimize the GUTZWILLER ENERGY FUNCTIONAL WITH RESPECT TO THE VARIATIONAL DENSITY MATRIX -+!
  !+-----------------------------------------------------------------------------------------------------+!
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
    select case(min_method)
    case('nlep')
       energy=gz_energy_recursive_nlep(optimized_vdm)
    case('cmin')
       energy=gz_energy_recursive_cmin(optimized_vdm)
    case('bryd')
       energy=gz_energy_broyden(optimized_vdm)
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




  !+- TO FIX the fixR routine!!!

  subroutine gz_optimization_vdm_Rhop(init_vdm,init_Rhop,opt_vdm,opt_Rhop) 
    real(8),dimension(Ns),intent(inout)    :: init_vdm
    complex(8),dimension(Ns),intent(inout) :: init_Rhop
    real(8),dimension(Ns),intent(out)      :: opt_vdm
    complex(8),dimension(Ns),intent(out)   :: opt_Rhop

    real(8),dimension(3*Ns) :: xmin
    real(8) :: Emin
    integer :: is,iter
    !
    do is=1,Ns
       xmin(is) = init_vdm(is)
       xmin(is+Ns) = dreal(init_Rhop(is))
       xmin(is+2*Ns) = dimag(init_Rhop(is))
    end do
    !    
    write(*,*) "XMIN"
    do is=1,3*Ns
       write(*,*) xmin(is)
    end do
    Emin=gz_energy_vdm_Rhop(xmin)
    write(*,*) Emin
    !call fmin_cg(xmin,gz_energy_vdm_Rhop,iter,Emin)
    !
  contains
    !
    function gz_energy_vdm_Rhop(x) result(GZ_energy)
      real(8),dimension(:) :: x
      real(8)              :: GZ_energy
      real(8),dimension(Ns) :: vdm
      complex(8),dimension(Ns,Ns) :: Rhop
      complex(8),dimension(Ns,Ns) :: slater_derivatives    
      real(8),dimension(Ns)           :: slater_lgr_multip,R_diag
      real(8),dimension(Ns,Ns,2) :: GZproj_lgr_multip  ! 
      real(8),dimension(Ns,Ns) :: GZproj_lgr_multip_  ! 
      real(8)                                :: E_Hstar,E_Hloc
      complex(8),dimension(nPhi)               :: GZvect_iter  ! GZ vector (during iterations)      
      integer :: is
      !
      if(size(x).ne.3*Ns) stop "gz_energy_vdm_Rhop/ wrong dimensions"
      !
      Rhop=zero
      do is=1,Ns
         vdm(is)=x(is)
         Rhop(is,is)=x(is+Ns) + xi*x(is+2*Ns)
         write(*,*) Rhop(is,is)
      end do
      !
      call slater_determinant_minimization_nlep(Rhop,vdm,E_Hstar,slater_lgr_multip,slater_derivatives,GZmin_verbose)       
      !
      call gz_projectors_minimization_fixR(slater_derivatives,vdm,Rhop,E_Hloc,GZvect_iter,GZproj_lgr_multip,GZmin_verbose)                   
      !
      GZ_energy=E_Hstar+E_Hloc
      !
    end function gz_energy_vdm_Rhop

  end subroutine gz_optimization_vdm_Rhop






  subroutine gz_optimization_vdm(init_vdm,optimized_vdm)
    real(8),dimension(1),intent(in)                :: init_vdm
    real(8),dimension(1),intent(out)               :: optimized_vdm
    integer :: iter,unit_vdm_opt,icall
    real(8) :: GZ_energy
    !
    optimized_vdm=init_vdm

    unit_vdm_opt=free_unit()
    open(unit_vdm_opt,file='vdm_optimization.out'); icall=0
    !
    call fmin_cg(optimized_vdm,gz_energy_vdm,iter,GZ_energy)


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


    !
  contains


    function gz_energy_vdm(x) result(GZ_energy)
      real(8),dimension(:) :: x
      real(8)              :: GZ_energy
      real(8),dimension(Ns) :: vdm
      complex(8),dimension(Ns,Ns) :: Rhop
      complex(8),dimension(Ns,Ns) :: slater_derivatives    
      real(8),dimension(Ns)           :: slater_lgr_multip,R_diag
      real(8),dimension(Ns,Ns,2) :: GZproj_lgr_multip  ! 
      real(8),dimension(Ns,Ns) :: GZproj_lgr_multip_  ! 
      real(8)                                :: E_Hstar,E_Hloc
      complex(8),dimension(nPhi)               :: GZvect_iter  ! GZ vector (during iterations)      
      integer :: is
      !
      !if(size(x).ne.Ns) stop "gz_energy_vdm/ wrong dimensions"
      !




      Rhop=zero
      do is=1,Ns
         vdm(is)=x(1)
      end do
      !
      select case(min_method)
      case('nlep')
         GZ_energy=gz_energy_recursive_nlep(vdm)
      case('cmin')
         GZ_energy=gz_energy_recursive_cmin(vdm)
      case('bryd')
         GZ_energy=gz_energy_broyden(vdm)
      end select
      !      
      icall = icall + 1
      write(unit_vdm_opt,'(20F18.10)') dble(icall),vdm,GZ_energy


    end function gz_energy_vdm

  end subroutine gz_optimization_vdm


  !+- TMP GALHAD
  subroutine gz_optimization_vdm_nlsq(init_vdm,optimized_vdm)
    real(8),dimension(Ns),intent(in)  :: init_vdm
    real(8),dimension(Ns),intent(out) :: optimized_vdm
    integer                           :: iter,unit_vdm_opt,icall
    real(8)                           :: GZ_energy
    !
    !LANCELOT VARIABLES
    !
    integer                           :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                           :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_
    real(8),allocatable               :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                           :: iunit,err_unit,ene_unit,Nsuccess    
    external                          :: energy_VS_vdm   

    optimized_vdm=init_vdm

    unit_vdm_opt=free_unit()
    open(unit_vdm_opt,file='vdm_optimization.out'); icall=0


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


    
    ! LANCELOT configuration parameters 
    !n_min       = Ns  ! number of minimization parameters  
    n_min       = Nvdm  ! number of minimization parameters  
    neq         = 0   ! number of equality constraints                   
    nin         = 0           ! number of in-equality constraints                   
    maxit       = 1000        ! maximum iteration number 
    gradtol     = 1.d-7       ! maximum norm of the gradient at convergence 
    feastol     = 1.d-7       ! maximum violation of parameters at convergence  
    print_level = lancelot_verbose           ! verbosity
    allocate(bl(n_min),bu(n_min),cx(neq+nin),y(neq+nin))
    bL = 3.d-1               ! lower bounds for minimization parameters
    bU = 1.d0-bL                 ! upper bounds for minimization parameters          
    !    
    call lancelot_simple(n_min,optimized_vdm,GZ_energy,exit_code,my_fun=energy_VS_vdm, &
         bl = bl, bu = bu,                                                                      &
         neq = neq, nin = nin,                                                                  &
         cx = cx, y = y, iters  = iter, maxit = maxit,                                          &
         gradtol = gradtol, feastol = feastol,                                                  &
         print_level = print_level )
    !+--------------------------------------------------------------------------------------+!    

  end subroutine gz_optimization_vdm_nlsq










  subroutine gz_optimization_simplex(simplex_init,optimized_vdm) 
    real(8),dimension(Ns+1,Ns),intent(inout) :: simplex_init
    real(8),dimension(Ns),intent(out)               :: optimized_vdm
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
    NP=Ns
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
       case('bryd')
          y(i_vertex)=gz_energy_broyden(p(i_vertex,:))
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
    case('bryd')
       call amoeba(p(1:MP,1:NP),y(1:MP),ftol,gz_energy_broyden,iter,amoeba_verbose)
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






  function gz_energy_broyden(n0)   result(GZ_energy)
    real(8),dimension(:),intent(in)        :: n0 !INPUT: Variational Density Matrix (VDM) (diagonal in istate)    
    real(8)                                :: GZ_energy !INPUT: Optimized GZ energy at fixed 
    real(8)                                :: GZ_energy_old,energy_err     ! Value of the GZ energy functional
    complex(8),dimension(Ns,Ns) :: R_init,Ropt        ! initial guess for the hopping renormalization matrix    
    complex(8),dimension(Ns,Ns)     :: slater_derivatives    
    real(8),dimension(Ns)           :: slater_lgr_multip,R_diag
    real(8),dimension(Ns,Ns)        :: GZproj_lgr_multip  ! 
    real(8)                                :: E_Hstar,E_Hloc
    complex(8),dimension(nPhi)               :: GZvect_iter  ! GZ vector (during iterations)
    !
    integer                                :: istate,iter,jstate,ifock,jfock,iphi,jphi,is
    integer                                :: unit
    logical                                :: bound

    !
    real(8),dimension(2*Ns) :: test_brIN,test_brOUT
    !
    real(8),dimension(2*Ns) :: R_broyden,test_f
    real(8),dimension(Ns) :: R_real
    real(8) :: Uin
    integer :: icall,info,i
    real(8) :: fmin
    !
    write(*,*) '********************'
    write(*,*) 'INPUT DENSITY',n0(:)
    bound=.false.
    do istate=1,Ns
       if(n0(istate).le.1.d-10.or.n0(istate).ge.1.d0-1.d-10) bound=.true.
    end do
    !
    if(.not.bound) then
       !

       !<BAUSTELLE
       R_init=0.d0 
       do istate=1,Ns
          R_init(istate,istate)=Rseed
       end do
       !

       !+- at some point I should update the routine such that all the matrix Rhop is considered
       !   in the broyden solution of the fixed point problem...
       do is=1,Ns
          R_broyden(is) = Rseed 
          R_broyden(is+Ns) = 0.d0
       end do
       !BAUSTELLE>

       icall=0

       !
       call fzero_broyden(root_function,R_broyden)
       !
       Ropt=zero
       Rseed=0.01d0
       do is=1,Ns
          Ropt(is,is) = R_broyden(is)+xi*R_broyden(is+Ns)
          R_diag(is) = dreal(Ropt(is,is))
          !Rseed = Rseed + R_diag(is)/dble(Ns)
       end do

       !+- GET GROUND STATE ENERGY FOR THE OPTIMIZED RHOP -+!
       call slater_determinant_minimization_nlep(Ropt,n0,E_Hstar,slater_lgr_multip,slater_derivatives,GZmin_verbose)       
       !
       select case(lgr_method)
       case('amoeba')
          call gz_projectors_minimization_nlep(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,GZmin_verbose)                   
       case('fsolve')
          call gz_projectors_minimization_nlep_fsolve(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,GZmin_verbose)                   
       end select
       !

       GZ_energy = E_Hstar + E_Hloc       

       test_f = root_function(R_broyden)
       write(*,*) 'ICALL',icall,test_f
       write(*,*) 'BROYDEN ENERGY',GZ_energy
       !


       !

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
    !<TMP    
    !write(*,*) 'broyden optimized energy',GZ_energy
    !TMP>
  contains
    !

    function root_function(Rhop) result(f)
      real(8),dimension(:),intent(in) :: Rhop
      real(8),dimension(size(Rhop))   :: f
      complex(8),dimension(Ns,Ns)     :: Rmatrix     
      complex(8),dimension(Ns,Ns)     :: slater_derivatives    
      complex(8),dimension(Ns,Ns)     :: Rnew ! hopping matrix renormalization (during iterations)
      real(8),dimension(Ns)           :: slater_lgr_multip,R_diag
      real(8),dimension(Ns,Ns)        :: GZproj_lgr_multip  ! 
      real(8)                         :: E_Hstar,E_Hloc
      complex(8),dimension(nPhi)      :: GZvect  ! GZ vector (during iterations)
      integer :: is
      !
      if(size(Rhop).ne.2*Ns) stop "root function/wrong dimensions in Rhop"
      !
      icall = icall+1
      Rmatrix=zero
      do is=1,Ns
         Rmatrix(is,is) = Rhop(is) +xi*Rhop(is+Ns)
      end do
      call slater_determinant_minimization_nlep(Rmatrix,n0,E_Hstar,slater_lgr_multip,slater_derivatives,GZmin_verbose)
      !+----------------------------+!
      !+- GZproj STEP MINIMIZATION -+!
      !+----------------------------+!    
      select case(lgr_method)
      case('amoeba')
         call gz_projectors_minimization_nlep(slater_derivatives,n0,E_Hloc,GZvect,GZproj_lgr_multip,GZmin_verbose)   
      case('fsolve')
         call gz_projectors_minimization_nlep_fsolve(slater_derivatives,n0,E_Hloc,GZvect,GZproj_lgr_multip,GZmin_verbose)   
      end select
      !
      Rnew=hopping_renormalization_normal(GZvect,n0)
      !      

      !

      do is=1,Ns
         f(is) = dreal(Rnew(is,is)-Rmatrix(is,is))
         f(is+Ns) = dimag(Rnew(is,is)-Rmatrix(is,is))         
      end do

      !<DEBUG
      write(44,*) icall,f(1),dreal(Rnew(1,1))
      !DEBUG>
      !
    end function root_function




  end function gz_energy_broyden












  function gz_energy_recursive_nlep(n0)   result(GZ_energy)
    real(8),dimension(:),intent(in)        :: n0 !INPUT: Variational Density Matrix (VDM) (diagonal in istate)    
    real(8)                                :: GZ_energy !INPUT: Optimized GZ energy at fixed 
    real(8)                                :: GZ_energy_old,energy_err     ! Value of the GZ energy functional
    complex(8),dimension(Ns,Ns) :: R_init        ! initial guess for the hopping renormalization matrix    
    complex(8),dimension(Ns,Ns) :: slater_derivatives    
    complex(8),dimension(Ns,Ns) :: R_iter,R_old ! hopping matrix renormalization (during iterations)
    real(8),dimension(Ns)           :: slater_lgr_multip,R_diag
    real(8),dimension(Ns,Ns) :: GZproj_lgr_multip  ! 
    real(8)                                :: E_Hstar,E_Hloc
    complex(8),dimension(nPhi)               :: GZvect_iter  ! GZ vector (during iterations)
    !
    integer                                :: istate,iter,jstate,ifock,jfock,iphi,jphi,is
    integer                                :: unit
    logical                                :: bound


    !
    write(*,*) '********************'
    write(*,*) 'INPUT DENSITY',n0(:)
    bound=.false.
    do istate=1,Ns
       if(n0(istate).lt.1.d-11.or.n0(istate).gt.1.d0-1.d-11) bound=.true.
    end do
    !
    if(.not.bound) then
       
       !+- initialize Rhop accordin to a given wanted symmetry
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

          !<DEBUG
          !write(*,*) R_iter
          !DEBUG>
          

          call slater_determinant_minimization_nlep(R_iter,n0,E_Hstar,slater_lgr_multip,slater_derivatives,GZmin_verbose)       

          !<DEBUG
          ! slater_derivatives=zero
          ! do is=1,Ns
          !    slater_derivatives(is,is) = 2.d0*R_iter(is,is)*e0test
          !    !write(*,*) slater_derivatives(is,is),R_iter(is,is),e0test
          ! end do
          !stop
          !DEBUG>

          !+----------------------------+!
          !+- GZproj STEP MINIMIZATION -+!
          !+----------------------------+!    

          select case(lgr_method)
          case('amoeba')
             call gz_projectors_minimization_nlep(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,GZmin_verbose)   
          case('fsolve')
             call gz_projectors_minimization_nlep_fsolve(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,GZmin_verbose)   
          end select
          
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
          !<DEBUG
          !write(33,*) r_diag(1)
          !DEBUG>

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

    real(8),dimension(Ns)              :: slater_lgr_multip
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
    real(8),dimension(Ns)           :: slater_lgr_multip
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

END MODULE GZ_OPTIMIZED_ENERGY




subroutine energy_VS_vdm(x,GZ_energy,i) 
  USE GZ_VARS_GLOBAL
  USE GZ_OPTIMIZED_ENERGY
  !USE GZ_PROJECTORS
  USE GZ_MATRIX_BASIS 
  implicit none

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
  !if(size(x).ne.Ns) stop "gz_energy_vdm/ wrong dimensions"  
  !size(x) == Nmin
  !
  do is=1,Ns
     imap = vdm_map(is)
     vdm(is) = x(imap)
  end do
  !
  Rhop=zero
  ! do is=1,Ns
  !    vdm(is)=x(is)
  ! end do
  !
  select case(min_method)
  case('nlep')
     GZ_energy=gz_energy_recursive_nlep(vdm)
  case('cmin')
     GZ_energy=gz_energy_recursive_cmin(vdm)
  case('bryd')
     GZ_energy=gz_energy_broyden(vdm)
  end select
  !      
end subroutine energy_VS_vdm
