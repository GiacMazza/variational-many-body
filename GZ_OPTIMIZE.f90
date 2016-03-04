MODULE GZ_OPTIMIZED_ENERGY
  USE SCIFOR
  !
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_EFFECTIVE_HOPPINGS
  USE GZ_LOCAL_FOCK  
  !
  USE LANCELOT_simple_double
  !
  USE GZ_MATRIX_BASIS
  USE GZ_ENERGY_MINIMIZATION
  USE MIN_AMOEBA
  implicit none
  private
  !

  !+- OPTIMIZATION ROUTINES w/ respect the VariationalDensityMatrix (VDM) -+!
  public :: gz_optimization_vdm_simplex   ! Simplex method (AMOEBA.f90)  !
  public :: gz_optimization_vdm_nlsq      !+- constrained NonLinearLeastSquare method (GALHAD)

  !+- OPTIMIZATION considering VDM and Renormalization_matrices as free parameters -+!
  public :: gz_optimization_vdm_Rhop 
  !public :: gz_optimization_vdm_Rhop_
  !
CONTAINS
  !+-----------------------------------------------------------------------------------------------------+!
  !+- PURPOSE: Minimize the GUTZWILLER ENERGY FUNCTIONAL WITH RESPECT TO THE VARIATIONAL DENSITY MATRIX -+!
  !+-----------------------------------------------------------------------------------------------------+!
  include 'gz_optimization.f90'
  !
  subroutine gz_optimization_vdm_Rhop(init_Rhop,init_vdm) 
    complex(8),dimension(Ns,Ns) :: init_Rhop
    real(8),dimension(Ns,Ns) :: init_vdm


    real(8),dimension(Ns,Ns) :: init_lgr
    real(8),dimension(Ns) :: vdm
    complex(8),dimension(Ns,Ns) :: Rhop
    real(8) :: delta,tmp_ene
    integer :: iter    !
    !
    integer                           :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                           :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_
    real(8),allocatable               :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                           :: iunit,err_unit,ene_unit,Nsuccess  
    real(8),dimension(:),allocatable  :: xmin
    integer :: Nopt,iopt,Nslater_lgr,NRhop,Nproj_lgr
    integer :: is,js,imap

    Nslater_lgr = Nopt_diag + Nopt_odiag
    NRhop = Nopt_diag + Nopt_odiag
    Nproj_lgr = Nopt_odiag
    !
    Nopt = Nslater_lgr + 2*NRhop + Nproj_lgr
    allocate(xmin(Nopt))    

    !
    !+- initialize slater_density_constraints starting from the initial guess for the variational density_matrix
    !   
    do is=1,Ns
       vdm(is) = init_vdm(is,is)
    end do
    write(*,*) 'te perdi qua'
    call  slater_minimization_lgr(init_Rhop,vdm,tmp_ene,init_lgr)
    !
    iopt=0
    do is=1,Ns
       do js=1,Ns
          imap = opt_map(is,js)
          if(imap.gt.0) then
             xmin(imap) = init_lgr(is,js)
             xmin(imap+Nslater_lgr) = dreal(init_Rhop(is,js))
             xmin(imap+Nslater_lgr+NRhop) = dimag(init_Rhop(is,js))
             if(is.ne.js) xmin(iopt+Nslater_lgr+2*NRhop) = 0.d0
          end if
       end do
    end do
    n_min       = Nopt      ! number of minimization parameters  
    neq         = 0         ! number of equality constraints                   
    nin         = 0         ! number of in-equality constraints                   
    maxit       = 100       ! maximum iteration number 
    gradtol     = 1.d-7     ! maximum norm of the gradient at convergence 
    feastol     = 1.d-7     ! maximum violation of parameters at convergence  
    print_level = lancelot_verbose           ! verbosity
    allocate(bl(n_min),bu(n_min),cx(neq+nin),y(neq+nin))
    !    
    bL=0.d0
    bU=0.d0
    do is=1,Nslater_lgr  !+- 1,...,Nslater_lgr
       bL(is) = -10.d0               ! lower bounds for minimization parameters
       bU(is) = 10.d0                ! upper bounds for minimization parameters          
    end do
    do is=1,NRhop !+- 1,....,NRhop Real + Imag
       !
       bL(is+Nslater_lgr) =  -1.d0               ! lower bounds for minimization parameters
       bU(is+Nslater_lgr) =  1.d0                 ! upper bounds for minimization parameters          
       !
       bL(is+Nslater_lgr+NRhop) = 0.d0               ! lower bounds for minimization parameters
       bU(is+Nslater_lgr+NRhop) = 0.d0                 ! upper bounds for minimization parameters    
    end do
    do is=1,Nproj_lgr
       !           !+- 1,....,Nproj_lgr
       bL(is+Nslater_lgr+2*NRhop) = 0.d0               ! lower bounds for minimization parameters
       bU(is+Nslater_lgr+2*NRhop) = 0.d0                 ! upper bounds for minimization parameters
    end do

    delta=1.d0
    call lancelot_simple(n_min,xmin,delta,exit_code,my_fun=R_VDM_free_opt_function, &
         bl = bl, bu = bu,                                                                      &
         neq = neq, nin = nin,                                                                  &
         cx = cx, y = y, iters  = iter, maxit = maxit,                                          &
         gradtol = gradtol, feastol = feastol,                                                  &
         print_level = print_level )
    !+--------------------------------------------------------------------------------------+!    
    write(*,*) 'DELTA OPTIMIZATION',delta,exit_code,xmin
    optimization_flag=.true.
    if(allocated(GZ_vector)) deallocate(GZ_vector)
    allocate(GZ_vector(Nphi))
    call R_VDM_free_opt_function(xmin,delta)
  end subroutine gz_optimization_vdm_Rhop
  !
  !+---------------------------------------------------------------------------+!
  !+- CONSTRAINED MINIMIZATION WITH RESPECT TO THE VARIATIONAL DENSITY MATRIX -+!
  !+---------------------------------------------------------------------------+!
  !
  subroutine gz_optimization_vdm_nlsq(init_vdm,optimized_vdm)
    real(8),dimension(Nopt_diag),intent(in)  :: init_vdm
    real(8),dimension(Nopt_diag),intent(out) :: optimized_vdm
    integer                             :: iter,unit_vdm_opt,icall
    real(8)                             :: GZ_energy
    !
    !LANCELOT VARIABLES
    !
    integer                           :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                           :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_
    real(8),allocatable               :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                           :: iunit,err_unit,ene_unit,Nsuccess    

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
    n_min       = Nopt_diag  ! number of minimization parameters  
    neq         = 0   ! number of equality constraints                   
    nin         = 0           ! number of in-equality constraints                   
    maxit       = 1000        ! maximum iteration number 
    gradtol     = 1.d-7       ! maximum norm of the gradient at convergence 
    feastol     = 1.d-7       ! maximum violation of parameters at convergence  
    print_level = lancelot_verbose           ! verbosity
    allocate(bl(n_min),bu(n_min),cx(neq+nin),y(neq+nin))
    bL = 0.d0               ! lower bounds for minimization parameters
    bU = 1.d0-bL                 ! upper bounds for minimization parameters          
    !    
    call lancelot_simple(n_min,optimized_vdm,GZ_energy,exit_code,my_fun=gz_get_energy_vdm, &
         bl = bl, bu = bu,                                                                      &
         neq = neq, nin = nin,                                                                  &
         cx = cx, y = y, iters  = iter, maxit = maxit,                                          &
         gradtol = gradtol, feastol = feastol,                                                  &
         print_level = print_level )
    !+--------------------------------------------------------------------------------------+!    

    optimization_flag=.true.
    allocate(GZ_vector(Nphi))
    call gz_get_energy_vdm(optimized_vdm,GZ_energy)

  end subroutine gz_optimization_vdm_nlsq








  !+---------------------------------------------------------------------+!
  !+- SIMPLEX MINIMIZATION W/ RESPECT TO THE VARIATIONAL DENSITY MATRIX -+!
  !+---------------------------------------------------------------------+!
  subroutine gz_optimization_vdm_simplex(simplex_init,optimized_vdm) 
    real(8),dimension(Ns+1,Ns),intent(inout) :: simplex_init
    real(8),dimension(Ns),intent(out)               :: optimized_vdm
    !real(8),intent(out)                                    :: optimized_energy    
    !+- amoeba_variables-+!
    real(8),allocatable,dimension(:,:)                     :: p
    real(8),allocatable,dimension(:)                       :: y
    real(8)                                                :: ftol
    integer                                                :: np,mp,i_vertex,j_vertex,i_dim,iter
    !integer,allocatable,dimension(:)                       :: idum
    !integer                                                :: tmp_dum
    !real(8)                                                :: rnd,tot_dens
    real(8)                                                :: GZ_energy
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
       y(i_vertex) = gz_energy_vdm(p(i_vertex,:))
       write(amoeba_unit,*) p(i_vertex,:),y(i_vertex)
    end do
    !
    ftol=amoeba_min_tol
    call amoeba(p(1:MP,1:NP),y(1:MP),ftol,gz_energy_vdm,iter,amoeba_verbose)
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
    !
    optimization_flag=.true.
    if(.not.allocated(GZ_vector)) allocate(GZ_vector(Nphi))
    GZ_energy = gz_energy_vdm(optimized_vdm)
    !
  end subroutine gz_optimization_vdm_simplex












END MODULE GZ_OPTIMIZED_ENERGY
