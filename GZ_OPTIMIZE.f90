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
  public :: gz_optimization_vdm_simplex  !+- simplex minimization w/ respect to the vdm
  !public :: gz_optimization_simplex_vdm_R  !+- simplex minimization w/ respect to the vdm and R (to be CODED)
  !
  public :: gz_optimization_vdm_nlsq !+- constrained optimization w/ respect to the vdm  
  public :: gz_optimization_vdm_Rhop_nlsq !+- constrained optimization w/ respect to the vdm & Rhop (to be CODED)
  !
CONTAINS
  !+-----------------------------------------------------------------------------------------------------+!
  !+- PURPOSE: Minimize the GUTZWILLER ENERGY FUNCTIONAL WITH RESPECT TO THE VARIATIONAL DENSITY MATRIX -+!
  !+-----------------------------------------------------------------------------------------------------+!

  
  !+- TO FIX the fixR routine!!!
  subroutine gz_optimization_vdm_Rhop_nlsq(init_vdm,init_Rhop,opt_vdm,opt_Rhop) 
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
!      call slater_determinant_minimization_nlep(Rhop,vdm,E_Hstar,slater_lgr_multip,slater_derivatives,GZmin_verbose)       
      !
!      call gz_projectors_minimization_fixR(slater_derivatives,vdm,Rhop,E_Hloc,GZvect_iter,GZproj_lgr_multip,GZmin_verbose)                   
      !
      GZ_energy=E_Hstar+E_Hloc
      !
    end function gz_energy_vdm_Rhop
  end subroutine gz_optimization_vdm_Rhop_nlsq
  

  !+---------------------------------------------------------------------------+!
  !+- CONSTRAINED MINIMIZATION WITH RESPECT TO THE VARIATIONAL DENSITY MATRIX -+!
  !+---------------------------------------------------------------------------+!
  subroutine gz_optimization_vdm_nlsq(init_vdm,optimized_vdm)
    real(8),dimension(Nvdm),intent(in)  :: init_vdm
    real(8),dimension(Nvdm),intent(out) :: optimized_vdm
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
    call lancelot_simple(n_min,optimized_vdm,GZ_energy,exit_code,my_fun=gz_get_energy_vdm, &
         bl = bl, bu = bu,                                                                      &
         neq = neq, nin = nin,                                                                  &
         cx = cx, y = y, iters  = iter, maxit = maxit,                                          &
         gradtol = gradtol, feastol = feastol,                                                  &
         print_level = print_level )
    !+--------------------------------------------------------------------------------------+!    

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
       y(i_vertex) = gz_energy_vdm(p(i_vertex,:))
       write(amoeba_unit,*) p(i_vertex,:),y(i_vertex)
    end do
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
  end subroutine gz_optimization_vdm_simplex












END MODULE GZ_OPTIMIZED_ENERGY
