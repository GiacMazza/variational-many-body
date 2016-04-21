
  !
  subroutine gz_optimization_vdm_Rhop_ghld(init_Rhop,init_vdm) 
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
    maxit       = 300       ! maximum iteration number 
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
    if(GZmin_verbose) then
       write(*,*) 'DELTA OPTIMIZATION',delta,exit_code,xmin,cx
    end if
    optimization_flag=.true.
    if(allocated(GZ_vector)) deallocate(GZ_vector)
    allocate(GZ_vector(Nphi))
    call R_VDM_free_opt_function(xmin,delta)
  end subroutine gz_optimization_vdm_Rhop_ghld






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
    real(8),dimension(2*Ns) :: R_broyden
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

       do is=1,Ns
          R_broyden(is) = dreal(R_init(is,is))
          R_real(is) = dreal(R_init(is,is))
          R_broyden(is+Ns) = dimag(R_init(is,is))
       end do
       !BAUSTELLE>

       icall=0

       ! Rseed=0.11d0
       ! Uin =-0.1d0
       ! do i=1,1
       !    Rseed = Rseed - 0.01
       !    !call fzero_broyden(root_functionU,Uin)
       !    write(*,*) "root Uin",root_functionU(Uin)
       !    write(*,*) "root 10.d0",root_functionU(10.d0)
       !    !stop
       !    Uin=fzero_brentq(root_functionU,Uin,10.d0)
       !    write(55,*) Rseed,Uin
       ! end do
       !
       ! stop
       !
       !
       call fzero_broyden(root_function,R_broyden)
       !call fixed_point_sub(R_broyden,root_function,xtol=1.d-6)
       !call fmin_cg(R_broyden,min_function,iter,fmin)
       write(*,*) 'ICALL',icall
       !
       Ropt=zero
       Rseed=0.01d0
       do is=1,Ns
          Ropt(is,is) = R_broyden(is)+xi*R_broyden(is+Ns)
          R_diag(is) = dreal(Ropt(is,is))
          !Rseed = Rseed + R_diag(is)/dble(Ns)
       end do

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
       !



       !<BAUSTELLE
       ! write(opt_GZ_unit,*) n0
       ! write(opt_GZ_unit,*)
       ! !
       ! do iphi=1,Nphi
       !    write(opt_GZ_unit,*) GZvect_iter(iphi)
       ! end do
       ! !
       ! write(opt_GZ_unit,*)
       ! write(opt_GZ_unit,*)
       ! write(opt_energy_unit,*) n0,GZ_energy,E_Hloc,E_Hstar
       ! write(opt_rhop_unit,*) n0,R_diag(1:Ns)
       !BAUSTELLE>


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
    write(*,*) 'broyden optimized energy',GZ_energy
    !TMP>
  contains
    !
    
    function root_functionU(Uin) result(f)
      real(8),intent(in) :: Uin
      !real(8),dimension(:) :: Rhop
      real(8)   :: f
      complex(8),dimension(Ns,Ns)     :: Rmatrix     
      complex(8),dimension(Ns,Ns)     :: slater_derivatives    
      complex(8),dimension(Ns,Ns,Lk) :: slater_matrix_el    
      complex(8),dimension(Ns,Ns)     :: Rnew ! hopping matrix renormalization (during iterations)
      real(8),dimension(Ns)           :: slater_lgr_multip,R_diag
      real(8),dimension(Ns,Ns)        :: GZproj_lgr_multip  ! 
      real(8)                         :: E_Hstar,E_Hloc,GZ_energy
      complex(8),dimension(nPhi)      :: GZvect  ! GZ vector (during iterations)
      integer :: is
      !
      !
      icall = icall+1
      Rmatrix=zero
      do is=1,Ns
         Rmatrix(is,is) = Rseed
      end do

      Uloc(1)=Uin
      Uloc(2)=Uin
      Ust=Uloc(1)
      call build_local_hamiltonian
      phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
      phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)
      call slater_determinant_minimization_nlep(Rmatrix,n0,E_Hstar,slater_lgr_multip,slater_derivatives,iverbose=.false.)       

      !<DEBUG
      ! write(*,*) 'slater derivatives ok'
      ! slater_derivatives=zero
      ! do is=1,Ns
      !    slater_derivatives(is,is) = 2.d0*Rmatrix(is,is)*e0test
      !    write(*,*) slater_derivatives(is,is)
      ! end do
      !DEBUG>

      select case(lgr_method)
      case('amoeba')
         call gz_projectors_minimization_nlep(slater_derivatives,n0,E_Hloc,GZvect,GZproj_lgr_multip,iverbose=.false.)                   
      case('fsolve')
         call gz_projectors_minimization_nlep_fsolve(slater_derivatives,n0,E_Hloc,GZvect,GZproj_lgr_multip,iverbose=.false.)                   
      end select


      !stop


      ! call slater_determinant_minimization_cmin(Rmatrix,n0,E_Hstar,slater_lgr_multip,slater_matrix_el,iverbose=.false.)       
      ! GZvect=1.d0/sqrt(dble(Nphi))
      ! call gz_projectors_minimization_cmin(slater_matrix_el,n0,GZvect,GZ_energy,GZproj_lgr_multip,.true.)


      !
      Rnew=hopping_renormalization_normal(GZvect,n0)
      !      
      f=0.d0
      do is=1,Ns
         !
         !f = f + dreal(Rnew(is,is)-Rmatrix(is,is))
         !f = f + dimag(Rnew(is,is)-Rmatrix(is,is))**2.d0
         !f(is+Ns) = dimag(Rnew(is,is)-Rmatrix(is,is))         
         !
         !<DEBUG
         !write(*,*) f(is),f(is+Ns),Rmatrix(is,is),Rnew(is,is)
         !DEBUG>
      end do
      f = f + dreal(Rnew(1,1)-Rmatrix(1,1))
      write(*,*)  Uin,Rnew(1,1),Rmatrix(1,1)
      ! write(*,*) "---------"
      ! write(*,*) "---------"
      ! write(*,*) "---------"
      !
    end function root_functionU
    
    

    
    function root_function(Rhop) result(f)
      real(8),dimension(:),intent(in) :: Rhop
      !real(8),dimension(:) :: Rhop
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
      write(*,*) 'entering root finding'

      if(size(Rhop).ne.2*Ns) stop "root function/wrong dimensions in Rhop"
      !
      icall = icall+1
      Rmatrix=zero
      do is=1,Ns
         Rmatrix(is,is) = Rhop(is) +xi*Rhop(is+Ns)
      end do
      call slater_determinant_minimization_nlep(Rmatrix,n0,E_Hstar,slater_lgr_multip,slater_derivatives,iverbose=.true.)       
      !+----------------------------+!
      !+- GZproj STEP MINIMIZATION -+!
      !+----------------------------+!    
      select case(lgr_method)
      case('amoeba')
         call gz_projectors_minimization_nlep(slater_derivatives,n0,E_Hloc,GZvect,GZproj_lgr_multip,iverbose=.true.)   
      case('fsolve')
         call gz_projectors_minimization_nlep_fsolve(slater_derivatives,n0,E_Hloc,GZvect,GZproj_lgr_multip,iverbose=.true.)   
      end select
      !
      Rnew=hopping_renormalization_normal(GZvect,n0)
      !      

      !<DEBUG
      write(*,*) "calling root funciton"
      do is=1,Ns
         f(is) = dreal(Rnew(is,is)-Rmatrix(is,is))
         f(is+Ns) = dimag(Rnew(is,is)-Rmatrix(is,is))         
         !
         ! f(is) = dreal(Rnew(is,is))
         ! f(is+Ns) = dimag(Rnew(is,is))

         ! f(is) = dreal(Rnew(is,is)-Rmatrix(is,is))
         ! f(is+Ns) = dimag(Rnew(is,is)-Rmatrix(is,is))         
         !
         !<DEBUG
         write(*,*) f(is),f(is+Ns),Rmatrix(is,is),Rnew(is,is)
         !write(*,*) f(is),Rmatrix(is,is),Rnew(is,is)
         !DEBUG>
      end do
      write(44,*) icall,f(1),dreal(Rnew(1,1))
      write(*,*) "---------"
      write(*,*) "---------"
      write(*,*) "---------"
      !DEBUG>
      !
    end function root_function



    function min_function(Rhop) result(f)
      !real(8),dimension(:),intent(in) :: Rhop
      real(8),dimension(:) :: Rhop
      real(8)   :: f
      complex(8),dimension(Ns,Ns)     :: Rmatrix     
      complex(8),dimension(Ns,Ns)     :: slater_derivatives    
      complex(8),dimension(Ns,Ns)     :: Rnew ! hopping matrix renormalization (during iterations)
      real(8),dimension(Ns)           :: slater_lgr_multip,R_diag
      real(8),dimension(Ns,Ns)        :: GZproj_lgr_multip  ! 
      real(8)                         :: E_Hstar,E_Hloc
      complex(8),dimension(nPhi)      :: GZvect  ! GZ vector (during iterations)
      integer :: is
      !
      write(*,*) 'entering root finding'

      if(size(Rhop).ne.2*Ns) stop "root function/wrong dimensions in Rhop"
      !
      icall = icall+1
      Rmatrix=zero
      do is=1,Ns
         Rmatrix(is,is) = Rhop(is) +xi*Rhop(is+Ns)
      end do
      call slater_determinant_minimization_nlep(Rmatrix,n0,E_Hstar,slater_lgr_multip,slater_derivatives,iverbose=.true.)       
      !+----------------------------+!
      !+- GZproj STEP MINIMIZATION -+!
      !+----------------------------+!    
      call gz_projectors_minimization_nlep(slater_derivatives,n0,E_Hloc,GZvect,GZproj_lgr_multip,iverbose=.true.)                   
      !
      Rnew=hopping_renormalization_normal(GZvect,n0)
      !      

      write(*,*) "calling min funciton"
      f=0.d0
      do is=1,Ns
         ! f(is) = dreal(Rnew(is,is)-Rmatrix(is,is))
         ! f(is+Ns) = dimag(Rnew(is,is)-Rmatrix(is,is))         
         !
         !f(is) = abs(dreal(Rnew(is,is)))
         f = f + dreal(Rnew(is,is)-Rmatrix(is,is))**2.d0         
         f = f + dimag(Rnew(is,is)-Rmatrix(is,is))**2.d0
         ! f(is) = dreal(Rnew(is,is)-Rmatrix(is,is))
         ! f(is+Ns) = dimag(Rnew(is,is)-Rmatrix(is,is))         
         !
         !<DEBUG
         !write(*,*) f(is),f(is+Ns),Rmatrix(is,is),Rnew(is,is)
         !write(*,*) f(is),Rmatrix(is,is),Rnew(is,is)
         !DEBUG>
      end do
      write(44,*) icall,f
      write(*,*) "---------"
      write(*,*) "---------"
      write(*,*) "---------"
      !
    end function min_function



  end function gz_energy_broyden





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


  subroutine gz_get_energy_vdm_Rhop(x,GZ_energy,i) 
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
    integer :: is,js,imap
    !
    Rhop=zero
    do is=1,Ns
       imap = vdm_map(is)
       vdm(is) = x(imap)
       do js=1,Ns
          imap = vdm_c_map(is,js)
          if(imap.gt.0) Rhop(is,js) = x(imap+Nvdm) + xi*x(imap+2*Nvdm)
       end do
    end do
    !

    write(*,*) vdm
    write(*,*) Rhop
    !stop
  end subroutine gz_get_energy_vdm_Rhop





  function gz_energy_vdm_Rhop(vdm,Rhop) result(GZ_energy)
    real(8),dimension(Ns) :: vdm
    complex(8),dimension(Ns,Ns) :: Rhop
    real(8)              :: GZ_energy
    real(8),dimension(Ns,Ns)           :: slater_lgr_multip,R_diag
    real(8),dimension(Ns,Ns,2) :: GZproj_lgr_multip  ! 
    real(8),dimension(Ns,Ns) :: GZproj_lgr_multip_  ! 
    real(8)                                :: E_Hstar,E_Hloc
    complex(8),dimension(nPhi)               :: GZvect_iter  ! GZ vector (during iterations)      
    integer :: is
    !
    !call slater_determinant_minimization(Rhop,vdm,E_Hstar,slater_lgr_multip,GZmin_verbose)          
    !
    !call gz_projectors_minimization_fixR(vdm,Rhop,E_Hloc,GZvect_iter,GZproj_lgr_multip(:,:,1),GZproj_lgr_multip(:,:,2),GZmin_verbose)            

    !
    GZ_energy=E_Hstar+E_Hloc
    !
  end function gz_energy_vdm_Rhop























!# SYMMETRY STUFF
function get_mask_equal_values(xValues,mask,eps) result(Ns)
  !
  real(8),dimension(:)               :: xValues
  integer,dimension(:,:),allocatable :: mask
  real(8),optional                   :: eps
  integer                            :: Ns
  real(8)                            :: eps_
  real(8),dimension(:,:),allocatable :: tmp_mask
  integer,dimension(:),allocatable   :: search_index
  integer                            :: Nx,ix,jx,is
  real(8) :: deltaX
  !
  eps_=1.d-8
  if(present(eps)) eps_=eps
  if(allocated(mask)) deallocate(mask)
  Nx=size(xValues)
  allocate(tmp_mask(Nx,Nx),search_index(Nx))  
  tmp_mask=0
  search_index=0
  NS=0
  do ix=1,Nx
     if(search_index(ix).ge.0) then
        NS = NS + 1
        do jx=1,Nx
           deltaX = abs(Xvalues(jx)-Xvalues(ix))
           if(deltaX.lt.eps_) then
              search_index(jx) = -1
              tmp_mask(NS,jx) = 1
           end if
        end do
     end if
  end do
  allocate(mask(NS,Nx))
  do is=1,NS
     mask(is,:)=tmp_mask(is,:)
  end do
  deallocate(tmp_mask)  
  !
end function get_mask_equal_values



subroutine group_statesNS(eigenvalues)
  real(8),dimension(:,:)               :: eigenvalues
  real(8),dimension(:,:),allocatable               :: indep_eigen_
  real(8),dimension(:,:,:),allocatable               :: indep_eigen
  integer,dimension(:,:),allocatable   :: eigen_label_mask
  integer,dimension(:,:),allocatable   :: map_group_
  integer,dimension(:),allocatable     :: map_group,mult,mult_
  integer,dimension(:,:),allocatable   :: group_mask_
  integer,dimension(:,:),allocatable   :: final_mask
  integer,dimension(:,:,:),allocatable :: group_mask
  real(8),dimension(:,:),allocatable   :: group_eigenvalues
  real(8),dimension(:),allocatable   :: NS_multiplets,NSSz_multiplets
  integer,dimension(2)                 :: count_equal_states,Ngroup
  integer                              :: igroup,i,j,k
  integer                              :: M,N
  integer                              :: imap,imap_,Ntmp,tmp,imult
  logical                              :: tmp_flag
  !
  N=size(eigenvalues,1)
  M=size(eigenvalues,2)
  !
  allocate(eigen_label_mask(M,2))
  eigen_label_mask=0
  !
  !+- S2 eigenvalues -+!
  eigen_label_mask(1,1)=1
  eigen_label_mask(5,1)=1
  !+- Sz eigenvalues -+!
  eigen_label_mask(1,2)=1
  eigen_label_mask(2,2)=1
  eigen_label_mask(5,2)=1
  !
  !
  Ngroup=0
  do i=1,M
     Ngroup(1:2) = Ngroup(1:2) + eigen_label_mask(i,1:2)
  end do
  !


  write(*,*) 'CREATE MASKS'
  
  allocate(group_mask(N,N,2)); group_mask=0
  allocate(indep_eigen(N,M,2)); indep_eigen=0.d0
  do j=1,2
     allocate(group_eigenvalues(N,Ngroup(j)))
     igroup=0
     do i=1,M     
        if(eigen_label_mask(i,j).eq.1) then
           igroup=igroup+1
           group_eigenvalues(:,igroup) = eigenvalues(:,i)
        end if
     end do
     write(*,*) j
     if(.not.allocated(group_mask_))allocate(group_mask_(N,N))
     count_equal_states(j) = get_mask_equal_eigenvalues(group_eigenvalues,group_mask_,indep_eigen_)  
     do i=1,count_equal_states(j)
        group_mask(i,:,j)  = group_mask_(i,:)        
        indep_eigen(i,1:Ngroup(j),j) = indep_eigen_(i,1:Ngroup(j))
     end do
     !<TEST
     do i=1,N
        write(*,*) i,group_mask(1:count_equal_states(j),i,j)
     end do
     write(*,*)
     !TEST>
     deallocate(group_eigenvalues)
  end do  
  !
  !allocate(map_group_(N,2),map_group(N))

  write(*,*) 'CREATE MAPS'

  allocate(map_group_(N,2))
  do j=1,2     
     imap=0
     do i=1,count_equal_states(j)
        do k=1,N
           if(group_mask(i,k,j).eq.1) then
              imap=imap+1
              map_group_(imap,j) = k
           end if
        end do
     end do
     !<TEST
     do i=1,N
        write(*,'(5F7.2,I4)') eigenvalues(map_group_(i,j),:),map_group_(i,j)
     end do
     write(*,*)     
     !TEST>
  end do


  !
  ! allocate(NS_mulitplets(Ngroup(1)),NSSz_mulitplets(Ngroup(2)))  
  ! tmp_flag=.true.
  ! j=1
  ! do i=1,Ngroup(j)
  !    tmp=0
  !    do k=1,N
  !       tmp = tmp + group_mask(i,k,j)
  !       if(group_mask(i,k,j).eq.1.and.tmp_flag) then
  !          NS_multiplets
  !       end if
  !    end do
     
  ! end do

  
  write(*,*) 'STORE MULTIPLETS'


  !
  allocate(map_group(N),mult_(N))
  imap=0
  mult_=0
  imult=0
  do i=1,count_equal_states(1)
     do j=1,count_equal_states(2)
        tmp=0
        tmp_flag=.true.
        do k=1,N
           tmp = tmp + group_mask(i,k,1)*group_mask(j,k,2)
           if(group_mask(i,k,1)*group_mask(j,k,2).eq.1) then
              if(tmp_flag) then
                 imap_=k; tmp_flag=.false.
              end if
              imap=imap+1
              map_group(imap) = k                            
           end if
        end do
        !
        if(tmp.ne.0) then
           imult=imult+1
           !write(*,*) i,j,indep_eigen(imult,1:Ngroup(2),2),tmp,imap_
           !write(*,*) i,j,eigenvalues(map_group(imap_),1),eigenvalues(map_group(imap_),2),eigenvalues(map_group(imap_),5),tmp,imap_
           write(*,*) eigenvalues(imap_,1),eigenvalues(imap_,2),eigenvalues(imap_,5),tmp
           mult_(imult)=tmp
        end if
        !
     end do
  end do

  allocate(mult(imult))
  do i=1,imult
     mult(imult) = mult_(imult)
  end do

  !<TEST
  do i=1,N
     write(*,'(5F7.2,2I4)') eigenvalues(map_group(i),:),map_group(i)
  end do
  write(*,*)     
  !TEST>




end subroutine group_statesNS




subroutine basisirr_reps  
  !+- SPIN AND ORBITAL ANGULAR MOMENTUM OPERATOR -+!
  complex(8),dimension(nFock,nFock,3)     :: Svec
  complex(8),dimension(nFock,nFock)       :: S2
  !
  complex(8),dimension(nFock,nFock,3)     :: isoSvec
  complex(8),dimension(nFock,nFock)       :: isoS2
  !
  complex(8),dimension(2,2,3)             :: sigma_pauli
  complex(8),dimension(3,3,3)             :: levi_civita
  complex(8),dimension(nFock,nFock)       :: Splus,Sminus
  complex(8),dimension(nFock,nFock)       :: tmp_matrix,test,test_
  !
  complex(8),dimension(:,:,:),allocatable :: test_joint_diag
  complex(8),dimension(:,:),allocatable   :: test_jointV,jointV
  complex(8),dimension(:),allocatable   :: tmp_vec
  real(8),dimension(:,:),allocatable      :: test_jointD
  !
  real(8),dimension(nFock,nFock)          :: S2diag,Ntest
  real(8),dimension(nFock)                :: tmp
  real(8),dimension(nFock)                :: S2eigen,SZeigen
  real(8),dimension(nFock)                :: Svalue,MSvalue
  real(8),dimension(nFock,2)              :: spin_state
  real(8),dimension(2)                    :: tmp_spin_state
  integer,dimension(nFock)                :: MS,search_index
  real(8),dimension(:),allocatable        :: get_Svalue
  integer,dimension(:),allocatable        :: count_NSstates
  integer,dimension(:,:),allocatable      :: irreducible_states,tmp_irreducible_states,SZ_states
  integer                                 :: i,j,k,iorb,jorb,ispin,jspin,istate,jstate,is,iss,jfock,ifock
  real(8)                                 :: storeS,tmp_sz,deltaS,tmp_test,modV,modV_

  integer                                 :: map,NS,NMS
  integer,dimension(nFock) :: ker_map


  integer :: imap,jmap,imin,imax,jmin,jmax,ii,jj,dim_irr

  type(local_multiplets),dimension(:),allocatable :: mult_list
  type(intarray),dimension(:),allocatable :: irr_reps,irr_reps_
  integer :: Nirr_reps,jtmp,Nineq_reps
  integer,dimension(:,:),allocatable :: equ_reps,equ_reps_


  integer,dimension(:,:),allocatable :: eigen_labels
  logical :: ker_flag


  !+- build sigma pauli and levi-civita tensor -+!
  sigma_pauli=0.d0
  !
  sigma_pauli(1,2,1) = 1.d0!one
  sigma_pauli(2,1,1) = 1.d0!one
  !
  sigma_pauli(1,2,2) = -xi
  sigma_pauli(2,1,2) =  xi
  !
  sigma_pauli(1,1,3) = 1.d0!one
  sigma_pauli(2,2,3) = -1.d0!*one
  !
  !
  levi_civita=0.d0
  levi_civita(1,2,3) =  1.d0
  levi_civita(1,3,2) = -1.d0
  !
  levi_civita(2,3,1) =  1.d0
  levi_civita(3,2,1) = -1.d0
  !
  levi_civita(3,1,2) =  1.d0
  levi_civita(3,2,1) = -1.d0
  !

  S2=zero
  do i=1,3
     Svec(:,:,i)=0.d0
     do iorb=1,Norb
        do ispin=1,2
           do jspin=1,2
              istate=index(ispin,iorb)
              jstate=index(jspin,iorb)
              Svec(:,:,i) = Svec(:,:,i) + &
                   0.5d0*sigma_pauli(ispin,jspin,i)*matmul(CC(istate,:,:),CA(jstate,:,:))
           end do
        end do
     end do
     S2 = S2 + matmul(Svec(:,:,i),Svec(:,:,i))
  end do

  Splus  = Svec(:,:,1) + xi*Svec(:,:,2)
  Sminus = Svec(:,:,1) - xi*Svec(:,:,2)

  Ntest=0.d0
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        Ntest = Ntest + matmul(CC(istate,:,:),CA(istate,:,:))
     end do
  end do

  isoS2=0.d0
  do i=1,3
     isoSvec(:,:,i)=0.d0
     select case(Norb)
     case(1)        
        forall(ifock=1:nFock) isoSvec(ifock,ifock,i) = 1.d0
     case(2)
        do iorb=1,Norb
           do jorb=1,Norb
              do ispin=1,2
                 istate=index(ispin,iorb)
                 jstate=index(ispin,jorb)
                 isoSvec(:,:,i) = isoSvec(:,:,i) + &
                      0.5d0*sigma_pauli(iorb,jorb,i)*matmul(CC(istate,:,:),CA(jstate,:,:))
              end do
           end do
        end do
     case(3) 
        do iorb=1,Norb
           do jorb=1,Norb
              do ispin=1,2
                 istate=index(ispin,iorb)
                 jstate=index(ispin,jorb)
                 isoSvec(:,:,i) = isoSvec(:,:,i) + &
                      xi*levi_civita(i,iorb,jorb)*matmul(CC(istate,:,:),CA(jstate,:,:))
              end do
           end do
        end do
     end select
     isoS2 = isoS2 + matmul(isoSvec(:,:,i),isoSvec(:,:,i))
  end do

  !check isoS2 commutes with S2
  test=matmul(S2,isoS2)
  test_=matmul(isoS2,S2)
  do ifock=1,nFock
     do jfock=1,nFock
        if(abs(test(ifock,jfock)-test_(ifock,jfock)).gt.1.d-8) write(*,*) ifock,jfock
     end do
  end do
  !stop

  !< TEST JACOBI JOINT DIAGONALIZATION 
  ! tmp_matrix = Svec(:,:,2)
  ! !call matrix_diagonalize(tmp_matrix,S2eigen,'V','U')
  ! call matrix_diagonalize(tmp_matrix,S2eigen)
  ! S2diag=0.d0
  ! do i=1,nFock
  !    S2diag(i,i)=S2eigen(i)
  !    write(*,*) i,S2eigen(i),Ntest(i,i)
  ! end do

  allocate(test_joint_diag(nFock,nFock,3),test_jointV(nFock,nFock),test_jointD(nFock,3))
  test_joint_diag(:,:,1)=S2
  test_joint_diag(:,:,2)=Svec(:,:,3)
  test_joint_diag(:,:,3)=Ntest

  ! allocate(test_joint_diag(nFock,nFock,5),test_jointV(nFock,nFock),test_jointD(nFock,5))
  ! test_joint_diag(:,:,1)=S2
  ! test_joint_diag(:,:,2)=Svec(:,:,3)
  ! test_joint_diag(:,:,3)=isoS2
  ! test_joint_diag(:,:,4)=isoSvec(:,:,3)
  ! test_joint_diag(:,:,5)=Ntest


  write(*,*)
  call simultaneous_diag(test_joint_diag,test_jointV,test_jointD,eps=1.d-10) 



  allocate(eigen_labels(3,1)); eigen_labels = 0
  eigen_labels(1,1) = 1
  eigen_labels(3,1) = 1
  !
  ! eigen_labels(1,2) = 1
  ! eigen_labels(2,2) = 1
  ! eigen_labels(3,2) = 1
  !
  call get_multiplets_list(test_jointD,eigen_labels,mult_list)

  !+- I obtained the basis for irreducible representation of total-spin rotations -+!
  allocate(jointV(nFock,nFock))
  allocate(tmp_vec(nFock))
  ifock=0
  ker_map = 0
  do i=1,mult_list(1)%N_mult



     do j=1,mult_list(1)%Nequiv_mult(i)

        map = mult_list(1)%Maps(i)%index(j)
        !apply S+
        tmp_vec = test_jointV(:,map)

        ! write(*,*) map
        ! write(*,'(20F7.2)') dreal(tmp_vec)
        ! write(*,'(20F7.2)') dimag(tmp_vec)

        tmp_vec = matmul(Splus,tmp_vec)
        modV=0.d0
        do k=1,nFock
           modV=modV+tmp_vec(k)*conjg(tmp_vec(k))
        end do
        modV=sqrt(modV)

        if(abs(modV).lt.1.d-10) then
           ker_map(map) = 1
        end if

     end do
  end do

  ifock=0  
  Nirr_reps=0
  allocate(irr_reps_(nFock))
  do i=1,nFock
     allocate(irr_reps_(i)%index(4))
  end do
  allocate(equ_reps_(nFock,nFock)); equ_reps_=0
  jtmp=0
  do jj=1,mult_list(1)%N_mult
     do ii=1,mult_list(1)%Nequiv_mult(jj)
        !
        i=mult_list(1)%Maps(jj)%index(ii)        
        if(ker_map(i).eq.1) then
           !
           tmp_vec = test_jointV(:,i)        
           modV=0.d0
           do k=1,nFock
              modV=modV+tmp_vec(k)*conjg(tmp_vec(k))
           end do
           modV=sqrt(modV)                   
           !
           imin = ifock+1
           !
           dim_irr=0
           do while(modV.gt.1.d-10) 
              !
              dim_irr=dim_irr+1
              ifock = ifock + 1
              jointV(:,ifock) = tmp_vec/modV
              tmp_vec = matmul(Sminus,tmp_vec)
              !
              modV=0.d0
              do k=1,nFock
                 modV=modV+tmp_vec(k)*conjg(tmp_vec(k))
              end do
              modV=sqrt(modV)           
           end do
           imax = ifock
           j=mult_list(1)%inv_map(i)
           !
           Nirr_reps = Nirr_reps+1
           !
           if(j.eq.jtmp) then
              equ_reps_(Nirr_reps,Nineq_reps) = 1
           else
              Nineq_reps = Nineq_reps+1
              equ_reps_(Nirr_reps,Nineq_reps) = 1
           end if
           jtmp=j
           !
           irr_reps_(Nirr_reps)%index(1) = imin 
           irr_reps_(Nirr_reps)%index(2) = imax
           irr_reps_(Nirr_reps)%index(3) = dim_irr
           irr_reps_(Nirr_reps)%index(4) = mult_list(1)%inv_map(i)
           !
        end if
     end do
  end do
  write(*,*) Nirr_reps,Nineq_reps

  allocate(irr_reps(Nirr_reps))
  do i=1,Nirr_reps
     allocate(irr_reps(i)%index(4))
  end do
  allocate(equ_reps(Nirr_reps,Nineq_reps))
  !
  do i=1,Nirr_reps
     irr_reps(i)%index(:) = irr_reps_(i)%index(:)
     do j=1,Nineq_reps
        equ_reps(i,j) = equ_reps_(i,j)
     end do
  end do

  write(*,*) 'IRREDUCIBLE REPS'
  do i=1,Nirr_reps
     write(*,*) irr_reps(i)%index(:)
  end do

  write(*,*) 'EQUIVALENT IRREDUCIBLE REPS'
  do i=1,Nirr_reps
     write(*,*) equ_reps(i,:)
  end do


end subroutine basisirr_reps
