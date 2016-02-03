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
