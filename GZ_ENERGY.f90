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
  public :: gz_get_energy_vdm,gz_get_energy_vdm_Rhop
  !
  public :: get_gz_ground_state_estimation !+- get the GZground state estimation at fixed n0
  !public :: get_gz_ground_state_estimation_vdm_R !+- get the GZground state estimation at fixed n0 and Rhop (TO BE CODED) ...and overload the two..
  !
  public :: gz_energy_vdm,gz_energy_vdm_Rhop
  !
contains
  !
  include "slater_min_routines.f90"
  include "gzproj_min_routines.f90"
  !
  !+- basically a post-processing routine... -+!
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
    call slater_determinant_minimization(Rhop,vdm,E_Hstar,slater_lgr_multip,GZmin_verbose)       
    !
    call gz_projectors_minimization_fixR(vdm,Rhop,E_Hloc,GZvect_iter,GZproj_lgr_multip(:,:,1),GZproj_lgr_multip(:,:,2),GZmin_verbose)                   
    !
    GZ_energy=E_Hstar+E_Hloc
    !
  end function gz_energy_vdm_Rhop



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
    GZ_energy = gz_energy_vdm_Rhop(vdm,Rhop)
    !      
  end subroutine gz_get_energy_vdm_Rhop










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













