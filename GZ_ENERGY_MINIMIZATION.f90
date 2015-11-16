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
  USE GZ_AUX_FUNX
  USE GZ_PROJECTORS
  !
  implicit none
  private
  !
  public :: gz_projectors_minimization_nlep,slater_determinant_minimization_nlep
  public :: slater_determinant_minimization_cmin,gz_projectors_minimization_cmin
  !
  real(8),dimension(:),allocatable,public :: vdm
  !
contains
  !
  include 'energy_functional_aux_routines.f90'
  !  
  !+-------------------------------------------------+!
  !+- NON-LINEAR EIGENVALUE PROBLEM (nlep) ROUTINES -+!
  !+-------------------------------------------------+!
  ! slater
  subroutine slater_determinant_minimization_nlep(Rhop,n0_target,Estar,lgr_multip,slater_derivatives,iverbose) 
    real(8),dimension(state_dim,state_dim),intent(in)  :: Rhop         !input:  renrmalization matrix
    real(8),dimension(state_dim),intent(in)            :: n0_target    !input:  variational density matrix
    real(8),intent(out)                                :: Estar        !output: Slater Deter GS energy
    real(8),dimension(state_dim),intent(out)           :: lgr_multip   !output: Slater Deter lagrange multipliers
    real(8),dimension(state_dim,state_dim),intent(out) :: slater_derivatives !output: Slater Deter GS energy derivatives
    logical,optional                                   :: iverbose     !input:  Verbosity level
    !
    real(8),dimension(state_dim)                       :: lgr
    real(8),dimension(state_dim,state_dim)             :: Hk,test_dens
    real(8),dimension(state_dim)                       :: ek,tmp_lgr,err_dens
    integer                                            :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb
    real(8),dimension(state_dim)                       :: lgr_multip_vec
    logical                                            :: iverbose_

    iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose

    lgr=0.d0
    call fsolve(get_delta_local_density_matrix_diag,lgr,tol=1.d-12,info=info)
    call store_slater_ground_state(Rhop,lgr,Estar,slater_derivatives)
    lgr_multip=lgr

    if(iverbose_) then
       write(*,*)
       write(*,*) "Slater Determinant: Lagrange Multipliers - diagonal form -"
       write(*,'(20F18.10)') lgr_multip
       write(*,*) "Slater Determinant: Variational density error"
       err_dens=get_delta_local_density_matrix_diag(lgr_multip)
       write(*,"(20F18.10)") err_dens
       write(*,*) "Slater Determinant: Ground State energy"
       write(*,"(20F18.10)") Estar
       write(*,*)
    end if
    !
  contains    
    !
    include 'self_minimization_slater_routines.f90'
    !
  end subroutine slater_determinant_minimization_nlep
  ! gz_projectors
  subroutine gz_projectors_minimization_nlep(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip,iverbose)
    real(8),dimension(state_dim,state_dim),intent(in) :: slater_derivatives !input:  Slater Deter GZ energy derivatives
    real(8),dimension(state_dim),intent(in)           :: n0_target          !input:  Variational density matrix
    real(8),dimension(state_dim) :: lgr
    real(8),dimension(state_dim,state_dim),intent(out):: lgr_multip         !output: GZprojectors Lagrange Multipliers -diagonal-
    real(8),dimension(nFock)                          :: GZvect   !output: GZvector
    real(8)                                           :: E_Hloc   !output: optimized local energy
    real(8)                                           :: lgr_symm(1)
    logical,optional                                  :: iverbose
    real(8),dimension(state_dim)            :: err_dens
    real(8) :: delta_out
    logical                                 :: iverbose_
    integer                                           :: info,istate,i,j

    iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose    
    !
    do istate=1,state_dim
       lgr(istate)=(n0_target(istate)-0.5d0)*2.d0       
    end do
    call fsolve(get_delta_proj_variational_density_diag,lgr,tol=1.d-15,info=info)    
    !
    lgr_multip=0.d0
    do istate=1,state_dim
       lgr_multip(istate,istate)=lgr(istate)
    end do
    call get_GZproj_ground_state(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect)
    !
    if(iverbose_) then
       write(*,*)
       write(*,*) "GZ projectors: Lagrange Parameters -diagonal case-",info
       write(*,'(10F18.10)') lgr(1:state_dim)
       err_dens=get_delta_proj_variational_density_diag(lgr)
       write(*,*) "GZ projectors: Variational density matrix error"
       write(*,'(10F18.10)') err_dens
       write(*,*) "GZ projectors: Optimized Local Energy"
       write(*,'(10F18.10)') E_Hloc
       write(*,*)
    end if
    !
  contains
    !
    include 'self_minimization_GZproj_routines.f90'
    !
  end subroutine gz_projectors_minimization_nlep


  !+--------------------------------------------+!
  !+- CONSTRAINED MINIMIZATION (cmin) ROUTINES -+!
  !+--------------------------------------------+!
  ! slater
  subroutine slater_determinant_minimization_cmin(Rhop,n0_target,Estar,lgr_multip,slater_matrix_el,iverbose) 
    real(8),dimension(state_dim,state_dim),intent(in)  :: Rhop         !input:  renrmalization matrix
    real(8),dimension(state_dim),intent(in)            :: n0_target    !input:  variational density matrix
    real(8),intent(out)                                :: Estar        !output: Slater Deter GS energy
    real(8),dimension(state_dim),intent(out)           :: lgr_multip   !output: Slater Deter lagrange multipliers
    real(8),dimension(state_dim,state_dim,Lk),intent(out) :: slater_matrix_el !output: Slater Deter GS energy derivatives
    logical,optional                                   :: iverbose     !input:  Verbosity level
    !
    real(8),dimension(state_dim)                       :: lgr
    real(8),dimension(state_dim,state_dim)             :: Hk,test_dens
    real(8),dimension(state_dim)                       :: ek,tmp_lgr,err_dens
    integer                                            :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb
    real(8),dimension(state_dim)                       :: lgr_multip_vec
    logical                                            :: iverbose_
    !
    iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose
    !
    lgr=0.d0
    call fsolve(get_delta_local_density_matrix_diag,lgr,tol=1.d-12,info=info)
    call store_slater_ground_state_cmin(Rhop,lgr,Estar,slater_matrix_el)
    lgr_multip=lgr
    !
    if(iverbose_) then
       write(*,*)
       write(*,*) "Slater Determinant: Lagrange Multipliers - diagonal form -"
       write(*,'(20F18.10)') lgr_multip
       write(*,*) "Slater Determinant: Variational density error"
       err_dens=get_delta_local_density_matrix_diag(lgr_multip)
       write(*,"(20F18.10)") err_dens
       write(*,*) "Slater Determinant: Ground State energy"
       write(*,"(20F18.10)") Estar
       write(*,*)
    end if
    !
  contains    
    !
    include 'self_minimization_slater_routines.f90'
    !
  end subroutine slater_determinant_minimization_cmin
  ! gz_projectors (lancelot subroutine)
  subroutine gz_projectors_minimization_cmin(slater_matrix_el,n0,GZvect_indep,GZenergy,GZproj_lgr_multip,GZmin_verbose)
    real(8),dimension(state_dim,state_dim,Lk),intent(in) :: slater_matrix_el
    real(8),dimension(state_dim),intent(in)              :: n0
    real(8),dimension(nFock_indep),intent(inout)         :: GZvect_indep
    real(8),intent(inout)                                :: GZenergy
    real(8),dimension(state_dim,state_dim),intent(inout) :: GZproj_lgr_multip
    logical                                              :: GZmin_verbose
    integer                                              :: iter,iter_max,istate,i,iter_
    integer                                              :: iorb,ispin,i_ind,ifock
    !
    !LANCELOT VARIABLES
    !
    integer                                              :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                                              :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_
    real(8),allocatable                                  :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                                              :: iunit,err_unit,ene_unit,Nsuccess    
    external  :: energy_GZproj_functional   

    if(allocated(slater_matrix_elements)) deallocate(slater_matrix_elements)
    if(allocated(vdm)) deallocate(vdm)
    allocate(slater_matrix_elements(state_dim,state_dim,Lk))
    slater_matrix_elements=slater_matrix_el
    allocate(vdm(state_dim)); vdm=n0
    !
    ! LANCELOT configuration parameters 
    n_min       = nFock_indep ! number of minimization parameters
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
    call lancelot_simple(n_min,GZvect_indep,GZenergy,exit_code,my_fun=energy_GZproj_functional, &
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
          write(GZmin_unit_,*) GZvect_indep(i)
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
       write(GZmin_unit_,*) 'Eiteration',Ephi
    end if
    deallocate(bl,bu,cx,y)
    deallocate(slater_matrix_elements)
  end subroutine gz_projectors_minimization_cmin
  !
END MODULE GZ_ENERGY_MINIMIZATION




subroutine energy_GZproj_functional(x,f,i)
  USE GZ_VARS_GLOBAL
  USE GZ_ENERGY_MINIMIZATION
  USE GZ_PROJECTORS
  implicit none
  !+- routine variables -+!
  real(8), intent(in)           :: x(:)
  real(8), intent(out)          :: f
  integer, intent(in), optional :: i
  real(8),allocatable           :: phi_(:)
  real(8),allocatable           :: niorb(:)
  real(8)                       :: Estar
  real(8),dimension(state_dim)  :: local_dens_
  real(8),dimension(state_dim,state_dim)  :: Rhop    
  real(8),dimension(state_dim,state_dim)  :: Hk
  !
  integer                       :: iorb,ispin,istate,jstate,ik,ifock,jfock,jorb,jspin
  !
  allocate(phi_(nFock))
  allocate(niorb(Norb))
  !
  do ifock=1,nFock
     phi_(ifock)=x(full2indep_fock(ifock))     
  end do
  !
  do iorb=1,Norb
     niorb(iorb) =  0.d0
     do ispin=1,2
        istate=index(ispin,iorb)
        niorb(iorb) = niorb(iorb) + vdm(istate)
     end do
  end do
  !
  Rhop=Rhop_matrix(phi_,vdm)
  !
  if (.not.present(i)) then
     !+- FREE ENERGY ESTIMATION WITHIN GZ APPROX -+!
     Estar=0.d0
     do ik=1,Lk
        Hk=0.d0
        do iorb=1,Norb
           do ispin=1,2
              do jorb=1,Norb
                 do jspin=1,2
                    istate=index(ispin,iorb)
                    jstate=index(jspin,jorb)               
                    ! build up the hopping hamiltonian !
                    if(ispin.eq.jspin) then
                       if(iorb.eq.jorb) then
                          Hk(istate,jstate)=epsik(ik)
                       else
                          Hk(istate,jstate)=hybik(ik)
                       end if
                    end if
                 end do
              end do
           end do
        end do
        ! hopping renormalization !
        Hk=matmul(Hk,Rhop)
        Hk=matmul(Rhop,Hk)
        do istate=1,state_dim
           do jstate=1,state_dim
              Estar = Estar + slater_matrix_elements(istate,jstate,ik)*wtk(ik)*Hk(istate,jstate)
           end do
        end do
     end do
     !
     f=Estar
     do ifock=1,nFock
        do jfock=1,nFock
           f=f+phi_(ifock)*phi_traces_basis_Hloc(ifock,jfock)*phi_(jfock)
        end do
     end do
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
              f=f+gz_local_diag(phi_,dens(istate,:,:))
           end do
           f=f-niorb(iorb)
        case(2)
           f=0.d0
           do ifock=1,nFock
              f=f+phi_(ifock)*phi_(ifock)
           end do
           f=f-1.d0
        end select
     case(2)
        select case(i)
        case(1)
           f=0.d0
           iorb=1
           do ispin=1,2
              istate=index(ispin,iorb)
              f=f+gz_local_diag(phi_,dens(istate,:,:))
           end do
           f=f-niorb(iorb)
        case(2)
           f=0.d0
           iorb=2
           do ispin=1,2
              istate=index(ispin,iorb)
              f=f+gz_local_diag(phi_,dens(istate,:,:))
           end do
           f=f-niorb(iorb)
        case(3)
           f=0.d0
           do ifock=1,nFock
              f=f+phi_(ifock)*phi_(ifock)
           end do
           f=f-1.d0
        end select        
     end select
  end if
end subroutine energy_GZproj_functional






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

