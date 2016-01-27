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
  USE MIN_AMOEBA
  !
  USE MATRIX_SPARSE
  !
  implicit none
  private
  !
  public :: gz_projectors_minimization_nlep,slater_determinant_minimization_nlep
  public :: slater_determinant_minimization_cmin,gz_projectors_minimization_cmin


  public :: gz_projectors_minimization_fixR
  public :: gz_projectors_minimization_nlep_obs

  public :: free_gz_projectors_init
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
    complex(8),dimension(Ns,Ns),intent(in)  :: Rhop         !input:  renrmalization matrix
    real(8),dimension(Ns),intent(in)     :: n0_target    !input:  variational density matrix
    real(8),intent(out)                  :: Estar        !output: Slater Deter GS energy
    real(8),dimension(Ns),intent(out)    :: lgr_multip   !output: Slater Deter lagrange multipliers
    complex(8),dimension(Ns,Ns),intent(out) :: slater_derivatives !output: Slater Deter GS energy derivatives
    logical,optional                     :: iverbose     !input:  Verbosity level
    !
    real(8),dimension(Ns)                :: lgr
    complex(8),dimension(Ns,Ns)             :: Hk
    real(8),dimension(Ns)                :: ek,tmp_lgr,err_dens
    integer                              :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb
    real(8),dimension(Ns)                :: lgr_multip_vec
    logical                              :: iverbose_

    iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose
    !
    lgr=0.d0
    call fsolve(get_delta_local_density_matrix_diag,lgr,tol=1.d-12,info=info)
    call store_slater_ground_state(Rhop,lgr,Estar,slater_derivatives)
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
  end subroutine slater_determinant_minimization_nlep
  ! gz_projectors
  subroutine gz_projectors_minimization_nlep_obs(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip,iverbose)
    complex(8),dimension(Ns,Ns),intent(in)  :: slater_derivatives !input:  Slater Deter GZ energy derivatives
    real(8),dimension(Ns),intent(in)     :: n0_target          !input:  Variational density matrix
    real(8),dimension(Ns)                :: lgr
    real(8),dimension(Ns*Ns)             :: lgr_full
    real(8),dimension(Ns,Ns),intent(out) :: lgr_multip         !output: GZprojectors Lagrange Multipliers -diagonal-
    complex(8),dimension(nPhi)              :: GZvect   !output: GZvector
    real(8)                              :: E_Hloc,Emin   !output: optimized local energy
    real(8)                              :: lgr_symm(1)
    logical,optional                     :: iverbose
    real(8),dimension(Ns)                :: err_dens
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
    do istate=1,Ns
       !lgr(istate)=(n0_target(istate)-0.5d0)*2.d0
       lgr(istate)=(0.5d0-n0_target(istate))*2.d0
    end do
    lgr=0.d0
    
    
    !<DEBUG
    write(*,*) '!+- DEBUG -+!'
    write(*,*) get_delta_proj_variational_density_diag(lgr)
    write(*,*) '!+---------+!'
    !stop
    !DEBUG>
    !
    call fsolve(get_delta_proj_variational_density_diag,lgr,tol=1.d-15,info=info)    
    !
    err_dens=get_delta_proj_variational_density_diag(lgr)
    off_constraints=0.d0
    do is=1,Ns
       !  off_constraints = off_constraints + abs(err_dens(is))**2.d0
    end do
    !<TMP TEST
    !off_constraints = 1.d-4
    !TMP TEST>
    if(off_constraints.gt.1.d-6) then 
       !if the fsolve didn't succed try with amoeba
       !minimizing the deviation from the variational target density no_target
       NP=NS
       MP=NP+1
       allocate(y(MP),p(MP,NP))
       !+- initialize simplex -+!
       p(1,:) = 0.d0
       do i=1,NS
          do j=1,Ns
             p(i+1,j) = (0.5d0-n0_target(j))*dble(i)
          end do
       end do
       !
       do i_vertex=1,MP
          y(i_vertex)=deltaVDM_proj(p(i_vertex,:))
       end do
       call amoeba_er(p(1:MP,1:NP),y(1:MP),1.d-10,deltaVDM_proj,iter,verbose=.false.)
       lgr=p(1,:) 
    end if
    !
    lgr_multip=0.d0
    do istate=1,Ns
       lgr_multip(istate,istate)=lgr(istate)
    end do
    call get_GZproj_ground_state(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect)
    !
    if(iverbose_) then
       write(*,*)
       write(*,*) "GZ projectors: Lagrange Parameters -diagonal case-"
       write(*,'(10F18.10)') lgr(1:Ns)
       !
       err_dens=get_delta_proj_variational_density_diag(lgr)
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
    include 'self_minimization_GZproj_routines.f90'
    !
  end subroutine gz_projectors_minimization_nlep_obs






  !+- AS IT IS THIS ROUTINE IS WRONG, INDEED IN ORDER TO IMPLEMENT THE FIX-R MINIMIZATION NO SLATER DERIVATIVES HAVE TO BE CONSIDERED....I'LL COME BACK LATER ON THIS... -+!
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










  subroutine gz_projectors_minimization_nlep(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip,iverbose)
    complex(8),dimension(Ns,Ns),intent(in)  :: slater_derivatives !input:  Slater Deter GZ energy derivatives
    real(8),dimension(Ns),intent(in)     :: n0_target          !input:  Variational density matrix
    real(8),dimension(Ns)                :: lgr
    real(8),dimension(Ns*Ns)             :: lgr_full
    real(8),dimension(Ns,Ns),intent(out) :: lgr_multip         !output: GZprojectors Lagrange Multipliers -diagonal-
    complex(8),dimension(nPhi)              :: GZvect   !output: GZvector
    real(8)                              :: E_Hloc,Emin   !output: optimized local energy
    real(8)                              :: lgr_symm(1)
    logical,optional                     :: iverbose
    real(8),dimension(Ns)                :: err_dens
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
    do istate=1,Ns
       lgr(istate)=(0.5d0-n0_target(istate))*2.d0
    end do
    !
    !+- LAGRANGE PARAMETER FIXING:
    !Simplex minimization minimizing of the deviation from the variational target density no_target
    NP=NS
    MP=NP+1
    allocate(y(MP),p(MP,NP))
    !+- initialize simplex -+!
    p(1,:) = 0.d0
    do i=1,NS
       do j=1,Ns
          p(i+1,j) = (0.5d0-n0_target(j))*dble(i)
       end do
    end do
    !
    do i_vertex=1,MP
       y(i_vertex)=deltaVDM_proj(p(i_vertex,:))
    end do
    call amoeba_er(p(1:MP,1:NP),y(1:MP),1.d-10,deltaVDM_proj,iter,verbose=.false.)
    lgr=p(1,:); deallocate(y,p)
 
    !
    lgr_multip=0.d0
    do istate=1,Ns
       lgr_multip(istate,istate)=lgr(istate)
    end do
    call get_GZproj_ground_state(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect)
    !
    if(iverbose_) then
       write(*,*)
       write(*,*) "GZ projectors: Lagrange Parameters -diagonal case-"
       write(*,'(10F18.10)') lgr(1:Ns)
       !
       err_dens=get_delta_proj_variational_density_diag(lgr)
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
    include 'self_minimization_GZproj_routines.f90'
    !
  end subroutine gz_projectors_minimization_nlep
  !
  subroutine free_gz_projectors_init(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip,iverbose)
    complex(8),dimension(Ns,Ns),intent(in)  :: slater_derivatives !input:  Slater Deter GZ energy derivatives
    real(8),dimension(Ns),intent(in)     :: n0_target          !input:  Variational density matrix
    real(8),dimension(Ns)                :: lgr
    real(8),dimension(Ns,Ns),intent(out) :: lgr_multip         !output: GZprojectors Lagrange Multipliers -diagonal-
    complex(8),dimension(Nphi)              :: GZvect   !output: GZvector
    real(8)                              :: E_Hloc   !output: optimized local energy
    real(8)                              :: lgr_symm(1)
    logical,optional                     :: iverbose
    real(8),dimension(Ns)                :: err_dens
    real(8)                              :: delta_out
    logical                              :: iverbose_
    integer                              :: info,istate,i,j

    iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose    
    !
    do istate=1,Ns
       lgr(istate)=(n0_target(istate)-0.5d0)*2.d0       
    end do
    call fsolve(get_delta_free_proj_variational_density_diag,lgr,tol=1.d-15,info=info)    
    !
    lgr_multip=0.d0
    do istate=1,Ns
       lgr_multip(istate,istate)=lgr(istate)
    end do
    call get_GZproj_free_ground_state(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect)
    !
    if(iverbose_) then
       write(*,*)
       write(*,*) "GZ projectors: Lagrange Parameters -diagonal case-",info
       write(*,'(10F18.10)') lgr(1:Ns)
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
  end subroutine free_gz_projectors_init





  !+--------------------------------------------+!
  !+- CONSTRAINED MINIMIZATION (cmin) ROUTINES -+!
  !+--------------------------------------------+!
  ! slater
  subroutine slater_determinant_minimization_cmin(Rhop,n0_target,Estar,lgr_multip,slater_matrix_el,iverbose) 
    complex(8),dimension(Ns,Ns),intent(in)     :: Rhop         !input:  renrmalization matrix
    real(8),dimension(Ns),intent(in)        :: n0_target    !input:  variational density matrix
    real(8),intent(out)                     :: Estar        !output: Slater Deter GS energy
    real(8),dimension(Ns),intent(out)       :: lgr_multip   !output: Slater Deter lagrange multipliers
    complex(8),dimension(Ns,Ns,Lk),intent(out) :: slater_matrix_el !output: Slater Deter GS energy derivatives
    logical,optional                        :: iverbose     !input:  Verbosity level
    !
    real(8),dimension(Ns)                   :: lgr
    complex(8),dimension(Ns,Ns)                :: Hk,test_dens
    real(8),dimension(Ns)                   :: ek,tmp_lgr,err_dens
    integer                                 :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb
    real(8),dimension(Ns)                   :: lgr_multip_vec
    
    real(8),dimension(Norb)                 :: lgr_orb,test_dens_orb
    logical                                 :: iverbose_
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
    complex(8),dimension(Ns,Ns,Lk),intent(in) :: slater_matrix_el
    real(8),dimension(Ns),intent(in)       :: n0
    complex(8),dimension(Nphi),intent(inout)  :: GZvect_indep    
    real(8),dimension(Nphi)   :: GZvect_indep_
    real(8),intent(inout)                  :: GZenergy
    real(8),dimension(Ns,Ns),intent(inout) :: GZproj_lgr_multip
    logical                                :: GZmin_verbose
    integer                                :: iter,iter_max,istate,i,iter_
    integer                                :: iorb,ispin,i_ind,ifock
    !
    !LANCELOT VARIABLES
    !
    integer                                :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                                :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_
    real(8),allocatable                    :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                                :: iunit,err_unit,ene_unit,Nsuccess    
    external                               :: energy_GZproj_functional   

    if(allocated(slater_matrix_elements)) deallocate(slater_matrix_elements)
    if(allocated(vdm)) deallocate(vdm)
    allocate(slater_matrix_elements(Ns,Ns,Lk))
    slater_matrix_elements=slater_matrix_el
    allocate(vdm(Ns)); vdm=n0
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
  end subroutine gz_projectors_minimization_cmin
  !
END MODULE GZ_ENERGY_MINIMIZATION




subroutine energy_GZproj_functional(x,f,i)
  USE GZ_VARS_GLOBAL
  USE GZ_ENERGY_MINIMIZATION
  !USE GZ_PROJECTORS
  USE GZ_MATRIX_BASIS 
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

