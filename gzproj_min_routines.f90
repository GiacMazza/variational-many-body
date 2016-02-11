!+- OBSOLETE ROUTINE TO BE REMOVED

subroutine gz_projectors_minimization_fixR(n0_target,R_target,E_Hloc,GZvect,lgr_multip_dens,lgr_multip_Rhop,iverbose)
  real(8),dimension(Ns),intent(in)     :: n0_target          !input:  Variational density matrix
  complex(8),dimension(Ns,Ns),intent(in)  :: R_target          !input:  Renormalization matrices
  real(8),dimension(2*Nvdm_c)                :: lgr
  real(8),dimension(Ns,Ns)                :: lgr_multip_dens
  complex(8),dimension(Ns,Ns)                :: lgr_multip_Rhop
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
  integer                              :: info,istate,jstate,i,j,iter
  !+- amoeba_variables-+!
  real(8),allocatable,dimension(:,:)   :: p
  real(8),allocatable,dimension(:)     :: y
  real(8)                              :: ftol,delta
  integer                              :: np,mp,i_vertex,j_vertex,i_dim,is,imap

  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose    
  !
  lgr=0.d0
  call fmin_cg(lgr,constraints_deviation,iter,delta)
  !
  lgr_multip_dens=0.d0
  do istate=1,Ns
     do jstate=1,Ns
        imap = vdm_c_map(istate,jstate)
        if(imap.gt.0) then
           lgr_multip_dens(istate,jstate)=lgr(imap)
           lgr_multip_Rhop(istate,jstate)=lgr(imap+Nvdm_c)  
        end if
     end do
  end do
  !
  call get_GZproj_ground_state_fixR(n0_target,lgr_multip_dens,lgr_multip_Rhop,E_Hloc,GZvect)
  !
  if(iverbose_) then
     write(*,*)
     write(*,*) "GZ projectors: Lagrange Parameters -diagonal case-"
     write(*,'(10F18.10)') lgr
     !
     !err_dens=constraints_deviation_amb(lgr)
     !err_dens_tmp=get_delta_proj_variational_density_diag(lgr_tmp)
     write(*,*) "GZ projectors: Variational density matrix error"
     write(*,'(10F18.10)') delta
     !
     write(*,*) "GZ projectors: Optimized Local Energy"
     write(*,'(10F18.10)') E_Hloc
     write(*,*)
  end if
  !
contains
  !
  function constraints_deviation(lm_) result(delta)
    real(8),dimension(:)                   :: lm_
    real(8)           :: delta
    real(8),dimension(Ns,Ns) :: lm_dens,lm_Rhop
    real(8),dimension(Ns,Ns) :: proj_variational_density
    complex(8),dimension(Ns,Ns) :: proj_Rhop
    complex(8),dimension(Nphi,Nphi)          :: H_projectors,H_tmp
    real(8),dimension(Nphi)               :: H_eigens
    real(8) :: tmp_ene

    integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
    integer :: iphi,jphi,imap
    !
    lm_dens=0.d0
    lm_Rhop=0.d0
    do istate=1,Ns
       do jstate=1,Ns
          imap = vdm_c_map(istate,jstate)
          if(imap.gt.0) then
             lm_dens(istate,jstate)=lm_(imap)
             lm_Rhop(istate,jstate)=lm_(imap+Nvdm_c)
          end if
       end do
    end do
    !

    !<DEBUG
    ! write(*,*) 'lm_dens'
    ! write(*,*) lm_dens
    ! write(*,*) 'n0_target'
    ! write(*,*) n0_target
    ! write(*,*) 'R_target'
    ! write(*,*) R_target
    ! write(*,*) 'lm_Rhop'
    ! write(*,*) lm_Rhop
    ! !DEBUG>
    !stop

    proj_variational_density=0.d0
    !+- build up the local H_projectors -+!
    H_projectors=zero
    H_projectors=phi_traces_basis_Hloc
    do istate=1,Ns
       do jstate=1,Ns
          !
          H_projectors = H_projectors + lm_dens(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
          H_projectors = H_projectors + lm_Rhop(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)
       end do
    end do
    !  
    call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
    !
    proj_variational_density=0.d0
    proj_Rhop=zero
    !
    do istate=1,Ns
       do jstate=1,Ns
          proj_variational_density(istate,jstate) = &
               trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))

          proj_Rhop(istate,jstate) = &
               trace_phi_basis(H_projectors(:,1),phi_traces_basis_Rhop(istate,jstate,:,:))


          proj_Rhop(istate,jstate) = proj_Rhop(istate,jstate) - & 
               R_target(istate,jstate)*sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))  !+-CHECK WELL
       end do
       proj_variational_density(istate,istate) = proj_variational_density(istate,istate) - n0_target(istate)      
    end do
    !
    delta=0.d0
    do istate=1,Ns     
       do jstate=1,Ns
          delta = delta + abs(proj_variational_density(istate,jstate))**2.d0
          delta = delta + abs(proj_Rhop(istate,jstate))**2.d0
       end do
    end do
    write(*,*) 'delta',delta,lm_
  end function constraints_deviation


  !
end subroutine gz_projectors_minimization_fixR





subroutine gz_projectors_minimization_nlep(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip,iverbose)
  complex(8),dimension(Ns,Ns),intent(in)  :: slater_derivatives !input:  Slater Deter GZ energy derivatives
  real(8),dimension(Ns),intent(in)     :: n0_target          !input:  Variational density matrix
  real(8),dimension(Nvdm)                :: lgr
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
  integer                              :: info,istate,jstate,i,j,iter
  !+- amoeba_variables-+!
  real(8),allocatable,dimension(:,:)   :: p
  real(8),allocatable,dimension(:)     :: y
  real(8)                              :: ftol,delta
  integer                              :: np,mp,i_vertex,j_vertex,i_dim,is
  integer                              :: imap

  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose    
  !

  lgr=0.d0
  call fmin_cg(lgr,get_delta_proj_variational_density,iter,delta)

  lgr_multip=0.d0
  do istate=1,Ns
     do jstate=1,Ns
        imap = vdm_c_map(istate,jstate)
        if(imap.gt.0) lgr_multip(istate,jstate)=lgr(imap)
     end do
  end do
  call get_GZproj_ground_state(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect)
  !
  if(iverbose_) then
     write(*,*)
     write(*,*) "GZ projectors: Lagrange Parameters -diagonal case-"
     do is=1,Ns
        write(*,'(20F18.10)') lgr_multip(is,:)
     end do
     !
     write(*,*) "GZ projectors: Variational density matrix error"
     write(*,'(10F18.10)') delta
     !
     write(*,*) "GZ projectors: Optimized Local Energy"
     write(*,'(10F18.10)') E_Hloc
     write(*,*)
  end if
  !
contains
  !
  function get_delta_proj_variational_density(lm_) result(delta)
    real(8),dimension(:)                   :: lm_
    real(8)           :: delta
    real(8),dimension(Ns,Ns) :: lm
    real(8),dimension(Ns,Ns) :: delta_proj_variational_density,proj_variational_density
    complex(8),dimension(Nphi,Nphi)          :: H_projectors,H_tmp
    real(8),dimension(Nphi)               :: H_eigens
    real(8) :: tmp_ene

    integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
    integer :: iphi,jphi,imap
    !
    lm=0.d0
    do istate=1,Ns
       do jstate=1,Ns
          imap = vdm_c_map(istate,jstate)
          if(imap.gt.0) lm(istate,jstate)=lm_(imap)
       end do
    end do
    !
    proj_variational_density=0.d0
    !+- build up the local H_projectors -+!
    H_projectors=zero
    H_projectors=phi_traces_basis_Hloc
    do istate=1,Ns
       do jstate=1,Ns
          H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
          H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
       end do
    end do
    !  
    call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
    !
    do istate=1,Ns
       do jstate=1,Ns
          !
          proj_variational_density(istate,jstate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
          !
       end do
    end do
    delta_proj_variational_density = proj_variational_density
    do istate=1,Ns
       delta_proj_variational_density(istate,istate) = & 
            delta_proj_variational_density(istate,istate) - n0_target(istate) 
    end do
    delta=0.d0
    do istate=1,Ns     
       do jstate=1,Ns
          delta = delta + abs(delta_proj_variational_density(istate,jstate))**2.d0
       end do
    end do
  end function get_delta_proj_variational_density
  !  include 'self_minimization_GZproj_routines.f90'
  !
end subroutine gz_projectors_minimization_nlep
!
subroutine gz_projectors_minimization_nlep_fsolve(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip,iverbose)
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
  lgr=0.d0        
  !
  call fsolve(get_delta_proj_variational_density_diag,lgr,tol=1.d-15,info=info)    
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
  !+- consider only the diagonal matrix lm(istate,istate) -+!
  function get_delta_proj_variational_density_diag(lm_) result(delta_proj_variational_density_vec)
    real(8),dimension(:)                   :: lm_
    real(8),dimension(Ns)           :: delta_proj_variational_density_vec
    real(8),dimension(Ns,Ns) :: lm
    real(8),dimension(Ns,Ns) :: delta_proj_variational_density,proj_variational_density
    complex(8),dimension(Nphi,Nphi)          :: H_projectors,H_tmp
    real(8),dimension(Nphi)               :: H_eigens
    real(8) :: tmp_ene

    integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
    integer :: iphi,jphi
    !
    !+- TEST SPARSE MATRIX
    type(sparse_matrix) :: sparseH
    logical :: test_sparse
    integer :: isparse
    !
    lm=0.d0
    do istate=1,Ns
       lm(istate,istate) = lm_(istate)
    end do
    !
    proj_variational_density=0.d0
    !+- build up the local H_projectors -+!
    H_projectors=zero
    H_projectors=phi_traces_basis_Hloc
    do istate=1,Ns
       do jstate=1,Ns
          H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
          H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
       end do
    end do
    !  
    call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
    !
    do istate=1,Ns
       do jstate=1,Ns
          !
          proj_variational_density(istate,jstate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
          !
       end do
    end do
    do istate=1,Ns     
       delta_proj_variational_density_vec(istate) = proj_variational_density(istate,istate) - n0_target(istate)      
    end do
  end function get_delta_proj_variational_density_diag


  !include 'self_minimization_GZproj_routines.f90'
  !
end subroutine gz_projectors_minimization_nlep_fsolve
!
!
subroutine free_gz_projectors_init(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip,iverbose)
  complex(8),dimension(Ns,Ns),intent(in)  :: slater_derivatives !input:  Slater Deter GZ energy derivatives
  real(8),dimension(Ns),intent(in)     :: n0_target          !input:  Variational density matrix
  real(8),dimension(Nvdm_c)                :: lgr
  real(8),dimension(Ns,Ns),intent(out) :: lgr_multip         !output: GZprojectors Lagrange Multipliers -diagonal-
  complex(8),dimension(Nphi)              :: GZvect   !output: GZvector
  real(8)                              :: E_Hloc   !output: optimized local energy
  real(8)                              :: lgr_symm(1)
  logical,optional                     :: iverbose
  real(8),dimension(Ns)                :: err_dens
  real(8)                              :: delta_out,delta
  logical                              :: iverbose_
  integer                              :: info,istate,jstate,i,j,imap,iter

  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose    
  !
  lgr=0.d0
  call fmin_cg(lgr,get_delta_proj_variational_density,iter,delta)
  !
  lgr_multip=0.d0
  do istate=1,Ns
     do jstate=1,Ns
        imap = vdm_c_map(istate,jstate)
        if(imap.gt.0) lgr_multip(istate,jstate)=lgr(imap)
     end do
  end do
  call get_GZproj_free_ground_state(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect)
  !
  if(iverbose_) then
     write(*,*)
     write(*,*) "GZ projectors: Lagrange Parameters -diagonal case-",info
     do istate=1,Ns
        write(*,'(10F18.10)') lgr_multip(istate,:)
     end do
     write(*,*) "GZ projectors: Variational density matrix error"
     write(*,'(10F18.10)') delta
     write(*,*) "GZ projectors: Optimized Local Energy"
     write(*,'(10F18.10)') E_Hloc
     write(*,*)
  end if
  !
contains
  !
  !include 'self_minimization_GZproj_routines.f90'
  function get_delta_proj_variational_density(lm_) result(delta)
    real(8),dimension(:)                   :: lm_
    real(8)           :: delta
    real(8),dimension(Ns,Ns) :: lm
    real(8),dimension(Ns,Ns) :: delta_proj_variational_density,proj_variational_density
    complex(8),dimension(Nphi,Nphi)          :: H_projectors,H_tmp
    real(8),dimension(Nphi)               :: H_eigens
    real(8) :: tmp_ene

    integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
    integer :: iphi,jphi,imap
    !
    lm=0.d0
    do istate=1,Ns
       do jstate=1,Ns
          imap = vdm_c_map(istate,jstate)
          if(imap.gt.0) lm(istate,jstate)=lm_(imap)
       end do
    end do
    !
    proj_variational_density=0.d0
    !+- build up the local H_projectors -+!
    H_projectors=zero
    H_projectors=phi_traces_basis_Hloc
    do istate=1,Ns
       do jstate=1,Ns
          H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
          H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
       end do
    end do
    !  
    call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
    !
    do istate=1,Ns
       do jstate=1,Ns
          !
          proj_variational_density(istate,jstate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
          !
       end do
    end do
    delta_proj_variational_density = proj_variational_density
    do istate=1,Ns
       delta_proj_variational_density(istate,istate) = & 
            delta_proj_variational_density(istate,istate) - n0_target(istate) 
    end do
    delta=0.d0
    do istate=1,Ns     
       do jstate=1,Ns
          delta = delta + abs(delta_proj_variational_density(istate,jstate))**2.d0
       end do
    end do
  end function get_delta_proj_variational_density
  !
end subroutine free_gz_projectors_init
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








subroutine get_GZproj_ground_state(n0,slater_derivatives,lgr_multip,E_Hloc,GZvect) 
  real(8),dimension(Ns)           :: n0
  complex(8),dimension(Ns,Ns) :: slater_derivatives
  real(8),dimension(Ns,Ns) :: lgr_multip
  real(8)                                :: E_Hloc
  complex(8),dimension(Nphi)               :: GZvect
  complex(8),dimension(Nphi,Nphi)         :: H_projectors
  real(8),dimension(Nphi)               :: H_eigens  
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                                :: iphi,jphi
  !
  !+- build up the local H_projectors -+!
  H_projectors=phi_traces_basis_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + &
             slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0(jstate)*(1.d0-n0(jstate)))
        H_projectors = H_projectors + lgr_multip(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
     end do
  end do
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !
  GZvect=H_projectors(1:Nphi,1)
  !
  E_Hloc=0.d0
  do iphi=1,Nphi
     do jphi=1,Nphi
        E_Hloc=E_Hloc+GZvect(iphi)*phi_traces_basis_Hloc(iphi,jphi)*GZvect(jphi)
     end do
  end do
  !
end subroutine get_GZproj_ground_state





subroutine get_GZproj_ground_state_(n0,dens_lgr,Rhop_lgr,E_Hloc,GZvect) 
  real(8),dimension(Ns)           :: n0
  complex(8),dimension(Ns,Ns) :: Rhop_lgr
  real(8),dimension(Ns,Ns) :: dens_lgr
  real(8)                                :: E_Hloc
  complex(8),dimension(Nphi)               :: GZvect
  complex(8),dimension(Nphi,Nphi)         :: H_projectors
  real(8),dimension(Nphi)               :: H_eigens  
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                                :: iphi,jphi
  !
  !+- build up the local H_projectors -+!
  H_projectors=phi_traces_basis_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + dens_lgr(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
        H_projectors = H_projectors + Rhop_lgr(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)
     end do
  end do

  !<DEBUG
  write(*,*) 'STO QUA'
  !write(*,*) lgr_multip
  !DEBUG>

  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !
  GZvect=H_projectors(1:Nphi,1)
  !
  E_Hloc=0.d0
  do iphi=1,Nphi
     do jphi=1,Nphi
        E_Hloc=E_Hloc+GZvect(iphi)*phi_traces_basis_Hloc(iphi,jphi)*GZvect(jphi)
     end do
  end do
  !
end subroutine get_GZproj_ground_state_




subroutine get_GZproj_free_ground_state(n0,slater_derivatives,lgr_multip,E_Hloc,GZvect) 
  real(8),dimension(Ns)           :: n0
  complex(8),dimension(Ns,Ns) :: slater_derivatives
  real(8),dimension(Ns,Ns) :: lgr_multip
  real(8)                                :: E_Hloc
  complex(8),dimension(nphi)               :: GZvect
  complex(8),dimension(Nphi,Nphi)         :: H_projectors
  real(8),dimension(Nphi)               :: H_eigens
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                                :: iphi,jphi
  !
  !+- build up the local H_projectors -+!
  H_projectors=phi_traces_basis_free_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + &
             slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0(jstate)*(1.d0-n0(jstate)))
        H_projectors = H_projectors + lgr_multip(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
     end do
  end do
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !
  GZvect=H_projectors(1:Nphi,1)
  !
  E_Hloc=0.d0
  do iphi=1,Nphi
     do jPhi=1,Nphi
        E_Hloc=E_Hloc+GZvect(iphi)*phi_traces_basis_Hloc(iphi,jphi)*GZvect(jphi)
     end do
  end do
end subroutine get_GZproj_free_ground_state





subroutine get_GZproj_ground_state_fixR(n0,lgr_multip_dens,lgr_multip_Rhop,E_Hloc,GZvect) 
  real(8),dimension(Ns)           :: n0
  real(8),dimension(Ns,Ns) :: lgr_multip_dens
  complex(8),dimension(Ns,Ns) :: lgr_multip_Rhop
  real(8)                                :: E_Hloc
  complex(8),dimension(Nphi)               :: GZvect
  complex(8),dimension(Nphi,Nphi)         :: H_projectors
  real(8),dimension(Nphi)               :: H_eigens  
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                                :: iphi,jphi
  !
  !+- build up the local H_projectors -+!
  H_projectors=phi_traces_basis_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + lgr_multip_dens(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
        H_projectors = H_projectors + lgr_multip_Rhop(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)
     end do
  end do
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !
  GZvect=H_projectors(1:Nphi,1)
  !
  E_Hloc=E_Hloc+trace_phi_basis(GZvect,phi_traces_basis_Hloc)
  !  
end subroutine get_GZproj_ground_state_fixR
