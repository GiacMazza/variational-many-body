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





subroutine get_GZproj_ground_state_fixR(n0,slater_derivatives,lgr_multip,E_Hloc,GZvect) 
  real(8),dimension(Ns)           :: n0
  complex(8),dimension(Ns,Ns) :: slater_derivatives
  real(8),dimension(Ns,Ns,2) :: lgr_multip
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
        H_projectors = H_projectors + lgr_multip(istate,jstate,1)*phi_traces_basis_dens(istate,jstate,:,:)
        !        H_projectors = H_projectors + lgr_multip(istate,jstate,2)*phi_traces_basis_Rhop(istate,jstate,:,:)
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
