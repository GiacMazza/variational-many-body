subroutine gz_projectors_minimization_nlep(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip,ifree,iverbose)
  complex(8),dimension(Ns,Ns),intent(in)  :: slater_derivatives !input:  Slater Deter GZ energy derivatives
  real(8),dimension(Ns),intent(in)     :: n0_target          !input:  Variational density matrix
  real(8),dimension(:),allocatable                :: lgr
  real(8),dimension(Ns*Ns)             :: lgr_full
  real(8),dimension(Ns,Ns),intent(out) :: lgr_multip         !output: GZprojectors Lagrange Multipliers -diagonal-
  complex(8),dimension(nPhi)              :: GZvect   !output: GZvector
  real(8)                              :: E_Hloc,Emin   !output: optimized local energy
  real(8)                              :: lgr_symm(1)
  logical,optional                     :: iverbose
  logical,optional                     :: ifree
  real(8),dimension(Ns)                :: err_dens
  real(8),dimension(Ns*Ns)             :: err_dens_full
  real(8)                              :: off_constraints
  real(8),allocatable,dimension(:)                              :: delta_out
  logical                              :: iverbose_,ifree_
  integer                              :: info,istate,jstate,i,j,iter
  !+- amoeba_variables-+!
  real(8),allocatable,dimension(:,:)   :: p
  real(8),allocatable,dimension(:)     :: y
  real(8)                              :: ftol,delta
  integer                              :: np,mp,i_vertex,j_vertex,i_dim,is
  integer                              :: imap

  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose    
  ifree_=.false.;if(present(ifree)) ifree_=ifree
  allocate(lgr(Nopt_diag+Nopt_odiag));  lgr=0.d0
  allocate(delta_out(Nopt_diag+Nopt_odiag))
  !
  !+- HERE Nopt = 2*Nvdm_NC_opt; allocate(lgr(Nopt))
  !
  lgr = 0.d0!lgr_init_gzproj
  select case(lgr_method)
  case('CG_min')
     call fmin_cg(lgr,get_delta_proj_variational_density,iter,delta)
  case('f_zero')
     call fsolve(fix_density,lgr,tol=1.d-10,info=iter)
     !<DEBUG
     !lgr = 0.d0
     !DEBUG>
     delta_out=fix_density(lgr)
     write(*,*) delta_out
     delta=0.d0
     do is=1,Nopt_diag+Nopt_odiag
        delta = delta + delta_out(is)**2.d0
     end do
  end select
  !
  lgr_init_gzproj = lgr
  !
  !
  !+- HERE ---> call stride_v2m(lgr,lgr_multip)
  !
  lgr_multip=0.d0
  do istate=1,Ns
     do jstate=1,Ns
        imap = opt_map(istate,jstate)
        if(imap.gt.0) lgr_multip(istate,jstate)=lgr(imap)
     end do
  end do
  !
  call gz_proj_minimization_fixed_lgr(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect,free_flag=ifree_)
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
  !<DEBUG
  !stop
  !DEBUG>
  !
contains
  !
  function get_delta_proj_variational_density(lm_) result(delta)
    real(8),dimension(:)            :: lm_
    real(8)                         :: delta
    real(8),dimension(Ns,Ns)        :: lm
    real(8),dimension(Ns,Ns)        :: delta_proj_variational_density,proj_variational_density
    complex(8),dimension(Nphi,Nphi) :: H_projectors,H_tmp
    real(8),dimension(Nphi)         :: H_eigens
    real(8)                         :: tmp_ene

    integer                         :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
    integer                         :: iphi,jphi,imap
    !
    lm=0.d0
    do istate=1,Ns
       do jstate=1,Ns
          imap = opt_map(istate,jstate)
          if(imap.gt.0) lm(istate,jstate)=lm_(imap)
       end do
    end do
    !
    !
    !+- HERE use stride_v2m
    !
    proj_variational_density=0.d0
    !+- build up the local H_projectors -+!
    !+- put all this story inside a routine 
    ! H_projectors=zero
    ! H_projectors=phi_traces_basis_Hloc
    ! if(ifree_) H_projectors = phi_traces_basis_free_Hloc
    ! do istate=1,Ns
    !    do jstate=1,Ns
    !       H_projectors = H_projectors + 2.d0*slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))          
    !       H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
    !    end do
    ! end do
    !
    call build_H_GZproj(H_projectors,slater_derivatives,n0_target,lm,ifree_)
    !  
    !call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
    call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
    !
    do istate=1,Ns
       do jstate=1,Ns
          !
          proj_variational_density(istate,jstate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
          !write(*,*) proj_variational_density(istate,jstate)
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
    !write(*,*) lm_,delta
  end function get_delta_proj_variational_density


  function fix_Density(lm_) result(delta)
    real(8),dimension(:)         :: lm_
    real(8),dimension(size(lm_)) :: delta
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
          imap = opt_map(istate,jstate)
          if(imap.gt.0) lm(istate,jstate)=lm_(imap)
       end do
    end do
    !
    !+- here use stride_v2m
    !

    !
    proj_variational_density=0.d0
    !+- build up the local H_projectors -+!
    !
    call build_H_GZproj(H_projectors,slater_derivatives,n0_target,lm,ifree_)
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
    !
    delta = 0.d0
    do istate=1,Ns
       do jstate=1,Ns
          imap = opt_map(istate,jstate)   !+- here we have a problem...  !-----> SOLUTION: use stride and then check for it!
          if(imap.gt.0) delta(imap) = delta(imap) + delta_proj_variational_density(istate,jstate)**2.d0
       end do
    end do
    !
  end function fix_Density
  !
end subroutine gz_projectors_minimization_nlep
!
subroutine gz_projectors_minimization_lgr_(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip,ifree,iverbose)
  complex(8),dimension(Ns,Ns),intent(in)  :: slater_derivatives !input:  Slater Deter GZ energy derivatives
  real(8),dimension(Ns),intent(in)     :: n0_target          !input:  Variational density matrix
  real(8),dimension(:),allocatable                :: lgr
  complex(8),dimension(Ns,Ns),intent(out) :: lgr_multip         !output: GZprojectors Lagrange Multipliers -diagonal-
  complex(8),dimension(nPhi)              :: GZvect   !output: GZvector
  real(8)                              :: E_Hloc,Emin   !output: optimized local energy
  logical,optional                     :: iverbose
  logical,optional                     :: ifree
  real(8),dimension(Ns)                :: err_dens
  real(8)                              :: off_constraints
  real(8),allocatable,dimension(:)                              :: delta_out
  logical                              :: iverbose_,ifree_
  integer                              :: info,istate,jstate,i,j,iter,i0
  !+- amoeba_variables-+!
  real(8),allocatable,dimension(:,:)   :: p
  real(8),allocatable,dimension(:)     :: y
  real(8)                              :: ftol,delta
  integer                              :: np,mp,i_vertex,j_vertex,i_dim,is
  integer                              :: imap,jmap  
  integer                              :: Nopt
  complex(8),dimension(:),allocatable  :: lgr_cmplx

  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose    
  ifree_=.false.;if(present(ifree)) ifree_=ifree
  !
  Nopt = 2*Nvdm_NC_opt; allocate(lgr(Nopt)); allocate(delta_out(Nopt))  
  lgr=0.d0  
  select case(lgr_method)
  case('CG_min')
     call fmin_cg(lgr,get_delta_proj_variational_density,iter,delta)
  case('f_zero')
     call fsolve(fix_density,lgr,tol=1.d-10,info=iter)
     delta_out=fix_density(lgr)
     delta=0.d0
     do is=1,2*Nopt_lgr
        delta = delta + delta_out(is)**2.d0
     end do
  end select
  !
  allocate(lgr_cmplx(Nvdm_NC_opt))
  do i=1,Nvdm_NC_opt
     lgr_cmplx(i) = lgr(i)+xi*lgr(i+Nvdm_NC_opt)
  end do
  call vdm_NC_stride_v2m(lgr_cmplx,lgr_multip)
  !
  call gz_proj_minimization_fixed_lgr_(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect,free_flag=ifree_)
  !
  if(iverbose_) then
     write(*,*)
     write(*,*) "GZ projectors: Lagrange Parameters -diagonal case-"
     do is=1,Ns
        write(*,'(20F8.4)') dreal(lgr_multip(is,:)),dimag(lgr_multip(is,:))
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
    complex(8),dimension(:),allocatable :: lm_cmplx,delta_cmplx
    complex(8),dimension(Ns,Ns) :: lm
    complex(8),dimension(Ns,Ns) :: delta_proj_variational_density,proj_variational_density
    complex(8),dimension(Nphi,Nphi)          :: H_projectors,H_tmp
    complex(8),dimension(Nphi)          :: proj_gs
    real(8),dimension(Nphi)               :: H_eigens
    real(8) :: tmp_ene,tmp_eigen

    integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
    integer :: iphi,jphi,imap,jmap,js,is,i,i0,Ndegen,Nopt_cmplx
    !
    allocate(lm_cmplx(Nvdm_NC_opt))
    do i=1,Nvdm_NC_opt
       lm_cmplx(i) = lm_(i)+xi*lm_(i+Nvdm_NC_opt)
    end do
    call vdm_NC_stride_v2m(lm_cmplx,lm)    
    deallocate(lm_cmplx)
    !
    proj_variational_density=0.d0
    !+- build up the local H_projectors -+!
    call build_H_GZproj_(H_projectors,slater_derivatives,n0_target,lm,ifree_)
    !  
    call matrix_diagonalize(H_projectors,H_eigens)         
    !
    proj_gs = H_projectors(:,1)
    !
    do is=1,Ns
       do js=1,Ns
          proj_variational_density(is,js) = trace_phi_basis(proj_gs,phi_traces_basis_dens(is,js,:,:))
       end do
    end do
    delta_proj_variational_density = proj_variational_density
    do istate=1,Ns
       delta_proj_variational_density(istate,istate) = & 
            delta_proj_variational_density(istate,istate) - n0_target(istate) 
    end do
    !
    delta=0.d0
    do istate=1,Ns     
       do jstate=1,Ns
          delta = delta + delta_proj_variational_density(istate,jstate)*conjg(delta_proj_variational_density(istate,jstate))
       end do
    end do

  end function get_delta_proj_variational_density
  !
  function fix_density(lm_) result(delta)
    real(8),dimension(:)         :: lm_
    real(8),dimension(size(lm_)) :: delta
    complex(8),dimension(:),allocatable :: lm_cmplx,delta_cmplx
    complex(8),dimension(Ns,Ns) :: lm
    complex(8),dimension(Ns,Ns) :: delta_proj_variational_density,proj_variational_density
    complex(8),dimension(Nphi,Nphi)          :: H_projectors,H_tmp
    complex(8),dimension(Nphi)          :: proj_gs
    real(8),dimension(Nphi)               :: H_eigens
    real(8) :: tmp_ene,tmp_eigen

    integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
    integer :: iphi,jphi,imap,jmap,js,is,i,i0,Ndegen,Nopt_cmplx
    !
    allocate(lm_cmplx(Nvdm_NC_opt))
    do i=1,Nvdm_NC_opt
       lm_cmplx(i) = lm_(i)+xi*lm_(i+Nvdm_NC_opt)
    end do
    call vdm_NC_stride_v2m(lm_cmplx,lm)    
    deallocate(lm_cmplx)
    !
    proj_variational_density=0.d0
    !+- build up the local H_projectors -+!
    call build_H_GZproj_(H_projectors,slater_derivatives,n0_target,lm,ifree_)
    !  
    call matrix_diagonalize(H_projectors,H_eigens)         
    !
    proj_gs = H_projectors(:,1)
    !
    do is=1,Ns
       do js=1,Ns
          proj_variational_density(is,js) = trace_phi_basis(proj_gs,phi_traces_basis_dens(is,js,:,:))
       end do
    end do
    delta_proj_variational_density = proj_variational_density
    do istate=1,Ns
       delta_proj_variational_density(istate,istate) = & 
            delta_proj_variational_density(istate,istate) - n0_target(istate) 
    end do
    !    
    delta=0.d0
    allocate(delta_cmplx(Nvdm_NC_opt))
    call vdm_NC_stride_m2v(delta_proj_variational_density,delta_cmplx)
    do i=1,Nvdm_NC_opt
       delta(i) = dreal(delta_cmplx(i))
       delta(i+Nvdm_NC_opt) = dimag(delta_cmplx(i))
    end do
    deallocate(delta_cmplx)
  end function fix_density
  !
end subroutine gz_projectors_minimization_lgr_
!
subroutine gz_projectors_minimization_cmin(slater_matrix_el,n0,GZenergy,GZvect_indep,GZproj_lgr_multip,iverbose)
  complex(8),dimension(Ns,Ns,Lk),intent(in) :: slater_matrix_el
  real(8),dimension(Ns),intent(in)       :: n0
  complex(8),dimension(Nphi),intent(inout)  :: GZvect_indep    
  real(8),dimension(Nphi)   :: GZvect_indep_
  real(8),intent(inout)                  :: GZenergy
  real(8),dimension(Ns,Ns),intent(inout) :: GZproj_lgr_multip
  logical,optional                       :: iverbose
  logical                       :: iverbose_
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

  real(8) :: tmp_test_constraint
  !
  iverbose_=.false. ; if(present(iverbose)) iverbose_=iverbose
  !
  vdm=n0
  !
  ! LANCELOT configuration parameters 
  n_min       = Nphi ! number of minimization parameters
  neq         = Ns*Ns+1!Nopt_diag+Nopt_odiag + 1    ! number of equality constraints         
  nin         = 0           ! number of in-equality constraints                   
  maxit       = 1000        ! maximum iteration number 
  gradtol     = 1.d-7       ! maximum norm of the gradient at convergence 
  feastol     = 1.d-7       ! maximum violation of parameters at convergence  
  print_level = lancelot_verbose           ! verbosity
  allocate(bl(n_min),bu(n_min),cx(neq+nin),y(neq+nin))
  bL = 0.d0               ! lower bounds for minimization parameters
  bU = 1.d0-bL                 ! upper bounds for minimization parameters          
  if(iverbose_) then
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
     GZvect_indep = 1.d0/sqrt(dble(n_min))
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
  !<DEBUG
  ! write(*,'(20F7.3)') cx  
  ! tmp_test_constraint=0.d0
  ! do i=1,n_min
  !    tmp_test_constraint = tmp_test_constraint + GZvect_indep_(i)*GZvect_indep_(i)
  ! end do
  ! write(*,*) tmp_test_constraint
  ! write(*,*)
  !DEBUG>
  GZproj_lgr_multip=0.d0
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        GZproj_lgr_multip(istate,istate)=cx(iorb)
     end do
  end do
  !
  if(iverbose_) then
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
    integer                      :: is,js,imap
    logical :: c_flag
    !
    f=1.d0
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
    Rhop=hopping_renormalization_normal(phi_,vdm)
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
                Estar = Estar + slater_matrix_el(istate,jstate,ik)*wtk(ik)*Hk(istate,jstate)
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
       imap = 0
       c_flag=.false.
       do is=1,Ns
          do js=1,Ns
             if(is.eq.js) then
                !
                f=0.d0
                f = f + trace_phi_basis(phi_,phi_traces_basis_dens(is,js,:,:))
                f = f - vdm(is)
                !
             else
                !
                f = 0.d0
                f = f + trace_phi_basis(phi_,phi_traces_basis_dens(is,js,:,:))
                !
             end if
          end do
       end do
       !
       if(i.eq.Nopt_diag+Nopt_odiag+1) then
          f=0.d0
          do iphi=1,Nphi
             f = f + conjg(phi_(iphi))*phi_(iphi)
          end do
          !
          f=f-1.d0
          !
       end if
    end if
  end subroutine energy_GZproj_functional


end subroutine gz_projectors_minimization_cmin
!





!
subroutine get_GZproj_ground_state(n0,slater_derivatives,lgr_multip,E_Hloc,GZvect,free_flag) 
  real(8),dimension(Ns)           :: n0
  complex(8),dimension(Ns,Ns) :: slater_derivatives
  real(8),dimension(Ns,Ns) :: lgr_multip
  real(8)                                :: E_Hloc
  complex(8),dimension(Nphi)               :: GZvect
  logical,optional :: free_flag
  logical :: free_flag_
  complex(8),dimension(Nphi,Nphi)         :: H_projectors
  real(8),dimension(Nphi)               :: H_eigens  
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                                :: iphi,jphi,is

  real(8),dimension(Ns,Ns) :: proj_variational_density
  !
  !+- build up the local H_projectors -+!
  free_flag_=.false.; if(present(free_flag)) free_flag_=free_flag
  !
  call build_H_GZproj(H_projectors,slater_derivatives,n0,lgr_multip)
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !
  GZvect=H_projectors(1:Nphi,1)
  !
  E_Hloc=trace_phi_basis(GZvect,phi_traces_basis_Hloc)
  !
  !<DEBUG
  write(*,*) 'tmp_test GZproj_variational_density'
  do istate=1,Ns
     do jstate=1,Ns
        !
        proj_variational_density(istate,jstate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
        !
     end do
     write(*,*) proj_variational_density(istate,:)
  end do
  write(*,*)
  do istate=1,Nphi
     write(*,*) H_eigens(istate),H_projectors(istate,1)
  end do
  do is=1,Ns
     write(*,*) slater_derivatives(is,:)
  end do
  !DEBUG>
  !
end subroutine get_GZproj_ground_state


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
        H_projectors = H_projectors + lgr_multip_dens(istate,jstate)*phi_traces_basis_dens_hc(istate,jstate,:,:)
        H_projectors = H_projectors - lgr_multip_Rhop(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)
        H_projectors = H_projectors - conjg(lgr_multip_Rhop(istate,jstate))*phi_traces_basis_Rhop_hc(istate,jstate,:,:)
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





subroutine gz_proj_minimization_fixed_lgr_(n0,slater_derivatives,lgr_multip,E_Hloc,GZvect,free_flag) 
  real(8),dimension(Ns)           :: n0
  complex(8),dimension(Ns,Ns) :: slater_derivatives
  complex(8),dimension(Ns,Ns) :: lgr_multip
  real(8)                                :: E_Hloc
  complex(8),dimension(Nphi)               :: GZvect
  logical,optional :: free_flag
  logical :: free_flag_
  complex(8),dimension(Nphi,Nphi)         :: H_projectors
  real(8),dimension(Nphi)               :: H_eigens  
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                                :: iphi,jphi,is

  real(8),dimension(Ns,Ns) :: proj_variational_density
  !
  !+- build up the local H_projectors -+!
  free_flag_=.false.; if(present(free_flag)) free_flag_=free_flag
  !
  call build_H_GZproj_(H_projectors,slater_derivatives,n0,lgr_multip)
  !
  call matrix_diagonalize(H_projectors,H_eigens)         
  !
  GZvect=H_projectors(1:Nphi,1)
  !
  E_Hloc=trace_phi_basis(GZvect,phi_traces_basis_Hloc)
  !
end subroutine gz_proj_minimization_fixed_lgr_


subroutine gz_proj_minimization_fixed_Rhop(n0,lgr_multip,lgr_multip_Rhop,E_Hloc,GZvect,free_flag) 
  real(8),dimension(Ns)           :: n0
  complex(8),dimension(Ns,Ns)   :: lgr_multip
  complex(8),dimension(Ns,Ns)     :: lgr_multip_Rhop
  real(8)                         :: E_Hloc
  complex(8),dimension(Nphi)      :: GZvect
  logical,optional                :: free_flag
  logical                         :: free_flag_
  complex(8),dimension(Nphi,Nphi) :: H_projectors
  real(8),dimension(Nphi)         :: H_eigens  
  integer                         :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                         :: iphi,jphi,is,js
  real(8),dimension(Ns,Ns)        :: proj_variational_density
  !
  !+- build up the local H_projectors -+!
  free_flag_=.false.; if(present(free_flag)) free_flag_=free_flag
  !
  H_projectors=zero
  H_projectors = phi_traces_basis_free_Hloc
  if(.not.free_flag_) H_projectors = phi_traces_basis_Hloc
  do is=1,Ns
     do js=1,Ns
        ! !        
        H_projectors = H_projectors + lgr_multip(is,js)*phi_traces_basis_dens(is,js,:,:)
        H_projectors = H_projectors + conjg(lgr_multip(is,js))*phi_traces_basis_dens_hc(is,js,:,:)
        ! !
        H_projectors = H_projectors - lgr_multip_Rhop(is,js)*phi_traces_basis_Rhop(is,js,:,:)
        H_projectors = H_projectors - conjg(lgr_multip_Rhop(is,js))*phi_traces_basis_Rhop_hc(is,js,:,:)
     end do
  end do
  !
  call matrix_diagonalize(H_projectors,H_eigens)
  !
  GZvect=H_projectors(1:Nphi,1)
  E_Hloc=trace_phi_basis(GZvect,phi_traces_basis_Hloc)
  !
end subroutine gz_proj_minimization_fixed_Rhop











!########################################################################!
!########################################################################!
!########################################################################!
!########################################################################!
!                                                                        !
!                           SUPERCONDUCTING ROUTINES                     !
!                                                                        !
!########################################################################!
!########################################################################!
!########################################################################!
!########################################################################!
subroutine gz_projectors_minimization_nlep_superc(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip,ifree,iverbose)
  complex(8),dimension(2,Ns,Ns),intent(in)  :: slater_derivatives !input:  Slater Deter GZ energy derivatives
  real(8),dimension(Ns),intent(in)     :: n0_target          !input:  Variational density matrix
  real(8),dimension(:),allocatable                :: lgr,lgr_normal,lgr_anomalous
  real(8),dimension(Ns*Ns)             :: lgr_full
  complex(8),dimension(2,Ns,Ns),intent(out) :: lgr_multip         !output: GZprojectors Lagrange Multipliers -diagonal-
  complex(8),dimension(nPhi)              :: GZvect   !output: GZvector
  real(8)                              :: E_Hloc,Emin   !output: optimized local energy
  real(8)                              :: lgr_symm(1)
  logical,optional                     :: iverbose
  logical,optional                     :: ifree
  real(8),dimension(Ns)                :: err_dens
  real(8),dimension(Ns*Ns)             :: err_dens_full
  real(8)                              :: off_constraints
  real(8),allocatable,dimension(:)                              :: delta_out
  logical                              :: iverbose_,ifree_
  integer                              :: info,istate,jstate,i,j,iter,i0
  !+- amoeba_variables-+!
  real(8),allocatable,dimension(:,:)   :: p
  real(8),allocatable,dimension(:)     :: y
  real(8)                              :: ftol,delta
  integer                              :: np,mp,i_vertex,j_vertex,i_dim,is
  integer                              :: imap,jmap  
  integer                              :: Nopt
  complex(8),dimension(:),allocatable  :: lgr_cmplx

  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose    
  ifree_=.false.;if(present(ifree)) ifree_=ifree
  ! allocate(lgr(Nopt_lgr));  lgr=0.d0
  ! allocate(delta_out(Nopt_lgr))
  !
  ! allocate(lgr(2*Nopt_lgr));  lgr=0.d0
  ! allocate(delta_out(2*Nopt_lgr))

  ! allocate(lgr_normal(2*Nopt_normal))
  ! allocate(lgr_anomalous(2*Nopt_anomalous))
  !
  !
  !lgr = 0.d0!lgr_init_gzproj

  !<FUCKING HARD DEBUGGING
  ! lgr(1)=-0.
  ! do i=1,60
  !    lgr(1)=lgr(1)-0.01
  !
  !    delta_out = fix_density_(lgr)
  !    delta=0.d0
  !    do is=1,2*Nopt_lgr
  !       delta = delta + delta_out(is)**2.d0
  !    end do
  !    write(300,*) lgr(1),delta
  !    delta = get_delta_proj_variational_density_(lgr)
  !    write(200,*) lgr(1),delta
  ! end do
  !  stop
  !FUCKING HARD DEBUGGING>
  !+- HERE something like Nopt = 2*Nvdm_NC_opt + 2*Nvdm_AC_opt; allocate(lgr(Nopt))  
  !
  Nopt = 2*Nvdm_NC_opt + 2*Nvdm_AC_opt; allocate(lgr(Nopt))
  lgr=0.d0  
  select case(lgr_method)
  case('CG_min')
     call fmin_cg(lgr,get_delta_proj_variational_density,iter,delta)
  case('f_zero')
     call fsolve(fix_density,lgr,tol=1.d-10,info=iter)
     delta_out=fix_density(lgr)
     !write(*,*) delta_out
     delta=0.d0
     do is=1,2*Nopt_lgr
        delta = delta + delta_out(is)**2.d0
     end do
  end select
  !
  !lgr_init_gzproj = lgr
  !
  !
  !+- HERE something like call vdm_NC_stride_v2m(lgr(1:2*Nvdm_NC_opt,lgr_multip(1,:,:)))    OF COURSE CHANGE ALL THIS STUFF INSIDE THE ROUTINES fix_density and get_delta!!!!!!
  !                       call vdm_NC_stride_v2m(lgr(2*Nvdm_NC_opt+1:Nopt,lgr_multip(2,:,:)))
  !
  ! lgr_multip=0.d0
  ! do istate=1,Ns
  !    do jstate=1,Ns
  !       imap = opt_map(istate,jstate)
  !       jmap = opt_map_anomalous(istate,jstate)
  !       if(imap.gt.0) lgr_multip(1,istate,jstate) = lgr(imap) + xi*lgr(imap+Nopt_lgr)
  !       if(jmap.gt.0) lgr_multip(2,istate,jstate) = lgr(jmap+Nopt_normal) + xi*lgr(jmap+Nopt_normal+Nopt_lgr)
  !    end do
  ! end do
  !

  allocate(lgr_cmplx(Nvdm_NC_opt))
  do i=1,Nvdm_NC_opt
     lgr_cmplx(i) = lgr(i)+xi*lgr(i+Nvdm_NC_opt)
  end do
  call vdm_NC_stride_v2m(lgr_cmplx,lgr_multip(1,:,:))
  i0=2*Nvdm_NC_opt
  deallocate(lgr_cmplx)
  allocate(lgr_cmplx(Nvdm_AC_opt))  
  do i=1,Nvdm_AC_opt
     lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_AC_opt)
  end do
  call vdm_AC_stride_v2m(lgr_cmplx,lgr_multip(2,:,:))

  !
  ! call vdm_NC_stride_v2m(lgr(1:2*Nvdm_NC_opt,lgr_multip(1,:,:)))    
  ! call vdm_AC_stride_v2m(lgr(2*Nvdm_NC_opt+1:Nopt,lgr_multip(2,:,:)))
  !
  ! do i=1,Nopt_normal
  !    lgr_normal(i)=lgr(i)
  !    lgr_normal(i+Nopt_normal)=lgr(i+Nopt_lgr)       
  ! end do
  ! do i=1,Nopt_anomalous
  !    lgr_anomalous(i)=lgr(i+Nopt_normal)
  !    lgr_anomalous(i+Nopt_anomalous)=lgr(i+Nopt_normal+Nopt_lgr)       
  ! end do
  ! lgr_multip(1,:,:) = symmetry_stride_vec_mat(lgr_normal)
  ! lgr_multip(2,:,:) = symmetry_stride_vec_mat_(lgr_anomalous)
  !
  call gz_proj_minimization_fixed_lgr_superc(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect,free_flag=ifree_)
  !
  if(iverbose_) then
     write(*,*)
     write(*,*) "GZ projectors: Lagrange Parameters -diagonal case-"
     do is=1,Ns
        write(*,'(20F8.4)') dreal(lgr_multip(1,is,:)),dimag(lgr_multip(1,is,:))
     end do
     write(*,*)
     do is=1,Ns
        write(*,'(20F8.4)') dreal(lgr_multip(2,is,:)),dimag(lgr_multip(2,is,:))
     end do
     !
     write(*,*) "GZ projectors: Variational density matrix error"
     write(*,'(10F18.10)') delta!,delta_out
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
          imap = opt_map(istate,jstate)
          if(imap.gt.0) lm(istate,jstate)=lm_(imap)
       end do
    end do
    !
    proj_variational_density=0.d0
    !
    call build_H_GZproj(H_projectors,slater_derivatives,n0_target,lm,ifree_)
    !  
    call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
    !
    do istate=1,Ns
       do jstate=1,Ns
          !
          proj_variational_density(istate,jstate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
          !write(*,*) proj_variational_density(istate,jstate)
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
    !write(*,*) lm_,delta
  end function get_delta_proj_variational_density


  function fix_density(lm_) result(delta)
    real(8),dimension(:)         :: lm_
    real(8),dimension(size(lm_)) :: delta
    complex(8),dimension(:),allocatable :: lm_cmplx,delta_cmplx
    complex(8),dimension(2,Ns,Ns) :: lm
    complex(8),dimension(2,Ns,Ns) :: delta_proj_variational_density,proj_variational_density
    complex(8),dimension(Nphi,Nphi)          :: H_projectors,H_tmp
    complex(8),dimension(Nphi)          :: proj_gs
    real(8),dimension(Nphi)               :: H_eigens
    real(8) :: tmp_ene,tmp_eigen

    integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
    integer :: iphi,jphi,imap,jmap,js,is,i,i0,Ndegen,Nopt_cmplx
    !
    !
    !    if(size(lm_).ne.2*Nopt_lgr) stop 'wrong lgr_parameters size @ gz_projectors_minimization_nlep_superc/fix_density'   
    ! Nopt_cmplx=Nopt/2
    ! do i=1,Nopt_cmplx
    !    lm_cmplx(i) = lm_(i) + xi*lm_(i+Nopt_cmplx)
    ! end do
    !
    !
    allocate(lm_cmplx(Nvdm_NC_opt))
    do i=1,Nvdm_NC_opt
       lm_cmplx(i) = lm_(i)+xi*lm_(i+Nvdm_NC_opt)
    end do
    call vdm_NC_stride_v2m(lm_cmplx,lm(1,:,:))    
    deallocate(lm_cmplx)
    i0=2*Nvdm_NC_opt
    allocate(lm_cmplx(Nvdm_AC_opt))
    do i=1,Nvdm_AC_opt
       lm_cmplx(i) = lm_(i0+i)+xi*lm_(i0+i+Nvdm_AC_opt)
    end do
    call vdm_AC_stride_v2m(lm_cmplx,lm(2,:,:))    
    deallocate(lm_cmplx)
    !
    ! lm=0.d0
    ! do istate=1,Ns
    !    do jstate=1,Ns
    !       imap = opt_map(istate,jstate)
    !       jmap = opt_map_anomalous(istate,jstate)
    !       if(imap.gt.0) lm(1,istate,jstate)=lm_cmplx(imap)
    !       if(jmap.gt.0) lm(2,istate,jstate)=lm_cmplx(jmap+Nopt_normal)
    !    end do
    ! end do
    !

    proj_variational_density=0.d0
    !+- build up the local H_projectors -+!
    call build_H_GZproj_superc(H_projectors,slater_derivatives,n0_target,lm,ifree_)
    !  
    !<FUCKING HARD DEBUGGING
    ! write(250,*) lm_(1)
    ! do iphi=1,Nphi
    !    do jphi=1,Nphi
    !       write(250,*) H_projectors(iphi,jphi),iphi,jphi
    !    end do
    !    write(250,*)
    ! end do
    !FUCKING HARD DEBUGGING>
    call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
    !<TEST_DEGENERATE
    ! tmp_eigen=H_eigens(1)
    ! do iphi=2,Nphi
    !    if(abs(tmp_eigen-H_eigens(iphi)).gt.1.d-10) exit       
    !    tmp_eigen=H_eigens(iphi)
    ! end do
    ! Ndegen=iphi-1
    ! proj_gs=zero
    ! do iphi=1,1
    !    proj_gs = proj_gs + H_projectors(:,iphi)/sqrt(dble(Ndegen))
    ! end do
    ! write(*,*) 'Ndegen!!',Ndegen
    !TEST_DEGENERATE>
    proj_gs = H_projectors(:,1)
    !
    do is=1,Ns
       do js=1,Ns
          !
          proj_variational_density(1,is,js) = trace_phi_basis(proj_gs,phi_traces_basis_dens(is,js,:,:))
          proj_variational_density(2,is,js) = trace_phi_basis(proj_gs,phi_traces_basis_dens_anomalous(is,js,:,:))
          !
          !          write(*,*) proj_variational_density(1,is,js)
       end do
    end do
    delta_proj_variational_density = proj_variational_density
    do istate=1,Ns
       delta_proj_variational_density(1,istate,istate) = & 
            delta_proj_variational_density(1,istate,istate) - n0_target(istate) 
    end do
    !    
    delta=0.d0
    allocate(delta_cmplx(Nvdm_NC_opt))
    call vdm_NC_stride_m2v(delta_proj_variational_density(1,:,:),delta_cmplx)
    do i=1,Nvdm_NC_opt
       delta(i) = dreal(delta_cmplx(i))
       delta(i+Nvdm_NC_opt) = dimag(delta_cmplx(i))
    end do
    deallocate(delta_cmplx)
    i0 = 2*Nvdm_NC_opt
    allocate(delta_cmplx(Nvdm_AC_opt))
    call vdm_AC_stride_m2v(delta_proj_variational_density(2,:,:),delta_cmplx)
    do i=1,Nvdm_AC_opt
       delta(i0+i) = dreal(delta_cmplx(i))
       delta(i0+i+Nvdm_AC_opt) = dimag(delta_cmplx(i))       
    end do
    deallocate(delta_cmplx)
    !
  end function fix_density
  !
end subroutine gz_projectors_minimization_nlep_superc




subroutine gz_proj_minimization_fixed_lgr_superc(n0,slater_derivatives,lgr_multip,E_Hloc,GZvect,free_flag) 
  real(8),dimension(Ns)           :: n0
  complex(8),dimension(2,Ns,Ns) :: slater_derivatives
  complex(8),dimension(2,Ns,Ns) :: lgr_multip
  real(8)                                :: E_Hloc
  complex(8),dimension(Nphi)               :: GZvect
  logical,optional :: free_flag
  logical :: free_flag_
  complex(8),dimension(Nphi,Nphi)         :: H_projectors
  real(8),dimension(Nphi)               :: H_eigens  
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                                :: iphi,jphi,is

  real(8),dimension(Ns,Ns) :: proj_variational_density
  !
  !+- build up the local H_projectors -+!
  free_flag_=.false.; if(present(free_flag)) free_flag_=free_flag
  !
  call build_H_GZproj_superc(H_projectors,slater_derivatives,n0,lgr_multip)
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !
  GZvect=H_projectors(1:Nphi,1)
  !
  E_Hloc=trace_phi_basis(GZvect,phi_traces_basis_Hloc)
  !
  !<DEBUG
  ! write(*,*) 'tmp_test GZproj_variational_density'
  ! do istate=1,Ns
  !    do jstate=1,Ns
  !       !
  !       proj_variational_density(istate,jstate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
  !       !
  !    end do
  !    write(*,*) proj_variational_density(istate,:)
  ! end do
  ! write(*,*)
  ! do istate=1,Nphi
  !    write(*,*) H_eigens(istate),H_projectors(istate,1)
  ! end do
  ! do is=1,Ns
  !    write(*,*) slater_derivatives(1,is,:)
  ! end do
  !DEBUG>
  !
end subroutine gz_proj_minimization_fixed_lgr_superc


!subroutine get_GZproj_ground_state_fixR(n0,lgr_multip_dens,lgr_multip_Rhop,E_Hloc,GZvect) 

subroutine gz_proj_minimization_fixed_lgr_superc_(n0,lgr_multip,lgr_multip_Rhop,lgr_multip_Qhop,E_Hloc,GZvect,free_flag) 
  real(8),dimension(Ns)           :: n0
  complex(8),dimension(2,Ns,Ns)   :: lgr_multip
  complex(8),dimension(Ns,Ns)     :: lgr_multip_Rhop,lgr_multip_Qhop
  real(8)                         :: E_Hloc
  complex(8),dimension(Nphi)      :: GZvect
  logical,optional                :: free_flag
  logical                         :: free_flag_
  complex(8),dimension(Nphi,Nphi) :: H_projectors
  real(8),dimension(Nphi)         :: H_eigens  
  integer                         :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                         :: iphi,jphi,is,js
  real(8),dimension(Ns,Ns)        :: proj_variational_density
  !
  !+- build up the local H_projectors -+!
  free_flag_=.false.; if(present(free_flag)) free_flag_=free_flag
  !
  H_projectors=zero
  H_projectors = phi_traces_basis_free_Hloc
  if(.not.free_flag_) H_projectors = phi_traces_basis_Hloc
  do is=1,Ns
     do js=1,Ns
        ! !        
        H_projectors = H_projectors + dreal(lgr_multip(1,is,js))*phi_traces_basis_dens(is,js,:,:)
        H_projectors = H_projectors + dreal(lgr_multip(1,is,js))*phi_traces_basis_dens_hc(is,js,:,:)
        ! !
        H_projectors = H_projectors + lgr_multip(2,is,js)*phi_traces_basis_dens_anomalous(is,js,:,:)
        H_projectors = H_projectors + conjg(lgr_multip(2,is,js))*phi_traces_basis_dens_anomalous_hc(is,js,:,:)
        ! !
        H_projectors = H_projectors - lgr_multip_Rhop(is,js)*phi_traces_basis_Rhop(is,js,:,:)
        H_projectors = H_projectors - conjg(lgr_multip_Rhop(is,js))*phi_traces_basis_Rhop_hc(is,js,:,:)
        ! !
        H_projectors = H_projectors - lgr_multip_Qhop(is,js)*phi_traces_basis_Qhop(is,js,:,:)
        H_projectors = H_projectors - conjg(lgr_multip_Qhop(is,js))*phi_traces_basis_Qhop_hc(is,js,:,:)
        ! !
     end do
  end do
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')
  !
  GZvect=H_projectors(1:Nphi,1)
  E_Hloc=trace_phi_basis(GZvect,phi_traces_basis_Hloc)
  !
end subroutine gz_proj_minimization_fixed_lgr_superc_






subroutine build_H_GZproj(H_projectors,slater_derivatives,n0,lgr_multip,ifree)
  complex(8),dimension(Nphi,Nphi) :: H_projectors
  complex(8),dimension(Ns,Ns) :: slater_derivatives
  real(8),dimension(Ns) :: n0
  real(8),dimension(Ns,Ns) :: lgr_multip
  logical,optional :: ifree
  logical :: ifree_
  integer :: is,js
  !
  ifree_=.false.;if(present(ifree)) ifree_=ifree  
  H_projectors=zero
  H_projectors=phi_traces_basis_Hloc
  if(ifree_) H_projectors = phi_traces_basis_free_Hloc  
  do is=1,Ns
     do js=1,Ns
        H_projectors = H_projectors + slater_derivatives(is,js)*phi_traces_basis_Rhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))          
        H_projectors = H_projectors + conjg(slater_derivatives(is,js))*phi_traces_basis_Rhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))          
        H_projectors = H_projectors + lgr_multip(is,js)*phi_traces_basis_dens(is,js,:,:)
        H_projectors = H_projectors + lgr_multip(is,js)*phi_traces_basis_dens_hc(is,js,:,:)
     end do
  end do
  !
end subroutine build_H_GZproj



subroutine build_H_GZproj_(H_projectors,slater_derivatives,n0,lgr_multip,ifree)
  complex(8),dimension(Nphi,Nphi) :: H_projectors
  complex(8),dimension(Ns,Ns) :: slater_derivatives
  real(8),dimension(Ns) :: n0
  complex(8),dimension(Ns,Ns) :: lgr_multip
  logical,optional :: ifree
  logical :: ifree_
  integer :: is,js
  !
  ifree_=.false.;if(present(ifree)) ifree_=ifree  
  H_projectors=zero
  H_projectors=phi_traces_basis_Hloc
  if(ifree_) H_projectors = phi_traces_basis_free_Hloc  
!  write(*,*) 'TMP_SLATER_DERIV'
  do is=1,Ns
     do js=1,Ns
        H_projectors = H_projectors + slater_derivatives(is,js)*phi_traces_basis_Rhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))          
        H_projectors = H_projectors + conjg(slater_derivatives(is,js))*phi_traces_basis_Rhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))          
        H_projectors = H_projectors + lgr_multip(is,js)*phi_traces_basis_dens(is,js,:,:)
        H_projectors = H_projectors + conjg(lgr_multip(is,js))*phi_traces_basis_dens_hc(is,js,:,:)
     end do
!     write(*,*) slater_derivatives(is,:)
  end do
  !
end subroutine build_H_GZproj_


subroutine build_H_GZproj_superc(H_projectors,slater_derivatives,n0,lgr_multip,ifree)
  complex(8),dimension(Nphi,Nphi) :: H_projectors
  complex(8),dimension(2,Ns,Ns) :: slater_derivatives
  real(8),dimension(Ns) :: n0
  complex(8),dimension(2,Ns,Ns) :: lgr_multip
  logical,optional :: ifree
  logical :: ifree_
  integer :: is,js
  !
  ifree_=.false.;if(present(ifree)) ifree_=ifree  
  H_projectors=zero
  H_projectors=phi_traces_basis_Hloc
  if(ifree_) H_projectors = phi_traces_basis_free_Hloc  
  do is=1,Ns
     do js=1,Ns
        !
        H_projectors = H_projectors + slater_derivatives(1,is,js)*phi_traces_basis_Rhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))          
        H_projectors = H_projectors + conjg(slater_derivatives(1,is,js))*phi_traces_basis_Rhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))          
        !
        H_projectors = H_projectors + slater_derivatives(2,is,js)*phi_traces_basis_Qhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))          
        H_projectors = H_projectors + conjg(slater_derivatives(2,is,js))*phi_traces_basis_Qhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))          
        ! !        
        H_projectors = H_projectors + lgr_multip(1,is,js)*phi_traces_basis_dens(is,js,:,:)
        H_projectors = H_projectors + conjg(lgr_multip(1,is,js))*phi_traces_basis_dens_hc(is,js,:,:)
        ! !
        H_projectors = H_projectors + lgr_multip(2,is,js)*phi_traces_basis_dens_anomalous(is,js,:,:)
        H_projectors = H_projectors + conjg(lgr_multip(2,is,js))*phi_traces_basis_dens_anomalous_hc(is,js,:,:)
        !
     end do
  end do
  !
end subroutine build_H_GZproj_superc




