subroutine slater_minimization_lgr(Rhop,n0_target,Estar,lgr_multip,n0_out,slater_derivatives,slater_matrix_el,iverbose) 
  complex(8),dimension(Ns,Ns),intent(in)  :: Rhop         !input:  renrmalization matrix
  real(8),dimension(Ns),intent(in)        :: n0_target    !input:  variational density matrix
  real(8),intent(out)                     :: Estar        !output: Slater Deter GS energy
  real(8),dimension(Ns,Ns),intent(out)    :: lgr_multip   !output: Slater Deter lagrange multipliers
  real(8),dimension(Ns,Ns),optional :: n0_out
  complex(8),dimension(Ns,Ns),optional    :: slater_derivatives
  complex(8),dimension(Ns,Ns,Lk),optional :: slater_matrix_el
  !
  logical,optional                     :: iverbose     !input:  Verbosity level
  !
  complex(8),dimension(Ns,Ns)    :: slater_derivatives_
  complex(8),dimension(Ns,Ns,Lk) :: slater_matrix_el_
  real(8),dimension(Ns,Ns) :: n0_out_
  !
  real(8),dimension(:),allocatable                :: lgr,delta_out     !+- real indeendent lgr_vector -+!


  complex(8),dimension(Ns,Ns)             :: Hk
  real(8),dimension(Ns)                :: tmp_lgr,err_dens
  integer                              :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb,iter,imap,is,js
  real(8),dimension(Ns)                :: lgr_multip_vec
  logical                              :: iverbose_
  real(8)  :: delta
  !
  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose
  !    
  allocate(lgr(Nopt_diag+Nopt_odiag));lgr=0.d0
  allocate(delta_out(Nopt_diag+Nopt_odiag))

  !+- HERE CHANGE Nopt = 2*Nvdm_NC_opt; allocate(Nvdm_NC_opt)

  !
  lgr = 0.d0!lgr_init_slater
  select case(lgr_method)
  case('CG_min')
     call fmin_cg(lgr,get_delta_local_density_matrix,iter,delta,itmax=20)
  case('f_zero')
     call fsolve(fix_density,lgr,tol=1.d-10,info=iter)
     delta_out=fix_density(lgr)
     delta=0.d0
     do is=1,Nopt_diag+Nopt_odiag
        delta = delta + delta_out(is)**2.d0
     end do
  end select
  lgr_init_slater=lgr
  !
  lgr_multip=0.d0
  do istate=1,Ns
     do jstate=1,Ns
        imap = opt_map(istate,jstate)
        if(imap.gt.0) lgr_multip(istate,jstate)=lgr(imap)
     end do
  end do
  !
  !+- HERE dump ----> call stride_v2m(lgr,lgr_multip)
  !
  call slater_minimization_fixed_lgr(Rhop,lgr_multip,Estar,n0_out_,slater_derivatives_,slater_matrix_el_)
  !
  if(present(n0_out)) n0_out=n0_out_
  if(present(slater_derivatives)) slater_derivatives = slater_derivatives_
  if(present(slater_matrix_el)) slater_matrix_el=slater_matrix_el_
  !
  if(iverbose_) then
     write(*,*)
     write(*,*) "Slater Determinant: input Rhop"
     do is=1,Ns
        write(*,'(20F7.3)') dreal(Rhop(is,:)),dimag(Rhop(is,:))
     end do
     write(*,*)
     write(*,*) "Slater Determinant: Lagrange Multipliers"
     do is=1,Ns
        write(*,'(20F18.10)') lgr_multip(is,:)
     end do
     write(*,*) "Slater Determinant: Variational density error"
     write(*,"(20F18.10)") delta
     write(*,*) "Slater Determinant: Ground State energy"
     write(*,"(20F18.10)") Estar
     write(*,*)
  end if
  !
contains    
  !
  function get_delta_local_density_matrix(lm_) result(delta)
    real(8),dimension(:)   :: lm_
    real(8)                :: delta
    real(8),dimension(Ns)  :: delta_local_density_matrix_vec
    real(8),dimension(Ns*Ns)  :: delta_local_density_matrix_vec_
    real(8),dimension(Ns,Ns)  :: lm
    real(8),dimension(Ns,Ns)  :: delta_local_density_matrix,local_density_matrix
    real(8),dimension(Ns,Ns) :: Hk,tmp
    real(8),dimension(Ns)          :: ek
    integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,imap
    !
    lm=0.d0
    do istate=1,Ns
       do jstate=1,Ns
          imap = opt_map(istate,jstate)
          if(imap.gt.0) lm(istate,jstate)=lm_(imap)
       end do
    end do
    !+- HERE dump ---> stride_v2m  (...in turn change to complex lgr-parameters) -+!
    !
    local_density_matrix=0.d0
    !
    do ik=1,Lk     
       Hk=0.d0
       ek=0.d0
       ! hopping renormalization !
       Hk=matmul(Hk_tb(:,:,ik),Rhop)
       Hk=matmul(Rhop,Hk)
       ! add Lagrange multipliers !
       Hk=Hk+lm                     
       ! diagonalize hamiltonian !
       call  matrix_diagonalize(Hk,ek,'V','L')
       !compute local density matrix
       do istate=1,Ns
          do jstate=1,Ns
             do kstate=1,Ns
                !
                local_density_matrix(istate,jstate) = &
                     local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*Hk(istate,kstate)*Hk(jstate,kstate)*wtk(ik)
                !
             end do
          end do
       end do
    end do
    ! return variation of local density matrix with respect to the target values
    delta_local_density_matrix = local_density_matrix
    do istate=1,Ns
       delta_local_density_matrix(istate,istate) = delta_local_density_matrix(istate,istate) - n0_target(istate)      
    end do
    delta=0.d0
    do istate=1,Ns
       do jstate=1,Ns
          delta = delta + abs(delta_local_density_matrix(istate,jstate))**2.d0
       end do
    end do
    !
    write(*,*) lm_
    write(*,*) 'delta',delta
  end function get_delta_local_density_matrix



  function fix_density(lm_) result(delta)
    real(8),dimension(:)   :: lm_
    real(8),dimension(size(lm_))      :: delta
    real(8),dimension(Ns)  :: delta_local_density_matrix_vec
    real(8),dimension(Ns*Ns)  :: delta_local_density_matrix_vec_
    real(8),dimension(Ns,Ns)  :: lm
    real(8),dimension(Ns,Ns)  :: delta_local_density_matrix,local_density_matrix
    real(8),dimension(Ns,Ns) :: Hk,tmp
    real(8),dimension(Ns)          :: ek
    integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,imap
    !
    complex(8),dimension(Ns,Ns) :: Rhop_dag
    lm=0.d0
    do istate=1,Ns
       do jstate=1,Ns
          imap = opt_map(istate,jstate)
          if(imap.gt.0) lm(istate,jstate)=lm_(imap)
          Rhop_dag(istate,jstate) = conjg(Rhop(jstate,istate))
       end do
    end do
    !+- HERE use strides instead of opt_map
    !
    local_density_matrix=0.d0
    !
    do ik=1,Lk
       Hk=0.d0
       ek=0.d0
       ! hopping renormalization !
       Hk=matmul(Hk_tb(:,:,ik),Rhop)
       Hk=matmul(Rhop_dag,Hk)
       ! add Lagrange multipliers !
       Hk=Hk+lm                     
       ! diagonalize hamiltonian !
       call  matrix_diagonalize(Hk,ek,'V','L')
       !compute local density matrix
       do istate=1,Ns
          do jstate=1,Ns
             do kstate=1,Ns
                !
                local_density_matrix(istate,jstate) = &
                     local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*Hk(istate,kstate)*Hk(jstate,kstate)*wtk(ik)
                !
             end do
          end do
       end do
    end do
    ! return variation of local density matrix with respect to the target values
    delta_local_density_matrix = local_density_matrix
    do istate=1,Ns
       delta_local_density_matrix(istate,istate) = delta_local_density_matrix(istate,istate) - n0_target(istate)      
    end do
    !
    !+- here use stride and AND, in turn, insert a check for the stride!!!
    !
    delta= 0.d0
    do istate=1,Ns
       do jstate=1,Ns
          imap = opt_map(istate,jstate)
          if(imap.gt.0) then
             delta(imap) = delta(imap) + delta_local_density_matrix(istate,jstate)**2.d0
             ! if(istate.eq.jstate) then
             !    delta(imap) = local_density_matrix(istate,jstate) - &
             !         n0_target(istate)
             ! else
             !    delta(imap) = local_density_matrix(istate,jstate)
             ! end if
          end if
       end do
    end do
  end function fix_density
end subroutine slater_minimization_lgr
!
subroutine slater_minimization_fixed_lgr(Rhop,lm,Estar,n0,slater_derivatives,slater_matrix_el)
  complex(8),dimension(Ns,Ns),intent(in)  :: Rhop
  real(8),dimension(Ns,Ns),intent(in)     :: lm
  real(8)                                 :: Estar
  ! optional inputs
  real(8),dimension(Ns,Ns),optional       :: n0
  complex(8),dimension(Ns,Ns),optional    :: slater_derivatives
  complex(8),dimension(Ns,Ns,Lk),optional :: slater_matrix_el
  real(8),dimension(Ns,Ns)                :: n0_
  complex(8),dimension(Ns,Ns)             :: slater_derivatives_
  complex(8),dimension(Ns,Ns,Lk)          :: slater_matrix_el_
  complex(8),dimension(Ns,Ns)             :: Hk,tmp,Hk_bare,Hstar
  real(8),dimension(Ns)                   :: ek
  integer                                 :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,is,js,ks,kks
  !
  complex(8),dimension(Ns,Ns)             :: Rhop_dag  
  !
  Estar=0.d0
  slater_derivatives_=zero
  n0_ = 0.d0
  slater_matrix_el_ = zero
  !
  do istate=1,Ns
     do jstate=1,Ns
        Rhop_dag(istate,jstate) = conjg(Rhop(jstate,istate))
     end do
  end do

  do ik=1,Lk
     !
     Hk=0.d0
     ek=0.d0
     Hk_bare=Hk_tb(:,:,ik)
     ! hopping renormalization !
     Hk=matmul(Hk_tb(:,:,ik),Rhop)
     Hk=matmul(Rhop_dag,Hk)  
     Hstar=Hk
     Hk=Hk+lm
     !
     !<TMP
     ! do is=1,Ns
     !    do js=1,Ns
     !       Hk(is,js) = Hk(is,js) + lm(js,is)
     !    end do
     ! end do
     !TMP>
     !
     ! diagonalize hamiltonian !
     call  matrix_diagonalize(Hk,ek,'V','L')
     ! store slater determinant matrix elements
     do iorb=1,Norb
        do ispin=1,2
           do jorb=1,Norb
              do jspin=1,2
                 istate=index(ispin,iorb)
                 jstate=index(jspin,jorb)               
                 tmp(istate,jstate)=0.d0
                 do kstate=1,Ns
                    tmp(istate,jstate) = tmp(istate,jstate) + conjg(Hk(istate,kstate))*fermi(ek(kstate),beta)*Hk(jstate,kstate)
                    
                    Estar = Estar + Hstar(istate,jstate)*conjg(Hk(istate,kstate))*Hk(jstate,kstate)*fermi(ek(kstate),beta)*wtk(ik)
                    n0_(istate,jstate) = n0_(istate,jstate) + conjg(Hk(istate,kstate))*fermi(ek(kstate),beta)*Hk(jstate,kstate)*wtk(ik)                    
                    slater_matrix_el_(istate,jstate,ik) = slater_matrix_el_(istate,jstate,ik) + &
                         conjg(Hk(istate,kstate))*Hk(jstate,kstate)*fermi(ek(kstate),beta)
                 end do
              end do
           end do
        end do
     end do


     tmp = slater_matrix_el_(:,:,ik)     
     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns
                 !
                 slater_derivatives_(is,js) = slater_derivatives_(is,js) + conjg(Rhop(kks,ks))*Hk_bare(kks,is)*tmp(ks,js)*wtk(ik)
                 
              end do
           end do
        end do
     end do

     ! do istate=1,Ns
     !    do jstate=1,Ns
     !       !+-  store slater ground state derivatives -+!
     !       do is=1,Ns
     !          do js=1,Ns
     !             do ks=1,Ns
     !                slater_derivatives_(istate,jstate) = slater_derivatives_(istate,jstate) + &
     !                     Hk(jstate,is)*fermi(ek(is),beta)*conjg(Hk(js,is))*Rhop_dag(js,ks)*Hk_bare(ks,istate)*wtk(ik)
     !             end do
     !          end do
     !       end do           
     !    end do
     ! end do
     !
     ! tmp=matmul(Rhop_dag,tmp)
     ! tmp=matmul(Hk_bare,tmp)             
     ! do istate=1,Ns
     !    do jstate=1,Ns             
     !       slater_derivatives_(istate,jstate) = &
     !            slater_derivatives_(istate,jstate) + 2.d0*tmp(istate,jstate)*wtk(ik)
     !    end do
     ! end do
     
  end do
  if(present(n0)) then
     n0=n0_
  end if
  if(present(slater_derivatives)) then
     slater_derivatives=slater_derivatives_
  end if
  if(present(slater_matrix_el)) then
     slater_matrix_el=slater_matrix_el_
  end if

  !<DEBUG
  write(*,*) 
  write(*,*) 'tmp test'
  do istate=1,Ns
     write(*,*) slater_derivatives_(istate,:)
  end do
  !DEBUG>

  !
end subroutine slater_minimization_fixed_lgr






















subroutine slater_minimization_lgr_(Rhop,n0_target,Estar,lgr_multip,n0_out,slater_derivatives,slater_matrix_el,iverbose) 
  complex(8),dimension(Ns,Ns),intent(in)  :: Rhop         !input:  renrmalization matrix
  real(8),dimension(Ns),intent(in)        :: n0_target    !input:  variational density matrix
  real(8),intent(out)                     :: Estar        !output: Slater Deter GS energy
  complex(8),dimension(Ns,Ns),intent(out)    :: lgr_multip   !output: Slater Deter lagrange multipliers
  complex(8),dimension(Ns,Ns),optional :: n0_out
  complex(8),dimension(Ns,Ns),optional    :: slater_derivatives
  complex(8),dimension(Ns,Ns,Lk),optional :: slater_matrix_el
  !
  logical,optional                     :: iverbose     !input:  Verbosity level
  !
  complex(8),dimension(Ns,Ns)    :: slater_derivatives_
  complex(8),dimension(Ns,Ns,Lk) :: slater_matrix_el_
  complex(8),dimension(Ns,Ns) :: n0_out_
  !
  real(8),dimension(:),allocatable                :: lgr,delta_out     !+- real indeendent lgr_vector -+!
  complex(8),dimension(:),allocatable   :: lgr_cmplx
  !
  complex(8),dimension(Ns,Ns)             :: Hk
  real(8),dimension(Ns)                :: tmp_lgr,err_dens
  integer                              :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb,iter,imap,is,js,i,Nopt
  real(8),dimension(Ns)                :: lgr_multip_vec
  logical                              :: iverbose_
  real(8)  :: delta
  !
  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose
  !
  Nopt=2*Nvdm_NC_opt; allocate(lgr(Nopt));allocate(delta_out(Nopt))
  !
  lgr=0.d0  
  select case(lgr_method)
  case('CG_min')
     call fmin_cg(lgr,get_delta_local_density_matrix,iter,delta,itmax=20)
  case('f_zero')
     call fsolve(fix_density,lgr,tol=1.d-10,info=iter)
     delta_out=fix_density(lgr)
     delta=0.d0
     do is=1,Nopt
        delta = delta + delta_out(is)**2.d0
     end do
  end select
  !
  allocate(lgr_cmplx(Nvdm_NC_opt))
  do i=1,Nvdm_NC_opt
     lgr_cmplx(i) = lgr(i)+xi*lgr(i+Nvdm_NC_opt)
  end do
  call vdm_NC_stride_v2m(lgr_cmplx,lgr_multip)
  deallocate(lgr_cmplx)
  !
  call slater_minimization_fixed_lgr_(Rhop,lgr_multip,Estar,n0_out_,slater_derivatives_,slater_matrix_el_)
  !
  if(present(n0_out)) n0_out=n0_out_
  if(present(slater_derivatives)) slater_derivatives = slater_derivatives_
  if(present(slater_matrix_el)) slater_matrix_el=slater_matrix_el_
  !
  if(iverbose_) then
     write(*,*)
     write(*,*) "Slater Determinant: input Rhop"
     do is=1,Ns
        write(*,'(20F7.3)') dreal(Rhop(is,:)),dimag(Rhop(is,:))
     end do
     write(*,*)
     write(*,*) "Slater Determinant: Lagrange Multipliers"
     do is=1,Ns
        write(*,'(20F18.10)') lgr_multip(is,:)
     end do
     write(*,*) "Slater Determinant: Variational density error"
     write(*,"(20F18.10)") delta
     write(*,*) "Slater Determinant: Ground State energy"
     write(*,"(20F18.10)") Estar
     write(*,*)
  end if
  !
contains    
  !
  function get_delta_local_density_matrix(lm_) result(delta)
    real(8),dimension(:)                :: lm_
    real(8)                             :: delta
    complex(8),dimension(Ns,Ns)         :: lm
    complex(8),dimension(Ns,Ns)         :: delta_local_density_matrix,local_density_matrix
    complex(8),dimension(Ns,Ns)         :: Hk,tmp
    complex(8),dimension(:),allocatable :: lm_cmplx
    real(8),dimension(Ns)               :: ek
    integer                             :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,imap
    complex(8),dimension(Ns,Ns)         :: Rhop_dag
    !
    allocate(lm_cmplx(Nvdm_NC_opt))
    do i=1,Nvdm_NC_opt
       lm_cmplx(i) = lm_(i)+xi*lm_(i+Nvdm_NC_opt)
    end do
    call vdm_NC_stride_v2m(lm_cmplx,lm)    
    deallocate(lm_cmplx)
    !
    do is=1,Ns
       do js=1,Ns
          Rhop_dag(is,js) = conjg(Rhop(js,is))
       end do
    end do
    !
    local_density_matrix=0.d0
    !
    do ik=1,Lk     
       Hk=0.d0
       ek=0.d0
       ! hopping renormalization !
       Hk=matmul(Hk_tb(:,:,ik),Rhop)
       Hk=matmul(Rhop_dag,Hk)
       ! add Lagrange multipliers !
       Hk=Hk+lm                     
       ! diagonalize hamiltonian !
       call  matrix_diagonalize(Hk,ek)
       !compute local density matrix
       do istate=1,Ns
          do jstate=1,Ns
             do kstate=1,Ns
                !
                local_density_matrix(istate,jstate) = &
                     local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*conjg(Hk(istate,kstate))*Hk(jstate,kstate)*wtk(ik)
                !
             end do
          end do
       end do
    end do
    ! return variation of local density matrix with respect to the target values
    delta_local_density_matrix = local_density_matrix
    do istate=1,Ns
       delta_local_density_matrix(istate,istate) = delta_local_density_matrix(istate,istate) - n0_target(istate)      
    end do
    delta=0.d0
    do istate=1,Ns
       do jstate=1,Ns
          delta = delta + abs(delta_local_density_matrix(istate,jstate))**2.d0
       end do
    end do
    !
    write(*,*) lm_
    write(*,*) 'delta No_SLATER',delta
  end function get_delta_local_density_matrix
  !
  function fix_density(lm_) result(delta)
    real(8),dimension(:)   :: lm_
    real(8),dimension(size(lm_))      :: delta
    complex(8),dimension(Ns,Ns)  :: lm
    complex(8),dimension(Ns,Ns)  :: delta_local_density_matrix,local_density_matrix
    complex(8),dimension(Ns,Ns)  :: delta_local_density_check
    complex(8),dimension(:),allocatable :: lm_cmplx,delta_cmplx
    complex(8),dimension(Ns,Ns) :: Hk,tmp
    real(8),dimension(Ns)          :: ek
    real(8)                        :: check
    integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,imap
    !
    complex(8),dimension(Ns,Ns) :: Rhop_dag
    !
    allocate(lm_cmplx(Nvdm_NC_opt))
    do i=1,Nvdm_NC_opt
       lm_cmplx(i) = lm_(i)+xi*lm_(i+Nvdm_NC_opt)
    end do
    call vdm_NC_stride_v2m(lm_cmplx,lm)    
    deallocate(lm_cmplx)
    !
    do istate=1,Ns
       do jstate=1,Ns
          Rhop_dag(istate,jstate) = conjg(Rhop(jstate,istate))
       end do
    end do
    !
    local_density_matrix=0.d0
    do ik=1,Lk
       Hk=0.d0
       ek=0.d0
       ! hopping renormalization !
       Hk=matmul(Hk_tb(:,:,ik),Rhop)
       Hk=matmul(Rhop_dag,Hk)
       ! add Lagrange multipliers !
       Hk=Hk+lm                     
       ! diagonalize hamiltonian !
       call  matrix_diagonalize(Hk,ek)
       !compute local density matrix
       do istate=1,Ns
          do jstate=1,Ns
             do kstate=1,Ns
                local_density_matrix(istate,jstate) = &
                     local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*conjg(Hk(istate,kstate))*Hk(jstate,kstate)*wtk(ik)
             end do
          end do
       end do
    end do
    ! return variation of local density matrix with respect to the target values
    delta_local_density_matrix = local_density_matrix
    do istate=1,Ns
       delta_local_density_matrix(istate,istate) = delta_local_density_matrix(istate,istate) - n0_target(istate)      
    end do
    !
    delta=0.d0
    allocate(delta_cmplx(Nvdm_NC_opt))
    call vdm_NC_stride_m2v(delta_local_density_matrix,delta_cmplx)
    do i=1,Nvdm_NC_opt
       delta(i) = dreal(delta_cmplx(i))
       delta(i+Nvdm_NC_opt) = dimag(delta_cmplx(i))
    end do
    !+- check if the stride is compatible with the results -+!
    call vdm_NC_stride_v2m(delta_cmplx,delta_local_density_check)
    check = 0.d0
    do is=1,Ns
       do js=1,Ns
          check = check + & 
               (delta_local_density_check(is,js)-delta_local_density_matrix(is,js))*conjg(delta_local_density_check(is,js)-delta_local_density_matrix(is,js))
       end do
    end do
    if(check.gt.1.d-10) stop "CHECK STRIDES @ fix_density_normal"
    write(*,*) delta
  end function fix_density
end subroutine slater_minimization_lgr_
!
subroutine slater_minimization_fixed_lgr_(Rhop,lm,Estar,n0,slater_derivatives,slater_matrix_el)
  complex(8),dimension(Ns,Ns),intent(in)  :: Rhop
  complex(8),dimension(Ns,Ns),intent(in)  :: lm
  real(8)                                 :: Estar
  ! optional inputs
  complex(8),dimension(Ns,Ns),optional    :: n0 
  complex(8),dimension(Ns,Ns),optional    :: slater_derivatives
  complex(8),dimension(Ns,Ns,Lk),optional :: slater_matrix_el
  complex(8),dimension(Ns,Ns)                :: n0_
  complex(8),dimension(Ns,Ns)             :: slater_derivatives_
  complex(8),dimension(Ns,Ns,Lk)          :: slater_matrix_el_
  complex(8),dimension(Ns,Ns)             :: Hk,tmp,Hk_bare,Hstar
  real(8),dimension(Ns)                   :: ek
  integer                                 :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,is,js,ks,kks
  !
  complex(8),dimension(Ns,Ns)             :: Rhop_dag  
  !
  do istate=1,Ns
     do jstate=1,Ns
        Rhop_dag(istate,jstate) = conjg(Rhop(jstate,istate))
     end do
  end do
  !
  Estar=0.d0
  slater_derivatives_=zero
  n0_ = 0.d0
  slater_matrix_el_ = zero
  do ik=1,Lk
     !
     Hk=0.d0
     ek=0.d0
     Hk_bare=Hk_tb(:,:,ik)
     ! hopping renormalization !
     Hk=matmul(Hk_tb(:,:,ik),Rhop)
     Hk=matmul(Rhop_dag,Hk)  
     Hstar=Hk
     Hk=Hk+lm
     !
     call  matrix_diagonalize(Hk,ek)
     ! store slater determinant matrix elements
     do istate=1,Ns
        do jstate=1,Ns
           tmp(istate,jstate)=0.d0
           do kstate=1,Ns
              tmp(istate,jstate) = tmp(istate,jstate) + conjg(Hk(istate,kstate))*fermi(ek(kstate),beta)*Hk(jstate,kstate)
              Estar = Estar + Hstar(istate,jstate)*conjg(Hk(istate,kstate))*Hk(jstate,kstate)*fermi(ek(kstate),beta)*wtk(ik)
              n0_(istate,jstate) = n0_(istate,jstate) + conjg(Hk(istate,kstate))*fermi(ek(kstate),beta)*Hk(jstate,kstate)*wtk(ik)                    
              slater_matrix_el_(istate,jstate,ik) = slater_matrix_el_(istate,jstate,ik) + &
                   conjg(Hk(istate,kstate))*Hk(jstate,kstate)*fermi(ek(kstate),beta)
           end do
        end do
     end do
     !
     tmp = slater_matrix_el_(:,:,ik)     
     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns
                 slater_derivatives_(is,js) = slater_derivatives_(is,js) + conjg(Rhop(kks,ks))*Hk_bare(kks,is)*tmp(ks,js)*wtk(ik)
              end do
           end do
        end do
     end do
  end do
  if(present(n0)) then
     n0=n0_
  end if
  if(present(slater_derivatives)) then
     slater_derivatives=slater_derivatives_
  end if
  if(present(slater_matrix_el)) then
     slater_matrix_el=slater_matrix_el_
  end if
  !
end subroutine slater_minimization_fixed_lgr_





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





subroutine slater_minimization_lgr_superc(Rhop,Qhop,n0_target,Estar,lgr_multip,n0_out,slater_derivatives,slater_matrix_el,iverbose) 
  !
  complex(8),dimension(Ns,Ns),intent(in)      :: Rhop         !input:  normal    hopping renormalization matrix
  complex(8),dimension(Ns,Ns),intent(in)      :: Qhop         !input:  anomalous hopping renormalization matrix 
  !
  real(8),dimension(Ns),intent(in)            :: n0_target    !input:  variational density matrix in the natural basis
  real(8),intent(out)                         :: Estar        !output: Slater Deter GS energy
  complex(8),dimension(2,Ns,Ns),intent(out)   :: lgr_multip   !output: Slater Deter lagrange multipliers
  complex(8),dimension(2,Ns,Ns),optional      :: n0_out             !output: Slater Deter ground state density matrix
  complex(8),dimension(2,Ns,Ns),optional      :: slater_derivatives !output: Slater Deter ground state renormalization derivatives
  complex(8),dimension(2*Ns,2*Ns,Lk),optional :: slater_matrix_el   !output: Slater Deter ground state matrix elements
  !
  logical,optional                            :: iverbose     !input:  Verbosity level
  !
  complex(8),dimension(2,Ns,Ns)               :: slater_derivatives_
  complex(8),dimension(2*Ns,2*Ns,Lk)          :: slater_matrix_el_
  complex(8),dimension(2,Ns,Ns)               :: n0_out_
  !
  real(8),dimension(:),allocatable            :: lgr
  complex(8),dimension(:),allocatable         :: lgr_cmplx
  real(8),dimension(:),allocatable            ::   delta_out     !+- real indeendent lgr_vector -+!

  complex(8),dimension(Ns,Ns)                 :: Rhop_dag         
  complex(8),dimension(Ns,Ns)                 :: Qhop_dag         
  !
  complex(8),dimension(Ns,Ns)             :: Hk
  real(8),dimension(Ns)                :: tmp_lgr,err_dens
  integer                              :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb,iter,imap,jmap,is,js,i,i0
  real(8),dimension(Ns)                :: lgr_multip_vec
  logical                              :: iverbose_
  real(8)  :: delta
  integer :: Nopt
  !
  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose
  !    
  Nopt=2*Nvdm_NC_opt + 2*Nvdm_AC_opt;allocate(lgr(Nopt));allocate(delta_out(Nopt))
  !
  lgr=0.d0
  select case(lgr_method)
  case('CG_min')
     call fmin_cg(lgr,get_delta_local_density_matrix,iter,delta,itmax=20)
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
  call vdm_NC_stride_v2m(lgr_cmplx,lgr_multip(1,:,:))
  deallocate(lgr_cmplx)
  i0=2*Nvdm_NC_opt
  allocate(lgr_cmplx(Nvdm_AC_opt))  
  do i=1,Nvdm_AC_opt
     lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_AC_opt)
  end do
  call vdm_AC_stride_v2m(lgr_cmplx,lgr_multip(2,:,:))  
  !
  call slater_minimization_fixed_lgr_superc(Rhop,Qhop,lgr_multip,Estar,n0_out_,slater_derivatives_,slater_matrix_el_)
  !
  if(present(n0_out)) n0_out=n0_out_
  if(present(slater_derivatives)) slater_derivatives = slater_derivatives_
  if(present(slater_matrix_el)) slater_matrix_el=slater_matrix_el_
  !
  if(iverbose_) then
     write(*,*)
     write(*,*) "Slater Determinant: input Rhop"
     do is=1,Ns
        write(*,'(20F7.3)') dreal(Rhop(is,:)),dimag(Rhop(is,:))
     end do
     write(*,*)
     write(*,*)
     write(*,*) "Slater Determinant: input Qhop"
     do is=1,Ns
        write(*,'(20F7.3)') dreal(Qhop(is,:)),dimag(Qhop(is,:))
     end do
     write(*,*)

     write(*,*) "Slater Determinant: Lagrange Multipliers"
     do is=1,Ns
        write(*,'(20F7.3)') dreal(lgr_multip(1,is,:)),dimag(lgr_multip(1,is,:))
     end do
     write(*,*)
     do is=1,Ns
        write(*,'(20F7.3)') dreal(lgr_multip(2,is,:)),dimag(lgr_multip(2,is,:))
     end do
     write(*,*) "Slater Determinant: Variational density error"
     write(*,"(20F18.10)") delta
     write(*,*) "Slater Determinant: Ground State energy"
     write(*,"(20F18.10)") Estar
     write(*,*)
  end if
  !
contains    
  !

  function get_delta_local_density_matrix(lm_) result(delta)
    real(8),dimension(:)                :: lm_
    complex(8),dimension(:),allocatable :: lm_cmplx,delta_cmplx
    complex(8),dimension(2,Ns,Ns)       :: lm
    real(8)                             :: delta
    complex(8),dimension(2,Ns,Ns)       :: delta_local_density_matrix,local_density_matrix
    complex(8),dimension(2*Ns,2*Ns)     :: Hk
    complex(8),dimension(Ns,Ns)         :: Hk_tmp
    real(8),dimension(2*Ns)             :: ek
    real(8)                             :: nqp
    integer                             :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,imap,jmap,i,i0
    integer                             :: is,js,ks
    !
    ! !+- reconstruct complex lgr paramters -+!
    ! do i=1,Nopt_lgr
    !    lm_cmplx(i) = lm_(i) + xi*lm_(i+Nopt_lgr)
    ! end do
    ! !+- reconstruct full optimization parameter matrices -+!
    ! lm=0.d0
    ! do istate=1,Ns
    !    do jstate=1,Ns
    !       imap = opt_map(istate,jstate)
    !       if(imap.gt.0) lm(1,istate,jstate)=lm_cmplx(imap)
    !       jmap = opt_map_anomalous(istate,jstate)
    !       if(jmap.gt.0) then
    !          lm(2,istate,jstate) = lm_cmplx(jmap+Nopt_normal)          
    !          !lm(jstate+Ns,istate) = conjg(lm_cmplx(jmap+Nopt_normal))
    !       end if
    !    end do
    ! end do
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
    local_density_matrix=0.d0
    !
    do istate=1,Ns
       do jstate=1,Ns
          Rhop_dag(istate,jstate) = conjg(Rhop(jstate,istate))
          Qhop_dag(istate,jstate) = conjg(Qhop(jstate,istate))
       end do
    end do
    !
    do ik=1,Lk     
       Hk=0.d0
       ek=0.d0
       !
       ! hopping renormalization !
       !
       Hk_tmp=matmul(Hk_tb(:,:,ik),Rhop)
       Hk_tmp=matmul(Rhop_dag,Hk_tmp)
       Hk(1:Ns,1:Ns) = Hk_tmp 
       !
       Hk_tmp=matmul(Hk_tb(:,:,ik),Qhop)
       Hk_tmp=matmul(Rhop_dag,Hk_tmp)
       Hk(1:Ns,Ns+1:2*Ns) = Hk_tmp 
       !
       Hk_tmp=matmul(Hk_tb(:,:,ik),Rhop)
       Hk_tmp=matmul(Qhop_dag,Hk_tmp)
       Hk(Ns+1:2*Ns,1:Ns) = Hk_tmp
       !
       Hk_tmp=matmul(Hk_tb(:,:,ik),Qhop)
       Hk_tmp=matmul(Qhop_dag,Hk_tmp)
       Hk(Ns+1:2*Ns,Ns+1:2*Ns) = Hk_tmp
       !
       ! add Lagrange multipliers !
       Hk(1:Ns,1:Ns)=Hk(1:Ns,1:Ns)+lm(1,:,:)
       Hk(1:Ns,Ns+1:2*Ns)=Hk(1:Ns,1:Ns)+lm(2,:,:)
       do is=1,Ns
          do js=1,Ns
             Hk(is+Ns,js) = Hk(is+Ns,js) + conjg(lm(2,js,is))
          end do
       end do
       ! diagonalize hamiltonian !
       call  matrix_diagonalize(Hk,ek,'V','L')
       !compute local density matrix [Normal and Anomalous part]       
       do is=1,Ns
          do js=1,Ns
             !
             do ks=1,Ns
                nqp = fermi(ek(ks)-ek(ks+Ns),beta)                
                !
                local_density_matrix(1,is,js) = local_density_matrix(1,is,js) + &
                     conjg(Hk(is,ks))*Hk(js,ks)*nqp*wtk(ik) + conjg(Hk(is,ks+Ns))*Hk(js,ks+Ns)*(1.d0-nqp)*wtk(ik)
                !
                local_density_matrix(2,istate,jstate) = local_density_matrix(2,istate,jstate) + &
                     conjg(Hk(is+Ns,ks))*Hk(js,ks)*nqp*wtk(ik) + conjg(Hk(is+Ns,ks+Ns))*Hk(js,ks+Ns)*(1.d0-nqp)*wtk(ik)
             end do
             !
          end do
       end do

    end do
    !+- un attimo crucial point...
    ! return variation of local density matrix with respect to the target values
    delta_local_density_matrix = local_density_matrix
    do istate=1,Ns
       delta_local_density_matrix(1,istate,istate) = delta_local_density_matrix(1,istate,istate) - n0_target(istate)      
    end do
    !
    delta=0.d0
    do istate=1,Ns
       do jstate=1,Ns
          delta = delta + delta_local_density_matrix(1,istate,jstate)*conjg(delta_local_density_matrix(1,istate,jstate))
          delta = delta + delta_local_density_matrix(2,istate,jstate)*conjg(delta_local_density_matrix(2,istate,jstate))
       end do
    end do
    !
  end function get_delta_local_density_matrix
  !
  function fix_density(lm_) result(delta)
    real(8),dimension(:)                :: lm_
    real(8),dimension(size(lm_))        :: delta
    complex(8),dimension(:),allocatable :: lm_cmplx,delta_cmplx
    complex(8),dimension(2,Ns,Ns)       :: lm   !+- this may be also (2,Ns,Ns)
    complex(8),dimension(2,Ns,Ns)       :: delta_local_density_matrix,local_density_matrix
    complex(8),dimension(2*Ns,2*Ns)     :: Hk
    complex(8),dimension(Ns,Ns)         :: Hk_tmp
    real(8),dimension(2*Ns)             :: ek
    integer                             :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,imap,jmap,i,i0
    integer                             :: is,js,ks
    real(8)                             :: nqp
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
    ! do i=1,Nopt_lgr
    !    lm_cmplx(i) = lm_(i) + xi*lm_(i+Nopt_lgr)
    ! end do
    ! !+- reconstruct full optimization parameter matrices -+!
    ! ! lm=0.d0
    ! ! do istate=1,Ns
    ! !    do jstate=1,Ns
    ! !       imap = opt_map(istate,jstate)
    ! !       if(imap.gt.0) lm(istate,jstate)=lm_cmplx(imap)
    ! !       jmap = opt_map_anomalous(istate,jstate)
    ! !       if(jmap.gt.0) then
    ! !          lm(istate,jstate+Ns) = lm_cmplx(jmap+Nopt_normal)          
    ! !          lm(jstate+Ns,istate) = conjg(lm_cmplx(jmap+Nopt_normal))
    ! !       end if
    ! !    end do
    ! ! end do
    ! lm=0.d0
    ! do istate=1,Ns
    !    do jstate=1,Ns
    !       imap = opt_map(istate,jstate)
    !       if(imap.gt.0) lm(1,istate,jstate)=lm_cmplx(imap)
    !       jmap = opt_map_anomalous(istate,jstate)
    !       if(jmap.gt.0) then
    !          lm(2,istate,jstate) = lm_cmplx(jmap+Nopt_normal)          
    !       end if
    !    end do
    ! end do
    !
    !
    local_density_matrix=0.d0
    !
    do istate=1,Ns
       do jstate=1,Ns
          Rhop_dag(istate,jstate) = conjg(Rhop(jstate,istate))
          Qhop_dag(istate,jstate) = conjg(Qhop(jstate,istate))
       end do
    end do
    !
    do ik=1,Lk     
       Hk=0.d0
       ek=0.d0
       !
       ! hopping renormalization !
       !
       Hk_tmp=matmul(Hk_tb(:,:,ik),Rhop)
       Hk_tmp=matmul(Rhop_dag,Hk_tmp)
       Hk(1:Ns,1:Ns) = Hk_tmp 
       !
       Hk_tmp=matmul(Hk_tb(:,:,ik),Qhop)
       Hk_tmp=matmul(Rhop_dag,Hk_tmp)
       Hk(1:Ns,Ns+1:2*Ns) = Hk_tmp 
       !
       Hk_tmp=matmul(Hk_tb(:,:,ik),Rhop)
       Hk_tmp=matmul(Qhop_dag,Hk_tmp)
       Hk(Ns+1:2*Ns,1:Ns) = Hk_tmp
       !
       Hk_tmp=matmul(Hk_tb(:,:,ik),Qhop)
       Hk_tmp=matmul(Qhop_dag,Hk_tmp)
       Hk(Ns+1:2*Ns,Ns+1:2*Ns) = Hk_tmp
       !
       ! add Lagrange multipliers !
       Hk(1:Ns,1:Ns)=Hk(1:Ns,1:Ns)+lm(1,:,:)
       Hk(1:Ns,Ns+1:2*Ns)=Hk(1:Ns,1:Ns)+lm(2,:,:)
       do is=1,Ns
          do js=1,Ns
             Hk(is+Ns,js) = Hk(is+Ns,js) + conjg(lm(2,js,is))
          end do
       end do
       ! diagonalize hamiltonian !
       call  matrix_diagonalize(Hk,ek,'V','L')
       !compute local density matrix [Normal and Anomalous part]       
       do is=1,Ns
          do js=1,Ns
             !
             do ks=1,Ns
                nqp = fermi(ek(ks)-ek(ks+Ns),beta)                
                !
                local_density_matrix(1,is,js) = local_density_matrix(1,is,js) + &
                     conjg(Hk(is,ks))*Hk(js,ks)*nqp*wtk(ik) + conjg(Hk(is,ks+Ns))*Hk(js,ks+Ns)*(1.d0-nqp)*wtk(ik)
                !
                local_density_matrix(2,is,js) = local_density_matrix(2,is,js) + &
                     conjg(Hk(is+Ns,ks))*Hk(js,ks)*nqp*wtk(ik) + conjg(Hk(is+Ns,ks+Ns))*Hk(js,ks+Ns)*(1.d0-nqp)*wtk(ik)
             end do
             !
          end do
       end do
       !
    end do
    ! return variation of local density matrix with respect to the target values
    delta_local_density_matrix = local_density_matrix
    do istate=1,Ns
       delta_local_density_matrix(1,istate,istate) = delta_local_density_matrix(1,istate,istate) - n0_target(istate)      
    end do
    write(*,*) '--> inside_fix_density <--'
    do istate=1,Ns
       write(*,'(12F8.4)') local_density_matrix(1,istate,:),n0_target(istate)
    end do
    write(*,*)
    write(*,*) 'lm_',lm_
    write(*,*)
    do istate=1,Ns
       write(*,'(12F8.4)') lm(1,istate,1:Ns)
    end do
    ! write(*,*)
    ! do istate=1,Ns
    !    write(*,'(12F8.4)') delta_local_density_matrix(1,istate,:)
    ! end do
    ! write(*,*)
    ! do istate=1,Ns
    !    write(*,'(12F8.4)') dimag(delta_local_density_matrix(1,istate,:))
    !    !write(*,'(12F8.4)') dimag(lm(istate,:))
    ! end do
    ! write(*,*)
    ! do istate=1,Ns
    !    write(*,'(12F8.4)') dimag(delta_local_density_matrix(2,istate,:))
    !    !write(*,'(12F8.4)') dreal(lm(istate,:))
    ! end do
    ! write(*,*)
    ! write(*,*) 'DEBUG TEST LOCAL'
    ! do istate=1,2*Ns
    !    write(*,'(12F8.4)') dreal(local_density_matrix(istate,:))
    !    !write(*,'(12F8.4)') dreal(lm(istate,:))
    ! end do
    ! write(*,*)
    ! do istate=1,2*Ns
    !    write(*,'(12F8.4)') dimag(local_density_matrix(istate,:))
    !    !write(*,'(12F8.4)') dimag(lm(istate,:))
    ! end do
    ! write(*,*)

    delta=0.d0
    allocate(delta_cmplx(Nvdm_NC_opt))
    call vdm_NC_stride_m2v(delta_local_density_matrix(1,:,:),delta_cmplx)
    do i=1,Nvdm_NC_opt
       delta(i) = dreal(delta_cmplx(i))
       delta(i+Nvdm_NC_opt) = dimag(delta_cmplx(i))
    end do
    deallocate(delta_cmplx)
    i0 = 2*Nvdm_NC_opt
    allocate(delta_cmplx(Nvdm_AC_opt))
    call vdm_AC_stride_m2v(delta_local_density_matrix(2,:,:),delta_cmplx)
    do i=1,Nvdm_AC_opt
       delta(i0+i) = dreal(delta_cmplx(i))
       delta(i0+i+Nvdm_AC_opt) = dimag(delta_cmplx(i))       
    end do
    deallocate(delta_cmplx)
    write(*,*) delta
    
    !
    ! delta = 0.d0
    ! do istate=1,Ns
    !    do jstate=1,Ns
    !       !
    !       imap = opt_map(istate,jstate)
    !       jmap = opt_map_anomalous(istate,jstate)          
    !       !

    !       if(imap.gt.0) then
    !          delta(imap) = delta(imap) + dreal(delta_local_density_matrix(1,istate,jstate))**2.d0          
    !          delta(imap+Nopt_lgr) = delta(imap+Nopt_lgr) + dimag(delta_local_density_matrix(1,istate,jstate))**2.d0          
    !       end if

    !       if(jmap.gt.0) then 
    !          delta(jmap+Nopt_normal) = delta(jmap+Nopt_normal) +  dreal(delta_local_density_matrix(2,istate,jstate))**2.d0
    !          delta(jmap+Nopt_normal+Nopt_lgr) = delta(jmap+Nopt_normal+Nopt_lgr) +  dimag(delta_local_density_matrix(2,istate,jstate))**2.d0
    !       end if
    !       !
    !    end do
    ! end do
    !
  end function fix_density


end subroutine slater_minimization_lgr_superc




subroutine slater_minimization_fixed_lgr_superc(Rhop,Qhop,lm,Estar,n0,slater_derivatives,slater_matrix_el)     
  complex(8),dimension(Ns,Ns),intent(in)      :: Rhop
  complex(8),dimension(Ns,Ns),intent(in)      :: Qhop  

  complex(8),dimension(2,Ns,Ns),intent(in)    :: lm  
  real(8)                                     :: Estar
  ! optional inputs
  complex(8),dimension(2,Ns,Ns),optional      :: n0
  complex(8),dimension(2,Ns,Ns),optional      :: slater_derivatives
  complex(8),dimension(2*Ns,2*Ns,Lk),optional :: slater_matrix_el
  complex(8),dimension(2,Ns,Ns)               :: n0_
  complex(8),dimension(2*Ns,2*Ns)             :: n0_tmp
  complex(8),dimension(2,Ns,Ns)               :: slater_derivatives_
  complex(8),dimension(2*Ns,2*Ns,Lk)          :: slater_matrix_el_

  complex(8),dimension(Ns,Ns)                 :: Rhop_dag
  complex(8),dimension(Ns,Ns)                 :: Qhop_dag  
  complex(8),dimension(2*Ns,2*Ns)             :: Hk,Hstar,tmp
  complex(8),dimension(Ns,Ns)                 :: Hk_tmp
  real(8),dimension(2*Ns)                     :: ek

  complex(8),dimension(Ns,Ns)                 :: Hk_bare
  integer                                     :: iorb,jorb,ispin,jspin,istate,jstate,kstate,kkstate,ik
  integer                                     ::is,js,ks,kks
  real(8)                                     :: nqp
  !
  Estar=0.d0
  slater_derivatives_=zero
  n0_ = 0.d0
  n0_tmp=zero
  slater_matrix_el_ = zero
  !
  do istate=1,Ns
     do jstate=1,Ns
        Rhop_dag(istate,jstate) = conjg(Rhop(jstate,istate))
        Qhop_dag(istate,jstate) = conjg(Qhop(jstate,istate))
     end do
  end do
  !
  do ik=1,Lk     
     Hk_bare=Hk_tb(:,:,ik)
     Hk=0.d0
     ek=0.d0
     !
     ! hopping renormalization !
     !
     Hk_tmp=matmul(Hk_tb(:,:,ik),Rhop)
     Hk_tmp=matmul(Rhop_dag,Hk_tmp)
     Hk(1:Ns,1:Ns) = Hk_tmp 
     !
     Hk_tmp=matmul(Hk_tb(:,:,ik),Qhop)
     Hk_tmp=matmul(Rhop_dag,Hk_tmp)
     Hk(1:Ns,Ns+1:2*Ns) = Hk_tmp 
     !
     Hk_tmp=matmul(Hk_tb(:,:,ik),Rhop)
     Hk_tmp=matmul(Qhop_dag,Hk_tmp)
     Hk(Ns+1:2*Ns,1:Ns) = Hk_tmp
     !
     Hk_tmp=matmul(Hk_tb(:,:,ik),Qhop)
     Hk_tmp=matmul(Qhop_dag,Hk_tmp)
     Hk(Ns+1:2*Ns,Ns+1:2*Ns) = Hk_tmp
     !
     ! add Lagrange multipliers !
     Hstar = Hk
     !Hk = Hk+lm
     Hk(1:Ns,1:Ns)=Hk(1:Ns,1:Ns)+lm(1,:,:)
     Hk(1:Ns,Ns+1:2*Ns)=Hk(1:Ns,Ns+1:2*Ns)+lm(2,:,:)  
     do is=1,Ns
        do js=1,Ns
           Hk(is+Ns,js) = Hk(is+Ns,js) + conjg(lm(2,js,is))
        end do
     end do
     ! diagonalize hamiltonian !
     call  matrix_diagonalize(Hk,ek) 
     !compute local density matrix (now 2Ns x 2Ns)
     do is = 1,2*Ns
        do js = 1,2*Ns           
           do ks = 1,Ns
              nqp = fermi(ek(ks)-ek(ks+Ns),beta)
              !
              Estar = Estar + &
                   conjg(Hk(is,ks))*Hk(js,ks)*Hstar(is,js)*nqp*wtk(ik) + &
                   conjg(Hk(is,ks+Ns))*Hk(js,ks+Ns)*Hstar(is,js)*(1.d0-nqp)*wtk(ik)
              !
              n0_tmp(is,js) = n0_tmp(is,js) + & 
                   conjg(Hk(is,ks))*Hk(js,ks)*nqp*wtk(ik) + &
                   conjg(Hk(is,ks+Ns))*Hk(js,ks+Ns)*(1.d0-nqp)*wtk(ik)
              !
              slater_matrix_el_(is,js,ik) = slater_matrix_el_(is,js,ik) + &
                   conjg(Hk(is,ks))*Hk(js,ks)*nqp + &
                   conjg(Hk(is,ks+Ns))*Hk(js,ks+Ns)*(1.d0-nqp)
           end do
        end do
     end do
     !
     tmp = slater_matrix_el_(:,:,ik)     
     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns
                 !
                 slater_derivatives_(1,is,js) = slater_derivatives_(1,is,js) + &
                      (conjg(Rhop(kks,ks))*Hk_bare(kks,is)*tmp(ks,js) + conjg(Qhop(kks,ks))*Hk_bare(kks,is)*tmp(ks+Ns,js))*wtk(ik)
                 !
                 slater_derivatives_(2,is,js) = slater_derivatives_(2,is,js) + &
                      (conjg(Rhop(kks,ks))*Hk_bare(kks,is)*tmp(ks,js+Ns) + conjg(Qhop(kks,ks))*Hk_bare(kks,is)*tmp(ks+Ns,js+Ns))*wtk(ik)
                 !
              end do
           end do
        end do
     end do
     !
  end do
  n0_(1,:,:)=n0_tmp(1:Ns,1:Ns)
  n0_(2,:,:)=n0_tmp(1:Ns,1+Ns:2*Ns)
  if(present(n0)) then
     n0=n0_
  end if
  if(present(slater_derivatives)) then
     slater_derivatives=slater_derivatives_
  end if
  if(present(slater_matrix_el)) then
     slater_matrix_el=slater_matrix_el_
  end if
  !
  !<DEBUG
  ! write(*,*) 'DEBUG slater_derivatives'
  ! do istate=1,Ns
  !    write(*,*) dreal(slater_derivatives_(1,istate,:)),dimag(slater_derivatives_(1,istate,:))
  ! end do
  ! write(*,*)
  ! do istate=1,Ns
  !    write(*,*) dreal(slater_derivatives_(2,istate,:)),dimag(slater_derivatives_(2,istate,:))
  ! end do
  !DEBUG>


  !<DEBUG
  write(*,*) 'N0_SLATER'
  do istate=1,2*Ns
     write(*,'(20F8.4)') dreal(n0_tmp(istate,:)),dimag(n0_tmp(istate,:))
  end do
  write(*,*)
  do istate=1,Ns
     write(*,'(20F8.4)') dreal(lm(1,istate,:)),dimag(lm(1,istate,:))
  end do
  write(*,*)
  do istate=1,Ns
     write(*,'(20F8.4)') dreal(lm(2,istate,:)),dimag(lm(2,istate,:))
  end do
  !stop
  !DEBUG>

  


  ! do istate=1,2*Ns
  !    write(*,'(10F8.4)') dreal(n0(istate,:))
  ! end do
  ! write(*,*)
  ! do istate=1,2*Ns
  !    write(*,'(10F8.4)') dimag(n0(istate,:))
  ! end do
  ! write(*,*)
  ! do istate=1,2*Ns
  !    write(*,'(10F8.4)') dreal(lm(istate,:))
  ! end do
  ! write(*,*)
  ! do istate=1,2*Ns
  !    write(*,'(10F8.4)') dimag(lm(istate,:))
  ! end do
  !
  ! write(*,*) size(slater_matrix_el,1),size(slater_matrix_el,2),size(slater_matrix_el,3)
  ! write(*,*) size(slater_matrix_el_,1),size(slater_matrix_el_,2),size(slater_matrix_el_,3)
  !
  !<DEBUG
  !write(*,*) 'Ma QUA TI PIANTI?!?!?!'
  !DEBUG>

end subroutine slater_minimization_fixed_lgr_superc
