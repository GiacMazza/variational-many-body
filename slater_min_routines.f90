subroutine slater_minimization_lgr(Rhop,n0_target,Estar,lgr_multip,n0_out,slater_derivatives,slater_matrix_el,iverbose) 
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
    
    ! write(*,*) 'delta_CMPLX'
    ! write(*,*) delta_cmplx
    
    ! check = 0.d0
    ! do is=1,Ns
    !    do js=1,Ns
    !       check = check + & 
    !            (delta_local_density_check(is,js)-delta_local_density_matrix(is,js))*conjg(delta_local_density_check(is,js)-delta_local_density_matrix(is,js))
    !    end do
    !    !write(*,*) dreal(delta_local_density_check(is,:))
    !    write(*,*) dreal(delta_local_density_matrix(is,:))
    !    write(*,*)
    ! end do
    ! write(*,*)
    ! write(*,*)
    ! if(check.gt.1.d-10) stop "CHECK STRIDES @ fix_density_normal"
    write(*,*) delta
  end function fix_density
end subroutine slater_minimization_lgr
!
subroutine slater_minimization_fixed_lgr(Rhop,lm,Estar,n0,slater_derivatives,slater_matrix_el,store)
  complex(8),dimension(Ns,Ns),intent(in)  :: Rhop
  complex(8),dimension(Ns,Ns),intent(in)  :: lm
  real(8)                                 :: Estar
  ! optional inputs
  complex(8),dimension(Ns,Ns),optional    :: n0 
  complex(8),dimension(Ns,Ns),optional    :: slater_derivatives
  complex(8),dimension(Ns,Ns,Lk),optional :: slater_matrix_el
  logical,optional                        :: store
  logical                                 :: store_
  complex(8),dimension(Ns,Ns)             :: n0_
  complex(8),dimension(Ns,Ns)             :: slater_derivatives_
  complex(8),dimension(Ns,Ns,Lk)          :: slater_matrix_el_
  complex(8),dimension(Ns,Ns)             :: Hk,tmp,Hk_bare,Hstar
  real(8),dimension(Ns)                   :: ek
  real(8),dimension(Ns,Lk)                   :: ek_store
  integer                                 :: i,iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,is,js,ks,kks
  integer                                 :: unit_store_slater_el,unit_store_slater_ek,istore,unit_store_qp
  complex(8),dimension(Ns*Ns)             :: tmp_matrix_el
  !
  complex(8),dimension(Ns,Ns)             :: Rhop_dag  
  character(len=20)                        :: store_file_suffix
  character(len=17)                        :: store_file

  complex(8),dimension(Ns,Ns,lw)          :: qp_gloc
  complex(8)                :: tmp_gk,iw
  real(8) :: w
  
  !
  Estar=0.d0
  slater_derivatives_=zero
  n0_ = 0.d0
  slater_matrix_el_ = zero
  !
  store_=.false.; if(present(store)) store_=store
  if(store_) then
     unit_store_slater_el = free_unit()
     open(unit_store_slater_el,file='store_slater_determinant_ground_state_el.out')
  end if
  do is=1,Ns
     do js=1,Ns
        Rhop_dag(is,js) = conjg(Rhop(js,is))
     end do
  end do
  !
  if(store_) qp_gloc = zero
  !
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
     ek_store(:,ik) = ek
     !
     ! store slater determinant matrix elements
     istore=0; tmp_matrix_el=zero
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
           istore=istore+1
           tmp_matrix_el(istore) = slater_matrix_el_(istate,jstate,ik)

           if(store_) then !+- compute the single-particle greens function for the QP spectrum      
              do i=1,lw
                 do kstate=1,Ns
                    w=wr(i)
                    iw=cmplx(w,0.05d0)
                    tmp_gk = 1.d0/(iw-ek(kstate))
                    qp_gloc(istate,jstate,i) = qp_gloc(istate,jstate,i) + Hk(istate,kstate)*conjg(Hk(jstate,kstate))*tmp_gk*wtk(ik)
                 enddo
              end do
           end if

        end do
     end do
     !
     if(store_) then
        do is=1,Ns*Ns
           write(unit_store_slater_el,'(10(F18.10))') tmp_matrix_el(is)
        end do
        write(unit_store_slater_el,'(10(F18.10))') 
     end if
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
  if(store_) then
     do is=1,Ns
        do js=is,Ns
           unit_store_qp = free_unit()
           store_file_suffix="_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))
           open(unit_store_qp,file="QP_GLOC_realw"//reg(store_file_suffix)//".data")
           !
           do i=1,lw
              write(unit_store_qp,'(5(F18.10))') wr(i),-dimag(qp_gloc(is,js,i))/pi,dreal(qp_gloc(is,js,i)),dimag(qp_gloc(is,js,i))
           end do
           close(unit_store_qp)
        end do
     end do
  end if
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
end subroutine slater_minimization_fixed_lgr





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
  complex(8),dimension(2,Ns,Ns,Lk),optional :: slater_matrix_el   !output: Slater Deter ground state matrix elements
  !
  logical,optional                            :: iverbose     !input:  Verbosity level
  !
  complex(8),dimension(2,Ns,Ns)               :: slater_derivatives_
  complex(8),dimension(2,Ns,Ns,Lk)          :: slater_matrix_el_
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
     do is=1,2*Nopt
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
    !
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
    !< TMP
    ! write(*,*) '--> inside_fix_density <--'
    ! do istate=1,Ns
    !    write(*,'(12F8.4)') local_density_matrix(1,istate,:),n0_target(istate)
    ! end do
    ! write(*,*)
    ! write(*,*) 'lm_',lm_
    ! write(*,*)
    ! do istate=1,Ns
    !    write(*,'(12F8.4)') lm(1,istate,1:Ns)
    ! end do
    !TMP>
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
    !write(*,*) delta
  end function fix_density


end subroutine slater_minimization_lgr_superc
!
subroutine slater_minimization_fixed_lgr_superc(Rhop,Qhop,lm,Estar,n0,slater_derivatives,slater_matrix_el,store)     
  complex(8),dimension(Ns,Ns),intent(in)      :: Rhop
  complex(8),dimension(Ns,Ns),intent(in)      :: Qhop  
  complex(8),dimension(2,Ns,Ns),intent(in)    :: lm  
  real(8)                                     :: Estar
  ! optional inputs
  complex(8),dimension(2,Ns,Ns),optional      :: n0
  complex(8),dimension(2,Ns,Ns),optional      :: slater_derivatives
  complex(8),dimension(2,Ns,Ns,Lk),optional :: slater_matrix_el
  logical,optional                            :: store
  complex(8),dimension(2,Ns,Ns)               :: n0_
  complex(8),dimension(2*Ns,2*Ns)             :: n0_tmp
  complex(8),dimension(2,Ns,Ns)               :: slater_derivatives_
  complex(8),dimension(2*Ns,2*Ns,Lk)          :: slater_superc_matrix_el_
  complex(8),dimension(2,Ns,Ns,Lk)            :: slater_matrix_el_
  logical                                     :: store_
  !
  complex(8),dimension(Ns,Ns)                 :: Rhop_dag
  complex(8),dimension(Ns,Ns)                 :: Qhop_dag  
  complex(8),dimension(2*Ns,2*Ns)             :: Hk,Hstar,tmp
  complex(8),dimension(Ns,Ns)                 :: Hk_tmp
  real(8),dimension(2*Ns)                     :: ek
  real(8),dimension(Ns)                       :: eps_ik
  !
  complex(8),dimension(Ns,Ns)                 :: Hk_bare
  integer                                     :: iorb,jorb,ispin,jspin,istate,jstate,kstate,kkstate,ik
  integer                                     ::is,js,ks,kks
  real(8)                                     :: nqp
  !
  integer                                 :: unit_store_slater_el,unit_store_slater_ek,istore,i,unit_store_qp
  complex(8),dimension(2*Ns*2*Ns)             :: tmp_matrix_el
  !
  complex(8),dimension(Ns,Ns,lw)          :: qp_gloc
  complex(8)                :: tmp_gk,iw
  real(8) :: w
  character(len=20)                        :: store_file_suffix
  !

  Estar=0.d0
  slater_derivatives_=zero
  n0_ = 0.d0
  n0_tmp=zero
  slater_superc_matrix_el_ = zero
  !
  do istate=1,Ns
     do jstate=1,Ns
        Rhop_dag(istate,jstate) = conjg(Rhop(jstate,istate))
        Qhop_dag(istate,jstate) = conjg(Qhop(jstate,istate))
     end do
  end do
  !
  store_=.false.; if(present(store)) store_=store
  if(store_) then
     unit_store_slater_el = free_unit()
     open(unit_store_slater_el,file='store_slater_determinant_ground_state_el.out')
  end if
  !
  !
  if(store_) qp_gloc=zero
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
     do is=1,Ns
        eps_ik(is) = ek(is)-ek(is+Ns)
     end do
     !compute local density matrix (now 2Ns x 2Ns)
     istore=0; tmp_matrix_el=zero
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
              slater_superc_matrix_el_(is,js,ik) = slater_superc_matrix_el_(is,js,ik) + &
                   conjg(Hk(is,ks))*Hk(js,ks)*nqp + &
                   conjg(Hk(is,ks+Ns))*Hk(js,ks+Ns)*(1.d0-nqp)
           end do
           istore=istore+1
           tmp_matrix_el(istore) = slater_superc_matrix_el_(istate,jstate,ik)
        end do
     end do
     !
     if(store_) then
        do is=1,Ns
           do js=1,Ns
              if(store_) then !+- compute the single-particle greens function for the QP spectrum      
                 do i=1,lw
                    do kstate=1,Ns
                       w=wr(i)
                       iw=cmplx(w,0.01d0)
                       tmp_gk = 1.d0/(iw-(ek(kstate)-ek(kstate+Ns)))
                       qp_gloc(is,js,i) = qp_gloc(is,js,i) + Hk(is,kstate)*conjg(Hk(js,kstate))*tmp_gk*wtk(ik)
                       tmp_gk = 1.d0/(iw+(ek(kstate)-ek(kstate+Ns)))
                       qp_gloc(is,js,i) = qp_gloc(is,js,i) + Hk(is,kstate+Ns)*conjg(Hk(js,kstate+Ns))*tmp_gk*wtk(ik)
                    enddo
                 end do
              end if
           end do
        end do
     end if
     !
     if(store_) then
        do is=1,4*Ns*Ns
           write(unit_store_slater_el,'(10(F18.10))') tmp_matrix_el(is)
        end do
        write(unit_store_slater_el,'(10(F18.10))') 
     end if
     !
     tmp = slater_superc_matrix_el_(:,:,ik)     
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
  !
  slater_matrix_el_(1,1:Ns,1:Ns,:) = slater_superc_matrix_el_(1:Ns,1:Ns,:)
  slater_matrix_el_(2,1:Ns,1:Ns,:) = slater_superc_matrix_el_(Ns+1:2*Ns,1:Ns,:)
  !
  n0_(1,:,:)=n0_tmp(1:Ns,1:Ns)
  n0_(2,:,:)=n0_tmp(1:Ns,1+Ns:2*Ns)
  !
  if(store_) then
     do is=1,Ns
        do js=is,Ns
           unit_store_qp = free_unit()
           store_file_suffix="_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))
           open(unit_store_qp,file="QP_GLOC_realw"//reg(store_file_suffix)//".data")
           !
           do i=1,lw
              write(unit_store_qp,'(5(F18.10))') wr(i),-dimag(qp_gloc(is,js,i))/pi,dreal(qp_gloc(is,js,i)),dimag(qp_gloc(is,js,i))
           end do
           close(unit_store_qp)
        end do
     end do
  end if
  !
  !
  !
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
end subroutine slater_minimization_fixed_lgr_superc
