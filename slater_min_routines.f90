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
  call slater_minimization_fixed_lgr(Rhop,lgr_multip,Estar,n0_out_,slater_derivatives_,slater_matrix_el_)

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
    !
    ! delta_local_density_matrix = local_density_matrix
    ! do istate=1,Ns
    !    delta_local_density_matrix(istate,istate) = delta_local_density_matrix(istate,istate) - n0_target(istate)      
    ! end do
    ! delta=0.d0
    ! do istate=1,Ns
    !    do jstate=1,Ns
    !       delta = delta + abs(delta_local_density_matrix(istate,jstate))**2.d0
    !    end do
    ! end do
    ! !
  end function fix_density


end subroutine slater_minimization_lgr
!
subroutine slater_minimization_fixed_lgr(Rhop,lm,Estar,n0,slater_derivatives,slater_matrix_el)     
  complex(8),dimension(Ns,Ns),intent(in) :: Rhop
  real(8),dimension(Ns,Ns),intent(in)           :: lm
  real(8)                                           :: Estar
  ! optional inputs
  real(8),dimension(Ns,Ns),optional            :: n0
  complex(8),dimension(Ns,Ns),optional            :: slater_derivatives
  complex(8),dimension(Ns,Ns,Lk),optional            :: slater_matrix_el
  real(8),dimension(Ns,Ns) :: n0_
  complex(8),dimension(Ns,Ns)  :: slater_derivatives_
  complex(8),dimension(Ns,Ns,Lk)            :: slater_matrix_el_
  complex(8),dimension(Ns,Ns)            :: Hk,tmp,Hk_bare,Hstar
  real(8),dimension(Ns)                      :: ek
  integer                                           :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,is,js,ks

  complex(8),dimension(Ns,Ns) :: Rhop_dag
  
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

     do istate=1,Ns
        do jstate=1,Ns
           !+-  store slater ground state derivatives -+!
           do is=1,Ns
              do js=1,Ns
                 do ks=1,Ns
                    slater_derivatives_(istate,jstate) = slater_derivatives_(istate,jstate) + &
                         Hk(jstate,is)*fermi(ek(is),beta)*conjg(Hk(js,is))*Rhop_dag(js,ks)*Hk_bare(ks,istate)*wtk(ik)
                 end do
              end do
           end do           
        end do
     end do
     !
     ! tmp=matmul(Rhop_dag,tmp)
     ! tmp=matmul(Hk_bare,tmp)             
     ! do istate=1,Ns
     !    do jstate=1,Ns             
     !       slater_derivatives_(istate,jstate) = &
     !            slater_derivatives_(istate,jstate) + 1.d0*tmp(istate,jstate)*wtk(ik)
     !    end do
     ! end do
     !
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
     write(*,*) n0_(istate,:)
  end do
  !DEBUG>

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
  complex(8),dimension(Ns,Ns),intent(in)    :: Rhop         !input:  normal    hopping renormalization matrix
  complex(8),dimension(Ns,Ns),intent(in)    :: Qhop         !input:  anomalous hopping renormalization matrix 
  !
  real(8),dimension(Ns),intent(in)          :: n0_target    !input:  variational density matrix in the natural basis
  real(8),intent(out)                       :: Estar        !output: Slater Deter GS energy
  complex(8),dimension(2*Ns,2*Ns),intent(out) :: lgr_multip   !output: Slater Deter lagrange multipliers
                                                            !+- this rather should be (2,Ns,Ns)
  complex(8),dimension(2*Ns,2*Ns),optional         :: n0_out             !output: Slater Deter ground state density matrix
  complex(8),dimension(2,Ns,Ns),optional    :: slater_derivatives !output: Slater Deter ground state renormalization derivatives
  complex(8),dimension(2*Ns,2*Ns,Lk),optional :: slater_matrix_el   !output: Slater Deter ground state matrix elements
  !
  logical,optional                          :: iverbose     !input:  Verbosity level
  !
  complex(8),dimension(2,Ns,Ns)             :: slater_derivatives_
  complex(8),dimension(2*Ns,2*Ns,Lk)          :: slater_matrix_el_
  complex(8),dimension(2*Ns,2*Ns)                  :: n0_out_
  !
  real(8),dimension(:),allocatable                :: lgr
  real(8),dimension(:),allocatable ::   delta_out     !+- real indeendent lgr_vector -+!

  complex(8),dimension(Ns,Ns)    :: Rhop_dag         !input:  normal    hopping renormalization matrix
  complex(8),dimension(Ns,Ns)    :: Qhop_dag         !input:  anomalous hopping renormalization matrix 


  complex(8),dimension(Ns,Ns)             :: Hk
  real(8),dimension(Ns)                :: tmp_lgr,err_dens
  integer                              :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb,iter,imap,jmap,is,js
  real(8),dimension(Ns)                :: lgr_multip_vec
  logical                              :: iverbose_
  real(8)  :: delta
  !
  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose
  !    
  !<DEBUG
  ! write(*,*) 'Rhop'
  ! do istate=1,Ns    
  !    write(*,*) dreal(Rhop(istate,:))
  ! end do
  ! write(*,*) 'Qhop'
  ! do istate=1,Ns    
  !    write(*,*) dreal(Qhop(istate,:))
  ! end do
  !DEBUG>

  
  allocate(lgr(2*Nopt_lgr))
  lgr=0.d0
  !
  select case(lgr_method)
  case('CG_min')
     call fmin_cg(lgr,get_delta_local_density_matrix,iter,delta,itmax=20)
  case('f_zero')
     call fsolve(fix_density,lgr,tol=1.d-10,info=iter)
     delta_out=fix_density(lgr)
     write(*,*) delta_out
     write(*,*) n0_target
     delta=0.d0
     do is=1,2*Nopt_lgr
        delta = delta + delta_out(is)**2.d0
     end do
  end select
  
  !
  lgr_multip=0.d0
  do istate=1,Ns
     do jstate=1,Ns
        imap = opt_map(istate,jstate)
        jmap = opt_map_anomalous(istate,jstate)
        if(imap.gt.0) lgr_multip(istate,jstate) = lgr(imap) + xi*lgr(imap+Nopt_lgr)
        if(jmap.gt.0) then
           lgr_multip(istate,jstate+Ns) = lgr(jmap+Nopt_normal) + xi*lgr(jmap+Nopt_normal+Nopt_lgr)
           lgr_multip(jstate+Ns,istate) = lgr(jmap+Nopt_normal) - xi*lgr(jmap+Nopt_normal+Nopt_lgr)
        end if
     end do
  end do
  !
  call slater_minimization_fixed_lgr_superc(Rhop,Qhop,lgr_multip,Estar,n0_out_,slater_derivatives_,slater_matrix_el_)

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
     do is=1,2*Ns
        write(*,'(20F7.3)') dreal(lgr_multip(is,:))
     end do
     write(*,*)
     do is=1,2*Ns
        write(*,'(20F7.3)') dimag(lgr_multip(is,:))
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
    real(8),dimension(:)          :: lm_
    complex(8),dimension(Nopt_lgr)       :: lm_cmplx
    complex(8),dimension(2*Ns,2*Ns) :: lm
    real(8)                :: delta

    complex(8),dimension(2,Ns,Ns)  :: delta_local_density_matrix,local_density_matrix
    complex(8),dimension(2*Ns,2*Ns) :: Hk
    complex(8),dimension(Ns,Ns) :: Hk_tmp
    real(8),dimension(2*Ns)          :: ek
    real(8)                          :: nqp
    integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,imap,jmap,i
    integer :: is,js,ks
    !
    !+- divide lgr_optimization parameters in normal and anomalous (--> complex <--)
    !   Nopt_lgr = Nopt_normal + Nopt_anomalous

    !+- reconstruct complex lgr paramters -+!
    do i=1,Nopt_lgr
       lm_cmplx(i) = lm_(i) + xi*lm_(i+Nopt_lgr)
    end do
    !+- reconstruct full optimization parameter matrices -+!
    lm=0.d0
    do istate=1,Ns
       do jstate=1,Ns
          imap = opt_map(istate,jstate)
          if(imap.gt.0) lm(istate,jstate)=lm_cmplx(imap)
          jmap = opt_map_anomalous(istate,jstate)
          if(jmap.gt.0) then
             lm(istate,jstate+Ns) = lm_cmplx(jmap+Nopt_normal)          
             lm(jstate+Ns,istate) = conjg(lm_cmplx(jmap+Nopt_normal))
          end if
       end do
    end do
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
       Hk=Hk+lm                     
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
    real(8),dimension(:)          :: lm_
    real(8),dimension(size(lm_))      :: delta
    complex(8),dimension(Nopt_lgr)       :: lm_cmplx
    complex(8),dimension(2*Ns,2*Ns) :: lm   !+- this may be also (2,Ns,Ns)

    complex(8),dimension(2,Ns,Ns)  :: delta_local_density_matrix,local_density_matrix
    complex(8),dimension(2*Ns,2*Ns) :: Hk
    complex(8),dimension(Ns,Ns) :: Hk_tmp
    real(8),dimension(2*Ns)          :: ek
    integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,imap,jmap,i
    integer                        :: is,js,ks
    real(8) :: nqp
    !
    !+- divide lgr_optimization parameters in normal and anomalous (--> complex <--)
    !   Nopt_lgr = Nopt_normal + Nopt_anomalous

    !+- reconstruct complex lgr paramters -+!

    !write(*,*) 'QUA DENTRO CI ARRIVO'

    do i=1,Nopt_lgr
       lm_cmplx(i) = lm_(i) + xi*lm_(i+Nopt_lgr)
    end do
    !+- reconstruct full optimization parameter matrices -+!
    lm=0.d0
    do istate=1,Ns
       do jstate=1,Ns
          imap = opt_map(istate,jstate)
          if(imap.gt.0) lm(istate,jstate)=lm_cmplx(imap)
          jmap = opt_map_anomalous(istate,jstate)
          if(jmap.gt.0) then
             lm(istate,jstate+Ns) = lm_cmplx(jmap+Nopt_normal)          
             lm(jstate+Ns,istate) = conjg(lm_cmplx(jmap+Nopt_normal))
          end if
       end do
    end do
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
       Hk=Hk+lm                     
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

    ! write(*,*) 'DEBUG TEST LM'
    ! do istate=1,2*Ns
    !    !write(*,'(12F8.4)') dreal(delta_local_density_matrix(istate,:))
    !    write(*,'(12F8.4)') dreal(lm(istate,:))
    ! end do
    ! write(*,*)
    ! do istate=1,2*Ns
    !    !write(*,'(12F8.4)') dreal(delta_local_density_matrix(istate,:))
    !    write(*,'(12F8.4)') dimag(lm(istate,:))
    ! end do
    ! write(*,*)
    
    delta_local_density_matrix = local_density_matrix
    do istate=1,Ns
       delta_local_density_matrix(1,istate,istate) = delta_local_density_matrix(1,istate,istate) - n0_target(istate)      
    end do
    
    ! write(*,*) 'DEBUG TEST DELTA_LOCAL'
    ! do istate=1,Ns
    !    write(*,'(12F8.4)') dreal(delta_local_density_matrix(1,istate,:))
    !    !write(*,'(12F8.4)') dreal(lm(istate,:))
    ! end do
    ! write(*,*)
    ! do istate=1,Ns
    !    write(*,'(12F8.4)') dreal(delta_local_density_matrix(2,istate,:))
    !    !write(*,'(12F8.4)') dreal(lm(istate,:))
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



    delta = 0.d0
    do istate=1,Ns
       do jstate=1,Ns
          !
          imap = opt_map(istate,jstate)
          jmap = opt_map_anomalous(istate,jstate)          
          !

          if(imap.gt.0) then
             !write(*,*) 'imap',imap,istate,jstate
             delta(imap) = delta(imap) + dreal(delta_local_density_matrix(1,istate,jstate))**2.d0          
             delta(imap+Nopt_lgr) = delta(imap+Nopt_lgr) + dimag(delta_local_density_matrix(1,istate,jstate))**2.d0          
             !             write(*,*) imap,dreal(delta_local_density_matrix(istate,jstate)),delta(imap)
             !             write(*,*) imap+Nopt_lgr
          end if
          
          if(jmap.gt.0) then 
             !write(*,*) 'jmap',jmap,istate,jstate
             delta(jmap+Nopt_normal) = delta(jmap+Nopt_normal) +  dreal(delta_local_density_matrix(2,istate,jstate))**2.d0
             delta(jmap+Nopt_normal+Nopt_lgr) = delta(jmap+Nopt_normal+Nopt_lgr) +  dimag(delta_local_density_matrix(2,istate,jstate))**2.d0
             !write(*,*) jmap+Nopt_normal
             !write(*,*) jmap+Nopt_lgr+Nopt_normal
          end if
          !
       end do
    end do
    !
    ! write(*,*) 'lm_',lm_
    ! write(*,'(12F8.4)') delta(:)
    ! write(*,'(20F5.2)') delta_local_density_matrix(1,:,:)
    ! write(*,'(20F5.2)') delta_local_density_matrix(2,:,:)
    ! write(*,*)
    ! stop
  end function fix_density
  

end subroutine slater_minimization_lgr_superc




subroutine slater_minimization_fixed_lgr_superc(Rhop,Qhop,lm,Estar,n0,slater_derivatives,slater_matrix_el)     
  complex(8),dimension(Ns,Ns),intent(in)     :: Rhop
  complex(8),dimension(Ns,Ns),intent(in)     :: Qhop  

  complex(8),dimension(2*Ns,2*Ns),intent(in) :: lm  
  real(8)                                    :: Estar
  ! optional inputs
  complex(8),dimension(2*Ns,2*Ns),optional          :: n0
  complex(8),dimension(2,Ns,Ns),optional            :: slater_derivatives
  complex(8),dimension(2*Ns,2*Ns,Lk),optional            :: slater_matrix_el
  complex(8),dimension(2*Ns,2*Ns) :: n0_
  complex(8),dimension(2,Ns,Ns)  :: slater_derivatives_
  complex(8),dimension(2*Ns,2*Ns,Lk)            :: slater_matrix_el_

  complex(8),dimension(Ns,Ns)     :: Rhop_dag
  complex(8),dimension(Ns,Ns)     :: Qhop_dag  
  complex(8),dimension(2*Ns,2*Ns) :: Hk,Hstar,tmp
  complex(8),dimension(Ns,Ns) :: Hk_tmp
  real(8),dimension(2*Ns)          :: ek

  complex(8),dimension(Ns,Ns)            :: Hk_bare
  integer                                           :: iorb,jorb,ispin,jspin,istate,jstate,kstate,kkstate,ik
  integer ::is,js,ks,kks
  real(8) :: nqp
  !
  Estar=0.d0
  slater_derivatives_=zero
  n0_ = 0.d0
  slater_matrix_el_ = zero
  !
  do istate=1,Ns
     do jstate=1,Ns
        Rhop_dag(istate,jstate) = conjg(Rhop(jstate,istate))
        Qhop_dag(istate,jstate) = conjg(Qhop(jstate,istate))
     end do
  end do
  

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
     Hk = Hk+lm                     
     ! diagonalize hamiltonian !
     call  matrix_diagonalize(Hk,ek,'V','L')
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
              n0_(is,js) = n0_(is,js) + & 
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
