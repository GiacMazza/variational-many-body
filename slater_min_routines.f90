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
  real(8),dimension(:),allocatable                :: lgr     !+- real indeendent lgr_vector -+!


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
  lgr=-0.5
  do is=1,Ns
     do js=1,Ns
        imap = opt_map(is,js)
        if(imap.gt.0) then
           if(is.eq.js) then
              if(n0_target(is).le.1.d-4)  lgr(imap) = Wband
              if(n0_target(is).ge.1.d0-1.d-4)  lgr(imap) =  -Wband
           end if
        end if
     end do
  end do
  !<DEBUG
  !DEBUG>
  lgr = lgr_init_slater
  call fmin_cgminimize(lgr,get_delta_local_density_matrix,iter,delta,itmax=20)
  lgr_init_slater=lgr
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
     write(*,*) "Slater Determinant: Lagrange Multipliers - OK -"
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
          !imap = vdm_c_map(istate,jstate)
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
  end function get_delta_local_density_matrix
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
  integer                                           :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  !
  Estar=0.d0
  slater_derivatives_=zero
  n0_ = 0.d0
  slater_matrix_el_ = zero
  !

 
  do ik=1,Lk
     !
     Hk=0.d0
     ek=0.d0
     Hk_bare=Hk_tb(:,:,ik)
     ! hopping renormalization !
     Hk=matmul(Hk_tb(:,:,ik),Rhop)
     Hk=matmul(Rhop,Hk)
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
                    tmp(istate,jstate) = tmp(istate,jstate) + Hk(istate,kstate)*fermi(ek(kstate),beta)*conjg(Hk(jstate,kstate))
                    Estar = Estar + Hk(istate,kstate)*Hk(jstate,kstate)*Hstar(istate,jstate)*fermi(ek(kstate),beta)*wtk(ik)
                    n0_(istate,jstate) = n0_(istate,jstate) + Hk(istate,kstate)*fermi(ek(kstate),beta)*conjg(Hk(jstate,kstate))*wtk(ik)                    
                    slater_matrix_el_(istate,jstate,ik) = slater_matrix_el_(istate,jstate,ik) + &
                         Hk(istate,kstate)*Hk(jstate,kstate)*fermi(ek(kstate),beta)
                 end do
              end do
           end do
        end do
     end do
     ! store slater ground state derivatives
     tmp=matmul(Rhop,tmp)
     tmp=matmul(Hk_bare,tmp)             
     do istate=1,Ns
        do jstate=1,Ns             
           slater_derivatives_(istate,jstate) = &
                slater_derivatives_(istate,jstate) + 2.d0*tmp(istate,jstate)*wtk(ik)
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
end subroutine slater_minimization_fixed_lgr



