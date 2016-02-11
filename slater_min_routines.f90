subroutine slater_determinant_minimization(Rhop,n0_target,Estar,lgr_multip,iverbose) 
  complex(8),dimension(Ns,Ns),intent(in)  :: Rhop         !input:  renrmalization matrix
  real(8),dimension(Ns),intent(in)        :: n0_target    !input:  variational density matrix
  real(8),intent(out)                  :: Estar        !output: Slater Deter GS energy
  real(8),dimension(Ns,Ns),intent(out)    :: lgr_multip   !output: Slater Deter lagrange multipliers
  real(8),dimension(Nvdm_c)                :: lgr     !+- real indeendent lgr_vector -+!


  logical,optional                     :: iverbose     !input:  Verbosity level
  !

  complex(8),dimension(Ns,Ns)             :: Hk
  real(8),dimension(Ns)                :: tmp_lgr,err_dens
  integer                              :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb,iter,imap,is
  real(8),dimension(Ns)                :: lgr_multip_vec
  logical                              :: iverbose_
  real(8)  :: delta

  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose
  !    
  lgr=0.d0
  call fmin_cg(lgr,get_delta_local_density_matrix,iter,delta)
  lgr_multip=0.d0
  do istate=1,Ns
     do jstate=1,Ns
        imap = vdm_c_map(istate,jstate)
        if(imap.gt.0) lgr_multip(istate,jstate)=lgr(imap)
     end do
  end do  
  call store_slater_ground_state(Rhop,lgr_multip,Estar)
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
          imap = vdm_c_map(istate,jstate)
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
  end function get_delta_local_density_matrix  
  !include 'self_minimization_slater_routines.f90'
  !
end subroutine slater_determinant_minimization







subroutine slater_determinant_minimization_nlep(Rhop,n0_target,Estar,lgr_multip,slater_derivatives,iverbose) 
  complex(8),dimension(Ns,Ns),intent(in)  :: Rhop         !input:  renrmalization matrix
  real(8),dimension(Ns),intent(in)     :: n0_target    !input:  variational density matrix
  real(8),intent(out)                  :: Estar        !output: Slater Deter GS energy
  real(8),dimension(Ns,Ns),intent(out)    :: lgr_multip   !output: Slater Deter lagrange multipliers
  real(8),dimension(Nvdm_c)                :: lgr     !+- real indeendent lgr_vector -+!

  complex(8),dimension(Ns,Ns),intent(out) :: slater_derivatives !output: Slater Deter GS energy derivatives
  logical,optional                     :: iverbose     !input:  Verbosity level
  !

  complex(8),dimension(Ns,Ns)             :: Hk
  real(8),dimension(Ns)                :: tmp_lgr,err_dens
  integer                              :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb,iter,imap,is
  real(8),dimension(Ns)                :: lgr_multip_vec
  logical                              :: iverbose_
  real(8)  :: delta

  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose
  !    
  lgr=0.d0
  call fmin_cg(lgr,get_delta_local_density_matrix,iter,delta)
  !call fsolve(get_delta_local_density_matrix_diag,lgr,tol=1.d-12,info=info)    
  lgr_multip=0.d0
  do istate=1,Ns
     do jstate=1,Ns
        imap = vdm_c_map(istate,jstate)
        if(imap.gt.0) lgr_multip(istate,jstate)=lgr(imap)
     end do
  end do
  !
  call store_slater_ground_state(Rhop,lgr_multip,Estar,slater_derivatives)
  !
  if(iverbose_) then
     write(*,*)
     write(*,*) "Slater Determinant: Lagrange Multipliers - diagonal form -"
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
          imap = vdm_c_map(istate,jstate)
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
  end function get_delta_local_density_matrix
  !include 'self_minimization_slater_routines.f90'
  !
end subroutine slater_determinant_minimization_nlep




!+--------------------------------------------+!
!+- CONSTRAINED MINIMIZATION (cmin) ROUTINES -+!
!+--------------------------------------------+!
! slater
subroutine slater_determinant_minimization_cmin(Rhop,n0_target,Estar,lgr_multip,slater_matrix_el,iverbose) 
  complex(8),dimension(Ns,Ns),intent(in)     :: Rhop         !input:  renrmalization matrix
  real(8),dimension(Ns),intent(in)        :: n0_target    !input:  variational density matrix
  real(8),intent(out)                     :: Estar        !output: Slater Deter GS energy
  real(8),dimension(Ns,Ns),intent(out)       :: lgr_multip   !output: Slater Deter lagrange multipliers
  complex(8),dimension(Ns,Ns,Lk),intent(out) :: slater_matrix_el !output: Slater Deter GS energy derivatives
  logical,optional                        :: iverbose     !input:  Verbosity level
  !
  real(8),dimension(Nvdm)                   :: lgr
  complex(8),dimension(Ns,Ns)                :: Hk,test_dens
  real(8),dimension(Ns)                   :: ek,tmp_lgr,err_dens
  integer                                 :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb,imap,iter
  real(8),dimension(Ns)                   :: lgr_multip_vec

  real(8),dimension(Norb)                 :: lgr_orb,test_dens_orb
  real(8) :: delta
  logical                                 :: iverbose_
  !
  iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose
  !
  ! lgr=0.d0    
  ! call fsolve(get_delta_local_density_matrix_diag,lgr,tol=1.d-12,info=info)    
  ! call store_slater_ground_state_cmin(Rhop,lgr,Estar,slater_matrix_el)


  lgr=0.d0
  call fmin_cg(lgr,get_delta_local_density_matrix,iter,delta)
  !call fsolve(get_delta_local_density_matrix_diag,lgr,tol=1.d-12,info=info)    
  lgr_multip=0.d0
  do istate=1,Ns
     do jstate=1,Ns
        imap = vdm_c_map(istate,jstate)
        if(imap.gt.0) lgr_multip(istate,jstate)=lgr(imap)
     end do
  end do
  call store_slater_ground_state_cmin(Rhop,lgr_multip,Estar,slater_matrix_el)    
  !
  if(iverbose_) then
     write(*,*)
     write(*,*) "Slater Determinant: Lagrange Multipliers - diagonal form -"
     do istate=1,Ns
        write(*,'(20F18.10)') lgr_multip(istate,:)
     end do
     write(*,*) "Slater Determinant: Variational density error"
     !err_dens=get_delta_local_density_matrix_diag(lgr_multip)
     write(*,"(20F18.10)") delta
     write(*,*) "Slater Determinant: Ground State energy"
     write(*,"(20F18.10)") Estar
     write(*,*)
  end if
  !
contains    
  !
  !include 'self_minimization_slater_routines.f90'
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
          imap = vdm_c_map(istate,jstate)
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
  end function get_delta_local_density_matrix

end subroutine slater_determinant_minimization_cmin



subroutine get_slater_ground_state(Rhop,lm,Estar,n0,slater_derivatives)     
  complex(8),dimension(Ns,Ns),intent(in) :: Rhop
  real(8),dimension(Ns,Ns),intent(in)           :: lm
  real(8)                                           :: Estar
  complex(8),dimension(Ns,Ns)            :: slater_derivatives
  real(8),dimension(Ns,Ns) :: n0
  complex(8),dimension(Ns,Ns)  :: slater_derivatives_
  complex(8),dimension(Ns,Ns)            :: Hk,tmp,Hk_bare,Hstar
  real(8),dimension(Ns)                      :: ek
  integer                                           :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  Estar=0.d0
  slater_derivatives_=0.d0
  n0 = 0.d0
  do ik=1,Lk
     !
     Hk=0.d0
     ek=0.d0
     !
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
                    n0(istate,jstate) = n0(istate,jstate) + Hk(istate,kstate)*fermi(ek(kstate),beta)*conjg(Hk(jstate,kstate))*wtk(ik)
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
  slater_derivatives=slater_derivatives_
end subroutine get_slater_ground_state



subroutine store_slater_ground_state(Rhop,lm,Estar,slater_derivatives)     
  complex(8),dimension(Ns,Ns),intent(in) :: Rhop
  real(8),dimension(Ns,Ns),intent(in)           :: lm
  real(8)                                           :: Estar
  complex(8),dimension(Ns,Ns),optional            :: slater_derivatives
  complex(8),dimension(Ns,Ns)  :: slater_derivatives_
  complex(8),dimension(Ns,Ns)            :: Hk,tmp,Hk_bare,Hstar
  real(8),dimension(Ns)                      :: ek
  integer                                           :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  Estar=0.d0
  slater_derivatives_=0.d0
  do ik=1,Lk
     !
     Hk=0.d0
     ek=0.d0
     !
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
  if(present(slater_derivatives)) slater_derivatives=slater_derivatives_
  !<DEBUG
  ! write(*,*) 'slater derivatives probably wrong!'
  ! do istate=1,Ns
  !    write(*,'(6F18.10)') slater_derivatives(istate,1:Ns)
  ! end do
  !DEBUG>
end subroutine store_slater_ground_state





subroutine store_slater_ground_state_cmin(Rhop,lm,Estar,slater_matrix_el)     
  complex(8),dimension(Ns,Ns),intent(in) :: Rhop
  real(8),dimension(Ns,Ns),intent(in)           :: lm
  real(8)                                           :: Estar
  complex(8),dimension(Ns,Ns,Lk)            :: slater_matrix_el
  complex(8),dimension(Ns,Ns)            :: Hk,tmp,Hk_bare,Hstar
  real(8),dimension(Ns)                      :: ek
  integer                                           :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  Estar=0.d0
  slater_matrix_el=0.d0
  do ik=1,Lk
     Hk=0.d0
     ek=0.d0
     !
     Hk_bare=Hk_tb(:,:,ik)
     ! hopping renormalization !
     Hk=matmul(Hk_tb(:,:,ik),Rhop)
     Hk=matmul(Rhop,Hk)
     Hstar=Hk
     ! add Lagrange multipliers !
     Hk = Hk + lm
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
                 slater_matrix_el(istate,jstate,ik)=0.d0
                 do kstate=1,Ns
                    Estar = Estar + Hk(istate,kstate)*Hk(jstate,kstate)*Hstar(istate,jstate)*fermi(ek(kstate),beta)*wtk(ik)
                    slater_matrix_el(istate,jstate,ik) = slater_matrix_el(istate,jstate,ik) + &
                         Hk(istate,kstate)*Hk(jstate,kstate)*fermi(ek(kstate),beta)

                    ! Estar = Estar + Hk(istate,kstate)*Hk(jstate,kstate)*Hstar(istate,jstate)*heaviside(-1.d0*ek(kstate))*wtk(ik)
                    ! slater_matrix_el(istate,jstate,ik) = slater_matrix_el(istate,jstate,ik) + &
                    !      Hk(istate,kstate)*Hk(jstate,kstate)*heaviside(-1.d0*ek(kstate))

                 end do
              end do
           end do
        end do
     end do
  end do
end subroutine store_slater_ground_state_cmin
