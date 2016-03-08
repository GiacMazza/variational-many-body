function get_delta_local_density_matrix_diag(lm_) result(delta_local_density_matrix_vec)
  real(8),dimension(:)  :: lm_
  real(8),dimension(state_dim)  :: delta_local_density_matrix_vec
  real(8),dimension(state_dim*state_dim)  :: delta_local_density_matrix_vec_
  real(8),dimension(state_dim,state_dim)  :: lm
  real(8),dimension(state_dim,state_dim)  :: delta_local_density_matrix,local_density_matrix
  real(8),dimension(state_dim,state_dim) :: Hk,tmp
  real(8),dimension(state_dim)          :: ek
  integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  !
  lm=0.d0
  do istate=1,state_dim
     lm(istate,istate)=lm_(istate)
  end do
  local_density_matrix=0.d0
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
     do istate=1,state_dim
        do jstate=1,state_dim
           do kstate=1,state_dim
              local_density_matrix(istate,jstate) = &
                   local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*Hk(istate,kstate)*Hk(jstate,kstate)*wtk(ik)
           end do
        end do
     end do
  end do
  ! return variation of local density matrix with respect to the target values
  delta_local_density_matrix = local_density_matrix
  do istate=1,state_dim
     delta_local_density_matrix(istate,istate) = delta_local_density_matrix(istate,istate) - n0_target(istate)      
  end do
  do istate=1,state_dim
     delta_local_density_matrix_vec(istate)=delta_local_density_matrix(istate,istate)
  end do
end function get_delta_local_density_matrix_diag
!
function get_delta_local_density_matrix_full(lm_) result(delta_local_density_matrix_vec)
  real(8),dimension(:)  :: lm_
  real(8),dimension(state_dim,state_dim)  :: lm_full
  real(8),dimension(state_dim*state_dim)  :: delta_local_density_matrix_vec
  real(8),dimension(state_dim,state_dim)  :: lm
  real(8),dimension(state_dim,state_dim)  :: delta_local_density_matrix,local_density_matrix
  real(8),dimension(state_dim,state_dim) :: Hk,tmp
  real(8),dimension(state_dim)          :: ek
  integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  !
  call vec2mat_stride(lm_,lm)
  local_density_matrix=0.d0
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
     do istate=1,state_dim
        do jstate=1,state_dim
           do kstate=1,state_dim
              local_density_matrix(istate,jstate) = &
                   local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*Hk(istate,kstate)*Hk(jstate,kstate)*wtk(ik)                        
           end do
        end do
     end do
  end do
  ! return variation of local density matrix with respect to the target values
  delta_local_density_matrix = local_density_matrix
  do istate=1,state_dim
     delta_local_density_matrix(istate,istate) = delta_local_density_matrix(istate,istate) - n0_target(istate)      
  end do
  call mat2vec_stride(delta_local_density_matrix,delta_local_density_matrix_vec)
end function get_delta_local_density_matrix_full
!



function get_local_density_matrix_diag(lm_) result(local_density_matrix)
  real(8),dimension(state_dim) :: lm_
  real(8),dimension(state_dim,state_dim) :: lm
  real(8),dimension(state_dim,state_dim) :: local_density_matrix
  real(8),dimension(state_dim,state_dim) :: Hk
  real(8),dimension(state_dim)           :: ek
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  lm=0.d0
  do istate=1,state_dim
     lm(istate,istate)=lm_(istate)
  end do
  local_density_matrix=0.d0
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
     do iorb=1,Norb
        do ispin=1,2
           do jorb=1,Norb
              do jspin=1,2
                 istate=index(ispin,iorb)
                 jstate=index(jspin,jorb)               
                 do kstate=1,state_dim
                    local_density_matrix(istate,jstate) = &
                         local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*Hk(istate,kstate)*Hk(jstate,kstate)*wtk(ik)                        
                 end do
              end do
           end do
        end do
     end do
  end do
end function get_local_density_matrix_diag
!
function get_local_density_matrix_full(lm) result(local_density_matrix)
  real(8),dimension(state_dim,state_dim)  :: lm
  real(8),dimension(state_dim,state_dim)  :: local_density_matrix
  real(8),dimension(state_dim,state_dim) :: Hk
  real(8),dimension(state_dim)          :: ek
  integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  local_density_matrix=0.d0
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
     do iorb=1,Norb
        do ispin=1,2
           do jorb=1,Norb
              do jspin=1,2
                 istate=index(ispin,iorb)
                 jstate=index(jspin,jorb)               
                 do kstate=1,state_dim
                    local_density_matrix(istate,jstate) = &
                         local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*Hk(istate,kstate)*Hk(jstate,kstate)*wtk(ik)                        
                 end do
              end do
           end do
        end do
     end do
  end do
end function get_local_density_matrix_full
!
!ORBITAL DIAGONAL TEST
!
function get_delta_local_density_matrix_orb(lm_) result(delta_local_density_matrix_vec)
  real(8),dimension(:)  :: lm_
  real(8),dimension(Norb)  :: delta_local_density_matrix_vec
  real(8),dimension(state_dim*state_dim)  :: delta_local_density_matrix_vec_
  real(8),dimension(state_dim,state_dim)  :: lm
  real(8),dimension(state_dim,state_dim)  :: delta_local_density_matrix,local_density_matrix
  real(8),dimension(state_dim,state_dim) :: Hk,tmp
  real(8),dimension(state_dim)          :: ek
  integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  !
  lm=0.d0
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        lm(istate,istate)= lm_(iorb)
     end do
  end do
  !
  local_density_matrix=0.d0
  do ik=1,Lk
     Hk=0.d0
     ek=0.d0
     ! hopping renormalization !
     Hk=matmul(Hk_tb(:,:,ik),Rhop)
     Hk=matmul(Rhop,Hk)
     ! add Lagrange multipliers !
     Hk=Hk-lm                     
     ! diagonalize hamiltonian !
     call  matrix_diagonalize(Hk,ek,'V','L')
     !compute local density matrix
     do istate=1,state_dim
        do jstate=1,state_dim
           do kstate=1,state_dim
              local_density_matrix(istate,jstate) = &
                   local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*Hk(istate,kstate)*Hk(jstate,kstate)*wtk(ik)
           end do
        end do
     end do
  end do
  ! return variation of local density matrix with respect to the target values
  delta_local_density_matrix = local_density_matrix
  do istate=1,state_dim
     delta_local_density_matrix(istate,istate) = delta_local_density_matrix(istate,istate) - n0_target(istate)      
  end do
  do iorb=1,Norb
     delta_local_density_matrix_vec(iorb)=0.d0
     do ispin=1,2
        istate=index(ispin,iorb)
        delta_local_density_matrix_vec(iorb)=&
             delta_local_density_matrix_vec(iorb)+delta_local_density_matrix(istate,istate)
     end do
  end do
end function get_delta_local_density_matrix_orb







