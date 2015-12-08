subroutine get_GZproj_ground_state(n0,slater_derivatives,lgr_multip,E_Hloc,GZvect) 
  real(8),dimension(state_dim)           :: n0
  real(8),dimension(state_dim,state_dim) :: slater_derivatives,lgr_multip
  real(8)                                :: E_Hloc
  real(8),dimension(nFock)               :: GZvect
  !
  real(8),dimension(nFock,nFock)         :: H_projectors
  real(8),dimension(nFock)               :: H_eigens
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  !
  !+- build up the local H_projectors -+!
  H_projectors=phi_traces_basis_Hloc
  do istate=1,state_dim
     do jstate=1,state_dim
        H_projectors = H_projectors + &
             slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0(jstate)*(1.d0-n0(jstate)))
        H_projectors = H_projectors + lgr_multip(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
     end do
  end do
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  GZvect=H_projectors(1:nFock,1)
  !
  E_Hloc=0.d0
  do ifock=1,nFock
     do jfock=1,nFock
        E_Hloc=E_Hloc+GZvect(ifock)*phi_traces_basis_Hloc(ifock,jfock)*GZvect(jfock)
     end do
  end do
end subroutine get_GZproj_ground_state



subroutine get_GZproj_free_ground_state(n0,slater_derivatives,lgr_multip,E_Hloc,GZvect) 
  real(8),dimension(state_dim)           :: n0
  real(8),dimension(state_dim,state_dim) :: slater_derivatives,lgr_multip
  real(8)                                :: E_Hloc
  real(8),dimension(nFock)               :: GZvect
  !
  real(8),dimension(nFock,nFock)         :: H_projectors
  real(8),dimension(nFock)               :: H_eigens
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  !
  !+- build up the local H_projectors -+!
  H_projectors=phi_traces_basis_free_Hloc
  do istate=1,state_dim
     do jstate=1,state_dim
        H_projectors = H_projectors + &
             slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0(jstate)*(1.d0-n0(jstate)))
        H_projectors = H_projectors + lgr_multip(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
     end do
  end do
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  GZvect=H_projectors(1:nFock,1)
  !
  E_Hloc=0.d0
  do ifock=1,nFock
     do jfock=1,nFock
        E_Hloc=E_Hloc+GZvect(ifock)*phi_traces_basis_Hloc(ifock,jfock)*GZvect(jfock)
     end do
  end do
end subroutine get_GZproj_free_ground_state





subroutine store_slater_ground_state(Rhop,lm,Estar,slater_derivatives)     
  real(8),dimension(state_dim,state_dim),intent(in) :: Rhop
  real(8),dimension(state_dim),intent(in)           :: lm
  real(8)                                           :: Estar
  real(8),dimension(state_dim,state_dim)            :: slater_derivatives
  real(8),dimension(state_dim,state_dim)            :: Hk,tmp,Hk_bare,Hstar
  real(8),dimension(state_dim)                      :: ek
  integer                                           :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  Estar=0.d0
  slater_derivatives=0.d0
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
     do istate=1,state_dim
        Hk(istate,istate)=Hk(istate,istate)+lm(istate)             
     end do
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
                 do kstate=1,state_dim
                    tmp(istate,jstate) = tmp(istate,jstate) + Hk(istate,kstate)*fermi(ek(kstate),beta)*Hk(jstate,kstate)
                    Estar = Estar + Hk(istate,kstate)*Hk(jstate,kstate)*Hstar(istate,jstate)*fermi(ek(kstate),beta)*wtk(ik)
                 end do
              end do
           end do
        end do
     end do
     ! store slater ground state derivatives
     tmp=matmul(Rhop,tmp)
     tmp=matmul(Hk_bare,tmp)             
     do istate=1,state_dim
        do jstate=1,state_dim             
           slater_derivatives(istate,jstate) = &
                slater_derivatives(istate,jstate) + 2.d0*tmp(istate,jstate)*wtk(ik)
        end do
     end do
  end do
end subroutine store_slater_ground_state




subroutine store_slater_ground_state_cmin(Rhop,lm,Estar,slater_matrix_el)     
  real(8),dimension(state_dim,state_dim),intent(in) :: Rhop
  real(8),dimension(state_dim),intent(in)           :: lm
  real(8)                                           :: Estar
  real(8),dimension(state_dim,state_dim,Lk)            :: slater_matrix_el
  real(8),dimension(state_dim,state_dim)            :: Hk,tmp,Hk_bare,Hstar
  real(8),dimension(state_dim)                      :: ek
  integer                                           :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  Estar=0.d0
  slater_matrix_el=0.d0
  do ik=1,Lk
     Hk=0.d0
     ek=0.d0
     ! do iorb=1,Norb
     !    do ispin=1,2
     !       do jorb=1,Norb
     !          do jspin=1,2
     !             istate=index(ispin,iorb)
     !             jstate=index(jspin,jorb)               
     !             ! build up the hopping hamiltonian !
     !             if(ispin.eq.jspin) then
     !                if(iorb.eq.jorb) then
     !                   Hk(istate,jstate)=epsik(ik)
     !                else
     !                   Hk(istate,jstate)=hybik(ik)
     !                end if
     !             end if
     !          end do
     !       end do
     !    end do
     ! end do
     Hk_bare=Hk_tb(:,:,ik)
     ! hopping renormalization !
     Hk=matmul(Hk_tb(:,:,ik),Rhop)
     Hk=matmul(Rhop,Hk)
     Hstar=Hk
     ! add Lagrange multipliers !
     do istate=1,state_dim
        Hk(istate,istate)=Hk(istate,istate)+lm(istate)             
     end do
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
                 do kstate=1,state_dim
                    Estar = Estar + Hk(istate,kstate)*Hk(jstate,kstate)*Hstar(istate,jstate)*fermi(ek(kstate),beta)*wtk(ik)
                    slater_matrix_el(istate,jstate,ik) = slater_matrix_el(istate,jstate,ik) + &
                         Hk(istate,kstate)*Hk(jstate,kstate)*fermi(ek(kstate),beta)
                 end do
              end do
           end do
        end do
     end do
  end do
end subroutine store_slater_ground_state_cmin



