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






subroutine store_slater_ground_state(Rhop,lm,Estar,slater_derivatives)     
  complex(8),dimension(Ns,Ns),intent(in) :: Rhop
  real(8),dimension(Ns,Ns),intent(in)           :: lm
  real(8)                                           :: Estar
  complex(8),dimension(Ns,Ns)            :: slater_derivatives
  complex(8),dimension(Ns,Ns)            :: Hk,tmp,Hk_bare,Hstar
  real(8),dimension(Ns)                      :: ek
  integer                                           :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
  Estar=0.d0
  slater_derivatives=0.d0
  !<DEBUG
  write(*,*) size(lm,1),size(lm,2),size(Hk,1),size(Hk,2)
  write(*,*) Rhop
  !DEBUG>
  
  do ik=1,Lk
     !write(*,*) ik     
     Hk=0.d0
     ek=0.d0
     !
     Hk_bare=Hk_tb(:,:,ik)
     ! hopping renormalization !
     Hk=matmul(Hk_tb(:,:,ik),Rhop)
     Hk=matmul(Rhop,Hk)
     Hstar=Hk
     !<DEBUG     
     ! write(*,*) 'IK',ik
     ! write(*,*) Rhop
     ! write(*,*) 'IK',ik
     !DEBUG>
     ! add Lagrange multipliers !
     
     !do istate=1,Ns
     !        Hk(istate,istate)=Hk(istate,istate)+lm(istate,istate)             
     !end do

!     Hk=Hk+lm

     ! diagonalize hamiltonian !
     call  matrix_diagonalize(Hk,ek,'V','L')

     !<DEBUG     
     write(*,*) 'IK',ik
     !DEBUG>

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

                    ! tmp(istate,jstate) = tmp(istate,jstate) + Hk(istate,kstate)*heaviside(-1.d0*ek(kstate))*conjg(Hk(jstate,kstate))
                    ! Estar = Estar + Hk(istate,kstate)*Hk(jstate,kstate)*Hstar(istate,jstate)*heaviside(-1.d0*ek(kstate))*wtk(ik)

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
           slater_derivatives(istate,jstate) = &
                slater_derivatives(istate,jstate) + 2.d0*tmp(istate,jstate)*wtk(ik)
        end do
     end do
  end do

  !<DEBUG
  ! write(*,*) 'slater derivatives probably wrong!'
  ! do istate=1,Ns
  !    write(*,'(6F18.10)') slater_derivatives(istate,1:Ns)
  ! end do
  !DEBUG>

end subroutine store_slater_ground_state




subroutine store_slater_ground_state_cmin(Rhop,lm,Estar,slater_matrix_el)     
  complex(8),dimension(Ns,Ns),intent(in) :: Rhop
  real(8),dimension(Ns),intent(in)           :: lm
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
     do istate=1,Ns
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



