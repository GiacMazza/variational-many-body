function constraints_deviation(lm_) result(delta_constraints)
  real(8),dimension(:)            :: lm_
  real(8),dimension(2*Ns)         :: delta_constraints
  real(8),dimension(Ns,Ns)        :: lm_dens,lm_Rhop
  real(8),dimension(Ns,Ns)        :: delta_constraints_dens,proj_variational_density
  real(8),dimension(Ns,Ns)        :: delta_constraints_Rhop
  complex(8),dimension(Ns,Ns)        :: proj_Rhop
  !
  complex(8),dimension(Nphi,Nphi) :: H_projectors,H_tmp
  real(8),dimension(Nphi)         :: H_eigens
  real(8)                         :: tmp_ene
  !
  integer                         :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                         :: iphi,jphi
  !

  !
  lm_dens=0.d0
  lm_Rhop=0.d0
  do istate=1,Ns
     lm_dens(istate,istate) = lm_(istate)
     lm_Rhop(istate,istate) = lm_(istate+Ns)
  end do
  write(*,*) lm_(1:Ns)
  write(*,*) lm_(Ns+1:2*Ns)

  !
  !+- build up the local H_projectors -+!
  H_projectors=zero
  H_projectors=phi_traces_basis_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
        H_projectors = H_projectors + lm_dens(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
        ! H_projectors = H_projectors + 0.5d0*lm_Rhop(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)
        ! H_projectors = H_projectors + 0.5d0*conjg(lm_Rhop(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:))
     end do
  end do
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !

  proj_variational_density=0.d0
  proj_Rhop=zero
  !
  do istate=1,Ns
     do jstate=1,Ns
        proj_variational_density(istate,jstate) = &
             trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
        proj_Rhop(istate,jstate) = &
             trace_phi_basis(H_projectors(:,1),phi_traces_basis_Rhop(istate,jstate,:,:))
        !write(*,*) proj_variational_density(istate,jstate)
     end do
  end do
  !
  delta_constraints=0.d0
  do istate=1,Ns     
     delta_constraints(istate) = proj_variational_density(istate,istate) - n0_target(istate)      
     !delta_constraints(istate+Ns) = abs(proj_Rhop(istate,istate) - R_target(istate,istate)*sqrt(n0_target(istate)*(1.d0-n0_target(istate))))
     !write(*,*) delta_constraints(istate+Ns)
  end do
  !
  write(*,*) 
  do is=1,2*Ns
     write(*,*) delta_constraints(is),lm_(is)
  end do
  write(*,*) 
  !stop
  !
end function constraints_deviation




function constraints_deviation_tmp(lm_) result(delta_constraints)
  real(8),dimension(:)            :: lm_
  real(8),dimension(Ns)           :: delta_constraints
  real(8),dimension(Ns,Ns)        :: lm_dens,lm_Rhop
  real(8),dimension(Ns,Ns)        :: delta_constraints_dens,proj_variational_density
  real(8),dimension(Ns,Ns)        :: delta_constraints_Rhop
  complex(8),dimension(Ns,Ns)     :: proj_Rhop
  !
  complex(8),dimension(Nphi,Nphi) :: H_projectors,H_tmp
  real(8),dimension(Nphi)         :: H_eigens
  real(8)                         :: tmp_ene
  !
  integer                         :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                         :: iphi,jphi
  !

  !
  lm_dens=0.d0
  lm_Rhop=0.d0
  do istate=1,Ns
     lm_dens(istate,istate) = lm_(1)
!     lm_Rhop(istate,istate) = lm_(istate+Ns)
  end do
  !write(*,*) lm_(1:Ns)
  !write(*,*) lm_(Ns+1:2*Ns)

  !
  !+- build up the local H_projectors -+!
  H_projectors=zero
  H_projectors=phi_traces_basis_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
        H_projectors = H_projectors + lm_dens(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
        ! H_projectors = H_projectors + 0.5d0*lm_Rhop(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)
        ! H_projectors = H_projectors + 0.5d0*conjg(lm_Rhop(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:))
     end do
  end do
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !

  proj_variational_density=0.d0
  proj_Rhop=zero
  !
  do istate=1,Ns
     do jstate=1,Ns
        proj_variational_density(istate,jstate) = &
             trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
        proj_Rhop(istate,jstate) = &
             trace_phi_basis(H_projectors(:,1),phi_traces_basis_Rhop(istate,jstate,:,:))
        !write(*,*) proj_variational_density(istate,jstate)
     end do
  end do
  !
  delta_constraints=0.d0
  do istate=1,Ns     
     delta_constraints(istate) = proj_variational_density(istate,istate) - n0_target(istate)      
     !delta_constraints(istate+Ns) = abs(proj_Rhop(istate,istate) - R_target(istate,istate)*sqrt(n0_target(istate)*(1.d0-n0_target(istate))))
     !write(*,*) delta_constraints(istate+Ns)
  end do
  !
  write(*,*) 
  do is=1,Ns
     write(*,*) delta_constraints(is),lm_(1)
  end do
  write(*,*) 
  !stop
  !
end function constraints_deviation_tmp














function constraints_deviation_amb(lm_) result(delta_constraints)
  real(8),dimension(:),intent(in)            :: lm_
  real(8)         :: delta_constraints
  real(8),dimension(Ns,Ns)        :: lm_dens,lm_Rhop
  real(8),dimension(Ns,Ns)        :: delta_constraints_dens,proj_variational_density
  real(8),dimension(Ns,Ns)        :: delta_constraints_Rhop
  complex(8),dimension(Ns,Ns)        :: proj_Rhop,Rtmp
  !
  complex(8),dimension(Nphi,Nphi) :: H_projectors,H_tmp
  real(8),dimension(Nphi)         :: H_eigens
  real(8)                         :: tmp_ene
  !
  integer                         :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer                         :: iphi,jphi
  !

  !
  lm_dens=0.d0
  lm_Rhop=0.d0
  do istate=1,Ns
     lm_dens(istate,istate) = lm_(istate)
     !lm_Rhop(istate,istate) = lm_(istate+Ns)
  end do
  !write(*,*) lm_(Ns+1:2*Ns)
  !Rtmp=zero
  


  !
  !+- build up the local H_projectors -+!
  H_projectors=zero
  H_projectors=phi_traces_basis_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
        H_projectors = H_projectors + lm_dens(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
        !H_projectors = H_projectors + lm_Rhop(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)
     end do
  end do
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !

  proj_variational_density=0.d0
  proj_Rhop=zero
  do istate=1,Ns
     do jstate=1,Ns
        proj_variational_density(istate,jstate) = &
             trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
        !proj_Rhop(istate,jstate) = &
        !             trace_phi_basis(H_projectors(:,1),phi_traces_basis_Rhop(istate,jstate,:,:))
        !write(*,*) proj_variational_density(istate,jstate)
     end do
  end do
  !
  delta_constraints=0.d0
  do istate=1,Ns
     do jstate=1,Ns
        if(istate.eq.jstate) then
           delta_constraints = delta_constraints + &
                abs(proj_variational_density(istate,jstate)-n0_target(istate))**2.d0
        else
           delta_constraints = delta_constraints + &
                abs(proj_variational_density(istate,jstate))**2.d0
        end if
        !delta_constraints = delta_constraints + &
!             abs(proj_Rhop(istate,jstate)-R_target(istate,jstate)*sqrt((n0_target(jstate))*(1.d0-n0_target(jstate))))**2.d0
   
        !if(istate.ne.jstate) write(*,*) proj_Rhop(istate,jstate)
     end do
  end do

  !
  ! write(*,*) 
  ! do is=1,2*Ns
  !write(*,*) delta_constraints,lm_(1:Ns)
  !
end function constraints_deviation_amb









