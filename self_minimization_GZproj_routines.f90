!+- consider only the diagonal matrix lm(istate,istate) -+!
function get_delta_proj_variational_density_diag(lm_) result(delta_proj_variational_density_vec)
  real(8),dimension(:)                   :: lm_
  real(8),dimension(Ns)           :: delta_proj_variational_density_vec
  real(8),dimension(Ns,Ns) :: lm
  real(8),dimension(Ns,Ns) :: delta_proj_variational_density,proj_variational_density
  real(8),dimension(Nphi,Nphi)          :: H_projectors,H_tmp
  real(8),dimension(Nphi)               :: H_eigens

  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer :: iphi,jphi

  !<HARD DEBUGGING
  real(8),dimension(Ns)   :: tmp
  integer :: kstate

  complex(8),dimension(Nphi,Nphi,1) :: tmp_matrix
  complex(8),dimension(Nphi,Nphi) :: tmpV
  real(8),dimension(Nphi,1) :: tmp_eigens
  real(8) :: eps=1.d-10
  !HARD DEBUGGING>


  !
  lm=0.d0
  do istate=1,Ns
     lm(istate,istate) = lm_(istate)
  end do
  !
  proj_variational_density=0.d0
  !+- build up the local H_projectors -+!
  H_projectors=phi_traces_basis_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
        H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
     end do
  end do
  !
  
  ! tmp_matrix(:,:,1) = H_projectors
  ! call simultaneous_diag(tmp_matrix,tmpV,tmp_eigens,eps)

  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !H_projectors(:,1) = tmpV(:,Nphi)
  

  !
  do istate=1,Ns
     do jstate=1,Ns
        !
        proj_variational_density(istate,jstate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
        !
     end do
  end do
  do istate=1,Ns     
     delta_proj_variational_density_vec(istate) = proj_variational_density(istate,istate) - n0_target(istate)      
  end do
end function get_delta_proj_variational_density_diag
!
function get_proj_variational_density_diag(lm_) result(proj_variational_density)
  real(8),dimension(:)           :: lm_
  real(8),dimension(Ns)   :: proj_variational_density
  real(8),dimension(Ns,Ns)   :: lm
  real(8),dimension(Nphi,Nphi) :: H_projectors
  real(8),dimension(Nphi)       :: H_eigens
  integer                        :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer :: iphi,jphi

  !<HARD DEBUGGING
  real(8),dimension(Ns)   :: tmp
  !HARD DEBUGGING>

  !
  lm=0.d0
  do istate=1,Ns
     lm(istate,istate)=lm_(istate)
  end do
  !
  H_projectors=phi_traces_basis_Hloc   ! build up the local H_projectors
  proj_variational_density=0.d0  
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
        H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
     end do
  end do
  !
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !
  do istate=1,Ns
     !
     proj_variational_density(istate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,istate,:,:))
     !
  end do

  !<HARD DEBUGGING
  ! do istate=1,Ns
  !    !
  !    tmp(istate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,istate,:,:))
  !    !
  ! end do
  ! write(*,*) tmp(:)
  !HARD DEBUGGING>



end function get_proj_variational_density_diag



!+- Full Lagrange multiplier matrix (Ns x Ns) ...not really tested... -+!
function get_delta_proj_variational_density_full(lm_) result(delta_proj_variational_density_vec)
  real(8),dimension(:)                   :: lm_
  real(8),dimension(Ns*Ns) :: delta_proj_variational_density_vec
  real(8),dimension(Ns,Ns) :: lm
  real(8),dimension(Ns,Ns) :: delta_proj_variational_density,proj_variational_density
  real(8),dimension(Nphi,Nphi)         :: H_projectors
  real(8),dimension(Nphi)               :: H_eigens
  !
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer :: iphi,jphi
  !
  call vec2mat_stride(lm_,lm)
  proj_variational_density=0.d0
  !+- build up the local H_projectors -+!
  H_projectors=phi_traces_basis_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(istate)*(1.d0-n0_target(istate)))
        H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
     end do
  end do
  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  do istate=1,Ns
     do jstate=1,Ns
        !
        proj_variational_density(istate,jstate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))        
        !
     end do
  end do
  delta_proj_variational_density=proj_variational_density
  do istate=1,Ns
     delta_proj_variational_density(istate,istate) = delta_proj_variational_density(istate,istate) - n0_target(istate)      
  end do
  call mat2vec_stride(delta_proj_variational_density,delta_proj_variational_density_vec)
end function get_delta_proj_variational_density_full
!
function get_proj_variational_density_full(lm_) result(proj_variational_density)
  real(8),dimension(:)                   :: lm_
  real(8),dimension(Ns,Ns) :: proj_variational_density
  real(8),dimension(Ns,Ns) :: lm
  real(8),dimension(Nphi,Nphi)         :: H_projectors
  real(8),dimension(Nphi)               :: H_eigens
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer :: iphi,jphi
  !
  call vec2mat_stride(lm_,lm)
  proj_variational_density=0.d0
  !+- build up the local H_projectors -+!
  H_projectors=phi_traces_basis_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
        H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
     end do
  end do
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
end function get_proj_variational_density_full

!+- consider only the diagonal matrix lm(istate,istate) -+!
function get_delta_free_proj_variational_density_diag(lm_) result(delta_proj_variational_density_vec)
  real(8),dimension(:)                   :: lm_
  real(8),dimension(Ns)           :: delta_proj_variational_density_vec
  real(8),dimension(Ns,Ns) :: lm
  real(8),dimension(Ns,Ns) :: delta_proj_variational_density,proj_variational_density
  real(8),dimension(Nphi,Nphi)         :: H_projectors
  real(8),dimension(Nphi)               :: H_eigens
  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer :: iphi,jphi
  !
  lm=0.d0
  do istate=1,Ns
     lm(istate,istate) = lm_(istate)
  end do
  !
  proj_variational_density=0.d0
  !+- build up the local H_projectors -+!

  H_projectors=phi_traces_basis_free_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
        H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
     end do
  end do
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
  do istate=1,Ns     
     delta_proj_variational_density_vec(istate) = proj_variational_density(istate,istate) - n0_target(istate)      
  end do
end function get_delta_free_proj_variational_density_diag

















function proj_energy(lm_) result(E_Hloc)
  real(8),dimension(:),intent(in)                   :: lm_
  real(8)           :: E_Hloc
  real(8),dimension(Ns,Ns) :: lm
  real(8),dimension(Ns,Ns) :: delta_proj_variational_density,proj_variational_density
  real(8),dimension(Nphi,Nphi)          :: H_projectors,H_tmp
  real(8),dimension(Nphi)               :: H_eigens
  real(8),dimension(Nphi)               :: GZvect

  integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
  integer :: iphi,jphi

  !
  lm=0.d0
  do istate=1,Ns
     lm(istate,istate) = lm_(istate)
  end do
  !
  E_Hloc=0.d0
  !+- build up the local H_projectors -+!
  H_projectors=phi_traces_basis_Hloc
  do istate=1,Ns
     do jstate=1,Ns
        H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
        H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
     end do
  end do
  !

  !write(*,*) lm_
  ! do istate=1,Nphi
  !    write(*,'(20F6.2)') H_projectors(istate,:)
  ! end do


  call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
  !
  GZvect=H_projectors(1:Nphi,1)
  ! E_Hloc = trace_phi_basis(GZvect,phi_traces_basis_Hloc)

  ! do istate=1,Ns
  !    do jstate=1,Ns
  !       !
  !       proj_variational_density(istate,jstate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
  !       E_Hloc = E_Hloc + lm(istate,jstate)*proj_variational_density(istate,jstate)

  !       if(istate.eq.jstate) E_Hloc = E_Hloc - lm(istate,jstate)*n0_target(istate)
  !       !
  !    end do
  ! end do
  ! E_Hloc = -E_Hloc

  E_Hloc = 0.d0

  do istate=1,Ns
     do jstate=1,Ns
        !
        proj_variational_density(istate,jstate) = trace_phi_basis(H_projectors(:,1),phi_traces_basis_dens(istate,jstate,:,:))
        if(istate.eq.jstate) then
           E_Hloc = E_Hloc + abs(proj_variational_density(istate,jstate)-n0_target(istate))**2.d0
        else
           E_Hloc = E_Hloc + abs(proj_variational_density(istate,jstate))**2.d0
        end if

        !if(istate.eq.jstate) E_Hloc = E_Hloc - lm(istate,jstate)*n0_target(istate)
        !
     end do
  end do


end function proj_energy













