MODULE GZ_TWO_SITES_FOCK
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE SF_IOTOOLS
  implicit none
  private
  !  
  public :: initialize_two_sites_fock_space       
  !public :: build_local_hamiltonian
  !
CONTAINS
  !
  include 'two_sites_fock_space_algebra.f90'
  include 'local_fock_space_observables.f90'
  !
  subroutine initialize_two_sites_fock_space
    integer :: iorb,jorb,ispin,jspin,ifock,jfock,isite,ii
    integer,dimension(:),allocatable   :: Fock,ivec
    !
    Ns = 2*Norb*Nsite
    NFock = 2**Ns
    !
    allocate(ivec(Ns))
    allocate(Fock(Nfock))
    do ifock=1,NFock
       Fock(ifock)=ifock
       call bdecomp(Fock(ifock),ivec)
       call get_state_number(ivec,jfock)
       if(jfock/=ifock) stop "error initialization Fock space"
       write(*,*) '|',ivec(1:Norb),',',ivec(1+Norb:2*Norb),',',ivec(1+2*Norb:3*Norb),',',ivec(1+3*Norb:4*Norb),'>'
    end do
    !
    !+- Allocate and initialize stride -+! 
    allocate(i_ios(2,Norb,Nsite))
    do ispin=1,2
       ii = (ispin-1)*Norb*Nsite
       do iorb=1,Norb
          do isite=1,Nsite
             i_ios(ispin,iorb,isite)=ii+iorb+(isite-1)*Nsite
          enddo
       enddo
    end do
    !
    call build_two_sites_fock_algebra
    !
    !+- here restrict the Hilbert space to the two-particle sector-+
    ! call build_local_hamiltonian
    ! !
    ! call build_local_observables
    !
    call fock_to_2p_hilbert
    !
  end subroutine initialize_two_sites_fock_space

  subroutine fock_to_2p_hilbert
    integer :: iip,ifock,in
    integer,dimension(:),allocatable   :: ivec

    Nh_2p = binomial(Ns,2)
    write(*,*) Nh_2p    
    allocate(ivec(Ns))
    allocate(ifk_to_i2p(nFock),i2p_to_ifk(Nh_2p))
    !
    ifk_to_i2p=0
    iip=0    
    do ifock=1,nFock
       call bdecomp(ifock,ivec)
       in=sum(ivec(1:Ns))       
       if(in.eq.2) then
          iip=iip+1
          i2p_to_ifk(iip) = ifock
          ifk_to_i2p(ifock) = iip
       end if
    end do
    write(*,*) i2p_to_ifk
    write(*,*) ifk_to_i2p    
    do iip=1,Nh_2p       
       ifock=i2p_to_ifk(iip)       
       call bdecomp(ifock,ivec)
       write(*,*) '|',ivec(1:Norb),',',ivec(1+Norb:2*Norb),',',ivec(1+2*Norb:3*Norb),',',ivec(1+3*Norb:4*Norb),'>'       
    end do

    allocate(CC_(Ns,nh_2p,nh_2p),CA_(Ns,nh_2p,nh_2p))
    do is=1,Ns
       do iip=1,nh_2p
          do jjp=1,nh_2p
             ifock=i2p_to_ifk(iip)
             jfock=i2p_to_ifk(jjp)
             CC_(is,iip,jjp) = CC(is,ifock,jfock)
             CA_(is,iip,jjp) = CA(is,ifock,jfock)
          end do
       end do
    end do
    
  end subroutine fock_to_2p_hilbert


  ! !
  subroutine build_two_sites_hamiltonian(elevels)
    integer :: iorb,jorb,ispin,jspin,istate,jstate,is_up,is_dn,js_up,js_dn,ifock,is
    real(8),dimension(Ns) :: elevels
    complex(8),dimension(nh_2p,nh_2p) :: tmpH
    integer :: iip,jjp
    !
    allocate(two_sites_Hamiltonian(nh_2p,nh_2p),tmpH(nh_2p,nh_2p))    
    do is=1,Ns
       two_sites_Hamiltonian=two_sites_Hamiltonian+elevels(is)*matmul(CC_(is,:,:),CA_(is,:,:)) 
    end do
    !
    tmph=0.d0
    do isite=1,Nsite
       do iorb=1,Norb
          is=i_ios(1,iorb,isite)
          js=i_ios(2,iorb,isite)          
          tmph=matmul(CC_(is,:,:),CA_(is,:,:))
          tmph=matmul(CA_(js,:,:),tmph)
          tmph=matmul(CC_(js,:,:),tmph)
          two_sites_Hamiltonian=two_sites_Hamiltonian+U*tmpH
       end do
    end do

    tmph=0.d0
    do ispin=1,2
       do jspin=1,2
          do iorb=1,Norb
             do jorb=1,Norb
                is=i_ios(ispin,iorb,1)
                js=i_ios(jspin,jorb,2)                          
                tmph=matmul(CC_(is,:,:),CA_(is,:,:))
                tmph=matmul(CA_(js,:,:),tmph)
                tmph=matmul(CC_(js,:,:),tmph)
             end do
          end do
       end do
    end do
    
  end subroutine build_two_sites_hamiltonian
  !   !+- energy of the atomic levels -+!
  !   allocate(atomic_energy_levels(Ns))
  !   select case(Norb)
  !   case(2)
  !      do iorb=1,Norb
  !         do ispin=1,2
  !            istate=index(ispin,iorb)
  !            atomic_energy_levels(istate) = Cfield*0.5d0
  !            if(iorb.eq.2) atomic_energy_levels(istate) = -Cfield*0.5d0
  !         end do
  !      end do
  !   case default
  !      atomic_energy_levels=0.d0        
  !   end select


  !   do istate=1,Ns
  !      state_dens(istate,:,:) = matmul(cc(istate,:,:),ca(istate,:,:))
  !   end do
  !   !+- FREE PART OF THE LOCAL HAMILTONIAN -+!
  !   local_hamiltonian=0.d0
  !   do istate=1,Ns
  !      local_hamiltonian = local_hamiltonian + atomic_energy_levels(istate)*state_dens(istate,:,:)
  !      local_hamiltonian = local_hamiltonian - xmu*state_dens(istate,:,:)
  !   end do
  !   local_hamiltonian_free=local_hamiltonian
  !   !+- LOCAL INTERACTION HAMILTONIAN -+!

  !   !+- INTRA-ORBITAL density-density  -+!
  !   do iorb=1,Norb
  !      is_up = index(1,iorb)
  !      is_dn = index(2,iorb)
  !      local_hamiltonian = local_hamiltonian + &
  !           Uloc(iorb)*matmul(state_dens(is_up,:,:),state_dens(is_dn,:,:))
  !   end do
    
  !   if(Norb.gt.1) then
  !      !+- INTER-ORBITAL density-density opposite spins -+!
  !      do iorb=1,Norb
  !         do jorb=iorb+1,Norb
  !            !
  !            is_up = index(1,iorb)
  !            is_dn = index(2,iorb)
  !            !
  !            js_up = index(1,jorb)
  !            js_dn = index(2,jorb)
  !            !
  !            local_hamiltonian = local_hamiltonian + &
  !                 Ust*matmul(state_dens(is_up,:,:),state_dens(js_dn,:,:))+ &
  !                 Ust*matmul(state_dens(is_dn,:,:),state_dens(js_up,:,:))
  !         end do
  !      end do
  !      !+- INTER-ORBITAL density-density parallel spins -+!
  !      do iorb=1,Norb
  !         do jorb=iorb+1,Norb
  !            !
  !            is_up = index(1,iorb)
  !            is_dn = index(2,iorb)
  !            !
  !            js_up = index(1,jorb)
  !            js_dn = index(2,jorb)
  !            !             
  !            local_hamiltonian = local_hamiltonian + &
  !                 (Ust-Jh)*matmul(state_dens(is_up,:,:),state_dens(js_up,:,:)) + &
  !                 (Ust-Jh)*matmul(state_dens(is_dn,:,:),state_dens(js_dn,:,:))
  !         end do
  !      end do
  !      !+- SPIN-FLIP interaction -+!
  !      do iorb=1,Norb
  !         do jorb=1,Norb
  !            !
  !            is_up = index(1,iorb)
  !            is_dn = index(2,iorb)
  !            !
  !            js_up = index(1,jorb)
  !            js_dn = index(2,jorb)
  !            !             
  !            if(iorb.ne.jorb) then
  !               tmp=matmul(CC(js_dn,:,:),CA(js_up,:,:))
  !               tmp=matmul(CA(is_dn,:,:),tmp)
  !               tmp=matmul(CC(is_up,:,:),tmp)
  !               local_hamiltonian = local_hamiltonian - Jsf*tmp
  !            end if
  !         end do
  !      end do
  !      !+- PAIR-HOPPING interaction -+!
  !      do iorb=1,Norb
  !         do jorb=1,Norb
  !            !
  !            is_up = index(1,iorb)
  !            is_dn = index(2,iorb)
  !            !
  !            js_up = index(1,jorb)
  !            js_dn = index(2,jorb)
  !            !             
  !            if(iorb.ne.jorb) then
  !               tmp=matmul(CA(js_dn,:,:),CA(js_up,:,:))
  !               tmp=matmul(CC(is_dn,:,:),tmp)
  !               tmp=matmul(CC(is_up,:,:),tmp)
  !               local_hamiltonian = local_hamiltonian + Jph*tmp
  !            end if
  !         end do
  !      end do
  !   end if
  !   !+- CHEMICAL POTENTIAL FOR PH CONDITION -+!
  !   !mu_ph = Uloc(1)*0.5d0 + dble(Norb-1)*0.5d0*(2.d0*Ust-Jh)
  !   do iorb=1,Norb
  !      do ispin=1,2
  !         is=index(ispin,iorb)
  !         mu_ph_(is) = Uloc(iorb)*0.5d0 + dble(Norb-1)*0.5d0*(2.d0*Ust-Jh)
  !      end do
  !   end do
  !   do istate=1,Ns       
  !      local_hamiltonian = local_hamiltonian - mu_ph_(is)*state_dens(istate,:,:)
  !   end do
  !   !
  ! end subroutine build_local_hamiltonian
  ! !

  ! subroutine build_local_observables  !+---> forse piu' corretto chiamarli local operators...

  !   allocate(op_dens(Ns,nFock,nFock))
  !   allocate(op_local_dens(Ns,Ns,nFock,nFock))
    
  !   !allocate(op_docc(Norb,nFock,nFock)) !+- TO BE REM
  !   !allocate(op_dens_dens_orb(Norb,Norb,nFock,nFock)) !!+- TO BE REM
  !   allocate(op_dens_dens(Ns,Ns,nFock,nFock))
  !   allocate(op_spin_flip(Norb,Norb,nFock,nFock))
  !   allocate(op_pair_hopping(Norb,Norb,nFock,nFock))
  !   allocate(op_sc_order(Ns,Ns,nFock,nFock))

  !   allocate(op_local_dens_anomalous(Ns,Ns,nFock,nFock))
    
  !   allocate(op_spin2(nFock,nFock))
  !   allocate(op_spinZ(nFock,nFock))
  !   allocate(op_isospin2(nFock,nFock))
  !   allocate(op_isospinZ(nFock,nFock))
    


  !   !    Uhubbard=density_density_interaction(CC,CA)    
  !   !dens_dens_interaction=rotationally_invariant_density_density(CC,CA)  !HERE MAY ADD SINGLET SPLITTING TERMS, SPIN FLIPS, PAIR HOPPINGS, etc...
  !   !op_docc                 = local_doubly(CC,CA)
    
  !   !op_dens                 = local_density(CC,CA) !+- TO BE REMOVED
    
  !   op_local_dens           = local_density_matrix(CC,CA)
  !   op_local_dens_anomalous = local_density_matrix_anomalous(CC,CA)
  !   !
  !   !op_dens_dens_orb        = local_density_density_orb(CC,CA) !+- TO BE REMOVED    
  !   op_dens_dens            = local_density_density(CC,CA)
  !   op_spin_flip            = local_spin_flip(CC,CA)
  !   op_pair_hopping         = local_pair_hopping(CC,CA)
  !   op_sc_order             = local_sc_order(CC,CA)
  !   !
  !   op_spin2                = local_spin2(CC,CA)
  !   op_spinZ                = local_spinZ(CC,CA)
  !   op_isospin2             = local_isospin2(CC,CA)
  !   op_isospinZ             = local_isospinZ(CC,CA)
  !   !
  ! end subroutine build_local_observables



  ! !


  ! function sz_rotate(fock_in) result(fock_out)
  !   integer                      :: fock_in,fock_out
  !   integer,dimension(Ns) :: state_in,state_out    
  !   integer                      :: iorb,istate
  !   call bdecomp(fock_in,state_in)    
  !   do iorb=1,Norb
  !      state_out(iorb) = state_in(iorb+Norb)
  !      state_out(iorb+Norb) = state_in(iorb)
  !   end do
  !   fock_out=1
  !   do istate=0,Ns-1
  !      fock_out = fock_out + state_out(istate+1)*2**istate
  !   end do
  ! end function sz_rotate


  ! ! subroutine get_spin_indep_states
  ! !   integer :: i_ind,i,iorb,istate
  ! !   integer :: tmp_search(nFock),tmp_target(nFock)
  ! !   integer :: ifock,isymm
  ! !   integer :: check_maps
  ! !   integer :: test_vec(Ns)
  ! !   !+- get independent states under sz rotation symmetry -+!        
  ! !   tmp_search=0
  ! !   i_ind=0
  ! !   do ifock=1,nFock       
  ! !      tmp_target(ifock)=sz_rotate(ifock)              
  ! !      if(tmp_search(ifock).ge.0) then
  ! !         i_ind=i_ind+1
  ! !         tmp_search(ifock)=ifock
  ! !         if(tmp_target(ifock).ne.ifock) tmp_search(tmp_target(ifock)) = -1
  ! !      end if
  ! !   end do
  ! !   nFock_indep=i_ind
  ! !   allocate(fock_indep(nFock_indep),full2indep_fock(nFock),indep2full_fock(nFock_indep,2))
  ! !   i_ind=0
  ! !   do i=1,nFock
  ! !      if(tmp_search(i).ge.0) then
  ! !         i_ind=i_ind+1
  ! !         fock_indep(i_ind) = tmp_search(i)
  ! !      end if
  ! !   end do
  ! !   do i_ind=1,nFock_indep
  ! !      full2indep_fock(fock_indep(i_ind))=i_ind       
  ! !      full2indep_fock(tmp_target(fock_indep(i_ind)))=i_ind
  ! !   end do
  ! !   do i_ind=1,nFock_indep       
  ! !      indep2full_fock(i_ind,1) = fock_indep(i_ind)
  ! !      indep2full_fock(i_ind,2) = tmp_target(fock_indep(i_ind))
  ! !   end do
  ! !   !+- check maps +-!
  ! !   do i_ind=1,nFock_indep
  ! !      do isymm=1,2
  ! !         check_maps=indep2full_fock(i_ind,isymm)
  ! !         if(i_ind /= full2indep_fock(check_maps)) stop "WRONG MAPS"
  ! !      end do
  ! !   end do
  ! ! end subroutine get_spin_indep_states

  ! ! Rotationally invariant Hubbard interaction !
  ! function rotationally_invariant_density_density(cc,ca) result(Oi)
  !   real(8),dimension(Ns,nFock,nFock) :: cc,ca
  !   real(8),dimension(nFock,nFock) :: Oi,Id
  !   real(8),dimension(nFock,nFock) :: ni
  !   integer                        :: i,ispin,iorb,istate
  !   !
  !   Id=0.d0                 
  !   do i=1,nFock
  !      Id(i,i)=1.d0*Norb
  !   end do
  !   !
  !   ni=0.d0    
  !   do ispin=1,2
  !      do iorb=1,Norb
  !         istate=index(ispin,iorb)
  !         ni=ni+matmul(cc(istate,:,:),ca(istate,:,:))
  !      end do
  !   end do
  !   !
  !   Oi=matmul(ni-Id,ni-Id)
  ! end function Rotationally_invariant_density_density



  ! function density_density_interaction(cc,ca) result(Oi)
  !   real(8),dimension(Ns,nFock,nFock) :: cc,ca
  !   real(8),dimension(nFock,nFock) :: Oi,Id
  !   real(8),dimension(nFock,nFock) :: ni,nj,n_up,n_dn
  !   integer                        :: i,ispin,iorb,istate,jorb
  !   !
  !   Id=0.d0                 
  !   do i=1,nFock
  !      Id(i,i)=1.d0
  !   end do
  !   !
  !   ni=0.d0    

  !   Oi=0.d0
  !   do iorb=1,Norb
  !      ispin=1
  !      istate=index(ispin,iorb)
  !      n_up=matmul(cc(istate,:,:),ca(istate,:,:))
  !      ispin=2
  !      istate=index(ispin,iorb)
  !      n_dn=matmul(cc(istate,:,:),ca(istate,:,:))
  !      Oi = Oi + matmul(n_up,n_dn)
  !   end do

  !   do iorb=1,Norb
  !      ni=0.d0
  !      do ispin=1,2
  !         istate=index(ispin,iorb)
  !         ni=ni+matmul(cc(istate,:,:),ca(istate,:,:))
  !      end do
  !      do jorb=1,iorb-1
  !         nj=0.d0
  !         do ispin=1,2
  !            istate=index(ispin,jorb)
  !            nj=nj+matmul(cc(istate,:,:),ca(istate,:,:))
  !         end do
  !         Oi = Oi + matmul(ni,nj)
  !      end do
  !   end do
  ! end function density_density_interaction
  !
END MODULE GZ_TWO_SITES_FOCK
