MODULE GZ_LOCAL_HAMILTONIAN
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_AUX_FUNX
  USE SF_LINALG
  USE SF_IOTOOLS
  USE MATRIX_SPARSE
  implicit none
  private


  public :: get_local_hamiltonian_trace
  public :: get_local_hamiltonian_HLOCphi


CONTAINS
  !
  subroutine get_local_hamiltonian_trace(atomic_energy_levels_)    
    real(8),dimension(Ns),optional :: atomic_energy_levels_
    real(8),dimension(Ns) :: atomic_energy_levels
    integer :: iorb,jorb,ispin,jspin,istate,jstate,is_up,is_dn,js_up,js_dn,ifock,is
    ! real(8),dimension(Ns,nFock,nFock) :: state_dens
    ! real(8),dimension(nFock,nFock) :: tmp
    real(8) :: mu_ph
    !
    if(allocated(phi_traces_basis_Hloc)) deallocate(phi_traces_basis_Hloc)
    if(allocated(phi_traces_basis_free_Hloc)) deallocate(phi_traces_basis_free_Hloc)
    allocate(phi_traces_basis_Hloc(Nphi,Nphi));phi_traces_basis_Hloc = zero
    allocate(phi_traces_basis_free_Hloc(Nphi,Nphi));phi_traces_basis_free_Hloc = zero
    !
    !
    atomic_energy_levels=0.d0;if(present(atomic_energy_levels_)) atomic_energy_levels=atomic_energy_levels_
    !
    !
    write(*,*) 'local_hamiltonian',Uloc,Ust,Jh,Jsf,Jph,atomic_energy_levels
    !+- FREE PART OF THE LOCAL HAMILTONIAN -+!
    do istate=1,Ns
       phi_traces_basis_free_Hloc = phi_traces_basis_free_Hloc + (atomic_energy_levels(istate)-xmu)*phi_traces_basis_local_dens(istate,istate,:,:)
       ! local_hamiltonian = local_hamiltonian + atomic_energy_levels(istate)*state_dens(istate,:,:)
       ! local_hamiltonian = local_hamiltonian - xmu*state_dens(istate,:,:)
    end do
    phi_traces_basis_Hloc = phi_traces_basis_free_Hloc   
    !+- LOCAL INTERACTION HAMILTONIAN -+!

    !+- INTRA-ORBITAL density-density  -+!
    do iorb=1,Norb
       is_up = index(1,iorb)
       is_dn = index(2,iorb)
       phi_traces_basis_Hloc =  phi_traces_basis_Hloc + &
            Uloc(iorb)*phi_traces_basis_dens_dens(is_up,is_dn,:,:)       
    end do

    if(Norb.gt.1) then
       !+- INTER-ORBITAL density-density opposite spins -+!
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             !
             is_up = index(1,iorb)
             is_dn = index(2,iorb)
             !
             js_up = index(1,jorb)
             js_dn = index(2,jorb)
             !
             phi_traces_basis_Hloc =  phi_traces_basis_Hloc + &
                  Ust*phi_traces_basis_dens_dens(is_up,js_dn,:,:) + &
                  Ust*phi_traces_basis_dens_dens(is_dn,js_up,:,:)
          end do
       end do
       !+- INTER-ORBITAL density-density parallel spins -+!
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             !
             is_up = index(1,iorb)
             is_dn = index(2,iorb)
             !
             js_up = index(1,jorb)
             js_dn = index(2,jorb)
             !             
             phi_traces_basis_Hloc =  phi_traces_basis_Hloc + &
                  (Ust-Jh)*phi_traces_basis_dens_dens(is_up,js_up,:,:) + &
                  (Ust-Jh)*phi_traces_basis_dens_dens(is_dn,js_dn,:,:)
          end do
       end do
       !+- SPIN-FLIP interaction -+!
       do iorb=1,Norb
          do jorb=1,Norb
             !
             is_up = index(1,iorb)
             is_dn = index(2,iorb)
             !
             js_up = index(1,jorb)
             js_dn = index(2,jorb)
             !             
             if(iorb.ne.jorb) then
                phi_traces_basis_Hloc =  phi_traces_basis_Hloc - Jsf*phi_traces_basis_spin_flip(iorb,jorb,:,:)                     
             end if
          end do
       end do
       !+- PAIR-HOPPING interaction -+!
       do iorb=1,Norb
          do jorb=1,Norb
             !
             is_up = index(1,iorb)
             is_dn = index(2,iorb)
             !
             js_up = index(1,jorb)
             js_dn = index(2,jorb)
             !             
             if(iorb.ne.jorb) then
                phi_traces_basis_Hloc =  phi_traces_basis_Hloc + Jph*phi_traces_basis_pair_hopping(iorb,jorb,:,:)
             end if
             !
          end do
       end do
    end if
    !+- CHEMICAL POTENTIAL FOR PH CONDITION -+!
    mu_ph = Uloc(1)*0.5d0 + dble(Norb-1)*0.5d0*(2.d0*Ust-Jh) !
    do iorb=1,Norb
       do ispin=1,2
          istate=index(ispin,iorb)          
          phi_traces_basis_Hloc =  phi_traces_basis_Hloc - mu_ph*phi_traces_basis_local_dens(istate,istate,:,:)
       end do
    end do
    !
  end subroutine get_local_hamiltonian_trace




  function get_local_hamiltonian_HLOCphi(gzproj_in,atomic_energy_levels_) result(gzproj_out)
    complex(8),dimension(Nphi) :: gzproj_in,gzproj_out
    type(sparse_matrix_csr_z) :: htmp
    real(8),dimension(Ns),optional :: atomic_energy_levels_
    real(8),dimension(Ns) :: atomic_energy_levels
    integer :: iorb,jorb,ispin,jspin,istate,jstate,is_up,is_dn,js_up,js_dn,ifock,is
    ! real(8),dimension(Ns,nFock,nFock) :: state_dens
    ! real(8),dimension(nFock,nFock) :: tmp
    real(8) :: mu_ph
    !
    !
    atomic_energy_levels=0.d0;if(present(atomic_energy_levels_)) atomic_energy_levels=atomic_energy_levels_
    !
    gzproj_out = zero
    !
    !+- FREE PART OF THE LOCAL HAMILTONIAN -+!
    do istate=1,Ns
       htmp = sp_scalar_matrix_csr(phi_spTraces_basis_local_dens(istate,istate),(atomic_energy_levels(istate)-xmu))
       gzproj_out = gzproj_out + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj_in)
       call sp_delete_matrix(htmp)
       !    phi_traces_basis_free_Hloc = phi_traces_basis_free_Hloc + (atomic_energy_levels(istate)-xmu)*phi_traces_basis_local_dens(istate,istate,:,:)
    end do
    ! phi_traces_basis_Hloc = phi_traces_basis_free_Hloc
    
    !+- LOCAL INTERACTION HAMILTONIAN -+!

    !+- INTRA-ORBITAL density-density  -+!
    do iorb=1,Norb
       is_up = index(1,iorb)
       is_dn = index(2,iorb)
       !
       htmp = sp_scalar_matrix_csr(phi_spTraces_basis_dens_dens(is_up,is_dn),Uloc(iorb))
       gzproj_out = gzproj_out + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj_in)
       call sp_delete_matrix(htmp)
       !
       ! phi_traces_basis_Hloc =  phi_traces_basis_Hloc + &
       !      Uloc(iorb)*phi_traces_basis_dens_dens(is_up,is_dn,:,:)       
    end do
    
    if(Norb.gt.1) then
       !+- INTER-ORBITAL density-density opposite spins -+!
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             !
             is_up = index(1,iorb)
             is_dn = index(2,iorb)
             !
             js_up = index(1,jorb)
             js_dn = index(2,jorb)
             !
             htmp = sp_scalar_matrix_csr(phi_spTraces_basis_dens_dens(is_up,js_dn),Ust)
             gzproj_out = gzproj_out + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj_in)
             call sp_delete_matrix(htmp)
             !
             htmp = sp_scalar_matrix_csr(phi_spTraces_basis_dens_dens(is_dn,js_up),Ust)
             gzproj_out = gzproj_out + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj_in)
             call sp_delete_matrix(htmp)

                 
             ! phi_traces_basis_Hloc =  phi_traces_basis_Hloc + &
             !      Ust*phi_traces_basis_dens_dens(is_up,js_dn,:,:) + &
             !      Ust*phi_traces_basis_dens_dens(is_dn,js_up,:,:)
          end do
       end do
       !+- INTER-ORBITAL density-density parallel spins -+!
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             !
             is_up = index(1,iorb)
             is_dn = index(2,iorb)
             !
             js_up = index(1,jorb)
             js_dn = index(2,jorb)
             !
             htmp = sp_scalar_matrix_csr(phi_spTraces_basis_dens_dens(is_up,js_up),Ust-Jh)
             gzproj_out = gzproj_out + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj_in)
             call sp_delete_matrix(htmp)

             htmp = sp_scalar_matrix_csr(phi_spTraces_basis_dens_dens(is_dn,js_dn),Ust-Jh)
             gzproj_out = gzproj_out + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj_in)
             call sp_delete_matrix(htmp)

             
             ! phi_traces_basis_Hloc =  phi_traces_basis_Hloc + &
             !      (Ust-Jh)*phi_traces_basis_dens_dens(is_up,js_up,:,:) + &
             !      (Ust-Jh)*phi_traces_basis_dens_dens(is_dn,js_dn,:,:)
             
          end do
       end do
       !+- SPIN-FLIP interaction -+!
       do iorb=1,Norb
          do jorb=1,Norb
             !
             is_up = index(1,iorb)
             is_dn = index(2,iorb)
             !
             js_up = index(1,jorb)
             js_dn = index(2,jorb)
             !             
             if(iorb.ne.jorb) then
                htmp = sp_scalar_matrix_csr(phi_spTraces_basis_spin_flip(iorb,jorb),-1.d0*Jsf)
                gzproj_out = gzproj_out + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj_in)
                call sp_delete_matrix(htmp)
   
                
                ! phi_traces_basis_Hloc =  phi_traces_basis_Hloc - Jsf*phi_traces_basis_spin_flip(iorb,jorb,:,:)
             end if
          end do
       end do
       !+- PAIR-HOPPING interaction -+!
       do iorb=1,Norb
          do jorb=1,Norb
             !
             is_up = index(1,iorb)
             is_dn = index(2,iorb)
             !
             js_up = index(1,jorb)
             js_dn = index(2,jorb)
             !             
             if(iorb.ne.jorb) then
                htmp = sp_scalar_matrix_csr(phi_spTraces_basis_pair_hopping(iorb,jorb),Jph)
                gzproj_out = gzproj_out + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj_in)
                call sp_delete_matrix(htmp)

                
                ! phi_traces_basis_Hloc =  phi_traces_basis_Hloc + Jph*phi_traces_basis_pair_hopping(iorb,jorb,:,:)
             end if
             !
          end do
       end do
    end if
    !+- CHEMICAL POTENTIAL FOR PH CONDITION -+!

    mu_ph = Uloc(1)*0.5d0 + dble(Norb-1)*0.5d0*(2.d0*Ust-Jh) !
    do iorb=1,Norb
       do ispin=1,2
          istate=index(ispin,iorb)
          !
          htmp = sp_scalar_matrix_csr(phi_spTraces_basis_local_dens(istate,istate),-1.d0*mu_ph)
          gzproj_out = gzproj_out + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj_in)
          call sp_delete_matrix(htmp)
          !
          ! phi_traces_basis_Hloc =  phi_traces_basis_Hloc - mu_ph*phi_traces_basis_local_dens(istate,istate,:,:)
          !
       end do
    end do
    !
  end function get_local_hamiltonian_HLOCphi




  



END MODULE GZ_LOCAL_HAMILTONIAN
