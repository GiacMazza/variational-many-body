MODULE GZ_MATRIX_BASIS
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_AUX_FUNX
  USE SF_LINALG
  implicit none
  private
  !
  public :: init_variational_matrices
  public :: trace_phi_basis

  public :: enforce_su2_rotational_symmetry
  !
CONTAINS

  include 'symmetries.f90'
  
  subroutine init_variational_matrices
    !
    Nphi = get_dimension_phi_basis()
    !
    call build_matrix_basis
    !
    call build_traces_matrix_basis
    !
  end subroutine init_variational_matrices


  function get_dimension_phi_basis() result(dim_phi)
    integer :: dim_phi
    dim_phi = nFock_indep
  end function get_dimension_phi_basis


  !+- PORCATA: TO BE UPDATED....
  subroutine build_matrix_basis 
    integer :: iphi,ifock,jfock,isymm
    real(8),dimension(nFock,nFock) :: Id
    real(8),dimension(:,:),allocatable :: test_trace
    !
    Id=0.d0; 
    forall(ifock=1:nFock) Id(ifock,ifock)=1.d0
    !
    allocate(phi_basis(Nphi,nFock,nFock),phi_basis_dag(Nphi,nFock,nFock))    
    phi_basis=0.d0
    !
    do iphi=1,Nphi       
       do isymm=1,2
          ifock = indep2full_fock(iphi,isymm)
          phi_basis(iphi,ifock,ifock) = &
               phi_basis(iphi,ifock,ifock) + 1.d0
       end do
    end do
    !
    do iphi=1,Nphi
       do ifock=1,nFock
          do jfock=1,nFock
             phi_basis_dag(iphi,jfock,ifock)=phi_basis(iphi,ifock,jfock)
          end do
       end do
    end do
    !
    allocate(test_trace(Nphi,Nphi))
    test_trace=get_traces_basis_phiOphi(Id)
    !
    do iphi=1,Nphi
       phi_basis(iphi,:,:)     = phi_basis(iphi,:,:)/sqrt(test_trace(iphi,iphi))
       phi_basis_dag(iphi,:,:) = phi_basis_dag(iphi,:,:)/sqrt(test_trace(iphi,iphi))
    end do
    !<SAFE_CHECK
    do iphi=1,Nphi
       write(30,*) iphi,indep2full_fock(iphi,:)
       do ifock=1,Nfock
          write(30,'(20F6.3)') phi_basis(iphi,ifock,:)
       end do
    end do
    !SAFE_CHECK>
    ! test_trace=get_traces_basis_phiOphi(Id)
    ! do iphi=1,Nphi
    !    write(*,'(20F6.3)') test_trace(iphi,:)
    ! end do
    ! write(*,*) "HERE",Nphi
  end subroutine build_matrix_basis


  function trace_phi_basis(phi_vect,phi_trace) result(trace)
    real(8),dimension(Nphi) :: phi_vect
    real(8),dimension(Nphi,Nphi) :: phi_trace
    real(8) :: trace
    integer :: iphi,jphi
    trace=0.d0
    do iphi=1,Nphi
       do jphi=1,Nphi
          trace = trace + phi_vect(iphi)*phi_vect(jphi)*phi_trace(iphi,jphi)
       end do
    end do
  end function trace_phi_basis




  subroutine build_traces_matrix_basis
    integer :: is,js,iorb,jorb
    !+- let's assume that the opertaors in Fock space have been already built


    !+- local operators Tr(Phi+ Oi Phi) -+!
    allocate(phi_traces_basis_Hloc(Nphi,Nphi))
    allocate(phi_traces_basis_free_Hloc(Nphi,Nphi))
    allocate(phi_traces_basis_local_dens(Ns,Ns,Nphi,Nphi))
    allocate(phi_traces_basis_dens_dens_orb(Norb,Norb,Nphi,Nphi))
    allocate(phi_traces_basis_docc_orb(Norb,Nphi,Nphi))    
    allocate(phi_traces_basis_spin_flip(Norb,Norb,Nphi,Nphi))
    allocate(phi_traces_basis_pair_hopping(Norb,Norb,Nphi,Nphi))

    phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
    phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)

    do is=1,Ns
       do js=1,Ns
          phi_traces_basis_local_dens(is,js,:,:) = &
               get_traces_basis_phiOphi(local_dens(is,js,:,:))
       end do
    end do
    !
    do iorb=1,Norb
       do jorb=1,Norb
          phi_traces_basis_dens_dens_orb(iorb,jorb,:,:) = &
               get_traces_basis_phiOphi(dens_dens_orb(iorb,jorb,:,:))
          phi_traces_basis_spin_flip(iorb,jorb,:,:) = &
               get_traces_basis_phiOphi(spin_flip(iorb,jorb,:,:))
          phi_traces_basis_pair_hopping(iorb,jorb,:,:) = &
               get_traces_basis_phiOphi(pair_hopping(iorb,jorb,:,:))
       end do
       phi_traces_basis_docc_orb(iorb,:,:) = &
            get_traces_basis_phiOphi(docc(iorb,:,:))
    end do
    !

    !+- density constraints Tr(Phi+ Phi C_i) C_i=density matrix -+!
    allocate(phi_traces_basis_dens(Ns,Ns,Nphi,Nphi))
    allocate(phi_traces_basis_Cdens(Ns,Ns,Nphi,Nphi))

    do is=1,Ns
       do js=1,Ns
          phi_traces_basis_dens(is,js,:,:) = &
               get_traces_basis_phiphiO_s(local_dens(is,js,:,:))
          phi_traces_basis_Cdens(is,js,:,:) = &
               get_traces_basis_phiphiO(local_dens(is,js,:,:))
       end do
    end do

    !+- Hoppings -+!
    allocate(phi_traces_basis_Rhop(Ns,Ns,Nphi,Nphi))
    do is=1,Ns
       do js=1,Ns
          !phi_traces_basis_Rhop(is,js,:,:) = get_traces_basis_phiAphiB(CA(is,:,:),CC(js,:,:)) 
          phi_traces_basis_Rhop(is,js,:,:) = get_traces_basis_phiAphiB_s(CA(is,:,:),CC(js,:,:)) 
       end do
    end do

  end subroutine build_traces_matrix_basis






  !                                                         !
  !HERE  BUILD ALL THE TRACES MATRICES FOR LOCAL OBSERVABLES! 
  !                                                         !

  function get_traces_basis_phiOphi(Oi) result(trace_matrix)
    real(8),dimension(nFock,nFock) :: Oi !local_operator whose traces should be computed
    real(8),dimension(Nphi,Nphi)   :: trace_matrix
    real(8),dimension(nFock,nFock) :: tmp
    integer                        :: kfock,iphi,jphi
    !
    trace_matrix=0.d0                
    do iphi=1,Nphi
       do jphi=1,Nphi
          tmp=0.d0
          tmp=matmul(Oi,phi_basis(jphi,:,:))
          tmp=matmul(phi_basis_dag(iphi,:,:),tmp)
          do kfock=1,nFock
             trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + tmp(kfock,kfock)                  
          end do
       end do
    end do
    !
  end function get_traces_basis_phiOphi
  !
  function get_traces_basis_phiphiO(Oi) result(trace_matrix)
    real(8),dimension(nFock,nFock) :: Oi !local_operator whose traces should be computed
    real(8),dimension(Nphi,Nphi)   :: trace_matrix
    real(8),dimension(nFock,nFock) :: tmp
    integer                        :: kfock,iphi,jphi
    !
    trace_matrix=0.d0                
    do iphi=1,Nphi
       do jphi=1,Nphi
          tmp=0.d0
          tmp=matmul(phi_basis(jphi,:,:),Oi)
          tmp=matmul(phi_basis_dag(iphi,:,:),tmp)
          do kfock=1,nFock
             trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + tmp(kfock,kfock)                  
          end do
       end do
    end do
    !
  end function get_traces_basis_phiphiO

  function get_traces_basis_phiphiO_s(Oi) result(trace_matrix)
    real(8),dimension(nFock,nFock) :: Oi !local_operator whose traces should be computed
    real(8),dimension(Nphi,Nphi)   :: trace_matrix
    real(8),dimension(nFock,nFock) :: tmp
    integer                        :: kfock,iphi,jphi
    !
    trace_matrix=0.d0                
    do iphi=1,Nphi
       do jphi=1,Nphi
          tmp=0.d0
          tmp=matmul(phi_basis(jphi,:,:),Oi)
          tmp=matmul(phi_basis_dag(iphi,:,:),tmp)
          do kfock=1,nFock
             trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + 0.5d0*tmp(kfock,kfock)                  
          end do
          tmp=0.d0
          tmp=matmul(phi_basis(iphi,:,:),Oi)
          tmp=matmul(phi_basis_dag(jphi,:,:),tmp)
          do kfock=1,nFock
             trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + 0.5d0*tmp(kfock,kfock)                  
          end do
       end do
    end do
    !
  end function get_traces_basis_phiphiO_s


  function get_traces_basis_phiAphiB(A,B) result(trace_matrix)
    real(8),dimension(nFock,nFock) :: A,B 
    real(8),dimension(Nphi,Nphi)   :: trace_matrix
    real(8),dimension(nFock,nFock) :: tmp
    integer                        :: kfock,iphi,jphi
    !
    trace_matrix=0.d0                
    do iphi=1,Nphi
       do jphi=1,Nphi
          tmp=0.d0
          tmp=matmul(phi_basis(jphi,:,:),B)
          tmp=matmul(A,tmp)
          tmp=matmul(phi_basis_dag(iphi,:,:),tmp)
          do kfock=1,nFock
             trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + tmp(kfock,kfock)                  
          end do
       end do
    end do
    !
  end function get_traces_basis_phiAphiB



  function get_traces_basis_phiAphiB_s(A,B) result(trace_matrix)
    real(8),dimension(nFock,nFock) :: A,B 
    real(8),dimension(Nphi,Nphi)   :: trace_matrix
    real(8),dimension(nFock,nFock) :: tmp
    integer                        :: kfock,iphi,jphi
    !
    trace_matrix=0.d0                
    do iphi=1,Nphi
       do jphi=1,Nphi
          tmp=0.d0
          tmp=matmul(phi_basis(jphi,:,:),B)
          tmp=matmul(A,tmp)
          tmp=matmul(phi_basis_dag(iphi,:,:),tmp)
          do kfock=1,nFock
             trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + 0.5d0*tmp(kfock,kfock)
          end do
          tmp=0.d0
          tmp=matmul(phi_basis(iphi,:,:),B)
          tmp=matmul(A,tmp)
          tmp=matmul(phi_basis_dag(jphi,:,:),tmp)
          do kfock=1,nFock
             trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + 0.5d0*tmp(kfock,kfock)
          end do
       end do
    end do
    !
  end function get_traces_basis_phiAphiB_s
  
  !
END MODULE GZ_MATRIX_BASIS
