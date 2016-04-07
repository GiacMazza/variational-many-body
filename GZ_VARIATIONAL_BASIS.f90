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
  !
  public :: get_traces_basis_phiOphi
  !
  public :: symmetry_stride_vec_mat,symmetry_stride_mat_vec
  public :: symmetry_stride_vec_mat_,symmetry_stride_mat_vec_
  !
  type intarray
     private
     integer,dimension(:),allocatable :: index
  end type intarray
  type local_multiplets
     private
     integer :: N_mult
     integer :: size_mult
     real(8),dimension(:,:),allocatable :: mult
     integer,dimension(:),allocatable :: Nequiv_mult
     type(intarray),dimension(:),allocatable :: Maps
     integer,dimension(:),allocatable :: inv_map
  end type local_multiplets
  !
CONTAINS
  !
  include 'symmetries.f90'
  !
  subroutine init_variational_matrices
    !
    Nphi = get_dimension_phi_basis()
    call build_traces_matrix_basis
    !
  end subroutine init_variational_matrices



  !< VERY VERY TEMPORARAY ROUTINES
  function symmetry_stride_vec_mat(lgr) result(lgr_matrix)
    real(8),dimension(:) :: lgr
    complex(8),dimension(:),allocatable ::  lgr_cmplx
    integer :: dim
    complex(8),dimension(Ns,Ns) :: lgr_matrix

    integer :: i,ilgr,iorb,jorb,ispin,istate,jstate

    dim = Norb*(Norb+1)/2
    if(size(lgr).ne.2*dim) stop "what the fuck vec_mat!"
    allocate(lgr_cmplx(dim))
    do i=1,dim
       lgr_cmplx(i) = lgr(i) + xi*lgr(i+dim)
    end do
    !
    lgr_matrix=zero
    ilgr=0
    do iorb=1,Norb
       do jorb=1,iorb
          ilgr = ilgr + 1
          do ispin=1,2
             !
             istate=index(ispin,iorb)
             jstate=index(ispin,jorb)
             !
             lgr_matrix(istate,jstate) = lgr_cmplx(ilgr) 
             if(iorb.ne.jorb) lgr_matrix(jstate,istate) = conjg(lgr_matrix(istate,jstate))
             !
          end do
       end do
    end do
    !
    ! dim = Norb*(Norb+1)/2
    ! if(size(lgr).ne.2*dim) stop "what the fuck vec_mat!"
    ! allocate(lgr_cmplx(dim))
    ! do i=1,dim
    !    lgr_cmplx(i) = lgr(i) + xi*lgr(i+dim)
    ! end do
    !
  end function symmetry_stride_vec_mat




  subroutine symmetry_stride_mat_vec(mat,vec) 
    complex(8),dimension(Ns,Ns) :: mat
    real(8),dimension(:) :: vec
    integer :: dim
    integer ::i,ilgr,iorb,jorb,ispin,istate,jstate
    !
    dim = Norb*(Norb+1)/2
    if(size(vec).ne.2*dim) stop "what the fuck mat_vec!"
    !
    ilgr=0
    do iorb=1,Norb
       do jorb=1,iorb
          ilgr = ilgr + 1
          do ispin=1,2
             !
             istate=index(ispin,iorb)
             jstate=index(ispin,jorb)
             !
             vec(ilgr) = dreal(mat(istate,jstate))
             vec(ilgr+dim) = dimag(mat(istate,jstate))
             !
          end do
       end do
    end do
    !
  end subroutine symmetry_stride_mat_vec



  !+- anomalous spin-singlet channel -+! 
  function symmetry_stride_vec_mat_(lgr) result(lgr_matrix)
    real(8),dimension(:) :: lgr
    complex(8),dimension(:),allocatable ::  lgr_cmplx
    integer :: dim
    complex(8),dimension(Ns,Ns) :: lgr_matrix

    integer ::i,ilgr,iorb,jorb,ispin,jspin,istate,jstate

    dim = Norb*(Norb+1)/2
    if(size(lgr).ne.2*dim) stop "what the fuck vec_mat_!"
    allocate(lgr_cmplx(dim))
    do i=1,dim
       lgr_cmplx(i) = lgr(i) + xi*lgr(i+dim)
    end do
    !
    ilgr=0
    lgr_matrix=zero
    do iorb=1,Norb
       do jorb=1,iorb
          ilgr = ilgr + 1
          !
          do ispin=1,2
             do jspin=ispin+1,2
                istate=index(ispin,iorb)
                jstate=index(jspin,jorb)
                write(*,*) istate,jstate
                lgr_matrix(istate,jstate) =  lgr_cmplx(ilgr) 
                lgr_matrix(jstate,istate) = -lgr_cmplx(ilgr) 
             end do
          end do
       end do
    end do
    !
  end function symmetry_stride_vec_mat_




  subroutine symmetry_stride_mat_vec_(mat,vec) 
    complex(8),dimension(Ns,Ns) :: mat
    real(8),dimension(:) :: vec
    integer :: dim
    integer ::i,ilgr,iorb,jorb,ispin,jspin,istate,jstate
    !
    dim = Norb*(Norb+1)/2
    if(size(vec).ne.2*dim) stop "what the fuck mat_vec_!"

    ilgr=0
    do iorb=1,Norb
       do jorb=1,iorb
          ilgr = ilgr + 1
          !
          do ispin=1,2
             do jspin=ispin+1,2
                istate=index(ispin,iorb)
                jstate=index(jspin,jorb)

                vec(ilgr) = dreal(mat(istate,jstate)) 
                vec(ilgr+dim) = dimag(mat(istate,jstate))
             end do
          end do
          !          
       end do
    end do
    !
  end subroutine symmetry_stride_mat_vec_
  !
  !< VERY VERY TEMPORARY ROUTINES
  !
  function get_dimension_phi_basis() result(dim_phi)
    integer                                         :: dim_phi
    integer,dimension(:,:),allocatable              :: irr_reps
    integer,dimension(:,:),allocatable              :: equ_reps
    complex(8),dimension(nFock,nFock)               :: Virr_reps
    integer                                         :: Nirr_reps,Nineq_reps,i,j,Neq,ifock,jfock,iphi,jphi

    complex(8),dimension(:,:,:),allocatable         :: phi_irr
    complex(8),dimension(:,:,:),allocatable         :: phi_fock

    complex(8),dimension(nFock,nFock)               :: tmp_matrix
    complex(8)                                      :: tmp_trace

    real(8),dimension(nFock,nFock)                  :: Id
    real(8),dimension(:,:),allocatable              :: test_trace
    !
    type(local_multiplets),dimension(:),allocatable :: mult_list    
    !
    real(8)                                         :: tmp_check

    Id=0.d0; 
    forall(ifock=1:nFock) Id(ifock,ifock)=1.d0

    ! call basis_SZ_irr_reps_test(mult_list,Virr_reps)
    ! call get_matrix_basis_irr_reps_test(mult_list,phi_irr)
    ! stop

    select case(wf_symmetry)
    case(0)
       call basis_O1xSU2_irr_reps(irr_reps,equ_reps,Virr_reps)
    case(1)
       call basis_O1cXSU2sXSU2c_irr_reps(irr_reps,equ_reps,Virr_reps)
    case(2)
       call basis_O1cXSU2sXisoZ_irr_reps(irr_reps,equ_reps,Virr_reps)
    case(3)
       call basis_SU2sXSU2c_irr_reps(irr_reps,equ_reps,Virr_reps)
    case(4)
       call basis_SU2_irr_reps(irr_reps,equ_reps,Virr_reps)
    case(5)
       call basis_SZ_irr_reps(irr_reps,equ_reps,Virr_reps)
    case(6)
       call basis_SU2sXisoZ_irr_reps(irr_reps,equ_reps,Virr_reps)
    end select
    !
    Nirr_reps=size(irr_reps,1)
    Nineq_reps=size(equ_reps,2)
    !
    write(*,*) Nirr_reps,Nineq_reps
    !
    write(*,*) 'IRREDUCIBLE REPS'
    do i=1,Nirr_reps
       write(*,*) irr_reps(i,:)
    end do
    !
    write(*,*) 'EQUIVALENT IRREDUCIBLE REPS'
    do i=1,Nirr_reps
       write(*,*) equ_reps(i,:)
    end do
    !
    call get_matrix_basis_irr_reps(irr_reps,equ_reps,phi_irr)
    call get_matrix_basis_original_fock(phi_irr,phi_fock,Virr_reps)
    !



    dim_phi=size(phi_fock,1)
    !
    !dim_phi = nFock*nFock
    Nphi=dim_phi   
    write(*,*) "NPHI",Nphi
    !
    allocate(phi_basis(dim_phi,nFock,nFock))
    allocate(phi_basis_dag(dim_phi,nFock,nFock))

    phi_basis = phi_fock


    !< tmp TEST
    ! phi_basis=0.d0
    ! iphi=0
    ! do ifock=1,nFock
    !    do jfock=1,nFock
    !       iphi=iphi+1
    !       phi_basis(iphi,ifock,jfock) = 1.d0
    !    end do
    ! end do
    ! END TMP_TEST>

    do iphi=1,Nphi
       do ifock=1,nFock
          do jfock=1,nFock
             phi_basis_dag(iphi,jfock,ifock)=conjg(phi_basis(iphi,ifock,jfock))
          end do
       end do
    end do
    ! 

    allocate(test_trace(Nphi,Nphi))
    test_trace=get_traces_basis_phiOphi(Id)

    !+- safe check on the matrix orthogonality
    do iphi=1,Nphi
       do jphi=1,Nphi
          if(abs(test_trace(iphi,jphi)).gt.1.d-8.and.iphi.ne.jphi) &
               stop "Variational basis not Trace-orthogonal"
       end do
    end do

    do iphi=1,Nphi
       phi_basis(iphi,:,:)     = phi_basis(iphi,:,:)/sqrt(test_trace(iphi,iphi))
       phi_basis_dag(iphi,:,:) = phi_basis_dag(iphi,:,:)/sqrt(test_trace(iphi,iphi))
    end do

    !<TMP_CHECK
    ! tmp_check=0.d0
    ! do iphi=1,Nphi
    !    do ifock=1,nFock
    !       do jfock=1,nFock
    !          tmp_check = tmp_check + dimag(phi_basis(iphi,ifock,jfock))**2.d0
    !       end do
    !    end do
    ! end do
    ! write(*,*) 'TMP_CHECK',tmp_check
    ! stop
    !TMP_CHECK>

    !stop "Variational basis  TRACE-ORTHOGONAL"

    !
    ! do iphi=1,Nphi
    !    write(30,*) iphi!,indep2full_fock(iphi,:)
    !    do ifock=1,Nfock
    !       write(30,'(20F6.3)') phi_basis(iphi,ifock,:)
    !    end do
    ! end do
    ! !SAFE_CHECK>
    ! test_trace=get_traces_basis_phiOphi(Id)
    ! do iphi=1,Nphi
    !    write(*,'(20F6.3)') test_trace(iphi,:)
    ! end do
    !
    !write(*,*) "NPHI",Nphi




    !stop

    !dim_phi = nFock_indep
    !dim_phi=nFock
    !< check full phi_matrix
    !dim_phi = nFock*nfock
    !>

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
    ! do iphi=1,Nphi       
    !    do isymm=1,2
    !       ifock = indep2full_fock(iphi,isymm)
    !       phi_basis(iphi,ifock,ifock) = &
    !            phi_basis(iphi,ifock,ifock) + 1.d0
    !    end do
    ! end do
    !

    !
    ! do ifock=1,nFock
    !    phi_basis(ifock,ifock,ifock) = 1.d0
    ! end do
    !

    phi_basis = 0.d0
    do ifock=1,nFock
       do jfock=1,nFock
          iphi = (ifock-1)*nFock + jfock
          phi_basis(iphi,ifock,jfock) = 1.d0
       end do
    end do
    ! !
    do iphi=1,Nphi
       do ifock=1,nFock
          do jfock=1,nFock
             phi_basis_dag(iphi,jfock,ifock)=conjg(phi_basis(iphi,ifock,jfock))
          end do
       end do
    end do
    ! !
    allocate(test_trace(Nphi,Nphi))
    test_trace=get_traces_basis_phiOphi(Id)
    !
    do iphi=1,Nphi
       phi_basis(iphi,:,:)     = phi_basis(iphi,:,:)/sqrt(test_trace(iphi,iphi))
       phi_basis_dag(iphi,:,:) = phi_basis_dag(iphi,:,:)/sqrt(test_trace(iphi,iphi))
    end do
    !<SAFE_CHECK
    do iphi=1,Nphi
       write(30,*) iphi!,indep2full_fock(iphi,:)
       do ifock=1,Nfock
          write(30,'(20F6.3)') phi_basis(iphi,ifock,:)
       end do
    end do
    !SAFE_CHECK>
    test_trace=get_traces_basis_phiOphi(Id)
    do iphi=1,Nphi
       write(*,'(20F6.3)') test_trace(iphi,:)
    end do
    write(*,*) "HERE",Nphi

    !    stop
  end subroutine build_matrix_basis


  function trace_phi_basis(phi_vect,phi_trace) result(trace)
    complex(8),dimension(Nphi) :: phi_vect
    complex(8),dimension(Nphi,Nphi) :: phi_trace
    complex(8) :: trace
    integer :: iphi,jphi
    trace=0.d0
    do iphi=1,Nphi
       do jphi=1,Nphi
          trace = trace + conjg(phi_vect(iphi))*phi_vect(jphi)*phi_trace(iphi,jphi)
       end do
    end do
  end function trace_phi_basis




  subroutine build_traces_matrix_basis
    integer :: is,js,iorb,jorb,iphi,jphi
    !+- let's assume that the opertaors in Fock space have been already built



    !+- local operators Tr(Phi+ Oi Phi) -+!
    allocate(phi_traces_basis_Hloc(Nphi,Nphi))
    allocate(phi_traces_basis_free_Hloc(Nphi,Nphi))
    allocate(phi_traces_basis_local_dens(Ns,Ns,Nphi,Nphi))
    allocate(phi_traces_basis_dens_dens_orb(Norb,Norb,Nphi,Nphi))
    allocate(phi_traces_basis_docc_orb(Norb,Nphi,Nphi))    
    allocate(phi_traces_basis_spin_flip(Norb,Norb,Nphi,Nphi))
    allocate(phi_traces_basis_pair_hopping(Norb,Norb,Nphi,Nphi))
    allocate(phi_traces_basis_sc_order(Ns,Ns,Nphi,Nphi))
    !
    allocate(phi_traces_basis_spin2(Nphi,Nphi))
    allocate(phi_traces_basis_spinZ(Nphi,Nphi))
    allocate(phi_traces_basis_isospin2(Nphi,Nphi))
    allocate(phi_traces_basis_isospinZ(Nphi,Nphi))

    phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
    phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)
    do is=1,Ns
       do js=1,Ns
          phi_traces_basis_local_dens(is,js,:,:) = &
               get_traces_basis_phiOphi(op_local_dens(is,js,:,:))
          write(*,*) is,js
       end do
    end do

    !
    do iorb=1,Norb
       do jorb=1,Norb
          phi_traces_basis_dens_dens_orb(iorb,jorb,:,:) = &
               get_traces_basis_phiOphi(op_dens_dens_orb(iorb,jorb,:,:))
          !
          phi_traces_basis_spin_flip(iorb,jorb,:,:) = &
               get_traces_basis_phiOphi(op_spin_flip(iorb,jorb,:,:))
          !
          phi_traces_basis_pair_hopping(iorb,jorb,:,:) = &
               get_traces_basis_phiOphi(op_pair_hopping(iorb,jorb,:,:))
       end do
       phi_traces_basis_docc_orb(iorb,:,:) = &
            get_traces_basis_phiOphi(op_docc(iorb,:,:))

    end do
    !

    !+- density constraints Tr(Phi+ Phi C_i) C_i=density matrix -+!
    allocate(phi_traces_basis_dens(Ns,Ns,Nphi,Nphi),phi_traces_basis_dens_hc(Ns,Ns,Nphi,Nphi))
    do is=1,Ns
       do js=1,Ns
          phi_traces_basis_dens(is,js,:,:) = &    !+- probably is more correct to call this global variable !+- phi_traces_basis_vdm
               get_traces_basis_phiphiO(op_local_dens(is,js,:,:))

          ! phi_traces_basis_dens(is,js,:,:) = &    !+- probably is more correct to call this global variable !+- phi_traces_basis_vdm
          !      get_traces_basis_phiphiO_s(op_local_dens(is,js,:,:))

       end do
    end do
    do iphi=1,Nphi
       do jphi=1,Nphi
          phi_traces_basis_dens_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_dens(:,:,jphi,iphi))
       end do
    end do
    !+- anomalous part -+!
    allocate(phi_traces_basis_dens_anomalous(Ns,Ns,Nphi,Nphi))
    allocate(phi_traces_basis_dens_anomalous_hc(Ns,Ns,Nphi,Nphi))
    do is=1,Ns
       do js=1,Ns
          phi_traces_basis_dens_anomalous(is,js,:,:) = &    !+- probably is more correct to call this global variable !+- phi_traces_basis_vdm
               get_traces_basis_phiphiO(op_local_dens_anomalous(is,js,:,:))
       end do
    end do
    do iphi=1,Nphi
       do jphi=1,Nphi
          phi_traces_basis_dens_anomalous_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_dens_anomalous(:,:,jphi,iphi))
       end do
    end do





    !<DEBUG
    do iphi=1,Nphi
       do jphi=1,Nphi
          write(907,'(10F18.10)') dble(iphi),dble(jphi),phi_traces_basis_dens(1,1,iphi,jphi),phi_traces_basis_dens(2,2,iphi,jphi),phi_traces_basis_dens(3,3,iphi,jphi),phi_traces_basis_dens(4,4,iphi,jphi)
       end do
    end do
    do iphi=1,nFock
       do jphi=1,nFock
          write(908,'(10F18.10)') op_local_dens(1,1,iphi,jphi),op_local_dens(2,2,iphi,jphi),op_local_dens(3,3,iphi,jphi),op_local_dens(4,4,iphi,jphi)
       end do
    end do
    do iphi=1,Nphi
       do jphi=1,Nphi
          write(909,'(10F18.10)') phi_traces_basis_local_dens(1,1,iphi,jphi),phi_traces_basis_local_dens(2,2,iphi,jphi),phi_traces_basis_local_dens(3,3,iphi,jphi),phi_traces_basis_local_dens(4,4,iphi,jphi)
       end do
    end do
    !    stop
    !DEBUG>



    !+- Hoppings -+!
    allocate(phi_traces_basis_Rhop(Ns,Ns,Nphi,Nphi),phi_traces_basis_Rhop_hc(Ns,Ns,Nphi,Nphi))
    do is=1,Ns
       do js=1,Ns
          !phi_traces_basis_Rhop(is,js,:,:) = get_traces_basis_phiAphiB_s(CA(is,:,:),CC(js,:,:)) 
          phi_traces_basis_Rhop(is,js,:,:) = get_traces_basis_phiAphiB(CA(is,:,:),CC(js,:,:)) 
       end do
    end do
    do iphi=1,Nphi
       do jphi=1,Nphi
          phi_traces_basis_Rhop_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_Rhop(:,:,jphi,iphi))
       end do
    end do
    !
    if(gz_superc) then
       allocate(phi_traces_basis_Qhop(Ns,Ns,Nphi,Nphi),phi_traces_basis_Qhop_hc(Ns,Ns,Nphi,Nphi))
       do is=1,Ns
          do js=1,Ns
             phi_traces_basis_Qhop(is,js,:,:) = get_traces_basis_phiAphiB(CA(is,:,:),CA(js,:,:)) 
          end do
       end do
       do iphi=1,Nphi
          do jphi=1,Nphi
             phi_traces_basis_Qhop_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_Qhop(:,:,jphi,iphi))
          end do
       end do
    end if
    !
    do is=1,Ns
       do js=1,Ns
          phi_traces_basis_sc_order(is,js,:,:) = get_traces_basis_phiOphi(op_sc_order(is,js,:,:)) 
       end do
    end do
    !
    phi_traces_basis_spin2 = get_traces_basis_phiOphi_z(op_spin2) 
    phi_traces_basis_spinZ = get_traces_basis_phiOphi_z(op_spinZ) 
    phi_traces_basis_isospin2 = get_traces_basis_phiOphi_z(op_isospin2) 
    phi_traces_basis_isospinZ = get_traces_basis_phiOphi_z(op_isospinZ) 
    !
  end subroutine build_traces_matrix_basis






  !                                                         !
  !HERE  BUILD ALL THE TRACES MATRICES FOR LOCAL OBSERVABLES! 
  !                                                         !

  function get_traces_basis_phiOphi(Oi) result(trace_matrix)
    real(8),dimension(nFock,nFock) :: Oi !local_operator whose traces should be computed
    complex(8),dimension(Nphi,Nphi)   :: trace_matrix
    complex(8),dimension(nFock,nFock) :: tmp
    integer                        :: kfock,ifock,jfock,iphi,jphi
    !
    trace_matrix=0.d0                
    do iphi=1,Nphi
       do jphi=1,Nphi
          tmp=0.d0
          !
          do ifock=1,nFock
             do jfock=1,nFock
                do kfock=1,nFock
                   trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + &
                        conjg(phi_basis(iphi,ifock,kfock))*phi_basis(jphi,jfock,kfock)*Oi(ifock,jfock)
                end do
             end do
          end do
          !
       end do
    end do
    !
  end function get_traces_basis_phiOphi

  function get_traces_basis_phiOphi_z(Oi) result(trace_matrix)
    complex(8),dimension(nFock,nFock) :: Oi !local_operator whose traces should be computed
    complex(8),dimension(Nphi,Nphi)   :: trace_matrix
    complex(8),dimension(nFock,nFock) :: tmp
    integer                        :: kfock,ifock,jfock,iphi,jphi
    !
    trace_matrix=0.d0                
    do iphi=1,Nphi
       do jphi=1,Nphi
          tmp=0.d0
          !
          do ifock=1,nFock
             do jfock=1,nFock
                do kfock=1,nFock
                   trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + &
                        conjg(phi_basis(iphi,ifock,kfock))*phi_basis(jphi,jfock,kfock)*Oi(ifock,jfock)
                end do
             end do
          end do
          !
       end do
    end do
    !
  end function get_traces_basis_phiOphi_z




  !
  function get_traces_basis_phiphiO(Oi) result(trace_matrix)
    real(8),dimension(nFock,nFock) :: Oi !local_operator whose traces should be computed
    complex(8),dimension(Nphi,Nphi)   :: trace_matrix
    complex(8),dimension(nFock,nFock) :: tmp
    integer                        :: kfock,iphi,jphi
    !
    trace_matrix=zero                
    do iphi=1,Nphi
       do jphi=1,Nphi
          tmp=zero
          tmp=matmul(phi_basis(jphi,:,:),Oi)
          tmp=matmul(phi_basis_dag(iphi,:,:),tmp)
          do kfock=1,nFock
             trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + tmp(kfock,kfock)                  
          end do
       end do
    end do
    !
  end function get_traces_basis_phiphiO
  !
  function get_traces_basis_phiAphiB(A,B) result(trace_matrix)
    real(8),dimension(nFock,nFock)    :: A,B 
    complex(8),dimension(Nphi,Nphi)   :: trace_matrix
    complex(8),dimension(nFock,nFock) :: tmp
    integer                           :: kfock,iphi,jphi
    !
    trace_matrix=zero
    do iphi=1,Nphi
       do jphi=1,Nphi
          tmp=zero
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
  !


  !+- SYMMETRIC OBSOLETE ROUTINES

  function get_traces_basis_phiAphiB_s(A,B) result(trace_matrix)
    real(8),dimension(nFock,nFock)    :: A,B 
    complex(8),dimension(Nphi,Nphi)   :: trace_matrix
    complex(8),dimension(nFock,nFock) :: tmp
    integer                           :: kfock,iphi,jphi
    !
    trace_matrix=zero
    do iphi=1,Nphi
       do jphi=1,Nphi
          tmp=zero
          tmp=matmul(phi_basis(jphi,:,:),B)
          tmp=matmul(A,tmp)
          tmp=matmul(phi_basis_dag(iphi,:,:),tmp)
          do kfock=1,nFock
             trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + 0.5d0*tmp(kfock,kfock)                  
          end do
          tmp=zero
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
  function get_traces_basis_phiphiO_s(Oi) result(trace_matrix)
    real(8),dimension(nFock,nFock) :: Oi !local_operator whose traces should be computed
    complex(8),dimension(Nphi,Nphi)   :: trace_matrix
    complex(8),dimension(nFock,nFock) :: tmp
    integer                        :: kfock,iphi,jphi
    !
    trace_matrix=zero                
    do iphi=1,Nphi
       do jphi=1,Nphi
          tmp=zero
          tmp=matmul(phi_basis(jphi,:,:),Oi)
          tmp=matmul(phi_basis_dag(iphi,:,:),tmp)
          do kfock=1,nFock
             trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + 0.5d0*tmp(kfock,kfock)                  
          end do
          tmp=zero
          tmp=matmul(phi_basis(iphi,:,:),Oi)
          tmp=matmul(phi_basis_dag(jphi,:,:),tmp)
          do kfock=1,nFock
             trace_matrix(iphi,jphi) = trace_matrix(iphi,jphi) + 0.5d0*tmp(kfock,kfock)                  
          end do          
       end do
    end do
    !
  end function get_traces_basis_phiphiO_s



END MODULE GZ_MATRIX_BASIS
