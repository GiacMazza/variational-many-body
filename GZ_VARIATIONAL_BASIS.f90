MODULE GZ_MATRIX_BASIS
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_AUX_FUNX
  USE SF_LINALG
  USE SF_IOTOOLS
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
  subroutine init_variational_matrices(symm,store_dir_,read_dir_)
    integer          :: symm
    logical          :: store,iread
    logical          :: store_,iread_
    character(len=200),optional :: read_dir_
    character(len=200),optional :: store_dir_
    character(len=200) :: read_dir
    character(len=200) :: store_dir
    !
    if(present(read_dir_)) then
       iread_=.true.
       if(.not.gz_superc) then
          read_dir=trim(read_dir_)//'NORB'//reg(txtfy(Norb))//'SYMM'//reg(txtfy(symm))//'_normal/'
       else
          read_dir=trim(read_dir_)//'NORB'//reg(txtfy(Norb))//'SYMM'//reg(txtfy(symm))//'_superc/'
       end if
    end if
    !
    store_=.false.
    if(present(store_dir_)) then
       store_=.true.
       if(.not.gz_superc) then
          store_dir=trim(store_dir_)//'NORB'//reg(txtfy(Norb))//'SYMM'//reg(txtfy(symm))//'_normal/'
       else
          store_dir=trim(store_dir_)//'NORB'//reg(txtfy(Norb))//'SYMM'//reg(txtfy(symm))//'_superc/'
       end if
    end if
    !
    Nphi = get_dimension_phi_basis(symm)    
    !
    if(.not.iread_) then
       call build_traces_matrix_basis    
    else
       call read_traces_matrix_basis(read_dir)    
    end if
    !
    ! call system('mkdir -v TMP_TEST_DIR')
    ! call system('mv *.gutz TMP_TEST_DIR')
    !
    ! if(.not.gz_superc) then
    !    store_dir='/homepmc/giacomo.mazza/etc_local/GZ_basis/NORB'//reg(txtfy(Norb))//'SYMM'//reg(txtfy(symm))//'_normal/'
    ! else
    !    store_dir='/homepmc/giacomo.mazza/etc_local/GZ_basis/NORB'//reg(txtfy(Norb))//'SYMM'//reg(txtfy(symm))//'_superc/'
    ! end if
    if(store_) call store_traces_matrix_basis(store_dir)
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
  function get_dimension_phi_basis(symm_index) result(dim_phi)
    integer                                         :: symm_index
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
    
    select case(symm_index)
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
    !dim_phi = nFock
    Nphi=dim_phi   
    write(*,*) "NPHI",Nphi
    !
    allocate(phi_basis(dim_phi,nFock,nFock))
    allocate(phi_basis_dag(dim_phi,nFock,nFock))

    phi_basis = phi_fock


    !< TMP_TEST_DIAGONAL_BASIS
    ! phi_basis=0.d0
    ! iphi=0
    ! do ifock=1,nFock
    !    iphi=iphi+1
    !    phi_basis(iphi,iphi,iphi) = 1.d0
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

    ! allocate(test_trace(Nphi,Nphi))
    ! test_trace=get_traces_basis_phiOphi(Id)

    ! !+- safe check on the matrix orthogonality
    ! do iphi=1,Nphi
    !    do jphi=1,Nphi
    !       if(abs(test_trace(iphi,jphi)).gt.1.d-8.and.iphi.ne.jphi) &
    !            stop "Variational basis not Trace-orthogonal"
    !    end do
    ! end do

    ! do iphi=1,Nphi
    !    phi_basis(iphi,:,:)     = phi_basis(iphi,:,:)/sqrt(test_trace(iphi,iphi))
    !    phi_basis_dag(iphi,:,:) = phi_basis_dag(iphi,:,:)/sqrt(test_trace(iphi,iphi))
    ! end do

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






  subroutine store_traces_matrix_basis(store_dir)
    integer          :: is,js,iorb,jorb,iphi,jphi,unit_store
    logical          :: store
    real(8)          :: x
    character(len=200) :: file_name
    character(len=200) :: store_dir

    store_dir=trim(store_dir)
    call system('mkdir -p '//store_dir)

    do is=1,Ns
       do js=1,Ns
          unit_store=free_unit()
          file_name=reg(store_dir)//"phi_traces_basis_local_dens_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          open(unit_store,file=file_name)
          do iphi=1,Nphi
             do jphi=1,Nphi
                x=conjg(phi_traces_basis_local_dens(is,js,iphi,jphi))*phi_traces_basis_local_dens(is,js,iphi,jphi)
                if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_local_dens(is,js,iphi,jphi)
             end do
          end do
          close(unit_store)
          unit_store=free_unit()
          file_name=reg(store_dir)//"phi_traces_basis_dens_dens_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          open(unit_store,file=file_name)
          do iphi=1,Nphi
             do jphi=1,Nphi
                x=conjg(phi_traces_basis_dens_dens(is,js,iphi,jphi))*phi_traces_basis_dens_dens(is,js,iphi,jphi)
                if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_dens_dens(is,js,iphi,jphi)
             end do
          end do
          close(unit_store)
       end do
    end do

    !
    do iorb=1,Norb
       do jorb=1,Norb
          !
          ! unit_store=free_unit()
          ! file_name=reg(store_dir)//"phi_traces_basis_dens_dens_orb_iorb"//reg(txtfy(iorb))//"_jorb"//reg(txtfy(jorb))//".gutz"             
          ! open(unit_store,file=file_name)
          ! do iphi=1,Nphi
          !    do jphi=1,Nphi
          !       x=conjg(phi_traces_basis_dens_dens_orb(iorb,jorb,iphi,jphi))*phi_traces_basis_dens_dens_orb(iorb,jorb,iphi,jphi)
          !       if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_dens_dens_orb(iorb,jorb,iphi,jphi)
          !    end do
          ! end do
          ! close(unit_store)
          !
          !+- Orbital spin-flip -+!
          unit_store=free_unit()
          file_name=reg(store_dir)//"phi_traces_basis_spin_flip_iorb"//reg(txtfy(iorb))//"_jorb"//reg(txtfy(jorb))//".gutz"             
          open(unit_store,file=file_name)
          do iphi=1,Nphi
             do jphi=1,Nphi
                x=conjg(phi_traces_basis_spin_flip(iorb,jorb,iphi,jphi))*phi_traces_basis_spin_flip(iorb,jorb,iphi,jphi)
                if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_spin_flip(iorb,jorb,iphi,jphi)
             end do
          end do
          close(unit_store)
          !

          !+- Orbital pair-hopping -+!
          unit_store=free_unit()
          file_name=reg(store_dir)//"phi_traces_basis_pair_hopping_iorb"//reg(txtfy(iorb))//"_jorb"//reg(txtfy(jorb))//".gutz"             
          open(unit_store,file=file_name)
          do iphi=1,Nphi
             do jphi=1,Nphi
                x=conjg(phi_traces_basis_pair_hopping(iorb,jorb,iphi,jphi))*phi_traces_basis_pair_hopping(iorb,jorb,iphi,jphi)
                if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_pair_hopping(iorb,jorb,iphi,jphi)
             end do
          end do
          close(unit_store)
          !          
       end do
       !+- orbital double occupancy -+!
       ! unit_store=free_unit()
       ! file_name=reg(store_dir)//"phi_traces_basis_docc_orb_iorb"//reg(txtfy(iorb))//".gutz"             
       ! open(unit_store,file=file_name)
       ! do iphi=1,Nphi
       !    do jphi=1,Nphi
       !       x=conjg(phi_traces_basis_docc_orb(iorb,iphi,jphi))*phi_traces_basis_docc_orb(iorb,iphi,jphi)
       !       if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_docc_orb(iorb,iphi,jphi)
       !    end do
       ! end do
       ! close(unit_store)
    end do
    !

    !+- density constraints Tr(Phi+ Phi C_i) C_i=density matrix -+!
    do is=1,Ns
       do js=1,Ns
          unit_store=free_unit()
          file_name=reg(store_dir)//"phi_traces_basis_dens_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          open(unit_store,file=file_name)
          do iphi=1,Nphi
             do jphi=1,Nphi
                x=conjg(phi_traces_basis_dens(is,js,iphi,jphi))*phi_traces_basis_dens(is,js,iphi,jphi)
                if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_dens(is,js,iphi,jphi)
             end do
          end do
          close(unit_store)
       end do
    end do

    !+- anomalous part -+!
    do is=1,Ns
       do js=1,Ns
          unit_store=free_unit()
          file_name=reg(store_dir)//"phi_traces_basis_dens_anomalous_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          open(unit_store,file=file_name)
          do iphi=1,Nphi
             do jphi=1,Nphi
                x=conjg(phi_traces_basis_dens_anomalous(is,js,iphi,jphi))*phi_traces_basis_dens_anomalous(is,js,iphi,jphi)
                if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_dens_anomalous(is,js,iphi,jphi)
             end do
          end do
          close(unit_store)
       end do
    end do
    ! do iphi=1,Nphi
    !    do jphi=1,Nphi
    !       phi_traces_basis_dens_anomalous_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_dens_anomalous(:,:,jphi,iphi))
    !    end do
    ! end do

    !+- Hoppings -+!
    do is=1,Ns
       do js=1,Ns
          unit_store=free_unit()
          file_name=reg(store_dir)//"phi_traces_basis_Rhop_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          open(unit_store,file=file_name)
          do iphi=1,Nphi
             do jphi=1,Nphi
                x=conjg(phi_traces_basis_Rhop(is,js,iphi,jphi))*phi_traces_basis_Rhop(is,js,iphi,jphi)
                if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_Rhop(is,js,iphi,jphi)
             end do
          end do
          close(unit_store)
       end do
    end do
    ! do iphi=1,Nphi
    !    do jphi=1,Nphi
    !       phi_traces_basis_Rhop_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_Rhop(:,:,jphi,iphi))
    !    end do
    ! end do
    !
    if(gz_superc) then
       do is=1,Ns
          do js=1,Ns
             unit_store=free_unit()
             file_name=reg(store_dir)//"phi_traces_basis_Qhop_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
             open(unit_store,file=file_name)
             do iphi=1,Nphi
                do jphi=1,Nphi
                   x=conjg(phi_traces_basis_Qhop(is,js,iphi,jphi))*phi_traces_basis_Qhop(is,js,iphi,jphi)
                   if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_Qhop(is,js,iphi,jphi)
                end do
             end do
             close(unit_store)
          end do
       end do
       ! do iphi=1,Nphi
       !    do jphi=1,Nphi
       !       phi_traces_basis_Qhop_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_Qhop(:,:,jphi,iphi))
       !    end do
       ! end do
    end if
    !
    do is=1,Ns
       do js=1,Ns
          unit_store=free_unit()
          file_name=reg(store_dir)//"phi_traces_basis_sc_order_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          open(unit_store,file=file_name)
          do iphi=1,Nphi
             do jphi=1,Nphi
                x=conjg(phi_traces_basis_sc_order(is,js,iphi,jphi))*phi_traces_basis_sc_order(is,js,iphi,jphi)
                if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_sc_order(is,js,iphi,jphi)
             end do
          end do
          close(unit_store)
       end do
    end do
    !
    unit_store=free_unit()
    file_name=reg(store_dir)//"phi_traces_basis_spin2.gutz"
    open(unit_store,file=file_name)
    do iphi=1,Nphi
       do jphi=1,Nphi
          x=conjg(phi_traces_basis_spin2(iphi,jphi))*phi_traces_basis_spin2(iphi,jphi)
          if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_spin2(iphi,jphi)
       end do
    end do
    close(unit_store)
    unit_store=free_unit()
    file_name=reg(store_dir)//"phi_traces_basis_spinZ.gutz"
    open(unit_store,file=file_name)
    do iphi=1,Nphi
       do jphi=1,Nphi
          x=conjg(phi_traces_basis_spinZ(iphi,jphi))*phi_traces_basis_spinZ(iphi,jphi)
          if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_spinZ(iphi,jphi)
       end do
    end do
    close(unit_store)
    unit_store=free_unit()
    file_name=reg(store_dir)//"phi_traces_basis_isospin2.gutz"
    open(unit_store,file=file_name)
    do iphi=1,Nphi
       do jphi=1,Nphi
          x=conjg(phi_traces_basis_isospin2(iphi,jphi))*phi_traces_basis_isospin2(iphi,jphi)
          if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_isospin2(iphi,jphi)
       end do
    end do
    close(unit_store)
    unit_store=free_unit()
    file_name=reg(store_dir)//"phi_traces_basis_isospinZ.gutz"
    open(unit_store,file=file_name)
    do iphi=1,Nphi
       do jphi=1,Nphi
          x=conjg(phi_traces_basis_isospinZ(iphi,jphi))*phi_traces_basis_isospinZ(iphi,jphi)
          if(x.gt.1.d-12) write(unit_store,'(2I4,2F18.10)') iphi,jphi,phi_traces_basis_isospinZ(iphi,jphi)
       end do
    end do
    close(unit_store)
    !
  end subroutine store_traces_matrix_basis





  subroutine read_traces_matrix_basis(read_dir)
    integer          :: is,js,iorb,jorb,iphi,jphi,unit_store,ios
    integer :: ifile
    logical          :: store,read_file
    real(8)          :: x,re_x,im_x
    character(len=200) :: file_name
    character(len=200) :: read_dir
    integer :: flen

    read_dir=trim(read_dir)

    allocate(phi_traces_basis_local_dens(Ns,Ns,Nphi,Nphi));phi_traces_basis_local_dens=zero
    !allocate(phi_traces_basis_dens_dens_orb(Norb,Norb,Nphi,Nphi));phi_traces_basis_dens_dens_orb=zero
    !allocate(phi_traces_basis_docc_orb(Norb,Nphi,Nphi));phi_traces_basis_docc_orb=zero
    allocate(phi_traces_basis_dens_dens(Ns,Ns,Nphi,Nphi));phi_traces_basis_dens_dens=zero
    allocate(phi_traces_basis_spin_flip(Norb,Norb,Nphi,Nphi));phi_traces_basis_spin_flip=zero
    allocate(phi_traces_basis_pair_hopping(Norb,Norb,Nphi,Nphi));phi_traces_basis_pair_hopping=zero
    allocate(phi_traces_basis_sc_order(Ns,Ns,Nphi,Nphi));phi_traces_basis_sc_order=zero
    !
    allocate(phi_traces_basis_spin2(Nphi,Nphi));phi_traces_basis_spin2=zero
    allocate(phi_traces_basis_spinZ(Nphi,Nphi));phi_traces_basis_spinZ=zero
    allocate(phi_traces_basis_isospin2(Nphi,Nphi));phi_traces_basis_isospin2=zero
    allocate(phi_traces_basis_isospinZ(Nphi,Nphi));phi_traces_basis_isospinZ=zero
    !
    write(*,*) "READING PHI TRACES"
    !
    write(*,*)
    write(*,*) 
    do is=1,Ns
       do js=1,Ns
          !
          write(*,*) "READING local_dens",is,js
          !
          file_name=reg(read_dir)//"phi_traces_basis_local_dens_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          inquire(file=file_name,exist=read_file)
          if(read_file) then
             !
             unit_store=free_unit()
             open(unit_store,file=file_name,status='old')
             flen=0
             do
                read (unit_store,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(unit_store)          
             open(unit_store,file=file_name,status='old')
             do ifile=1,flen
                read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
                phi_traces_basis_local_dens(is,js,iphi,jphi) = re_x+xi*im_x
             end do
             close(unit_store)
             !
          else
             write(*,*) 'file',file_name,'not stored'
             stop 
          end if

          write(*,*)
          write(*,*) "READING local_dens_dens",is,js

          file_name=reg(read_dir)//"phi_traces_basis_dens_dens_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          inquire(file=file_name,exist=read_file)
          if(read_file) then
             !
             unit_store=free_unit()
             open(unit_store,file=file_name,status='old')
             flen=0
             do
                read (unit_store,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(unit_store)          
             open(unit_store,file=file_name,status='old')
             do ifile=1,flen
                read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
                phi_traces_basis_dens_dens(is,js,iphi,jphi) = re_x+xi*im_x
             end do
             close(unit_store)
             !
          else
             write(*,*) 'file',file_name,'not stored'
             stop 
          end if
       end do
    end do






    !
    do iorb=1,Norb
       do jorb=1,Norb
          !
          !
          !+- Orbital spin-flip -+!
          file_name=reg(read_dir)//"phi_traces_basis_spin_flip_iorb"//reg(txtfy(iorb))//"_jorb"//reg(txtfy(jorb))//".gutz"             
          write(*,*)
          write(*,*) "READING spin_flip",iorb,jorb          
          inquire(file=file_name,exist=read_file)
          if(read_file) then
             !
             unit_store=free_unit()
             open(unit_store,file=file_name,status='old')
             flen=0
             do
                read (unit_store,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(unit_store)          
             open(unit_store,file=file_name,status='old')
             do ifile=1,flen
                read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
                phi_traces_basis_spin_flip(iorb,jorb,iphi,jphi) = re_x+xi*im_x
             end do
             close(unit_store)
             !
          else
             write(*,*) 'file',file_name,'not stored'
             stop 
          end if
          !

          !+- Orbital pair-hopping -+!
          write(*,*)
          write(*,*) "READING pair hopping",iorb,jorb          
          file_name=reg(read_dir)//"phi_traces_basis_pair_hopping_iorb"//reg(txtfy(iorb))//"_jorb"//reg(txtfy(jorb))//".gutz"             
          inquire(file=file_name,exist=read_file)
          if(read_file) then
             !
             unit_store=free_unit()
             open(unit_store,file=file_name,status='old')
             flen=0
             do
                read (unit_store,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(unit_store)          
             open(unit_store,file=file_name,status='old')
             do ifile=1,flen
                read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
                phi_traces_basis_pair_hopping(iorb,jorb,iphi,jphi) = re_x+xi*im_x
             end do
             close(unit_store)
             !
          else
             write(*,*) 'file',file_name,'not stored'
             stop 
          end if
          !          
       end do
       ! !+- orbital double occupancy -+!
       ! file_name=reg(read_dir)//"phi_traces_basis_docc_orb_iorb"//reg(txtfy(iorb))//".gutz"             
       ! inquire(file=file_name,exist=read_file)
       ! if(read_file) then
       !    !
       !    unit_store=free_unit()
       !    open(unit_store,file=file_name,status='old')
       !    flen=0
       !    do
       !       read (unit_store,*,iostat=ios) x
       !       if (ios/=0) exit     
       !       flen=flen+1
       !    end do
       !    close(unit_store)          
       !    open(unit_store,file=file_name,status='old')
       !    do ifile=1,flen
       !       read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
       !       phi_traces_basis_docc_orb(iorb,iphi,jphi) = re_x+xi*im_x
       !    end do
       !    close(unit_store)
       !    !
       ! else
       !    write(*,*) 'file',file_name,'not stored'
       !    stop 
       ! end if
    end do
    !

    !+- density constraints Tr(Phi+ Phi C_i) C_i=density matrix -+!
    allocate(phi_traces_basis_dens(Ns,Ns,Nphi,Nphi),phi_traces_basis_dens_hc(Ns,Ns,Nphi,Nphi))
    phi_traces_basis_dens=zero
    phi_traces_basis_dens_hc=zero
    do is=1,Ns
       do js=1,Ns
          write(*,*)
          write(*,*) "READING density constraint",is,js          
          file_name=reg(read_dir)//"phi_traces_basis_dens_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          inquire(file=file_name,exist=read_file)
          if(read_file) then
             !
             unit_store=free_unit()
             open(unit_store,file=file_name,status='old')
             flen=0
             do
                read (unit_store,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(unit_store)          
             open(unit_store,file=file_name,status='old')
             do ifile=1,flen
                read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
                phi_traces_basis_dens(is,js,iphi,jphi) = re_x+xi*im_x
             end do
             close(unit_store)
             !
          else
             write(*,*) 'file',file_name,'not stored'
             stop 
          end if
       end do
    end do
    do iphi=1,Nphi
       do jphi=1,Nphi
          phi_traces_basis_dens_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_dens(:,:,jphi,iphi))
       end do
    end do
    !
    !+- anomalous part -+!
    allocate(phi_traces_basis_dens_anomalous(Ns,Ns,Nphi,Nphi));phi_traces_basis_dens_anomalous=zero
    allocate(phi_traces_basis_dens_anomalous_hc(Ns,Ns,Nphi,Nphi));phi_traces_basis_dens_anomalous_hc=zero
    !
    do is=1,Ns
       do js=1,Ns
          write(*,*)
          write(*,*) "READING density constraint anomalous",is,js          
          file_name=reg(read_dir)//"phi_traces_basis_dens_anomalous_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          inquire(file=file_name,exist=read_file)
          if(read_file) then
             !
             unit_store=free_unit()
             open(unit_store,file=file_name,status='old')
             flen=0
             do
                read (unit_store,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(unit_store)          
             open(unit_store,file=file_name,status='old')
             do ifile=1,flen
                read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
                phi_traces_basis_dens_anomalous(is,js,iphi,jphi) = re_x+xi*im_x
             end do
             close(unit_store)
             !
          else
             write(*,*) 'file',file_name,'not stored'
             stop 
          end if
       end do
    end do
    do iphi=1,Nphi
       do jphi=1,Nphi
          phi_traces_basis_dens_anomalous_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_dens_anomalous(:,:,jphi,iphi))
       end do
    end do

    !+- Hoppings -+!
    allocate(phi_traces_basis_Rhop(Ns,Ns,Nphi,Nphi),phi_traces_basis_Rhop_hc(Ns,Ns,Nphi,Nphi))
    phi_traces_basis_Rhop=zero;phi_traces_basis_Rhop_hc=zero
    do is=1,Ns
       do js=1,Ns
          write(*,*)
          write(*,*) "READING Rhop",is,js          

          file_name=reg(read_dir)//"phi_traces_basis_Rhop_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          inquire(file=file_name,exist=read_file)
          if(read_file) then
             !
             unit_store=free_unit()
             open(unit_store,file=file_name,status='old')
             flen=0
             do
                read (unit_store,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(unit_store)          
             open(unit_store,file=file_name,status='old')
             do ifile=1,flen
                read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
                phi_traces_basis_Rhop(is,js,iphi,jphi) = re_x+xi*im_x
             end do
             close(unit_store)
             !
          else
             write(*,*) 'file',file_name,'not stored'
             stop 
          end if
       end do
    end do
    do iphi=1,Nphi
       do jphi=1,Nphi
          phi_traces_basis_Rhop_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_Rhop(:,:,jphi,iphi))
       end do
    end do
    ! !
    if(gz_superc) then
       allocate(phi_traces_basis_Qhop(Ns,Ns,Nphi,Nphi),phi_traces_basis_Qhop_hc(Ns,Ns,Nphi,Nphi))
       phi_traces_basis_Qhop=zero;phi_traces_basis_Qhop_hc=zero
       do is=1,Ns
          do js=1,Ns

             write(*,*)
             write(*,*) "READING Qhop",is,js          


             file_name=reg(read_dir)//"phi_traces_basis_Qhop_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
             inquire(file=file_name,exist=read_file)
             if(read_file) then
                !
                unit_store=free_unit()
                open(unit_store,file=file_name,status='old')
                flen=0
                do
                   read (unit_store,*,iostat=ios) x
                   if (ios/=0) exit     
                   flen=flen+1
                end do
                close(unit_store)          
                open(unit_store,file=file_name,status='old')
                do ifile=1,flen
                   read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
                   phi_traces_basis_Qhop(is,js,iphi,jphi) = re_x+xi*im_x
                end do
                close(unit_store)
                !
             else
                write(*,*) 'file',file_name,'not stored'
                stop 
             end if
          end do
       end do
       do iphi=1,Nphi
          do jphi=1,Nphi
             phi_traces_basis_Qhop_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_Qhop(:,:,jphi,iphi))
          end do
       end do
    end if
    ! !
    do is=1,Ns
       do js=1,Ns

          write(*,*)
          write(*,*) "READING SC_order_param",is,js          


          file_name=reg(read_dir)//"phi_traces_basis_sc_order_is"//reg(txtfy(is))//"_js"//reg(txtfy(js))//".gutz"             
          inquire(file=file_name,exist=read_file)
          if(read_file) then
             !
             unit_store=free_unit()
             open(unit_store,file=file_name,status='old')
             flen=0
             do
                read (unit_store,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(unit_store)          
             open(unit_store,file=file_name,status='old')
             do ifile=1,flen
                read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
                phi_traces_basis_sc_order(is,js,iphi,jphi) = re_x+xi*im_x
             end do
             close(unit_store)
             !
          else
             write(*,*) 'file',file_name,'not stored'
             stop 
          end if
       end do
    end do
    ! 

    write(*,*)
    write(*,*) "READING ANGULAR MOMENTA"


    file_name=reg(read_dir)//"phi_traces_basis_spin2.gutz"    
    inquire(file=file_name,exist=read_file)
    if(read_file) then
       !
       unit_store=free_unit()
       open(unit_store,file=file_name,status='old')
       flen=0
       do
          read (unit_store,*,iostat=ios) x
          if (ios/=0) exit     
          flen=flen+1
       end do
       close(unit_store)          
       open(unit_store,file=file_name,status='old')
       do ifile=1,flen
          read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
          phi_traces_basis_spin2(iphi,jphi) = re_x+xi*im_x
       end do
       close(unit_store)
       !
    else
       write(*,*) 'file',file_name,'not stored'
       stop 
    end if
    !
    file_name=reg(read_dir)//"phi_traces_basis_spinZ.gutz"
    inquire(file=file_name,exist=read_file)
    if(read_file) then
       !
       unit_store=free_unit()
       open(unit_store,file=file_name,status='old')
       flen=0
       do
          read (unit_store,*,iostat=ios) x
          if (ios/=0) exit     
          flen=flen+1
       end do
       close(unit_store)          
       open(unit_store,file=file_name,status='old')
       do ifile=1,flen
          read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
          phi_traces_basis_spinZ(iphi,jphi) = re_x+xi*im_x
       end do
       close(unit_store)
       !
    else
       write(*,*) 'file',file_name,'not stored'
       stop 
    end if
    !
    file_name=reg(read_dir)//"phi_traces_basis_isospin2.gutz"
    inquire(file=file_name,exist=read_file)
    if(read_file) then
       !
       unit_store=free_unit()
       open(unit_store,file=file_name,status='old')
       flen=0
       do
          read (unit_store,*,iostat=ios) x
          if (ios/=0) exit     
          flen=flen+1
       end do
       close(unit_store)          
       open(unit_store,file=file_name,status='old')
       do ifile=1,flen
          read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
          phi_traces_basis_isospin2(iphi,jphi) = re_x+xi*im_x
       end do
       close(unit_store)
       !
    else
       write(*,*) 'file',file_name,'not stored'
       stop 
    end if
    !+- 
    file_name=reg(read_dir)//"phi_traces_basis_isospinZ.gutz"
    inquire(file=file_name,exist=read_file)
    if(read_file) then
       !
       unit_store=free_unit()
       open(unit_store,file=file_name,status='old')
       flen=0
       do
          read (unit_store,*,iostat=ios) x
          if (ios/=0) exit     
          flen=flen+1
       end do
       close(unit_store)          
       open(unit_store,file=file_name,status='old')
       do ifile=1,flen
          read(unit_store,'(2I4,2F18.10)') iphi,jphi,re_x,im_x
          phi_traces_basis_isospinZ(iphi,jphi) = re_x+xi*im_x
       end do
       close(unit_store)
       !
    else
       write(*,*) 'file',file_name,'not stored'
       stop 
    end if
    !
  end subroutine read_traces_matrix_basis
  


  



  subroutine build_traces_matrix_basis
    integer          :: is,js,iorb,jorb,iphi,jphi,unit_store
    real(8)          :: x
    character(len=50) :: file_name
    !+- let's assume that the opertaors in Fock space have been already built

    !+- local operators Tr(Phi+ Oi Phi) -+!
    allocate(phi_traces_basis_local_dens(Ns,Ns,Nphi,Nphi))
    !+- this new guy contains all phi_traces_basis_ndens_ndens
    allocate(phi_traces_basis_dens_dens(Ns,Ns,Nphi,Nphi))
    !
    allocate(phi_traces_basis_spin_flip(Norb,Norb,Nphi,Nphi))
    allocate(phi_traces_basis_pair_hopping(Norb,Norb,Nphi,Nphi))
    allocate(phi_traces_basis_sc_order(Ns,Ns,Nphi,Nphi))
    !
    allocate(phi_traces_basis_spin2(Nphi,Nphi))
    allocate(phi_traces_basis_spinZ(Nphi,Nphi))
    allocate(phi_traces_basis_isospin2(Nphi,Nphi))
    allocate(phi_traces_basis_isospinZ(Nphi,Nphi))

    write(*,*) 'BUILDING PHI TRACES: local density matrix'
    do is=1,Ns
       do js=1,Ns
          write(*,*) 'IS',is,'JS',js
          phi_traces_basis_local_dens(is,js,:,:) = &
               get_traces_basis_phiOphi(op_local_dens(is,js,:,:))
       end do
    end do
    !
    write(*,*) 'BUILDING PHI TRACES: local density-density'
    write(*,*) 
    do is=1,Ns
       do js=1,Ns
          write(*,*) 'IS',is,'JS',js
          phi_traces_basis_dens_dens(is,js,:,:) = &
               get_traces_basis_phiOphi(op_dens_dens(is,js,:,:))
       end do
    end do
    !
    write(*,*) 'BUILDING PHI TRACES: local spin-flip and pair-hopping'
    do iorb=1,Norb
       do jorb=1,Norb
          !+- Orbital spin-flip -+!
          write(*,*) 'ORBITAL SPIN-FLIP','IORB',iorb,'JORB',jorb
          phi_traces_basis_spin_flip(iorb,jorb,:,:) = &
               get_traces_basis_phiOphi(op_spin_flip(iorb,jorb,:,:))
          !+- Orbital pair-hopping -+!
          write(*,*) 'ORBITAL PAIR-HOPPING','IORB',iorb,'JORB',jorb
          phi_traces_basis_pair_hopping(iorb,jorb,:,:) = &
               get_traces_basis_phiOphi(op_pair_hopping(iorb,jorb,:,:))
       end do
    end do
    !+- density constraints Tr(Phi+ Phi C_i) C_i=density matrix -+!
    write(*,*) 'BUILDING PHI TRACES: density constraints'
    allocate(phi_traces_basis_dens(Ns,Ns,Nphi,Nphi),phi_traces_basis_dens_hc(Ns,Ns,Nphi,Nphi))
    do is=1,Ns
       do js=1,Ns
          write(*,*) 'IS',is,'JS',js
          phi_traces_basis_dens(is,js,:,:) = &    
               get_traces_basis_phiphiO(op_local_dens(is,js,:,:))
       end do
    end do
    do iphi=1,Nphi
       do jphi=1,Nphi
          phi_traces_basis_dens_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_dens(:,:,jphi,iphi))
       end do
    end do
    !

    !+- anomalous part -+!
    allocate(phi_traces_basis_dens_anomalous(Ns,Ns,Nphi,Nphi))
    allocate(phi_traces_basis_dens_anomalous_hc(Ns,Ns,Nphi,Nphi))
    write(*,*) 'BUILDING PHI TRACES: anomalous density constraints'
    do is=1,Ns
       do js=1,Ns
          write(*,*) 'IS',is,'JS',js
          phi_traces_basis_dens_anomalous(is,js,:,:) = &    !+- probably is more correct to call this global variable !+- phi_traces_basis_vdm
               get_traces_basis_phiphiO(op_local_dens_anomalous(is,js,:,:))
       end do
    end do
    do iphi=1,Nphi
       do jphi=1,Nphi
          phi_traces_basis_dens_anomalous_hc(:,:,iphi,jphi) = conjg(phi_traces_basis_dens_anomalous(:,:,jphi,iphi))
       end do
    end do

    !+- Hoppings -+!
    write(*,*) 'BUILDING PHI TRACES: normal hopping renormalization'
    allocate(phi_traces_basis_Rhop(Ns,Ns,Nphi,Nphi),phi_traces_basis_Rhop_hc(Ns,Ns,Nphi,Nphi))
    do is=1,Ns
       do js=1,Ns
          write(*,*) 'IS',is,'JS',js
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
       write(*,*) 'BUILDING PHI TRACES: anomalous hopping renormalization'
       allocate(phi_traces_basis_Qhop(Ns,Ns,Nphi,Nphi),phi_traces_basis_Qhop_hc(Ns,Ns,Nphi,Nphi))
       do is=1,Ns
          do js=1,Ns
             write(*,*) 'IS',is,'JS',js
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
    write(*,*) 'BUILDING PHI TRACES: local order parameter'
    do is=1,Ns
       do js=1,Ns
          write(*,*) 'IS',is,'JS',js
          phi_traces_basis_sc_order(is,js,:,:) = get_traces_basis_phiOphi(op_sc_order(is,js,:,:)) 
       end do
    end do
    !
    write(*,*) 'BUILDING PHI TRACES: local SPIN2'
    phi_traces_basis_spin2 = get_traces_basis_phiOphi_z(op_spin2) 
    write(*,*) 'BUILDING PHI TRACES: local SPIN_z'
    phi_traces_basis_spinZ = get_traces_basis_phiOphi_z(op_spinZ) 
    write(*,*) 'BUILDING PHI TRACES: local ISOSPIN2'
    phi_traces_basis_isospin2 = get_traces_basis_phiOphi_z(op_isospin2) 
    write(*,*) 'BUILDING PHI TRACES: local ISOSPIN_z'
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
