program GUTZ_mb
  USE SCIFOR
  !
  USE DMFT_MISC
  USE SF_PARSE_INPUT
  !
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_LOCAL_HAMILTONIAN
  USE GZ_OPTIMIZED_ENERGY
  USE GZ_ENERGY_MINIMIZATION
  USE GZ_EFFECTIVE_HOPPINGS
  !
  USE GZ_MATRIX_BASIS
  !
  implicit none
  !
  !+- hamiltonian details -+!
  integer                            :: ispin,iorb,i,istate,jstate,ifock,jorb
  integer,dimension(:),allocatable   :: fock_vec
  complex(8),dimension(:),allocatable               :: init_vec
  real(8),dimension(:),allocatable   :: variational_density_natural
  real(8),dimension(:),allocatable   :: vdm_init,vdm_out
  complex(8),dimension(:),allocatable   :: R_init,R_out
  complex(8),dimension(:),allocatable   :: Rhop_init,Rhop_out
  complex(8),dimension(:,:),allocatable   :: Rhop_init_matrix
  real(8),dimension(:,:),allocatable :: variational_density_natural_simplex
  real(8),dimension(:,:),allocatable :: variational_density_matrix
  integer                            :: out_unit,iter
  integer                            :: lattice ! 2=square;3=cubic
  real(8),dimension(:),allocatable :: epsik,hybik
  integer :: Nx,is
  !
  real(8) :: tmp_emin,Uiter,tmp_ene,orb_pol,Jh_ratio,Jiter
  complex(8),dimension(:,:),allocatable :: slater_lgr_init,gzproj_lgr_init
  integer :: Nopt
  character(len=5) :: dir_suffix
  character(len=6) :: dir_iter
  real(8),dimension(:),allocatable :: dump_seed
  integer :: expected_flen,flen,unit
  logical :: seed_file

  !
  character(len=6) :: sweep !sweepU,sweepJ
  real(8) :: sweep_start,sweep_stop,sweep_step
  integer ::  Nsweep,iseed
  !
  character(len=200) :: store_dir,read_dir
  real(8),dimension(:),allocatable :: energy_levels
  !
  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=10)
  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"Wband","inputGZ.conf",default=1.d0)
  call parse_input_variable(sweep,"SWEEP","inputGZ.conf",default='sweepJ')  
  call parse_input_variable(sweep_start,"SWEEP_START","inputGZ.conf",default=-0.4d0)  
  call parse_input_variable(sweep_stop,"SWEEP_STOP","inputGZ.conf",default=0.d0)  
  call parse_input_variable(sweep_step,"SWEEP_STEP","inputGZ.conf",default=0.05d0)  
  call parse_input_variable(read_dir,"READ_DIR","inputGZ.conf",default='~/etc_local/GZ_basis/')
  call parse_input_variable(store_dir,"STORE_DIR","inputGZ.conf",default='./READ_PHI_TRACES/')  
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")

  if(Norb.eq.1.and.wf_symmetry.eq.1) then
     write(*,*) 'WARNING THE O(1) x SU(2)c x ORBITAL_ROTATION = O(1) x SU(2)c for the Norb=1 case!'
     wf_symmetry=0
  end if


  !NOTE: ON HUNDS COUPLINGS: (see arXiv 1207.3033)
  !
  !NORB=3 ROTATIONAL INVARIANT HAMILTONIAN       :: Jsf=Jh, Jph=U-Ust-J   (NO relation between Ust and U)
  !       FULLY ROTATIONAL INVARIANT HAMILTONIAN :: Jsf=Jh, Jph=J, Ust = U - 2J   
  !
  !NORB=2 FULLY ROTATIONAL INVARIANCE            :: Jsf=Jh, Jph=0, Ust = U - J  
  !       PARTIAL O(2) INVARIANCE                :: Jsf=Jh, Jph=0, Ust = U
  !       ONLY SPIN SU(2) INVARIANCE             :: Jsf=Jh, Jph=J, Ust = U - 2J  
  ! Jh = Jh*Uloc(1)
  ! Jsf = Jh
  ! Jph = Jh
  ! Ust = Uloc(1)-2.d0*Jh
  !
  call initialize_local_fock_space
  !
  !  
  call init_variational_matrices(wf_symmetry,store_dir_=store_dir,read_dir_=read_dir)  
  allocate(energy_levels(Ns))
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        energy_levels(istate) = Cfield*0.5d0
        if(iorb.eq.2) energy_levels(istate) = -Cfield*0.5d0
     end do
  end do
  !
  call build_lattice_model

  !  
  allocate(variational_density_natural_simplex(Ns+1,Ns))
  allocate(variational_density_natural(Ns))

  !


  !
  NRhop_opt=2;   Rhop_stride_v2m => Rhop_vec2mat; Rhop_stride_m2v => Rhop_mat2vec 
  Nvdm_NC_opt=2; vdm_NC_stride_v2m => vdm_NC_vec2mat ; vdm_NC_stride_m2v => vdm_NC_mat2vec
  Nvdm_NCoff_opt=0; vdm_NCoff_stride_v2m => vdm_NCoff_vec2mat ; vdm_NCoff_stride_m2v => vdm_NCoff_mat2vec
  Nopt = NRhop_opt + Nvdm_NC_opt + Nvdm_NCoff_opt
  Nopt = 2*Nopt
  !
  allocate(Rhop_init_matrix(Ns,Ns)); 
  allocate(slater_lgr_init(Ns,Ns),gzproj_lgr_init(Ns,Ns))  
  !
  Nopt_reduced = 1 + 1 + 1 + 1
  stride_zeros_orig2red => stride2reduced
  stride_zeros_red2orig => stride2orig
  !
  !
  !
  expected_flen=Nopt
  inquire(file="Rn0_root_seed.conf",exist=seed_file)
  if(seed_file) then
     flen=file_length("Rn0_root_seed.conf")
     unit=free_unit()
     open(unit,file="Rn0_root_seed.conf")
     write(*,*) 'reading root seed from file RQn0_root_seed.conf'
     if(flen.eq.expected_flen) then
        allocate(dump_seed(flen))
        !+- read from file -+!
        do i=1,flen
           read(unit,*) dump_seed(i)
        end do
        !+------------------+!        
        call dump2mats(dump_seed,Rhop_init_matrix,slater_lgr_init,gzproj_lgr_init)        
     else
        write(*,*) 'RQn0_root_seed.conf in the wrong form',flen,expected_flen
        write(*,*) 'please check your initialization file for the root finding'
        stop
     end if
     !
  else
     !
     call init_Rhop_seed(Rhop_init_matrix)
     !
     slater_lgr_init=0.d0
     slater_lgr_init(1,1) =  Cfield*0.5d0
     slater_lgr_init(2,2) = -Cfield*0.5d0
     gzproj_lgr_init=zero
     !
  end if

  do is=1,Ns
     write(*,*) Rhop_init_matrix(is,:)
  end do
  

  !
  allocate(vdm_init(Ns))
  !call initialize_variational_density(vdm_init)
  vdm_init=0.5d0
  allocate(variational_density_matrix(Ns,Ns)); variational_density_matrix=0.d0
  do is=1,Ns
     variational_density_matrix(is,is) = vdm_init(is)
  end do
  !
  select case(sweep)
  case('sweepJ')
     !+- sweep JHund -+!
     Nsweep = abs(sweep_start-sweep_stop)/sweep_step
     Jiter = sweep_start
     do i=1,35
        !
        Jh = 0.d0
        Jsf = 0.d0
        Jph = Jh
        !
        call get_local_hamiltonian_trace(energy_levels)
        !
        write(dir_suffix,'(F4.2)') Jiter
        dir_iter="J"//trim(dir_suffix)
        call system('mkdir -v '//dir_iter)     
        !
        call gz_optimization_vdm_Rhop_reduced(Rhop_init_matrix,slater_lgr_init,gzproj_lgr_init)
        call get_gz_ground_state(GZ_vector)
        call print_output
        call system('cp * '//dir_iter)
        call system('rm *.out *.data fort* ')
        Jiter = Jiter + sweep_step
     end do
  case('sweepU')
     Nsweep = abs(sweep_start-sweep_stop)/abs(sweep_step)
     Uiter=sweep_start
     do i=1,Nsweep
        do iorb=1,Norb
           Uloc(iorb) = Uiter
        end do
        Ust = Uiter
        !
        call get_local_hamiltonian_trace(energy_levels)
        !
        write(dir_suffix,'(F4.2)') Uiter
        dir_iter="U"//trim(dir_suffix)
        call system('mkdir -v '//dir_iter)     
        !
        call gz_optimization_vdm_Rhop_reduced(Rhop_init_matrix,slater_lgr_init,gzproj_lgr_init)
        call get_gz_ground_state(GZ_vector)
        slater_lgr_init=0.d0
        call print_output
        call system('cp * '//dir_iter)
        call system('rm *.out *.data fort* ')
        !
        Uiter = Uiter + sweep_step
     end do
  end select
  !
CONTAINS
  !
  subroutine build_lattice_model  
    !
    integer :: ix,iy,iz,ik,Nk
    real(8),allocatable,dimension(:)   :: kx
    real(8)                            :: ts,test_k,kx_,ky_,kz_,wini,wfin,de,test_n1,test_n2
    !

    Lk=Nx
    allocate(epsik(Lk),wtk(Lk),hybik(Lk))

    wini=-Wband/2.d0
    wfin= Wband/2.d0
    epsik=linspace(wini,wfin,Lk,mesh=de)
    !
    test_k=0.d0
    test_n1=0.d0;test_n2=0.d0
    do ix=1,Lk
       wtk(ix)=4.d0/Wband/pi*sqrt(1.d0-(2.d0*epsik(ix)/Wband)**2.d0)*de
       !wtk(ix) = 1.d0/Wband*de
       if(ix==1.or.ix==Lk) wtk(ix)=0.d0
       test_n1=test_n1+wtk(ix)*fermi(epsik(ix)+Cfield*0.5d0,beta)
       test_n2=test_n2+wtk(ix)*fermi(epsik(ix)-Cfield*0.5d0,beta)
       write(77,*) epsik(ix),wtk(ix)
    end do
    hybik=0.d0
    !write(*,*) test_n1,test_n2,Cfield; stop


    ! allocate(kx(Nx))
    ! kx = linspace(0.d0,pi,Nx,.true.,.true.)
    ! Lk=Nx*Nx*Nx
    ! allocate(epsik(Lk),wtk(Lk),hybik(Lk))
    ! ik=0
    ! do ix=1,Nx
    !    do iy=1,Nx
    !       do iz=1,Nx
    !          ik=ik+1
    !          !kx_=dble(ix)/dble(Nx)*pi
    !          epsik(ik) = -2.d0/6.d0*(cos(kx(ix))+cos(kx(iy))+cos(kx(iz))) 
    !          hybik(ik) = 0.1d0/6.d0*(cos(kx(ix))-cos(kx(iy)))*cos(kx(iz)) 
    !          wtk(ik) = 1.d0/dble(Lk)
    !       end do
    !    end do
    ! end do

    call get_free_dos(epsik,wtk,file='DOS_free.kgrid')
    !stop

    allocate(Hk_tb(Ns,Ns,Lk))    
    Hk_tb=0.d0
    do ik=1,Lk
       do iorb=1,Norb
          do jorb=1,Norb
             do ispin=1,2
                istate=index(ispin,iorb)
                jstate=index(ispin,jorb)
                Hk_tb(istate,istate,ik) = epsik(ik)
                if(iorb.eq.jorb)  then
                   Hk_tb(istate,jstate,ik) = epsik(ik)
                else
                   Hk_tb(istate,jstate,ik) = hybik(ik)
                end if
             end do
          end do
       end do
    end do

    !<EXTREMA RATIO TEST
    ! e0test=0.d0
    ! do ik=1,Lk
    !    e0test = e0test + fermi_zero(epsik(ik),0.d0)*epsik(ik)*wtk(ik)
    ! end do
    !EXTREMA RATIO TEST>



  end subroutine build_lattice_model



  subroutine print_output(vdm_simplex,vdm_opt)
    real(8),dimension(Ns+1,Ns),optional :: vdm_simplex
    real(8),dimension(Nvdm_NC_opt-Nvdm_NCoff_opt),optional :: vdm_opt
    integer :: out_unit,istate,iorb,iphi,ifock,jfock
    integer,dimension(Ns) :: fock_state
    complex(8),dimension(Ns) :: tmp
    real(8) :: deltani,delta_tmp,vdm_tmp

    real(8),dimension(nFock,nFock) :: test_full_phi


    out_unit=free_unit()
    open(out_unit,file='optimized_projectors.data')
    test_full_phi=0.d0
    do iphi=1,Nphi
       !
       test_full_phi = test_full_phi + GZ_vector(iphi)*phi_basis(iphi,:,:)
       write(out_unit,*) GZ_vector(iphi)
    end do
    write(out_unit,*) '!+-----------------------------+!'
    write(out_unit,*) '!+-----------------------------+!'
    write(out_unit,*) '!+-----------------------------+!'
    do ifock=1,nFock
       do jfock=1,nFock
          write(out_unit,*) test_full_phi(ifock,jfock),ifock,jfock
       end do
    end do
    close(out_unit)    

    !
    out_unit=free_unit()
    open(out_unit,file='optimized_internal_energy.data')
    write(out_unit,'(5F18.10)') GZ_opt_energy,GZ_opt_kinetic,GZ_opt_Eloc
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_variational_density_matrix.data')
    write(out_unit,*) 'NORMAL VDM'
    do istate=1,Ns
       tmp(istate)=GZ_opt_VDM(istate,istate)
       write(out_unit,'(20F18.10)') GZ_opt_VDM(istate,:)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') tmp(1:Ns)
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_Rhop_matrix.data')
    do istate=1,Ns
       tmp(istate)=GZ_opt_Rhop(istate,istate)
       write(out_unit,'(20F18.10)') GZ_opt_Rhop(istate,1:Ns)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') tmp(1:Ns)
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_density.data')
    do istate=1,Ns
       write(out_unit,'(20F18.10)') gz_dens_matrix(istate,1:Ns)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') gz_dens(1:Ns)    
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='local_density_density.data')
    do is=1,Ns
       write(out_unit,'(20F18.10)') gz_dens_dens(is,:)
    end do
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='local_angular_momenta.data')
    write(out_unit,'(20F18.10)') gz_spin2,gz_spinZ,gz_isospin2,gz_isospinZ
    close(out_unit)
    !
    if(present(vdm_simplex)) then
       out_unit=free_unit()
       open(out_unit,file='vdm_simplex.restart')
       do jstate=1,Ns+1
          if(jstate.le.Ns) then
             do istate=1,Ns
                write(out_unit,'(20F18.10)') vdm_simplex(jstate,istate)
             end do
             if(jstate.le.Ns) write(out_unit,*)  'x'
          else
             do istate=1,Ns
                deltani=vdm_simplex(jstate,istate)-0.5
                if(deltani.gt.0.d0) then
                   delta_tmp=0.9999-vdm_simplex(jstate,istate)
                   vdm_tmp=vdm_simplex(jstate,istate)+delta_tmp*0.1
                   write(out_unit,'(20F18.10)') vdm_tmp
                else
                   delta_tmp=vdm_simplex(jstate,istate)-0.0001
                   vdm_tmp=vdm_simplex(jstate,istate)-delta_tmp*0.1
                   write(out_unit,'(20F18.10)') vdm_tmp
                end if
             end do
          end if
       end do
       close(out_unit)
    end if
    !
    if(present(vdm_opt)) then
       out_unit=free_unit()
       open(out_unit,file='vdm_seed.restart')
       do istate=1,Nvdm_NC_opt-Nvdm_NCoff_opt
          write(out_unit,'(10F18.10)')  vdm_opt(istate)
       end do
       close(out_unit)
    end if

  end subroutine print_output





  subroutine stride2reduced(x_orig,x_reduced)
    real(8),dimension(:) :: x_orig
    real(8),dimension(:) :: x_reduced
    if(size(x_orig).ne.Nopt) stop "error orig @ stride2reduced"
    if(size(x_reduced).ne.Nopt_reduced) stop "error reduced @ stride2reduced"
    !
    x_reduced(1) = x_orig(1)
    x_reduced(2) = x_orig(2)
    x_reduced(3) = x_orig(5)
    x_reduced(4) = x_orig(6)
    !
  end subroutine stride2reduced

  subroutine stride2orig(x_reduced,x_orig)
    real(8),dimension(:) :: x_orig
    real(8),dimension(:) :: x_reduced
    if(size(x_orig).ne.Nopt) stop "error orig @ stride2reduced"
    if(size(x_reduced).ne.Nopt_reduced) stop "error reduced @ stride2reduced"
    !
    x_orig=0.d0
    x_orig(1) = x_reduced(1)
    x_orig(2) = x_reduced(2)
    x_orig(5) = x_reduced(3)
    x_orig(6) = x_reduced(4)
    !
  end subroutine stride2orig












  subroutine Rhop_vec2mat(Rhop_indep,Rhop_mat)
    complex(8),dimension(:)   :: Rhop_indep
    complex(8),dimension(:,:) :: Rhop_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(Rhop_mat,1).ne.size(Rhop_mat,2)) stop "wrong stride"
    if(size(Rhop_mat,1).ne.Ns) stop "wrong stride"
    if(size(Rhop_indep).ne.NRhop_opt) stop "wrong stride!"    
    Rhop_mat = zero
    do iorb=1,Norb
       do ispin=1,2
          is=index(ispin,iorb)
          Rhop_mat(is,is) = Rhop_indep(iorb)
       end do
    end do
    !
  end subroutine Rhop_vec2mat
  subroutine Rhop_mat2vec(Rhop_mat,Rhop_indep)
    complex(8),dimension(:,:) :: Rhop_mat
    complex(8),dimension(:)   :: Rhop_indep
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    complex(8) :: test_stride
    real(8) :: test
    if(size(Rhop_mat,1).ne.size(Rhop_mat,2)) stop "wrong stride"
    if(size(Rhop_mat,1).ne.Ns) stop "wrong stride"
    if(size(Rhop_indep).ne.NRhop_opt) stop "wrong stride!"    
    !
    ispin=1;iorb=1;is=index(ispin,iorb)
    Rhop_indep(1)=Rhop_mat(is,is)
    ispin=1;iorb=2;is=index(ispin,iorb)
    Rhop_indep(2)=Rhop_mat(is,is)
    !
  end subroutine Rhop_mat2vec



  subroutine vdm_NC_vec2mat(vdm_NC_indep,vdm_NC_mat)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NC_opt) stop "wrong stride!"    
    !
    vdm_NC_mat = zero
    do iorb=1,Norb
       do ispin=1,2
          is=index(ispin,iorb)
          vdm_NC_mat(is,is) = vdm_NC_indep(iorb)
       end do
    end do
    Nopt_odiag = 0
    !
  end subroutine vdm_NC_vec2mat
  subroutine vdm_NC_mat2vec(vdm_NC_mat,vdm_NC_indep)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NC_opt) stop "wrong stride!"    
    !
    ispin=1;iorb=1;is=index(ispin,iorb)
    vdm_NC_indep(1)=vdm_NC_mat(is,is)
    ispin=1;iorb=2;is=index(ispin,iorb)
    vdm_NC_indep(2)=vdm_NC_mat(is,is)
    !
  end subroutine vdm_NC_mat2vec



  subroutine vdm_NCoff_vec2mat(vdm_NC_indep,vdm_NC_mat)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NCoff_opt) stop "wrong stride!"    
    !
    vdm_NC_mat = zero
    !
  end subroutine vdm_NCoff_vec2mat
  subroutine vdm_NCoff_mat2vec(vdm_NC_mat,vdm_NC_indep)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NCoff_opt) stop "wrong stride!"    
    !
    vdm_NC_indep = zero
    !
  end subroutine vdm_NCoff_mat2vec



end program GUTZ_mb



!AMOEBA TEST


