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
  integer                            :: ispin,jspin,iorb,i,j,istate,jstate,ifock,jorb
  integer,dimension(:),allocatable   :: fock_vec
  complex(8),dimension(:),allocatable               :: init_vec
  real(8),dimension(:),allocatable   :: variational_density_natural
  real(8),dimension(:),allocatable   :: vdm_init,vdm_out
  complex(8),dimension(:),allocatable   :: R_out
  complex(8),dimension(:),allocatable   :: Rhop_init,Rhop_out
  complex(8),dimension(:,:),allocatable   :: Rhop_init_matrix
  real(8),dimension(:,:),allocatable :: variational_density_natural_simplex
  real(8),dimension(:,:),allocatable :: variational_density_matrix
  integer                            :: out_unit,iter,iloop
  integer                            :: lattice ! 2=square;3=cubic
  real(8),dimension(:),allocatable :: epsik,hybik
  integer :: Nx,is,js,imap,jmap,Nopt
  !
  real(8) :: tmp_emin,Uiter,tmp_ene,orb_pol,tmp_real,Jh_ratio,Jiter
  real(8),dimension(:),allocatable :: dump_seed
  integer :: expected_flen,flen,unit
  logical :: seed_file


  character(len=5) :: dir_suffix
  character(len=6) :: dir_iter
  !
  complex(8),dimension(:,:,:),allocatable :: slater_lgr_init,gzproj_lgr_init
  complex(8),dimension(:,:),allocatable   :: R_init,Q_init

  complex(8),dimension(:,:,:),allocatable :: slater_lgr_init_,gzproj_lgr_init_
  complex(8),dimension(:,:),allocatable   :: R_init_,Q_init_

  complex(8),dimension(1) :: tmpQ
  character(len=6) :: sweep !sweepU,sweepJ
  real(8) :: sweep_start,sweep_stop,sweep_step
  integer ::  Nsweep,iseed
  real(8),dimension(:),allocatable :: x_reseed
  real(8) :: nread,Ntest,ndelta,nerr
  real(8) :: ndelta_
  logical :: converged_mu
  !

  !
  character(len=200) :: store_dir,read_dir
  real(8),dimension(:),allocatable :: energy_levels
  !

  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=10)
  call parse_input_variable(Wband,"Wband","inputGZ.conf",default=1.d0)
  call parse_input_variable(sweep_start,"SWEEP_START","inputGZ.conf",default=-0.4d0)  
  call parse_input_variable(sweep_stop,"SWEEP_STOP","inputGZ.conf",default=0.d0)  
  call parse_input_variable(sweep_step,"SWEEP_STEP","inputGZ.conf",default=0.05d0)  
  call parse_input_variable(read_dir,"READ_DIR","inputGZ.conf",default='~/etc_local/GZ_basis/')
  call parse_input_variable(store_dir,"STORE_DIR","inputGZ.conf",default='./READ_PHI_TRACES/')
  call parse_input_variable(nread,"NREAD","inputGZ.conf",default=0.d0)
  call parse_input_variable(nerr,"NERR","inputGZ.conf",default=1.d-7)
  call parse_input_variable(ndelta,"NDELTA","inputGZ.conf",default=1.d-4)
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")

  Norb=1
  wf_symmetry=4
  Cfield=0.d0
  !
  call initialize_local_fock_space
  !
  call init_variational_matrices(wf_symmetry,read_dir_=read_dir)  
  !  call init_variational_matrices(wf_symmetry)  
  allocate(energy_levels(Ns)); energy_levels=0.d0
  !  
  call build_lattice_model
  !
  NRhop_opt=1;   Rhop_stride_v2m => Rhop_vec2mat; Rhop_stride_m2v => Rhop_mat2vec 
  NQhop_opt=1;   Qhop_stride_v2m => Qhop_vec2mat; Qhop_stride_m2v => Qhop_mat2vec
  Nvdm_NC_opt=1; vdm_NC_stride_v2m => vdm_NC_vec2mat ; vdm_NC_stride_m2v => vdm_NC_mat2vec
  Nvdm_NCoff_opt=0; vdm_NCoff_stride_v2m => vdm_NCoff_vec2mat ; vdm_NCoff_stride_m2v => vdm_NCoff_mat2vec
  Nvdm_AC_opt=1; vdm_AC_stride_v2m => vdm_AC_vec2mat ; vdm_AC_stride_m2v => vdm_AC_mat2vec
  !
  Nopt = NRhop_opt + NQhop_opt + Nvdm_NC_opt + Nvdm_NCoff_opt + 2*Nvdm_AC_opt
  Nopt = 2*Nopt
  Nopt_reduced = 1 + 1 + 1 + 1 + 1
  !
  stride_zeros_orig2red => stride2reduced
  stride_zeros_red2orig => stride2orig
  !
  allocate(R_init(Ns,Ns),Q_init(Ns,Ns))
  allocate(slater_lgr_init(2,Ns,Ns),gzproj_lgr_init(2,Ns,Ns))
  slater_lgr_init=0.d0
  gzproj_lgr_init=0.d0
  !
  expected_flen=Nopt
  inquire(file="RQn0_root_seed.conf",exist=seed_file)
  if(seed_file) then
     flen=file_length("RQn0_root_seed.conf")
     unit=free_unit()
     open(unit,file="RQn0_root_seed.conf")
     write(*,*) 'reading root seed from file RQn0_root_seed.conf'
     if(flen.eq.expected_flen) then
        allocate(dump_seed(flen))
        !+- read from file -+!
        do i=1,flen
           read(unit,*) dump_seed(i)
        end do
        !+------------------+!        
        call dump2mats_superc(dump_seed,R_init,Q_init,slater_lgr_init,gzproj_lgr_init)        
     else
        write(*,*) 'RQn0_root_seed.conf in the wrong form',flen,expected_flen
        write(*,*) 'please check your initialization file for the root finding'
        stop
     end if
     !
  else
     !
     call init_Rhop_seed(R_init)
     call init_Qhop_seed(Q_init)
     !
     slater_lgr_init(1,:,:)=0.3d0
     slater_lgr_init(2,:,:)=0.d0
     gzproj_lgr_init(1,:,:)=0.3d0
     gzproj_lgr_init(2,:,:)=0.d0  
     !
  end if
  allocate(x_reseed(2*Nopt))



  Nsweep = abs(sweep_start-sweep_stop)/abs(sweep_step)
  Uiter=sweep_start
  do i=1,Nsweep
     !
     do iorb=1,Norb
        Uloc(iorb) = Uiter
     end do
     write(dir_suffix,'(F5.2)') Uiter
     dir_iter="U"//trim(dir_suffix)
     call system('mkdir -v '//dir_iter)     
     !
     if(nread/=0.d0) then
        converged_mu=.false.
        iloop=0
        ndelta_=ndelta
        do while(.not.converged_mu.AND.iloop<1000)
           iloop = iloop + 1
           call get_local_hamiltonian_trace
           unit=free_unit()
           open(unit,file='local_hamiltonian_parameters.out')
           write(unit,*) 'Uloc',Uloc
           write(unit,*) 'Ust',Ust
           write(unit,*) 'Jh',Jh
           write(unit,*) 'Jsf',Jsf
           write(unit,*) 'Jph',Jph
           write(unit,*) 'xmu',xmu
           close(unit)           
           call gz_optimization_vdm_Rhop_superc_reduced(R_init,Q_init,slater_lgr_init,gzproj_lgr_init)
           call get_gz_ground_state_superc(GZ_vector)  
           Ntest=sum(gz_dens(1:Ns))
           converged_mu=.true.
           call search_chemical_potential(xmu,Ntest,converged_mu,nread,ndelta,nerr)
        end do
        ndelta=ndelta_
     else
        call get_local_hamiltonian_trace
        unit=free_unit()
        open(unit,file='local_hamiltonian_parameters.out')
        write(unit,*) 'Uloc',Uloc
        write(unit,*) 'Ust',Ust
        write(unit,*) 'Jh',Jh
        write(unit,*) 'Jsf',Jsf
        write(unit,*) 'Jph',Jph
        write(unit,*) 'xmu',xmu
        close(unit)           
        call gz_optimization_vdm_Rhop_superc_reduced(R_init,Q_init,slater_lgr_init,gzproj_lgr_init)
        call get_gz_ground_state_superc(GZ_vector)  
     end if
     !
     call print_output_superc
     !
     call system('cp * '//dir_iter)
     call system('rm *.out *.data fort* ')
     Uiter = Uiter + sweep_step
  end do
  !
CONTAINS
  !
  subroutine build_lattice_model  
    !
    integer :: ix,iy,iz,ik,Nk
    real(8),allocatable,dimension(:)   :: kx
    real(8)                            :: ts,test_k,kx_,ky_,kz_,wini,wfin,de,test_n1,test_n2
    !

    ! Lk=Nx
    ! allocate(epsik(Lk),wtk(Lk),hybik(Lk))
    ! wini = -5.d0
    ! wfin =5.d0
    ! epsik=linspace(wini,wfin,Lk,mesh=de)
    ! !
    ! test_k=0.d0
    ! test_n1=0.d0;test_n2=0.d0
    ! do ix=1,Lk
    !    ! wtk(ix)=4.d0/Wband/pi*sqrt(1.d0-(2.d0*epsik(ix)/Wband)**2.d0)*de
    !    ! !wtk(ix) = 1.d0/Wband*de
    !    ! if(ix==1.or.ix==Lk) wtk(ix)=0.d0
    !    if(epsik(ix).gt.-Wband/2.d0.and.epsik(ix).lt.Wband/2) then
    !       wtk(ix)=4.d0/Wband/pi*sqrt(1.d0-(2.d0*epsik(ix)/Wband)**2.d0)*de
    !    else
    !       wtk(ix) = 0.d0
    !    end if
    !    test_n1=test_n1+wtk(ix)*fermi(epsik(ix)+Cfield*0.5d0,beta)
    !    test_n2=test_n2+wtk(ix)*fermi(epsik(ix)-Cfield*0.5d0,beta)
    !    write(77,*) epsik(ix),wtk(ix)
    ! end do
    ! hybik=0.d0
    !write(*,*) test_n1,test_n2,Cfield; stop
    !
    allocate(kx(Nx))
    kx = linspace(0.d0,pi,Nx,.true.,.true.)
    Lk=Nx*Nx*Nx
    allocate(epsik(Lk),wtk(Lk),hybik(Lk))
    ik=0
    do ix=1,Nx
       do iy=1,Nx
          do iz=1,Nx
             ik=ik+1
             epsik(ik) = -2.d0/6.d0*(cos(kx(ix))+cos(kx(iy))+cos(kx(iz))) 
             hybik(ik) =  0.d0/6.d0*(cos(kx(ix))-cos(kx(iy)))*cos(kx(iz)) 
             wtk(ik) = 1.d0/dble(Lk)
          end do
       end do
    end do
    !
    call get_free_dos(epsik,wtk,file='DOS_free.kgrid')
    !
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





  subroutine print_output(vdm_simplex)
    real(8),dimension(Ns+1,Ns),optional :: vdm_simplex
    integer :: out_unit,istate,iorb,iphi,ifock,jfock
    integer,dimension(Ns) :: fock_state
    real(8),dimension(Ns) :: tmp
    real(8) :: deltani,delta_tmp,vdm_tmp

    real(8),dimension(nFock,nFock) :: test_full_phi

    out_unit=free_unit()
    open(out_unit,file='optimized_projectors.data')

    !+- CHANGE THE NAME OF GZ_opt_projector_diag -+!
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
    !out_unit=free_unit()
    open(out_unit,file='optimized_variational_density_matrix.data')
    !write(out_unit,'(20F18.10)') variational_density_natural(1:Ns)
    do istate=1,Ns
       tmp(istate)=GZ_opt_VDM(istate,istate)
       write(out_unit,'(20F18.10)') GZ_opt_VDM(istate,:)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') tmp(1:Ns)
    close(out_unit)
    !
    !out_unit=free_unit()
    open(out_unit,file='optimized_Rhop_matrix.data')
    do istate=1,Ns
       tmp(istate)=GZ_opt_Rhop(istate,istate)
       write(out_unit,'(20F18.10)') GZ_opt_Rhop(istate,1:Ns)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') tmp(1:Ns)
    close(out_unit)
    !
    !out_unit=free_unit()
    open(out_unit,file='optimized_density.data')
    do istate=1,Ns
       write(out_unit,'(20F18.10)') gz_dens_matrix(istate,1:Ns)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') gz_dens(1:Ns)    
    close(out_unit)
    !
    !out_unit=free_unit()
    open(out_unit,file='local_density_density.data')
    do is=1,Ns
       write(out_unit,'(20F18.10)') gz_dens_dens(is,:)
    end do
    close(out_unit)
    !
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
  end subroutine print_output






  subroutine print_output_superc(vdm_simplex)
    real(8),dimension(Ns+1,Ns),optional :: vdm_simplex
    integer :: out_unit,istate,iorb,iphi,ifock,jfock,is,js,ik
    integer,dimension(Ns) :: fock_state
    complex(8),dimension(Ns) :: tmp
    real(8) :: deltani,delta_tmp,vdm_tmp

    real(8),dimension(nFock,nFock) :: test_full_phi

    out_unit=free_unit()
    open(out_unit,file='optimized_projectors.data')

    !+- CHANGE THE NAME OF GZ_opt_projector_diag -+!
    test_full_phi=0.d0
    do iphi=1,Nphi
       !
       test_full_phi = test_full_phi + GZ_vector(iphi)*phi_basis(iphi,:,:)
       write(out_unit,'(2F18.10)') GZ_vector(iphi)
    end do
    close(out_unit)    
    open(out_unit,file='optimized_phi_matrix.data')    
    do ifock=1,nFock
       do jfock=1,nFock
          write(out_unit,*) test_full_phi(ifock,jfock),ifock,jfock
       end do
    end do
    close(out_unit)
    !
    !+- STORE SLATER DETERMINANT -+!
    do is=1,Ns
       do js=1,Ns
          out_unit=free_unit()
          open(out_unit,file='optimized_slater_normal_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data')
          do ik=1,Lk
             write(out_unit,'(2F18.10)') dreal(GZ_opt_slater_superc(1,is,js,ik)),dimag(GZ_opt_slater_superc(1,is,js,ik))
          end do
          close(out_unit)
       end do
    end do
    do is=1,Ns
       do js=1,Ns
          out_unit=free_unit()
          open(out_unit,file='optimized_slater_anomalous_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data')
          do ik=1,Lk
             write(out_unit,'(2F18.10)') dreal(GZ_opt_slater_superc(2,is,js,ik)),dimag(GZ_opt_slater_superc(2,is,js,ik))
          end do
          close(out_unit)
       end do
    end do
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
       tmp(istate)=GZ_opt_VDM_superc(1,istate,istate)
       write(out_unit,'(20F18.10)') dreal(GZ_opt_VDM_superc(1,istate,:)),dimag(GZ_opt_VDM_superc(1,istate,:))
    end do
    write(out_unit,*)
    write(out_unit,*) 'ANOMALOUS VDM'
    do istate=1,Ns
       write(out_unit,'(20F18.10)') dreal(GZ_opt_VDM_superc(2,istate,:)),dimag(GZ_opt_VDM_superc(2,istate,:))
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
    open(out_unit,file='optimized_Qhop_matrix.data')
    do istate=1,Ns
       write(out_unit,'(20F18.10)') GZ_opt_Qhop(istate,1:Ns)
    end do
    ! write(out_unit,*) ! on the last line store the diagonal elements
    ! write(out_unit,'(20F18.10)') tmp(1:Ns)
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
    !out_unit=free_unit()
    open(out_unit,file='local_sc_order.data')
    do is=1,Ns
       write(out_unit,'(20F18.10)') dreal(gz_sc_order(is,:)),dimag(gz_sc_order(is,:))
    end do
    close(out_unit)


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
  end subroutine print_output_superc




  subroutine stride2reduced(x_orig,x_reduced)
    real(8),dimension(:) :: x_orig
    real(8),dimension(:) :: x_reduced
    if(size(x_orig).ne.Nopt) stop "error orig @ stride2reduced"
    if(size(x_reduced).ne.Nopt_reduced) stop "error reduced @ stride2reduced"
    !+- pecionata -+!
    x_reduced(1) = x_orig(1)
    x_reduced(2) = x_orig(3)
    x_reduced(3) = x_orig(5)
    x_reduced(4) = x_orig(7)
    x_reduced(5) = x_orig(9)
  end subroutine stride2reduced

  subroutine stride2orig(x_reduced,x_orig)
    real(8),dimension(:) :: x_orig
    real(8),dimension(:) :: x_reduced
    if(size(x_orig).ne.Nopt) stop "error orig @ stride2reduced"
    if(size(x_reduced).ne.Nopt_reduced) stop "error reduced @ stride2reduced"
    !+- pecionata -+!
    x_orig=0.d0
    x_orig(1) = x_reduced(1)
    x_orig(3) = x_reduced(2)
    x_orig(5) = x_reduced(3)
    x_orig(7) = x_reduced(4)
    x_orig(9) = x_reduced(5)
  end subroutine stride2orig









  !+- STRIDES -+!
  ! subroutine Rhop_vec2mat(Rhop_indep,Rhop_mat)
  !   complex(8),dimension(:)   :: Rhop_indep
  !   complex(8),dimension(:,:) :: Rhop_mat
  !   integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
  !   if(size(Rhop_mat,1).ne.size(Rhop_mat,2)) stop "wrong stride"
  !   if(size(Rhop_mat,1).ne.Ns) stop "wrong stride"
  !   if(size(Rhop_indep).ne.NRhop_opt) stop "wrong stride!"    
  !   Rhop_mat = zero
  !   do iorb=1,Norb
  !      do ispin=1,2
  !         is=index(ispin,iorb)
  !         Rhop_mat(is,is) = Rhop_indep(iorb)
  !      end do
  !   end do
  !   !
  ! end subroutine Rhop_vec2mat
  ! subroutine Rhop_mat2vec(Rhop_mat,Rhop_indep)
  !   complex(8),dimension(:,:) :: Rhop_mat
  !   complex(8),dimension(:)   :: Rhop_indep
  !   integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
  !   complex(8) :: test_stride
  !   real(8) :: test
  !   if(size(Rhop_mat,1).ne.size(Rhop_mat,2)) stop "wrong stride"
  !   if(size(Rhop_mat,1).ne.Ns) stop "wrong stride"
  !   if(size(Rhop_indep).ne.NRhop_opt) stop "wrong stride!"    
  !   !
  !   do iorb=1,Norb
  !      ispin=1
  !      is=index(ispin,iorb)
  !      Rhop_indep(iorb)=Rhop_mat(is,is)
  !   end do
  !   !
  ! end subroutine Rhop_mat2vec

  ! subroutine Qhop_vec2mat(Qhop_indep,Qhop_mat)
  !   complex(8),dimension(:)   :: Qhop_indep
  !   complex(8),dimension(:,:) :: Qhop_mat
  !   integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
  !   if(size(Qhop_mat,1).ne.size(Qhop_mat,2)) stop "wrong stride"
  !   if(size(Qhop_mat,1).ne.Ns) stop "wrong stride"
  !   if(size(Qhop_indep).ne.NQhop_opt) stop "wrong stride!"    

  !   !+- to understand why the SU2xisoZ allows only inter-orbital SC ???? -+!
  !   Qhop_mat = zero
  !   do iorb=1,Norb
  !      do jorb=1,Norb
  !         do ispin=1,2
  !            jspin=3-ispin
  !            is=index(ispin,iorb)
  !            js=index(jspin,jorb)
  !            !             write(*,*) is,js,size(Qhop_mat,1),size(Qhop_mat,2)
  !            if(iorb.ne.jorb) then
  !               Qhop_mat(is,js) = (-1.d0)**dble(jspin)*Qhop_indep(1)
  !            else
  !               Qhop_mat(is,js) = zero
  !            end if
  !         end do
  !      end do
  !   end do
  ! end subroutine Qhop_vec2mat
  ! subroutine Qhop_mat2vec(Qhop_mat,Qhop_indep)
  !   complex(8),dimension(:)   :: Qhop_indep
  !   complex(8),dimension(:,:) :: Qhop_mat
  !   integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
  !   if(size(Qhop_mat,1).ne.size(Qhop_mat,2)) stop "wrong stride"
  !   if(size(Qhop_mat,1).ne.Ns) stop "wrong stride"
  !   if(size(Qhop_indep).ne.NQhop_opt) stop "wrong stride!"    
  !   !
  !   iorb=1;jorb=2;ispin=1;jspin=2
  !   is=index(ispin,iorb)
  !   js=index(jspin,jorb)
  !   Qhop_indep(1) = Qhop_mat(is,js)
  ! end subroutine Qhop_mat2vec



  ! subroutine vdm_NC_vec2mat(vdm_NC_indep,vdm_NC_mat)
  !   complex(8),dimension(:)   :: vdm_NC_indep
  !   complex(8),dimension(:,:) :: vdm_NC_mat
  !   integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
  !   if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
  !   if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
  !   if(size(vdm_NC_indep).ne.Nvdm_NC_opt) stop "wrong stride!"    
  !   !
  !   vdm_NC_mat = zero
  !   do iorb=1,Norb
  !      do ispin=1,2
  !         is=index(ispin,iorb)
  !         vdm_NC_mat(is,is) = vdm_NC_indep(iorb)
  !      end do
  !   end do
  !   Nopt_odiag = 0
  !   !
  ! end subroutine vdm_NC_vec2mat
  ! subroutine vdm_NC_mat2vec(vdm_NC_mat,vdm_NC_indep)
  !   complex(8),dimension(:)   :: vdm_NC_indep
  !   complex(8),dimension(:,:) :: vdm_NC_mat
  !   integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
  !   if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
  !   if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
  !   if(size(vdm_NC_indep).ne.Nvdm_NC_opt) stop "wrong stride!"    
  !   !
  !   ispin=1
  !   do iorb=1,Norb
  !      is=index(ispin,iorb)       
  !      vdm_NC_indep(iorb) = vdm_NC_mat(is,is)
  !   end do
  !   !
  ! end subroutine vdm_NC_mat2vec



  ! subroutine vdm_NCoff_vec2mat(vdm_NC_indep,vdm_NC_mat)
  !   complex(8),dimension(:)   :: vdm_NC_indep
  !   complex(8),dimension(:,:) :: vdm_NC_mat
  !   integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
  !   if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
  !   if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
  !   if(size(vdm_NC_indep).ne.Nvdm_NCoff_opt) stop "wrong stride!"    
  !   !
  !   vdm_NC_mat = zero
  !   !
  ! end subroutine vdm_NCoff_vec2mat
  ! subroutine vdm_NCoff_mat2vec(vdm_NC_mat,vdm_NC_indep)
  !   complex(8),dimension(:)   :: vdm_NC_indep
  !   complex(8),dimension(:,:) :: vdm_NC_mat
  !   integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
  !   if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
  !   if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
  !   if(size(vdm_NC_indep).ne.Nvdm_NCoff_opt) stop "wrong stride!"    
  !   !
  !   vdm_NC_indep = zero
  !   !
  ! end subroutine vdm_NCoff_mat2vec





  ! !
  ! subroutine vdm_AC_vec2mat(vdm_AC_indep,vdm_AC_mat)
  !   complex(8),dimension(:)   :: vdm_AC_indep
  !   complex(8),dimension(:,:) :: vdm_AC_mat
  !   integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
  !   if(size(vdm_AC_mat,1).ne.size(vdm_AC_mat,2)) stop "wrong stride"
  !   if(size(vdm_AC_mat,1).ne.Ns) stop "wrong stride"
  !   if(size(vdm_AC_indep).ne.Nvdm_AC_opt) stop "wrong stride!"    
  !   !
  !   vdm_AC_mat = zero
  !   do iorb=1,Norb
  !      do jorb=1,Norb
  !         do ispin=1,2
  !            jspin=3-ispin
  !            is=index(ispin,iorb)
  !            js=index(jspin,jorb)
  !            if(iorb.ne.jorb) then
  !               vdm_AC_mat(is,js) = (-1.d0)**dble(jspin)*vdm_AC_indep(1)
  !            else
  !               vdm_AC_mat(is,js) = zero
  !            end if
  !         end do
  !      end do
  !   end do
  !   !
  ! end subroutine vdm_AC_vec2mat
  ! subroutine vdm_AC_mat2vec(vdm_AC_mat,vdm_AC_indep)
  !   complex(8),dimension(:)   :: vdm_AC_indep
  !   complex(8),dimension(:,:) :: vdm_AC_mat
  !   integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
  !   if(size(vdm_AC_mat,1).ne.size(vdm_AC_mat,2)) stop "wrong stride"
  !   if(size(vdm_AC_mat,1).ne.Ns) stop "wrong stride"
  !   if(size(vdm_AC_indep).ne.Nvdm_AC_opt) stop "wrong stride!"    
  !   !
  !   iorb=1;jorb=2;ispin=1;jspin=2
  !   is=index(ispin,iorb)
  !   js=index(jspin,jorb)
  !   vdm_AC_indep(1) = vdm_AC_mat(is,js)
  !   !
  ! end subroutine vdm_AC_mat2vec













  subroutine Rhop_vec2mat(Rhop_indep,Rhop_mat)
    complex(8),dimension(:)   :: Rhop_indep
    complex(8),dimension(:,:) :: Rhop_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(Rhop_mat,1).ne.size(Rhop_mat,2)) stop "wrong stride"
    if(size(Rhop_mat,1).ne.Ns) stop "wrong stride"
    if(size(Rhop_indep).ne.NRhop_opt) stop "wrong stride!"    
    Rhop_mat = zero
    do is=1,Ns
       Rhop_mat(is,is) = Rhop_indep(1)
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
    Rhop_indep(1)=Rhop_mat(1,1)
    !
  end subroutine Rhop_mat2vec

  subroutine Qhop_vec2mat(Qhop_indep,Qhop_mat)
    complex(8),dimension(:)   :: Qhop_indep
    complex(8),dimension(:,:) :: Qhop_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(Qhop_mat,1).ne.size(Qhop_mat,2)) stop "wrong stride"
    if(size(Qhop_mat,1).ne.Ns) stop "wrong stride"
    if(size(Qhop_indep).ne.NQhop_opt) stop "wrong stride!"    
    Qhop_mat = zero
    do iorb=1,Norb
       do jorb=1,Norb
          do ispin=1,2
             jspin=3-ispin
             is=index(ispin,iorb)
             js=index(jspin,jorb)
             if(iorb.eq.jorb) then
                Qhop_mat(is,js) = (-1.d0)**dble(jspin)*Qhop_indep(1)
             else
                Qhop_mat(is,js) = zero
             end if
          end do
       end do
    end do
  end subroutine Qhop_vec2mat
  subroutine Qhop_mat2vec(Qhop_mat,Qhop_indep)
    complex(8),dimension(:)   :: Qhop_indep
    complex(8),dimension(:,:) :: Qhop_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(Qhop_mat,1).ne.size(Qhop_mat,2)) stop "wrong stride"
    if(size(Qhop_mat,1).ne.Ns) stop "wrong stride"
    if(size(Qhop_indep).ne.NQhop_opt) stop "wrong stride!"    
    !
    iorb=1;jorb=1;ispin=1;jspin=2
    is=index(ispin,iorb)
    js=index(jspin,jorb)
    Qhop_indep(1) = Qhop_mat(is,js)
  end subroutine Qhop_mat2vec



  subroutine vdm_NC_vec2mat(vdm_NC_indep,vdm_NC_mat)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NC_opt) stop "wrong stride!"    
    !
    vdm_NC_mat = zero
    do is=1,Ns       
       vdm_NC_mat(is,is) = vdm_NC_indep(1)
    end do
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
    vdm_NC_indep(1) = vdm_NC_mat(1,1)
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





  !
  subroutine vdm_AC_vec2mat(vdm_AC_indep,vdm_AC_mat)
    complex(8),dimension(:)   :: vdm_AC_indep
    complex(8),dimension(:,:) :: vdm_AC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_AC_mat,1).ne.size(vdm_AC_mat,2)) stop "wrong stride"
    if(size(vdm_AC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_AC_indep).ne.Nvdm_AC_opt) stop "wrong stride!"    
    !
    vdm_AC_mat = zero
    do iorb=1,Norb
       do jorb=1,Norb
          do ispin=1,2
             jspin=3-ispin
             is=index(ispin,iorb)
             js=index(jspin,jorb)
             if(iorb.eq.jorb) then
                vdm_AC_mat(is,js) = (-1.d0)**dble(jspin)*vdm_AC_indep(1)
             else
                vdm_AC_mat(is,js) = zero
             end if
          end do
       end do
    end do
    !
  end subroutine vdm_AC_vec2mat
  subroutine vdm_AC_mat2vec(vdm_AC_mat,vdm_AC_indep)
    complex(8),dimension(:)   :: vdm_AC_indep
    complex(8),dimension(:,:) :: vdm_AC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_AC_mat,1).ne.size(vdm_AC_mat,2)) stop "wrong stride"
    if(size(vdm_AC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_AC_indep).ne.Nvdm_AC_opt) stop "wrong stride!"    
    !
    iorb=1;jorb=1;ispin=1;jspin=2
    is=index(ispin,iorb)
    js=index(jspin,jorb)
    vdm_AC_indep(1) = vdm_AC_mat(is,js)
    !
  end subroutine vdm_AC_mat2vec







end program GUTZ_mb



!AMOEBA TEST


