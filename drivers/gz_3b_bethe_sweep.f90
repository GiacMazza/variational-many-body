program GUTZ_mb
  USE SCIFOR
  !
  USE DMFT_MISC
  USE SF_PARSE_INPUT
  !
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
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
  complex(8),dimension(:),allocatable   :: R_out,vdm_tmp
  complex(8),dimension(:),allocatable   :: Rhop_init,Rhop_out
  complex(8),dimension(:,:),allocatable   :: Rhop_init_matrix,tmp_vdm_mat
  real(8),dimension(:,:),allocatable :: variational_density_natural_simplex
  real(8),dimension(:,:),allocatable :: variational_density_matrix
  integer                            :: out_unit,iter
  integer                            :: lattice ! 2=square;3=cubic
  real(8),dimension(:),allocatable :: epsik,hybik
  integer :: Nx,is,js,imap,jmap,Nopt
  !
  real(8) :: tmp_emin,Uiter,tmp_ene,orb_pol,tmp_real,Jh_ratio,Jiter,VDMiter
  real(8),dimension(:),allocatable :: dump_seed
  integer :: expected_flen,flen,unit
  logical :: seed_file


  character(len=5) :: dir_suffix
  character(len=8) :: dir_iter
  !
  complex(8),dimension(:,:,:),allocatable :: slater_lgr_init,gzproj_lgr_init
  complex(8),dimension(:,:),allocatable   :: R_init,Q_init

  complex(8),dimension(:,:,:),allocatable :: slater_lgr_init_,gzproj_lgr_init_
  complex(8),dimension(:,:),allocatable   :: R_init_,Q_init_
  complex(8),dimension(1) :: tmpQ

  character(len=11) :: task !e_sweep_vdm,min_simplex,min_galahad,nRQfree_min
  character(len=6) :: sweep !sweepU,sweepJ
  real(8) :: sweep_start,sweep_stop,sweep_step
  integer ::  Nsweep,iseed
  


  !
  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=1000)
  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"Wband","inputGZ.conf",default=1.d0)
  call parse_input_variable(lattice,"LAT_DIMENSION","inputGZ.conf",default=3)  
  call parse_input_variable(task,"TASK","inputGZ.conf",default="min_galahad")  
  call parse_input_variable(sweep,"SWEEP","inputGZ.conf",default='sweepJ')  
  call parse_input_variable(sweep_start,"SWEEP_START","inputGZ.conf",default=-0.4d0)  
  call parse_input_variable(sweep_stop,"SWEEP_STOP","inputGZ.conf",default=0.d0)  
  call parse_input_variable(sweep_step,"SWEEP_STEP","inputGZ.conf",default=0.05d0)  
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")
  !
  if(Norb.eq.1.and.wf_symmetry.eq.1) then
     write(*,*) 'WARNING THE O(1) x SU(2)c x ORBITAL_ROTATION = O(1) x SU(2)c for the Norb=1 case!'
     wf_symmetry=0
  end if
  !
  !NOTE: ON HUNDS COUPLINGS:
  !NORB=3 RATATIONAL INVARIANT HAMILTONIAN       :: Jsf=Jh, Jph=U-Ust-J   (NO relation between Ust and U)
  !       FULLY ROTATIONAL INVARIANT HAMILTONIAN :: Jsf=Jh, Jph=J, Ust = U - 2J   
  ! Jh = Jh*Uloc(1)
  
  !
  Jsf = Jh
  Jph = Jh
  Ust = Uloc(1)-2.d0*Jh
  !
  call initialize_local_fock_space
  !
  call init_variational_matrices
  !
  call build_lattice_model
  !
  NRhop_opt=1;   Rhop_stride_v2m => Rhop_vec2mat; Rhop_stride_m2v => Rhop_mat2vec 
  !  NQhop_opt=1;   Qhop_stride_v2m => Qhop_vec2mat; Qhop_stride_m2v => Qhop_mat2vec
  Nvdm_NC_opt=1; vdm_NC_stride_v2m => vdm_NC_vec2mat ; vdm_NC_stride_m2v => vdm_NC_mat2vec
  Nvdm_NCoff_opt=0; vdm_NCoff_stride_v2m => vdm_NCoff_vec2mat ; vdm_NCoff_stride_m2v => vdm_NCoff_mat2vec
  !  Nvdm_AC_opt=1; vdm_AC_stride_v2m => vdm_AC_vec2mat ; vdm_AC_stride_m2v => vdm_AC_mat2vec
  Nopt = NRhop_opt  + Nvdm_NC_opt + Nvdm_NCoff_opt 
  !
  select case(task)
  case("min_galahad")
     !+- nlsq_minimization using the galahad routine -+!     
     allocate(vdm_init(Nvdm_NC_opt-Nvdm_NCoff_opt),vdm_out(Nvdm_NC_opt-Nvdm_NCoff_opt))
     !     
     call initialize_variational_density(vdm_init)
     !
     select case(sweep)
     case('sweepJ')
        !+- sweep JHund -+!
        Nsweep = abs(sweep_start-sweep_stop)/abs(sweep_step)
        Jiter = sweep_start
        do i=1,Nsweep
           !
           Jh = Jiter
           Jsf = Jh
           Jph = Jh
           Ust = Uloc(1)-2.d0*Jh
           call build_local_hamiltonian     
           phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
           phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)
           !
           write(dir_suffix,'(F4.2)') abs(Jiter)
           dir_iter="J"//trim(dir_suffix)
           call system('mkdir -v '//dir_iter)                
           !
           call gz_optimization_vdm_nlsq(vdm_init,vdm_out)  
           !
           call get_gz_ground_state(GZ_vector)
           !
           call print_output
           !
           call system('cp * '//dir_iter)
           call system('rm *.out *.data fort.* ')
           !
           vdm_init=vdm_out
           Jiter = Jiter + sweep_step
           !
        end do
     case('sweepU')
        Nsweep = abs(sweep_start-sweep_stop)/abs(sweep_step)        
        Uiter=sweep_start
        do i=1,Nsweep
           do iorb=1,Norb
              Uloc(iorb) = Uiter
           end do
           !
           Jsf = Jh
           Jph = Jh
           Ust = Uloc(1)-2.d0*Jh
           !
           call build_local_hamiltonian     
           phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
           phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)
           !
           write(dir_suffix,'(F4.2)') abs(Uiter)
           dir_iter="U"//trim(dir_suffix)
           call system('mkdir -v '//dir_iter)                
           !
           call gz_optimization_vdm_nlsq(vdm_init,vdm_out)  
           !
           call get_gz_ground_state(GZ_vector)
           !
           call print_output
           !
           call system('cp * '//dir_iter)
           call system('rm *.out *.data fort.* ')
           !
           vdm_init=vdm_out
           Uiter = Uiter + sweep_step
           !
        end do
     end select
           
  case("min_simplex")
     !+- simplex_minimization using the amoeab routine +-!
     !     
     allocate(variational_density_natural_simplex(Ns+1,Ns)); allocate(variational_density_natural(Ns))     
     call initialize_variational_density_simplex(variational_density_natural_simplex)

     select case(sweep)
     case('sweepJ')
        !+- sweep JHund -+!
        Nsweep = abs(sweep_start-sweep_stop)/abs(sweep_step)
        Jiter = sweep_start
        do i=1,Nsweep
           !
           Jh = Jiter
           Jsf = Jh
           Jph = Jh
           Ust = Uloc(1)-2.d0*Jh
           call build_local_hamiltonian     
           phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
           phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)
           !
           write(dir_suffix,'(F4.2)') abs(Jiter)
           dir_iter="J"//trim(dir_suffix)
           call system('mkdir -v '//dir_iter)                
           !
           call gz_optimization_vdm_simplex(variational_density_natural_simplex,variational_density_natural)  
           !
           call get_gz_ground_state(GZ_vector)
           !
           call print_output           
           !
           call system('cp * '//dir_iter)
           call system('rm *.out *.data fort.* ')
           !
           vdm_init=vdm_out
           Jiter = Jiter + sweep_step
           !
        end do
     case('sweepU')
        Nsweep = abs(sweep_start-sweep_stop)/abs(sweep_step)        
        Uiter=sweep_start
        do i=1,Nsweep
           do iorb=1,Norb
              Uloc(iorb) = Uiter
           end do
           !
           Jsf = Jh
           Jph = Jh
           Ust = Uloc(1)-2.d0*Jh
           !
           call build_local_hamiltonian     
           phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
           phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)
           !
           write(dir_suffix,'(F4.2)') abs(Uiter)
           dir_iter="U"//trim(dir_suffix)
           call system('mkdir -v '//dir_iter)                
           !
           call gz_optimization_vdm_simplex(variational_density_natural_simplex,variational_density_natural)  
           !
           call get_gz_ground_state(GZ_vector)
           !
           call print_output           
           !
           call system('cp * '//dir_iter)
           call system('rm *.out *.data fort.* ')
           !
           vdm_init=vdm_out
           Uiter = Uiter + sweep_step
           !
        end do
     end select
     !
  case("e_sweep_vdm")
     !
     !+- energy computation for a fixed Variational Density Matrix provided by the user (read independent vdm elements from file) -+!
     !
     allocate(vdm_tmp(Nvdm_NC_opt-Nvdm_NCoff_opt),vdm_init(Nvdm_NC_opt-Nvdm_NCoff_opt),tmp_vdm_mat(Ns,Ns))
     !
     call initialize_variational_density(vdm_init); vdm_tmp=vdm_init; deallocate(vdm_init)
     call vdm_NC_stride_v2m(vdm_tmp,tmp_vdm_mat); allocate(vdm_init(Ns))
     !
     do is=1,Ns
        vdm_init(is)=dreal(tmp_vdm_mat(is,is))
     end do
     !
     opt_energy_unit=free_unit()
     open(opt_energy_unit,file='GZ_OptEnergy_VS_vdm.out')
     opt_rhop_unit=free_unit()
     if(gz_superc) then
        open(opt_rhop_unit,file='GZ_OptRhop_VS_vdm.out')
        opt_qhop_unit=free_unit()
     end if
     open(opt_qhop_unit,file='GZ_OptQhop_VS_vdm.out')
     opt_GZ_unit=free_unit()
     open(opt_GZ_unit,file='GZ_OptProj_VS_vdm.out')
     if(GZmin_verbose) then
        GZmin_unit=free_unit()
        open(GZmin_unit,file='GZ_SelfCons_min_verbose.out')
        GZmin_unit_=free_unit()
        open(GZmin_unit_,file='GZ_proj_min.out')
     end if
     optimization_flag=.true.
     if(.not.allocated(GZ_vector)) allocate(GZ_vector(Nphi))

     Nsweep = abs(sweep_start-sweep_stop)/abs(sweep_step)
     vdm_init= sweep_start
     do i=1,Nsweep
        !
        write(dir_suffix,'(F4.2)') vdm_init(1)
        dir_iter="n"//trim(dir_suffix)
        call system('mkdir -v '//dir_iter)                
        tmp_emin = gz_energy_vdm(vdm_init);
        call get_gz_ground_state(GZ_vector)  
        call print_output        
        call system('cp * '//dir_iter)
        vdm_init = vdm_init + sweep_step
        !
     end do
     !
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
    !
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
    !          hybik(ik) = 0.d0/6.d0*(cos(kx(ix))-cos(kx(iy)))*cos(kx(iz)) 
    !          wtk(ik) = 1.d0/dble(Lk)
    !       end do
    !    end do
    ! end do
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
  end subroutine build_lattice_model





  subroutine print_output(vdm_simplex)
    real(8),dimension(Ns+1,Ns),optional :: vdm_simplex
    integer :: out_unit,istate,iorb,iphi,ifock,jfock
    integer,dimension(Ns) :: fock_state
    real(8),dimension(Ns) :: tmp
    real(8) :: deltani,delta_tmp,vdm_tmp

    complex(8),dimension(nFock,nFock) :: test_full_phi

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
    open(out_unit,file='orbital_double_occupancy.data')
    write(out_unit,'(20F18.10)') gz_docc(1:Norb)
    close(out_unit)
    !
    !out_unit=free_unit()
    open(out_unit,file='orbital_density_density.data')
    do iorb=1,Norb
       write(out_unit,'(20F18.10)') gz_dens_dens_orb(iorb,:)
    end do
    close(out_unit)
    !
    !out_unit=free_unit()
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






  subroutine print_output_superc(vdm_simplex,vdm_opt)
    real(8),dimension(Ns+1,Ns),optional :: vdm_simplex
    real(8),dimension(Nvdm_NC_opt-Nvdm_NCoff_opt),optional :: vdm_opt
    integer :: out_unit,istate,iorb,iphi,ifock,jfock
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
    open(out_unit,file='orbital_double_occupancy.data')
    write(out_unit,'(20F18.10)') gz_docc(1:Norb)
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='orbital_density_density.data')
    do iorb=1,Norb
       write(out_unit,'(20F18.10)') gz_dens_dens_orb(iorb,:)
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
    if(present(vdm_opt)) then
       out_unit=free_unit()
       open(out_unit,file='vdm_seed.restart')
       do istate=1,Nvdm_NC_opt-Nvdm_NCoff_opt
          write(out_unit,'(10F18.10)')  vdm_opt(istate)
       end do
       close(out_unit)
    end if

  end subroutine print_output_superc















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


