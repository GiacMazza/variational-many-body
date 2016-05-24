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
  USE GZ_ENERGY_MINIMIZATION
  USE GZ_OPTIMIZED_ENERGY
  !

  !
  USE GZ_MATRIX_BASIS
  !
  implicit none
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

  character(len=200) :: store_dir,read_dir
  real(8),dimension(:),allocatable :: energy_levels


  !
  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=1000)
  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"Wband","inputGZ.conf",default=1.d0)
  call parse_input_variable(task,"TASK","inputGZ.conf",default="min_galahad")  
  call parse_input_variable(read_dir,"READ_DIR","inputGZ.conf",default='~/etc_local/GZ_basis/')
  call parse_input_variable(store_dir,"STORE_DIR","inputGZ.conf",default='./READ_PHI_TRACES/')
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")

  if(Norb.eq.1.and.wf_symmetry.eq.1) then
     write(*,*) 'WARNING THE O(1) x SU(2)c x ORBITAL_ROTATION = O(1) x SU(2)c for the Norb=1 case!'
     wf_symmetry=0
  end if
  !
  !NOTE: ON HUNDS COUPLINGS:
  !NORB=3 RATATIONAL INVARIANT HAMILTONIAN       :: Jsf=Jh, Jph=U-Ust-J   (NO relation between Ust and U)
  !FULLY ROTATIONAL INVARIANT HAMILTONIAN        :: Jsf=Jh, Jph=J, Ust = U - 2J   
  !
  call initialize_local_fock_space
  !
  call init_variational_matrices(wf_symmetry,store_dir_=store_dir,read_dir_=read_dir)  
  !
  allocate(energy_levels(Ns))
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        energy_levels(istate) = Cfield*0.5d0
        if(iorb.eq.2) energy_levels(istate) = -Cfield*0.5d0
     end do
  end do
  call get_local_hamiltonian_trace(energy_levels)
  !
  call build_lattice_model
  !
  NRhop_opt=1;   Rhop_stride_v2m => Rhop_vec2mat; Rhop_stride_m2v => Rhop_mat2vec 
  Nvdm_NC_opt=1; vdm_NC_stride_v2m => vdm_NC_vec2mat ; vdm_NC_stride_m2v => vdm_NC_mat2vec
  Nvdm_NCoff_opt=0; vdm_NCoff_stride_v2m => vdm_NCoff_vec2mat ; vdm_NCoff_stride_m2v => vdm_NCoff_mat2vec
  !Nopt = NRhop_opt + NQhop_opt + Nvdm_NC_opt + Nvdm_NCoff_opt + 2*Nvdm_AC_opt
  !

  select case(task)
  case("min_galahad")
     !+- nlsq_minimization using the galahad routine -+!     
     allocate(vdm_init(Nvdm_NC_opt-Nvdm_NCoff_opt),vdm_out(Nvdm_NC_opt-Nvdm_NCoff_opt))
     !
     call initialize_variational_density(vdm_init)
     !
     call gz_optimization_vdm_nlsq(vdm_init,vdm_out)  
     !
     call get_gz_ground_state(GZ_vector)
     !
     call print_output(vdm_opt=vdm_out)
     !
  case("min_simplex")
     !+- simplex_minimization using the amoeab routine +-!
     !
     allocate(variational_density_natural_simplex(Ns+1,Ns)); allocate(variational_density_natural(Ns))     
     !
     call initialize_variational_density_simplex(variational_density_natural_simplex)
     !
     call gz_optimization_vdm_simplex(variational_density_natural_simplex,variational_density_natural)  
     !
     call get_gz_ground_state(GZ_vector)
     !
     call print_output(variational_density_natural_simplex)
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
     open(opt_rhop_unit,file='GZ_OptRhop_VS_vdm.out')
     opt_qhop_unit=free_unit()
     open(opt_GZ_unit,file='GZ_OptProj_VS_vdm.out')
     if(GZmin_verbose) then
        GZmin_unit=free_unit()
        open(GZmin_unit,file='GZ_SelfCons_min_verbose.out')
        GZmin_unit_=free_unit()
        open(GZmin_unit_,file='GZ_proj_min.out')
     end if
     optimization_flag=.true.
     if(.not.allocated(GZ_vector)) allocate(GZ_vector(Nphi))
     tmp_emin = gz_energy_vdm(vdm_init);
     call get_gz_ground_state(GZ_vector)    
     call print_output
     !
  end select


  !
CONTAINS
  !
  subroutine build_lattice_model  
    !
    integer :: ix,iy,iz,ik,Nk
    real(8),allocatable,dimension(:)   :: kx
    real(8)                            :: ts,test_k,kx_,ky_,kz_,wini,wfin,de
    !


    Lk=Nx
    allocate(epsik(Lk),wtk(Lk),hybik(Lk))

    wini=-Wband/2.d0
    wfin= Wband/2.d0
    epsik=linspace(wini,wfin,Lk,mesh=de)
    !
    test_k=0.d0
    do ix=1,Lk
       wtk(ix)=4.d0/Wband/pi*sqrt(1.d0-(2.d0*epsik(ix)/Wband)**2.d0)*de
       !wtk(ix) = 1.d0/Wband*de
       if(ix==1.or.ix==Lk) wtk(ix)=0.d0
       test_k=test_k+wtk(ix)
       write(77,*) epsik(ix),wtk(ix)
    end do
    hybik=0.d0
    write(*,*) test_k,de

    allocate(Hk_tb(Ns,Ns,Lk))

    Hk_tb=0.d0
    do ik=1,Lk
       do iorb=1,Norb
          do ispin=1,2
             istate=index(ispin,iorb)
             Hk_tb(istate,istate,ik) = epsik(ik)
          end do
       end do
    end do

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
    ! out_unit=free_unit()
    ! open(out_unit,file='orbital_double_occupancy.data')
    ! write(out_unit,'(20F18.10)') gz_docc(1:Norb)
    ! close(out_unit)
    ! !
    ! out_unit=free_unit()
    ! open(out_unit,file='orbital_density_density.data')
    ! do iorb=1,Norb
    !    write(out_unit,'(20F18.10)') gz_dens_dens_orb(iorb,:)
    ! end do
    ! close(out_unit)
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
  !   ispin=1;iorb=1;is=index(ispin,iorb)
  !   Rhop_indep(1)=Rhop_mat(is,is)
  !   ispin=1;iorb=2;is=index(ispin,iorb)
  !   Rhop_indep(2)=Rhop_mat(is,is)
  !   !
  ! end subroutine Rhop_mat2vec



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
  !   ispin=1;iorb=1;is=index(ispin,iorb)
  !   vdm_NC_indep(1)=vdm_NC_mat(is,is)
  !   ispin=1;iorb=2;is=index(ispin,iorb)
  !   vdm_NC_indep(2)=vdm_NC_mat(is,is)
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




end program GUTZ_mb



!AMOEBA TEST


