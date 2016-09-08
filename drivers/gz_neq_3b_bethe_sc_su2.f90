program GUTZ_mb
  USE SCIFOR
  !
  USE DMFT_MISC
  USE SF_PARSE_INPUT
  USE RK_IDE
  !
  USE GZ_AUX_FUNX
  USE GZ_neqAUX_FUNX
  USE GZ_DYNAMICS
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_LOCAL_HAMILTONIAN
  USE GZ_ENERGY_MINIMIZATION
  USE GZ_OPTIMIZED_ENERGY
  USE GZ_MATRIX_BASIS
  !
  implicit none
  real(8),dimension(:),allocatable :: epsik,hybik
  real(8) :: t
  integer :: Nx,out_unit,is,js,ik,it,itt,i,j,iorb,ispin
  integer :: nprint
  !
  character(len=200) :: store_dir,read_dir,read_optWF_dir
  complex(8),dimension(:,:,:,:),allocatable :: slater_init
  complex(8),dimension(:),allocatable     :: gz_proj_init
  !
  complex(8),dimension(:),allocatable     :: psi_t
  real(8),dimension(:,:),allocatable      :: Ut 
  real(8),dimension(:),allocatable      :: Jht
  real(8) :: r,s,tmpU
  real(8),dimension(3) :: s_orb
  !
  integer :: unit_neq_hloc
  integer :: unit_neq_local_dens
  integer :: unit_neq_local_dens_dens
  integer :: unit_neq_ene
  integer :: unit_neq_dens_constrSL
  integer :: unit_neq_dens_constrGZ
  integer :: unit_neq_dens_constrSLa
  integer :: unit_neq_dens_constrGZa
  integer :: unit_neq_constrU
  integer :: unit_neq_Rhop
  integer :: unit_neq_Qhop
  integer :: unit_neq_AngMom
  integer :: unit_neq_sc_order
  !
  !+- observables -+!
  complex(8),dimension(:),allocatable   :: Rhop
  complex(8),dimension(:,:),allocatable   :: Rhop_matrix,Qhop_matrix
  complex(8),dimension(:,:),allocatable :: local_density_matrix
  real(8),dimension(:,:),allocatable    :: local_dens_dens
  complex(8),dimension(:,:,:),allocatable :: dens_constrSL
  complex(8),dimension(:,:,:),allocatable :: dens_constrGZ
  real(8)                               :: unitary_constr
  real(8),dimension(4)                  :: local_angular_momenta
  real(8),dimension(3)                  :: energies
  complex(8),dimension(:,:),allocatable             :: sc_order
  !
  real(8),dimension(:),allocatable      :: dump_vect
  
  real(8) :: Uneq0,tStart_neqU,tRamp_neqU,tSin_neqU
  real(8),dimension(3) :: Uneq,dUneq
  real(8) :: Jhneq,Jhneq0,tStart_neqJ,tRamp_neqJ,tSin_neqJ,dJneq



  complex(8),dimension(:,:,:),allocatable :: slater_lgr_init,gzproj_lgr_init
  complex(8),dimension(:,:),allocatable   :: R_init,Q_init
  integer :: Nopt
  real(8),dimension(:),allocatable :: dump_seed
  integer :: expected_flen,flen,unit
  logical :: seed_file


  complex(8),dimension(:,:,:),allocatable :: td_lgr
  
  !
  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"WBAND","inputGZ.conf",default=2.d0)
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=1000)
  call parse_input_variable(read_dir,"READ_GZ_BASIS_DIR","inputGZ.conf",default='~/etc_local/GZ_basis/')
  call parse_input_variable(read_optWF_dir,"EQWF_DIR","inputGZ.conf",default='./')
  call parse_input_variable(store_dir,"STORE_GZ_BASIS_DIR","inputGZ.conf",default='./READ_PHI_TRACES/')
  call parse_input_variable(nprint,"NPRINT","inputGZ.conf",default=10)  
  

  call parse_input_variable(Uneq,"Uneq","inputGZ.conf",default=[0.d0,0.d0,0.d0])
  call parse_input_variable(Uneq0,"Uneq0","inputGZ.conf",default=0.d0) 
  call parse_input_variable(tStart_neqU,"TSTART_NEQU","inputGZ.conf",default=0.d0)
  call parse_input_variable(tRamp_neqU,"TRAMP_NEQU","inputGZ.conf",default=0.d0)  
  call parse_input_variable(tSin_neqU,"TSIN_NEQU","inputGZ.conf",default=0.5d0)
  call parse_input_variable(dUneq,"DUneq","inputGZ.conf",default=[0.d0,0.d0,0.d0]) 
  !
  call parse_input_variable(Jhneq,"Jhneq","inputGZ.conf",default=0.d0) 
  call parse_input_variable(Jhneq0,"Jhneq0","inputGZ.conf",default=0.d0) 
  call parse_input_variable(tStart_neqJ,"TSTART_NEQJ","inputGZ.conf",default=0.d0)
  call parse_input_variable(tRamp_neqJ,"TRAMP_NEQJ","inputGZ.conf",default=0.d0)  
  call parse_input_variable(tSin_neqJ,"TSIN_NEQJ","inputGZ.conf",default=0.5d0)
  call parse_input_variable(dJneq,"DJneq","inputGZ.conf",default=0.d0) 
  
  
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")
  if(Norb.eq.1.and.wf_symmetry.eq.1) then
     write(*,*) 'WARNING THE O(1) x SU(2)c x ORBITAL_ROTATION = O(1) x SU(2)c for the Norb=1 case!'
     wf_symmetry=0
  end if

  !
  call initialize_local_fock_space
  call init_variational_matrices(wf_symmetry,read_dir_=read_dir)    
  !
  call build_lattice_model; get_Hk_t => getHk
  allocate(eLevels(Ns)); eLevels=0.d0
  !

  NRhop_opt=2;   Rhop_stride_v2m => Rhop_vec2mat; Rhop_stride_m2v => Rhop_mat2vec 
  NQhop_opt=2;   Qhop_stride_v2m => Qhop_vec2mat; Qhop_stride_m2v => Qhop_mat2vec
  Nvdm_NC_opt=2; vdm_NC_stride_v2m => vdm_NC_vec2mat ; vdm_NC_stride_m2v => vdm_NC_mat2vec
  Nvdm_NCoff_opt=0; vdm_NCoff_stride_v2m => vdm_NCoff_vec2mat ; vdm_NCoff_stride_m2v => vdm_NCoff_mat2vec
  Nvdm_AC_opt=2; vdm_AC_stride_v2m => vdm_AC_vec2mat ; vdm_AC_stride_m2v => vdm_AC_mat2vec

  Nopt = NRhop_opt + NQhop_opt + Nvdm_NC_opt + Nvdm_NCoff_opt + 2*Nvdm_AC_opt
  Nopt = 2*Nopt
  Nopt_reduced = 2 + 2 + 2 + 2 + 2
  !
  stride_zeros_orig2red => stride2reduced
  stride_zeros_red2orig => stride2orig
  !
  allocate(Rgrid(Ns,Ns),Qgrid(Ns,Ns),Ngrid(Ns,Ns))
  Rgrid=.false.
  Ngrid=.false.
  do is=1,Ns
     Rgrid(is,is)=.true.
     Ngrid(is,js)=.true.
  end do
  Qgrid=.false.
  do iorb=1,Norb
     is=index(1,iorb)
     js=index(2,iorb)
     Qgrid(is,js) = .true.
  end do
  !





  !+- INITIALIZE TIME GRIDS -+!
  Nt_aux=2*Nt+1
  allocate(t_grid(Nt),t_grid_aux(Nt_aux))
  !
  t_grid = linspace(tstart,tstep*real(Nt-1,8),Nt)
  t_grid_aux = linspace(tstart,0.5d0*tstep*real(Nt_aux-1,8),Nt_aux)
  !

  allocate(Ut(3,Nt_aux))  
  do itt=1,Nt_aux
     t = t_grid_aux(itt) 
     !
     if(t.lt.tStart_neqU) then
        r=0.d0
     else
        if(t.lt.tStart_neqU+tRamp_neqU) then
           r = (1.d0 - 1.5d0*cos(pi*(t-tStart_neqU)/tRamp_neqU) + 0.5d0*(cos(pi*(t-tStart_neqU)/tRamp_neqU))**3)*0.5d0
        else
           r = 1.d0 
        end if
     end if
     !
!     if(t.lt.tStart_neqU+tRamp_neqU+tSin_neqU) then
     if(t.lt.tStart_neqU+tRamp_neqU) then
        s_orb = 1.d0
     else
        !      s_orb(:) = 1.d0 + dUneq(:)*dsin(2.d0*pi*t/tSin_neqU)
        s_orb(:) = 1.d0 + dUneq(:)*(1.d0-dcos(2.d0*pi*t/tSin_neqU))/2.d0
     end if
     !
     Ut(:,itt) = Uneq0 + r*(Uneq(:)*s_orb(:)-Uneq0)
  end do
  allocate(Jht(Nt_aux))  
  do itt=1,Nt_aux
     t = t_grid_aux(itt) 
     !
     if(t.lt.tStart_neqJ) then
        r=0.d0
     else
        if(t.lt.tStart_neqJ+tRamp_neqJ) then
           r = (1.d0 - 1.5d0*cos(pi*(t-tStart_neqJ)/tRamp_neqJ) + 0.5d0*(cos(pi*(t-tStart_neqJ)/tRamp_neqJ))**3)*0.5d0
        else
           r = 1.d0 
        end if
     end if
     !
     if(t.lt.tStart_neqU+tRamp_neqJ+tSin_neqJ) then
        s = 1.d0
     else
        s = 1.d0 + dJneq*dsin(2.d0*pi*t/tSin_neqJ)
     end if
     !
     Jht(itt) = Jhneq0 + r*(Jhneq*s-Jhneq0)
  end do
  !
  call setup_neq_hamiltonian(Uloc_t_=Ut,Jh_t_=Jht)
  !+- IMPOSE RELATIONS BETWEEN LOCAL INTERACTION PARAMETERS -+!
  !       NORB=3 RATATIONAL INVARIANT HAMILTONIAN       :: Jsf=Jh, Jph=U-Ust-J   (NO relation between Ust and U)
  !       FULLY ROTATIONAL INVARIANT HAMILTONIAN :: Jsf=Jh, Jph=J, Ust = U - 2J   
  !
  Jh_t  = Jh_t
  Jsf_t = Jh_t
  Jph_t = Jh_t
  Ust_t = Uloc_t(1,:)-2.d0*Jh_t ! be careful this has to be done using the first orbital -> we assume Ust is not changed during the excitation !
  !
  unit_neq_hloc = free_unit()
  open(unit_neq_hloc,file="neq_local_interaction_parameters.out")
  !
  do itt=1,Nt_aux
     write(unit_neq_hloc,*) t_grid_aux(itt),Uloc_t(:,itt),Jh_t(itt),Jsf_t(itt),Jph_t(itt),Ust_t(itt)
  end do
  close(unit_neq_hloc)


  !+- READ EQUILIBRIUM AND SETUP DYNAMICAL VECTOR -+!
  nDynamics = 2*Ns*Ns*Lk + Nphi
  allocate(psi_t(nDynamics))
  allocate(slater_init(2,Ns,Ns,Lk),gz_proj_init(Nphi))  
  !
  expected_flen=Nopt
  inquire(file="RQn0_root_seed.conf",exist=seed_file)
  if(seed_file) then
     flen=file_length("RQn0_root_seed.conf")
     unit=free_unit()
     open(unit,file="RQn0_root_seed.conf")
     write(*,*) 'reading equilibrium solution from file RQn0_root_seed.conf'
     if(flen.eq.expected_flen) then
        allocate(dump_seed(flen))
        !+- read from file -+!
        do i=1,flen
           read(unit,*) dump_seed(i)
        end do
        !+------------------+!
        !
        it = 1
        Uloc = Uloc_t(:,it)
        Ust = Ust_t(it)
        Jh = Jh_t(it)
        Jsf = Jsf_t(it)
        Jph = Jph_t(it)
        eLevels = eLevels_t(:,it)
        !
        call get_local_hamiltonian_trace
        !
        allocate(R_init(Ns,Ns),Q_init(Ns,Ns))
        allocate(slater_lgr_init(2,Ns,Ns),gzproj_lgr_init(2,Ns,Ns))
        slater_lgr_init=0.d0
        gzproj_lgr_init=0.d0
        call dump2mats_superc(dump_seed,R_init,Q_init,slater_lgr_init,gzproj_lgr_init)
        !
        call get_gz_optimized_vdm_Rhop_superc(R_init,Q_init,slater_lgr_init,gzproj_lgr_init)
        !
        call get_gz_ground_state_superc(GZ_vector)  
        !
        slater_init = GZ_opt_slater_superc
        gz_proj_init = GZ_vector
        allocate(td_lgr(2,Ns,Ns))
        td_lgr(1,:,:) = slater_lgr_init(2,:,:)
        td_lgr(2,:,:) = gzproj_lgr_init(2,:,:)
        call wfMatrix_superc_2_dynamicalVector(slater_init,gz_proj_init,psi_t)
     else
        write(*,*) 'RQn0_root_seed.conf in the wrong form',flen,expected_flen
        write(*,*) 'please check your file for the optimized wavefunction'
        stop
     end if
  else
     !
     allocate(td_lgr(2,Ns,Ns)); td_lgr=zero
     call read_optimized_variational_wf_superc(read_optWF_dir,slater_init,gz_proj_init,td_lgr(1,:,:),td_lgr(2,:,:))
     call wfMatrix_superc_2_dynamicalVector(slater_init,gz_proj_init,psi_t)
     it=1
     Uloc=Uloc_t(:,it)
     Ust =Ust_t(it)
     Jh=Jh_t(it)
     Jsf=Jsf_t(it)
     Jph=Jph_t(it)
     eLevels = eLevels_t(:,it)
     call get_local_hamiltonian_trace(eLevels)      
     !
  end if

  
        
  !
  call setup_neq_dynamics_superc
  !    
  unit_neq_Rhop = free_unit()
  open(unit_neq_Rhop,file='neq_Rhop_matrix.data')
  !
  unit_neq_Qhop = free_unit()
  open(unit_neq_Qhop,file='neq_Qhop_matrix.data')  
  !
  unit_neq_local_dens = free_unit()
  open(unit_neq_local_dens,file='neq_local_density_matrix.data')
  !
  unit_neq_local_dens_dens = free_unit()
  open(unit_neq_local_dens_dens,file='neq_local_dens_dens.data')
  !
  unit_neq_ene = free_unit()
  open(unit_neq_ene,file='neq_energy.data')
  !
  unit_neq_dens_constrSL = free_unit()
  open(unit_neq_dens_constrSL,file='neq_dens_constrSL.data')
  !
  unit_neq_dens_constrGZ = free_unit()
  open(unit_neq_dens_constrGZ,file='neq_dens_constrGZ.data')
  !
  unit_neq_dens_constrSLa = free_unit()
  open(unit_neq_dens_constrSLa,file='neq_dens_constrSLa.data')
  !
  unit_neq_dens_constrGZa = free_unit()
  open(unit_neq_dens_constrGZa,file='neq_dens_constrGZa.data')
  !
  unit_neq_constrU = free_unit()
  open(unit_neq_constrU,file='neq_constrU.data')
  !
  unit_neq_AngMom = free_unit()
  open(unit_neq_AngMom,file='neq_AngMom.data')
  !
  unit_neq_sc_order = free_unit()
  open(unit_neq_sc_order,file='neq_sc_order.data')

  allocate(Rhop(Ns));allocate(Rhop_matrix(Ns,Ns))
  allocate(Qhop_matrix(Ns,Ns))
  allocate(local_density_matrix(Ns,Ns))
  allocate(local_dens_dens(Ns,Ns))
  allocate(dens_constrSL(2,Ns,Ns))
  allocate(dens_constrGZ(2,Ns,Ns))  
  allocate(sc_order(Ns,Ns))
  allocate(dump_vect(Ns*Ns))

  !*) ACTUAL DYNAMICS 
  do it=1,Nt
     write(*,*) it,Nt
     !
     t=t_grid(it)
     !
     if(mod(it-1,nprint).eq.0) then        
        !
        call gz_neq_measure_superc_sp(psi_t,t)
        !
        do is=1,Ns
           call get_neq_Rhop(is,is,Rhop(is))
           do js=1,Ns
              call get_neq_Rhop(is,js,Rhop_matrix(is,js))              
              call get_neq_Qhop(is,js,Qhop_matrix(is,js))              
              call get_neq_local_dens(is,js,local_density_matrix(is,js))              
              call get_neq_local_dens_dens(is,js,local_dens_dens(is,js))              
              call get_neq_dens_constr_slater(is,js,dens_constrSL(1,is,js))
              call get_neq_dens_constr_gzproj(is,js,dens_constrGZ(1,is,js))
              call get_neq_dens_constrA_slater(is,js,dens_constrSL(2,is,js))
              call get_neq_dens_constrA_gzproj(is,js,dens_constrGZ(2,is,js))
              call get_neq_local_sc(is,js,sc_order(is,js))
           end do
        end do
        call get_neq_energies(energies)
        call get_neq_local_angular_momenta(local_angular_momenta)
        call get_neq_unitary_constr(unitary_constr)
        !
        call write_complex_matrix_grid(Rhop_matrix,unit_neq_Rhop,print_grid_Rhop,t)
        call write_complex_matrix_grid(Qhop_matrix,unit_neq_Qhop,print_grid_Qhop,t)
        call write_complex_matrix_grid(sc_order,unit_neq_sc_order,print_grid_SC,t)
        !
        call write_hermitean_matrix(local_density_matrix,unit_neq_local_dens,t)
        call write_hermitean_matrix(dens_constrSL(1,:,:),unit_neq_dens_constrSL,t)
        call write_hermitean_matrix(dens_constrGZ(1,:,:),unit_neq_dens_constrGZ,t)
        call write_hermitean_matrix(dens_constrSL(2,:,:),unit_neq_dens_constrSLa,t)
        call write_hermitean_matrix(dens_constrGZ(2,:,:),unit_neq_dens_constrGZa,t)
        call write_hermitean_matrix(local_density_matrix,unit_neq_local_dens,t)
        call write_symmetric_matrix(local_dens_dens,unit_neq_local_dens_dens,t)
        write(unit_neq_AngMom,'(10F18.10)') t,local_angular_momenta
        write(unit_neq_ene,'(10F18.10)') t,energies
        write(unit_neq_constrU,'(10F18.10)') t,unitary_constr
        !
     end if
     !
     call step_dynamics_td_lagrange_superc(nDynamics,tstep,t,psi_t,td_lgr,gz_equations_of_motion_superc_lgr_sp)
     !
  end do
  !
CONTAINS
  !
  subroutine build_lattice_model  
    implicit none
    !
    integer                          :: ix,iy,iz,ik,Nk,iorb,jorb,ispin,istate,jstate
    real(8),allocatable,dimension(:) :: kx
    real(8)                          :: ts,test_k,kx_,ky_,kz_,wini,wfin,de,n1,n2
    !
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
       if(ix==1.or.ix==Lk) wtk(ix)=0.d0
       test_k=test_k+wtk(ix)
       write(77,*) epsik(ix),wtk(ix)
    end do
    hybik=0.d0
    !
    call get_free_dos(epsik,wtk,file='DOS_free.kgrid')
    !
    allocate(Hk_tb(Ns,Ns,Lk))
    !
    Hk_tb=0.d0
    do ik=1,Lk
       do iorb=1,Norb
          do jorb=1,Norb
             do ispin=1,2
                istate=index(ispin,iorb)
                jstate=index(ispin,jorb)
                if(iorb.eq.jorb)  then
                   Hk_tb(istate,jstate,ik) = epsik(ik)
                else
                   Hk_tb(istate,jstate,ik) = hybik(ik)
                end if
             end do
          end do
       end do
    end do
    !
  end subroutine build_lattice_model


  subroutine getHk(Hk,ik,time)
    complex(8),dimension(:,:) :: Hk
    integer                   :: ik
    real(8)                   :: time
    if(size(Hk,1).ne.size(Hk,2)) stop "wrong dimenions in getHk"
    if(size(Hk,1).ne.Ns) stop "wrong dimenions in getHk"
    Hk = Hk_tb(:,:,ik)
  end subroutine getHk







  !+- STRIDES DEFINITION -+!







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
          if(iorb.eq.1) then
             Rhop_mat(is,is) = Rhop_indep(1)
          else
             Rhop_mat(is,is) = Rhop_indep(2)
          end if
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
    do iorb=1,Norb
       ispin=1
       is=index(ispin,iorb)
       if(iorb.eq.1) then
          Rhop_indep(1)=Rhop_mat(is,is)
       else
          Rhop_indep(2)=Rhop_mat(is,is)
       end if
    end do
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
                if(iorb.eq.1) then
                   Qhop_mat(is,js) = (-1.d0)**dble(jspin)*Qhop_indep(1)
                else
                   Qhop_mat(is,js) = (-1.d0)**dble(jspin)*Qhop_indep(2)
                end if
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
    !
    do iorb=1,Norb
       is=index(ispin,iorb)
       js=index(jspin,iorb)
       if(iorb.eq.1) then
          Qhop_indep(1) = Qhop_mat(is,js)
       else
          Qhop_indep(2) = Qhop_mat(is,js)
       end if
    end do
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
    do iorb=1,Norb
       do ispin=1,2
          is=index(ispin,iorb)
          if(iorb.eq.1) then
             vdm_NC_mat(is,is) = vdm_NC_indep(1)
          else
             vdm_NC_mat(is,is) = vdm_NC_indep(2)
          end if
       end do
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
    do iorb=1,Norb
       ispin=1
       is=index(ispin,iorb)
       if(iorb.eq.1) then
          vdm_NC_indep(1) = vdm_NC_mat(is,is)
       else
          vdm_NC_indep(2) = vdm_NC_mat(is,is)
       end if
    end do
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
                if(iorb.eq.1) then
                   vdm_AC_mat(is,js) = (-1.d0)**dble(jspin)*vdm_AC_indep(1)
                else
                   vdm_AC_mat(is,js) = (-1.d0)**dble(jspin)*vdm_AC_indep(2)
                end if
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
    do iorb=1,Norb
       is=index(ispin,iorb)
       js=index(jspin,iorb)
       if(iorb.eq.1) then
          vdm_AC_indep(1) = vdm_AC_mat(is,js)
       else
          vdm_AC_indep(2) = vdm_AC_mat(is,js)
       end if
    end do
    !
  end subroutine vdm_AC_mat2vec


  subroutine stride2reduced(x_orig,x_reduced)
    real(8),dimension(:) :: x_orig
    real(8),dimension(:) :: x_reduced
    if(size(x_orig).ne.Nopt) stop "error orig @ stride2reduced"
    if(size(x_reduced).ne.Nopt_reduced) stop "error reduced @ stride2reduced"
    !+- R
    x_reduced(1) = x_orig(1)
    x_reduced(2) = x_orig(2)
    !+- Q
    x_reduced(3) = x_orig(5)
    x_reduced(4) = x_orig(6)
    !+- LGR_slater normal
    x_reduced(5) = x_orig(9)
    x_reduced(6) = x_orig(10)
    !+- LGR_slater anomalous
    x_reduced(7) = x_orig(13)
    x_reduced(8) = x_orig(14)
    !+- GZproj slater anomalous
    x_reduced(9)  = x_orig(17)
    x_reduced(10) = x_orig(18)    
    !
  end subroutine stride2reduced
  subroutine stride2orig(x_reduced,x_orig)
    real(8),dimension(:) :: x_orig
    real(8),dimension(:) :: x_reduced
    if(size(x_orig).ne.Nopt) stop "error orig @ stride2reduced"
    if(size(x_reduced).ne.Nopt_reduced) stop "error reduced @ stride2reduced"
    x_orig=0.d0
    !+- R
    x_orig(1) = x_reduced(1)
    x_orig(2) = x_reduced(2)
    !+- Q
    x_orig(5) = x_reduced(3)
    x_orig(6) = x_reduced(4)
    !+- LGR slater 
    x_orig(9) = x_reduced(5)
    x_orig(10) = x_reduced(6)
    !+- LGR slater anomalous
    x_orig(13) = x_reduced(7)
    x_orig(14) = x_reduced(8)
    !+- GZproj sanomalous
    x_orig(17) = x_reduced(9)
    x_orig(18) = x_reduced(10)
  end subroutine stride2orig


  
end program GUTZ_mb



!AMOEBA TEST


