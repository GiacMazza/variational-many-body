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
  USE MATRIX_SPARSE
  !
  implicit none
  real(8),dimension(:),allocatable :: epsik,hybik
  real(8) :: t
  integer :: Nx,out_unit,is,js,ik,it,itt,i,j,iorb,jorb,ispin,unit,ifock,jfock,iphi,idyn
  integer :: nprint,nsave
  !
  character(len=200) :: store_dir,read_dir,read_optWF_dir
  complex(8),dimension(:,:,:,:),allocatable :: slater_init,slater_t
  complex(8),dimension(:,:),allocatable :: slater_normal
  complex(8),dimension(:,:),allocatable :: slater_anomalous
  complex(8),dimension(:),allocatable     :: gz_proj_init,gz_proj_t
  complex(8),dimension(:,:),allocatable     :: gz_proj_matrix_init
  !
  complex(8),dimension(:),allocatable     :: psi_t
  real(8),dimension(:,:),allocatable      :: Ut,CFt
  real(8),dimension(:),allocatable      :: Jht
  real(8) :: r,s,tmpU,rs !+- time_envelope
  real(8),dimension(3) :: s_orb
  real(8),dimension(:),allocatable :: scf
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
  integer :: unit_neq_weights
  integer :: unit_neq_pair_hopp
  integer :: unit_neq_spin_flip
  integer :: unit_save
  integer :: unit_lgrA_sl
  integer :: unit_lgrA_gz
  !
  !+- observables -+!
  complex(8),dimension(:),allocatable   :: Rhop
  complex(8),dimension(:,:),allocatable   :: Rhop_matrix,Qhop_matrix
  complex(8),dimension(:,:),allocatable :: local_density_matrix
  real(8),dimension(:,:),allocatable    :: local_dens_dens
  complex(8),dimension(:,:,:),allocatable :: dens_constrSL
  complex(8),dimension(:,:,:),allocatable :: dens_constrGZ
  complex(8),dimension(:,:),allocatable :: lgrA_constrSL
  complex(8),dimension(:,:),allocatable :: lgrA_constrGZ
  complex(8),dimension(:,:),allocatable :: pair_hopp
  complex(8),dimension(:,:),allocatable :: spin_flip
  real(8)                               :: unitary_constr
  real(8),dimension(4)                  :: local_angular_momenta
  real(8),dimension(3)                  :: energies
  complex(8),dimension(:,:),allocatable             :: sc_order
  !
  real(8),dimension(:),allocatable      :: dump_vect

  real(8) :: Uneq0,tStart_neqU,tRamp_neqU,tSin_neqU
  real(8) :: tenv,dU_tstop,dU_tcenter
  character(len=3) :: dUenv
  real(8),dimension(3) :: Uneq,dUneq
  real(8) :: Jhneq,Jhneq0,tStart_neqJ,tRamp_neqJ,tSin_neqJ,dJneq,Uasymm


  real(8) :: CFneq0,tStart_neqCF,tRamp_neqCF,tSin_neqCF
  real(8),dimension(3) :: CFneq,dCFneq


  complex(8),dimension(:,:,:),allocatable :: slater_lgr_init,gzproj_lgr_init
  complex(8),dimension(:,:),allocatable   :: R_init,Q_init
  integer :: Nopt
  real(8),dimension(:),allocatable :: dump_seed
  integer :: expected_flen,flen,cf_type
  logical :: seed_file
  logical :: read_full_phi
  character(len=2) :: tdLGR

  complex(8),dimension(:,:,:),allocatable :: td_lgr

  integer,allocatable,dimension(:,:) :: fock_states
  real(8),dimension(:),allocatable :: weights_fock


  integer,dimension(:),allocatable :: tmp_states,states210,ivec
  integer :: count_states,itest,nstate

  integer :: Ntmp,itmp,it0
  !
  real(8) :: tsave,x_re,x_im


  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"WBAND","inputGZ.conf",default=2.d0)
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=1000)
  call parse_input_variable(read_dir,"READ_GZ_BASIS_DIR","inputGZ.conf",default='~/etc_local/GZ_basis/')
  call parse_input_variable(read_optWF_dir,"EQWF_DIR","inputGZ.conf",default='./')
  call parse_input_variable(store_dir,"STORE_GZ_BASIS_DIR","inputGZ.conf",default='./READ_PHI_TRACES/')
  call parse_input_variable(nprint,"NPRINT","inputGZ.conf",default=10)  
  call parse_input_variable(nsave,"NSAVE","inputGZ.conf",default=2000)  
  !
  call parse_input_variable(Uneq,"Uneq","inputGZ.conf",default=[0.d0,0.d0,0.d0])
  call parse_input_variable(Uneq0,"Uneq0","inputGZ.conf",default=0.d0) 
  call parse_input_variable(tStart_neqU,"TSTART_NEQU","inputGZ.conf",default=0.d0)
  call parse_input_variable(tRamp_neqU,"TRAMP_NEQU","inputGZ.conf",default=0.d0)  
  call parse_input_variable(tenv,"TENV","inputGZ.conf",default=0.d0)  
  call parse_input_variable(dU_tstop,"TSTOP","inputGZ.conf",default=0.d0)  
  call parse_input_variable(dU_tcenter,"TCENTER","inputGZ.conf",default=0.d0)  

  call parse_input_variable(dUenv,"DUENV","inputGZ.conf",default='box')  
  call parse_input_variable(tSin_neqU,"TSIN_NEQU","inputGZ.conf",default=0.5d0)
  call parse_input_variable(Uasymm,"UASYMM","inputGZ.conf",default=0.d0) 
  !
  call parse_input_variable(CFneq0,"CFneq0","inputGZ.conf",default=0.d0) 
  call parse_input_variable(CFneq,"CFneq","inputGZ.conf",default=[0.d0,0.d0,0.d0]) 
  call parse_input_variable(dCFneq,"DCFneq","inputGZ.conf",default=[0.d0,0.d0,0.d0]) 
  call parse_input_variable(tStart_neqCF,"TSTART_NEQCF","inputGZ.conf",default=0.d0)
  call parse_input_variable(tRamp_neqCF,"TRAMP_NEQCF","inputGZ.conf",default=0.d0)  
  call parse_input_variable(tSin_neqCF,"TSIN_NEQCF","inputGZ.conf",default=0.5d0)
  call parse_input_variable(cf_type,"CF_TYPE","inputGZ.conf",default=0) 
  !
  call parse_input_variable(Jhneq,"Jhneq","inputGZ.conf",default=0.d0) 
  call parse_input_variable(Jhneq0,"Jhneq0","inputGZ.conf",default=0.d0) 
  call parse_input_variable(tStart_neqJ,"TSTART_NEQJ","inputGZ.conf",default=0.d0)
  call parse_input_variable(tRamp_neqJ,"TRAMP_NEQJ","inputGZ.conf",default=0.d0)  
  call parse_input_variable(tSin_neqJ,"TSIN_NEQJ","inputGZ.conf",default=0.5d0)
  call parse_input_variable(dJneq,"DJneq","inputGZ.conf",default=0.d0) 
  call parse_input_variable(tdLGR,"tdLGR","inputGZ.conf",default='sl')  !possibilitie sl/sg/no
  call parse_input_variable(read_full_phi,"read_phi","inputGZ.conf",default=.false.)

  call parse_input_variable(tsave,"TSAVE","inputGZ.conf",default=0.d0)
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")
  if(Norb.eq.1.and.wf_symmetry.eq.1) then
     write(*,*) 'WARNING THE O(1) x SU(2)c x ORBITAL_ROTATION = O(1) x SU(2)c for the Norb=1 case!'
     wf_symmetry=0
  end if
  !

  call initialize_local_fock_space




  allocate(tmp_states(nFock),ivec(Ns))
  count_states=0
  do ifock=1,NFock
     call bdecomp(ifock,ivec)
     nstate=sum(ivec)
     if(nstate.eq.Norb) then
        itest=0
        do iorb=1,Norb
           itest = itest + ivec(iorb)*ivec(iorb+Norb)
        end do
        if(itest.ne.0) then
           count_states = count_states + 1
           tmp_states(count_states) = ifock
        end if
     end if
  end do
  allocate(states210(count_states))
  states210=tmp_states(1:count_states)
  out_unit=free_unit()
  open(out_unit,file='210_states.info')
  do i=1,count_states
     call bdecomp(states210(i),ivec)
     write(out_unit,*) ivec,'   ',states210(i)
  end do

  call build_lattice_model; get_Hk_t => getHk
  allocate(eLevels(Ns)); eLevels=0.d0
  !
  !
  !
  call init_variational_matrices(wf_symmetry,read_dir_=read_dir)    
  !
  !
  !
  Nsl_normal_opt=2; sl_normal_stride_v2m => sl_normal_vec2mat; sl_normal_stride_m2v => sl_normal_mat2vec
  slNi_v2m => i2m_slN
  slAi_v2m => i2m_slA
  slNi_m2v => m2i_slN
  slAi_m2v => m2i_slA
  !
  Nsl_anomalous_opt=2; sl_anomalous_stride_v2m => sl_anomalous_vec2mat; sl_anomalous_stride_m2v => sl_anomalous_mat2vec
  NRhop_opt=2;   Rhop_stride_v2m => Rhop_vec2mat; Rhop_stride_m2v => Rhop_mat2vec
  NQhop_opt=2;   Qhop_stride_v2m => Qhop_vec2mat; Qhop_stride_m2v => Qhop_mat2vec
  Nvdm_NC_opt=2; vdm_NC_stride_v2m => vdm_NC_vec2mat ; vdm_NC_stride_m2v => vdm_NC_mat2vec
  Nvdm_NCoff_opt=0; vdm_NCoff_stride_v2m => vdm_NCoff_vec2mat ; vdm_NCoff_stride_m2v => vdm_NCoff_mat2vec


  Nvdm_AC_opt=3; vdm_AC_stride_v2m => vdm_AC_vec2mat ; vdm_AC_stride_m2v => vdm_AC_mat2vec
  call get_vdm_AC_mat_index
  !                                                                                                                                                                                  
  Nopt = NRhop_opt + NQhop_opt + Nvdm_NC_opt + Nvdm_NCoff_opt + 2*Nvdm_AC_opt
  Nopt = 2*Nopt
  Nopt_reduced = 2 + 2 + 2 + 2 + 2
  !                                                                                                                                                                                 
  stride_zeros_orig2red => stride2reduced
  stride_zeros_red2orig => stride2orig
  !                                    
  !
  allocate(Rgrid(Ns,Ns),Qgrid(Ns,Ns),Ngrid(Ns,Ns))
  Rgrid=.false.
  Ngrid=.false.
  do is=1,Ns
     Rgrid(is,is)=.true.
     Ngrid(is,is)=.true.  
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
  if(du_tstop==0.d0) dU_tstop = t_grid_aux(NT_aux) - tenv
  if(du_tcenter==0.d0) dU_tcenter = t_grid(NT/2) 

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

     select case(dUenv)
     case('box')
        if(t.lt.tStart_neqU+tRamp_neqU) then
           s_orb = 1.d0
        else                   
           if(t.ge.tenv) then
              !
              if(t.le.dU_tstop) then
                 rs = 1.d0
              else
                 if(t.le.dU_tstop+tenv) then
                    !+- ramp down
                    rs = (1.d0 - 1.5d0*cos(pi*(t-dU_tstop)/tenv) + &
                         0.5d0*(cos(pi*(t-dU_tstop)/tenv)**3))*0.5d0
                    rs = 1.d0 - rs
                 else
                    rs=0.d0
                 end if
              end if
              !
           else
              !
              rs = (1.d0 - 1.5d0*cos(pi*(t-tStart_neqU)/tenv) + 0.5d0*(cos(pi*(t-tStart_neqU)/tenv))**3)*0.5d0
              !
           end if
        end if
     case('exp')
        rs = exp(-(t-dU_tcenter)**2.d0/(2.d0*tenv**2.0))
     end select
     !s_orb(:) = 1.d0 + dUneq(:)*rs*(1.d0-dcos(2.d0*pi*t/tSin_neqU))/2.d0             
     dUneq=[2,1,1]*Uasymm*0.5d0
     s_orb(:) = 1.d0 + dUneq(:)*rs*(1.d0-dcos(2.d0*pi*t/tSin_neqU))/2.d0             
     ! dUneq=[1,0,1]*Uasymm
     ! s_orb(:) = s_orb(:) + dUneq(:)*rs*(1.d0+dcos(2.d0*pi*t/tSin_neqU))/2.d0             
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



  allocate(CFt(Ns,Nt_aux));allocate(scf(Ns))
  CFt=0.d0;scf=0.d0
  do itt=1,Nt_aux
     t = t_grid_aux(itt) 
     !
     if(t.lt.tStart_neqCF) then
        r=0.d0
     else
        if(t.lt.tStart_neqCF+tRamp_neqCF) then
           r = (1.d0 - 1.5d0*cos(pi*(t-tStart_neqCF)/tRamp_neqCF) + 0.5d0*(cos(pi*(t-tStart_neqCF)/tRamp_neqCF))**3)*0.5d0
        else
           r = 1.d0 
        end if
     end if
     !
     if(t.lt.tStart_neqCF+tRamp_neqCF) then
        !
        do iorb=1,Norb
           do ispin=1,2
              is=index(ispin,iorb)
              scf(is) = CFneq(iorb)
           end do
        end do
        !
     else
        !
        select case(cf_type)
        case(0)
           do iorb=1,Norb
              do ispin=1,2
                 is=index(ispin,iorb)
                 scf(is) = CFneq(iorb) + dCFneq(iorb)*dsin(2.d0*pi*t/tSin_neqCF)
              end do
           end do
        case(1)
           do iorb=1,Norb
              do ispin=1,2
                 is=index(ispin,iorb)
                 scf(is) = CFneq(iorb) + dCFneq(iorb)*(1.d0-dcos(2.d0*pi*t/tSin_neqCF))
              end do
           end do
        end select
     end if
     !     
     CFt(:,itt) = CFneq0 + r*(scf(:)-CFneq0)
     !
  end do
  !
  !
  call setup_neq_hamiltonian(Uloc_t_=Ut,Jh_t_=Jht,eLevels_t_=CFt)
  !+- IMPOSE RELATIONS BETWEEN LOCAL INTERACTION PARAMETERS -+!
  !       NORB=3 RATATIONAL INVARIANT HAMILTONIAN       :: Jsf=Jh, Jph=U-Ust-J   (NO relation between Ust and U)
  !       FULLY ROTATIONAL INVARIANT HAMILTONIAN :: Jsf=Jh, Jph=J, Ust = U - 2J   
  !
  Jh_t  = Jh_t
  Jsf_t = Jh_t
  Jph_t = Jh_t
  Ust_t = Uloc_t(1,1)-2.d0*Jh_t ! be careful this has to be done using the first orbital -> we assume Ust is not changed during the excitation !
  !
  unit_neq_hloc = free_unit()
  open(unit_neq_hloc,file="neq_local_interaction_parameters.out")
  !
  do itt=1,Nt_aux
     write(unit_neq_hloc,*) t_grid_aux(itt),Uloc_t(:,itt),Jh_t(itt),Jsf_t(itt),Jph_t(itt),Ust_t(itt),CFt(:,itt)
  end do
  close(unit_neq_hloc)
  !+- READ EQUILIBRIUM AND SETUP DYNAMICAL VECTOR -+!
  nDynamics = 2*Ns*Ns*Lk + Nphi
  allocate(psi_t(nDynamics))
  allocate(slater_init(2,Ns,Ns,Lk),gz_proj_init(Nphi),gz_proj_matrix_init(nFock,nFock))
  allocate(slater_normal(Nsl_normal_opt,Lk),slater_anomalous(Nsl_anomalous_opt,Lk))  
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
        call slater_full2reduced(slater_init,slater_normal,slater_anomalous)
        gz_proj_init = GZ_vector
        allocate(td_lgr(2,Ns,Ns))
        td_lgr(1,:,:) = zero!slater_lgr_init(2,:,:)
        td_lgr(2,:,:) = zero!gzproj_lgr_init(2,:,:)
        !
        call wfMatrix_superc_2_dynamicalVector(slater_init,gz_proj_init,psi_t)
        !
     else
        write(*,*) 'RQn0_root_seed.conf in the wrong form',flen,expected_flen
        write(*,*) 'please check your file for the optimized wavefunction'
        stop
     end if
  else
     !
     allocate(td_lgr(2,Ns,Ns)); td_lgr=zero

     if(read_full_phi) then
        call read_optimized_variational_wf(read_optWF_dir,slater_init,gz_proj_matrix_init)
        call get_phi_basis_decomposition(gz_proj_matrix_init,gz_proj_init)
     else
        call read_optimized_variational_wf(read_optWF_dir,slater_init,gz_proj_matrix_init)
        call get_phi_basis_decomposition(gz_proj_matrix_init,gz_proj_init)
        call read_optimized_variational_wf(read_optWF_dir,slater_init,gz_proj_init)
     end if


     out_unit=free_unit()
     open(out_unit,file='optimized_phi_matrix.data')
     do ifock=1,nFock
        do jfock=1,nFock
           write(out_unit,*) dreal(gz_proj_matrix_init(ifock,jfock)),dimag(gz_proj_matrix_init(ifock,jfock)),ifock,jfock
        end do
     end do
     close(out_unit)

     open(out_unit,file='optimized_projectors.data')
     do iphi=1,Nphi
        write(out_unit,'(2F18.10)') gz_proj_init(iphi)
     end do
     close(out_unit)
     call slater_full2reduced(slater_init,slater_normal,slater_anomalous)
     call wfMatrix_superc_2_dynamicalVector(slater_init,gz_proj_init,psi_t)

     if(tsave/=0.d0) then
        flen=file_length("save_dynamics.data")
        expected_flen = nDynamics
        unit=free_unit()
        open(unit,file="save_dynamics.data",status='old')
        write(*,*) 'reading saved solution'
        if(flen.eq.expected_flen) then
           !+- read from file -+!
           do i=1,flen
              read(unit,*) x_re,x_im
              psi_t(i) = x_re + xi*x_im
           end do
        else
           stop "save_dynamics.data wrong length"
        end if
        close(unit)
     end if
  end if

  it=1
  Uloc=Uloc_t(:,it)
  Ust =Ust_t(it)
  Jh=Jh_t(it)
  Jsf=Jsf_t(it)
  Jph=Jph_t(it)
  eLevels = eLevels_t(:,it)
  call get_local_hamiltonian_trace(eLevels)      
  !
  call setup_neq_dynamics_superc
  !    
  call open_data_files
  !
  allocate(Rhop(Ns));allocate(Rhop_matrix(Ns,Ns))
  allocate(Qhop_matrix(Ns,Ns))
  allocate(local_density_matrix(Ns,Ns))
  allocate(local_dens_dens(Ns,Ns))
  allocate(pair_hopp(Norb,Norb))
  allocate(spin_flip(Norb,Norb))
  allocate(dens_constrSL(2,Ns,Ns))
  allocate(dens_constrGZ(2,Ns,Ns))
  allocate(lgrA_constrSL(Ns,Ns))
  allocate(lgrA_constrGZ(Ns,Ns))
  allocate(sc_order(Ns,Ns))
  allocate(dump_vect(Ns*Ns))
  allocate(gz_proj_t(Nphi))
  !+- local fock weights -+!
  allocate(weights_fock(count_states))  
  !*) ACTUAL DYNAMICS 

  it0=t2it(tsave,tstep)

  do it=it0,Nt
     write(*,*) it,Nt
     !
     t=t_grid(it)
     !
     if(mod(it-1,nprint).eq.0) then        
        !
        call gz_neq_measure_superc_sp(psi_t,t,read_gzproj=gz_proj_t)
        !
        do is=1,Ns
           do js=1,Ns
              call get_neq_Rhop(is,js,Rhop_matrix(is,js))              
              call get_neq_Qhop(is,js,Qhop_matrix(is,js))              
              call get_neq_local_dens(is,js,local_density_matrix(is,js))              
              call get_neq_local_dens_dens(is,js,local_dens_dens(is,js))              
              call get_neq_dens_constr_slater(is,js,dens_constrSL(1,is,js))
              call get_neq_dens_constr_gzproj(is,js,dens_constrGZ(1,is,js))
              call get_neq_dens_constrA_slater(is,js,dens_constrSL(2,is,js))
              call get_neq_dens_constrA_gzproj(is,js,dens_constrGZ(2,is,js))
              call get_neq_lgrA_slater(is,js,lgrA_constrSL(is,js)) 
              call get_neq_lgrA_gzproj(is,js,lgrA_constrGZ(is,js))
              call get_neq_local_sc(is,js,sc_order(is,js))
           end do
        end do
        do iorb=1,Norb
           do jorb=1,Norb
              call get_neq_pair_hopp(iorb,jorb,pair_hopp(iorb,jorb))
              call get_neq_spin_flip(iorb,jorb,spin_flip(iorb,jorb))
           end do
        end do
        call get_neq_energies(energies)
        call get_neq_local_angular_momenta(local_angular_momenta)
        call get_neq_unitary_constr(unitary_constr)
        !
        !-> here it's useless to write all the matrix (only stride allowed one!!)
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
        call write_hermitean_matrix(lgrA_constrSL,unit_lgrA_sl,t)
        call write_hermitean_matrix(lgrA_constrGZ,unit_lgrA_gz,t)
        call write_hermitean_matrix(pair_hopp,unit_neq_pair_hopp,t)
        call write_hermitean_matrix(spin_flip,unit_neq_spin_flip,t)
        call write_symmetric_matrix(local_dens_dens,unit_neq_local_dens_dens,t)
        write(unit_neq_AngMom,'(10F18.10)') t,local_angular_momenta
        write(unit_neq_ene,'(10F18.10)') t,energies
        write(unit_neq_constrU,'(10F18.10)') t,unitary_constr

        !
        do ifock=1,count_states
           weights_fock(ifock)=weight_local_fock_states(states210(ifock),gz_proj_t)
        end do
        write(unit_neq_weights,'(20F18.10)') t,weights_fock(1:count_states)

        !
     end if
     !
     select case(tdlgr)
     case('sl')
        psi_t = RK4_step(nDynamics,4,tstep,t,psi_t,gz_eom_superc_lgrSL)
     case('sg')
        psi_t = RK4_step(nDynamics,4,tstep,t,psi_t,gz_eom_superc_lgrSLGZ)
     case('no')
        psi_t = RK4_step(nDynamics,4,tstep,t,psi_t,gz_equations_of_motion_superc_sp)
     end select
     !
     if(mod(it-1,nsave).eq.0) then
        unit_save=free_unit()
        open(unit_save,file='save_dynamics.data')
        do idyn=1,nDynamics
           write(unit_save,'(2F18.10)') psi_t(idyn)
        end do
        close(unit_save)
        open(unit_save,file='save_dynamics.info')
        write(unit_save,*) 'dynamics saved at time',t_grid(it)
        close(unit_save)
     end if
     !
  end do
  !
CONTAINS
  !

  subroutine open_data_files
    logical :: exist

    unit_neq_Rhop = free_unit()
    inquire(file="neq_Rhop_matrix.data", exist=exist)
    if (exist) then
       open(unit_neq_Rhop,file='neq_Rhop_matrix.data',status="old",position="append",action="write")
    else
       open(unit_neq_Rhop,file='neq_Rhop_matrix.data',status="new",action="write")
    end if
    !
    unit_neq_Qhop = free_unit()
    inquire(file="neq_Qhop_matrix.data", exist=exist)
    if (exist) then
       open(unit_neq_Qhop,file='neq_Qhop_matrix.data',status="old",position="append",action="write")
    else
       open(unit_neq_Qhop,file='neq_Qhop_matrix.data',status="new",action="write")
    end if
    !
    unit_neq_local_dens = free_unit()
    inquire(file="neq_local_density_matrix.data", exist=exist)
    if (exist) then
       open(unit_neq_local_dens,file='neq_local_density_matrix.data',status="old",position="append",action="write")
    else
       open(unit_neq_local_dens,file='neq_local_density_matrix.data',status="new",action="write")
    end if
    !
    unit_neq_local_dens_dens = free_unit()
    inquire(file="neq_local_dens_dens.data", exist=exist)
    if (exist) then
       open(unit_neq_local_dens_dens,file='neq_local_dens_dens.data',status="old",position="append",action="write")
    else
       open(unit_neq_local_dens_dens,file='neq_local_dens_dens.data',status="new",action="write")
    end if
    !
    unit_neq_pair_hopp = free_unit()
    inquire(file="neq_pair_hopp.data", exist=exist)
    if (exist) then
       open(unit_neq_pair_hopp,file='neq_pair_hopp.data',status="old",position="append",action="write")
    else
       open(unit_neq_pair_hopp,file='neq_pair_hopp.data',status="new",action="write")
    end if
    !
    unit_neq_spin_flip = free_unit()
    inquire(file="neq_spin_flip.data", exist=exist)
    if (exist) then
       open(unit_neq_spin_flip,file='neq_spin_flip.data',status="old",position="append",action="write")
    else
       open(unit_neq_spin_flip,file='neq_spin_flip.data',status="new",action="write")
    end if
    !
    unit_neq_ene = free_unit()
    inquire(file="neq_energy.data", exist=exist)
    if (exist) then
       open(unit_neq_ene,file='neq_energy.data',status="old",position="append",action="write")
    else
       open(unit_neq_ene,file='neq_energy.data',status="new",action="write")
    end if
    !
    unit_neq_dens_constrSL = free_unit()
    inquire(file="neq_dens_constrSL.data", exist=exist)
    if (exist) then
       open(unit_neq_dens_constrSL,file='neq_dens_constrSL.data',status="old",position="append",action="write")
    else
       open(unit_neq_dens_constrSL,file='neq_dens_constrSL.data',status="new",action="write")
    end if
    !
    unit_neq_dens_constrGZ = free_unit()
    inquire(file="neq_dens_constrGZ.data", exist=exist)
    if (exist) then
       open(unit_neq_dens_constrGZ,file='neq_dens_constrGZ.data',status="old",position="append",action="write")
    else
       open(unit_neq_dens_constrGZ,file='neq_dens_constrGZ.data',status="new",action="write")
    end if
    !
    unit_neq_dens_constrSLa = free_unit()
    inquire(file="neq_dens_constrSLa.data", exist=exist)
    if (exist) then
       open(unit_neq_dens_constrSLa,file='neq_dens_constrSLa.data',status="old",position="append",action="write")
    else
       open(unit_neq_dens_constrSLa,file='neq_dens_constrSLa.data',status="new",action="write")
    end if
    !open(unit_neq_dens_constrSLa,file='neq_dens_constrSLa.data')
    !
    unit_neq_dens_constrGZa = free_unit()
    inquire(file="neq_dens_constrGZa.data", exist=exist)
    if (exist) then
       open(unit_neq_dens_constrGZa,file='neq_dens_constrGZa.data',status="old",position="append",action="write")
    else
       open(unit_neq_dens_constrGZa,file='neq_dens_constrGZa.data',status="new",action="write")
    end if
    !open(unit_neq_dens_constrGZa,file='neq_dens_constrGZa.data')
    !
    unit_lgrA_sl = free_unit()
    inquire(file="neq_lgrA_constrSL.data", exist=exist)
    if (exist) then
       open(unit_lgrA_sl,file='neq_lgrA_constrSL.data',status="old",position="append",action="write")
    else
       open(unit_lgrA_sl,file='neq_lgrA_constrSL.data',status="new",action="write")
    end if
    !open(unit_lgrA_sl,file='neq_lgrA_constrSL.data')
    !
    unit_lgrA_gz = free_unit()
    inquire(file="neq_lgrA_constrGZ.data", exist=exist)
    if (exist) then
       open(unit_lgrA_gz,file='neq_lgrA_constrGZ.data',status="old",position="append",action="write")
    else
       open(unit_lgrA_gz,file='neq_lgrA_constrGZ.data',status="new",action="write")
    end if
    !open(unit_lgrA_gz,file='neq_lgrA_constrGZ.data')
    !
    unit_neq_constrU = free_unit()
    inquire(file="neq_constrU.data", exist=exist)
    if (exist) then
       open(unit_neq_constrU,file='neq_constrU.data',status="old",position="append",action="write")
    else
       open(unit_neq_constrU,file='neq_constrU.data',status="new",action="write")
    end if
    !open(unit_neq_constrU,file='neq_constrU.data')
    !
    unit_neq_AngMom = free_unit()
    inquire(file="neq_AngMom.data", exist=exist)
    if (exist) then
       open(unit_neq_AngMom,file='neq_AngMom.data',status="old",position="append",action="write")
    else
       open(unit_neq_AngMom,file='neq_AngMom.data',status="new",action="write")
    end if
    !open(unit_neq_AngMom,file='neq_AngMom.data')
    !
    unit_neq_sc_order = free_unit()
    inquire(file="neq_sc_order.data", exist=exist)
    if (exist) then
       open(unit_neq_sc_order,file='neq_sc_order.data',status="old",position="append",action="write")
    else
       open(unit_neq_sc_order,file='neq_sc_order.data',status="new",action="write")
    end if
    !open(unit_neq_sc_order,file='neq_sc_order.data')
    !
    unit_neq_weights = free_unit()
    inquire(file="neq_weights.data", exist=exist)
    if (exist) then
       open(unit_neq_weights,file='neq_weights.data',status="old",position="append",action="write")
    else
       open(unit_neq_weights,file='neq_weights.data',status="new",action="write")
    end if
    !open(unit_neq_weights,file='neq_weights.data')
    

  end subroutine open_data_files

  





  function setup_pointers(NRhop,NQhop,NvdmNC,NvdmNCo,NvdmAC) result(Nopt)
    implicit none
    integer :: Nrhop,NQhop,NvdmNC,NvdmNCo,NvdmAC
    integer :: Nopt
    !
    NRhop_opt=NRhop
    Rhop_stride_v2m => Rhop_vec2mat; Rhop_stride_m2v => Rhop_mat2vec 
    !
    NQhop_opt=NQhop
    Qhop_stride_v2m => Qhop_vec2mat; Qhop_stride_m2v => Qhop_mat2vec
    !
    Nvdm_NC_opt=NvdmNC
    vdm_NC_stride_v2m => vdm_NC_vec2mat ; vdm_NC_stride_m2v => vdm_NC_mat2vec
    !
    Nvdm_NCoff_opt=NvdmNCo 
    vdm_NCoff_stride_v2m => vdm_NCoff_vec2mat ; vdm_NCoff_stride_m2v => vdm_NCoff_mat2vec
    !
    Nvdm_AC_opt=NvdmAC
    vdm_AC_stride_v2m => vdm_AC_vec2mat ; vdm_AC_stride_m2v => vdm_AC_mat2vec
    !    
    Nopt = NRhop_opt + NQhop_opt + Nvdm_NC_opt + Nvdm_NCoff_opt + 2*Nvdm_AC_opt
    Nopt = 2*Nopt
    !
    Nopt_reduced = NRhop_opt + NQhop_opt + Nvdm_NC_opt + Nvdm_NCoff_opt + 2*Nvdm_AC_opt
    !
    stride_zeros_orig2red => stride2reduced
    stride_zeros_red2orig => stride2orig
    !
  end function setup_pointers





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
    implicit none
    complex(8),dimension(:,:) :: Hk
    integer                   :: ik
    real(8)                   :: time
    if(size(Hk,1).ne.size(Hk,2)) stop "wrong dimenions in getHk"
    if(size(Hk,1).ne.Ns) stop "wrong dimenions in getHk"
    Hk = Hk_tb(:,:,ik)
  end subroutine getHk







  !+- STRIDES DEFINITION -+!
  subroutine i2m_slN(iv,imI,imJ)
    implicit none
    integer :: iv
    integer :: imI,imJ
    integer :: iorb,jorb
    if(iv.gt.Nsl_normal_opt) stop "wrong stride!!"
    select case(iv)
    case(1)
       iorb=1
       imI=index(1,iorb)
       imJ=index(1,iorb)
    case(2)
       iorb=2
       imI=index(1,iorb)
       imJ=index(1,iorb)
    end select
  end subroutine i2m_slN
  subroutine i2m_slA(iv,imI,imJ)
    implicit none
    integer :: iv
    integer :: imI,imJ
    integer :: iorb,jorb
    if(iv.gt.Nsl_anomalous_opt) stop "wrong stride!!"
    select case(iv)
    case(1)
       iorb=1
       imI=index(1,iorb)
       imJ=index(2,iorb)
    case(2)
       iorb=2
       imI=index(1,iorb)
       imJ=index(2,iorb)
    end select

  end subroutine i2m_slA
  !
  subroutine m2i_slN(imI,imJ,iv)
    implicit none
    integer:: iv,imI,imJ
    integer :: iorb,jorb,ispin,jspin
    if(imI.gt.Ns.or.imJ.gt.Ns) stop "wrong stride!!"
    !
    ispin= (imI-1)/Norb+1    
    iorb = imI - (ispin-1)*Norb
    jspin= (imJ-1)/Norb+1    
    jorb = imJ - (jspin-1)*Norb
    if(ispin.ne.jspin) then
       iv=0
    else
       if(iorb.ne.jorb) then
          iv = 0
       else
          select case(iorb)
          case(1)
             iv = 1
          case(2)
             iv = 2
          case(3)
             iv = 2             
          end select
       end if
    end if
  end subroutine m2i_slN
  subroutine m2i_slA(iv,imI,imJ)
    implicit none
    integer:: iv,imI,imJ
    integer :: iorb,jorb,ispin,jspin
    if(imI.gt.Ns.or.imJ.gt.Ns) stop "wrong stride!!"
    !
    ispin= (imI-1)/Norb+1    
    iorb = imI - (ispin-1)*Norb
    jspin= (imJ-1)/Norb+1    
    jorb = imJ - (jspin-1)*Norb
    if(ispin.eq.jspin) then
       iv=0
    else
       if(iorb.ne.jorb) then
          iv = 0
       else
          select case(iorb)
          case(1)
             iv = 1
          case(2)
             iv = 2
          case(3)
             iv = 2             
          end select
       end if
    end if
  end subroutine m2i_slA

  subroutine sl_normal_vec2mat(slater_indep,slater_mat)
    implicit none
    complex(8),dimension(:)   :: slater_indep
    complex(8),dimension(:,:) :: slater_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(slater_mat,1).ne.size(slater_mat,2)) stop "wrong stride"
    if(size(slater_mat,1).ne.Ns) stop "wrong stride"
    if(size(slater_indep).ne.Nsl_normal_opt) stop "wrong stride!"    
    slater_mat = zero
    do iorb=1,Norb
       do ispin=1,2
          is=index(ispin,iorb)
          if(iorb.eq.1) then
             slater_mat(is,is) = slater_indep(1)
          else
             slater_mat(is,is) = slater_indep(2)
          end if
       end do
    end do
    !
  end subroutine sl_normal_vec2mat
  subroutine sl_normal_mat2vec(slater_mat,slater_indep)
    implicit none
    complex(8),dimension(:,:) :: slater_mat
    complex(8),dimension(:)   :: slater_indep
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    complex(8) :: test_stride
    real(8) :: test
    if(size(slater_mat,1).ne.size(slater_mat,2)) stop "wrong stride"
    if(size(slater_mat,1).ne.Ns) stop "wrong stride"
    if(size(slater_indep).ne.Nsl_normal_opt) stop "wrong stride!"    
    !
    do iorb=1,Norb
       ispin=1
       is=index(ispin,iorb)
       if(iorb.eq.1) then
          slater_indep(1)=slater_mat(is,is)
       else
          slater_indep(2)=slater_mat(is,is)
       end if
    end do
    !
  end subroutine sl_normal_mat2vec
  !
  subroutine sl_anomalous_vec2mat(slater_indep,slater_mat)
    implicit none
    complex(8),dimension(:)   :: slater_indep
    complex(8),dimension(:,:) :: slater_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(slater_mat,1).ne.size(slater_mat,2)) stop "wrong stride"
    if(size(slater_mat,1).ne.Ns) stop "wrong stride"
    if(size(slater_indep).ne.Nsl_anomalous_opt) stop "wrong stride!"    
    slater_mat = zero
    do iorb=1,Norb
       do jorb=1,Norb
          do ispin=1,2
             jspin=3-ispin
             is=index(ispin,iorb)
             js=index(jspin,jorb)
             if(iorb.eq.jorb) then
                if(iorb.eq.1) then
                   slater_mat(is,js) = (-1.d0)**dble(jspin)*slater_indep(1)
                else
                   slater_mat(is,js) = (-1.d0)**dble(jspin)*slater_indep(2)
                end if
             else
                slater_mat(is,js) = zero
             end if
          end do
       end do
    end do
  end subroutine sl_anomalous_vec2mat
  subroutine sl_anomalous_mat2vec(slater_mat,slater_indep)
    implicit none
    complex(8),dimension(:)   :: slater_indep
    complex(8),dimension(:,:) :: slater_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(slater_mat,1).ne.size(slater_mat,2)) stop "wrong stride"
    if(size(slater_mat,1).ne.Ns) stop "wrong stride"
    if(size(slater_indep).ne.Nsl_anomalous_opt) stop "wrong stride!"    
    !
    iorb=1;jorb=1;ispin=1;jspin=2
    !
    do iorb=1,Norb
       is=index(ispin,iorb)
       js=index(jspin,iorb)
       if(iorb.eq.1) then
          slater_indep(1) = slater_mat(is,js)
       else
          slater_indep(2) = slater_mat(is,js)
       end if
    end do
  end subroutine sl_anomalous_mat2vec




  subroutine Rhop_vec2mat(Rhop_indep,Rhop_mat)
    implicit none
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
    implicit none
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
    implicit none
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
    implicit none
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
    implicit none
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
    implicit none
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
    implicit none
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
    implicit none
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
    implicit none
    complex(8),dimension(:)   :: vdm_AC_indep
    complex(8),dimension(:,:) :: vdm_AC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin,iind
    if(size(vdm_AC_mat,1).ne.size(vdm_AC_mat,2)) stop "wrong stride"
    if(size(vdm_AC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_AC_indep).ne.Nvdm_AC_opt) stop "wrong stride!"    
    !
    vdm_AC_mat = zero
    ! do iorb=1,Norb
    !    do jorb=1,Norb
    !       do ispin=1,2
    !          jspin=3-ispin
    !          is=index(ispin,iorb)
    !          js=index(jspin,jorb)
    !          if(iorb.eq.jorb) then
    !             if(iorb.eq.1) then
    !                vdm_AC_mat(is,js) = (-1.d0)**dble(jspin)*vdm_AC_indep(1)
    !             else
    !                vdm_AC_mat(is,js) = (-1.d0)**dble(jspin)*vdm_AC_indep(2)
    !             end if
    !          else
    !             vdm_AC_mat(is,js) = zero
    !          end if
    !       end do
    !    end do
    ! end do
    iind=0
    do iorb=1,Norb
       iind=iind+1
       do ispin=1,2
          jspin=3-ispin
          is=index(ispin,iorb)
          js=index(jspin,iorb)
          vdm_AC_mat(is,js) = (-1.d0)**dble(jspin)*vdm_AC_indep(iind)
       end do
    end do    
    !
  end subroutine vdm_AC_vec2mat
  subroutine vdm_AC_mat2vec(vdm_AC_mat,vdm_AC_indep)
    implicit none
    complex(8),dimension(:)   :: vdm_AC_indep
    complex(8),dimension(:,:) :: vdm_AC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin,iind
    if(size(vdm_AC_mat,1).ne.size(vdm_AC_mat,2)) stop "wrong stride"
    if(size(vdm_AC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_AC_indep).ne.Nvdm_AC_opt) stop "wrong stride!"    
    !
    iorb=1;jorb=1;ispin=1;jspin=2
    iind=0
    do iorb=1,Norb
       is=index(ispin,iorb)
       js=index(jspin,iorb)
       iind=iind+1
       vdm_AC_indep(iind) = vdm_AC_mat(is,js)
    end do
    !
  end subroutine vdm_AC_mat2vec



  subroutine get_vdm_AC_mat_index
    integer :: iorb,jorb,ispin,jspin

    allocate(IS_vdmAC(Nvdm_AC_opt),JS_vdmAC(Nvdm_AC_opt))
    !
    !
    iorb=1;ispin=1
    jorb=1;jspin=2
    IS_vdmAC(1)=index(ispin,iorb)
    JS_vdmAC(1)=index(jspin,jorb)
    !
    !
    iorb=2;ispin=1
    jorb=2;jspin=2
    IS_vdmAC(2)=index(ispin,iorb)
    JS_vdmAC(2)=index(jspin,jorb)
    !
    !
    iorb=3;ispin=1
    jorb=3;jspin=2
    IS_vdmAC(3)=index(ispin,iorb)
    JS_vdmAC(3)=index(jspin,jorb)
  end subroutine get_vdm_AC_mat_index


  ! subroutine get_vdm_AC_vec_index
  !   integer :: iorb,jorb,ispin,jspin

  !   allocate(I_vdmAC(Ns,Ns)); I_vdmAC=0    
  !   do iorb=1,Norb
  !      do jorb=1,Norb
  !         do ispin=1,2
  !            jspin=3-ispin
  !            is=index(ispin,iorb)
  !            js=index(jspin,jorb)
  !            if(iorb.eq.jorb) then
  !               if(iorb.eq.1) then
  !                  I_vdmAC(is,js) = 1
  !               else
  !                  I_vdmAC(is,js) = 1                   vdm_AC_mat(is,js) = (-1.d0)**dble(jspin)*vdm_AC_indep(2)
  !               end if
  !            else
  !               vdm_AC_mat(is,js) = zero
  !            end if
  !         end do
  !      end do
  !   end do
  !   !
  ! end subroutine get_vdm_AC_vec_index





  subroutine stride2reduced(x_orig,x_reduced)
    implicit none
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
    implicit none
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


