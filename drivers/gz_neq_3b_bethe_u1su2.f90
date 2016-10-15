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
  integer :: Nx,out_unit,is,js,ik,it,itt,i,j,iorb,ispin,unit,ifock,jfock,iphi
  integer :: nprint
  !
  character(len=200) :: store_dir,read_dir,read_optWF_dir
  complex(8),dimension(:,:,:),allocatable :: slater_init,slater_t
  complex(8),dimension(:,:),allocatable :: slater_normal
  complex(8),dimension(:),allocatable     :: gz_proj_init,gz_proj_t
  complex(8),dimension(:,:),allocatable     :: gz_proj_matrix_init
  !
  complex(8),dimension(:),allocatable     :: psi_t
  real(8),dimension(:,:),allocatable      :: Ut,CFt
  real(8),dimension(:),allocatable      :: Jht
  real(8) :: r,s,tmpU
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
  !
  !+- observables -+!
  complex(8),dimension(:,:),allocatable   :: Rhop_matrix
  complex(8),dimension(:,:),allocatable :: local_density_matrix
  real(8),dimension(:,:),allocatable    :: local_dens_dens
  complex(8),dimension(:,:),allocatable :: dens_constrSL
  complex(8),dimension(:,:),allocatable :: dens_constrGZ
  real(8)                               :: unitary_constr
  real(8),dimension(4)                  :: local_angular_momenta
  real(8),dimension(3)                  :: energies
  complex(8),dimension(:,:),allocatable             :: sc_order
  !
  real(8),dimension(:),allocatable      :: dump_vect
  
  real(8) :: Uneq0,tStart_neqU,tRamp_neqU,tSin_neqU
  real(8),dimension(3) :: Uneq,dUneq
  real(8) :: Jhneq,Jhneq0,tStart_neqJ,tRamp_neqJ,tSin_neqJ,dJneq

  real(8) :: CFneq0,tStart_neqCF,tRamp_neqCF,tSin_neqCF
  real(8),dimension(3) :: CFneq,dCFneq
  complex(8),dimension(:,:),allocatable   :: R_init
  integer :: Nopt
  real(8),dimension(:),allocatable :: dump_seed
  integer :: expected_flen,flen,cf_type
  logical :: seed_file
  logical :: tdLGR

  complex(8),dimension(:,:,:),allocatable :: td_lgr

  integer,allocatable,dimension(:,:) :: fock_states
  real(8),dimension(:),allocatable :: weights_fock


  integer,dimension(:),allocatable :: tmp_states,states210,ivec
  integer :: count_states,itest,nstate
  complex(8),dimension(:,:),allocatable :: slater_lgr_init,gzproj_lgr_init
  integer :: Ntmp,itmp
  !



  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"WBAND","inputGZ.conf",default=2.d0)
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=1000)
  call parse_input_variable(read_dir,"READ_GZ_BASIS_DIR","inputGZ.conf",default='~/etc_local/GZ_basis/')
  call parse_input_variable(read_optWF_dir,"EQWF_DIR","inputGZ.conf",default='./')
  call parse_input_variable(store_dir,"STORE_GZ_BASIS_DIR","inputGZ.conf",default='./READ_PHI_TRACES/')
  call parse_input_variable(nprint,"NPRINT","inputGZ.conf",default=10)  
  !
  call parse_input_variable(Uneq,"Uneq","inputGZ.conf",default=[0.d0,0.d0,0.d0])
  call parse_input_variable(Uneq0,"Uneq0","inputGZ.conf",default=0.d0) 
  call parse_input_variable(tStart_neqU,"TSTART_NEQU","inputGZ.conf",default=0.d0)
  call parse_input_variable(tRamp_neqU,"TRAMP_NEQU","inputGZ.conf",default=0.d0)  
  call parse_input_variable(tSin_neqU,"TSIN_NEQU","inputGZ.conf",default=0.5d0)
  call parse_input_variable(dUneq,"DUneq","inputGZ.conf",default=[0.d0,0.d0,0.d0]) 
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
  !call parse_input_variable(tdLGR,"tdLGR","inputGZ.conf",default=.true.) 
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
  !
  !
  !
  call init_variational_matrices(wf_symmetry,read_dir_=read_dir)    
  !
  !
  !
  call build_lattice_model; get_Hk_t => getHk
  allocate(eLevels(Ns)); eLevels=0.d0
  !
  !
  !  
  NRhop_opt=2;   Rhop_stride_v2m => Rhop_vec2mat; Rhop_stride_m2v => Rhop_mat2vec
  Nvdm_NC_opt=2; vdm_NC_stride_v2m => vdm_NC_vec2mat ; vdm_NC_stride_m2v => vdm_NC_mat2vec
  Nvdm_NCoff_opt=0; vdm_NCoff_stride_v2m => vdm_NCoff_vec2mat ; vdm_NCoff_stride_m2v => vdm_NCoff_mat2vec
  !                                                                                                                                                                                  
  Nopt = NRhop_opt + Nvdm_NC_opt + Nvdm_NCoff_opt 
  Nopt = 2*Nopt
  Nopt_reduced = 2 + 2
  !
  stride_zeros_orig2red => stride2reduced
  stride_zeros_red2orig => stride2orig
  !                                    
  !
  allocate(Rgrid(Ns,Ns),Ngrid(Ns,Ns))
  Rgrid=.false.
  Ngrid=.false.
  do is=1,Ns
     Rgrid(is,is)=.true.
     Ngrid(is,is)=.true.  
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
     if(t.lt.tStart_neqU+tRamp_neqU) then
        s_orb = 1.d0
     else
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
  call setup_neq_hamiltonian(Uloc_t_=Ut,Jh_t_=Jht,eLevels_t_=CFt)
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
     write(unit_neq_hloc,*) t_grid_aux(itt),Uloc_t(:,itt),Jh_t(itt),Jsf_t(itt),Jph_t(itt),Ust_t(itt),CFt(:,itt)
  end do
  close(unit_neq_hloc)

  !+- READ EQUILIBRIUM AND SETUP DYNAMICAL VECTOR -+!
  nDynamics = Ns*Ns*Lk + Nphi
  allocate(psi_t(nDynamics))
  allocate(slater_init(Ns,Ns,Lk),gz_proj_init(Nphi),gz_proj_matrix_init(nFock,nFock))
  !
  expected_flen=Nopt
  inquire(file="Rn0_root_seed.conf",exist=seed_file)
  if(seed_file) then
     flen=file_length("Rn0_root_seed.conf")
     unit=free_unit()
     open(unit,file="Rn0_root_seed.conf")
     write(*,*) 'reading equilibrium solution from file Rn0_root_seed.conf'
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
        allocate(R_init(Ns,Ns))
        allocate(slater_lgr_init(Ns,Ns),gzproj_lgr_init(Ns,Ns))
        slater_lgr_init=0.d0
        gzproj_lgr_init=0.d0
        call dump2mats(dump_seed,R_init,slater_lgr_init,gzproj_lgr_init)
        !
        call get_gz_optimized_vdm_Rhop(R_init,slater_lgr_init,gzproj_lgr_init)
        !
        call get_gz_ground_state_superc(GZ_vector)  
        !
        slater_init = GZ_opt_slater
        gz_proj_init = GZ_vector
        call wfMatrix_2_dynamicalVector(slater_init,gz_proj_init,psi_t)
        !
     else
        write(*,*) 'Rn0_root_seed.conf in the wrong form',flen,expected_flen
        write(*,*) 'please check your file for the optimized wavefunction'
        stop
     end if
  else
     !
     call read_optimized_variational_wf(read_optWF_dir,slater_init,gz_proj_matrix_init)
     call get_phi_basis_decomposition(gz_proj_matrix_init,gz_proj_init)
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
    call wfMatrix_2_dynamicalVector(slater_init,gz_proj_init,psi_t)
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
  call setup_neq_dynamics
  !    
  unit_neq_Rhop = free_unit()
  open(unit_neq_Rhop,file='neq_Rhop_matrix.data')
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
  unit_neq_constrU = free_unit()
  open(unit_neq_constrU,file='neq_constrU.data')
  !
  unit_neq_AngMom = free_unit()
  open(unit_neq_AngMom,file='neq_AngMom.data')
  !
  unit_neq_sc_order = free_unit()
  open(unit_neq_sc_order,file='neq_sc_order.data')
  !
  unit_neq_weights = free_unit()
  open(unit_neq_weights,file='neq_weights.data')
  !
  allocate(Rhop_matrix(Ns,Ns))
  allocate(local_density_matrix(Ns,Ns))
  allocate(local_dens_dens(Ns,Ns))
  allocate(dens_constrSL(Ns,Ns))
  allocate(dens_constrGZ(Ns,Ns))  
  allocate(dump_vect(Ns*Ns))
  allocate(gz_proj_t(Nphi))
  allocate(weights_fock(count_states))
  
  
  
  !*) ACTUAL DYNAMICS 
  do it=1,Nt
     write(*,*) it,Nt
     !
     t=t_grid(it)
     !
     if(mod(it-1,nprint).eq.0) then        
        !
        call gz_neq_measure_sp(psi_t,t,read_gzproj=gz_proj_t)
        !
        do is=1,Ns
           do js=1,Ns
              call get_neq_Rhop(is,js,Rhop_matrix(is,js))              
              call get_neq_local_dens(is,js,local_density_matrix(is,js))              
              call get_neq_local_dens_dens(is,js,local_dens_dens(is,js))              
              call get_neq_dens_constr_slater(is,js,dens_constrSL(is,js))
              call get_neq_dens_constr_gzproj(is,js,dens_constrGZ(is,js))
           end do
        end do
        call get_neq_energies(energies)
        call get_neq_local_angular_momenta(local_angular_momenta)
        call get_neq_unitary_constr(unitary_constr)
        !
        call write_complex_matrix_grid(Rhop_matrix,unit_neq_Rhop,print_grid_Rhop,t)
        call write_hermitean_matrix(local_density_matrix,unit_neq_local_dens,t)
        call write_hermitean_matrix(dens_constrSL,unit_neq_dens_constrSL,t)
        call write_hermitean_matrix(dens_constrGZ,unit_neq_dens_constrGZ,t)
        call write_symmetric_matrix(local_dens_dens,unit_neq_local_dens_dens,t)
        write(unit_neq_AngMom,'(10F18.10)') t,local_angular_momenta
        write(unit_neq_ene,'(10F18.10)') t,energies
        write(unit_neq_constrU,'(10F18.10)') t,unitary_constr
        do ifock=1,count_states
           weights_fock(ifock)=weight_local_fock_states(states210(ifock),gz_proj_t)
        end do        
        write(unit_neq_weights,'(20F18.10)') t,weights_fock(1:count_states)
     end if
     !
     psi_t = RK4_step(nDynamics,4,tstep,t,psi_t,gz_equations_of_motion_sp)
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
  !
  subroutine getHk(Hk,ik,time)
    implicit none
    complex(8),dimension(:,:) :: Hk
    integer                   :: ik
    real(8)                   :: time
    if(size(Hk,1).ne.size(Hk,2)) stop "wrong dimenions in getHk"
    if(size(Hk,1).ne.Ns) stop "wrong dimenions in getHk"
    Hk = Hk_tb(:,:,ik)
  end subroutine getHk
  !
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
  !
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
  subroutine stride2reduced(x_orig,x_reduced)
    implicit none
    real(8),dimension(:) :: x_orig
    real(8),dimension(:) :: x_reduced
    if(size(x_orig).ne.Nopt) stop "error orig @ stride2reduced"
    if(size(x_reduced).ne.Nopt_reduced) stop "error reduced @ stride2reduced"
    !+- R
    x_reduced(1) = x_orig(1)
    x_reduced(2) = x_orig(2)
    !+- LGR_slater normal
    x_reduced(3) = x_orig(5)
    x_reduced(4) = x_orig(6)
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
    !+- LGR slater 
    x_orig(5) = x_reduced(3)
    x_orig(6) = x_reduced(4)
  end subroutine stride2orig


  
end program GUTZ_mb



!AMOEBA TEST


