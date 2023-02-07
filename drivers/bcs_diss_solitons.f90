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
  integer :: Nx,out_unit,is,js,ik,it,itt,i,iorb,ispin
  integer :: nprint
  !
  character(len=200) :: store_dir,read_dir,read_optWF_dir,read_finSC_dir
  complex(8),dimension(:,:,:,:),allocatable :: slater_init
  complex(8),dimension(:),allocatable     :: gz_proj_init
  !
  complex(8),dimension(:,:),allocatable :: bcs_wf
  !
  complex(8),dimension(:),allocatable     :: psi_t,psi_bcs_t,psi_bcs_check
  real(8),dimension(:,:),allocatable      :: Ut 
  real(8),dimension(:),allocatable      :: Jht
  real(8) :: r,s,tmpU,Ubcs,Ubcs0,Ubcsf
  !
  integer :: uio
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
  integer :: unit_neq_nqp
  integer :: unit_neq_bcs
  integer :: unit_proj
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
  complex(8),dimension(:),allocatable :: neq_gzproj
  real(8),dimension(:,:),allocatable :: nqp 
  !
  real(8),dimension(:),allocatable      :: dump_vect
  real(8) :: fin_sc_dir
  real(8) :: Uneq,Uneq0,tStart_neqU,tRamp_neqU,tSin_neqU,dUneq
  real(8) :: tStart_kdiss,tRamp_kdiss,tSin_kdiss
  real(8) :: Jhneq,Jhneq0,tStart_neqJ,tRamp_neqJ,tSin_neqJ,dJneq
  complex(8) :: bcs_sc_order,bcs_delta
  real(8) :: bcs_Kenergy,bcs_Uenergy,phiBCS,bcs_dens,bcs_energy
  real(8) :: sc_phase,bcs_sc_soliton,tmp_bcs,dot_bcs
  logical :: bcs_neq
  logical :: linear_ramp,trpz,flat_dos
  real(8) :: energy_init,delta_pm(2),delta_plus,delta_minus,delta_plus_save,period_save
  real(8) :: tmax,sc_max,esn,ecn,edn,ephi,ck,ce
  integer :: iter,isoliton
  logical :: soliton_solve,skip_bcs
  !
  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"WBAND","inputGZ.conf",default=2.d0)
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=1000)
  call parse_input_variable(read_dir,"READ_GZ_BASIS_DIR","inputGZ.conf",default='~/etc_local/GZ_basis/')
  call parse_input_variable(read_optWF_dir,"EQWF_DIR","inputGZ.conf",default='./')
  call parse_input_variable(read_finSC_dir,"FINSC_DIR","inputGZ.conf",default='./')
  call parse_input_variable(store_dir,"STORE_GZ_BASIS_DIR","inputGZ.conf",default='./READ_PHI_TRACES/')
  call parse_input_variable(nprint,"NPRINT","inputGZ.conf",default=10)  
  call parse_input_variable(bcs_neq,"BCS_NEQ","inputGZ.conf",default=.false.)  
  call parse_input_variable(flat_dos,"FLAT_DOS","inputGZ.conf",default=.false.)  

  call parse_input_variable(sc_seed,"sc_seed","inputGZ.conf",default=0.d0)  

  call parse_input_variable(linear_ramp,"LIN_RAMP","inputGZ.conf",default=.true.)  
  call parse_input_variable(trpz,"TRPZ","inputGZ.conf",default=.false.)  
  !
  call parse_input_variable(Uneq,"Uneq","inputGZ.conf",default=0.d0) 
  call parse_input_variable(Uneq0,"Uneq0","inputGZ.conf",default=0.d0)
  !
  call parse_input_variable(tStart_neqU,"TSTART_NEQU","inputGZ.conf",default=0.d0)
  call parse_input_variable(tRamp_neqU,"TRAMP_NEQU","inputGZ.conf",default=0.d0)
  !
  call parse_input_variable(tStart_kdiss,"TSTART_KDISS","inputGZ.conf",default=0.d0)
  call parse_input_variable(tRamp_kdiss,"TRAMP_KDISS","inputGZ.conf",default=0.d0)
  !
  call parse_input_variable(tSin_neqU,"TSIN_NEQU","inputGZ.conf",default=0.5d0)
  call parse_input_variable(tSin_kdiss,"TSIN_KDISS","inputGZ.conf",default=0.5d0)
  !
  call parse_input_variable(dUneq,"DUneq","inputGZ.conf",default=0.d0) 
  call parse_input_variable(sc_phase,"sc_phase","inputGZ.conf",default=0.d0) 
  call parse_input_variable(soliton_solve,"soliton_solve","inputGZ.conf",default=.false.) 

  !
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")
  !
  Norb=1
  wf_symmetry=4
  !
  call initialize_local_fock_space
  call build_lattice_model; get_Hk_t => getHk
  allocate(eLevels(Ns)); eLevels=0.d0
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
           if(linear_ramp) then
              r = (t-tStart_neqU)/tRamp_neqU
           else
              r = (1.d0 - 1.5d0*cos(pi*(t-tStart_neqU)/tRamp_neqU) + 0.5d0*(cos(pi*(t-tStart_neqU)/tRamp_neqU))**3)*0.5d0
           end if
        else
           r = 1.d0 
        end if
     end if
     !
     if(t.lt.tStart_neqU+tRamp_neqU+tSin_neqU) then
        s = 1.d0
     else
        s = 1.d0 + dUneq*dsin(2.d0*pi*t/tSin_neqU)
     end if
     !
     tmpU = Uneq0 + r*(Uneq*s-Uneq0)
     !
     Ut(:,itt) = Uneq0 + r*(Uneq*s-Uneq0)
  end do
  !

  Ubcs0=Uneq0
  Ubcsf=Uneq
  allocate(Ubcs_t(Nt_aux))
  unit_neq_hloc = free_unit()
  open(unit_neq_hloc,file="neq_Ubcs.out")
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
     if(t.lt.tStart_neqU+tRamp_neqU+tSin_neqU) then
        s = 1.d0
     else
        s = 1.d0 + dUneq*dsin(2.d0*pi*t/tSin_neqU)
     end if
     !
     tmpU = Ubcs0 + r*(Ubcsf*s-Ubcs0)
     !
     Ubcs_t(itt) = tmpU
     if(mod(itt-1,nprint).eq.0) then        
        write(unit_neq_hloc,'(2F18.10)') t,Ubcs_t(itt)
     end if
  end do
  close(unit_neq_hloc)
  !
  !
  !
  k_qp_diss=k_qp_diss*abs(Ubcsf)
  k_qp_loss=k_qp_loss*abs(Ubcsf)
  allocate(kdiss_t(Nt_aux),kpump_t(Nt_aux),kloss_t(Nt_aux))
  unit_neq_hloc = free_unit()
  open(unit_neq_hloc,file="neq_kdiss.out")
  do itt=1,Nt_aux
     t = t_grid_aux(itt) 
     !
     if(t.lt.tStart_kdiss) then
        r=0.d0
     else
        if(t.lt.tStart_kdiss+tRamp_kdiss) then
           r = (1.d0 - 1.5d0*cos(pi*(t-tStart_kdiss)/tRamp_kdiss) + 0.5d0*(cos(pi*(t-tStart_kdiss)/tRamp_kdiss))**3)*0.5d0
        else
           r = 1.d0 
        end if
     end if
     !
     if(t.lt.tStart_kdiss+tRamp_kdiss+tSin_kdiss) then
        s = 1.d0
     else
        s = 1.d0 + dUneq*dsin(2.d0*pi*t/tSin_kdiss)
     end if
     !
     kdiss_t(itt) = r*k_qp_diss
     kpump_t(itt) = r*k_qp_pump
     kloss_t(itt) = r*k_qp_loss
     if(mod(itt-1,nprint).eq.0) then        
        write(unit_neq_hloc,'(5F18.10)') t,kdiss_t(itt),kpump_t(itt),kloss_t(itt)
     end if
  end do
  close(unit_neq_hloc)
  !
  allocate(bcs_wf(3,Lk));bcs_wf=0.d0
  allocate(psi_bcs_t(3*Lk));psi_bcs_t=0.d0
  call init_BCS_wf(bcs_wf,Ubcsf,sc_phase,iprint=.true.)
  delta_pm(1) = 0.d0
  do ik=1,Lk
     delta_pm(1) = delta_pm(1) + 0.5d0*bcs_wf(1,ik)*wtk(ik)
  end do

  bcs_wf=0.d0
  call init_BCS_wf(bcs_wf,Ubcs0,sc_phase,eout=energy_init)
  energy_init=0d0
  delta_pm(2)  = 0.d0
  bcs_sc_order = 0d0
  do ik=1,Lk
     delta_pm(2) = delta_pm(2) + 0.5d0*bcs_wf(1,ik)*wtk(ik)
     energy_init = energy_init + epsik(ik)*bcs_wf(3,ik)*wtk(ik)
     bcs_sc_order = bcs_sc_order + 0.5d0*(bcs_wf(1,ik)+xi*bcs_wf(2,ik))*wtk(ik)
  end do
  !
  energy_init = energy_init - abs(Ubcsf)*abs(bcs_sc_order)**2d0
  write(*,*) delta_pm,energy_init
  !
  !
  !
  call fsolve(solitons_deltas,delta_pm,tol=1d-18,info=iter)

  delta_plus=max(delta_pm(1),delta_pm(2))
  delta_minus=min(delta_pm(1),delta_pm(2))
  !
  !
  call comelp (1.d0-delta_minus**2d0/delta_plus**2d0, ck, ce )
  !
  !
  uio=free_unit()
  open(uio,file="delta_solitons.out")  
  write(uio,*) delta_pm,solitons_deltas(delta_pm),ck,ce,delta_minus,delta_plus
  close(uio)


  call solitons_self_cons(delta_minus,delta_plus)
  ! delta_plus_save=delta_plus
  ! delta_pm=0d0
  ! !
  ! delta_pm(1)=fzero_brentq(solitons_deltas_diss,0d0,delta_plus_save)  
  ! ! call fsolve(,delta_pm(1),tol=1d-18,info=iter)
  ! delta_pm(2) = delta_plus_save
  
  ! uio=free_unit()
  ! open(uio,file="diss_delta_solitons.out")  
  ! write(uio,*) delta_pm
  ! close(uio)

  skip_bcs=.false.
  if(soliton_solve) then
     write(*,*) 'soliton solution w/out dynamics'
     skip_bcs=.true.
  end if
  
  if(.not.skip_bcs) then
     
     call BCSwf_2_dynamicalVector(bcs_wf,psi_bcs_t)  
     !
     !
     !
     allocate(psi_bcs_check(3*Lk))
     call BCSwf_2_dynamicalVector(bcs_wf,psi_bcs_check)  
     
     !
     it=1

     !+- this is all gz stuff -+!
     ! Uloc=Uloc_t(:,it)
     ! Ust =Ust_t(it)
     ! Jh=Jh_t(it)
     ! Jsf=Jsf_t(it)
     ! Jph=Jph_t(it)
     ! eLevels = eLevels_t(:,it)
     ! call get_local_hamiltonian_trace(eLevels)      
     ! !
     ! call setup_neq_dynamics_superc
     ! !    
     ! unit_neq_Rhop = free_unit()
     ! open(unit_neq_Rhop,file='neq_Rhop_matrix.data')
     ! !
     ! unit_neq_Qhop = free_unit()
     ! open(unit_neq_Qhop,file='neq_Qhop_matrix.data')  
     ! !
     ! unit_neq_local_dens = free_unit()
     ! open(unit_neq_local_dens,file='neq_local_density_matrix.data')
     ! !
     ! unit_neq_local_dens_dens = free_unit()
     ! open(unit_neq_local_dens_dens,file='neq_local_dens_dens.data')
     ! !
     ! unit_neq_ene = free_unit()
     ! open(unit_neq_ene,file='neq_energy.data')
     ! !
     ! unit_neq_dens_constrSL = free_unit()
     ! open(unit_neq_dens_constrSL,file='neq_dens_constrSL.data')
     ! !
     ! unit_neq_dens_constrGZ = free_unit()
     ! open(unit_neq_dens_constrGZ,file='neq_dens_constrGZ.data')
     ! !
     ! unit_neq_dens_constrSLa = free_unit()
     ! open(unit_neq_dens_constrSLa,file='neq_dens_constrSLa.data')
     ! !
     ! unit_neq_dens_constrGZa = free_unit()
     ! open(unit_neq_dens_constrGZa,file='neq_dens_constrGZa.data')
     ! !
     ! unit_neq_constrU = free_unit()
     ! open(unit_neq_constrU,file='neq_constrU.data')
     ! !
     ! unit_neq_AngMom = free_unit()
     ! open(unit_neq_AngMom,file='neq_AngMom.data')
     ! !
     ! unit_neq_sc_order = free_unit()
     ! open(unit_neq_sc_order,file='neq_sc_order.data')
     ! !
     ! unit_proj = free_unit()
     ! open(unit_proj,file='neq_proj.data')
     !
     unit_neq_bcs = free_unit()
     
     
     open(unit_neq_bcs,file='diss_bcs.data')
     !
     ! allocate(Rhop(Ns));allocate(Rhop_matrix(Ns,Ns))
     ! allocate(Qhop_matrix(Ns,Ns))
     ! allocate(local_density_matrix(Ns,Ns))
     ! allocate(local_dens_dens(Ns,Ns))
     ! allocate(dens_constrSL(2,Ns,Ns))
     ! allocate(dens_constrGZ(2,Ns,Ns))  
     ! allocate(sc_order(Ns,Ns))
     ! allocate(neq_gzproj(Nphi))
     ! allocate(nqp(Ns,Lk))
     ! allocate(dump_vect(Ns*Ns))
     
     !*) ACTUAL DYNAMICS (simple do loop measuring each nprint times)
     
     sc_max=0d0
     do it=1,Nt
        write(*,*) it,Nt
        !
        t=t_grid(it)
        !
        if(mod(it-1,nprint).eq.0) then        
           !
           ! call gz_neq_measure_superc(psi_t,t,neq_gzproj)
           ! !
           ! do is=1,Ns
           !    call get_neq_Rhop(is,is,Rhop(is))
           !    do js=1,Ns
           !       call get_neq_Rhop(is,js,Rhop_matrix(is,js))              
           !       call get_neq_Qhop(is,js,Qhop_matrix(is,js))              
           !       call get_neq_local_dens(is,js,local_density_matrix(is,js))              
           !       call get_neq_local_dens_dens(is,js,local_dens_dens(is,js))              
           !       call get_neq_dens_constr_slater(is,js,dens_constrSL(1,is,js))
           !       call get_neq_dens_constr_gzproj(is,js,dens_constrGZ(1,is,js))
           !       call get_neq_dens_constrA_slater(is,js,dens_constrSL(2,is,js))
           !       call get_neq_dens_constrA_gzproj(is,js,dens_constrGZ(2,is,js))
           !       call get_neq_local_sc(is,js,sc_order(is,js))
           !    end do
           !    ! do ik=1,Lk
           !    !    call get_neq_nqp(is,ik,nqp(is,ik))
           !    ! end do
           ! end do
           ! !
           ! call get_neq_energies(energies)
           ! call get_neq_local_angular_momenta(local_angular_momenta)
           ! call get_neq_unitary_constr(unitary_constr)
           ! !
           ! call write_complex_matrix_grid(Rhop_matrix,unit_neq_Rhop,print_grid_Rhop,t)
           ! call write_complex_matrix_grid(Qhop_matrix,unit_neq_Qhop,print_grid_Qhop,t)
           ! call write_complex_matrix_grid(sc_order,unit_neq_sc_order,print_grid_SC,t)
           ! !
           ! call write_hermitean_matrix(local_density_matrix,unit_neq_local_dens,t)
           ! call write_hermitean_matrix(dens_constrSL(1,:,:),unit_neq_dens_constrSL,t)
           ! call write_hermitean_matrix(dens_constrGZ(1,:,:),unit_neq_dens_constrGZ,t)
           ! call write_hermitean_matrix(dens_constrSL(2,:,:),unit_neq_dens_constrSLa,t)
           ! call write_hermitean_matrix(dens_constrGZ(2,:,:),unit_neq_dens_constrGZa,t)
           ! call write_hermitean_matrix(local_density_matrix,unit_neq_local_dens,t)
           ! call write_symmetric_matrix(local_dens_dens,unit_neq_local_dens_dens,t)
           ! write(unit_neq_AngMom,'(10F18.10)') t,local_angular_momenta
           ! write(unit_neq_ene,'(10F18.10)') t,energies
           ! write(unit_neq_constrU,'(10F18.10)') t,unitary_constr
           ! write(unit_proj,'(20F18.10)') t,neq_gzproj        
           !        
           !+- measure BCS -+!
           call dynamicalVector_2_BCSwf(psi_bcs_t,bcs_wf)
           bcs_sc_order = zero !<d+d+>
           bcs_Kenergy = zero
           bcs_delta=zero
           bcs_dens=zero
           do ik=1,Lk
              bcs_sc_order = bcs_sc_order + 0.5d0*(bcs_wf(1,ik)+xi*bcs_wf(2,ik))*wtk(ik)
              bcs_dens = bcs_dens + 0.5d0*(bcs_wf(3,ik)+1.d0)*wtk(ik)
              bcs_Kenergy = bcs_Kenergy + epsik(ik)*bcs_wf(3,ik)*wtk(ik)
           end do
           itt=t2it(t,tstep*0.5d0)
           bcs_Uenergy = Ubcs_t(itt)*abs(bcs_sc_order)**2.d0
           bcs_energy = bcs_Kenergy + bcs_Uenergy
           if(abs(bcs_sc_order).gt.sc_max) then
              sc_max=abs(bcs_sc_order)
              tmax=t
           end if
           !bcs_Uenergy = 2.d0*Ubcs_t(itt)*bcs_delta*conjg(bcs_delta)                
           write(unit_neq_bcs,'(30F18.10)') t,abs(bcs_sc_order),dreal(bcs_sc_order),dimag(bcs_sc_order),bcs_dens, bcs_Kenergy,bcs_Uenergy,bcs_energy, & 
                dreal(bcs_wf(1,Lk/2-10))**2d0+dreal(bcs_wf(2,Lk/2-10))**2d0+dreal(bcs_wf(3,Lk/2-10))**2d0, &
                dreal(bcs_wf(1,Lk/2-50))**2d0+dreal(bcs_wf(2,Lk/2-50))**2d0+dreal(bcs_wf(3,Lk/2-50))**2d0, &
                dreal(bcs_wf(1,Lk/2+10))**2d0+dreal(bcs_wf(2,Lk/2+10))**2d0+dreal(bcs_wf(3,Lk/2+10))**2d0, &
                dreal(bcs_wf(1,Lk/2+50))**2d0+dreal(bcs_wf(2,Lk/2+50))**2d0+dreal(bcs_wf(3,Lk/2+50))**2d0
           !
           !,bcs_Kenergy+bcs_Uenergy,bcs_Kenergy,bcs_Uenergy
           !     
           !
           ! call dynamicalVector_2_BCSwf(psi_bcs_check,bcs_wf)
           ! bcs_sc_order = zero !<d+d+>
           ! bcs_dens=zero
           ! do ik=1,Lk
           !    bcs_sc_order = bcs_sc_order + 0.5d0*(bcs_wf(1,ik)+xi*bcs_wf(2,ik))*wtk(ik)
           !    bcs_dens = bcs_dens + 0.5d0*(bcs_wf(3,ik)+1.d0)*wtk(ik)
           ! end do
           ! itt=t2it(t,tstep*0.5d0)
           ! bcs_Uenergy = 2.d0*Ubcs_t(itt)*bcs_delta*conjg(bcs_delta)        
           ! write(746,'(10F18.10)') t,dreal(bcs_sc_order),dimag(bcs_sc_order),bcs_dens
        end if
        psi_bcs_t = RK_step(3*Lk,4,tstep,t,psi_bcs_t,bcs_equations_of_motion)
     end do
     close(unit_neq_bcs)
  end if


  open(unit_neq_bcs,file='soliton_bcs.data')
  tmp_bcs=0d0
  do it=1,Nt
     t=t_grid(it)     
     call jelp((t-tmax)*delta_plus*abs(Ubcsf),1.d0-delta_minus**2d0/delta_plus**2d0,esn,ecn,edn,ephi)
     bcs_sc_soliton = delta_plus*edn
     
     call jelp((t-tmax)*delta_plus*abs(Ubcsf),1.d0-(delta_minus*0.1d0)**2d0/delta_plus**2d0,esn,ecn,edn,ephi)
     ecn=edn*delta_plus

     if(it.gt.1) then 
        dot_bcs=(bcs_sc_soliton-tmp_bcs)/(t_grid(it)-t_grid(it-1))
        write(unit_neq_bcs,'(10F18.10)') t,tmp_bcs,ecn
     end if
     tmp_bcs = bcs_sc_soliton
  end do
  close(unit_neq_bcs)







  !
CONTAINS
  

  subroutine solitons_self_cons(delta_minus,delta_plus)
    real(8),intent(inout) :: delta_minus,delta_plus
    real(8) :: period,ck,ce,delta_pm(2),dm_,dp_
    integer :: iself,is
    !
    !
    !
    
    

    !dm_=delta_minus;dp_=delta_plus
    dm_=0.0d0
    dp_=1d0

    do is=1,1
       isoliton=is
       call comelp (1.d0-dm_**2d0/dp_**2d0, ck, ce )
       period = 2d0*ck/abs(Ubcsf)/delta_plus    
       delta_minus=dm_;delta_plus=dp_
       do iself=1,100
          period_save = period
          
          delta_pm(1) = delta_minus
          delta_pm(2) = delta_plus
          call fsolve(solitons_deltas_diss,delta_pm,tol=1d-18,info=iter)
          delta_plus=max(delta_pm(1),delta_pm(2))
          delta_minus=min(delta_pm(1),delta_pm(2))
          
          call comelp (1.d0-delta_minus**2d0/delta_plus**2d0, ck, ce )
          period = 2d0*ck/abs(Ubcsf)/delta_plus
          
          write(400,*) period,abs(period-period_save),delta_pm,solitons_deltas_diss(delta_pm)
       end do
       write(400,*)
       write(500,*) is,period
    end do


    ! call comelp (1.d0-dm_**2d0/dp_**2d0, ck, ce )
    ! period = 2d0*ck/abs(Ubcsf)/delta_plus    
    ! delta_minus=dm_;delta_plus=dp_
    ! do iself=1,100
    !    period_save = period
       
    !    delta_pm(1) = delta_minus
    !    delta_pm(2) = delta_plus
    !    call fsolve(solitons_deltas_diss_,delta_pm,tol=1d-18,info=iter)
    !    delta_plus=max(delta_pm(1),delta_pm(2))
    !    delta_minus=min(delta_pm(1),delta_pm(2))

    !    call comelp (1.d0-delta_minus**2d0/delta_plus**2d0, ck, ce )
    !    period = 2d0*ck/abs(Ubcsf)/delta_plus
       
    !    write(500,*) period,abs(period-period_save),delta_pm,solitons_deltas_diss(delta_pm)
    ! end do

    !
    !
  end subroutine solitons_self_cons
  !
  !
  !
  function solitons_deltas_diss(deltas) result(self_cons)
    implicit none
    real(8),dimension(:),intent(in) :: deltas
    real(8),dimension(size(deltas)) :: self_cons
    real(8) :: delta_plus,delta_minus,denk,numk
    real(8) :: gamma_diss,gammak
    !    
    gamma_diss = period_save*(2d0*k_qp_loss)*(0.5d0 + dble(isoliton-1) )
    !
    if(size(deltas).ne.2) then
       write(*,*) 'size(deltas).ne.2'
       stop
    end if
    delta_minus=deltas(1);delta_plus=deltas(2)
    !
    self_cons=0d0
    do ik=1,Lk
       !
       denk = Ubcsf**2.d0*(delta_minus**2.d0+delta_plus**2.d0) + 4.d0*epsik(ik)**2.d0 
       denk = denk**2.d0
       denk = denk - 4.d0*Ubcsf**4.d0*delta_minus**2.d0*delta_plus**2.d0
       denk = denk**0.5d0
       !
       if(epsik(ik).lt.0) then
          gammak = 4d0*gamma_diss
       else
          gammak = 0d0*gamma_diss
       end if
       ! gammak = 2d0*gamma_diss
       self_cons(1) = self_cons(1) + 2.d0*epsik(ik)*sign(1.d0,epsik(ik))/denk*wtk(ik)*sqrt(1d0-gammak)
       !
       numk=0.5d0*(delta_plus**2.d0-delta_minus**2.d0)-2.d0*epsik(ik)**2.d0/Ubcsf**2.d0
       self_cons(2) = self_cons(2) + 2.d0*Ubcsf**2.d0*epsik(ik)*sign(1.d0,epsik(ik))/denk*numk*wtk(ik)*sqrt(1d0-gammak)
       !
    end do
    !
    !
    self_cons(1) = abs(Ubcsf)*self_cons(1)-1d0
    self_cons(2) = self_cons(2) - abs(Ubcsf)*delta_plus**2.d0-energy_init*(1d0-dble(isoliton-1)*gamma_diss)
    !
    ! 
  end function solitons_deltas_diss




  function solitons_deltas_diss_(deltas) result(self_cons)
    implicit none
    real(8),dimension(:) :: deltas
    real(8),dimension(size(deltas)) :: self_cons
    real(8) :: delta_plus,delta_minus,denk,numk
    real(8) :: gamma_diss,gammak
    !
    
    gamma_diss = period_save*2d0*k_qp_loss*1.5d0 

    if(size(deltas).ne.2) then
       write(*,*) 'size(deltas).ne.2'
       stop
    end if
    delta_minus=deltas(1);delta_plus=deltas(2)
    !
    self_cons=0d0
    do ik=1,Lk
       !
       denk = Ubcsf**2.d0*(delta_minus**2.d0+delta_plus**2.d0) + 4.d0*epsik(ik)**2.d0 
       denk = denk**2.d0
       denk = denk - 4.d0*Ubcsf**4.d0*delta_minus**2.d0*delta_plus**2.d0
       denk = denk**0.5d0
       !
       if(epsik(ik).lt.0) then
          gammak = 4d0*gamma_diss
       else
          gammak = 0.5d0*gamma_diss
       end if

       self_cons(1) = self_cons(1) + 2.d0*epsik(ik)*sign(1.d0,epsik(ik))/denk*wtk(ik)*sqrt(1d0-gammak)
       !
       numk=0.5d0*(delta_plus**2.d0-delta_minus**2.d0)-2.d0*epsik(ik)**2.d0/Ubcsf**2.d0
       self_cons(2) = self_cons(2) + 2.d0*Ubcsf**2.d0*epsik(ik)*sign(1.d0,epsik(ik))/denk*numk*wtk(ik)*sqrt(1d0-gammak)
       !
    end do
    !
    !
    self_cons(1) = abs(Ubcsf)*self_cons(1)-1d0
    self_cons(2) = self_cons(2) - abs(Ubcsf)*delta_plus**2.d0-energy_init*(1d0-1.d0*gamma_diss)
    !
    ! 
  end function solitons_deltas_diss_
  


  function solitons_deltas(deltas) result(self_cons)
    implicit none
    real(8),dimension(:),intent(in) :: deltas
    real(8),dimension(size(deltas)) :: self_cons
    real(8) :: delta_plus,delta_minus,denk,numk
    !
    if(size(deltas).ne.2) then
       write(*,*) 'size(deltas).ne.2'
       stop
    end if
    delta_minus=deltas(1);delta_plus=deltas(2)
    !
    self_cons=0d0
    do ik=1,Lk
       !
       denk = Ubcsf**2.d0*(delta_minus**2.d0+delta_plus**2.d0) + 4.d0*epsik(ik)**2.d0 
       denk = denk**2.d0
       denk = denk - 4.d0*Ubcsf**4.d0*delta_minus**2.d0*delta_plus**2.d0
       denk = denk**0.5d0
       !
       self_cons(1) = self_cons(1) + 2.d0*epsik(ik)*sign(1.d0,epsik(ik))/denk*wtk(ik)
       !
       numk=0.5d0*(delta_plus**2.d0-delta_minus**2.d0)-2.d0*epsik(ik)**2.d0/Ubcsf**2.d0
       self_cons(2) = self_cons(2) + 2.d0*Ubcsf**2.d0*epsik(ik)*sign(1.d0,epsik(ik))/denk*numk*wtk(ik)
       !
    end do
    self_cons(1) = abs(Ubcsf)*self_cons(1)-1d0
    self_cons(2) = self_cons(2) - abs(Ubcsf)*delta_plus**2.d0-energy_init
    !
    ! 
  end function solitons_deltas
    

  !
  subroutine build_lattice_model  
    implicit none
    
    integer                          :: ix,iy,iz,ik,Nk,iorb,jorb,ispin,istate,jstate
    real(8),allocatable,dimension(:) :: kx
    real(8)                          :: ts,test_k,kx_,ky_,kz_,wini,wfin,de,n1,n2
    !
    !
    Lk=Nx
    allocate(epsik(Lk),wtk(Lk),hybik(Lk))    
    wini=-Wband/2.d0*1.2d0
    wfin= Wband/2.d0*1.2d0
    ! wini=-5.d0
    ! wfin= 5.d0
    epsik=linspace(wini,wfin,Lk,mesh=de)
    !
    test_k=0.d0
    do ix=1,Lk
       if(epsik(ix).gt.-Wband/2.d0.and.epsik(ix).lt.Wband/2.d0) then
          wtk(ix)=4.d0/Wband/pi*sqrt(1.d0-(2.d0*epsik(ix)/Wband)**2.d0)*de
       else
          wtk(ix) = 0.d0
       end if
       if(flat_dos) then
          wtk(ix) = fermi(epsik(ix)-Wband/2d0,1000d0)*fermi(-epsik(ix)-Wband/2d0,1000d0)/Wband*de
       end if

       ! wtk(ix) = 1.d0/Wband*de
       if(ix==1.or.ix==Lk) wtk(ix)=wtk(ix)*0.5d0
       test_k=test_k+wtk(ix)
       write(77,*) epsik(ix),wtk(ix)
    end do
    hybik=0.d0
    wtk=wtk/test_k
    write(*,*) test_k,de
    !
    ! allocate(kx(Nx))
    ! kx = linspace(0.d0,pi,Nx,.true.,.true.)
    ! Lk=Nx*Nx*Nx
    ! allocate(epsik(Lk),wtk(Lk),hybik(Lk))
    ! ik=0
    ! test_k=0.d0;n1=0.d0;n2=0.d0
    ! do ix=1,Nx
    !    do iy=1,Nx
    !       do iz=1,Nx
    !          ik=ik+1
    !          !kx_=dble(ix)/dble(Nx)*pi
    !          epsik(ik) = -2.d0/6.d0*(cos(kx(ix))+cos(kx(iy))+cos(kx(iz))) 
    !          hybik(ik) = 1.d0/6.d0*(cos(kx(ix))-cos(kx(iy)))*cos(kx(iz)) 
    !          wtk(ik) = 1.d0/dble(Lk)
    !          n1=n1+fermi(epsik(ik)+cfield*0.5,beta)*wtk(ik)
    !          n2=n2+fermi(epsik(ik)-cfield*0.5,beta)*wtk(ik)
    !       end do
    !    end do
    ! end do
    ! !
    ! write(*,*) 'n1/n2'
    ! write(*,*) n1,n2,n1+n2
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







  subroutine print_output(vdm_simplex,vdm_opt)
    real(8),dimension(Ns+1,Ns),optional :: vdm_simplex
    real(8),dimension(Nvdm_NC_opt-Nvdm_NCoff_opt),optional :: vdm_opt
    integer :: out_unit,istate,jstate,iorb,iphi,ifock,jfock,is,js,ik
    integer,dimension(Ns) :: fock_state
    complex(8),dimension(Ns) :: tmp
    real(8) :: deltani,delta_tmp,vdm_tmp

    real(8),dimension(nFock,nFock) :: test_full_phi

    !+- STORE GZ PROJECTORS -+!
    out_unit=free_unit()
    open(out_unit,file='optimized_projectors.data')
    test_full_phi=0.d0
    do iphi=1,Nphi
       !
       test_full_phi = test_full_phi + GZ_vector(iphi)*phi_basis(iphi,:,:)
       write(out_unit,'(2F18.10)') GZ_vector(iphi)
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

    !+- STORE SLATER DETERMINANT -+!
    do is=1,Ns
       do js=1,Ns
          out_unit=free_unit()
          open(out_unit,file='optimized_slater_normal_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data')
          do ik=1,Lk
             write(out_unit,'(2F18.10)') dreal(GZ_opt_slater(is,js,ik)),dimag(GZ_opt_slater(is,js,ik))
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


  subroutine init_bcs_wf(bcs_wf,U,phase,eout,iprint)
    complex(8),dimension(:,:) :: bcs_wf
    real(8) :: U,bcs_sc_order,phase
    real(8) :: Ek
    real(8) :: sintk,costk,sinph,cosph
    complex(8) :: bcs_phi
    logical,optional :: iprint
    real(8),optional :: eout
    real(8)          :: eout_
    real(8)   :: sc_test
    real(8)  :: braket(2),rho0
    logical :: iprint_
    integer :: uio
    if(size(bcs_wf,1).ne.3) stop "error init bcs \sigma"
    if(size(bcs_wf,2).ne.Lk) stop "error init bcs Lk"
    !
    iprint_=.false.
    if(present(iprint)) iprint_=iprint
    !
    Ubcs=U
    !
    write(*,*) 'solving BCS'
    braket(1) = bcs_self_cons(1d-9)
    braket(2) = bcs_self_cons(1d0)    
    !
    if(braket(1)*braket(2).le.0d0) then
       bcs_sc_order=brentq(bcs_self_cons,1.d-9,1.d0)
    else
       if(U.lt.1d-8) then
          bcs_sc_order = 1e-10
       else
          bcs_sc_order = Wband/abs(Ubcs)*exp(-Wband/abs(Ubcs))
          write(*,*) 'approximate solution'
          stop

       end if
    end if
    !
    !
    bcs_phi = bcs_sc_order*exp(xi*phase*pi)
    sinph   = dimag(bcs_phi);cosph=dreal(bcs_phi)
    eout_   = 0d0
    !
    !
    sc_test=0d0
    do ik=1,Lk
       !
       Ek = sqrt(epsik(ik)**2.d0 + (bcs_sc_order*Ubcs)**2.d0)
       sintk=Ubcs/Ek
       costk=epsik(ik)/Ek
       !
       bcs_wf(1,ik) = -cosph*sintk!*tanh(beta*Ek*0.5d0)       
       bcs_wf(2,ik) = -sinph*sintk!*tanh(beta*Ek*0.5d0)
       bcs_wf(3,ik) = -costk!*tanh(beta*Ek*0.5d0)
       !
       eout_ = eout_ + epsik(ik)*bcs_wf(3,ik)*wtk(ik)
       !
    end do
    eout_ = eout_ - abs(Ubcs)*abs(bcs_phi)**2.d0
    if(present(eout)) eout=eout_
    
    if(iprint_) then
       uio=free_unit()
       open(unit=uio,file='bcs_equ.out')
       write(uio,'(10F18.10)') U,bcs_sc_order,U*bcs_sc_order
       write(uio,'(10F18.10)')
       do ik=1,Lk
          write(uio,'(10F18.10)') epsik(ik),bcs_wf(:,ik)           
       end do
       close(uio)
    end if
    ! write(*,*) 'tmp test; init bcs',bcs_sc_order,Ubcs,sc_test
    ! stop
    !
  end subroutine init_bcs_wf

  function bcs_self_cons(phi) result(x)
    real(8),intent(in) :: phi
    real(8) :: x
    real(8) :: Ek
    integer :: ik
    !
    x=0.d0
    do ik=1,Lk
       Ek = sqrt(epsik(ik)**2.d0 + (phi*Ubcs)**2.d0)
       x = x + wtk(ik)/Ek!*tanh(beta*Ek*0.5d0)
    end do
    x=x*abs(Ubcs)*0.5d0
    x=x-1.d0
    !
    ! write(*,*) 
    !
  end function bcs_self_cons

  function bcs_order_param(xU) result(phi)
    real(8),intent(in) :: xU
    real(8) :: phi
    !
    Ubcs=xU
    phi=brentq(bcs_self_cons,0.d0,1.d0)    
    !
  end function bcs_order_param

  subroutine getUbcs(phi,xu)
    real(8) :: phi,xu,tmp
    phiBCS=phi
    ! tmp=delta_bcs_order_param(-10.d0)
    ! write(*,*) phiBCS,tmp
    ! tmp=delta_bcs_order_param(-0.2d0)
    ! write(*,*) phiBCS,tmp
    xu = brentq(delta_bcs_order_param,-10.d0,-0.2d0)
    !write(*,*) xu
  end subroutine getUbcs


  function delta_bcs_order_param(xU) result(phi)
    real(8),intent(in) :: xU
    real(8) :: phi
    !
    Ubcs=xU
    phi=brentq(bcs_self_cons,0.d0,1.d0)    
    phi=phi-phiBCS
    !
  end function delta_bcs_order_param





  subroutine read_SC_order(read_dir,sc_order)
    character(len=200)             :: read_dir
    real(8)                        :: tmp,sc_order
    integer :: is,js,read_unit,flen,ios,i
    character(len=200) :: file_name
    logical :: check_file
    !
    file_name=reg(read_dir)//'local_sc_order.data'
    inquire(file=file_name,exist=check_file)
    if(check_file) then
       read_unit=free_unit()
       open(read_unit,file=file_name,status='old')
       flen=0
       do
          read (read_unit,*,iostat=ios) tmp
          if (ios/=0) exit     
          flen=flen+1
       end do
       close(read_unit)                    
       if(flen.ne.Ns) stop "READING SC ORDER PARAM: number of lines /= Ns"
       open(read_unit,file=file_name,status='old')
       do i=1,flen
          read(read_unit,'(2F18.10)') tmp
          if(i.eq.flen) sc_order = sqrt(tmp*tmp)
       end do
       close(read_unit)
    else
       write(*,*) 'FILE',file_name,'does not exist!!!'
    end if
    !
  end subroutine read_SC_order


end program GUTZ_mb



!AMOEBA TEST


