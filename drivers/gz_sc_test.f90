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
  complex(8),dimension(:),allocatable   :: R_init,R_out
  complex(8),dimension(:),allocatable   :: Rhop_init,Rhop_out
  complex(8),dimension(:,:),allocatable   :: Rhop_init_matrix
  real(8),dimension(:,:),allocatable :: variational_density_natural_simplex
  real(8),dimension(:,:),allocatable :: variational_density_matrix
  integer                            :: out_unit,iter
  integer                            :: lattice ! 2=square;3=cubic
  real(8),dimension(:),allocatable :: epsik,hybik
  integer :: Nx,is,js
  !
  real(8) :: tmp_emin,Uiter,tmp_ene,orb_pol,tmp_real

  character(len=5) :: dir_suffix
  character(len=6) :: dir_iter

  !
  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=10)
  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"Wband","inputGZ.conf",default=1.d0)
  call parse_input_variable(lattice,"LAT_DIMENSION","inputGZ.conf",default=3)  
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")

  if(Norb.eq.1.and.wf_symmetry.eq.1) then
     write(*,*) 'WARNING THE O(1) x SU(2)c x ORBITAL_ROTATION = O(1) x SU(2)c for the Norb=1 case!'
     wf_symmetry=0
  end if


  !NOTE: ON HUNDS COUPLINGS:
  !NORB=3 RATATIONAL INVARIANT HAMILTONIAN       :: Jsf=Jh, Jph=U-Ust-J   (NO relation between Ust and U)
  !       FULLY ROTATIONAL INVARIANT HAMILTONIAN :: Jsf=Jh, Jph=J, Ust = U - 2J   
  Jh = Jh*Uloc(1)
  Jsf = Jh
  Jph = Jh
  Ust = Uloc(1)-2.d0*Jh
  !
  call initialize_local_fock_space
  !
  !
  call init_variational_matrices
  ! do is=1,Ns
  !    do js=1,Ns        
  !       write(*,*) '!+----------------------------------------------+!'
  !       write(*,*) 'IS',is,'JS',js
  !       tmp_real=0.d0
  !       do i=1,Nphi
  !          do j=1,Nphi
  !             !tmp_real = tmp_real + phi_traces_basis_Rhop(is,js,i,j)*conjg(phi_traces_basis_Rhop(is,js,i,j))
  !             tmp_real = tmp_real + phi_traces_basis_sc_order(is,js,i,j)*conjg(phi_traces_basis_sc_order(is,js,i,j))
  !          end do
  !       end do
  !       write(*,*) tmp_real
  !       tmp_real=0.d0
  !       do i=1,Nphi
  !          do j=1,Nphi
  !             tmp_real = tmp_real + phi_traces_basis_Rhop(is,js,i,j)*conjg(phi_traces_basis_Rhop(is,js,i,j))
  !          end do
  !       end do
  !       write(*,*) tmp_real
  !       tmp_real=0.d0
  !       do i=1,Nphi
  !          do j=1,Nphi
  !             tmp_real = tmp_real + phi_traces_basis_Qhop(is,js,i,j)*conjg(phi_traces_basis_Qhop(is,js,i,j))
  !          end do
  !       end do
  !       write(*,*) tmp_real
  !       tmp_real=0.d0
  !       do i=1,nFock
  !          !           write(*,'(20F6.2)') op_sc_order(is,js,i,:)
  !          do j=1,nFock
  !             tmp_real = tmp_real + op_sc_order(is,js,i,j)*op_sc_order(is,js,i,j)
  !          end do
  !       end do
  !       !        write(*,*) tmp_real

  !       ! do i=1,Nphi
  !       !    
  !       ! end do
  !       ! write(*,*)
  !       ! write(*,*)

  !       ! write(*,*)
  !       ! do i=1,Nphi
  !       !    write(*,'(20F7.3)') dimag(phi_traces_basis_Rhop(is,js,i,:))
  !       ! end do
  !       ! write(*,*)
  !       ! write(*,*) 'xxx'
  !       ! write(*,*)
  !       ! do i=1,Nphi
  !       !    write(*,'(20F7.3)') dreal(phi_traces_basis_Qhop(is,js,i,:))
  !       ! end do
  !       ! write(*,*)
  !       ! do i=1,Nphi
  !       !    write(*,'(20F7.3)') dimag(phi_traces_basis_Qhop(is,js,i,:))
  !       ! end do
  !       ! write(*,*)
  !    end do
  ! end do
  !  
  allocate(variational_density_natural_simplex(Ns+1,Ns))
  allocate(variational_density_natural(Ns))
  call initialize_variational_density_simplex(variational_density_natural_simplex)
  call build_lattice_model
  !
  
  !+ maps for lagrange multipliers -+!
  
  !+- build maps for lagrange multipliers -+!
  Nopt_diag=1
  Nopt_odiag=0

  ! Nopt_diag=1
  ! Nopt_odiag=0
  
  
  Nopt_normal = 1!Nopt_diag + Nopt_odiag
  Nopt_anomalous = 2
  Nopt_lgr = Nopt_normal+Nopt_anomalous
  
  allocate(opt_map(Ns,Ns),opt_map_anomalous(Ns,Ns))
  opt_map = 0
  opt_map_anomalous = 0


  do ispin=1,2
     i=0
     do iorb=1,Norb
        do jorb=1,2
           is=index(ispin,iorb)
           js=index(ispin,jorb)           
           !i = i + 1
           !opt_map(is,js) = i
           if(iorb.eq.jorb) then
              opt_map(is,js) = 1!iorb
           else
              opt_map(is,js) = 0!iorb+jorb              
           end if
        end do
     end do
  end do

  do ispin=1,2
     jspin=3-ispin
     do iorb=1,Norb
        do jorb=1,Norb
           is=index(ispin,iorb)
           js=index(jspin,jorb)
           if(iorb.eq.jorb) then
              opt_map_anomalous(is,js) = 1
           else
              opt_map_anomalous(is,js) = 2
           end if
        end do
     end do
  end do
  !+-----------------------------------------+!

  write(*,*) 'fine codice mai'
  do istate=1,Ns
     write(*,*) opt_map(istate,:)
  end do
  write(*,*)
  do istate=1,Ns
     write(*,*) opt_map_anomalous(istate,:)
  end do
  !stop
  allocate(lgr_init_slater(Nopt_diag+Nopt_odiag))
  allocate(lgr_init_gzproj(Nopt_diag+Nopt_odiag))

  lgr_init_slater=  0.d0
  !
  lgr_init_gzproj=  0.d0

  variational_density_natural = variational_density_natural_simplex(1,1:Ns)
  tmp_emin = gz_energy_vdm(variational_density_natural)
  write(*,*) 'EMIN',tmp_emin

  stop


  call gz_optimization_vdm_simplex(variational_density_natural_simplex,variational_density_natural) 

  call get_gz_ground_state(GZ_vector)  
  call print_output



  !stop

  !

  !  tmp_emin = gz_energy_vdm_Rhop(variational_density_natural,Rhop_init_matrix)


  !  stop

  !allocate(vdm_init(1),vdm_out(1),R_init(1),R_out(1))

  ! vdm_init = variational_density_natural_simplex(Ns,1)
  ! R_init = 0.5d0
  !stop

  allocate(vdm_init(2),vdm_out(2)); 

  vdm_init(1) = 0.05
  vdm_init(2) = 0.95

  ! tmp_emin = gz_energy_vdm(variational_density_natural_simplex(1,:))
  ! write(*,*) 'TMP_EMIN',tmp_emin
  ! stop

  !call gz_optimization_vdm_nlsq(vdm_init,vdm_out)
  ! call get_gz_ground_state(GZ_vector)  
  ! call print_output
  ! stop




  orb_pol = 0.275d0
  do i=1,2
     do iorb=1,2
        do ispin=1,2
           is = index(ispin,iorb)

           select case(iorb)
           case(1)
              variational_density_natural_simplex(:,is) = 0.5d0 - orb_pol
           case(2)
              variational_density_natural_simplex(:,is) = 0.5d0 + orb_pol
           end select
        end do
     end do


     write(dir_suffix,'(F5.3)') orb_pol

     dir_iter="P"//trim(dir_suffix)

     write(*,*) dir_iter

     call system('mkdir -p '//dir_iter)     



     opt_energy_unit=free_unit()
     open(opt_energy_unit,file='GZ_OptEnergy_VS_vdm.out')
     opt_rhop_unit=free_unit()
     open(opt_rhop_unit,file='GZ_OptRhop_VS_vdm.out')
     opt_GZ_unit=free_unit()
     open(opt_GZ_unit,file='GZ_OptProj_VS_vdm.out')
     if(GZmin_verbose) then
        GZmin_unit=free_unit()
        open(GZmin_unit,file='GZ_SelfCons_min_verbose.out')
        GZmin_unit_=free_unit()
        open(GZmin_unit_,file='GZ_proj_min.out')
     end if
     !
     lgr_init_slater(1) =  0.5d0
     lgr_init_slater(2) = -0.5d0
     !
     lgr_init_gzproj(1) =  0.5d0
     lgr_init_gzproj(2) = -0.5d0
     !
     optimization_flag=.true.
     if(.not.allocated(GZ_vector)) allocate(GZ_vector(Nphi))     
     tmp_emin = gz_energy_vdm(variational_density_natural_simplex(1,:)); 


     !call gz_optimization_vdm_simplex(variational_density_natural_simplex,variational_density_natural)      
     !
     call get_gz_ground_state(GZ_vector)
     !
     call print_output

     call system('cp * '//dir_iter)

     call system('rm *.out *.data fort* ')

     ! Rseed=0.d0
     ! do is=1,Ns
     !    Rseed = Rseed + dreal(GZ_opt_Rhop(is,is))/dble(Ns)
     ! end do

     orb_pol = orb_pol + 0.025d0
  end do

  stop


  !call gz_optimization_vdm_Rhop(vdm_init,R_init,vdm_out,R_out)




  allocate(variational_density_matrix(Ns,Ns)); variational_density_matrix=0.d0
  do is=1,Ns
     variational_density_matrix(is,is) = variational_density_natural_simplex(3,is)
  end do
  !
  allocate(Rhop_init_matrix(Ns,Ns)); Rhop_init_matrix=zero
  do is=1,Ns
     Rhop_init_matrix(is,is) = 1.d0
  end do


  Uiter = -0.1

  do i=1,60

     Uiter = Uiter + 0.1
     Uloc(1) = Uiter
     Uloc(2) = Uiter
     Ust = Uiter

     call build_local_hamiltonian

     phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
     phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)



     write(dir_suffix,'(F4.2)') Uiter

     dir_iter="U"//trim(dir_suffix)

     write(*,*) dir_iter

     call system('mkdir '//dir_iter)     

     !
     call gz_optimization_vdm_Rhop(Rhop_init_matrix,variational_density_matrix)

     call get_gz_ground_state(GZ_vector)

     Rhop_init_matrix = GZ_opt_Rhop
     variational_density_matrix = GZ_opt_VDM

     call print_output

     call system('cp * '//dir_iter)

     call system('rm *.out *.data fort* ')

  end do




  write(*,*) "DAJE; quanto e' vera la madonna"



  !variational_density_natural_simplex(1,:)=0.5d0
  !
  !tmp_emin=gz_energy_broyden(variational_density_natural_simplex(1,:))  
  !
  !tmp_emin=gz_energy_recursive_nlep(variational_density_natural_simplex(1,:))
  stop
  !
  !<TEST MINIMIZATION
  ! allocate(vdm_init(Ns),vdm_out(Ns))
  ! allocate(Rhop_init(Ns),Rhop_out(Ns),init_vec(Nphi),Rhop_init_matrix(Ns,Ns))
  ! vdm_init=variational_density_natural_simplex(1,1:Ns)
  ! init_vec=1.d0/sqrt(dble(Nphi))
  ! Rhop_init_matrix=hopping_renormalization_normal(init_vec,vdm_init)  
  ! do is=1,NS
  !    Rhop_init(is) = Rhop_init_matrix(is,is)
  ! end do
  ! Rhop_init=one
  ! call gz_optimization_vdm_Rhop(vdm_init,Rhop_init,vdm_out,Rhop_out)
  ! stop
  !TEST MINIMIZATION
  !

  !
  call gz_optimization_simplex(variational_density_natural_simplex,variational_density_natural)    !

  !
  call get_gz_ground_state_estimation(variational_density_natural)
  !
  call print_output(variational_density_natural_simplex)
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

    ! wini=-Wband/2.d0
    ! wfin= Wband/2.d0
    ! epsik=linspace(wini,wfin,Lk,mesh=de)
    ! !
    ! test_k=0.d0
    ! test_n1=0.d0;test_n2=0.d0
    ! do ix=1,Lk
    !    wtk(ix)=4.d0/Wband/pi*sqrt(1.d0-(2.d0*epsik(ix)/Wband)**2.d0)*de
    !    !wtk(ix) = 1.d0/Wband*de
    !    if(ix==1.or.ix==Lk) wtk(ix)=0.d0
    !    test_n1=test_n1+wtk(ix)*fermi(epsik(ix)+Cfield*0.5d0,beta)
    !    test_n2=test_n2+wtk(ix)*fermi(epsik(ix)-Cfield*0.5d0,beta)
    !    write(77,*) epsik(ix),wtk(ix)
    ! end do
    ! hybik=0.d0
    !write(*,*) test_n1,test_n2,Cfield; stop


    allocate(kx(Nx))
    kx = linspace(0.d0,pi,Nx,.true.,.true.)
    Lk=Nx*Nx*Nx
    allocate(epsik(Lk),wtk(Lk),hybik(Lk))
    ik=0
    do ix=1,Nx
       do iy=1,Nx
          do iz=1,Nx
             ik=ik+1
             !kx_=dble(ix)/dble(Nx)*pi
             epsik(ik) = -2.d0/6.d0*(cos(kx(ix))+cos(kx(iy))+cos(kx(iz))) 
             hybik(ik) = 0.d0/6.d0*(cos(kx(ix))-cos(kx(iy)))*cos(kx(iz)) 
             wtk(ik) = 1.d0/dble(Lk)
          end do
       end do
    end do

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
       ! write(*,*) iphi
       ! do ifock=1,nFock
       !    write(*,'(20F7.1)') phi_basis(iphi,ifock,:)
       ! end do
       ! write(*,*)
       !
       test_full_phi = test_full_phi + GZ_vector(iphi)*phi_basis(iphi,:,:)
    end do
    do ifock=1,nFock
       do jfock=1,nFock
          write(out_unit,*) test_full_phi(ifock,jfock),ifock,jfock
       end do
    end do

    do iphi=1,Nphi
       write(*,*) GZ_vector(iphi)
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
    !write(out_unit,'(20F18.10)') variational_density_natural(1:Ns)
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


end program GUTZ_mb



!AMOEBA TEST


