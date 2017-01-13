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
  integer                            :: ispin,jspin,iorb,i,j,istate,jstate,ifock,jorb,is,js,jfock,dim,dim_s,ieigen
  integer,dimension(:),allocatable   :: fock_vec
  complex(8),dimension(:),allocatable               :: init_vec
  integer  :: Nsectors,isect,iout
  !
  integer,dimension(:),allocatable :: tmp_states,states210,ivec
  integer :: count_states,itest,nstate
  real(8),dimension(:),allocatable :: eigenv_sector
  real(8),dimension(:,:),allocatable :: hamiltonian_sector
  real(8),dimension(:),allocatable :: tmp_m0,tmp_mp,tmp_mm
  real(8),dimension(:,:),allocatable :: deltaH,tmpH
  real(8),dimension(:,:),allocatable :: S_up,T_dw,S_up_,T_up_
  logical :: exist
  !
  integer :: iflo,jflo,Nflo,m_up,m_dw
  integer,dimension(:),allocatable :: flo_sect
  !
  type sector
     integer :: dim_sect
     integer,dimension(:),allocatable :: states_sect
     real(8),dimension(:),allocatable :: E_sect
     real(8),dimension(:,:),allocatable :: V_sect
  end type sector
  type(sector),dimension(:),allocatable :: sectors
  type flo_spectrum
     integer :: dim
     integer :: dim_s
     integer :: dim_flo
     real(8),dimension(:),allocatable :: E_flo
     real(8),dimension(:,:),allocatable :: V_flo
  end type flo_spectrum
  type(flo_spectrum),dimension(:),allocatable :: flo_sectors
  
  integer,dimension(:),allocatable :: out_unit,out_unitV
  character(len=4) :: sect_suffix
  character(len=2) :: flo_suffix
  real(8) :: wflo,dU,check_flo
  real(8),dimension(:),allocatable :: tmp_gs_flo,gs_flo,gs_flo_new,gs_sect

  integer :: isfl,jsfl
  
  !+- PARSE INPUT DRIVER -+!
  ! call parse_input_variable(Nx,"Nx","inputGZ.conf",default=10)
  ! call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  ! call parse_input_variable(Wband,"Wband","inputGZ.conf",default=1.d0)
  ! call parse_input_variable(sweep,"SWEEP","inputGZ.conf",default='sweepJ')  
  ! call parse_input_variable(sweep_start,"SWEEP_START","inputGZ.conf",default=-0.4d0)  
  ! call parse_input_variable(sweep_stop,"SWEEP_STOP","inputGZ.conf",default=0.d0)  
  ! call parse_input_variable(sweep_step,"SWEEP_STEP","inputGZ.conf",default=0.05d0)  
  ! call parse_input_variable(read_dir,"READ_DIR","inputGZ.conf",default='~/etc_local/GZ_basis/')
  ! call parse_input_variable(store_dir,"STORE_DIR","inputGZ.conf",default='./READ_PHI_TRACES/')
  ! call parse_input_variable(nread,"NREAD","inputGZ.conf",default=0.d0)
  ! call parse_input_variable(nerr,"NERR","inputGZ.conf",default=1.d-7)
  call parse_input_variable(wflo,"Wfloquet","inputGZ.conf",default=0.5d0)
  call parse_input_variable(dU,"DU","inputGZ.conf",default=0.0d0)
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")

  if(Norb.eq.1.and.wf_symmetry.eq.1) then
     write(*,*) 'WARNING THE O(1) x SU(2)c x ORBITAL_ROTATION = O(1) x SU(2)c for the Norb=1 case!'
     wf_symmetry=0
  end if
  !
  !NOTE: ON HUNDS COUPLINGS:
  !NORB=3 RATATIONAL INVARIANT HAMILTONIAN       :: Jsfl=Jh, Jph=U-Ust-J   (NO relation between Ust and U)
  !       FULLY ROTATIONAL INVARIANT HAMILTONIAN :: Jsfl=Jh, Jph=J, Ust = U - 2J   
  Jsf = Jh
  Jph = Jh
  Ust = Uloc(1)-2.d0*Jh
  Uloc=Uloc(1)
  !
  !
  call initialize_local_fock_space_Nsector
  allocate(deltaH(nFock,nFock),tmpH(nFock,nFock))
  tmpH = local_hamiltonian
  Uloc(2) = Uloc(1)-dU
  Uloc(3) = Uloc(1)-dU
  call build_local_hamiltonian
  deltaH = -local_hamiltonian+tmpH
  !local_hamiltonian = tmpH
  !

  allocate(flo_sect(11))
  flo_sect= arange(-5,11,.true.)
  write(*,*) flo_Sect


  


  

 allocate(gs_flo(Nsectors),gs_flo_new(Nsectors),gs_sect(Nsectors)) 
 !
 gs_flo=0.d0
 gs_flo_new=0.d0
 !

 allocate(flo_sectors(Nsectors))
 
 do Nflo=1,31,2

    
    
    gs_flo_new=gs_flo
    
    do isect=1,Nsectors
       !       write(*,*) "diagonalising sector",isect
       dim_s = sectors(isect)%dim_sect
       dim=dim_s*Nflo
       flo_sectors(isect)%dim_s = dim_s
       flo_sectors(isect)%dim_flo = Nflo
       flo_sectors(isect)%dim = dim
       !
       allocate(flo_sectors(isect)%E_flo(dim))
       allocate(flo_sectors(isect)%V_flo(dim,dim))
       !
    end do

    do isect=1,Nsectors
       !
       dim_s = sectors(isect)%dim_sect
       dim=dim_s*Nflo
       !
       if(isect==1) write(*,*) dim,dim_s,Nflo
       allocate(hamiltonian_sector(dim,dim));hamiltonian_sector=0.d0
       allocate(eigenv_sector(dim))       
       do iflo=1,Nflo
          do is=1,dim_s
             do js=1,dim_s
                ifock=sectors(isect)%states_sect(is)
                jfock=sectors(isect)%states_sect(js)
                isfl = (iflo-1)*dim_s + is
                jsfl = (iflo-1)*dim_s + js
                hamiltonian_sector(isfl,jsfl) = hamiltonian_sector(isfl,jsfl) + &
                     local_hamiltonian(ifock,jfock) 
                if(is==js) then
                   hamiltonian_sector(isfl,jsfl) = hamiltonian_sector(isfl,jsfl) + dble(iflo-(Nflo+1)/2)*wflo
                end if
                !
                if(iflo.gt.1) then
                   jflo = iflo -1
                   isfl = (iflo-1)*dim_s + is
                   jsfl = (jflo-1)*dim_s + js
                   hamiltonian_sector(isfl,jsfl) = hamiltonian_sector(isfl,jsfl) +&
                        deltaH(ifock,jfock)*0.5d0 !+- "floquet hopping"
                end if
                !
                if(iflo.lt.Nflo) then
                   jflo = iflo +1
                   isfl = (iflo-1)*dim_s + is
                   jsfl = (jflo-1)*dim_s + js
                   hamiltonian_sector(isfl,jsfl) = hamiltonian_sector(isfl,jsfl) +&
                        deltaH(ifock,jfock)*0.5d0 !+- "floquet hopping"
                end if
             end do
          end do
          !
       end do
       call matrix_diagonalize(hamiltonian_sector,eigenv_sector)     
       !
       flo_sectors(isect)%E_flo = eigenv_sector
       gs_sect(isect) = eigenv_sector(1)
       
       
       deallocate(eigenv_sector,hamiltonian_sector)
       ! write(out_unitV(isect),'(100F18.10)')
    end do
    

    !+- check convergency on m=0 sector -+!    
    gs_flo = gs_sect + dble((Nflo+1)/2)*wflo    
    if(Nflo==1) then
       check_flo = 1.d0
    else
       check_flo=0.d0
       do isect=1,Nsectors
          check_flo = check_flo + (gs_flo(isect)-gs_flo_new(isect))**2.d0
       end do
    end if
    write(400,*) Nflo,check_flo
    write(401,'(10F18.10)') dble(Nflo),gs_flo
    write(402,'(10F18.10)') dble(Nflo),gs_flo_new
    write(403,'(10F18.10)') dble(Nflo),gs_sect!+dble((Nflo+1)/2)*wflo
    !
    if(check_flo.lt.1.d-20) then
       !
       allocate(out_unit(Nsectors*3))
       out_unit=free_units(Nsectors*3)       
       iout=0
       do isect=1,Nsectors
          do iflo=(Nflo+1)/2-1,(Nflo+1)/2+1
             iout=iout+1
             write(sect_suffix,'(I4.4)') isect-1          
             sect_suffix=trim(adjustr(sect_suffix))
             write(flo_suffix,'(I2.2)') iflo-(Nflo+1)/2+2
             flo_suffix=trim(adjustr(flo_suffix))
             inquire(file='spectrum_'//sect_suffix//"_floq_"//flo_suffix//".data", exist=exist)
             if (exist) then
                open(out_unit(iout),file='spectrum_'//sect_suffix//"_floq_"//flo_suffix//".data",status="old",position="append",action="write")
             else
                open(out_unit(iout),file='spectrum_'//sect_suffix//"_floq_"//flo_suffix//".data",status="new",action="write")
             end if
          end do
       end do
       !
       iout=0
       do isect=1,Nsectors
          dim_s = sectors(isect)%dim_sect
          dim=dim_s*Nflo
          do iflo=(Nflo+1)/2-1,(Nflo+1)/2+1
             iout=iout+1
             isfl = (iflo-1)*dim_s + 1
             jsfl = (iflo-1)*dim_s + dim_s
             write(out_unit(iout),'(400F18.10)') Uloc(1)-Uloc(2),flo_sectors(isect)%E_flo(isfl:jsfl)
          end do
       end do
       exit

    end if
    !
    do isect=1,Nsectors
       deallocate(flo_sectors(isect)%E_flo,flo_sectors(isect)%V_flo)
    end do
    !
 end do
 
 
 write(*,*) 'floquet_converged',Nflo 
 do isect=1,Nsectors
    close(out_unit(isect))
 end do
 
 
 
 stop





 ! out_unitV=free_units(Nsectors)
 ! do isect=1,Nsectors
 !    write(sect_suffix,'(I4.4)') isect-1
 !    sect_suffix=trim(adjustr(sect_suffix))
    
 !    inquire(file='Evect_'//sect_suffix//".data", exist=exist)
 !    if (exist) then
 !       open(out_unitV(isect),file='Evect_'//sect_suffix//".data",status="old",position="append",action="write")
 !    else
 !       open(out_unitV(isect),file='Evect_'//sect_suffix//".data",status="new",action="write")
 !    end if
 ! end do

 

 do isect=1,Nsectors
    write(*,*) "diagonalising sector",isect
    dim_s = sectors(isect)%dim_sect
    dim=dim_s*Nflo
    if(isect==1) write(*,*) dim,dim_s,Nflo
    allocate(hamiltonian_sector(dim,dim));hamiltonian_sector=0.d0
    allocate(eigenv_sector(dim))       
    do iflo=1,Nflo
       do is=1,dim_s
          do js=1,dim_s
             ifock=sectors(isect)%states_sect(is)
             jfock=sectors(isect)%states_sect(js)
             isfl = (iflo-1)*dim_s + is
             jsfl = (iflo-1)*dim_s + js
             hamiltonian_sector(isfl,jsfl) = hamiltonian_sector(isfl,jsfl) + &
                  local_hamiltonian(ifock,jfock) 
!             if(ifock==jfock) then
                hamiltonian_sector(isfl,jsfl) = hamiltonian_sector(isfl,jsfl) + dble(iflo-(Nflo+1)/2)*wflo
!             end if
             
             if(iflo.gt.1) then
                jflo = iflo -1
                isfl = (iflo-1)*dim_s + is
                jsfl = (jflo-1)*dim_s + js
                hamiltonian_sector(isfl,jsfl) = hamiltonian_sector(isfl,jsfl) +&
                     deltaH(ifock,jfock)*0.5d0 !+- "floquet hopping"
             end if
             !
             if(iflo.lt.Nflo) then
                jflo = iflo +1
                isfl = (iflo-1)*dim_s + is
                jsfl = (jflo-1)*dim_s + js
                hamiltonian_sector(isfl,jsfl) = hamiltonian_sector(isfl,jsfl) +&
                     deltaH(ifock,jfock)*0.5d0 !+- "floquet hopping"
             end if
          end do
       end do
       !
    end do
    call matrix_diagonalize(hamiltonian_sector,eigenv_sector)     
    !
    !
    do ieigen=1,dim_s
       allocate(tmp_m0(dim_s),tmp_mm(dim_s),tmp_mp(dim_s))
       iflo=(Nflo+1)/2 
       isfl = (iflo-1)*dim_s + 1
       jsfl = (iflo-1)*dim_s + dim_s
       tmp_m0(:) = hamiltonian_sector(isfl:jsfl,ieigen)
       tmp_mm=0.d0
       tmp_mp=0.d0
       if(Nflo.gt.1) then
          isfl = (iflo-2)*dim_s + 1
          jsfl = (iflo-2)*dim_s + dim_s
          tmp_mm = hamiltonian_sector(isfl:jsfl,ieigen)
          isfl = (iflo)*dim_s + 1
          jsfl = (iflo)*dim_s + dim_s
          tmp_mp = hamiltonian_sector(isfl:jsfl,ieigen)
       end if

       sectors(isect)%E_sect(ieigen) = 0.d0
       do is=1,dim_s
          do js=1,dim_s
             ifock=sectors(isect)%states_sect(is)
             jfock=sectors(isect)%states_sect(js)
!             if(ifock==jfock) then
                sectors(isect)%E_sect(ieigen) = sectors(isect)%E_sect(ieigen) + dble(iflo-(Nflo+1)/2)*wflo*tmp_m0(is)*tmp_m0(js)
!             end if
             sectors(isect)%E_sect(ieigen) = sectors(isect)%E_sect(ieigen) + local_hamiltonian(ifock,jfock)*tmp_m0(is)*tmp_m0(js)
             sectors(isect)%E_sect(ieigen) = sectors(isect)%E_sect(ieigen) + deltaH(ifock,jfock)*0.5d0*tmp_m0(is)*tmp_mm(js)
             sectors(isect)%E_sect(ieigen) = sectors(isect)%E_sect(ieigen) + deltaH(ifock,jfock)*0.5d0*tmp_m0(is)*tmp_mp(js)
          end do
       end do
       deallocate(tmp_m0,tmp_mm,tmp_mp)
    
    end do
    deallocate(eigenv_sector,hamiltonian_sector)
    !
    write(out_unit(isect),'(400F18.10)') sectors(isect)%E_sect(1:dim_s)
    !
 end do
 

 
 
 stop

 
  !
  stop
  

  ! allocate(tmp_states(nFock),ivec(Ns))
  ! count_states=0
  ! do ifock=1,NFock
  !    call bdecomp(ifock,ivec)
  !    nstate=sum(ivec)
  !    if(nstate.eq.Norb) then
  !       itest=0
  !       do iorb=1,Norb
  !          itest = itest + ivec(iorb)*ivec(iorb+Norb)
  !       end do
  !       if(itest.ne.0) then
  !          count_states = count_states + 1
  !          tmp_states(count_states) = ifock
  !       end if
  !    end if
  ! end do
  ! allocate(states210(count_states))
  ! states210=tmp_states(1:count_states)

  ! open(out_unit,file='210_states.info')
  ! do i=1,count_states
  !    call bdecomp(states210(i),ivec)
  !    write(out_unit,*) ivec,'   ',states210(i)
  ! end do
  ! close(out_unit)


  
  Jsf = Jh
  Jph = Jh
  Ust = Uloc(1)-2.d0*Jh
  call build_local_hamiltonian

!  allocate(local_eigenv(nFock))

  !
  !  open(out_unit,file='atomic_levels.info')
  !  write(out_unit,'(30F18.10)') Uloc(1),local_eigenv(1:20)
  
 ! stop



  !
  !
CONTAINS
  !

  subroutine initialize_local_fock_space_Nsector
    integer :: iorb,jorb,ispin,jspin,ifock,jfock,count,isect,Nstates,dim,is
    integer,dimension(:),allocatable   :: Fock,ivec,istates
    integer :: unit
    !
    !     type sector
    !    private
    !    integer :: dim_sect
    !    integer,dimension(:),allocatable :: states_sect
    !    integer,dimension(:),allocatable :: map_sect2fock     
    ! end type sector
    
    Ns = 2*Norb        
    NFock = 2**Ns
    Nsectors = Ns+1
    !
    allocate(ivec(Ns))
    allocate(Fock(Nfock))
    allocate(sectors(Nsectors),istates(Nsectors))
    do isect=1,Nsectors
       sectors(isect)%dim_sect=0
    end do

    unit=free_unit()
    open(unit,file="states_fock.info")
    !+- find the dimensions of each sector
    do ifock=1,NFock
       Fock(ifock)=ifock
       call bdecomp(Fock(ifock),ivec)
       call get_state_number(ivec,jfock)
       write(unit,'(10I3)') ivec(:),ifock,jfock
       count=sum(ivec(1:Ns))
       sectors(count+1)%dim_sect = sectors(count+1)%dim_sect + 1       
    end do

    do isect=1,Nsectors
       dim = sectors(isect)%dim_sect
       allocate(sectors(isect)%states_sect(dim))
       allocate(sectors(isect)%E_sect(dim))
       allocate(sectors(isect)%V_sect(dim,dim))
    end do
    istates=0
    do ifock=1,NFock
       Fock(ifock)=ifock
       call bdecomp(Fock(ifock),ivec)
       !call get_state_number(ivec,ifock)
       !       
       count=sum(ivec(1:Ns))
       istates(count+1)=istates(count+1)+1
       sectors(count+1)%states_sect(istates(count+1)) = ifock
    end do

    do isect=1,Nsectors
       write(unit,*) "number of particles",isect-1,"dimension",sectors(isect)%dim_sect
       Nstates=sectors(isect)%dim_sect
       do is=1,Nstates
          ifock = sectors(isect)%states_sect(is)
          call bdecomp(ifock,ivec)
          write(unit,*) ivec(:),'|',ifock,'>'
       end do
    end do
    !
    !+- Allocate and initialize stride -+! 
    allocate(index(2,Norb))
    do ispin=1,2
       do iorb=1,Norb
          index(ispin,iorb)=iorb+(ispin-1)*Norb
       enddo
    end do
    !
    call build_local_fock_algebra
    !
    call build_local_hamiltonian
    !
    call build_local_observables
    !
    ! TO BE REMOVED ONCE SU(2) AND GENERAL ROTATIONS SYMMETRIES ARE DIRECTELY IMPLEMENTED
    !call get_spin_indep_states
    !
  end subroutine initialize_local_fock_space_Nsector










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


