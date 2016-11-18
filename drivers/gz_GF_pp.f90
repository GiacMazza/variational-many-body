program GZ_GF
  USE SCIFOR
  !
  USE DMFT_MISC
  USE SF_PARSE_INPUT
  !
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_neqGREENS_FUNCTIONS
  implicit none
  

  integer :: Nx,Nw
  real(8) :: Wrange,dw
  character(len=200) :: read_neq_dir
  real(8),dimension(:),allocatable :: epsik,hybik
  complex(8),dimension(:,:,:),allocatable :: neq_Rhop
  complex(8),dimension(:,:,:,:),allocatable :: Gloc_ret_tt
  complex(8),dimension(:,:,:,:),allocatable :: Gloc_ret_tw,Gloc_ret_tw_
  complex(8),dimension(:,:),allocatable :: tmpG
  real(8),dimension(:),allocatable :: wre  
  integer :: it,jt,unit,iw,iti,jtj
  integer :: is
  
  call parse_input_variable(Wband,"WBAND","inputGZgz.conf",default=2.d0)
  call parse_input_variable(Nx,"Nx","inputGZgz.conf",default=1000)
  call parse_input_variable(read_neq_dir,"READ_NEQ_DIR","inputGZgz.conf",default='./')
  call parse_input_variable(Norb,"Norb","inputGZgz.conf",default=1)
  call parse_input_variable(Ntgf,"NTGF","inputGZgz.conf",default=100)
  call parse_input_variable(Nw,"Nw","inputGZgz.conf",default=200)
  call parse_input_variable(Wrange,"WRANGE","inputGZgz.conf",default=2.d0)
  call parse_input_variable(tstep,"TSTEP","inputGZgz.conf",default=1.d-2)
  call save_input_file("inputGZgz.conf")

  call initialize_local_fock_space    
  call build_lattice_model; get_Hk_t => getHk
  
  allocate(neq_rhop(2*Ntgf-1,Ns,Ns),t_grid(2*Ntgf-1))  
  call read_neq_Rhop(read_neq_dir,neq_Rhop)

  allocate(Gloc_ret_tt(Ntgf,Ntgf,Ns,Ns)) 
  call gz_get_Gloc_ret(neq_Rhop,Gloc_ret_tt)
  
  open(unit,file='Gloc_ret_tt.data')
  do it=1,Ntgf
     do jt=1,Ntgf
        iti = it + Ntgf - 1
        jtj = iti + 1 - jt
        write(unit,'(6F18.10)') t_grid(iti),t_grid(jtj),Gloc_ret_tt(it,jt,1,1)        
     end do
     write(unit,'(6F18.10)')
  end do
  close(unit)
  allocate(wre(Nw))
  wre=linspace(-wrange,wrange,Nw,mesh=dw)
  allocate(Gloc_ret_tw(Ntgf,Nw,Ns,Ns)); Gloc_ret_tw=zero
  allocate(Gloc_ret_tw_(Ntgf,Nw,Ns,Ns)); Gloc_ret_tw_=zero
  do is=1,Ns
     call get_relative_time_FT(Gloc_ret_tt(:,:,is,is),Gloc_ret_tw(:,:,is,is),wre)
  end do

  allocate(tmpG(Ns,Ns))
  Gloc_ret_tw_(1,:,:,:) = Gloc_ret_tw(1,:,:,:)

  do iw=1,Nw
     tmpG=zero
     do it=1,Ntgf-1
        tmpG=tmpG+Gloc_ret_tw(it,iw,:,:)*tstep*0.5d0
        tmpG=tmpG+Gloc_ret_tw(it+1,iw,:,:)*tstep*0.5d0
        iti = it + Ntgf -1
        Gloc_ret_tw_(it+1,iw,:,:) = tmpG/(t_grid(iti+1)-t_grid(Ntgf))
     end do
  end do
  
  
  open(unit,file='Gloc_ret_tw.data')
  do it=1,Ntgf
     do iw=1,Nw
        write(unit,'(6F18.10)') t_grid(it+Ntgf-1),wre(iw),Gloc_ret_tw(it,iw,1,1),Gloc_ret_tw_(it,iw,1,1)          
     end do
     write(unit,'(6F18.10)')
  end do
  close(unit)

  open(unit,file='Gloc_ret_tw_.data')
  do it=1,Ntgf
     do iw=1,Nw
        write(unit,'(10F18.10)') t_grid(it+Ntgf-1),wre(iw),Gloc_ret_tw(it,iw,1,1),Gloc_ret_tw_(it,iw,1,1)        
     end do
     write(unit,'(6F18.10)')
     write(unit,'(6F18.10)')
  end do
  close(unit)

contains
  
  subroutine read_neq_Rhop(read_neq_dir,neq_Rhop)
    implicit none
    character(len=200) :: read_neq_dir,file_name

    integer,dimension(:),allocatable :: ISr,JSr
    integer,dimension(:,:),allocatable :: grid_Rhop
    integer :: ir,i,ifile,it
    integer :: ios,unit
    logical :: read_check
    integer :: tmp_ir    
    real(8),dimension(:),allocatable :: dump_vect
    complex(8),dimension(2*Ntgf-1,Ns,Ns) :: neq_Rhop
    !
    read_neq_dir=trim(read_neq_dir)
    file_name=reg(read_neq_dir)//"grid_Rhop.info"
    inquire(file=file_name,exist=read_check)
    if(read_check) then     
       unit=free_unit()
       open(unit,file=file_name,status='old')
       ir=0
       do
          read (unit,*,iostat=ios) tmp_ir
          if (ios/=0) exit
          ir = ir + 1
       end do
       close(unit)
       open(unit,file=file_name,status='old')
       allocate(ISr(ir),JSr(ir))
       allocate(grid_Rhop(Ns,Ns));grid_Rhop=0
       do ifile=1,ir
          read(unit,*) tmp_ir,ISr(ifile),JSr(ifile)
          grid_Rhop(ISr(ifile),JSr(ifile))=1
       end do
       close(unit)
    else
       write(*,*) "file ",file_name,"  not present!"
       stop
    end if




    file_name=reg(read_neq_dir)//"neq_Rhop_matrix.data"


    inquire(file=file_name,exist=read_check)
    if(read_check) then     
       unit=free_unit()
       open(unit,file=file_name,status='old')
       !
       ios=0
       allocate(dump_vect(2*ir))     
       do it=1,Ntgf
          neq_Rhop(it+Ntgf,:,:) = zero
          read(unit,*,iostat=ios) t_grid(it+Ntgf-1),dump_vect(1:2*ir)
          !if(ios/=0.and.it.ne.2*Ntgf) stop "Ntgf > Nt/2"
          do i=1,ir
             neq_Rhop(it+Ntgf-1,ISr(i),JSr(i)) = dump_vect(i) + xi*dump_vect(i+ir)
          end do
       end do
       do it=1,Ntgf-1
          neq_Rhop(it,:,:) = neq_Rhop(Ntgf,:,:)
          t_grid(it) = -t_grid(2*Ntgf-it)
       end do
       close(unit)
    else
       write(*,*) "file ",file_name,"  not present!"
       stop
    end if
    !
    open(unit,file='neq_Rhop_read.data')
    do it=1,2*Ntgf-1
       !
       call write_complex_matrix_grid(neq_Rhop(it,:,:),unit,grid_Rhop,t_grid(it))
       !
    end do
    close(unit)


  end subroutine read_neq_Rhop


  

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
       if(epsik(ix).gt.-Wband/2.d0.and.epsik(ix).lt.Wband/2) then
          wtk(ix)=4.d0/Wband/pi*sqrt(1.d0-(2.d0*epsik(ix)/Wband)**2.d0)*de
       else
          wtk(ix) = 0.d0
       end if
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
    complex(8),dimension(:,:) :: Hk
    integer                   :: ik
    real(8)                   :: time
    if(size(Hk,1).ne.size(Hk,2)) stop "wrong dimenions in getHk"
    if(size(Hk,1).ne.Ns) stop "wrong dimenions in getHk"
    Hk = Hk_tb(:,:,ik)
  end subroutine getHk





  

end program GZ_GF
