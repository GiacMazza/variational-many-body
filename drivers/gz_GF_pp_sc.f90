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
  complex(8),dimension(:,:,:),allocatable :: neq_Rhop,neq_Qhop,neq_vdm,neq_dens
  complex(8),dimension(:,:),allocatable :: VSC_op
  complex(8),dimension(:,:,:,:),allocatable :: Gloc_ret_tt,Gloc_ret_tt_
  complex(8),dimension(:,:,:,:),allocatable :: Gloc_ret_tw,Gloc_ret_tw_
  complex(8),dimension(:,:),allocatable :: tmpG,Gloc_ret_tot_tt,Gloc_ret_tot_tw
  real(8),dimension(:),allocatable :: wre
  real(8),dimension(:,:),allocatable :: orig2nat_angle
  real(8) :: gap0
  integer :: it,jt,unit,iw,iti,jtj,uwrite
  integer :: is,iorb
  logical :: read_gloc
  real(8) :: ti,tf,dt_tmp
  complex(8),dimension(:),allocatable :: dumpGloc,dumpGloc_

  
  call parse_input_variable(Wband,"WBAND","inputGZgz.conf",default=2.d0)
  call parse_input_variable(Nx,"Nx","inputGZgz.conf",default=1000)
  call parse_input_variable(read_neq_dir,"READ_NEQ_DIR","inputGZgz.conf",default='./')
  call parse_input_variable(Norb,"Norb","inputGZgz.conf",default=1)
  call parse_input_variable(Ntgf,"NTGF","inputGZgz.conf",default=100)
  call parse_input_variable(Nt0,"NT0","inputGZgz.conf",default=100)
  call parse_input_variable(Nw,"Nw","inputGZgz.conf",default=200)
  call parse_input_variable(Wrange,"WRANGE","inputGZgz.conf",default=2.d0)
  call parse_input_variable(tstep,"TSTEP","inputGZgz.conf",default=1.d-2)
  call parse_input_variable(gap0,"gap0","inputGZgz.conf",default=0.d0)
  call parse_input_variable(read_gloc,"READ_GLOC","inputGZgz.conf",default=.false.)
  call save_input_file("inputGZgz.conf")

  call initialize_local_fock_space    
  call build_lattice_model; get_Hk_t => getHk

  Nttgf = Nt0+Ntgf-1
  
  allocate(neq_rhop(Nttgf,Ns,Ns),neq_qhop(Nttgf,Ns,Ns),t_grid(Nttgf))
  ti=-tstep*Nt0; tf=tstep*Ntgf
  t_grid=linspace(ti,tf,Nttgf,istart=.false.,iend=.false.,mesh=dt_tmp)
  if(abs(dt_tmp-tstep).gt.1.d-8) then
     write(*,*) dt_tmp,tstep,ti,tf
     stop "wrong time grid"
  end if
  write(*,*) t_grid(1),t_grid(Nt0),t_grid(Nttgf)
  call read_neq_Rhop_Qhop(read_neq_dir,neq_Rhop,neq_Qhop)
  !
  ! do it=1,Nttgf
  !    neq_Rhop(it,:,:) = neq_Rhop(Nt0,:,:)*(cos(gap0*t_grid(it))-xi*sin(gap0*t_grid(it)))
  !    neq_Qhop(it,:,:) = neq_Qhop(Nt0,:,:)*(cos(gap0*t_grid(it))+xi*sin(gap0*t_grid(it)))
  ! end do
  !
  allocate(Gloc_ret_tt(Ntgf,Ntgf,Ns,Ns),Gloc_ret_tt_(Ntgf,Ntgf,2*Ns,2*Ns)) 
  if(read_gloc) then
     call read_gloc_tt(read_neq_dir,Gloc_ret_tt)             
  else  
     call gz_get_Gloc_ret_superc(neq_Rhop,neq_Qhop,Gloc_ret_tt_)
     Gloc_ret_tt=Gloc_ret_tt_(:,:,1:Ns,1:Ns) 
  end if
  !
  !
  allocate(dumpGloc(Ns),dumpGloc_(Ns))
  !
  !
  unit=free_unit()
  open(unit,file='Gloc_ret_tt.data')
  do it=1,Ntgf
     do jt=1,Ntgf
        iti = it + Nt0 - 1
        jtj = iti + 1 - jt
        do is=1,Ns
           dumpGloc(is)=Gloc_ret_tt(it,jt,is,is)
        end do
        write(unit,'(20F18.10)') t_grid(iti),t_grid(jtj),dumpGloc(1:Ns)
     end do
     write(unit,'(6F18.10)')
  end do
  close(unit)



  open(unit,file='Gloc_ret_tt_.data')
  do it=1,Ntgf
     do jt=1,Ntgf
        iti = it + Nt0 - 1
        jtj = iti + 1 - jt
        do is=1,Ns
           dumpGloc(is)=Gloc_ret_tt(it,jt,is,is)
        end do
        write(unit,'(20F18.10)') t_grid(iti),t_grid(jtj),dumpGloc(1:Ns)
     end do
     write(unit,'(6F18.10)')
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
  !
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
        do is=1,Ns
           dumpGloc(is)=Gloc_ret_tw(it,iw,is,is)
           dumpGloc_(is)=Gloc_ret_tw_(it,iw,is,is)
        end do
        write(unit,'(30F18.10)') t_grid(it+Nt0-1),wre(iw),dumpGloc(1:Ns),dumpGloc_(1:Ns)
!        write(unit,'(20F18.10)') t_grid(it+Nt0-1),wre(iw),Gloc_ret_tw(it,iw,1,1),Gloc_ret_tw_(it,iw,1,1),Gloc_ret_tw(it,iw,2,2),Gloc_ret_tw_(it,iw,2,2)          
     end do
     write(unit,'(20F18.10)')
  end do
  close(unit)

  open(unit,file='Gloc_ret_tw_.data')
  do it=1,Ntgf
     do iw=1,Nw
        do is=1,Ns
           dumpGloc(is)=Gloc_ret_tw(it,iw,is,is)
           dumpGloc_(is)=Gloc_ret_tw_(it,iw,is,is)
        end do
        write(unit,'(30F18.10)') t_grid(it+Nt0-1),wre(iw),dumpGloc(1:Ns),dumpGloc_(1:Ns)
     end do
     write(unit,'(6F18.10)')
     write(unit,'(6F18.10)')
  end do
  close(unit)

contains

  subroutine read_neq_VDM(read_neq_dir,neq_vdm)
    implicit none
    character(len=200) :: read_neq_dir,file_name
    complex(8),dimension(2*Ntgf-1,Ns,Ns) :: neq_vdm
    integer :: ir,i,ifile,it,iq,is
    integer :: ios,unit
    logical :: read_check
    real(8),dimension(:),allocatable :: dump_vect
    !
    read_neq_dir=trim(read_neq_dir)   
    file_name=reg(read_neq_dir)//"neq_dens_constrGZ.data"
    inquire(file=file_name,exist=read_check)
    if(read_check) then     
       unit=free_unit()
       open(unit,file=file_name,status='old')
       !
       ios=0
       allocate(dump_vect(Ns))     
       do it=1,Ntgf
          neq_vdm(it+Ntgf-1,:,:) = zero
          read(unit,*,iostat=ios) t_grid(it+Ntgf-1),dump_vect(1:Ns)
          do is=1,Ns
             neq_vdm(it+Ntgf-1,is,is) = dump_vect(is) 
          end do
       end do
       do it=1,Ntgf-1
          t_grid(it) = -t_grid(2*Ntgf-it)
          neq_vdm(it,:,:) = neq_vdm(Ntgf,:,:)
       end do
       close(unit)
       deallocate(dump_vect)
    else
       write(*,*) "file ",file_name,"  not present!"
       stop
    end if
    !
    open(unit,file='neq_vdm_read.data')
    do it=1,2*Ntgf-1
       call write_hermitean_matrix(neq_vdm(it,:,:),unit,t_grid(it))
    end do
    close(unit)
    !
  end subroutine read_neq_VDM




  subroutine read_gloc_tt(read_neq_dir,neq_Gloc)
    implicit none
    character(len=200) :: read_neq_dir,file_name
    complex(8),dimension(Ntgf,Ntgf,Ns,Ns) :: neq_Gloc
    integer :: ir,i,ifile,it,iq,is
    integer :: ios,unit
    logical :: read_check
    real(8),dimension(:),allocatable :: dump_vect
    real(8) :: t,tt
    !
    write(*,*) "reading GLOC"
    read_neq_dir=trim(read_neq_dir)   
    file_name=reg(read_neq_dir)//"Gloc_ret_tt.data"
    inquire(file=file_name,exist=read_check)
    if(read_check) then     
       unit=free_unit()
       open(unit,file=file_name,status='old')
       !
       ios=0
       allocate(dump_vect(2*Ns))     
       do it=1,Ntgf
          do jt=1,Ntgf              
             read(unit,*,iostat=ios) t,tt,dump_vect(1:2*Ns)
             do is=1,Ns
                neq_gloc(it,jt,is,is) = dump_vect(2*is-1)+xi*dump_vect(2*is)
             end do
             ! neq_gloc(it,jt,2,2) = dump_vect(3)+xi*dump_vect(4)
             ! neq_gloc(it,jt,3,3) = dump_vect(3)+xi*dump_vect(4)
             ! neq_gloc(it,jt,4,4) = dump_vect(1)+xi*dump_vect(2)
             ! neq_gloc(it,jt,5,5) = dump_vect(3)+xi*dump_vect(4)
             ! neq_gloc(it,jt,6,6) = dump_vect(3)+xi*dump_vect(4)
          end do
       end do
       close(unit)
       deallocate(dump_vect)
    else
       write(*,*) "file ",file_name,"  not present!"
       stop
    end if    
    !
  end subroutine read_gloc_tt



  subroutine read_neq_dens(read_neq_dir,neq_dens)
    implicit none
    character(len=200) :: read_neq_dir,file_name
    complex(8),dimension(2*Ntgf-1,Ns,Ns) :: neq_dens
    integer :: ir,i,ifile,it,iq,is
    integer :: ios,unit
    logical :: read_check
    real(8),dimension(:),allocatable :: dump_vect
    !
    read_neq_dir=trim(read_neq_dir)   
    file_name=reg(read_neq_dir)//"neq_local_density_matrix.data"
    inquire(file=file_name,exist=read_check)
    if(read_check) then     
       unit=free_unit()
       open(unit,file=file_name,status='old')
       !
       ios=0
       allocate(dump_vect(Ns))     
       do it=1,Ntgf
          neq_dens(it+Ntgf-1,:,:) = zero
          read(unit,*,iostat=ios) t_grid(it+Ntgf-1),dump_vect(1:Ns)
          do is=1,Ns
             neq_dens(it+Ntgf-1,is,is) = dump_vect(is) 
          end do
       end do
       do it=1,Ntgf-1
          t_grid(it) = -t_grid(2*Ntgf-it)
          neq_dens(it,:,:) = neq_dens(Ntgf,:,:)
       end do
       close(unit)
       deallocate(dump_vect)
    else
       write(*,*) "file ",file_name,"  not present!"
       stop
    end if
    !
    open(unit,file='neq_density_read.data')
    do it=1,2*Ntgf-1
       call write_hermitean_matrix(neq_dens(it,:,:),unit,t_grid(it))
    end do
    close(unit)
    !
  end subroutine read_neq_dens

  

  
  subroutine read_neq_Rhop_Qhop(read_neq_dir,neq_Rhop,neq_Qhop)
    implicit none
    character(len=200) :: read_neq_dir,file_name

    integer,dimension(:),allocatable :: ISr,JSr,ISq,JSq
    integer,dimension(:,:),allocatable :: grid_Rhop,grid_Qhop
    integer :: ir,i,ifile,it,iq
    integer :: ios,unit
    logical :: read_check
    integer :: tmp_ir,tmp_iq
    real(8) :: tmpt
    real(8),dimension(:),allocatable :: dump_vect
    complex(8),dimension(Nttgf,Ns,Ns) :: neq_Rhop,neq_Qhop
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
       neq_Rhop = zero                 
       do it=1,Ntgf
          read(unit,*,iostat=ios) tmpt,dump_vect(1:2*ir) !HERE
          do i=1,ir
             neq_Rhop(it+Nt0-1,ISr(i),JSr(i)) = dump_vect(i) + xi*dump_vect(i+ir)
          end do
          if(abs(tmpt-t_grid(it+Nt0-1)).gt.1.d-8) stop "file to read and given tstep not compatibles"
       end do
       do it=1,Nt0-1
          !          t_grid(it) = -t_grid(2*Ntgf-it)
          neq_Rhop(it,:,:) = neq_Rhop(Nt0,:,:)*cos(gap0*t_grid(it))
          neq_Rhop(it,:,:) = neq_Rhop(it,:,:)-neq_Rhop(Nt0,:,:)*xi*sin(gap0*t_grid(it))
          !
          ! neq_Rhop(it,:,:) = neq_Rhop(Ntgf,:,:)
          ! t_grid(it) = -t_grid(2*Ntgf-it)
       end do

       !<TMP
       ! do it=1,2*Ntgf-1
       !    neq_Rhop(it,:,:) = neq_Rhop(Ntgf,:,:)*(cos(gap0*t_grid(it))-xi*sin(gap0*t_grid(it)))
       !    neq_Qhop(it,:,:) = neq_Qhop(Ntgf,:,:)*(cos(gap0*t_grid(it))-xi*sin(gap0*t_grid(it)))
       ! end do       
       !TMP>
       
       close(unit)
       deallocate(dump_vect)
    else
       write(*,*) "file ",file_name,"  not present!"
       stop
    end if
    !
    open(unit,file='neq_Rhop_read.data')
    do it=1,Nttgf
       !
       call write_complex_matrix_grid(neq_Rhop(it,:,:),unit,grid_Rhop,t_grid(it))
       !
    end do
    close(unit)



    file_name=reg(read_neq_dir)//"grid_Qhop.info"
    inquire(file=file_name,exist=read_check)
    if(read_check) then     
       unit=free_unit()
       open(unit,file=file_name,status='old')
       iq=0
       do
          read (unit,*,iostat=ios) tmp_iq
          if (ios/=0) exit
          iq = iq + 1
       end do
       close(unit)
       open(unit,file=file_name,status='old')
       allocate(ISq(iq),JSq(iq))
       allocate(grid_Qhop(Ns,Ns));grid_Qhop=0
       do ifile=1,iq
          read(unit,*) tmp_iq,ISq(ifile),JSq(ifile)
          grid_Qhop(ISq(ifile),JSq(ifile))=1
       end do
       close(unit)
    else
       write(*,*) "file ",file_name,"  not present!"
       stop
    end if
    
    file_name=reg(read_neq_dir)//"neq_Qhop_matrix.data"
    inquire(file=file_name,exist=read_check)
    if(read_check) then     
       unit=free_unit()
       open(unit,file=file_name,status='old')
       !
       ios=0
       allocate(dump_vect(2*iq))
       neq_Qhop = zero
       do it=1,Ntgf
          read(unit,*,iostat=ios) tmpt,dump_vect(1:2*iq)
          if(abs(tmpt-t_grid(it+Nt0-1)).gt.1.d-8) stop "file to read and given tstep not compatibles"
          do i=1,iq
             neq_Qhop(it+Nt0-1,ISq(i),JSq(i)) = dump_vect(i) + xi*dump_vect(i+iq)
          end do
       end do
       do it=1,Nt0-1
          !t_grid(it) = -t_grid(2*Ntgf-it)
          neq_Qhop(it,:,:) = neq_Qhop(Nt0,:,:)*cos(gap0*t_grid(it))
          neq_Qhop(it,:,:) = neq_Qhop(it,:,:)+neq_Qhop(Nt0,:,:)*xi*sin(gap0*t_grid(it))
       end do
       close(unit)
    else
       write(*,*) "file ",file_name,"  not present!"
       stop
    end if
    !
    open(unit,file='neq_Qhop_read.data')
    do it=1,Nttgf
       !
       call write_complex_matrix_grid(neq_Qhop(it,:,:),unit,grid_Qhop,t_grid(it))
       !
    end do
    close(unit)
    deallocate(dump_vect)


    
  end subroutine read_neq_Rhop_Qhop


  

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
       !write(77,*) epsik(ix),wtk(ix)
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
