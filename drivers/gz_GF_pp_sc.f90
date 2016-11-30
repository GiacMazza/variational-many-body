program GZ_GF
  USE SCIFOR
  !
  USE DMFT_MISC
  USE SF_PARSE_INPUT
  !
  USE GZ_AUX_FUNX
  USE GZ_neqAUX_FUNX
  USE GZ_DYNAMICS
  USE RK_IDE
  !
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_neqGREENS_FUNCTIONS

  USE MPI

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
  logical :: read_gloc,no_dynamics,add_lgrA  
  real(8) :: ti,tf,dt_tmp
  complex(8),dimension(:),allocatable :: dumpGloc,dumpGloc_
  complex(8),dimension(:,:,:),allocatable :: neq_lgrA
  real(8) :: deps

  !+- START MPI -+!
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)




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
  call parse_input_variable(no_dynamics,"NO_DYN","inputGZgz.conf",default=.false.)
  call parse_input_variable(add_lgrA,"ADD_LGR","inputGZgz.conf",default=.false.)
  call parse_input_variable(deps,"DEPS","inputGZgz.conf",default=0.01d0)

  if(mpiID==0) call save_input_file("inputGZgz.conf")

  call initialize_local_fock_space    
  call build_lattice_model; get_Hk_t => getHk

  Nttgf = Nt0+Ntgf-1

  allocate(neq_rhop(Nttgf,Ns,Ns),neq_qhop(Nttgf,Ns,Ns),neq_lgrA(Nttgf,Ns,Ns),t_grid(Nttgf))
  ti=-tstep*Nt0; tf=tstep*Ntgf
  t_grid=linspace(ti,tf,Nttgf,istart=.false.,iend=.false.,mesh=dt_tmp)
  if(abs(dt_tmp-tstep).gt.1.d-8) then
     write(*,*) dt_tmp,tstep,ti,tf
     stop "wrong time grid"
  end if
  write(*,*) t_grid(1),t_grid(Nt0),t_grid(Nttgf)
  !
  if(mpiID==0) call read_neq_Rhop_Qhop(read_neq_dir,neq_Rhop,neq_Qhop)
  call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)  
  call MPI_BCAST(neq_Rhop,Nttgf*Ns*Ns,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIerr)
  call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)  
  call MPI_BCAST(neq_Qhop,Nttgf*Ns*Ns,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIerr)
  call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
  
  if(mpiID==0) call read_neq_lgr(read_neq_dir,neq_lgrA)
  call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)  
  call MPI_BCAST(neq_lgrA,Nttgf*Ns*Ns,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,MPIerr)
  call MPI_BARRIER(MPI_COMM_WORLD,MPIerr)
  !
  !
  !
  if(no_dynamics) then
     do it=1,Nttgf
        neq_Rhop(it,:,:) = neq_Rhop(Nt0,:,:)*(cos(gap0*t_grid(it))-xi*sin(gap0*t_grid(it)))
        neq_Qhop(it,:,:) = neq_Qhop(Nt0,:,:)*(cos(gap0*t_grid(it))+xi*sin(gap0*t_grid(it)))
        neq_lgrA(it,:,:) = zero!neq_lgrA(Nt0,:,:)
     end do
  end if
  !
  !
  !
  allocate(Gloc_ret_tt(Ntgf,Ntgf,Ns,Ns),Gloc_ret_tt_(Ntgf,Ntgf,2*Ns,2*Ns)) 
  if(read_gloc) then
     call read_gloc_tt(read_neq_dir,Gloc_ret_tt)             
  else
     if(add_lgrA) then
        !call gz_get_Gloc_ret_superc_mpi(neq_Rhop,neq_Qhop,Gloc_ret_tt_,neq_lgrA)
        !call gz_get_Gloc_ret_superc_mpi(neq_Rhop,neq_Qhop,Gloc_ret_tt)
        call gz_get_Gloc_ret_superc_mpik(neq_Rhop,neq_Qhop,Gloc_ret_tt,neq_lgrA)
     else
        call gz_get_Gloc_ret_superc_mpik(neq_Rhop,neq_Qhop,Gloc_ret_tt)
        !call gz_get_Gloc_ret_superc_diag_hk(neq_Rhop,neq_Qhop,Gloc_ret_tt_)
     end if
     !Gloc_ret_tt=Gloc_ret_tt_(:,:,1:Ns,1:Ns) 
  end if
  !
  !
  !
  allocate(dumpGloc(Ns),dumpGloc_(Ns))
  !
  !
  if(mpiID==0) then
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
     !
     open(unit,file='Gloc_ret_tt_.data')
     do it=1,Ntgf
        do jt=1,Ntgf
           iti = it + Nt0 - 1
           jtj = iti + 1 - jt
           do is=1,Ns
              dumpGloc(is)=Gloc_ret_tt(it,jt,is,is)
           end do
           write(unit,'(20F18.10)') t_grid(iti),t_grid(jtj),Gloc_ret_tt(it,jt,1,1)!dumpGloc(1:Ns)
        end do
        write(unit,'(6F18.10)')
        write(unit,'(6F18.10)')
     end do
     close(unit)
  end if


  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  allocate(wre(Nw))
  wre=linspace(-wrange,wrange,Nw,mesh=dw)
  allocate(Gloc_ret_tw(Ntgf,Nw,Ns,Ns)); Gloc_ret_tw=zero
  allocate(Gloc_ret_tw_(Ntgf,Nw,Ns,Ns)); Gloc_ret_tw_=zero
  do is=1,Ns
     call get_relative_time_FT(Gloc_ret_tt(:,:,is,is),Gloc_ret_tw(:,:,is,is),wre,deps)
  end do
  !  if(mpiID==0) write(*,*) Gloc_ret_tw(:,:,1,1)
  !  stop
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

  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  if(mpiID==0) then
     unit=free_unit()
     open(unit,file='Gloc_ret_tw.data')
     do it=1,Ntgf
        do iw=1,Nw
           do is=1,Ns
              dumpGloc(is)=Gloc_ret_tw(it,iw,is,is)
              dumpGloc_(is)=Gloc_ret_tw_(it,iw,is,is)
           end do
           write(unit,'(30F18.10)') t_grid(it+Nt0-1),wre(iw),dumpGloc(1:Ns),dumpGloc_(1:Ns)  
        end do
        write(*,*) it,iw,Ntgf!,Nw
        write(unit,'(20F18.10)')
     end do

     ! do it=1,Ntgf
     !    do iw=1,Nw
     !       do is=1,Ns
     !          dumpGloc(is)=zero!Gloc_ret_tw(it,iw,is,is)
     !       end do
     !       write(unit,'(20F18.10)') t_grid(it+Nt0-1),wre(iw),Gloc_ret_tw(it,iw,1,1)
     !    end do
     !    write(unit,'(6F18.10)')
     ! end do
     close(unit)
     !
     unit=free_unit()
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

     ! do it=1,Ntgf
     !    do jt=1,Ntgf
     !       iti = it + Nt0 - 1
     !       jtj = iti + 1 - jt
     !       do is=1,Ns
     !          dumpGloc(is)=Gloc_ret_tt(it,jt,is,is)
     !       end do
     ! ! !       write(unit,'(30F18.10)') t_grid(it+Nt0-1),wre(iw),dumpGloc(1:Ns),dumpGloc_(1:Ns)  
     !       write(unit,'(20F18.10)') t_grid(iti),t_grid(jtj),Gloc_ret_tt(it,jt,1,1)!dumpGloc(1:Ns)
     !    end do
     !    write(unit,'(6F18.10)')
     !    write(unit,'(6F18.10)')
     ! end do
     ! close(unit)

     ! unit=free_unit()
     ! write(*,*) 'qui',mpiID,unit
     ! open(unit,file='Gloc_ret_tw.data')
     ! ! do it=1,Ntgf
     ! !    do iw=1,Nw
     ! !       do is=1,Ns
     ! !          dumpGloc(is)=Gloc_ret_tw(it,iw,is,is)
     ! !          dumpGloc_(is)=Gloc_ret_tw_(it,iw,is,is)
     ! !       end do
     ! !       write(unit,'(30F18.10)') t_grid(it+Nt0-1),wre(iw),dumpGloc(1:Ns),dumpGloc_(1:Ns)  
     ! !    end do
     ! !    write(*,*) it,iw,Ntgf!,Nw
     ! !    write(unit,'(20F18.10)')
     ! ! end do
     ! close(unit)
     ! write(*,*) 'qui',mpiID,unit
     ! open(unit,file='Gloc_ret_tw_.data')
     ! do it=1,Ntgf
     !    do iw=1,Nw
     !       do is=1,Ns
     !          dumpGloc(is)=Gloc_ret_tw(it,iw,is,is)
     !          dumpGloc_(is)=Gloc_ret_tw_(it,iw,is,is)
     !       end do
     !       write(unit,'(30F18.10)') t_grid(it+Nt0-1),wre(iw),dumpGloc(1:Ns),dumpGloc_(1:Ns)
     !    end do
     !    write(unit,'(6F18.10)')
     !    write(unit,'(6F18.10)')
     ! end do
     ! close(unit)
  end if
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  call MPI_FINALIZE(mpiERR)
  !
contains


  ! subroutine get_neq_lgrAC(read_dir,neq_Rhop,neq_Qhop,slater_lgrA)
  !   implicit none
  !   character(len=200)  :: read_dir
  !   complex(8),dimension(2,Ns,Ns,Lk) :: slater_init
  !   complex(8),dimension(Nttgf,Ns,Ns) :: neq_Rhop,neq_Qhop
  !   complex(8),dimension(:,:,:),allocatable :: neq_Rhop_,neq_Qhop_
  !   complex(8),dimension(:),allocatable :: psi_t,psi_tmp
  !   complex(8),dimension(:,:,:),allocatable :: slater_lgrA
  !   complex(8),dimension(:),allocatable :: lgr_cmplx
  !   real(8) :: t,ti,tf
  !   integer :: it,is,js,unit
  !   !
  !   Nt_aux=2*Nttgf-1
  !   allocate(neq_Rhop_(Nt_aux,Ns,Ns),neq_Qhop_(Nt_aux,Ns,Ns))
  !   allocate(t_grid_aux(Nt_aux));
  !   ti=t_grid(1); tf=t_grid(Nttgf)
  !   t_grid_aux=linspace(ti,tf,Nt_aux,istart=.true.,iend=.true.)
  !   !t_grid_aux = linspace(tstart,0.5d0*tstep*real(Nt_aux-1,8),Nt_aux)

  !   allocate(slater_lgrA(Nttgf,Ns,Ns))

  !   do it=0,Nttgf-1
  !      !
  !      neq_Rhop_(2*it+1,:,:) = neq_Rhop(it+1,:,:)
  !      if(it.lt.Nttgf-1) then
  !         neq_Rhop_(2*it+2,:,:) = neq_Rhop(it+1,:,:)*0.5d0+neq_Rhop(it+2,:,:)*0.5d0
  !      end if
  !      !
  !      neq_Qhop_(2*it+1,:,:) = neq_Qhop(it+1,:,:)
  !      if(it.lt.Nttgf-1) then
  !         neq_Qhop_(2*it+2,:,:) = neq_Qhop(it+1,:,:)*0.5d0+neq_Qhop(it+2,:,:)*0.5d0
  !      end if
  !      !
  !   end do
  !   !
  !   nDynamics = 2*Ns*Ns*Lk 
  !   allocate(psi_t(nDynamics),psi_tmp(nDynamics))
  !   !
  !   call read_optimized_variational_wf_slater_superc(read_dir,slater_init)
  !   call wfMatrix_superc_2_dynamicalSlater(slater_init,psi_t)
  !   call setup_neq_slater_dynamics_superc(neq_Rhop_,neq_Qhop_)
  !   !
  !   Nvdm_AC_opt=2; vdm_AC_stride_v2m => vdm_AC_vec2mat ; vdm_AC_stride_m2v => vdm_AC_mat2vec
  !   allocate(lgr_cmplx(Nvdm_AC_opt))
  !   !
  !   unit=free_unit()
  !   open(unit,file='neq_lgrA.data')
  !   do it=1,Nttgf
  !      t=t_grid(it)
  !      write(*,*) it,t
  !      psi_tmp = gz_eom_slater_superc_lgr(t,psi_t,nDynamics)
  !      do is=1,Ns
  !         do js=1,Ns
  !            call get_neq_lgrA_slater(is,js,slater_lgrA(it,is,js))
  !         end do
  !      end do
  !      !
  !      call vdm_AC_stride_m2v(slater_lgrA(it,:,:),lgr_cmplx)
  !      write(unit,'(10F18.10)')  t_grid(it),lgr_cmplx
  !      !
  !   end do
  !   close(unit)
  ! end subroutine get_neq_lgrAC



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











  subroutine read_neq_lgr(read_neq_dir,neq_lgrA)
    implicit none
    character(len=200) :: read_neq_dir,file_name
    integer,dimension(:),allocatable :: ISr,JSr,ISq,JSq
    integer,dimension(:,:),allocatable :: grid_Rhop,grid_Qhop
    integer :: ir,i,ifile,it,iq,is,js
    integer :: ios,unit
    logical :: read_check
    integer :: tmp_ir,tmp_iq
    real(8) :: tmpt
    real(8),dimension(:),allocatable :: dump_vect
    complex(8),dimension(Nttgf,Ns,Ns) :: neq_lgrA
    !
    read_neq_dir=trim(read_neq_dir)    

    file_name=reg(read_neq_dir)//"neq_lgrA_constrSL.data"
    inquire(file=file_name,exist=read_check)
    if(read_check) then     
       unit=free_unit()
       open(unit,file=file_name,status='old')
       !
       ios=0
       ir=Ns*Ns
       iq=Ns*(Ns-1)/2
       allocate(dump_vect(ir)) !+- lenght of the array used to dump file
       neq_lgrA = zero
       do it=1,Ntgf
          read(unit,*,iostat=ios) tmpt,dump_vect(1:ir) !HERE
          i=0
          do is=1,Ns
             i=i+1
             !+-diagonal [only real]
             neq_lgrA(it+Nt0-1,is,is) = dump_vect(i) 
          end do
          do is=1,Ns
             do js=is+1,Ns 
                i=i+1
                !+- upper diagonal 
                neq_lgrA(it+Nt0-1,is,js) = dump_vect(i) + xi*dump_vect(i+iq)
             end do
          end do
          do is=1,Ns
             do js=1,is-1
                !+- hermitean part
                neq_lgrA(it+Nt0-1,is,js) = conjg(neq_lgrA(it+Nt0-1,js,is))
             end do
          end do
          if(abs(tmpt-t_grid(it+Nt0-1)).gt.1.d-8) stop "file to read and given tstep not compatibles"
       end do
       !+- continuation to the trivial dynamics case
       do it=1,Nt0-1
          neq_lgrA(it,:,:) = neq_lgrA(Nt0,:,:)
       end do
       close(unit)
       deallocate(dump_vect)
    else
       write(*,*) "file ",file_name,"  not present!"
       stop
    end if
    !
    open(unit,file='neq_lgrA_constrSL_read.data')
    do it=1,Nttgf
       !
       call write_hermitean_matrix(neq_lgrA(it,:,:),unit,t_grid(it))
       !
    end do
    close(unit)

  end subroutine read_neq_lgr










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




  subroutine vdm_AC_vec2mat(vdm_AC_indep,vdm_AC_mat)
    implicit none
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
    implicit none
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




end program GZ_GF
