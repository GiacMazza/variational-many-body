MODULE GZ_neqAUX_FUNX
  USE GZ_VARS_GLOBAL
  USE SF_IOTOOLS
  USE SF_LINALG
  implicit none
  private
  !

  interface read_optimized_variational_wf
     module procedure read_optimized_variational_wf_normal,read_optimized_variational_wf_normal_,read_optimized_variational_wf_superc,read_optimized_variational_wf_superc_
  end interface read_optimized_variational_wf

  public :: read_optimized_variational_wf!,read_optimized_variational_wf_superc
  public :: wfMatrix_2_dynamicalVector,dynamicalVector_2_wfMatrix
  public :: wfMatrix_superc_2_dynamicalVector,dynamicalVector_2_wfMatrix_superc
  public :: wfMatrix_superc_2_dynamicalVector_,dynamicalVector_2_wfMatrix_superc_
  public :: slater_full2reduced
  public :: BCSwf_2_dynamicalVector,dynamicalVector_2_BCSwf
  public :: t2it
  !
CONTAINS
  !
  subroutine read_optimized_variational_wf_normal(read_dir,slater_determinant,phi_vector)
    character(len=200)             :: read_dir
    complex(8),dimension(Nphi)     :: phi_vector
    complex(8),dimension(Ns,Ns,Lk) :: slater_determinant    
    real(8)                        :: tmp_re,tmp_im,x
    integer :: is,js,read_unit,flen,ios,i
    character(len=200) :: file_name
    logical :: check_file
    !
    read_dir=trim(read_dir)
    do is=1,Ns
       do js=1,Ns

          file_name=reg(read_dir)//'optimized_slater_normal_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data'
          inquire(file=file_name,exist=check_file)
          if(check_file) then
             read_unit = free_unit()
             open(read_unit,file=file_name,status='old')
             flen=0
             do
                read (read_unit,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(read_unit)                    
             if(flen.ne.Lk) stop "READING OPTIMIZED SLATER DETRMINANT: number of lines /= Lk"
             open(read_unit,file=file_name,status='old')
             write(*,*) "Reading from file",file_name
             do i=1,flen
                read(read_unit,'(2F18.10)') tmp_re,tmp_im
                slater_determinant(is,js,i) = tmp_re + xi*tmp_im
             end do
             close(read_unit)
          else
             write(*,*) "FILE",file_name,"does not exist!!!"
             stop 
          end if
       end do
    end do
    !
    file_name=reg(read_dir)//'optimized_projectors.data'
    inquire(file=file_name,exist=check_file)
    if(check_file) then
       read_unit=free_unit()
       open(read_unit,file=file_name,status='old')
       flen=0
       do
          read (read_unit,*,iostat=ios) tmp_re
          if (ios/=0) exit     
          flen=flen+1
       end do
       close(read_unit)                    
       if(flen.ne.Nphi) stop "READING OPTIMIZED GZ PROJECTORS: number of lines /= Nphi"
       open(read_unit,file=file_name,status='old')
       do i=1,flen
          read(read_unit,'(2F18.10)') tmp_re,tmp_im
          phi_vector(i) = tmp_re + xi*tmp_im
       end do
       close(read_unit)
    else
       write(*,*) 'FILE',file_name,'does not exist!!!'
    end if
  end subroutine read_optimized_variational_wf_normal



  subroutine read_optimized_variational_wf_normal_(read_dir,slater_determinant,phi_matrix)
    character(len=200)             :: read_dir
    complex(8),dimension(nFock,nFock)     :: phi_matrix
    complex(8),dimension(Ns,Ns,Lk) :: slater_determinant    
    real(8)                        :: tmp_re,tmp_im,x
    integer :: is,js,read_unit,flen,ios,i,tmp_if,tmp_jf
    character(len=200) :: file_name
    logical :: check_file
    !
    read_dir=trim(read_dir)
    do is=1,Ns
       do js=1,Ns

          file_name=reg(read_dir)//'optimized_slater_normal_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data'
          inquire(file=file_name,exist=check_file)
          if(check_file) then
             read_unit = free_unit()
             open(read_unit,file=file_name,status='old')
             flen=0
             do
                read (read_unit,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(read_unit)                    
             if(flen.ne.Lk) stop "READING OPTIMIZED SLATER DETRMINANT: number of lines /= Lk"
             open(read_unit,file=file_name,status='old')
             write(*,*) "Reading from file",file_name
             do i=1,flen
                read(read_unit,'(2F18.10)') tmp_re,tmp_im
                slater_determinant(is,js,i) = tmp_re + xi*tmp_im
             end do
             close(read_unit)
          else
             write(*,*) "FILE",file_name,"does not exist!!!"
             stop 
          end if
       end do
    end do
    !
    file_name=reg(read_dir)//'optimized_phi_matrix.data'
    inquire(file=file_name,exist=check_file)
    if(check_file) then
       read_unit=free_unit()
       open(read_unit,file=file_name,status='old')
       flen=0
       do
          read (read_unit,*,iostat=ios) tmp_re
          if (ios/=0) exit     
          flen=flen+1
       end do
       close(read_unit)                    
       if(flen.ne.nFock*nFock) stop "READING OPTIMIZED GZ PHI-PROJECTORS: number of lines /= nFock*nFock"
       open(read_unit,file=file_name,status='old')
       phi_matrix=zero
       do i=1,flen
          read(read_unit,'(F18.10,2I2)') tmp_re,tmp_if,tmp_jf
          tmp_im=0.d0
          phi_matrix(tmp_if,tmp_jf) = tmp_re + xi*tmp_im
       end do
       close(read_unit)
    else
       write(*,*) 'FILE',file_name,'does not exist!!!'
    end if
  end subroutine read_optimized_variational_wf_normal_




  subroutine read_optimized_variational_wf_superc(read_dir,slater_determinant,phi_vector,slater_lgr_A,proj_lgr_A)
    character(len=200)             :: read_dir
    complex(8),dimension(Nphi)     :: phi_vector
    complex(8),dimension(2,Ns,Ns,Lk) :: slater_determinant    
    complex(8),dimension(Ns,Ns),optional :: slater_lgr_A,proj_lgr_A
    real(8)                        :: tmp_re,tmp_im,x
    integer :: is,js,read_unit,flen,ios,i
    character(len=200) :: file_name
    logical :: check_file
    !
    read_dir=trim(read_dir)
    do is=1,Ns
       do js=1,Ns
          file_name=reg(read_dir)//'optimized_slater_normal_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data'
          inquire(file=file_name,exist=check_file)
          if(check_file) then
             read_unit = free_unit()
             open(read_unit,file=file_name,status='old')
             flen=0
             do
                read (read_unit,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(read_unit)                    
             if(flen.ne.Lk) stop "READING OPTIMIZED SLATER DETRMINANT NORMAL: number of lines /= Lk"
             open(read_unit,file=file_name,status='old')
             write(*,*) "Reading from file",file_name
             do i=1,flen
                read(read_unit,'(2F18.10)') tmp_re,tmp_im
                slater_determinant(1,is,js,i) = tmp_re + xi*tmp_im
             end do
             close(read_unit)
          else
             write(*,*) "FILE",file_name,"does not exist!!!"
             stop 
          end if
       end do
    end do
    !
    do is=1,Ns
       do js=1,Ns
          file_name=reg(read_dir)//'optimized_slater_anomalous_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data'
          inquire(file=file_name,exist=check_file)
          if(check_file) then
             read_unit = free_unit()
             open(read_unit,file=file_name,status='old')
             flen=0
             do
                read (read_unit,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(read_unit)                    
             if(flen.ne.Lk) stop "READING OPTIMIZED SLATER DETRMINANT ANOMALOUS: number of lines /= Lk"
             open(read_unit,file=file_name,status='old')
             write(*,*) "Reading from file",file_name
             do i=1,flen
                read(read_unit,'(2F18.10)') tmp_re,tmp_im
                slater_determinant(2,is,js,i) = tmp_re + xi*tmp_im
             end do
             close(read_unit)
          else
             write(*,*) "FILE",file_name,"does not exist!!!"
             stop 
          end if
       end do
    end do
    !
    file_name=reg(read_dir)//'optimized_projectors.data'
    inquire(file=file_name,exist=check_file)
    if(check_file) then
       read_unit=free_unit()
       open(read_unit,file=file_name,status='old')
       write(*,*) "Reading from file",file_name
       flen=0
       do
          read (read_unit,*,iostat=ios) tmp_re
          if (ios/=0) exit     
          flen=flen+1
       end do
       close(read_unit)                    
       if(flen.ne.Nphi) stop "READING OPTIMIZED GZ PROJECTORS: number of lines /= Nphi"
       open(read_unit,file=file_name,status='old')
       do i=1,flen
          read(read_unit,'(2F18.10)') tmp_re,tmp_im
          phi_vector(i) = tmp_re + xi*tmp_im
       end do
       close(read_unit)
    else
       write(*,*) 'FILE',file_name,'does not exist!!!'
    end if


    if(present(slater_lgr_A)) then
       !
       do is=1,Ns
          do js=1,Ns
             file_name=reg(read_dir)//'optimized_slaterLGR_anomalous_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data'
             inquire(file=file_name,exist=check_file)
             if(check_file) then
                read_unit = free_unit()
                open(read_unit,file=file_name,status='old')
                flen=0
                do
                   read (read_unit,*,iostat=ios) x
                   if (ios/=0) exit     
                   flen=flen+1
                end do
                close(read_unit)                    
                if(flen.ne.1) stop "READING OPTIMIZED SLATER DETRMINANT LGR ANOMALOUS: number of lines /= 1"
                open(read_unit,file=file_name,status='old')
                write(*,*) "Reading from file",file_name
                do i=1,flen
                   read(read_unit,'(2F18.10)') tmp_re,tmp_im
                   slater_lgr_A(is,js) = tmp_re + xi*tmp_im
                end do
                close(read_unit)
             else
                write(*,*) "FILE",file_name,"does not exist!!!"
                stop 
             end if
          end do
       end do
    end if


    if(present(proj_lgr_A)) then
       !
       do is=1,Ns
          do js=1,Ns
             file_name=reg(read_dir)//'optimized_gzprojLGR_anomalous_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data'
             inquire(file=file_name,exist=check_file)
             if(check_file) then
                read_unit = free_unit()
                open(read_unit,file=file_name,status='old')
                flen=0
                do
                   read (read_unit,*,iostat=ios) x
                   if (ios/=0) exit     
                   flen=flen+1
                end do
                close(read_unit)                    
                if(flen.ne.1) stop "READING OPTIMIZED GZ PROJ DETRMINANT LGR ANOMALOUS: number of lines /= 1"
                open(read_unit,file=file_name,status='old')
                write(*,*) "Reading from file",file_name
                do i=1,flen
                   read(read_unit,'(2F18.10)') tmp_re,tmp_im
                   proj_lgr_A(is,js) = tmp_re + xi*tmp_im
                end do
                close(read_unit)
             else
                write(*,*) "FILE",file_name,"does not exist!!!"
                stop 
             end if
          end do
       end do
    end if

  end subroutine read_optimized_variational_wf_superc






  subroutine read_optimized_variational_wf_superc_(read_dir,slater_determinant,phi_matrix,slater_lgr_A,proj_lgr_A)
    character(len=200)             :: read_dir
    complex(8),dimension(nFock,nFock)     :: phi_matrix
    complex(8),dimension(2,Ns,Ns,Lk) :: slater_determinant    
    complex(8),dimension(Ns,Ns),optional :: slater_lgr_A,proj_lgr_A
    real(8)                        :: tmp_re,tmp_im,x
    integer :: is,js,read_unit,flen,ios,i,tmp_if,tmp_jf
    character(len=200) :: file_name
    logical :: check_file
    !
    read_dir=trim(read_dir)
    do is=1,Ns
       do js=1,Ns
          file_name=reg(read_dir)//'optimized_slater_normal_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data'
          inquire(file=file_name,exist=check_file)
          if(check_file) then
             read_unit = free_unit()
             open(read_unit,file=file_name,status='old')
             flen=0
             do
                read (read_unit,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(read_unit)                    
             if(flen.ne.Lk) stop "READING OPTIMIZED SLATER DETRMINANT NORMAL: number of lines /= Lk"
             open(read_unit,file=file_name,status='old')
             write(*,*) "Reading from file",file_name
             do i=1,flen
                read(read_unit,'(2F18.10)') tmp_re,tmp_im
                slater_determinant(1,is,js,i) = tmp_re + xi*tmp_im
             end do
             close(read_unit)
          else
             write(*,*) "FILE",file_name,"does not exist!!!"
             stop 
          end if
       end do
    end do
    !
    do is=1,Ns
       do js=1,Ns
          file_name=reg(read_dir)//'optimized_slater_anomalous_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data'
          inquire(file=file_name,exist=check_file)
          if(check_file) then
             read_unit = free_unit()
             open(read_unit,file=file_name,status='old')
             flen=0
             do
                read (read_unit,*,iostat=ios) x
                if (ios/=0) exit     
                flen=flen+1
             end do
             close(read_unit)                    
             if(flen.ne.Lk) stop "READING OPTIMIZED SLATER DETRMINANT ANOMALOUS: number of lines /= Lk"
             open(read_unit,file=file_name,status='old')
             write(*,*) "Reading from file",file_name
             do i=1,flen
                read(read_unit,'(2F18.10)') tmp_re,tmp_im
                slater_determinant(2,is,js,i) = tmp_re + xi*tmp_im
             end do
             close(read_unit)
          else
             write(*,*) "FILE",file_name,"does not exist!!!"
             stop 
          end if
       end do
    end do
    !
    file_name=reg(read_dir)//'optimized_phi_matrix.data'
    inquire(file=file_name,exist=check_file)
    if(check_file) then
       read_unit=free_unit()
       open(read_unit,file=file_name,status='old')
       flen=0
       do
          read (read_unit,*,iostat=ios) tmp_re
          if (ios/=0) exit     
          flen=flen+1
       end do
       close(read_unit)                    
       if(flen.ne.nFock*nFock) stop "READING OPTIMIZED GZ PHI MATRIX: number of lines /= nFock*nFock"
       open(read_unit,file=file_name,status='old')
       write(*,*) "Reading from file",file_name
       phi_matrix=zero
       do i=1,flen
          read(read_unit,*) tmp_re,tmp_if,tmp_jf
          tmp_im=zero
          phi_matrix(tmp_if,tmp_jf) = tmp_re + xi*tmp_im
       end do
       close(read_unit)
    else
       write(*,*) 'FILE',file_name,'does not exist!!!'
    end if


    if(present(slater_lgr_A)) then
       !
       do is=1,Ns
          do js=1,Ns
             file_name=reg(read_dir)//'optimized_slaterLGR_anomalous_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data'
             inquire(file=file_name,exist=check_file)
             if(check_file) then
                read_unit = free_unit()
                open(read_unit,file=file_name,status='old')
                flen=0
                do
                   read (read_unit,*,iostat=ios) x
                   if (ios/=0) exit     
                   flen=flen+1
                end do
                close(read_unit)                    
                if(flen.ne.1) stop "READING OPTIMIZED SLATER DETRMINANT LGR ANOMALOUS: number of lines /= 1"
                open(read_unit,file=file_name,status='old')
                write(*,*) "Reading from file",file_name
                do i=1,flen
                   read(read_unit,'(2F18.10)') tmp_re,tmp_im
                   slater_lgr_A(is,js) = tmp_re + xi*tmp_im
                end do
                close(read_unit)
             else
                write(*,*) "FILE",file_name,"does not exist!!!"
                stop 
             end if
          end do
       end do
    end if


    if(present(proj_lgr_A)) then
       !
       do is=1,Ns
          do js=1,Ns
             file_name=reg(read_dir)//'optimized_gzprojLGR_anomalous_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data'
             inquire(file=file_name,exist=check_file)
             if(check_file) then
                read_unit = free_unit()
                open(read_unit,file=file_name,status='old')
                flen=0
                do
                   read (read_unit,*,iostat=ios) x
                   if (ios/=0) exit     
                   flen=flen+1
                end do
                close(read_unit)                    
                if(flen.ne.1) stop "READING OPTIMIZED GZ PROJ DETRMINANT LGR ANOMALOUS: number of lines /= 1"
                open(read_unit,file=file_name,status='old')
                write(*,*) "Reading from file",file_name
                do i=1,flen
                   read(read_unit,'(2F18.10)') tmp_re,tmp_im
                   proj_lgr_A(is,js) = tmp_re + xi*tmp_im
                end do
                close(read_unit)
             else
                write(*,*) "FILE",file_name,"does not exist!!!"
                stop 
             end if
          end do
       end do
    end if

  end subroutine read_optimized_variational_wf_superc_





  subroutine wfMatrix_2_dynamicalVector(slater,gzproj,dynVect)
    complex(8),dimension(Ns,Ns,Lk) :: slater
    complex(8),dimension(Nphi)      :: gzproj
    complex(8),dimension(nDynamics):: dynVect
    integer :: i,j,is,js,iphi,idyn,ik    
    idyn=0
    do ik=1,Lk
       do is=1,Ns
          do js=1,Ns
             idyn=idyn+1
             dynVect(idyn) = slater(is,js,ik)
          end do
       end do
    end do
    do iphi=1,Nphi
       idyn=idyn+1
       dynVect(idyn) = gzproj(iphi)
    end do
    if(idyn.ne.nDynamics) stop 'dimension problems wfMatrix_2_dynamicalVector'
  end subroutine wfMatrix_2_dynamicalVector
  !
  subroutine dynamicalVector_2_wfMatrix(dynVect,slater,gzproj)
    complex(8),dimension(Ns,Ns,Lk) :: slater
    complex(8),dimension(Nphi)      :: gzproj
    complex(8),dimension(nDynamics):: dynVect
    integer :: i,j,is,js,iphi,idyn,ik  
    idyn=0
    do ik=1,Lk
       do is=1,Ns
          do js=1,Ns
             idyn=idyn+1
             slater(is,js,ik) = dynVect(idyn)
          end do
       end do
    end do
    do iphi=1,Nphi
       idyn=idyn+1
       gzproj(iphi) = dynVect(idyn) 
    end do
    if(idyn.ne.nDynamics) stop 'dimension problems wfMatrix_2_dynamicalVector'
  end subroutine dynamicalVector_2_wfMatrix


  subroutine wfMatrix_superc_2_dynamicalVector(slater,gzproj,dynVect)
    complex(8),dimension(2,Ns,Ns,Lk) :: slater
    complex(8),dimension(Nphi)      :: gzproj
    complex(8),dimension(nDynamics):: dynVect
    integer :: i,j,is,js,iphi,idyn,ik    
    idyn=0
    do ik=1,Lk
       do is=1,Ns
          do js=1,Ns
             idyn=idyn+1
             dynVect(idyn) = slater(1,is,js,ik)
          end do
       end do
    end do
    do ik=1,Lk
       do is=1,Ns
          do js=1,Ns
             idyn=idyn+1
             dynVect(idyn) = slater(2,is,js,ik)
          end do
       end do
    end do
    do iphi=1,Nphi
       idyn=idyn+1
       dynVect(idyn) = gzproj(iphi)
    end do
    if(idyn.ne.nDynamics) stop 'dimension problems wfMatrix_superc_2_dynamicalVector'
  end subroutine wfMatrix_superc_2_dynamicalVector
  subroutine wfMatrix_superc_2_dynamicalVector_(slater_normal,slater_anomalous,gzproj,dynVect)
    complex(8),dimension(Nsl_normal_opt,Lk) :: slater_normal
    complex(8),dimension(Nsl_anomalous_opt,Lk) :: slater_anomalous
    complex(8),dimension(Nphi)      :: gzproj
    complex(8),dimension(nDynamics):: dynVect
    integer :: i,j,is,js,iphi,idyn,ik    
    idyn=0
    do ik=1,Lk
       do i=1,Nsl_normal_opt
          idyn = idyn + 1
          dynVect(idyn) = slater_normal(i,ik)
       end do
    end do
    do ik=1,Lk
       do i=1,Nsl_anomalous_opt
          idyn = idyn + 1
          dynVect(idyn) = slater_anomalous(i,ik)
       end do
    end do
    !
    do iphi=1,Nphi
       idyn=idyn+1
       dynVect(idyn) = gzproj(iphi)
    end do
    if(idyn.ne.nDynamics) stop 'dimension problems wfMatrix_superc_2_dynamicalVector_'
  end subroutine wfMatrix_superc_2_dynamicalVector_



  
  !
  subroutine dynamicalVector_2_wfMatrix_superc(dynVect,slater,gzproj)
    complex(8),dimension(2,Ns,Ns,Lk) :: slater
    complex(8),dimension(Nphi)      :: gzproj
    complex(8),dimension(nDynamics):: dynVect
    integer :: i,j,is,js,iphi,idyn,ik  
    idyn=0
    do ik=1,Lk
       do is=1,Ns
          do js=1,Ns
             idyn=idyn+1
             slater(1,is,js,ik) = dynVect(idyn)
          end do
       end do
    end do
    do ik=1,Lk
       do is=1,Ns
          do js=1,Ns
             idyn=idyn+1
             slater(2,is,js,ik) = dynVect(idyn)
          end do
       end do
    end do
    do iphi=1,Nphi
       idyn=idyn+1
       gzproj(iphi) = dynVect(idyn) 
    end do
    if(idyn.ne.nDynamics) stop 'dimension problems dynamicalVector_2_wfMatrix_superc'
  end subroutine dynamicalVector_2_wfMatrix_superc
  subroutine dynamicalVector_2_wfMatrix_superc_(dynVect,slater_normal,slater_anomalous,gzproj)
    complex(8),dimension(Nsl_normal_opt,Lk) :: slater_normal
    complex(8),dimension(Nsl_anomalous_opt,Lk) :: slater_anomalous
    complex(8),dimension(Nphi)      :: gzproj
    complex(8),dimension(nDynamics):: dynVect
    integer :: i,j,is,js,iphi,idyn,ik  
    idyn=0
    do ik=1,Lk
       do i=1,Nsl_normal_opt
          idyn=idyn+1
          slater_normal(i,ik) = dynVect(idyn)
       end do
    end do
    do ik=1,Lk
       do i=1,Nsl_anomalous_opt
          idyn=idyn+1
          slater_anomalous(i,ik) = dynVect(idyn)
       end do
    end do
    !
    do iphi=1,Nphi
       idyn=idyn+1
       gzproj(iphi) = dynVect(idyn) 
    end do
    !
    if(idyn.ne.nDynamics) stop 'dimension problems dynamicalVector_2_wfMatrix_superc_'
  end subroutine dynamicalVector_2_wfMatrix_superc_




  subroutine slater_full2reduced(slater,slater_normal,slater_anomalous)
    complex(8),dimension(2,Ns,Ns,Lk) :: slater
    complex(8),dimension(Nsl_normal_opt,Lk) :: slater_normal
    complex(8),dimension(Nsl_anomalous_opt,Lk) :: slater_anomalous
    integer :: ik
    do ik=1,Lk
       ! call sl_normal_stride_v2m(slater_normal(:,ik),slater(1,:,:,ik))
       ! call sl_anomalous_stride_v2m(slater_anomalous(:,ik),slater(2,:,:,ik))
       !
       call sl_normal_stride_m2v(slater(1,:,:,ik),slater_normal(:,ik))
       call sl_anomalous_stride_m2v(slater(2,:,:,ik),slater_anomalous(:,ik))
       !
    end do
  end subroutine slater_full2reduced

  




  


  subroutine BCSwf_2_dynamicalVector(bcs_wf,dynVect)
    complex(8),dimension(3,Lk) :: bcs_wf
    complex(8),dimension(3*LK):: dynVect
    integer :: i,j,is,js,iphi,idyn,ik    
    idyn=0
    do ik=1,Lk
       idyn=idyn+1
       dynVect(idyn) = bcs_wf(1,ik)
       dynVect(idyn+Lk) = bcs_wf(2,ik)
       dynVect(idyn+2*Lk) = bcs_wf(3,ik)
    end do
    if(idyn.ne.Lk) stop 'dimension problems BCSwf_2_dynamicalVector'
  end subroutine BCSwf_2_dynamicalVector
  !
  subroutine dynamicalVector_2_BCSwf(dynVect,bcs_wf)
    complex(8),dimension(3,Lk) :: bcs_wf
    complex(8),dimension(3*Lk):: dynVect
    integer :: i,j,is,js,iphi,idyn,ik  
    idyn=0
    do ik=1,Lk
       idyn=idyn+1
       bcs_wf(1,ik) = dynVect(idyn)
       bcs_wf(2,ik) = dynVect(idyn+Lk)
       bcs_wf(3,ik) = dynVect(idyn+2*Lk)
    end do
    if(idyn.ne.Lk) stop 'dimension problems dynamicalVector_2_BCSwf'
  end subroutine dynamicalVector_2_BCSwf




  function t2it(time,delta_t) result(it)
    real(8) :: time,delta_t
    integer :: it
    real(8) :: eps=1.d-8    
    it = (time+eps)/delta_t + 1
  end function t2it


END MODULE GZ_NEQAUX_FUNX
