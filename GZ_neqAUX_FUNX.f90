MODULE GZ_neqAUX_FUNX
  USE GZ_VARS_GLOBAL
  USE SF_IOTOOLS
  USE SF_LINALG
  implicit none
  private
  !
  public :: read_optimized_variational_wf
  public :: wfMatrix_2_dynamicalVector,dynamicalVector_2_wfMatrix
  public :: t2it
  !
CONTAINS
  !
  subroutine read_optimized_variational_wf(read_dir,slater_determinant,phi_vector)
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
  end subroutine read_optimized_variational_wf



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


  function t2it(time,delta_t) result(it)
    real(8) :: time,delta_t
    integer :: it
    real(8) :: eps=1.d-8    
    it = (time+eps)/delta_t + 1
  end function t2it


END MODULE GZ_NEQAUX_FUNX
