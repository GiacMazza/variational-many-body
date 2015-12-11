MODULE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE SF_IOTOOLS
  implicit none
  private
  !
  !  public :: initialize_local_fock_space
  !  public :: build_local_operators_fock_space
  !  public :: build_local_hamiltonian
  public :: bdecomp
  !  public :: get_spin_indep_states
  !public :: initialize_local_density
  public :: vec2mat_stride,mat2vec_stride
  public :: initialize_variational_density_simplex
  !
CONTAINS

  !+- Get binary decomposition of the state i into the configuration vector ivec -+!
  subroutine bdecomp(i,ivec)
    integer :: ivec(:)         
    integer :: l,i,Ntot
    logical :: busy
    Ntot=size(ivec)
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i-1,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end subroutine bdecomp
  !
  subroutine mat2vec_stride(mat,vec)
    real(8) :: mat(:,:)
    real(8) :: vec(:)
    integer :: n,m,vec_size,i,j,k
    n=size(mat,1)
    m=size(mat,2)
    vec_size=n*m
    if(vec_size.ne.size(vec)) stop "mat2vec vec_size.ne.size(vec)"
    do i=1,n
       do j=1,m
          k=(i-1)*m+j
          vec(k)=mat(i,j)
       end do
    end do
  end subroutine mat2vec_stride
  !
  subroutine vec2mat_stride(vec,mat)
    real(8) :: mat(:,:)
    real(8) :: vec(:)
    integer :: n,m,vec_size,i,j,k
    n=size(mat,1)
    m=size(mat,2)
    vec_size=n*m
    if(vec_size.ne.size(vec)) stop "vec2mat vec_size.ne.size(vec)"
    do k=1,vec_size
       i=(k-1)/m+1
       j=mod((k-1),m)+1
       mat(i,j)=vec(k)
    end do
  end subroutine vec2mat_stride
  
  
  !+- Get the input file initialization of the density symplex -+!
  subroutine initialize_variational_density_simplex(variational_density_simplex)
    real(8),dimension(Ns+1,Ns),intent(inout) :: variational_density_simplex
    logical                 :: IOfile
    integer                 :: unit,flen,i,j,expected_flen
    expected_flen=(Ns+1)*(Ns+1)-1
    inquire(file="vdm_simplex_seed.conf",exist=IOfile)
    if(IOfile) then
       flen=file_length("vdm_simplex_seed.conf")
       unit=free_unit()
       open(unit,file="vdm_simplex_seed.conf")
       write(*,*) 'reading denisty seed from file vdm_simplex_seed.conf'
       if(flen.eq.expected_flen) then
          !+- read from file -+!
          do i=1,Ns+1
             do j=1,Ns
                read(unit,*) variational_density_simplex(i,j)
             end do
             write(*,*) variational_density_simplex(i,:)
             if(i.le.Ns) read(unit,*)
          end do
       else
          write(*,*) 'vdm_simplex_seed.conf in the wrong form',flen,expected_flen
          write(*,*) 'please check your initialization file for the simplex'
          stop
       end if
    else
       write(*,*) 'vdm_simplex_seed.conf does not exist'
       write(*,*) 'please provide an initialization file for the simplex'
       stop
    end if
  end subroutine initialize_variational_density_simplex







END MODULE GZ_AUX_FUNX
