MODULE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE SF_IOTOOLS
  USE SF_LINALG
  implicit none
  private
  !
  public :: bdecomp
  public :: vec2mat_stride,mat2vec_stride
  public :: initialize_variational_density_simplex
  public :: initialize_variational_density  
  public :: init_Rhop_seed,init_Qhop_seed
  public :: get_free_single_paritcle_gf

  !+- at some point these two subroutines should be merged in SCIFOR -+!
  public :: simultaneous_diag
  public :: fixed_point_sub
  !
  public :: fermi_zero
  !
  public :: search_chemical_potential
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
       close(unit)
    else
       write(*,*) 'vdm_simplex_seed.conf does not exist'
       write(*,*) 'please provide an initialization file for the simplex'
       stop
    end if
  end subroutine initialize_variational_density_simplex






  subroutine initialize_variational_density(variational_density)
    real(8),dimension(:),intent(inout) :: variational_density
    logical                 :: IOfile
    integer                 :: unit,flen,i,j,expected_flen
    if(size(variational_density).ne.Nvdm_NC_opt-Nvdm_NCoff_opt) stop "initialize variational density /= Nvdm_NC_opt-Nvdm_NCoff_opt"
    expected_flen=size(variational_density)
    inquire(file="vdm_seed.conf",exist=IOfile)
    if(IOfile) then
       flen=file_length("vdm_seed.conf")
       unit=free_unit()
       open(unit,file="vdm_seed.conf")
       write(*,*) 'reading denisty seed from file vdm_seed.conf'
       if(flen.eq.expected_flen) then
          !+- read from file -+!
          do i=1,flen
             read(unit,*) variational_density(i)
          end do
       else
          write(*,*) 'vdm_seed.conf in the wrong form',flen,expected_flen
          write(*,*) 'variational density matrix will be initialized to the homogeneous case'
          variational_density = 0.5d0
       end if
       close(unit)
    else
       write(*,*) 'vdm_simplex_seed.conf does not exist'
       write(*,*) 'variational density matrix will be initialized to the homogeneous case'
       variational_density = 0.5d0
    end if
  end subroutine initialize_variational_density




  subroutine init_Rhop_seed(Rhop_init)
    complex(8),dimension(Ns,Ns) :: Rhop_init
    logical :: IOfile
    real(8) :: buffer_real,buffer_imag  
    integer :: expected_flen,flen,unit,is,js
    complex(8),dimension(NRhop_opt) :: tmp_init
    !
    inquire(file="Rhop_seed.conf",exist=IOfile)
    if(IOfile) then
       expected_flen=Ns*Ns
       flen=file_length("Rhop_seed.conf")
       unit=free_unit()
       open(unit,file="Rhop_seed.conf")
       write(*,*) 'reading Rhop_seed from file Rhop_seed.conf'
       if(flen.eq.expected_flen) then
          !+- read from file -+!
          do is=1,Ns
             do js=1,Ns
                read(unit,*) buffer_real,buffer_imag
                Rhop_init(is,js)= buffer_real+xi*buffer_imag
             end do
          end do
       else
          write(*,*) 'Rhop_seed.conf in the wrong form',flen,expected_flen
          write(*,*) 'Rhop will be initialized to the non-interacting case'
          !
          if(.not.gz_superc) then
             tmp_init = 1.d0
          else
             tmp_init = 1.d0/sqrt(2.d0)
          end if
          call Rhop_stride_v2m(tmp_init,Rhop_init)
          !
       end if
       close(unit)
    else
       !
       if(.not.gz_superc) then
          tmp_init = 1.d0
       else
          tmp_init = 1.d0/sqrt(2.d0)
       end if
       call Rhop_stride_v2m(tmp_init,Rhop_init)
       !
    end if
    !
    do is=1,Ns
       write(*,'(20F8.4)') dreal(Rhop_init(is,1:Ns)),dimag(Rhop_init(is,1:Ns))
    end do
    !
  end subroutine init_Rhop_seed




  subroutine init_Qhop_seed(Qhop_init)
    complex(8),dimension(Ns,Ns) :: Qhop_init
    logical :: IOfile
    real(8) :: buffer_real,buffer_imag  
    integer :: expected_flen,flen,unit,is,js
    complex(8),dimension(NRhop_opt) :: tmp_init
    !
    inquire(file="Qhop_seed.conf",exist=IOfile)
    if(IOfile) then
       expected_flen=Ns*Ns
       flen=file_length("Qhop_seed.conf")
       unit=free_unit()
       open(unit,file="Qhop_seed.conf")
       write(*,*) 'reading Qhop_seed from file Qhop_seed.conf'
       if(flen.eq.expected_flen) then
          !+- read from file -+!
          do is=1,Ns
             do js=1,Ns
                read(unit,*) buffer_real,buffer_imag
                Qhop_init(is,js)= buffer_real+xi*buffer_imag
             end do
          end do
       else
          write(*,*) 'Qhop_seed.conf in the wrong form',flen,expected_flen
          write(*,*) 'Qhop will be initialized to the non-interacting case'
          !
          tmp_init = 1.d0/sqrt(2.d0)
          call Qhop_stride_v2m(tmp_init,Qhop_init)
          !
       end if
       close(unit)
    else
       !
       tmp_init = 1.d0/sqrt(2.d0)
       call Qhop_stride_v2m(tmp_init,Qhop_init)
       !
    end if
    !
    do is=1,Ns
       write(*,'(20F8.4)') dreal(Qhop_init(is,1:Ns)),dimag(Qhop_init(is,1:Ns))
    end do
    !
  end subroutine init_Qhop_seed





















  subroutine search_chemical_potential(var,ntmp,converged,nread,ndelta,nerr,niter_)
    real(8),intent(inout) :: var
    real(8),intent(in)    :: ntmp
    logical,intent(inout) :: converged
    real(8),intent(in)    :: nread,nerr
    real(8),intent(inout) :: ndelta
    integer,optional      :: niter_
    integer               :: niter
    logical               :: bool
    real(8)               :: ndiff
    integer,save          :: count=0,totcount=0,i
    integer,save          :: nindex=0
    integer               :: nindex_old(3)
    real(8)               :: ndelta_old,nratio
    integer,save          :: nth_magnitude=-2,nth_magnitude_old=-2
    real(8),save          :: nth=1.d-2
    logical,save          :: ireduce=.true.
    integer               :: unit
    !
    niter=20;if(present(niter_)) niter=niter_
    !
    ndiff=ntmp-nread
    nratio = 0.5d0;!nratio = 1.d0/(6.d0/11.d0*pi)
    !
    !check actual value of the density *ntmp* with respect to goal value *nread*
    count=count+1
    totcount=totcount+1
    if(count>2)then
       do i=1,2
          nindex_old(i+1)=nindex_old(i)
       enddo
    endif
    nindex_old(1)=nindex
    !
    if(ndiff >= nth)then
       nindex=-1
    elseif(ndiff <= -nth)then
       nindex=1
    else
       nindex=0
    endif
    
    !
    write(*,*) 'cristo santo sto xmu searching'
    write(*,*) ndiff,nth,nindex
    !
    
    !
    ndelta_old=ndelta
    bool=nindex/=0.AND.( (nindex+nindex_old(1)==0).OR.(nindex+sum(nindex_old(:))==0) )
    !if(nindex_old(1)+nindex==0.AND.nindex/=0)then !avoid loop forth and back
    if(bool)then
       ndelta=ndelta_old*nratio !decreasing the step
    else
       ndelta=ndelta_old
    endif
    !
    if(ndelta_old<1.d-9)then
       ndelta_old=0.d0
       nindex=0
    endif
    !update chemical potential
    var=var+dble(nindex)*ndelta
    !
    write(*,*) 'var',var,nindex,ndelta
    !
    !xmu=xmu+dble(nindex)*ndelta
    !
    !Print information
    write(*,"(A,f16.9,A,f15.9)")"n    = ",ntmp," /",nread
    if(nindex>0)then
       write(*,"(A,es16.9,A)")"shift= ",nindex*ndelta," ==>"
    elseif(nindex<0)then
       write(*,"(A,es16.9,A)")"shift= ",nindex*ndelta," <=="
    else
       write(*,"(A,es16.9,A)")"shift= ",nindex*ndelta," == "
    endif
    write(*,"(A,f15.9)")"var  = ",var
    write(*,"(A,ES16.9,A,ES16.9)")"dn   = ",ndiff,"/",nth
    unit=free_unit()
    open(unit,file="GZsearch_mu_iteration.out",position="append")
    write(unit,*)var,ntmp,ndiff,nth
    close(unit)
    !
    !check convergence within actual threshold
    !if reduce is activetd
    !if density is in the actual threshold
    !if DMFT is converged
    !if threshold is larger than nerror (i.e. this is not last loop)
    bool=ireduce.AND.(abs(ndiff)<nth).AND.converged.AND.(nth>nerr)
    if(bool)then
       nth_magnitude_old=nth_magnitude        !save old threshold magnitude
       nth_magnitude=nth_magnitude_old-1      !decrease threshold magnitude || floor(log10(abs(ntmp-nread)))
       nth=max(nerr,10.d0**(nth_magnitude))   !set the new threshold 
       count=0                                !reset the counter
       converged=.false.                      !reset convergence
       ndelta=ndelta_old*nratio                  !reduce the delta step
       !
    endif
    !
    !if density is not converged set convergence to .false.
    if(abs(ntmp-nread)>nth)converged=.false.
    !
    !check convergence for this threshold
    !!---if smallest threshold-- NO MORE
    !if reduce is active (you reduced the treshold at least once)
    !if # iterations > max number
    !if not yet converged
    !set threshold back to the previous larger one.
    !bool=(nth==nerr).AND.ireduce.AND.(count>niter).AND.(.not.converged)
    bool=ireduce.AND.(count>niter).AND.(.not.converged)
    if(bool)then
       ireduce=.false.
       nth=10.d0**(nth_magnitude_old)
    endif
    !
    write(*,"(A,I5)")"count= ",count
    write(*,"(A,L2)"),"Converged=",converged
    print*,""
    !
  end subroutine search_chemical_potential






  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine get_free_single_paritcle_gf(ek,gk,wr,eps)
    real(8)             :: ek
    complex(8),dimension(:)          :: gk
    real(8),dimension(:)         :: wr
    real(8),optional :: eps
    real(8) :: eta
    integer                       :: i,ik,Lk,L
    real(8)                       :: w,dew
    complex(8)                    :: gf,iw
    !
    if(size(gk).ne.size(wr)) stop 'size error in get free single particle gf'
    L=size(gk)
    eta  = 0.01d0 ; if(present(eps))eta=eps
    dew=abs(wfin-wini)/real(L,8)
    do i=1,L
       !w  = wini + dble(i-1)*dew;
       w=wr(i)
       iw=cmplx(w,eta)
       gk(i) = 1.d0/(iw-ek)
    enddo
  end subroutine get_free_single_paritcle_gf







  subroutine simultaneous_diag(A,V,diagA,eps,iinit)
    implicit none
    complex(8),dimension(:,:,:) :: A
    complex(8),dimension(:,:)   :: V
    integer,optional :: iinit
    complex(8),dimension(:,:),allocatable   :: Rgivens
    complex(8),dimension(:,:),allocatable   :: Vinit
    real(8),dimension(:),allocatable   :: tmp_eigen
    complex(8),dimension(:,:,:),allocatable :: A_input
    complex(8),dimension(:,:,:),allocatable   :: tmp
    complex(8),dimension(:,:),allocatable   :: tmpV
    real(8),dimension(:,:)      :: diagA
    real(8),dimension(:,:),allocatable      :: tmp_diag
    real(8),dimension(3,3)      :: G
    complex(8),dimension(3)      :: Gvec
    complex(8),dimension(2,2)      :: tmp_rotate,tmpA,tmp_rotate_dag
    real(8),dimension(3) :: Geigen,Jacobi_vec
    real(8)                     :: eps,C,off
    complex(8)                  :: S
    integer                     :: N,M,NM,iinit_
    integer                     :: iloop,imax,i,j,k
    !
    real(8) :: min,max
    integer,dimension(:),allocatable :: sort,tmp_sort,ivec
    integer :: imin
    !
    integer                     :: ii,jj,kk,kkk
    logical :: flag_init
    !
    imax=200
    !
    if(size(A,1).ne.size(A,2)) stop "simultaneous diag: wrong dimensions" 
    N=size(A,1);M=size(A,3);NM=N*M
    if(size(V,1).ne.size(V,2)) stop "simultaneous diag: wrong dimensions" 
    if(size(V,1).ne.N) stop "simultaneous diag: wrong dimensions" 

    if(size(diagA,1).ne.N) stop "simultaneous diag: wrong dimensions"
    if(size(diagA,2).ne.M) stop "simultaneous diag: wrong dimensions"
    !

    allocate(A_input(N,N,M))
    A_input = A

    V=0.d0
    do i=1,N
       V(i,i) = one
    end do
    !
    iinit_=1
    if(present(iinit).and.iinit.le.M) then
       !
       iinit_=iinit
       allocate(Vinit(N,N)); Vinit=A(:,:,iinit_)
       allocate(tmp_eigen(N))
       call matrix_diagonalize(Vinit,tmp_eigen)       
       A=zero
       do k=1,M       
          do i=1,N
             do j=1,N
                A(i,j,k) = zero
                do ii=1,N
                   do jj=1,N
                      A(i,j,k) = A(i,j,k) + &
                           conjg(V(ii,i))*V(jj,j)*A_input(ii,jj,k)
                   end do
                end do
             end do
          end do
       end do
       !
    end if
    !
    allocate(tmp(N,N,M),Rgivens(N,N),tmpV(N,N));
    !
    imax=50
    do iloop=1,imax

       off = 0.d0
       do k=1,M
          off = off + off_diag(A(:,:,k))
       end do
       if(off.lt.eps) exit
       !
       do i=1,N
          do j=1,N           
             if(i.ne.j) then
                G=0.d0     
                do k=1,M
                   Gvec(1) = A(i,i,k) - A(j,j,k)
                   Gvec(2) = A(i,j,k) + A(j,i,k)
                   Gvec(3) = xi*(A(j,i,k) - A(i,j,k))
                   G = G + get_gmatrix(Gvec)
                end do
                !              
                call matrix_diagonalize(G,Geigen,'V','U')
                Jacobi_vec = G(:,3)
                !
                !Givens rotations
                if(Jacobi_vec(1).lt.0.d0) Jacobi_vec=-Jacobi_vec
                C = sqrt(Jacobi_vec(1)*0.5d0+0.5d0)
                S = Jacobi_vec(2)-xi*Jacobi_vec(3)
                S = 0.5d0*S/C
                !
                Rgivens=0.d0
                do k=1,N
                   Rgivens(k,k) = 1.d0
                end do
                !+- R A R^+ rotation -+!
                Rgivens(i,i) = C
                Rgivens(i,j) = conjg(S)
                Rgivens(j,i) = -S
                Rgivens(j,j) = C
                !
                do k=1,M              
                   !+-  compute rotations R A R^+  -+!
                   tmp(:,:,k)=A(:,:,k)              
                   do kk=1,N
                      tmp(kk,i,k) = A(kk,i,k)*conjg(Rgivens(i,i)) + A(kk,j,k)*conjg(Rgivens(i,j))
                      tmp(KK,j,k) = A(kk,i,k)*conjg(Rgivens(j,i)) + A(kk,j,k)*conjg(Rgivens(j,j))
                   end do
                   A(:,:,k)=tmp(:,:,k)
                   do kk=1,N
                      A(i,kk,k) =  Rgivens(i,i)*tmp(i,kk,k) + Rgivens(i,j)*tmp(j,kk,k)
                      A(j,kk,k) =  Rgivens(j,i)*tmp(i,kk,k) + Rgivens(j,j)*tmp(j,kk,k)
                   end do
                   !
                end do
                !
                !+- PAY ATTENTION: here I want to put the rotation in the standard form V^+ A V -+!
                tmpV = V
                do kk=1,N
                   V(kk,i) = tmpV(kk,i)*conjg(Rgivens(i,i)) + tmpV(kk,j)*conjg(Rgivens(i,j))
                   V(kk,j) = tmpV(kk,i)*conjg(Rgivens(j,i)) + tmpV(kk,j)*conjg(Rgivens(j,j))
                end do
             end if
          end do
       end do
       !
    end do
    !
    do k=1,M
       do i=1,N
          diagA(i,k)=A(i,i,k)
       end do
    end do
    !
    do i=1,N
       write(*,'(20F7.3)') diagA(i,:)
    end do
    !+- sort eigenvalues for k=1 -+!
    allocate(sort(N),tmp_sort(N)); forall(i=1:N) sort(i)=0; tmp_sort=sort  
    do i=1,N
       max = -10     
       do j=1,N
          if(tmp_sort(j).eq.0.and.diagA(j,1).gt.max) then
             max = diagA(j,1)
             imax = j
          end if
       end do
       tmp_sort(imax) = 1
       sort(i) = imax
    end do
    !

    !+- sort the transformation matrix according to the previous sorting -+!
    tmpV = V
    do i=1,N
       do j=1,N
          V(i,j) = tmpV(i,sort(j))
       end do
    end do
    !< SAFETY-TEST
    tmp=0.d0
    do k=1,M
       do i=1,N
          do j=1,N
             tmp(i,j,k) = 0.d0
             do ii=1,N
                do jj=1,N
                   tmp(i,j,k) = tmp(i,j,k) + conjg(V(ii,i))*V(jj,j)*A_input(ii,jj,k)
                end do
             end do
          end do
       end do
    end do
    !
    off = 0.d0
    do k=1,M
       off = off + off_diag(tmp(:,:,k))
    end do
    if(off.lt.eps) then
       do i=1,N
          diagA(i,:) = tmp(i,i,:)
       end do
    else
       stop "failures in the diagonalization"
    end if

    ! tmp_diag=diagA
    ! do i=1,N
    !    !diagA(i,:) = tmp_diag(sort(i),:)
    !    !write(*,'(I4,20F8.3)') sort(i),diagA(i,1),diagA(i,2),diagA(i,5)
    !    write(*,'(I4,20F8.3)') i,diagA(i,1),diagA(i,2),diagA(i,5)
    ! end do

    ! write(*,*)
    !
    !+- TO DO LIST -+!
    !+- match and group the equal entries of diagA
    !+- obtain map
    !+- sort the corresponding diagA according to one of the column values
    !  !
  contains
    !
    function get_Gmatrix(h) result(G)
      implicit none
      complex(8),dimension(:) :: h
      real(8),dimension(size(h),size(h)) :: G
      integer :: N,i,j
      N=size(h)
      ! if(size(G,1).ne.size(G,2)) stop 'wrong dimensions'
      ! if(size(G,1).ne.N) stop 'wrong dimensions'
      !
      do i=1,N
         do j=1,N
            G(i,j) = dreal(conjg(h(i))*h(j))
         end do
      end do
      !

    end function get_Gmatrix


    !
    function off_diag(A) result(off)
      complex(8),dimension(:,:) :: A
      real(8) :: off
      integer :: i,j,N
      N=size(A,1)
      !
      off = 0.d0
      do i=1,N
         do j=1,N
            if(i.ne.j) off = off + abs(A(i,j))**2.d0
         end do
      end do
      !
    end function off_diag
    !
  end subroutine simultaneous_diag



  function fermi_zero(x,b) result(f)
    real(8) :: x,f,b
    if(x.lt.0.d0) then
       f=1.d0
    else       
       f=0.d0
    end if
    if(x==0.d0) f=0.5d0
  end function fermi_zero


  subroutine fixed_point_sub(xin,func,itmax,xtol)
    implicit none
    real(8),dimension(:),intent(inout) :: xin
    integer,optional :: itmax
    real(8),optional :: xtol
    integer :: itmax_=200
    real(8) :: xtol_=1.d-8
    real(8) :: tmp_test
    integer :: dimX,i,iter
    real(8),dimension(size(xin)) :: P0,P1,P2,D,P,Pold,relerr
    real(8),dimension(:),allocatable :: Xarray
    real(8) :: Xscalar
    real(8) :: mix=0.1
    logical :: check
    interface 
       function func(x_)
         implicit none
         real(8),dimension(:),intent(in) :: x_
         real(8),dimension(size(x_))     :: func
       end function func
    end interface
    !
    if(present(xtol))  xtol_  = xtol
    if(present(itmax)) itmax_ = itmax
    dimX = size(xin)
    allocate(Xarray(dimX))
    Xarray = Xin
    P0 = Xarray
    Pold=P0
    do iter=1,itmax_
       !
       P1 = func(P0)
       P2 = func(P1)
       D = P2 - 2.d0*P1 + P0
       !
       check=.true.
       do i=1,dimX
          if(D(i) == 0.d0) then
             P(i) = P2(i)
          else
             P(i) = P0(i) - (P1(i)-P0(i))*(P1(i)-P0(i))/D(i)
          end if
          if(P0(i) == 0.d0) then
             relerr(i) = P(i) 
          else
             relerr(i) = (P(i)-P0(i))/P0(i)
          end if
          if(abs(relerr(i)).gt.xtol_) check=.false.
       end do
       if(check) then
          Xin=P
          return
       end if
       !
       P0=P
       Pold=P
       !
    end do
    write(*,*) "FIXED POINT:exceede number of iterations"
  end subroutine fixed_point_sub



END MODULE GZ_AUX_FUNX
