!+-------------------------------------------------------------------+
!PURPOSE  : input state |i> of the basis and calculates |j>=Cm|i>
!the sign of j has the phase convention
!m labels the sites
!+-------------------------------------------------------------------+
subroutine c(m,i,j,sgn)
  integer :: ib(Ns)
  integer :: i,j,m,km
  integer :: isg
  real(8) :: sgn
  call bdecomp(i,ib)
  if (ib(m)==0)then
     j=0
  else
     if(m==1)then
        j=i-1
     else
        km=sum(ib(1:m-1))
        isg=(-1)**km
        j=(i-2**(m-1))*isg
     endif
  endif
  sgn=dfloat(j)/dfloat(abs(j))
  if(j==0) sgn=0
  j=abs(j)
end subroutine c


subroutine cdg(m,i,j,sgn)
  integer :: ib(Ns)
  integer :: i,j,m,km
  integer :: isg
  real(8) :: sgn
  call bdecomp(i,ib)
  if (ib(m)==1)then
     j=0
  else
     if(m==1)then
        j=i+1
     else
        km=sum(ib(1:m-1))
        isg=(-1)**km
        j=(i+2**(m-1))*isg
     endif
  endif
  sgn=dfloat(j)/dfloat(abs(j))
  if(j==0) sgn=0
  j=abs(j)
end subroutine cdg


