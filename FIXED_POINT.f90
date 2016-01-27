MODULE FIXED_POINT
  implicit none

CONTAINS

  subroutine fixed_point_sub(xin,func,itmax,xtol)
    implicit none
    real(8),dimension(:),intent(inout) :: xin
    integer,optional :: itmax
    real(8),optional :: xtol
    integer :: itmax_=200
    real(8) :: xtol_=1.d-8
    real(8) :: tmp_test
    integer :: dimX,i,iter
    real(8),dimension(size(xin)) :: P0,P1,P2,D,P,relerr
    real(8),dimension(:),allocatable :: Xarray
    real(8) :: Xscalar
    logical :: check

    interface 
       function func(x_)
         implicit none
         real(8),dimension(:) :: x_
         real(8),dimension(size(x_))     :: func
       end function func
    end interface

    
    if(present(xtol))  xtol_  = xtol
    if(present(itmax)) itmax_ = itmax
    
    dimX = size(xin)
    allocate(Xarray(dimX))
    Xarray = Xin
    P0 = Xarray
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
       P0=P
    end do

    

  end subroutine fixed_point_sub
END MODULE FIXED_POINT
