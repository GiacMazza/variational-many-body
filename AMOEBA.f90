MODULE MIN_AMOEBA
  USE SF_IOTOOLS

CONTAINS

  SUBROUTINE amoeba(p,y,ftol,func,iter,verbose)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: iter
    REAL(8), INTENT(IN) :: ftol
    REAL(8), DIMENSION(:), INTENT(INOUT) :: y
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: p
    LOGICAL, OPTIONAL :: verbose
    integer           :: unit
    INTERFACE
       FUNCTION func(x)
         !USE nrtype
         IMPLICIT NONE
         REAL(8), DIMENSION(:), INTENT(IN) :: x
         REAL(8) :: func
       END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: ITMAX=200
    REAL(8), PARAMETER :: TINY=1.0d-10
    INTEGER :: ihi,ndim
    REAL(8), DIMENSION(size(p,2)) :: psum
    if(verbose) then
       unit=free_unit()
       open(unit,file='amoeba_verbose.out')
    end if
    call amoeba_private
    if(verbose) close(unit)
  CONTAINS
    !BL
    SUBROUTINE amoeba_private
      IMPLICIT NONE
      INTEGER :: i,ilo,inhi
      REAL(8) :: rtol,ysave,ytry,ytmp
      ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')

      iter=0
      psum(:)=sum(p(:,:),dim=1)
      do                  
         ilo=iminloc(y(:))
         ihi=imaxloc(y(:))
         ytmp=y(ihi)
         y(ihi)=y(ilo)
         inhi=imaxloc(y(:))
         y(ihi)=ytmp
         rtol=2.0d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
         if(verbose) write(unit,'(5F18.10)') dble(iter),rtol,y(ihi),y(ilo)
         if (rtol < ftol) then
            call swap_r(y(1),y(ilo))
            call swap_rv(p(1,:),p(ilo,:))                        
            RETURN
         end if
         if (iter >= ITMAX) then 
            call nrerror('ITMAX exceeded in amoeba')
            call swap_r(y(1),y(ilo))
            call swap_rv(p(1,:),p(ilo,:))
            RETURN
         end if
         ytry=amotry(-1.0d0)
         iter=iter+1

         if (ytry <= y(ilo)) then
            ytry=amotry(2.0d0)
            iter=iter+1
         else if (ytry >= y(inhi)) then
            ysave=y(ihi)
            ytry=amotry(0.5d0)
            iter=iter+1
            if (ytry >= ysave) then
               p(:,:)=0.5d0*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
               do i=1,ndim+1
                  if (i /= ilo) y(i)=func(p(i,:))
               end do
               iter=iter+ndim
               psum(:)=sum(p(:,:),dim=1)
            end if
         end if
      end do
    END SUBROUTINE amoeba_private
    !BL
    FUNCTION amotry(fac)
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: fac
      REAL(8) :: amotry
      REAL(8) :: fac1,fac2,ytry
      REAL(8), DIMENSION(size(p,2)) :: ptry
      fac1=(1.0d0-fac)/ndim
      fac2=fac1-fac
      ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
      ytry=func(ptry)
      if (ytry < y(ihi)) then
         y(ihi)=ytry
         psum(:)=psum(:)-p(ihi,:)+ptry(:)
         p(ihi,:)=ptry(:)
      end if
      amotry=ytry
    END FUNCTION amotry

  END SUBROUTINE amoeba



  SUBROUTINE amoeba_er(p,y,ftol,func,iter,verbose)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: iter
    REAL(8), INTENT(IN) :: ftol
    REAL(8), DIMENSION(:), INTENT(INOUT) :: y
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: p
    LOGICAL, OPTIONAL :: verbose
    integer           :: unit
    INTERFACE
       FUNCTION func(x)
         !USE nrtype
         IMPLICIT NONE
         REAL(8), DIMENSION(:), INTENT(IN) :: x
         REAL(8) :: func
       END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: ITMAX=400
    REAL(8), PARAMETER :: TINY=1.0d-10
    INTEGER :: ihi,ndim
    REAL(8), DIMENSION(size(p,2)) :: psum
    if(verbose) then
       unit=free_unit()
       open(unit,file='amoeba_ER_verbose.out')
    end if
    call amoeba_private
    if(verbose) close(unit)
  CONTAINS
    !BL
    SUBROUTINE amoeba_private
      IMPLICIT NONE
      INTEGER :: i,ilo,inhi
      REAL(8) :: rtol,ysave,ytry,ytmp
      ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')

      iter=0
      psum(:)=sum(p(:,:),dim=1)
      do                  
         ilo=iminloc(y(:))
         ihi=imaxloc(y(:))
         ytmp=y(ihi)
         y(ihi)=y(ilo)
         inhi=imaxloc(y(:))
         y(ihi)=ytmp
         rtol=2.0d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
         if(verbose) write(unit,'(5F18.10)') dble(iter),rtol,y(ihi),y(ilo)
         if (rtol < ftol) then
            call swap_r(y(1),y(ilo))
            call swap_rv(p(1,:),p(ilo,:))                        
            RETURN
         end if
         if (iter >= ITMAX) then 
            call nrerror('ITMAX exceeded in amoeba')
            call swap_r(y(1),y(ilo))
            call swap_rv(p(1,:),p(ilo,:))
            RETURN
         end if
         ytry=amotry(-1.0d0)
         iter=iter+1

         if (ytry <= y(ilo)) then
            ytry=amotry(2.0d0)
            iter=iter+1
         else if (ytry >= y(inhi)) then
            ysave=y(ihi)
            ytry=amotry(0.5d0)
            iter=iter+1
            if (ytry >= ysave) then
               p(:,:)=0.5d0*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
               do i=1,ndim+1
                  if (i /= ilo) y(i)=func(p(i,:))
               end do
               iter=iter+ndim
               psum(:)=sum(p(:,:),dim=1)
            end if
         end if
      end do
    END SUBROUTINE amoeba_private
    !BL
    FUNCTION amotry(fac)
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: fac
      REAL(8) :: amotry
      REAL(8) :: fac1,fac2,ytry
      REAL(8), DIMENSION(size(p,2)) :: ptry
      fac1=(1.0d0-fac)/ndim
      fac2=fac1-fac
      ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
      ytry=func(ptry)
      if (ytry < y(ihi)) then
         y(ihi)=ytry
         psum(:)=psum(:)-p(ihi,:)+ptry(:)
         p(ihi,:)=ptry(:)
      end if
      amotry=ytry
    END FUNCTION amotry

  END SUBROUTINE amoeba_er

  

  
  FUNCTION assert_eq(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq
    if (n1 == n2 .and. n2 == n3) then
       assert_eq=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq


  FUNCTION imaxloc(arr)
    REAL(8), DIMENSION(:), INTENT(IN) :: arr
    INTEGER :: imaxloc
    INTEGER, DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc=imax(1)
  END FUNCTION imaxloc


  FUNCTION iminloc(arr)
    REAL(8), DIMENSION(:), INTENT(IN) :: arr
    INTEGER, DIMENSION(1) :: imin
    INTEGER :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  END FUNCTION iminloc


  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    !STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror



  SUBROUTINE swap_r(a,b)
    REAL(8), INTENT(INOUT) :: a,b
    REAL(8) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r
  !BL
  SUBROUTINE swap_rv(a,b)
    REAL(8), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(8), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv
  

end MODULE MIN_AMOEBA
