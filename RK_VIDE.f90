MODULE RK_IDE
  implicit none
  private

  !RK coefficients
  real(8),dimension(:,:),allocatable :: A   
  real(8),dimension(:),allocatable   :: B,C

  complex(8),parameter :: Zi = (0.d0,1.d0)

  interface RK_step
     module procedure RK_IDE_step_d,RK_IDE_step_c,RK_step_z,RK_step_d
  end interface

  public :: RK_step
  
CONTAINS

  
  !+-------------------------+!
  !   REAL VALUED FUNCTIONS   !
  !+-------------------------+!
  FUNCTION RK_step_d(Nsys,mRK,h,t,y_old,funct) result(y_new)
    implicit none
    !INPUT
    integer                               ::  mRK      !RK order 
    integer                               ::  Nsys     !Number of equations
    integer                               ::  Nker     !Number of kernels
    real(8)                               ::  h        !time step
    real(8)                               ::  t        !time
    real(8),dimension(Nsys),intent(in)    ::  y_old    !y(t_n)

    interface
       ! function f(t,y(t))
       function funct(t,y,Nsys)
         implicit none
         integer  :: Nsys
         real(8),dimension(Nsys) :: y
         real(8)                 :: t
         real(8),dimension(Nsys) :: funct
       end function funct

    end interface

    !OUTPUT
    real(8),dimension(Nsys)               ::  y_new    !(t_{n+1})

    !internal variables
    integer                               :: irk
    integer                               :: iker
    integer                               :: j
    real(8),dimension(:,:),allocatable    :: Fn,Yn,Zn
    real(8),dimension(Nsys)               :: f


    if(mRK.ne.2.and.mRK.ne.4) then
       write(*,*) mRK,'RK order set to default value 2'
       mRK = 2
    end if

    !+- allocate and initialize stuff -+!
    allocate(A(mRK,mRK),B(mRK),C(mRK))
    allocate(Yn(Nsys,mRK))
    call RKcoeff(A,B,C,mRK)

    do j = 1,Nsys
       y_new(j) = y_old(j)
       Yn(j,:)  = y_old(j)
    end do

    !+- RUNGE-KUTTA step -+!
    do irk=1,mRK
       
       !+ Y(:,irk) +!
       do j=1,mRK
          if(A(irk,j).gt.1.d-5) then
             Yn(:,irk) = Yn(:,irk) + h*A(irk,j)*funct(t+C(j)*h,Yn(:,j),Nsys) 
          end if
       end do

       !+- step -+! 
       y_new(:) = y_new(:) + h*B(irk)*funct(t+C(irk)*h,Yn(:,irk),Nsys) 

    end do
    
    !+- deallocate -+!
    deallocate(A,B,C)
    deallocate(Yn)
    return
    
  END FUNCTION RK_step_d



  FUNCTION RK_step_z(Nsys,mRK,h,t,y_old,funct) result(y_new)
    implicit none
    !INPUT
    integer                                  ::  mRK      !RK order 
    integer                                  ::  Nsys     !Number of equations
    real(8)                                  ::  h        !time step
    real(8)                                  ::  t        !time
    complex(8),dimension(Nsys),intent(in)    ::  y_old    !y(t_n)
    
    interface
       ! function f(t,y(t))
       function funct(t,y,Nsys)
         implicit none
         integer                    :: Nsys
         complex(8),dimension(Nsys) :: y
         real(8)                    :: t
         complex(8),dimension(Nsys) :: funct
       end function funct

    end interface

    !OUTPUT
    complex(8),dimension(Nsys)               ::  y_new    !(t_{n+1})

    !internal variables
    integer                               :: irk
    integer                               :: iker
    integer                               :: j
    complex(8),dimension(:,:),allocatable    :: Fn,Yn,Zn
    complex(8),dimension(Nsys)               :: f

    if(mRK.ne.2.and.mRK.ne.4) then
       write(*,*) mRK,'RK order set to default value 2'
       mRK = 2
    end if

    !+- allocate and initialize stuff -+!
    allocate(A(mRK,mRK),B(mRK),C(mRK))
    allocate(Yn(Nsys,mRK))
    call RKcoeff(A,B,C,mRK)

    do j = 1,Nsys
       y_new(j) = y_old(j)
       Yn(j,:)  = y_old(j)
    end do

    !+- RUNGE-KUTTA step -+!
    do irk=1,mRK
       
       !+ Y(:,irk) +!
       do j=1,mRK
          if(A(irk,j).gt.1.d-5) then
             Yn(:,irk) = Yn(:,irk) + h*A(irk,j)*funct(t+C(j)*h,Yn(:,j),Nsys) 
          end if
       end do

       !+- step -+! 
       y_new(:) = y_new(:) + h*B(irk)*funct(t+C(irk)*h,Yn(:,irk),Nsys) 
    end do
    
    !+- deallocate -+!
    deallocate(A,B,C)
    deallocate(Yn)
    return
    
  END FUNCTION RK_step_z



  


  
  !+------------------------------------------------------------------+!
  !                                                                    !
  !  RK_VIDE (Volterra Integro-Differential Equations) evolution step  !
  !  dy/dt = f(t,y(t)) + int ds k(t,s,y(s))                            !
  !  spcial case of separable kernels                                  !
  !  k(t,t',y(t')) = kl(t)kr(t',y(t'))                                 !
  !                                                                    !
  !+------------------------------------------------------------------+!
  
  
  !+-------------------------+!
  !   REAL VALUED FUNCTIONS   !
  !+-------------------------+!
  FUNCTION RK_IDE_step_d(Nsys,Nker,mRK,h,t,y_old,funct,Lkernel,Rkernel,lag) result(y_new)
    implicit none
    !INPUT
    integer                               ::  mRK      !RK order 
    integer                               ::  Nsys     !Number of equations
    integer                               ::  Nker     !Number of kernels
    real(8)                               ::  h        !time step
    real(8)                               ::  t        !time
    real(8),dimension(Nsys),intent(in)    ::  y_old    !y(t_n)

    interface
       ! function f(t,y(t))
       function funct(t,y,Nsys)
         implicit none
         integer  :: Nsys
         real(8),dimension(Nsys) :: y
         real(8)                 :: t
         real(8),dimension(Nsys) :: funct
       end function funct

       ! Left kernel kl(t)
       function Lkernel(t,Nsys,Nker)
         implicit none
         integer                      :: Nsys
         integer                      :: Nker
         real(8)                      :: t
         real(8),dimension(Nsys,Nker) :: Lkernel
       end function Lkernel

       ! right kernel  kr(t,y(t))
       function Rkernel(t,y,Nsys,Nker)
         implicit none
         integer                      :: Nsys
         integer                      :: Nker
         real(8),dimension(Nsys)      :: y
         real(8)                      :: t
         real(8),dimension(Nsys,Nker) :: Rkernel
       end function Rkernel

    end interface
    real(8),dimension(Nsys,Nker),intent(inout) ::  lag      !"lag" term (Brunner)     

    !OUTPUT
    real(8),dimension(Nsys)               ::  y_new    !(t_{n+1})

    !internal variables
    integer                               :: irk
    integer                               :: iker
    integer                               :: j
    real(8),dimension(:,:),allocatable    :: Fn,Yn,Zn
    real(8),dimension(Nsys)               :: f
    real(8),dimension(Nsys,Nker)          :: Lk,Rk

    if(mRK.ne.2.and.mRK.ne.4) then
       write(*,*) mRK,'RK order set to default value 2'
       mRK = 2
    end if

    !+- allocate and initialize stuff -+!
    allocate(A(mRK,mRK),B(mRK),C(mRK))
    allocate(Yn(Nsys,mRK),Fn(Nsys,mRK),Zn(Nsys,mRK))
    call RKcoeff(A,B,C,mRK)

    do j = 1,Nsys
       y_new(j) = y_old(j)
       Yn(j,:)  = y_old(j)
       Zn(j,:)  = 0.d0
       Fn(j,:)  = 0.d0
    end do

    !+- RUNGE-KUTTA step -+!
    do irk=1,mRK
       
       !+ left kernels +!
       Lk = Lkernel(t+C(irk)*h,Nsys,Nker)
       
       !+ "lag" term +!
       do iker = 1,Nker
          Fn(:,irk) = Fn(:,irk) + Lk(:,iker)*lag(:,iker)
       end do
       
       !+ Y(:,irk) +!
       do j=1,mRK
          Yn(:,irk) = Yn(:,irk) + h*A(irk,j)*funct(t+C(j)*h,Yn(:,j),Nsys) + &
               h*A(irk,j)*Fn(:,j) + h*A(irk,j)*h*Zn(:,j)
       end do

       !+ increment function +!
       do j=1,mRK
          Lk = Lkernel(t+C(irk)*h,Nsys,Nker)
          Rk = Rkernel(t+C(j)*h,Yn(:,j),Nsys,Nker)
          do iker = 1,Nker
             Zn(:,irk) = Zn(:,irk) + A(irk,j)*Rk(:,iker)*Lk(:,iker)
          end do
       end do
       
       !+- step -+! 
       y_new(:) = y_new(:) + h*B(irk)*funct(t+C(irk)*h,Yn(:,irk),Nsys) + &
            h*B(irk)*Fn(:,irk) + h*B(irk)*h*Zn(:,irk)
    end do
    
    !+- update lag term -+!
    do irk =1,mRK
       Rk = Rkernel(t+C(irk)*h,Yn(:,irk),Nsys,Nker)
       do iker = 1,Nker
          lag(:,iker) = lag(:,iker) + B(irk)*Rk(:,iker)*h
       end do
    end do
    
    !+- deallocate -+!
    deallocate(A,B,C)
    deallocate(Yn,Fn,Zn)
    return
    
  END FUNCTION RK_IDE_step_d

  !+----------------------------+!
  !   COMPLEX VALUED FUNCTIONS   !
  !+----------------------------+!
  FUNCTION RK_IDE_step_c(Nsys,Nker,mRK,h,t,y_old,funct,Lkernel,Rkernel,lag) result(y_new)
    implicit none
    !INPUT
    integer                               ::  mRK      !RK order 
    integer                               ::  Nker     !Number of kernels
    integer                               ::  Nsys     !Number of eqs. in the sys
    real(8)                               ::  h        !time step
    real(8)                               ::  t        !time
    complex(8),dimension(Nsys),intent(in) ::  y_old    !y(t_n)
    
    interface
       ! function f(t,y(t))
       function funct(t,y,Nsys)
         implicit none
         integer  :: Nsys
         complex(8),dimension(Nsys) :: y
         real(8)                    :: t
         complex(8),dimension(Nsys) :: funct
       end function funct

       ! Left kernel kl(t)
       function Lkernel(t,Nsys,Nker)
         implicit none
         integer                         :: Nker
         integer                         :: Nsys
         real(8)                         :: t
         complex(8),dimension(Nsys,Nker) :: Lkernel
       end function Lkernel

       ! right kernel  kr(t,y(t))
       function Rkernel(t,y,Nsys,Nker)
         implicit none
         integer                         :: Nker
         integer                         :: Nsys
         complex(8),dimension(Nsys)      :: y
         real(8)                         :: t
         complex(8),dimension(Nsys,Nker) :: Rkernel
       end function Rkernel

    end interface
    complex(8),dimension(Nsys,Nker),intent(inout) ::  lag      !"lag" term (Brunner)     

    !OUTPUT
    complex(8),dimension(Nsys)               ::  y_new    !(t_{n+1})

    !internal variables
    integer                                  :: irk,iker
    integer                                  :: j
    complex(8),dimension(:,:),allocatable    :: Fn,Yn,Zn
    complex(8),dimension(Nsys,Nker)          :: Lk,Rk
    complex(8),dimension(Nsys)               :: f
    complex(8),dimension(Nsys)               :: test
    
    if(mRK.ne.2.and.mRK.ne.4) then
       write(*,*) mRK,'RK order set to default value 2'
       mRK = 2
    end if
    
    !+- allocate and initialize stuff -+!
    allocate(A(mRK,mRK),B(mRK),C(mRK))
    allocate(Yn(Nsys,mRK),Fn(Nsys,mRK),Zn(Nsys,mRK))
    call RKcoeff(A,B,C,mRK)

    do j = 1,Nsys
       y_new(j) = y_old(j)
       Yn(j,:)  = y_old(j)
       Zn(j,:)  = 0.d0
       Fn(j,:)  = 0.d0
    end do

    !+- RUNGE-KUTTA step -+!
    do irk=1,mRK
       ! left kernels
       Lk = Lkernel(t+C(irk)*h,Nsys,Nker)
       ! "lag" term
       do iker = 1,Nker
          Fn(:,irk) = Fn(:,irk) + Lk(:,iker)*lag(:,iker)
       end do
              
       ! Y(:,irk)
       do j=1,mRK
          Yn(:,irk) = Yn(:,irk) + h*A(irk,j)*funct(t+C(j)*h,Yn(:,j),Nsys) + &
               h*A(irk,j)*Fn(:,j) + h*A(irk,j)*h*Zn(:,j)
          if(dabs(A(irk,j)).gt.1.d-5) then
             test = funct(t+C(j)*h,Yn(1,j),Nsys) + Fn(:,j) + h*Zn(:,j)
             !             write(40,*) AIMAG(test(Nsys)),AIMAG(Yn(Nsys,j))
          end if
       end do


       
       !increment function
       do j=1,mRK
          Lk = Lkernel(t+C(irk)*h,Nsys,Nker)
          Rk = Rkernel(t+C(j)*h,Yn(:,j),Nsys,Nker)
          do iker = 1,Nker
             Zn(:,irk) = Zn(:,irk) + A(irk,j)*Rk(:,iker)*Lk(:,iker)
          end do
       end do
       !+- step -+! 
       y_new(:) = y_new(:) + h*B(irk)*funct(t+C(irk)*h,Yn(:,irk),Nsys) + &
            h*B(irk)*Fn(:,irk) + h*B(irk)*h*Zn(:,irk)
       
       test = funct(t+C(irk)*h,Yn(:,irk),Nsys) + Fn(:,irk) + h*Zn(:,irk)
       write(40,*) irk,AIMAG(test(Nsys)),t!,AIMAG(Yn(Nsys,irk)),AIMAG(Yn(Nsys,irk)+Zi)
    end do

    write(40,*)
    
    !+- update lag term -+!
    do irk =1,mRK
       Rk = Rkernel(t+C(irk)*h,Yn(:,irk),Nsys,Nker)
       do iker = 1,Nker
          lag(:,iker) = lag(:,iker) + B(irk)*Rk(:,iker)*h
       end do
    end do
    
    
    !+- deallocate -+!
    deallocate(A,B,C)
    deallocate(Yn,Fn,Zn)
    return
    
  END FUNCTION RK_IDE_step_c


  SUBROUTINE RKcoeff(A,B,C,mRK)
    integer,intent(in)                       :: mRK
    real(8),dimension(mRK,mRK),intent(inout) :: A
    real(8),dimension(mRK),intent(inout)     :: B,C
    
    A=0.d0
    B=0.d0
    C=0.d0
    
    if(mRK.eq.2) then
       A(2,1) = 0.5d0
       C(2) = 0.5d0
       B(2) = 1.d0
    end if
        
    if(mRK.eq.4) then
       A(2,1) = 0.5d0
       A(3,2) = 0.5d0
       A(4,3) = 1.d0
       
       C(2) = 0.5d0
       C(3) = 0.5d0
       C(4) = 1.d0
       
       B(1) = 1.d0/6.d0
       B(2) = 1.d0/3.d0
       B(3) = 1.d0/3.d0
       B(4) = 1.d0/6.d0
    end if
    
  END SUBROUTINE RKcoeff
  
END MODULE RK_IDE






