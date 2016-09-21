program GUTZ_mb
  USE SCIFOR
  !
  USE DMFT_MISC
  USE SF_PARSE_INPUT
  USE RK_IDE
  !
  USE GZ_AUX_FUNX
  USE GZ_neqAUX_FUNX
  USE GZ_DYNAMICS
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_LOCAL_HAMILTONIAN
  USE GZ_ENERGY_MINIMIZATION
  USE GZ_OPTIMIZED_ENERGY
  USE GZ_MATRIX_BASIS
  !
  implicit none
  real(8),dimension(2) :: psi_t
  real(8) :: Docc,phase,Zhop,xdop,ndens,r,tRamp,ene,Ut,Uf,t
  integer :: itt,it
  !
  call parse_input_variable(Uf,"Uf","inputGZ.conf",default=0.5d0)
  call parse_input_variable(xdop,"xdop","inputGZ.conf",default=0.0d0)
  call parse_input_variable(tramp,"TRAMP","inputGZ.conf",default=0.0d0)
  call parse_input_variable(tstep,"TSTEP","inputGZ.conf",default=0.005d0)
  call parse_input_variable(Nt,"Nt","inputGZ.conf",default=10000)
  
  
  call save_input_file("inputGZ.conf")
  

  

  !
  Docc=0.25d0*(1.d0-xdop)**2.d0
  phase= 0.d0
  
  ndens=1.d0-xdop
  !

  Nt_aux=2*Nt+1
  allocate(t_grid(Nt),t_grid_aux(Nt_aux))
  t_grid = linspace(0.d0,tstep*real(Nt-1,8),Nt)
  t_grid_aux = linspace(0.d0,0.5d0*tstep*real(Nt_aux-1,8),Nt_aux)

  do itt=1,Nt_aux
     !
     t = t_grid_aux(itt)
     !
     if(t.lt.tRamp) then
        r = (1.d0 - 1.5d0*cos(pi*t/tRamp) + 0.5d0*(cos(pi*t/tRamp))**3)*0.5d0
     else
        r = 1.d0
     end if
  end do
  !
  psi_t(1) = Docc
  psi_t(2) = 0.d0
  !

  
  open(unit=100,file='docc_phase_neq.data')
  do it=1,Nt
     !
     t=t_grid(it)
     !
     Docc=psi_t(1)
     phase=psi_t(2)
     !
     if(t.lt.tRamp) then
        r = (1.d0 - 1.5d0*cos(pi*t/tRamp) + 0.5d0*(cos(pi*t/tRamp))**3)*0.5d0
     else
        r = 1.d0
     end if
     !
     Ut = Uf*r

     Zhop = 2.d0*(ndens-2.d0*Docc)/(1.d0-xdop*xdop)*((sqrt(Docc+xdop)-sqrt(Docc))**2+4.d0*dcos(phase)**2.d0*sqrt(Docc**2+Docc*xdop))
     ene = -1.d0/8.d0*Zhop*(1.d0-xdop*xdop) + Ut*Docc+Ut*xdop*0.5        
     write(100,'(10F18.10)') t,Docc,Zhop,phase,ene
     !
     psi_t = RK_step(2,4,tstep,t,psi_t,docc_phase_eom)
     !
  end do
  !
CONTAINS

  function docc_phase_eom(time,y,Nsys) result(f)
    implicit none
    !inputs
    integer                                     :: Nsys ! nr of equations
    real(8)                                     :: time ! time variable
    real(8),dimension(Nsys)                  :: y    ! argument array
    real(8),dimension(Nsys)                  :: f    ! result
    !
    real(8) :: d,ph,Ut
    integer :: it

    !it = t2it(time,tstep*0.5d0)

    !
    if(time.lt.tRamp) then
       r = (1.d0 - 1.5d0*cos(pi*time/tRamp) + 0.5d0*(cos(pi*time/tRamp))**3)*0.5d0
    else
       r = 1.d0
    end if
    !
    Ut = Uf*r
    
    
    d=y(1);ph=y(2)
    
    f(1) = -2.d0*(1-xdop-2.d0*d)*sqrt(d*(d+xdop))*dcos(ph)*dsin(ph)
    !write(110,'(10F18.10)') f(1),d,ph,d*(d+xdop),xdop!,dcos(ph),dsin(ph),(1-xdop-2.d0*d)
    !
    f(2)=0.d0
    
    f(2) = 0.5d0*((sqrt(d)-sqrt(d+xdop))**2.d0+4.d0*dcos(ph)**2.d0*sqrt(d*(d+xdop)))
    f(2) = f(2) - 0.25d0*(1.d0-xdop-2.d0*d)*(2.d0*dcos(ph)**2.d0*(2.d0*d+xdop)-(sqrt(d)-sqrt(d+xdop))**2.d0)/sqrt(d*(d+xdop))
    f(2) = f(2) + Ut

    !f=f*0.5d0
    ! write(111,*) f(1),f(2),d,ph
    !
  end function docc_phase_eom

end program GUTZ_mb



!AMOEBA TEST


