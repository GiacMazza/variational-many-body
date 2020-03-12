function BCS_equations_of_motion(time,y,Nsys) result(f)
  implicit none
  
  integer                                     :: Nsys ! nr of equations                                                                                                                                   
  real(8)                                     :: time ! time variable                                                                                                                                      
  complex(8),dimension(Nsys)                  :: y    ! argument array                                                                                                                                     
  complex(8),dimension(Nsys)                  :: f    ! result 
  !
  complex(8),dimension(3,Lk)               :: bcsWF,bcsWF_dot
  complex(8) :: ekt
  complex(8) :: delta_t,phi_t
  complex(8),dimension(Lk) :: delta_tk
  real(8),dimension(Lk) :: n_tk
  complex(8),dimension(Ns,Ns)                 :: Hk
  complex(8) ::  nhh_dot 
  integer                                     :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  !
  real(8) :: n_t,tmp_delta
  real(8) :: lgr_t
  real(8) :: Sz_dot,nnsum
  complex(8) :: cmu

  !HERE write the GZ EQUATIONS OF MOTION
  if(Nsys.ne.3*LK) stop "wrong dimensions in the BCS equations of motion"
  !
  call dynamicalVector_2_BCSwf(y,BCSwf)
  !
  BCSwf_dot=zero
  !
  it = t2it(time,tstep*0.5d0)
  !
  delta_t = zero
  n_t=0.d0
  do ik=1,Lk
     delta_t = delta_t + 0.5d0*(bcsWF(1,ik)+xi*bcsWF(2,ik))*wtk(ik)
     delta_tk(ik) = 0.5d0*(bcsWF(1,ik)+xi*bcsWF(2,ik))
     n_t = n_t + 1.d0*(bcsWF(3,ik)+1.d0)*wtk(ik)
     n_tk(ik) = 0.5d0*(bcsWF(3,ik)+1.d0)
  end do
  phi_t = (Ubcs_t(it)+xi*kdiss_t(it))*delta_t !
  !delta_t = (Ubcs_t(it)+xi*kdiss_t(it))*delta_t !
  !

  Sz_dot=0.d0
  nnsum=0.d0
  cmu = Ubcs_t(it)*0.5d0*(1.d0-n_t)
  do ik=1,Lk
     call get_Hk_t(Hk,ik,time)
     ekt = Hk(1,1) + dreal(cmu)
     !
     bcsWF_dot(1,ik) = -2.d0*ekt*bcsWF(2,ik) + 2.d0*dimag(phi_t)*bcsWF(3,ik)
     bcsWF_dot(1,ik) = bcsWF_dot(1,ik) - kdiss_t(it)*n_t*bcsWF(1,ik) 
     bcsWF_dot(1,ik) = bcsWF_dot(1,ik) - 2.d0*kpump_t(it)*bcsWF(1,ik) 
     !
     bcsWF_dot(2,ik) =  2.d0*ekt*bcsWF(1,ik) - 2.d0*dreal(phi_t)*bcsWF(3,ik)
     bcsWF_dot(2,ik) = bcsWF_dot(2,ik) - kdiss_t(it)*n_t*bcsWF(2,ik) 
     bcsWF_dot(2,ik) = bcsWF_dot(2,ik) - 2.d0*kpump_t(it)*bcsWF(2,ik) 
     ! !
     bcsWF_dot(3,ik) =  2.d0*dreal(phi_t)*bcsWF(2,ik) - 2.d0*dimag(phi_t)*bcsWF(1,ik)
     bcsWF_dot(3,ik) = bcsWF_dot(3,ik) - kdiss_t(it)*n_t*(bcsWF(3,ik)+1.d0)
     bcsWF_dot(3,ik) = bcsWF_dot(3,ik) + 2.d0*kpump_t(it)*(1.d0-n_tk(ik))
     !

     !+ add here the non-hermitean part -+!
     !
     nhh_dot = -delta_tk(ik)*n_tk(ik) + delta_t*n_tk(ik)**2.d0-conjg(delta_t)*delta_tk(ik)
     !
     bcsWF_dot(1,ik) = bcsWF_dot(1,ik) - 2.d0*(1.d0-a_nhh)*kdiss_t(it)*2.d0*dreal(nhh_dot)
     bcsWF_dot(2,ik) = bcsWF_dot(2,ik) - 2.d0*(1.d0-a_nhh)*kdiss_t(it)*2.d0*dimag(nhh_dot)
     !
     nhh_dot = 0.5d0*n_t*(abs(delta_tk(ik))**2.d0-n_tk(ik)**2.d0)-2.d0*dreal(delta_t*conjg(delta_tk(ik))*n_tk(ik))
     bcsWF_dot(3,ik) = bcsWF_dot(3,ik) - 2.d0*(1.d0-a_nhh)*kdiss_t(it)*nhh_dot     
     !
     Sz_dot=Sz_dot+bcsWF_dot(3,ik)*wtk(ik)
     !     
     nnsum = nnsum + (1.d0-n_tk(ik))*n_tk(ik)*wtk(ik)
     !
  end do

  if(abs(nnsum).gt.1.d-12) then
     cmu = cmu -xi*Sz_dot/4.d0/nnsum
  end if
 
  !+- complex chemical potential -+!  
  if(diss_fixdens) then
     do ik=1,Lk
        !
        nhh_dot = -4.d0*dimag(cmu)*Delta_tk(ik)*n_tk(ik) + 2.d0*xi*conjg(cmu)*Delta_tk(ik)
        bcsWF_dot(1,ik) = bcsWF_dot(1,ik) + 2.d0*dreal(nhh_dot)
        bcsWF_dot(2,ik) = bcsWF_dot(2,ik) + 2.d0*dimag(nhh_dot)
        ! 
        nhh_dot = 2.d0*dimag(cmu)*n_tk(ik)*(1.d0-n_tk(ik))
        bcsWF_dot(3,ik) = bcsWF_dot(3,ik) + 2.d0*dreal(nhh_dot)
        !
        !+
     end do
     
  end if

  ! lgr_t=0.d0
  ! if(abs(dreal(delta_t)).gt.1.d-12) then
  !    lgr_t = +Sz_dot/4.d0/dreal(delta_t)
  ! end if
  !
  ! Sz_dot=0.d0
  ! tmp_delta=0.d0
  ! phi_t=0.d0
  ! phi_t = lgr_t*xi
  ! if(diss_fixdens) then
  !    do ik=1,Lk
  !       call get_Hk_t(Hk,ik,time)
  !       ekt = Hk(1,1) + Ubcs_t(it)*0.5d0*(1.d0-n_t)
  !       !
  !       bcsWF_dot(1,ik) = bcsWF_dot(1,ik) + 2.d0*dimag(phi_t)*bcsWF(3,ik)
  !       !
  !       bcsWF_dot(2,ik) = bcsWF_dot(2,ik) - 2.d0*dreal(phi_t)*bcsWF(3,ik)
  !       ! 
  !       bcsWF_dot(3,ik) = bcsWF_dot(3,ik) + 2.d0*dreal(phi_t)*bcsWF(2,ik) - 2.d0*dimag(phi_t)*bcsWF(1,ik)
  !       !
  !       !
  !       !Sz_dot=Sz_dot+bcsWF_dot(3,ik)*wtk(ik)
  !       !+
  !    end do
     
  ! end if

  ! !write(965,*) tmp_delta*0.5d0,dimag(delta_t)
  write(965,'(10F18.10)') time,cmu
  !
  !
  call BCSwf_2_dynamicalVector(bcsWF_dot,f)
  !
end function BCS_equations_of_motion




function BCS_eom(time,y,Nsys) result(f)
  implicit none               
  integer                                     :: Nsys ! nr of equations
  real(8)                                     :: time ! time variable                                                                                                                        
  complex(8),dimension(Nsys)                  :: y    ! argument array                                                                                                                                     
  complex(8),dimension(Nsys)                  :: f    ! result 
  !
  complex(8),dimension(3,Lk)               :: bcsWF,bcsWF_dot
  complex(8) :: ekt
  complex(8) :: delta_t,phi_t
  complex(8),dimension(Lk) :: delta_tk
  real(8),dimension(Lk) :: n_tk
  complex(8),dimension(Ns,Ns)                 :: Hk
  complex(8) ::  nhh_dot 
  complex(8) :: Delta_dot,nk_dot
  integer                                     :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  !
  real(8) :: n_t
  
  !HERE write the GZ EQUATIONS OF MOTION
  if(Nsys.ne.3*LK) stop "wrong dimensions in the BCS equations of motion"
  !
  call dynamicalVector_2_BCSwf(y,BCSwf)
  !
  BCSwf_dot=zero
  !
  it = t2it(time,tstep*0.5d0)
  !
  delta_t = zero
  n_t=0.d0
  do ik=1,Lk
     delta_t = delta_t + 0.5d0*(bcsWF(1,ik)+xi*bcsWF(2,ik))*wtk(ik)
     delta_tk(ik) = 0.5d0*(bcsWF(1,ik)+xi*bcsWF(2,ik))
     n_t = n_t +0.5d0*(bcsWF(3,ik)+1.d0)*wtk(ik)
     n_tk(ik) = 0.5d0*(bcsWF(3,ik)+1.d0)
  end do
  !
  phi_t = (Ubcs_t(it)+xi*kdiss_t(it))*delta_t 
  !
  do ik=1,Lk
     call get_Hk_t(Hk,ik,time)
     ekt = Hk(1,1)
     !
     !
     delta_dot=zero
     nk_dot=0.d0
     !
     !
     delta_dot = delta_dot + 2.d0*xi*ekt*delta_tk(ik) - 2.d0*kdiss_t(it)*n_t*Delta_tk(ik)   
     delta_dot = delta_dot - (xi*Ubcs_t(it)-kdiss_t(it))*delta_t*(2.d0*n_tk(ik)-1.d0)
     !
     nk_dot = nk_dot - 2.d0*kdiss_t(it)*n_t*n_tk(ik)
     nk_dot = nk_dot - (kdiss_t(it)+xi*Ubcs_t(it))*conjg(delta_t)*delta_tk(ik)
     nk_dot = nk_dot + (-kdiss_t(it)+xi*Ubcs_t(it))*delta_t*conjg(delta_tk(ik))
     !
     !
     bcsWF_dot(1,ik) = 2.d0*dreal(delta_dot)
     bcsWF_dot(2,ik) = 2.d0*dimag(delta_dot)
     bcsWF_dot(3,ik) = 2.d0*nk_dot
     !
     !
  end do
  !
  call BCSwf_2_dynamicalVector(bcsWF_dot,f)
  !
end function BCS_eom
