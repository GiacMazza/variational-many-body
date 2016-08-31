function gz_equations_of_motion(time,y,Nsys) result(f)
  implicit none
  !inputs                                                                                                                                                        
  integer                                     :: Nsys ! nr of equations
  real(8)                                     :: time ! time variable
  complex(8),dimension(Nsys)                  :: y    ! argument array
  complex(8),dimension(Nsys)                  :: f    ! result 
  !
  complex(8),dimension(Ns,Ns,Lk)              :: slater,slater_dot
  complex(8),dimension(Ns,Ns)                 :: Hk,tRR
  complex(8),dimension(Nphi)                  :: gzproj,gzproj_dot
  complex(8),dimension(Nphi,Nphi)             :: Hproj
  integer                                     :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  complex(8),dimension(Ns,Ns)                 :: tmpHk,Rhop,slater_derivatives,Rhop_hc
  complex(8),dimension(Ns,Ns)                 :: vdm_natural
  real(8),dimension(Ns)                       :: vdm_diag
  !

  !HERE write the GZ EQUATIONS OF MOTION
  if(Nsys.ne.nDynamics) stop "wrong dimensions in the GZ_equations_of_motion"
  !
  call dynamicalVector_2_wfMatrix(y,slater,gzproj)
  do is=1,Ns
     do js=1,Ns
        vdm_natural(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
     end do
     vdm_diag(is) = dreal(vdm_natural(is,is))
  end do
  Rhop = hopping_renormalization_normal(gzproj,vdm_diag)
  do is=1,Ns
     do js=1,Ns
        Rhop_hc(is,js) = conjg(Rhop(js,is))
     end do
  end do

  !
  slater_dot=zero
  gzproj_dot=zero
  !
  it = t2it(time,tstep*0.5d0)
  !
  slater_derivatives=zero
  do ik=1,Lk
     call get_Hk_t(Hk,ik,time)
     !
     tRR = matmul(Hk,Rhop)
     tRR = matmul(Rhop_hc,tRR)
     !
     ! do is=1,Ns
     !    do js=1,Ns
     !       tmpHk(is,js)=zero           
     !       do iis=1,Ns
     !          do jjs=1,Ns                 
     !             tmpHk(is,js) = tmpHk(is,js) + conjg(Rhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,js)
     !          end do
     !       end do
     !    end do
     ! end do
     !
     do is=1,Ns
        do js=1,Ns
           slater_dot(is,js,ik) = zero
           do ks=1,Ns
              slater_dot(is,js,ik) = slater_dot(is,js,ik) + tRR(js,ks)*slater(is,ks,ik)
              slater_dot(is,js,ik) = slater_dot(is,js,ik) - tRR(ks,is)*slater(ks,js,ik)
           end do
        end do
     end do
     !     
     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns
                 slater_derivatives(is,js) = slater_derivatives(is,js) + conjg(Rhop(kks,ks))*Hk(kks,is)*slater(ks,js,ik)*wtk(ik)
              end do
           end do
        end do
     end do
     !
  end do
  slater_dot = -xi*slater_dot

  !+- create HLOC

  Uloc=Uloc_t(:,it)
  Ust =Ust_t(it)
  Jh=Jh_t(it)
  Jsf=Jsf_t(it)
  Jph=Jph_t(it)
  eLevels = eLevels_t(:,it)
  write(777,*) time,Uloc
  call get_local_hamiltonian_trace(eLevels)  
  call build_neqH_GZproj(Hproj,slater_derivatives,Rhop,vdm_diag)
  !
  do iphi=1,Nphi
     gzproj_dot(iphi) = zero
     do jphi=1,Nphi
        gzproj_dot(iphi) = gzproj_dot(iphi) + Hproj(iphi,jphi)*gzproj(jphi)
     end do
  end do
  gzproj_dot = -xi*gzproj_dot
  !
  call wfMatrix_2_dynamicalVector(slater_dot,gzproj_dot,f)
  !
end function gz_equations_of_motion




function gz_equations_of_motion_superc(time,y,Nsys) result(f)
  implicit none
  !inputs
  integer                          :: Nsys ! nr of equations
  real(8)                          :: time ! time variable                                                                                                                                      
  complex(8),dimension(Nsys)       :: y    ! argument array                                                                                                                                     
  complex(8),dimension(Nsys)       :: f    ! result 
  !
  complex(8),dimension(2,Ns,Ns,Lk) :: slater_,slater_dot
  complex(8),dimension(3,Ns,Ns,Lk) :: slater
  complex(8),dimension(Nphi)       :: gzproj,gzproj_dot
  complex(8),dimension(Nphi,Nphi)  :: Hproj
  integer                          :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  complex(8),dimension(Ns,Ns)      :: tmpHk
  complex(8),dimension(Ns,Ns)      :: Hk
  complex(8),dimension(2*Ns,2*Ns)  :: tmp
  complex(8),dimension(2,Ns,Ns)    :: slater_derivatives
  complex(8),dimension(Ns,Ns)      :: Rhop,Qhop,Rhop_hc,Qhop_hc
  complex(8),dimension(Ns,Ns)      :: tRR,tRQ,tQR,tQQ
  complex(8),dimension(Ns,Ns)      :: vdm_natural
  real(8),dimension(Ns)            :: vdm_diag
  !

  !HERE write the GZ EQUATIONS OF MOTION
  if(Nsys.ne.nDynamics) stop "wrong dimensions in the GZ_equations_of_motion"
  !

  call dynamicalVector_2_wfMatrix_superc(y,slater_,gzproj)
  slater(1:2,:,:,:) = slater_(1:2,:,:,:)
  slater(3,:,:,:) = zero
  do is=1,Ns
     do js=1,Ns
        if(is.eq.js) slater(3,is,js,:) = 1.d0
        slater(3,is,js,:) = slater(3,is,js,:) - slater_(1,js,is,:)
     end do
  end do
  !
  do is=1,Ns
     do js=1,Ns
        vdm_natural(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
     end do
     vdm_diag(is) = dreal(vdm_natural(is,is))
  end do
  !
  Rhop = hopping_renormalization_normal(gzproj,vdm_diag)
  Qhop = hopping_renormalization_anomalous(gzproj,vdm_diag)
  do is=1,Ns
     do js=1,Ns
        Rhop_hc(is,js) = conjg(Rhop(js,is))
        Qhop_hc(is,js) = conjg(Qhop(js,is))
     end do
  end do
  !
  slater_dot=zero
  gzproj_dot=zero
  !
  it = t2it(time,tstep*0.5d0)
  !
  slater_derivatives = zero
  slater_dot = zero

  do ik=1,Lk
     call get_Hk_t(Hk,ik,time)
     !
     tRR = matmul(Hk,Rhop)
     tRR = matmul(Rhop_hc,tRR)
     !
     tRQ = matmul(Hk,Qhop)
     tRQ = matmul(Rhop_hc,tRQ)
     !
     tQR = matmul(Hk,Rhop)
     tQR = matmul(Qhop_hc,tQR)
     !
     tQQ = matmul(Hk,Qhop)
     tQQ = matmul(Qhop_hc,tQQ)
     !+-> to speed up a bit I've to work out on these matrix products <-+!
     do is=1,Ns
        do js=1,Ns           
           slater_dot(1,is,js,ik) = zero
           slater_dot(2,is,js,ik) = zero
           do ks=1,Ns
              !tRR_{js,ks}SL1(is,ks) - tRR_{ks,is}SL1(ks,js) 
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRR(js,ks)*slater(1,is,ks,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) - tRR(ks,is)*slater(1,ks,js,ik) !**
              !tRQ_{js,ks}SL2(is,ks) + tRQ_{ks,js}SL2(ks,is)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRQ(js,ks)*slater(2,is,ks,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRQ(ks,js)*slater(2,ks,is,ik) !**
              !tQR_{is,ks}CONJG_SL2(js,ks) + tQR_{ks,is}CONJG_SL2(ks,js)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQR(is,ks)*conjg(slater(2,js,ks,ik))
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQR(ks,is)*conjg(slater(2,ks,js,ik)) !**
              !tQQ_{is,ks}SL1(ks,js) - tQQ_{ks,js}SL1(is,ks)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQQ(is,ks)*slater(1,ks,js,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) - tQQ(ks,js)*slater(1,is,ks,ik) !**


              !tRR_{ks,js}SL2(ks,is) - tRR_{ks,is}SL2(ks,js)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tRR(ks,js)*slater(2,ks,is,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tRR(ks,is)*slater(2,ks,js,ik) !**                    
              !tQR_{js,ks}SL1(is,ks) - tQR_{ks,js}SL1(is,ks)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQR(js,ks)*slater(1,is,ks,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQR(ks,js)*slater(1,is,ks,ik) !**
              !tQR_{is,ks}SL3(ks,js) - tQR_{ks,is}SL3(ks,js)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQR(is,ks)*slater(3,ks,js,ik) 
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQR(ks,is)*slater(3,ks,js,ik) !**
              !tQQ_{is,ks}SL2(ks,js) - tQQ_{js,ks}SL2(ks,is)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQQ(is,ks)*slater(2,ks,js,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQQ(js,ks)*slater(2,ks,is,ik) !**
              !
           end do
           !tmp(1:Ns,1:Ns) = <c^+ c>
           tmp(is,js) = slater(1,is,js,ik)     
           !tmp(1:Ns,Ns+1:2*Ns) = <c^+ c^+>
           tmp(is,js+Ns) = slater(2,is,js,ik)     
           !tmp(Ns+1:2*Ns,1:Ns) = <c c>
           tmp(is+Ns,js) = conjg(slater(2,js,is,ik))     
           !tmp(Ns+1:2*Ns,Ns+1:2*Ns) = <c c^+>
           tmp(is+Ns,js+Ns) = slater(3,is,js,ik)     
           !
        end do
     end do
     !       
     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Rhop(kks,ks))*Hk(kks,is)*tmp(ks,js)*wtk(ik)    
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Qhop(kks,ks))*Hk(kks,is)*tmp(ks+Ns,js)*wtk(ik)
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Rhop(kks,ks))*Hk(kks,is)*tmp(ks,js+Ns)*wtk(ik) 
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Qhop(kks,ks))*Hk(kks,is)*tmp(ks+Ns,js+Ns)*wtk(ik)
                 !
              end do
           end do
        end do
     end do
     !
  end do
  !
  slater_dot = -xi*slater_dot
  !
   !+- create HLOC
  Uloc=Uloc_t(:,it)
  Ust =Ust_t(it)
  Jh=Jh_t(it)
  Jsf=Jsf_t(it)
  Jph=Jph_t(it)
  eLevels = eLevels_t(:,it)
  call get_local_hamiltonian_trace(eLevels)  
  call build_neqH_GZproj_superc(Hproj,slater_derivatives,Rhop,Qhop,vdm_diag)
  !
  do iphi=1,Nphi
     gzproj_dot(iphi) = zero
     do jphi=1,Nphi
        gzproj_dot(iphi) = gzproj_dot(iphi) - xi * Hproj(iphi,jphi)*gzproj(jphi)
     end do
  end do
  !
  
  call wfMatrix_superc_2_dynamicalVector(slater_dot,gzproj_dot,f)
  !
end function gz_equations_of_motion_superc





function gz_equations_of_motion_superc_sp(time,y,Nsys) result(f)
  implicit none
  !inputs
  integer                          :: Nsys ! nr of equations
  real(8)                          :: time ! time variable                                                                                                                                      
  complex(8),dimension(Nsys)       :: y    ! argument array                                                                                                                                     
  complex(8),dimension(Nsys)       :: f    ! result 
  !
  complex(8),dimension(2,Ns,Ns,Lk) :: slater_,slater_dot
  complex(8),dimension(3,Ns,Ns,Lk) :: slater
  complex(8),dimension(Nphi)       :: gzproj,gzproj_dot
  complex(8),dimension(Nphi,Nphi)  :: Hproj
  type(sparse_matrix_csr_z)        :: htmp
  complex(8)                       :: xtmp
  real(8)                          :: xtmp_
  integer                          :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  complex(8),dimension(Ns,Ns)      :: tmpHk
  complex(8),dimension(Ns,Ns)      :: Hk
  complex(8),dimension(2*Ns,2*Ns)  :: tmp
  complex(8),dimension(2,Ns,Ns)    :: slater_derivatives
  complex(8),dimension(Ns,Ns)      :: Rhop,Qhop,Rhop_hc,Qhop_hc
  complex(8),dimension(Ns,Ns)      :: tRR,tRQ,tQR,tQQ
  complex(8),dimension(Ns,Ns)      :: vdm_natural
  real(8),dimension(Ns)            :: vdm_diag,n0
  !

  !HERE write the GZ EQUATIONS OF MOTION
  if(Nsys.ne.nDynamics) stop "wrong dimensions in the GZ_equations_of_motion"
  !

  call dynamicalVector_2_wfMatrix_superc(y,slater_,gzproj)
  slater(1:2,:,:,:) = slater_(1:2,:,:,:)
  slater(3,:,:,:) = zero
  do is=1,Ns
     do js=1,Ns
        if(is.eq.js) slater(3,is,js,:) = 1.d0
        slater(3,is,js,:) = slater(3,is,js,:) - slater_(1,js,is,:)
     end do
  end do
  !
  do is=1,Ns
     do js=1,Ns
        vdm_natural(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
     end do
     vdm_diag(is) = dreal(vdm_natural(is,is))
  end do
  !
  Rhop = hopping_renormalization_normal_sp(gzproj,vdm_diag)
  Qhop = hopping_renormalization_anomalous_sp(gzproj,vdm_diag)
  n0=vdm_diag
  do is=1,Ns
     do js=1,Ns
        Rhop_hc(is,js) = conjg(Rhop(js,is))
        Qhop_hc(is,js) = conjg(Qhop(js,is))
     end do
  end do
  !
  slater_dot=zero
  gzproj_dot=zero
  !
  it = t2it(time,tstep*0.5d0)
  !
  slater_derivatives = zero
  slater_dot = zero
  !
  do ik=1,Lk
     call get_Hk_t(Hk,ik,time)
     !
     tRR = matmul(Hk,Rhop)
     tRR = matmul(Rhop_hc,tRR)
     !
     tRQ = matmul(Hk,Qhop)
     tRQ = matmul(Rhop_hc,tRQ)
     !
     tQR = matmul(Hk,Rhop)
     tQR = matmul(Qhop_hc,tQR)
     !
     tQQ = matmul(Hk,Qhop)
     tQQ = matmul(Qhop_hc,tQQ)
     !+-> to speed up a bit I've to work out on these matrix products <-+!
     do is=1,Ns
        do js=1,Ns           
           slater_dot(1,is,js,ik) = zero
           slater_dot(2,is,js,ik) = zero
           do ks=1,Ns
              !tRR_{js,ks}SL1(is,ks) - tRR_{ks,is}SL1(ks,js) 
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRR(js,ks)*slater(1,is,ks,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) - tRR(ks,is)*slater(1,ks,js,ik) !**
              !tRQ_{js,ks}SL2(is,ks) + tRQ_{ks,js}SL2(ks,is)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRQ(js,ks)*slater(2,is,ks,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRQ(ks,js)*slater(2,ks,is,ik) !**
              !tQR_{is,ks}CONJG_SL2(js,ks) + tQR_{ks,is}CONJG_SL2(ks,js)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQR(is,ks)*conjg(slater(2,js,ks,ik))
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQR(ks,is)*conjg(slater(2,ks,js,ik)) !**
              !tQQ_{is,ks}SL1(ks,js) - tQQ_{ks,js}SL1(is,ks)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQQ(is,ks)*slater(1,ks,js,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) - tQQ(ks,js)*slater(1,is,ks,ik) !**


              !tRR_{ks,js}SL2(ks,is) - tRR_{ks,is}SL2(ks,js)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tRR(ks,js)*slater(2,ks,is,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tRR(ks,is)*slater(2,ks,js,ik) !**                    
              !tQR_{js,ks}SL1(is,ks) - tQR_{ks,js}SL1(is,ks)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQR(js,ks)*slater(1,is,ks,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQR(ks,js)*slater(1,is,ks,ik) !**
              !tQR_{is,ks}SL3(ks,js) - tQR_{ks,is}SL3(ks,js)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQR(is,ks)*slater(3,ks,js,ik) 
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQR(ks,is)*slater(3,ks,js,ik) !**
              !tQQ_{is,ks}SL2(ks,js) - tQQ_{js,ks}SL2(ks,is)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQQ(is,ks)*slater(2,ks,js,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQQ(js,ks)*slater(2,ks,is,ik) !**
              !
           end do
           !tmp(1:Ns,1:Ns) = <c^+ c>
           tmp(is,js) = slater(1,is,js,ik)     
           !tmp(1:Ns,Ns+1:2*Ns) = <c^+ c^+>
           tmp(is,js+Ns) = slater(2,is,js,ik)     
           !tmp(Ns+1:2*Ns,1:Ns) = <c c>
           tmp(is+Ns,js) = conjg(slater(2,js,is,ik))     
           !tmp(Ns+1:2*Ns,Ns+1:2*Ns) = <c c^+>
           tmp(is+Ns,js+Ns) = slater(3,is,js,ik)     
           !
        end do
     end do
     !       
     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Rhop(kks,ks))*Hk(kks,is)*tmp(ks,js)*wtk(ik)    
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Qhop(kks,ks))*Hk(kks,is)*tmp(ks+Ns,js)*wtk(ik)
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Rhop(kks,ks))*Hk(kks,is)*tmp(ks,js+Ns)*wtk(ik) 
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Qhop(kks,ks))*Hk(kks,is)*tmp(ks+Ns,js+Ns)*wtk(ik)
                 !
              end do
           end do
        end do
     end do
     !
  end do
  !
  slater_dot = -xi*slater_dot
  !
   !+- create HLOC
  Uloc=Uloc_t(:,it)
  Ust =Ust_t(it)
  Jh=Jh_t(it)
  Jsf=Jsf_t(it)
  Jph=Jph_t(it)
  eLevels = eLevels_t(:,it)
  ! call get_local_hamiltonian_trace(eLevels)  
  ! call build_neqH_GZproj_superc(Hproj,slater_derivatives,Rhop,Qhop,vdm_diag)
  ! !
  ! do iphi=1,Nphi
  !    gzproj_dot(iphi) = zero
  !    do jphi=1,Nphi
  !       gzproj_dot(iphi) = gzproj_dot(iphi) - xi * Hproj(iphi,jphi)*gzproj(jphi)
  !    end do
  ! end do
  !
  gzproj_dot = get_local_hamiltonian_HLOCphi(gzproj,eLevels)
  do is=1,Ns
     do js=1,Ns
        !
        xtmp=slater_derivatives(1,is,js)/sqrt(n0(js)*(1.d0-n0(js)))
        xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
        if(xtmp_.gt.1.d-10) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Rhop(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(slater_derivatives(1,is,js))/sqrt(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Rhop_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !htmp = htmp + slater_derivatives(1,is,js)*phi_traces_basis_Rhop(is,js,iphi,jphi)/sqrt(n0(js)*(1.d0-n0(js)))
        !htmp = htmp + conjg(slater_derivatives(1,is,js))*phi_traces_basis_Rhop_hc(is,js,iphi,jphi)/sqrt(n0(js)*(1.d0-n0(js)))
        xtmp=slater_derivatives(1,is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
        xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
        if(xtmp_.gt.1.d-10) then           
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)        
           xtmp=conjg(slater_derivatives(1,is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !htmp = htmp + slater_derivatives(1,is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,iphi,jphi)
        !htmp = htmp + conjg(slater_derivatives(1,is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,iphi,jphi) !hc conjgate capra!        
        xtmp=slater_derivatives(2,is,js)/sqrt(n0(js)*(1.d0-n0(js)))
        xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
        if(xtmp_.gt.1.d-10) then           
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Qhop(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(slater_derivatives(2,is,js))/sqrt(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Qhop_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !htmp = htmp + slater_derivatives(2,is,js)*phi_traces_basis_Qhop(is,js,iphi,jphi)/sqrt(n0(js)*(1.d0-n0(js)))
        !htmp = htmp + conjg(slater_derivatives(2,is,js))*phi_traces_basis_Qhop_hc(is,js,iphi,jphi)/sqrt(n0(js)*(1.d0-n0(js)))
        xtmp=slater_derivatives(2,is,js)*Qhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
        xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
        if(xtmp_.gt.1.d-10) then           
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(slater_derivatives(2,is,js)*Qhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !htmp = htmp + slater_derivatives(2,is,js)*Qhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,iphi,jphi)
        !htmp = htmp + conjg(slater_derivatives(2,is,js)*Qhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,iphi,jphi) !+-> capra!
     end do
  end do
  !  
  gzproj_dot = -xi*gzproj_dot

  
  call wfMatrix_superc_2_dynamicalVector(slater_dot,gzproj_dot,f)
  !
end function gz_equations_of_motion_superc_sp








function gz_equations_of_motion_superc_lgr(time,y,Nsys) result(f)
  implicit none
  !inputs
  integer                          :: Nsys ! nr of equations
  real(8)                          :: time ! time variable                                                                                                                                      
  complex(8),dimension(Nsys)       :: y    ! argument array                                                                                                                                     
  complex(8),dimension(Nsys)       :: f    ! result 
  !
  complex(8),dimension(2,Ns,Ns,Lk) :: slater_,slater_dot
  complex(8),dimension(3,Ns,Ns,Lk) :: slater
  complex(8),dimension(Nphi)       :: gzproj,gzproj_dot
  complex(8),dimension(Nphi,Nphi)  :: Hproj
  integer                          :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  complex(8),dimension(Ns,Ns)      :: tmpHk
  complex(8),dimension(Ns,Ns)      :: Hk
  complex(8),dimension(2*Ns,2*Ns)  :: tmp
  complex(8),dimension(2,Ns,Ns)    :: slater_derivatives
  complex(8),dimension(Ns,Ns)      :: Rhop,Qhop,Rhop_hc,Qhop_hc
  complex(8),dimension(Ns,Ns)      :: tRR,tRQ,tQR,tQQ
  complex(8),dimension(Ns,Ns)      :: vdm_natural
  complex(8)                       :: htmp
  real(8),dimension(Ns)            :: vdm_diag,n0
  real(8) :: test_slater
  !

  !HERE write the GZ EQUATIONS OF MOTION
  if(Nsys.ne.nDynamics) stop "wrong dimensions in the GZ_equations_of_motion"
  !

  call dynamicalVector_2_wfMatrix_superc(y,slater_,gzproj)
  slater(1:2,:,:,:) = slater_(1:2,:,:,:)
  slater(3,:,:,:) = zero
  do is=1,Ns
     do js=1,Ns
        if(is.eq.js) slater(3,is,js,:) = 1.d0
        slater(3,is,js,:) = slater(3,is,js,:) - slater_(1,js,is,:)
     end do
  end do
  !
  do is=1,Ns
     do js=1,Ns
        vdm_natural(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
     end do
     vdm_diag(is) = dreal(vdm_natural(is,is))
  end do
  !
  Rhop = hopping_renormalization_normal(gzproj,vdm_diag)
  Qhop = hopping_renormalization_anomalous(gzproj,vdm_diag)
  do is=1,Ns
     do js=1,Ns
        Rhop_hc(is,js) = conjg(Rhop(js,is))
        Qhop_hc(is,js) = conjg(Qhop(js,is))
     end do
  end do
  !
  slater_dot=zero
  gzproj_dot=zero
  !
  it = t2it(time,tstep*0.5d0)
  !
  slater_derivatives = zero
  slater_dot = zero
  
  do ik=1,Lk
     call get_Hk_t(Hk,ik,time)
     !
     tRR = matmul(Hk,Rhop)
     tRR = matmul(Rhop_hc,tRR)
     !
     tRQ = matmul(Hk,Qhop)
     tRQ = matmul(Rhop_hc,tRQ)
     !
     tQR = matmul(Hk,Rhop)
     tQR = matmul(Qhop_hc,tQR)
     !+- add_lgr_multipliers
     do is=1,Ns
        do js=1,Ns
           tRQ(is,js) = tRQ(is,js) + neq_lgr(1,is,js)
           tQR(is,js) = tQR(is,js) + conjg(neq_lgr(1,js,is))
        end do
     end do
     !
     tQQ = matmul(Hk,Qhop)
     tQQ = matmul(Qhop_hc,tQQ)
     !
     do is=1,Ns
        do js=1,Ns           
           slater_dot(1,is,js,ik) = zero
           slater_dot(2,is,js,ik) = zero
           do ks=1,Ns
              !tRR_{js,ks}SL1(is,ks) - tRR_{ks,is}SL1(ks,js) 
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRR(js,ks)*slater(1,is,ks,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) - tRR(ks,is)*slater(1,ks,js,ik) !**
              !tRQ_{js,ks}SL2(is,ks) + tRQ_{ks,js}SL2(ks,is)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRQ(js,ks)*slater(2,is,ks,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRQ(ks,js)*slater(2,ks,is,ik) !**
              !tQR_{is,ks}CONJG_SL2(js,ks) + tQR_{ks,is}CONJG_SL2(ks,js)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQR(is,ks)*conjg(slater(2,js,ks,ik))
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQR(ks,is)*conjg(slater(2,ks,js,ik)) !**
              !tQQ_{is,ks}SL1(ks,js) - tQQ_{ks,js}SL1(is,ks)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQQ(is,ks)*slater(1,ks,js,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) - tQQ(ks,js)*slater(1,is,ks,ik) !**


              !tRR_{ks,js}SL2(ks,is) - tRR_{ks,is}SL2(ks,js)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tRR(ks,js)*slater(2,ks,is,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tRR(ks,is)*slater(2,ks,js,ik) !**                    
              !tQR_{js,ks}SL1(is,ks) - tQR_{ks,js}SL1(is,ks)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQR(js,ks)*slater(1,is,ks,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQR(ks,js)*slater(1,is,ks,ik) !**
              !tQR_{is,ks}SL3(ks,js) - tQR_{ks,is}SL3(ks,js)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQR(is,ks)*slater(3,ks,js,ik) 
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQR(ks,is)*slater(3,ks,js,ik) !**
              !tQQ_{is,ks}SL2(ks,js) - tQQ_{js,ks}SL2(ks,is)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQQ(is,ks)*slater(2,ks,js,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQQ(js,ks)*slater(2,ks,is,ik) !**
              !
           end do
           !tmp(1:Ns,1:Ns) = <c^+ c>
           tmp(is,js) = slater(1,is,js,ik)     
           !tmp(1:Ns,Ns+1:2*Ns) = <c^+ c^+>
           tmp(is,js+Ns) = slater(2,is,js,ik)     
           !tmp(Ns+1:2*Ns,1:Ns) = <c c>
           tmp(is+Ns,js) = conjg(slater(2,js,is,ik))     
           !tmp(Ns+1:2*Ns,Ns+1:2*Ns) = <c c^+>
           tmp(is+Ns,js+Ns) = slater(3,is,js,ik)     
           !
        end do
     end do
     !  

     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Rhop(kks,ks))*Hk(kks,is)*tmp(ks,js)*wtk(ik)    
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Qhop(kks,ks))*Hk(kks,is)*tmp(ks+Ns,js)*wtk(ik)
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Rhop(kks,ks))*Hk(kks,is)*tmp(ks,js+Ns)*wtk(ik) 
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Qhop(kks,ks))*Hk(kks,is)*tmp(ks+Ns,js+Ns)*wtk(ik)
                 !
              end do
           end do
        end do
     end do
     !
  end do
  !
  slater_dot = -xi*slater_dot
  !

  !+- create HLOC
  Uloc=Uloc_t(:,it)
  Ust =Ust_t(it)
  Jh=Jh_t(it)
  Jsf=Jsf_t(it)
  Jph=Jph_t(it)
  eLevels = eLevels_t(:,it)
  call get_local_hamiltonian_trace(eLevels)  
  call build_neqH_GZproj_superc_lgr(Hproj,slater_derivatives,Rhop,Qhop,vdm_diag)
  !
  do iphi=1,Nphi
     !
     gzproj_dot(iphi) = zero
     !
     do jphi=1,Nphi
        !
        gzproj_dot(iphi) = gzproj_dot(iphi) + Hproj(iphi,jphi) *gzproj(jphi)
     end do
  end do
  gzproj_dot = -xi*gzproj_dot
  ! 
  call wfMatrix_superc_2_dynamicalVector(slater_dot,gzproj_dot,f)
  !
end function gz_equations_of_motion_superc_lgr













function gz_equations_of_motion_superc_lgr_sp(time,y,Nsys) result(f)
  implicit none
  !inputs
  integer                          :: Nsys ! nr of equations
  real(8)                          :: time ! time variable                                                                                                                                      
  complex(8),dimension(Nsys)       :: y    ! argument array                                                                                                                                     
  complex(8),dimension(Nsys)       :: f    ! result 
  !
  complex(8),dimension(2,Ns,Ns,Lk) :: slater_,slater_dot
  complex(8),dimension(3,Ns,Ns,Lk) :: slater
  complex(8),dimension(Nphi)       :: gzproj,gzproj_dot
  complex(8),dimension(Nphi,Nphi)  :: Hproj
  type(sparse_matrix_csr_z)          :: htmp
  complex(8)                       :: xtmp
  real(8)                          :: xtmp_
  integer                          :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  complex(8),dimension(Ns,Ns)      :: tmpHk
  complex(8),dimension(Ns,Ns)      :: Hk
  complex(8),dimension(2*Ns,2*Ns)  :: tmp
  complex(8),dimension(2,Ns,Ns)    :: slater_derivatives
  complex(8),dimension(Ns,Ns)      :: Rhop,Qhop,Rhop_hc,Qhop_hc
  complex(8),dimension(Ns,Ns)      :: tRR,tRQ,tQR,tQQ
  complex(8),dimension(Ns,Ns)      :: vdm_natural
  real(8),dimension(Ns)            :: vdm_diag,n0
  real(8) :: test_slater
  !

  !HERE write the GZ EQUATIONS OF MOTION
  if(Nsys.ne.nDynamics) stop "wrong dimensions in the GZ_equations_of_motion"
  !

  call dynamicalVector_2_wfMatrix_superc(y,slater_,gzproj)
  slater(1:2,:,:,:) = slater_(1:2,:,:,:)
  slater(3,:,:,:) = zero
  do is=1,Ns
     do js=1,Ns
        if(is.eq.js) slater(3,is,js,:) = 1.d0
        slater(3,is,js,:) = slater(3,is,js,:) - slater_(1,js,is,:)
     end do
  end do
  !
  do is=1,Ns
     do js=1,Ns
        vdm_natural(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
     end do
     vdm_diag(is) = dreal(vdm_natural(is,is))
  end do
  !
  ! Rhop = hopping_renormalization_normal(gzproj,vdm_diag)
  ! Qhop = hopping_renormalization_anomalous(gzproj,vdm_diag)

  Rhop = hopping_renormalization_normal_sp(gzproj,vdm_diag)
  Qhop = hopping_renormalization_anomalous_sp(gzproj,vdm_diag)
  n0=vdm_diag

  do is=1,Ns
     do js=1,Ns
        Rhop_hc(is,js) = conjg(Rhop(js,is))
        Qhop_hc(is,js) = conjg(Qhop(js,is))
     end do
  end do
  !
  slater_dot=zero
  gzproj_dot=zero
  !
  it = t2it(time,tstep*0.5d0)
  !
  slater_derivatives = zero
  slater_dot = zero
  
  do ik=1,Lk
     call get_Hk_t(Hk,ik,time)
     !
     tRR = matmul(Hk,Rhop)
     tRR = matmul(Rhop_hc,tRR)
     !
     tRQ = matmul(Hk,Qhop)
     tRQ = matmul(Rhop_hc,tRQ)
     !
     tQR = matmul(Hk,Rhop)
     tQR = matmul(Qhop_hc,tQR)
     !+- add_lgr_multipliers
     do is=1,Ns
        do js=1,Ns
           tRQ(is,js) = tRQ(is,js) + neq_lgr(1,is,js)
           tQR(is,js) = tQR(is,js) + conjg(neq_lgr(1,js,is))
        end do
     end do
     !
     tQQ = matmul(Hk,Qhop)
     tQQ = matmul(Qhop_hc,tQQ)
     !
     do is=1,Ns
        do js=1,Ns           
           slater_dot(1,is,js,ik) = zero
           slater_dot(2,is,js,ik) = zero
           do ks=1,Ns
              !tRR_{js,ks}SL1(is,ks) - tRR_{ks,is}SL1(ks,js) 
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRR(js,ks)*slater(1,is,ks,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) - tRR(ks,is)*slater(1,ks,js,ik) !**
              !tRQ_{js,ks}SL2(is,ks) + tRQ_{ks,js}SL2(ks,is)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRQ(js,ks)*slater(2,is,ks,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tRQ(ks,js)*slater(2,ks,is,ik) !**
              !tQR_{is,ks}CONJG_SL2(js,ks) + tQR_{ks,is}CONJG_SL2(ks,js)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQR(is,ks)*conjg(slater(2,js,ks,ik))
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQR(ks,is)*conjg(slater(2,ks,js,ik)) !**
              !tQQ_{is,ks}SL1(ks,js) - tQQ_{ks,js}SL1(is,ks)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + tQQ(is,ks)*slater(1,ks,js,ik)
              slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) - tQQ(ks,js)*slater(1,is,ks,ik) !**


              !tRR_{ks,js}SL2(ks,is) - tRR_{ks,is}SL2(ks,js)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tRR(ks,js)*slater(2,ks,is,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tRR(ks,is)*slater(2,ks,js,ik) !**                    
              !tQR_{js,ks}SL1(is,ks) - tQR_{ks,js}SL1(is,ks)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQR(js,ks)*slater(1,is,ks,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQR(ks,js)*slater(1,is,ks,ik) !**
              !tQR_{is,ks}SL3(ks,js) - tQR_{ks,is}SL3(ks,js)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQR(is,ks)*slater(3,ks,js,ik) 
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQR(ks,is)*slater(3,ks,js,ik) !**
              !tQQ_{is,ks}SL2(ks,js) - tQQ_{js,ks}SL2(ks,is)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQQ(is,ks)*slater(2,ks,js,ik)
              slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQQ(js,ks)*slater(2,ks,is,ik) !**
              !
           end do
           !tmp(1:Ns,1:Ns) = <c^+ c>
           tmp(is,js) = slater(1,is,js,ik)     
           !tmp(1:Ns,Ns+1:2*Ns) = <c^+ c^+>
           tmp(is,js+Ns) = slater(2,is,js,ik)     
           !tmp(Ns+1:2*Ns,1:Ns) = <c c>
           tmp(is+Ns,js) = conjg(slater(2,js,is,ik))     
           !tmp(Ns+1:2*Ns,Ns+1:2*Ns) = <c c^+>
           tmp(is+Ns,js+Ns) = slater(3,is,js,ik)     
           !
        end do
     end do
     !  

     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Rhop(kks,ks))*Hk(kks,is)*tmp(ks,js)*wtk(ik)    
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Qhop(kks,ks))*Hk(kks,is)*tmp(ks+Ns,js)*wtk(ik)
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Rhop(kks,ks))*Hk(kks,is)*tmp(ks,js+Ns)*wtk(ik) 
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Qhop(kks,ks))*Hk(kks,is)*tmp(ks+Ns,js+Ns)*wtk(ik)
                 !
              end do
           end do
        end do
     end do
     !
  end do
  !
  slater_dot = -xi*slater_dot
  !

  !+- create HLOC
  Uloc=Uloc_t(:,it)
  Ust =Ust_t(it)
  Jh=Jh_t(it)
  Jsf=Jsf_t(it)
  Jph=Jph_t(it)
  eLevels = eLevels_t(:,it)
  ! call get_local_hamiltonian_trace(eLevels)  
  ! call build_neqH_GZproj_superc_lgr(Hproj,slater_derivatives,Rhop,Qhop,vdm_diag)
  !
  ! htmp = sp_scalar_matrix_csr(phi_spTraces_basis_local_dens(istate,istate),(atomic_energy_levels(istate)-xmu))
  ! gzproj_out = gzproj_out + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj_in)
  ! call sp_delete_matrix(htmp)

  gzproj_dot = get_local_hamiltonian_HLOCphi(gzproj,eLevels)
  do is=1,Ns
     do js=1,Ns
        !
        xtmp=slater_derivatives(1,is,js)/sqrt(n0(js)*(1.d0-n0(js)))
        xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
        if(xtmp_.gt.1.d-10) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Rhop(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(slater_derivatives(1,is,js))/sqrt(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Rhop_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !htmp = htmp + slater_derivatives(1,is,js)*phi_traces_basis_Rhop(is,js,iphi,jphi)/sqrt(n0(js)*(1.d0-n0(js)))
        !htmp = htmp + conjg(slater_derivatives(1,is,js))*phi_traces_basis_Rhop_hc(is,js,iphi,jphi)/sqrt(n0(js)*(1.d0-n0(js)))
        xtmp=slater_derivatives(1,is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
        xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
        if(xtmp_.gt.1.d-10) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)        
           xtmp=conjg(slater_derivatives(1,is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !htmp = htmp + slater_derivatives(1,is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,iphi,jphi)
        !htmp = htmp + conjg(slater_derivatives(1,is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,iphi,jphi) !hc conjgate capra!        
        xtmp=slater_derivatives(2,is,js)/sqrt(n0(js)*(1.d0-n0(js)))
        xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
        if(xtmp_.gt.1.d-10) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Qhop(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(slater_derivatives(2,is,js))/sqrt(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Qhop_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !htmp = htmp + slater_derivatives(2,is,js)*phi_traces_basis_Qhop(is,js,iphi,jphi)/sqrt(n0(js)*(1.d0-n0(js)))
        !htmp = htmp + conjg(slater_derivatives(2,is,js))*phi_traces_basis_Qhop_hc(is,js,iphi,jphi)/sqrt(n0(js)*(1.d0-n0(js)))
        xtmp=slater_derivatives(2,is,js)*Qhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
        xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
        if(xtmp_.gt.1.d-10) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(slater_derivatives(2,is,js)*Qhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !htmp = htmp + slater_derivatives(2,is,js)*Qhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,iphi,jphi)
        !htmp = htmp + conjg(slater_derivatives(2,is,js)*Qhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,iphi,jphi) !+-> capra!

        !+- > add time dependent lagrange multipliers <-!
        xtmp=neq_lgr(2,is,js)
        xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
        if(xtmp_.gt.1.d-10) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_anomalous(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(neq_lgr(2,is,js))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_anomalous_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !htmp = htmp + neq_lgr(2,is,js)*phi_traces_basis_dens_anomalous(is,js,iphi,jphi)
        !htmp = htmp + conjg(neq_lgr(2,is,js))*phi_traces_basis_dens_anomalous_hc(is,js,iphi,jphi)        
     end do
  end do
  !  
  gzproj_dot = -xi*gzproj_dot
  ! 
  call wfMatrix_superc_2_dynamicalVector(slater_dot,gzproj_dot,f)
  !
end function gz_equations_of_motion_superc_lgr_sp












subroutine build_neqH_GZproj(Hproj,slater_derivatives,Rhop,n0)
  complex(8),dimension(Nphi,Nphi) :: Hproj
  complex(8),dimension(Ns,Ns)     :: slater_derivatives,Rhop
  real(8),dimension(Ns)           :: n0 
  integer                         :: is,js
  Hproj=zero
  Hproj=phi_traces_basis_Hloc
  do is=1,Ns
     do js=1,Ns
        !
        Hproj = Hproj + slater_derivatives(is,js)*phi_traces_basis_Rhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
        Hproj = Hproj + conjg(slater_derivatives(is,js))*phi_traces_basis_Rhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))

        Hproj = Hproj + slater_derivatives(is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,:,:)
        Hproj = Hproj + conjg(slater_derivatives(is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,:,:)
        !
     end do
  end do
end subroutine build_neqH_GZproj




subroutine build_neqH_GZproj_superc(Hproj,slater_derivatives,Rhop,Qhop,n0)
  complex(8),dimension(Nphi,Nphi) :: Hproj
  complex(8),dimension(2,Ns,Ns)   :: slater_derivatives
  complex(8),dimension(Ns,Ns)     :: Rhop,Qhop
  real(8),dimension(Ns)           :: n0 
  integer                         :: is,js
  real(8)                         :: test_slater
  !
  Hproj=phi_traces_basis_Hloc
  do is=1,Ns
     do js=1,Ns
        !
        test_slater=conjg(slater_derivatives(1,is,js))*slater_derivatives(1,is,js)
        if(test_slater.gt.1.d-10) then
           Hproj = Hproj + slater_derivatives(1,is,js)*phi_traces_basis_Rhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
           Hproj = Hproj + conjg(slater_derivatives(1,is,js))*phi_traces_basis_Rhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
        !
           Hproj = Hproj + slater_derivatives(1,is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,:,:)
           Hproj = Hproj + conjg(slater_derivatives(1,is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,:,:)
        end if
        !
        !
        !
        test_slater=conjg(slater_derivatives(2,is,js))*slater_derivatives(2,is,js)
        if(test_slater.gt.1.d-10) then
           Hproj = Hproj + slater_derivatives(2,is,js)*phi_traces_basis_Qhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
           Hproj = Hproj + conjg(slater_derivatives(2,is,js))*phi_traces_basis_Qhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
           !
           Hproj = Hproj + slater_derivatives(2,is,js)*Qhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,:,:)
           Hproj = Hproj + conjg(slater_derivatives(2,is,js)*Qhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,:,:)
        end if
           !
     end do
  end do
end subroutine build_neqH_GZproj_superc



subroutine build_neqH_GZproj_superc_lgr(Hproj,slater_derivatives,Rhop,Qhop,n0)
  complex(8),dimension(Nphi,Nphi) :: Hproj
  complex(8),dimension(2,Ns,Ns)   :: slater_derivatives
  complex(8),dimension(Ns,Ns)     :: Rhop,Qhop
  real(8),dimension(Ns)           :: n0 
  integer                         :: is,js
  real(8)                         :: test_slater
  !
  Hproj=phi_traces_basis_Hloc
  do is=1,Ns
     do js=1,Ns
        !
        test_slater=sqrt(conjg(slater_derivatives(1,is,js))*slater_derivatives(1,is,js))
        if(test_slater.gt.1.d-10) then
           Hproj = Hproj + slater_derivatives(1,is,js)*phi_traces_basis_Rhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
           Hproj = Hproj + conjg(slater_derivatives(1,is,js))*phi_traces_basis_Rhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
           !
           Hproj = Hproj + slater_derivatives(1,is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,:,:)
           Hproj = Hproj + conjg(slater_derivatives(1,is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,:,:)
        end if
        test_slater=sqrt(conjg(slater_derivatives(2,is,js))*slater_derivatives(2,is,js))
        if(test_slater.gt.1.d-10) then
           !
           Hproj = Hproj + slater_derivatives(2,is,js)*phi_traces_basis_Qhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
           Hproj = Hproj + conjg(slater_derivatives(2,is,js))*phi_traces_basis_Qhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
           !
           Hproj = Hproj + slater_derivatives(2,is,js)*Qhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,:,:)
           Hproj = Hproj + conjg(slater_derivatives(2,is,js)*Qhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(js,js,:,:)
        end if
        Hproj = Hproj + neq_lgr(2,is,js)*phi_traces_basis_dens_anomalous(is,js,:,:)
        Hproj = Hproj + conjg(neq_lgr(2,is,js))*phi_traces_basis_dens_anomalous_hc(is,js,:,:)
        !
     end do
  end do
end subroutine build_neqH_GZproj_superc_lgr
