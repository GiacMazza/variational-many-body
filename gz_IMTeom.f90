!
function gz_imt_equations_of_motion(time,y,Nsys) result(f)
  implicit none
  !inputs                                                                                                                                   
  integer                                     :: Nsys ! nr of equations
  real(8)                                     :: time ! time variable
  complex(8),dimension(Nsys)                  :: y    ! argument array
  complex(8),dimension(Nsys)                  :: f    ! result 
  !
  complex(8),dimension(Ns,Ns,Lk)              :: slater,slater_dot,Hqp,Hqp_dot
  complex(8),dimension(Ns,Ns)                 :: Hk,tRR
  complex(8),dimension(Nphi)                  :: gzproj,gzproj_dot,tmp_gzproj_dot
  complex(8),dimension(Nphi,Nphi)             :: Hproj
  integer                                     :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  complex(8),dimension(Ns,Ns)                 :: tmpHk,Rhop,slater_derivatives,Rhop_hc
  real(8),dimension(Ns)                    :: tmp_eHk
  complex(8),dimension(Ns,Ns)                 :: vdm_natural
  real(8),dimension(Ns)                       :: vdm_diag,n0
  !
  type(sparse_matrix_csr_z)        :: htmp
  complex(8)                       :: xtmp
  !
  real(8) :: tmpU,tmpU_
  !
  it = imt_t2it(time,itstep*0.5d0)
  if(Nsys.ne.nDynamics) stop "wrong dimensions in the imt_GZ_equations_of_motion"
  call dynamicalVector_2_wfMatrix(y,Hqp,gzproj)
  gzproj_dot=zero
  Hqp_dot=zero
  !
  !+- get some shit
  do is=1,Ns
     do js=1,Ns
        vdm_natural(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
     end do
     vdm_diag(is) = dreal(vdm_natural(is,is))
  end do
  n0=vdm_diag
  Rhop = hopping_renormalization_normal(gzproj,vdm_diag)
  do is=1,Ns
     do js=1,Ns
        Rhop_hc(is,js) = conjg(Rhop(js,is))
     end do
  end do
  !
  !
  !
  slater_derivatives=zero
  do ik=1,Lk
     !
     Hk = matmul(Hk_tb(:,:,ik),Rhop)
     Hk = matmul(Rhop_hc,Hk)
     !
     Hk = Hk + imt_lgrNC
     !
     tmpHk = Hqp(:,:,ik) !+-> temporary H used to compute occupation of the slater states
     call matrix_diagonalize(tmpHk,tmp_eHk)       
     !
     do is=1,Ns
        do js=1,Ns
           slater(is,js,ik) = zero
           do ks=1,Ns
              slater(is,js,ik) = slater(is,js,ik) + conjg(tmpHk(is,ks))*tmpHk(js,ks)*fermi(tmp_eHk(ks),1.d0)
           end do
           !
           Hqp_dot(is,js,ik) = 2.0*Hk(is,js) 
           !
        end do
     end do
     Hk=Hk_tb(:,:,ik) !+- coglione!
     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns                 
                 slater_derivatives(is,js) = slater_derivatives(is,js) + conjg(Rhop(kks,ks))*Hk(kks,is)*slater(ks,js,ik)*wtk(ik)
              end do
           end do
        end do
     end do
  end do
  !
  eLevels = 0.d0
  !
  gzproj_dot = get_local_hamiltonian_HLOCphi(gzproj,eLevels)
  do is=1,Ns
     do js=1,Ns
        !
        xtmp=slater_derivatives(is,js)/sqrt(n0(js)*(1.d0-n0(js)))
        if(xtmp/=zero) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Rhop(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(xtmp)
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Rhop_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        xtmp=slater_derivatives(is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
        if(xtmp/=zero) then           
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(js,js),xtmp)   !js,js is right, capra!!
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)        
           xtmp=conjg(xtmp)
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !+- add lagrange parameters
        xtmp=-1.d0*imt_lgrNC(is,js)
        if(xtmp/=zero) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(xtmp)
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        ! !
        xtmp=-imt_lgr_local_dens
        if(xtmp/=zero.and.is.eq.js) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_local_dens(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !
     end do
  end do
  gzproj_dot=-gzproj_dot
  !
  !
  tmpU = 0.d0
  tmpU_= 0.d0
  do iphi=1,Nphi
     tmpU  = tmpU + gzproj(iphi)*conjg(gzproj(iphi))
     tmpU_ = tmpU_ + conjg(gzproj(iphi))*gzproj_dot(iphi)
  end do
  imt_lgrU = tmpU_/tmpU
  gzproj_dot = gzproj_dot - imt_lgrU*gzproj 
  !
  call wfMatrix_2_dynamicalVector(Hqp_dot,gzproj_dot,f)
  !
end function gz_imt_equations_of_motion




!+--- SUPERC ROUTINES ---+!
function gz_imt_equations_of_motion_superc(time,y,Nsys) result(f)
  implicit none
  !inputs                                                                                                                                   
  integer                          :: Nsys ! nr of equations
  real(8)                          :: time ! time variable
  complex(8),dimension(Nsys)       :: y    ! argument array
  complex(8),dimension(Nsys)       :: f    ! result 
  !
  complex(8),dimension(3,Ns,Ns,Lk) :: slater
  complex(8),dimension(3,Ns,Ns,Lk) :: Hqp,Hqp_dot
  complex(8),dimension(2,Ns,Ns) :: slater_derivatives
  complex(8),dimension(2*Ns,2*Ns)  :: Hk,tmpHqp,tmp
  complex(8),dimension(Nphi)       :: gzproj,gzproj_dot,tmp_gzproj_dot
  !
  integer                          :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  complex(8),dimension(Ns,Ns)      :: tmpHk,Rhop,Rhop_hc,Qhop,Qhop_hc,Hk0
  real(8),dimension(2*Ns)          :: tmp_eHk
  complex(8),dimension(Ns,Ns)      :: vdm_natural
  real(8),dimension(Ns)            :: vdm_diag,n0
  !
  type(sparse_matrix_csr_z)        :: htmp
  complex(8)                       :: xtmp,tmp_test
  !
  real(8)                          :: tmpU,tmpU_,nqp

  
  complex(8),dimension(2,Ns,Ns) :: tmp_lgr
  real(8) :: tmp_Estar
  complex(8),dimension(2,Ns,Ns,Lk) :: tmp_slater

  !
  it = imt_t2it(time,itstep*0.5d0)
  if(Nsys.ne.nDynamics) stop "wrong dimensions in the imt_GZ_equations_of_motion"
  !
  call imt_dynamicalVector_2_wfMatrix_superc(y,Hqp,gzproj)
  !
  gzproj_dot=zero
  Hqp_dot=zero
  !
  do is=1,Ns
     do js=1,Ns
        vdm_natural(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
     end do
     vdm_diag(is) = dreal(vdm_natural(is,is))
  end do
  !
  n0=vdm_diag
  Rhop=hopping_renormalization_normal_sp(gzproj,vdm_diag)
  Rhop_hc=conjg(transpose(Rhop))
  Qhop=hopping_renormalization_anomalous_sp(gzproj,vdm_diag)
  Qhop_hc=conjg(transpose(Qhop))
  !
  slater_derivatives=zero
  slater=zero
  do ik=1,Lk
     !
     tmpHk = matmul(Hk_tb(:,:,ik),Rhop)
     tmpHk = matmul(Rhop_hc,tmpHk)
     Hk(1:Ns,1:Ns) = tmpHk + imt_lgr_NC(1,:,:) !+ imt_lgr_NC(2,:,:)
     !     
     tmpHk = matmul(Hk_tb(:,:,ik),Qhop)
     tmpHk = matmul(Rhop_hc,tmpHk)
     Hk(1+Ns:2*Ns,1:Ns) = tmpHk + imt_lgr_AC(1,:,:) !+ 
     !
     tmpHk = matmul(Hk_tb(:,:,ik),Rhop)
     tmpHk = matmul(Qhop_hc,tmpHk)
     Hk(1:Ns,1+Ns:2*Ns) = tmpHk + conjg(transpose(imt_lgr_AC(1,:,:)))
     !
     tmpHk = matmul(Hk_tb(:,:,ik),Qhop)
     tmpHk = matmul(Qhop_hc,tmpHk)
     Hk(1+Ns:2*Ns,1+Ns:2*Ns) = tmpHk     
     !
     !
     Hqp_dot(1,:,:,ik) = 2.0*Hk(1:Ns,1:Ns)           
     Hqp_dot(2,:,:,ik) = 2.0*Hk(1:Ns,1+Ns:2*Ns)     !+- da dove viene questo meno qui??? -+!
!     Hqp_dot(4,:,:,ik) = conjg(transpose(Hqp_dot(2,:,:,ik)))!2.0*Hk(1+Ns:2*Ns,1:Ns)     
     Hqp_dot(3,:,:,ik) = 2.0*Hk(1+Ns:2*Ns,1+Ns:2*Ns) 
     !
     !
     tmpHqp(1:Ns,1:Ns) = Hqp(1,:,:,ik) 
     tmpHqp(1:Ns,1+Ns:2*Ns) = Hqp(2,:,:,ik) 
     tmpHqp(1+Ns:2*Ns,1:Ns) = conjg(transpose(Hqp(2,:,:,ik)))  ! Hqp(4,:,:,ik) !
     tmpHqp(1+Ns:2*Ns,1+Ns:2*Ns) = Hqp(3,:,:,ik)
     !
     !
     !
     !
     !write(*,*) ik,'before'
     !
     !
     !tmpHqp=Hk

     call matrix_diagonalize(tmpHqp,tmp_eHk)       

     ! if(ik.eq.10) then
     !    do is=1,2*Ns
     !       write(*,'(20F18.10)') dreal(tmpHqp(:,is))
     !    end do
     !    write(*,'(20F18.10)') 
     !    write(*,'(20F18.10)') 
     !    write(*,'(20F18.10)') 
     ! end if
     ! if(ik.eq.10) then
     !    do is=1,2*Ns
     !       write(*,'(20F18.10)') dreal(tmpHqp(:,is))
     !    end do
     !    write(*,'(20F18.10)') 
     !    write(*,'(20F18.10)') 
     !    write(*,'(20F18.10)') 
     !    ! do is=1,2*Ns
     !    !    do js=1,2*Ns
     !    !       tmp_test=0.d0
     !    !       do ks=1,2*Ns
     !    !          tmp_test = tmp_test + tmpHqp(is,ks)*conjg(tmpHqp(js,ks))
     !    !       end do
     !    !       write(*,*) tmp_test,is,js
     !    !    end do
     !    ! end do
     !    stop
     ! end if
     !write(*,*) ik,'after'
     !
     !
     do is=1,Ns
        do js=1,Ns
           !
           do ks=1,Ns
              !
              nqp = fermi(tmp_eHk(ks)-tmp_eHk(ks+Ns),1.d0)                
              !nqp = fermi(tmp_eHk(ks)-tmp_eHk(ks+Ns),beta)                
              !
              slater(1,is,js,ik) = slater(1,is,js,ik) + &    ! (*)
                   conjg(tmpHqp(is,ks))*tmpHqp(js,ks)*nqp + conjg(tmpHqp(is,ks+Ns))*tmpHqp(js,ks+Ns)*(1.d0-nqp)
              !
              slater(2,is,js,ik) = slater(2,is,js,ik) + &  ! (a)
                   conjg(tmpHqp(is+Ns,ks))*tmpHqp(js,ks)*nqp + conjg(tmpHqp(is+Ns,ks+Ns))*tmpHqp(js,ks+Ns)*(1.d0-nqp)
              !
              ! slater(2,is,js,ik) = slater(2,is,js,ik) + &    ! (b)
              !      conjg(tmpHqp(is,ks))*tmpHqp(js+Ns,ks)*nqp + conjg(tmpHqp(is,ks+Ns))*tmpHqp(js+Ns,ks+Ns)*(1.d0-nqp)
              !(a) and (b) not equivalent ???? !+- thy should be equivalent mmmm....
              !
           end do
           !
        end do
     end do
     !
     !do ik=1,Lk
     ! write(557,'(20F18.10)') dreal(slater(1,1,1,ik))
     ! write(558,'(20F18.10)') dreal(slater(2,1,2,ik))
     !end do
     do is=1,Ns
        do js=1,Ns
           if(is.eq.js)  slater(3,is,js,ik) = 1.d0
           slater(3,is,js,ik) = slater(3,is,js,ik) - slater(1,js,is,ik)           
           !           
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
     Hk0=Hk_tb(:,:,ik)      
     !+- slater derivatives
     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Rhop(kks,ks))*Hk0(kks,is)*tmp(ks,js)*wtk(ik)    
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Qhop(kks,ks))*Hk0(kks,is)*tmp(ks+Ns,js)*wtk(ik)
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Rhop(kks,ks))*Hk0(kks,is)*tmp(ks,js+Ns)*wtk(ik) 
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Qhop(kks,ks))*Hk0(kks,is)*tmp(ks+Ns,js+Ns)*wtk(ik)
                 !
              end do
           end do
        end do
     end do
     !+-
  end do
  !stop
  ! if(time==19.d0) then
  !    write(555,'(40F18.10)')
  !    write(555,'(40F18.10)')
  !    write(556,'(40F18.10)')
  !    write(556,'(40F18.10)')
  ! end if
  !
  eLevels = 0.d0
  !
  tmp_gzproj_dot=zero
  gzproj_dot = get_local_hamiltonian_HLOCphi(gzproj,eLevels)
  do is=1,Ns
     do js=1,Ns
        !
        xtmp=slater_derivatives(1,is,js)/sqrt(n0(js)*(1.d0-n0(js)))
        if(xtmp/=zero) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Rhop(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(slater_derivatives(1,is,js))/sqrt(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Rhop_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        xtmp=slater_derivatives(1,is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
        if(xtmp/=zero) then           
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)        
           xtmp=conjg(slater_derivatives(1,is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        xtmp=slater_derivatives(2,is,js)/sqrt(n0(js)*(1.d0-n0(js)))
        if(xtmp/=zero) then           
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Qhop(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(slater_derivatives(2,is,js))/sqrt(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Qhop_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        xtmp=slater_derivatives(2,is,js)*Qhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
        if(xtmp/=zero) then           
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(slater_derivatives(2,is,js)*Qhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(js,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !+- add lagrange parameters -+!
        xtmp=-1.d0*imt_lgr_NC(2,is,js)
        if(xtmp/=zero) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(xtmp)
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !
        ! xtmp=1.d0
        ! htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_anomalous(is,js),xtmp)
        ! write(660,'(20F18.10)') time,sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
        ! call sp_delete_matrix(htmp)
        ! stop
        !
        xtmp=1.d0*imt_lgr_AC(2,is,js)
        if(xtmp/=zero) then
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_anomalous(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
           xtmp=conjg(xtmp)
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_anomalous_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !
        xtmp=-imt_lgr_local_dens
        if(xtmp/=zero.and.is.eq.js) then           
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_local_dens(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !
     end do
  end do
  gzproj_dot=-gzproj_dot
  !
  tmpU = 0.d0
  tmpU_= 0.d0
  do iphi=1,Nphi
     tmpU  = tmpU + gzproj(iphi)*conjg(gzproj(iphi))
     tmpU_ = tmpU_ + conjg(gzproj(iphi))*gzproj_dot(iphi)
  end do
  imt_lgrU = tmpU_/tmpU
  gzproj_dot = gzproj_dot - imt_lgrU*gzproj 
  !
  call imt_wfMatrix_superc_2_dynamicalVector(Hqp_dot,gzproj_dot,f)
  !
end function gz_imt_equations_of_motion_superc



