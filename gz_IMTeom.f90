function gz_imt_equations_of_motion(time,y,Nsys) result(f)
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
  real(8),dimension(Ns)                    :: tmp_eHk
  complex(8),dimension(Ns,Ns)                 :: vdm_natural
  real(8),dimension(Ns)                       :: vdm_diag,n0
  !
  type(sparse_matrix_csr_z)        :: htmp
  complex(8)                       :: xtmp
  !
  !HERE write the GZ EQUATIONS OF MOTION

  !+- the main idea is:
  !   
  !   I have the integral up to time t0
  !   I approximate the remaining part of the integral as (time-t0)*f(time,y)  [this is the only tricky part]
  !   at the end I update the integral
  !
  it = imt_t2it(time,itstep*0.5d0)
  gzproj_dot=zero
  if(Nsys.ne.nDynamics) stop "wrong dimensions in the imt_GZ_equations_of_motion"
  gzproj = y
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
  slater_derivatives=zero
  do ik=1,Lk
     !
     tmpHk= gz_imt_qpH(ik,:,:) 
     Hk = matmul(Hk_tb(:,:,ik),Rhop)
     Hk = matmul(Rhop_hc,Hk)
     Hk = Hk + imt_lgrNC

     if(ik.eq.10) write(*,'(10F18.10)') tmpHk(1,1),Hk_tb(1,1,ik),time,imt_tstep
     tmpHk = tmpHk + 2.d0*(time-imt_tstep)*Hk !here add/or subtract the small increment integral !the two is important!!
     if(ik.eq.10) then
        write(*,*) tmpHk(1,1)
        write(*,*)
     end if

     !
     call matrix_diagonalize(tmpHk,tmp_eHk)       
     !
     do is=1,Ns
        do js=1,Ns
           slater(is,js,ik) = zero
           do ks=1,Ns
              slater(is,js,ik) = slater(is,js,ik) + conjg(tmpHk(is,ks))*tmpHk(js,ks)*fermi(tmp_eHk(ks),1.d0)
           end do
        end do
     end do
     !if(abs(time-999).lt.1.d-6) write(456,'(10F18.10)') slater(1,1,ik)
     !
     Hk=Hk-imt_lgrNC
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
  ! Uloc=Uloc_t(:,it)
  ! Ust =Ust_t(it)
  ! Jh=Jh_t(it)
  ! Jsf=Jsf_t(it)
  ! Jph=Jph_t(it)
  !  write(*,*) slater_derivatives
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
           !xtmp=conjg(slater_derivatives(is,js))/sqrt(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Rhop_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        xtmp=slater_derivatives(is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
        if(xtmp/=zero) then           
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)        
           xtmp=conjg(xtmp)
           !xtmp=conjg(slater_derivatives(is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
           htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(is,js),xtmp)
           gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
           call sp_delete_matrix(htmp)
        end if
        !+- add lagrange parameter
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
     end do
  end do
  f = -gzproj_dot
  !
end function gz_imt_equations_of_motion
!
function gz_imt_eom(time,y,Nsys) result(f)
  implicit none
  !inputs
  integer                                   :: Nsys ! nr of equations
  real(8)                                   :: time ! time variable
  complex(8),dimension(Nsys)                :: y    ! argument array
  complex(8),dimension(Nsys)                :: f    ! result 
  !
  !
  !
  real(8),dimension(:),allocatable          :: lgr,delta_out
  complex(8),dimension(:),allocatable       :: lgr_cmplx
  complex(8),dimension(Nphi)                :: gzproj
  integer                                   :: iter,Nopt,iphi
  integer                                   :: i,i0
  real(8)                                   :: delta
  real(8) :: tmpU,tmpU_
  !
  ! if(allocated(imt_lgrNC)) deallocate(imt_lgr)
  ! allocate(imt_lgrNC(Ns,Ns)); imt_lgrNC=zero
  !
  !HERE write the GZ EQUATIONS OF MOTION

  !   get the traces over the uncorrelated density distribution (as a function of the lagrange params)
  !   get the derivative of the Phi matrices, such that d/dt TR(Phi* Phi) = 0
  !   be in peace 

  ! modifications for the anomalous constriants: as before, with two more lgr params.



  if(Nsys.ne.nDynamics) stop "wrong dimensions in the GZ_equations_of_motion"
  !
  f = gz_imt_equations_of_motion(time,y,Nsys)
  
  gzproj = y
  !
  tmpU = 0.d0
  tmpU_= 0.d0
  do iphi=1,Nphi
     tmpU  = tmpU + gzproj(iphi)*conjg(gzproj(iphi))
     tmpU_ = tmpU_ + conjg(gzproj(iphi))*f(iphi)
  end do
  imt_lgrU = tmpU_/tmpU
  f = f - imt_lgrU*gzproj
  f=-f !?!?!?!
end function gz_imt_eom









subroutine update_qpH_integral(yold,ynew,lgrNC_old,lgrNC_new) 
  implicit none
  !inputs                                                                                                                                             
  complex(8),dimension(Nphi)                  :: yold
  complex(8),dimension(Nphi)                  :: ynew 
  complex(8),dimension(Nphi)                  :: gzproj
  complex(8),dimension(Ns,Ns)                 :: lgrNC_old,lgrNC_new
  !
  complex(8),dimension(Ns,Ns)                 :: Hk,tRR
  integer                                     :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  complex(8),dimension(Ns,Ns)                 :: tmpHk,Rhop,slater_derivatives,Rhop_hc
  real(8),dimension(Ns)                    :: tmp_eHk
  complex(8),dimension(Ns,Ns)                 :: vdm_natural
  real(8),dimension(Ns)                       :: vdm_diag,n0
  !
  type(sparse_matrix_csr_z)        :: htmp
  complex(8)                       :: xtmp
  !
  gzproj = yold
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
  do ik=1,Lk
     !
     Hk = matmul(Hk_tb(:,:,ik),Rhop)
     Hk = matmul(Rhop_hc,Hk)
     Hk = Hk + lgrNC_old
     gz_imt_qpH(ik,:,:) = gz_imt_qpH(ik,:,:) + 2.d0*Hk*itstep*0.5d0     
  end do



  gzproj = ynew
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
  do ik=1,Lk
     !
     Hk = matmul(Hk_tb(:,:,ik),Rhop)
     Hk = matmul(Rhop_hc,Hk)
     Hk = Hk + lgrNC_new
     gz_imt_qpH(ik,:,:) = gz_imt_qpH(ik,:,:) + 2.d0*Hk*itstep*0.5d0     
  end do

end subroutine update_qpH_integral



! function gz_IMTeom_superc(time,y,Nsys) result(f)
!   implicit none
!   !inputs
!   integer                                   :: Nsys ! nr of equations
!   real(8)                                   :: time ! time variable
!   complex(8),dimension(Nsys)                :: y    ! argument array
!   complex(8),dimension(Nsys)                :: f    ! result 
!   !
!   !
!   !
!   real(8),dimension(:),allocatable          :: lgr,delta_out
!   complex(8),dimension(:),allocatable       :: lgr_cmplx
!   complex(8),dimension(Nphi)                :: tmp_gzdot
!   complex(8),dimension(2,Nvdm_AC_opt,Ns,Ns) :: GZa_commAC
!   complex(8),dimension(Nvdm_AC_opt)         :: GZa_commH
!   complex(8),dimension(Ns,Ns)               :: SLa_constr_dot_not
!   complex(8),dimension(2*Ns,2*Ns)           :: SL_constr
!   integer                                   :: iter,Nopt,iphi
!   integer                                   :: i,i0
!   real(8)                                   :: delta
!   !
!   if(allocated(neq_lgr)) deallocate(neq_lgr)
!   allocate(neq_lgr(2,Ns,Ns)); neq_lgr=zero
!   !
!   !HERE write the GZ EQUATIONS OF MOTION
!   if(Nsys.ne.nDynamics) stop "wrong dimensions in the GZ_equations_of_motion"
!   !
!   call get_SLaGZa_not(y,SLa_constr_dot_not,SL_constr,GZa_commH,GZa_commAC)
!   !
!   Nopt=Nvdm_NC_opt
!   allocate(lgr_cmplx(Nopt))
!   allocate(lgr(2*Nopt+1));allocate(delta_out(2*Nopt+1))
!   !
!   !+- compute the derivative such that the derivative of the slater constraint is equal to zero -+!
!   call vdm_AC_stride_m2v(gz_imt_dens_lgr,lgr_cmplx)
!   do i=1,Nvdm_NC_opt
!      lgr(i) = dreal(lgr_cmplx(i))
!      lgr(i+Nvdm_NC_opt) = dimag(lgr_cmplx(i))
!   end do
!   lgr(2*Nopt+1) = gz_imt_unit_lgr
  


  
!   call fsolve(fix_anomalous_lgr_sl,lgr,tol=10.d-12,info=iter)  
!   delta_out = fix_anomalous_lgr_sl(lgr);  
!   write(*,*) "SL lgr fixed"
!   write(*,'(10F8.4)') delta_out,lgr
!   lgr_cmplx=zero
!   do i=1,Nvdm_AC_opt
!      lgr_cmplx(i) = lgr(i)+xi*lgr(i+Nvdm_AC_opt)
!   end do
!   call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(1,:,:))
!   gz_neq_dens_lgrA_slater = neq_lgr(1,:,:)
!   !
!   !
!   call vdm_AC_stride_m2v(gz_neq_dens_lgrA_gzproj,lgr_cmplx)
!   do i=1,Nvdm_AC_opt
!      lgr(i) = dreal(lgr_cmplx(i))
!      lgr(i+Nvdm_AC_opt) = dimag(lgr_cmplx(i))
!   end do
!   call fsolve(fix_anomalous_vdm_gz,lgr,tol=1.d-16,info=iter)
!   delta_out = fix_anomalous_vdm_gz(lgr);  
!   write(*,*) "GZ lgr fixed"
!   write(*,'(10F8.4)')  delta_out,lgr
!   lgr_cmplx=zero
!   do i=1,Nvdm_AC_opt
!      lgr_cmplx(i) = lgr(i)+xi*lgr(i+Nvdm_AC_opt)
!   end do
!   call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(2,:,:))  
!   gz_neq_dens_lgrA_gzproj = neq_lgr(2,:,:)
!   !
!   f = gz_equations_of_motion_superc_lgr_sp(time,y,Nsys)
!   !
! contains
!   !  
!   function fix_anomalous_vdm_gz(lgr) result(delta)
!     implicit none
!     real(8),dimension(:)                :: lgr
!     real(8),dimension(size(lgr))        :: delta
!     complex(8),dimension(:),allocatable :: lgr_cmplx,delta_cmplx
!     complex(8),dimension(Ns,Ns)         :: anomalous_constrGZ_dot,lgrGZ
!     real(8)                             :: tmp_test
!     complex(8),dimension(2,Ns,Ns,Lk)    :: slater_,slater_dot
!     complex(8),dimension(3,Ns,Ns,Lk)    :: slater
!     complex(8),dimension(Nphi)          :: gzproj,gzproj_dot,gztmp
!     complex(8),dimension(Nphi,Nphi)     :: Hproj
!     type(sparse_matrix_csr_z)           :: htmp
!     complex(8)                          :: xtmp
!     real(8)                             :: xtmp_
!     integer                             :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi,i0,i,iopt,jopt
!     complex(8),dimension(Ns,Ns)         :: tmpHk
!     complex(8),dimension(Ns,Ns)         :: Hk
!     complex(8),dimension(2*Ns,2*Ns)     :: tmp
!     complex(8),dimension(2,Ns,Ns)       :: slater_derivatives
!     complex(8),dimension(Ns,Ns)         :: Rhop,Qhop,Rhop_hc,Qhop_hc
!     complex(8),dimension(Ns,Ns)         :: tRR,tRQ,tQR,tQQ
!     complex(8),dimension(Ns,Ns)         :: vdm_natural
!     real(8),dimension(Ns)               :: vdm_diag,n0
!     real(8)                             :: test_slater
!     !+- dump slater_lgr_multipliers
!     allocate(lgr_cmplx(Nvdm_AC_opt),delta_cmplx(Nvdm_AC_opt))
!     do i=1,Nvdm_AC_opt
!        lgr_cmplx(i) = lgr(i) + xi*lgr(i+Nvdm_AC_opt)
!     end do
!     call vdm_AC_stride_v2m(lgr_cmplx,lgrGZ)  
!     !
!     do iopt=1,Nvdm_AC_opt
!        delta_cmplx(iopt) = GZa_commH(iopt)
!        do is=1,Ns
!           do js=1,Ns
!              delta_cmplx(iopt) = delta_cmplx(iopt) + &
!                   lgrGZ(is,js)*GZa_commAC(1,iopt,is,js) + &
!                   conjg(lgrGZ(is,js))*GZa_commAC(2,iopt,is,js) 
!           end do
!        end do
!     end do
!     !
!     do i=1,Nvdm_AC_opt
!        delta(i) = dreal(delta_cmplx(i))
!        delta(i+Nvdm_AC_opt) = dimag(delta_cmplx(i))
!     end do
!     deallocate(delta_cmplx)    
!     !
!     if(GZneq_verbose) then
!        write(*,*) 'GZ constraints'
!        write(*,'(20F18.10)') lgr
!        write(*,'(20F18.10)') delta
!     end if
!     !
!   end function fix_anomalous_vdm_gz
!   !
!   function fix_anomalous_lgr_sl(lgr) result(delta)
!     implicit none
!     real(8),dimension(:) :: lgr
!     real(8),dimension(size(lgr)) :: delta
!     complex(8),dimension(:),allocatable  :: lgr_cmplx,delta_cmplx
!     complex(8),dimension(Ns,Ns) :: anomalous_constrSL_dot,lgrSL
!     real(8) :: tmp_test
!     complex(8),dimension(2,Ns,Ns,Lk) :: slater_,slater_dot
!     complex(8),dimension(3,Ns,Ns,Lk) :: slater
!     complex(8),dimension(Nphi)       :: gzproj,gzproj_dot
!     complex(8),dimension(Nphi,Nphi)  :: Hproj
!     type(sparse_matrix_csr_z)          :: htmp
!     complex(8)                       :: xtmp
!     real(8)                          :: xtmp_
!     integer                          :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi,i0,i
!     complex(8),dimension(Ns,Ns)      :: anomalous_slater_constr_dot
!     complex(8),dimension(2,Ns,Ns)    :: slater_derivatives
!     complex(8),dimension(Ns,Ns)      :: Rhop,Qhop,Rhop_hc,Qhop_hc
!     complex(8),dimension(Ns,Ns)      :: tRR,tRQ,tQR,tQQ
!     complex(8),dimension(Ns,Ns)      :: vdm_natural
!     real(8),dimension(Ns)            :: vdm_diag,n0
!     real(8) :: test_slater
!     !
!     allocate(lgr_cmplx(Nvdm_AC_opt))
!     do i=1,Nvdm_AC_opt
!        lgr_cmplx(i) = lgr(i) + xi*lgr(i+Nvdm_AC_opt)
!     end do
!     call vdm_AC_stride_v2m(lgr_cmplx,lgrSL)
!     !
!     tRQ=lgrSL
!     tQR=conjg(transpose(lgrSL))
!     !
!     anomalous_slater_constr_dot = SLa_constr_dot_not
!     do is=1,Ns
!        do js=1,Ns          
!           !
!           do ks=1,Ns
!              !tQR_{js,ks}SL1(is,ks) - tQR_{ks,js}SL1(is,ks)
!              anomalous_slater_constr_dot(is,js) = anomalous_slater_constr_dot(is,js) + tQR(js,ks)*SL_constr(is,ks)
!              anomalous_slater_constr_dot(is,js) = anomalous_slater_constr_dot(is,js) - tQR(ks,js)*SL_constr(is,ks) !**
!              !tQR_{is,ks}SL3(ks,js) - tQR_{ks,is}SL3(ks,js)
!              anomalous_slater_constr_dot(is,js) = anomalous_slater_constr_dot(is,js) + tQR(is,ks)*SL_constr(ks+Ns,js+Ns) 
!              anomalous_slater_constr_dot(is,js) = anomalous_slater_constr_dot(is,js) - tQR(ks,is)*SL_constr(ks+Ns,js+Ns) !**  
!           end do
!           !
!        end do
!     end do
!     !
!     delta=0.d0
!     allocate(delta_cmplx(Nvdm_AC_opt))
!     call vdm_AC_stride_m2v(anomalous_slater_constr_dot,delta_cmplx)
!     do i=1,Nvdm_AC_opt
!        delta(i) = dreal(delta_cmplx(i))
!        delta(i+Nvdm_AC_opt) = dimag(delta_cmplx(i))
!     end do
!     deallocate(delta_cmplx)
!     !
!     if(GZneq_verbose) then
!        write(*,*) 'SL constraints'
!        write(*,'(20F18.10)') lgr
!        write(*,'(20F18.10)') delta
!     end if
!   end function fix_anomalous_lgr_sl
!   !
!   subroutine get_SLaGZa_not(y,anomalous_constrSL_dot,constrSL,GZa_commH,GZa_commAC)
!     implicit none
!     complex(8),dimension(nDynamics)               :: y
!     complex(8),dimension(Ns,Ns)                   :: anomalous_constrSL_dot,anomalous_constrGZ_dot
!     complex(8),dimension(2*Ns,2*Ns)               :: constrSL
!     complex(8),dimension(Nvdm_AC_opt)             :: GZa_commH
!     complex(8),dimension(2,Nvdm_AC_opt,Ns,Ns) :: GZa_commAC
!     complex(8),dimension(2,Nvdm_AC_opt,Nvdm_AC_opt) :: GZa_commAC_
!     real(8)                                       :: tmp_test
!     complex(8),dimension(2,Ns,Ns,Lk)              :: slater_,slater_dot
!     complex(8),dimension(3,Ns,Ns,Lk)              :: slater
!     complex(8),dimension(Nphi)                    :: gzproj,gzproj_dot
!     complex(8),dimension(Nphi,Nphi)               :: Hproj
!     type(sparse_matrix_csr_z)                     :: htmp
!     complex(8)                                    :: xtmp
!     real(8)                                       :: xtmp_
!     integer                                       :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi,i0,i,iopt,jopt
!     complex(8),dimension(Ns,Ns)                   :: tmpHk
!     complex(8),dimension(Ns,Ns)                   :: Hk
!     complex(8),dimension(2*Ns,2*Ns)               :: tmp,HK_full,dotK
!     complex(8),dimension(2,Ns,Ns)                 :: slater_derivatives
!     complex(8),dimension(Ns,Ns)                   :: Rhop,Qhop,Rhop_hc,Qhop_hc
!     complex(8),dimension(Ns,Ns)                   :: tRR,tRQ,tQR,tQQ
!     complex(8),dimension(Ns,Ns)                   :: vdm_natural
!     real(8),dimension(Ns)                         :: vdm_diag,n0
!     real(8)                                       :: test_slater
!     complex(8),dimension(Nvdm_AC_opt,Nphi)             :: wc,wc_
!     complex(8),dimension(Ns,Ns,Nphi)             :: w_ac,w_ac_
!     !
!     !
!     call dynamicalVector_2_wfMatrix_superc(y,slater_,gzproj)
!     slater(1:2,:,:,:) = slater_(1:2,:,:,:)
!     slater(3,:,:,:) = zero
!     do is=1,Ns
!        do js=1,Ns
!           if(is.eq.js) slater(3,is,js,:) = 1.d0
!           slater(3,is,js,:) = slater(3,is,js,:) - slater_(1,js,is,:)
!        end do
!     end do
!     !
!     do is=1,Ns
!        do js=1,Ns
!           vdm_natural(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
!        end do
!        vdm_diag(is) = dreal(vdm_natural(is,is))
!     end do
!     !
!     Rhop = hopping_renormalization_normal_sp(gzproj,vdm_diag)
!     Qhop = hopping_renormalization_anomalous_sp(gzproj,vdm_diag)
!     n0=vdm_diag
!     !
!     do is=1,Ns
!        do js=1,Ns
!           Rhop_hc(is,js) = conjg(Rhop(js,is))
!           Qhop_hc(is,js) = conjg(Qhop(js,is))
!        end do
!     end do
!     !
!     anomalous_constrSL_dot=zero
!     constrSL=zero
!     !
!     it = im_t2it(time,tstep*0.5d0)
!     !
!     slater_derivatives = zero

!     do ik=1,Lk
!        call get_Hk_t(Hk,ik,time)
!        !
!        tRR = matmul(Hk,Rhop)
!        tRR = matmul(Rhop_hc,tRR)
!        !
!        tRQ = matmul(Hk,Qhop)
!        tRQ = matmul(Rhop_hc,tRQ)
!        !
!        tQR = matmul(Hk,Rhop)
!        tQR = matmul(Qhop_hc,tQR)
!        ! !+- add_lgr_multipliers
!        !
!        tQQ = matmul(Hk,Qhop)
!        tQQ = matmul(Qhop_hc,tQQ)
!        !       
!        constrSL(1:Ns,1:Ns) = constrSL(1:Ns,1:Ns)  + slater(1,:,:,ik)*wtk(ik)
!        constrSL(1:Ns,1+Ns:2*Ns) = constrSL(1:Ns,1+Ns:2*Ns)  + slater(2,:,:,ik)*wtk(ik)
!        constrSL(1+Ns:2*Ns,1:Ns) = constrSL(1+Ns:2*Ns,1:Ns)  + conjg(transpose(slater(2,:,:,ik))*wtk(ik))
!        constrSL(1+Ns:2*Ns,1+Ns:2*Ns) = constrSL(1+Ns:2*Ns,1+Ns:2*Ns)  + slater(3,:,:,ik)*wtk(ik)
!        !
!        do is=1,Ns
!           do js=1,Ns                   
!              slater_dot(2,is,js,ik) = zero
!              do ks=1,Ns
!                 !tRR_{ks,js}SL2(ks,is) - tRR_{ks,is}SL2(ks,js)
!                 slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tRR(ks,js)*slater(2,ks,is,ik)
!                 slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tRR(ks,is)*slater(2,ks,js,ik) !**                    
!                 !tQR_{js,ks}SL1(is,ks) - tQR_{ks,js}SL1(is,ks)
!                 slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQR(js,ks)*slater(1,is,ks,ik)
!                 slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQR(ks,js)*slater(1,is,ks,ik) !**
!                 !tQR_{is,ks}SL3(ks,js) - tQR_{ks,is}SL3(ks,js)
!                 slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQR(is,ks)*slater(3,ks,js,ik) 
!                 slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQR(ks,is)*slater(3,ks,js,ik) !**
!                 !tQQ_{is,ks}SL2(ks,js) - tQQ_{js,ks}SL2(ks,is)
!                 slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + tQQ(is,ks)*slater(2,ks,js,ik)
!                 slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - tQQ(js,ks)*slater(2,ks,is,ik) !**
!                 !
!              end do
!              anomalous_constrSL_dot(is,js) = anomalous_constrSL_dot(is,js) + slater_dot(2,is,js,ik)*wtk(ik)
!           end do
!        end do

!        do is=1,Ns
!           do js=1,Ns
!              do ks=1,Ns
!                 do kks=1,Ns
!                    !
!                    slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
!                         conjg(Rhop(kks,ks))*Hk(kks,is)*slater(1,ks,js,ik)*wtk(ik)    
!                    slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
!                         conjg(Qhop(kks,ks))*Hk(kks,is)*conjg(slater(2,js,ks,ik))*wtk(ik)
!                    slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
!                         conjg(Rhop(kks,ks))*Hk(kks,is)*slater(2,ks,js,ik)*wtk(ik) 
!                    slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
!                         conjg(Qhop(kks,ks))*Hk(kks,is)*slater(3,ks,js,ik)*wtk(ik)
!                 end do
!              end do
!           end do
!        end do
!     end do
!     !
!     Uloc=Uloc_t(:,it)
!     Ust =Ust_t(it)
!     Jh=Jh_t(it)
!     Jsf=Jsf_t(it)
!     Jph=Jph_t(it)
!     eLevels = eLevels_t(:,it)
!     !
!     gzproj_dot = get_local_hamiltonian_HLOCphi(gzproj,eLevels)
!     do is=1,Ns
!        do js=1,Ns
!           !
!           xtmp=slater_derivatives(1,is,js)/sqrt(n0(js)*(1.d0-n0(js)))
!           xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
!           if(xtmp/=zero) then
!              htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Rhop(is,js),xtmp)
!              gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
!              call sp_delete_matrix(htmp)
!              xtmp=conjg(slater_derivatives(1,is,js))/sqrt(n0(js)*(1.d0-n0(js)))
!              htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Rhop_hc(is,js),xtmp)
!              gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
!              call sp_delete_matrix(htmp)
!           end if
!           xtmp=slater_derivatives(1,is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
!           xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
!           if(xtmp/=zero) then
!              htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(js,js),xtmp)
!              gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
!              call sp_delete_matrix(htmp)        
!              xtmp=conjg(slater_derivatives(1,is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
!              htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(js,js),xtmp)
!              gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
!              call sp_delete_matrix(htmp)
!           end if
!           xtmp=slater_derivatives(2,is,js)/sqrt(n0(js)*(1.d0-n0(js)))
!           xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
!           if(xtmp/=zero) then
!              htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Qhop(is,js),xtmp)
!              gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
!              call sp_delete_matrix(htmp)
!              xtmp=conjg(slater_derivatives(2,is,js))/sqrt(n0(js)*(1.d0-n0(js)))
!              htmp=sp_scalar_matrix_csr(phi_spTraces_basis_Qhop_hc(is,js),xtmp)
!              gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
!              call sp_delete_matrix(htmp)
!           end if
!           xtmp=slater_derivatives(2,is,js)*Qhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
!           xtmp_=dreal(xtmp)**2.d0+dimag(xtmp)**2.d0
!           if(xtmp/=zero) then
!              htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens(js,js),xtmp)
!              gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
!              call sp_delete_matrix(htmp)
!              xtmp=conjg(slater_derivatives(2,is,js)*Qhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))
!              htmp=sp_scalar_matrix_csr(phi_spTraces_basis_dens_hc(js,js),xtmp)
!              gzproj_dot = gzproj_dot + sp_matrix_vector_product_csr_z(Nphi,htmp,gzproj)
!              call sp_delete_matrix(htmp)
!           end if
!        end do
!     end do
!     !
!     do is=1,Ns
!        do js=1,Ns
!           w_ac(is,js,:) = sp_matrix_vector_product_csr_z(Nphi,phi_spTraces_basis_dens_anomalous(is,js),gzproj)
!           w_ac_(is,js,:) = sp_matrix_vector_product_csr_z(Nphi,phi_spTraces_basis_dens_anomalous_hc(is,js),gzproj)
!        end do
!     end do
!     !
!     do iopt=1,Nvdm_AC_opt
!        iis=IS_vdmAC(iopt)
!        jjs=JS_vdmAC(iopt)
!        GZa_commH(iopt) = zero       
!        do is=1,Ns
!           do js=1,Ns
!              GZa_commAC(:,iopt,is,js) = zero
!              do iphi=1,Nphi
!                 GZa_commAC(1,iopt,is,js) = GZa_commAC(1,iopt,is,js) + &
!                      conjg(w_ac_(iis,jjs,iphi))*w_ac(is,js,iphi) - conjg(w_ac_(is,js,iphi))*w_ac(iis,jjs,iphi)
!                 GZa_commAC(2,iopt,is,js) = GZa_commAC(2,iopt,is,js) + &
!                      conjg(w_ac_(iis,jjs,iphi))*w_ac_(is,js,iphi) - conjg(w_ac(is,js,iphi))*w_ac(iis,jjs,iphi)
!              end do
!           end do
!        end do
!        do iphi=1,Nphi
!           GZa_commH(iopt) = GZa_commH(iopt) +&
!                conjg(w_ac_(iis,jjs,iphi))*gzproj_dot(iphi) - conjg(gzproj_dot(iphi))*w_ac(iis,jjs,iphi)
!        end do
!     end do
!   end subroutine get_SLaGZa_not
!   !
! end function gz_IMTeom_superc
