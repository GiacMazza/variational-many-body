function R_VDM_free_zeros_(x)  result(Fout)   !+- change this name
  real(8),dimension(:)                :: x      
  !
  real(8),dimension(size(x))          :: Fout
  real(8),dimension(size(x))          :: tmp_x
  real(8)                             :: E_Hstar,E_Hloc
  complex(8),dimension(Nphi)          :: GZvect
  !
  complex(8),dimension(Ns,Ns)         :: Rhop,Rhop_GZproj
  complex(8),dimension(Ns,Ns)       :: slater_derivatives
  complex(8),dimension(Ns,Ns)       :: slater_lgr_multip
  complex(8),dimension(Ns,Ns)       :: proj_lgr_multip
  complex(8),dimension(Ns,Ns)         :: Rhop_lgr_multip
  !
  complex(8),dimension(Ns,Ns)       :: n0_slater,n0_GZproj
  real(8),dimension(Ns)               :: n0_diag
  real(8)                             :: tmp,delta_out
  complex(8)                          :: tmpR
  integer                             :: imap,is,js,i,ix,i0
  integer                             :: Nslater_lgr,NQhop,NRhop,Nproj_lgr,Nopt
  !
  complex(8),dimension(:),allocatable :: dump_input
  complex(8),dimension(:),allocatable :: dump_output
  !
  complex(8),dimension(Ns,Ns)         :: Rhop_out
  complex(8),dimension(Ns,Ns)       :: slater_out,GZproj_out
  logical                           :: test_n
  !
  NRhop = 2*NRhop_opt
  !
  Nslater_lgr = 2*Nvdm_NC_opt
  Nproj_lgr = 2*Nvdm_NCoff_opt
  !
  Nopt = NRhop + Nslater_lgr + Nproj_lgr
  if(size(x).ne.Nopt) stop "R_VDM_free_zeros_: size input vector not equal to Nopt."
  !
  call dump2mats(x,Rhop,slater_lgr_multip,proj_lgr_multip)
  !  
  call slater_minimization_fixed_lgr(Rhop,slater_lgr_multip,E_Hstar,n0_slater,slater_derivatives)  
  !
  test_n=.true.
  !  
  do is=1,Ns
     n0_diag(is) = n0_slater(is,is)
     if(n0_diag(is).gt.1.d0-1.d-8.or.n0_diag(is).lt.1.d-8) test_n=.false.
  end do
  !

  if(test_n) then

     do is=1,Ns
        do js=1,Ns
           Rhop_lgr_multip(is,js) = -1.d0*slater_derivatives(is,js)/sqrt(n0_slater(js,js)*(1.d0-n0_slater(js,js)))
        end do
     end do
     !
     do is=1,Ns
        tmp = 0.d0
        do js=1,Ns
           tmp = tmp + dreal(Rhop_lgr_multip(is,js)*Rhop(is,js))
        end do
        proj_lgr_multip(is,is) = -0.5d0*slater_lgr_multip(is,is) + &  
             0.5d0*(1.d0-2.d0*n0_slater(is,is))/sqrt(n0_slater(is,is)*(1.d0-n0_slater(is,is)))*tmp 
     end do
     !
     ! write(*,*) 'prima',n0_diag
     ! write(*,*) slater_lgr_multip
     ! write(*,*) slater_lgr_multip
     call gz_proj_minimization_fixed_lgr_hop(n0_diag,proj_lgr_multip,Rhop_lgr_multip,E_Hloc,GZvect)   
     !write(*,*) 'dopo' 
     !
     n0_GZproj = 0.d0
     do is=1,Ns
        do js=1,Ns
           n0_GZproj(is,js)   = trace_phi_basis(GZvect,phi_traces_basis_dens(is,js,:,:))
           Rhop_GZProj(is,js) = trace_phi_basis(GZvect,phi_traces_basis_Rhop(is,js,:,:))
        end do
     end do
     !  
     do is=1,Ns
        do js=1,Ns
           Rhop_out(is,js) = Rhop_GZproj(is,js)-Rhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js)))
           if(is.eq.js) then
              slater_out(is,js) = n0_slater(is,js) - n0_GZproj(is,js)
           else
              slater_out(is,js) = n0_slater(is,js)
           end if
           GZproj_out(is,js) = n0_GZproj(is,js)
        end do
     end do
     !
     call dump2vec(tmp_x,Rhop_out,slater_out,GZproj_out)
     !
     Fout=tmp_x
     delta_out=0.d0
     do is=1,Nopt
        delta_out = delta_out + Fout(is)*Fout(is)
     end do

  else

     Fout=1
  end if

  !
  if(GZmin_verbose) then
     write(root_unit,'(30F18.10)') delta_out,Fout(1:Nopt)
     write(xmin_unit,'(30F18.10)') x(1:Nopt)
  end if
  !
  if(optimization_flag) then
     !+- store final informations to global variables -+!                   
     GZ_vector      = GZvect
     GZ_opt_energy  = E_Hstar+E_Hloc
     GZ_opt_kinetic = E_Hstar
     GZ_opt_Eloc    = E_Hloc
     !+- here I should modify this ?? porca troia quante robbe da modificare...
     GZ_opt_slater_lgr = slater_lgr_multip
     GZ_opt_proj_lgr = proj_lgr_multip
     !
  end if
  !
end function R_VDM_free_zeros_
!












function R_VDM_free_zeros__(x)  result(Fout)   !+- change this name
  real(8),dimension(:)                :: x      
  real(8),dimension(size(x))          :: Fout
  !
  integer                             :: Nslater_lgr,NQhop,NRhop,Nproj_lgr,Nopt
  real(8),dimension(:),allocatable    :: x_orig,Fout_orig
  !
  NRhop = 2*NRhop_opt
  !
  Nslater_lgr = 2*Nvdm_NC_opt
  Nproj_lgr = 2*Nvdm_NCoff_opt
  !
  Nopt = NRhop + Nslater_lgr + Nproj_lgr
  if(size(x).ne.Nopt_reduced) stop "R_VDM_free_zeros__: size input vector not equal to Nopt_reduced."
  !
  allocate(x_orig(Nopt),Fout_orig(Nopt))
  call stride_zeros_red2orig(x,x_orig)
  Fout_orig = R_VDM_free_zeros_(x_orig)
  call stride_zeros_orig2red(Fout_orig,Fout)
  !
end function R_VDM_free_zeros__
!







!+- SUPERCONDUCTING ROUTINES -+!


!
!
function R_Q_VDM_free_zeros_superc(x)  result(Fout)
  real(8),dimension(:)                :: x      
  !
  real(8),dimension(size(x))          :: Fout
  real(8),dimension(size(x))          :: tmp_x
  real(8)                             :: E_Hstar,E_Hloc
  complex(8),dimension(Nphi)          :: GZvect
  !
  complex(8),dimension(Ns,Ns)         :: Rhop,Rhop_GZproj
  complex(8),dimension(Ns,Ns)         :: Qhop,Qhop_GZproj
  complex(8),dimension(2,Ns,Ns)       :: slater_derivatives
  complex(8),dimension(2,Ns,Ns)       :: slater_lgr_multip
  complex(8),dimension(2,Ns,Ns)       :: proj_lgr_multip
  complex(8),dimension(Ns,Ns)         :: Rhop_lgr_multip
  complex(8),dimension(Ns,Ns)         :: Qhop_lgr_multip
  !
  complex(8),dimension(2,Ns,Ns)       :: n0_slater,n0_GZproj
  real(8),dimension(Ns)               :: n0_diag
  real(8)                             :: tmp,delta_out
  complex(8)                          :: tmpR
  integer                             :: imap,is,js,i,ix,i0
  integer                             :: Nslater_lgr,NQhop,NRhop,Nproj_lgr,Nopt
  !
  complex(8),dimension(:),allocatable :: dump_input
  complex(8),dimension(:),allocatable :: dump_output
  !
  complex(8),dimension(Ns,Ns)         :: Rhop_out,Qhop_out
  complex(8),dimension(2,Ns,Ns)       :: slater_out,GZproj_out
  !
  NRhop = 2*NRhop_opt
  NQhop = 2*NQhop_opt
  !
  Nslater_lgr = 2*(Nvdm_NC_opt+Nvdm_AC_opt)
  Nproj_lgr = 2*(Nvdm_NCoff_opt+Nvdm_AC_opt) 
  !
  Nopt = NRhop + NQhop + Nslater_lgr + Nproj_lgr
  if(size(x).ne.Nopt) stop "R_Q_VDM_free_seros_superc: size input vector not equal to Nopt."
  !
  call dump2mats_superc(x,Rhop,Qhop,slater_lgr_multip,proj_lgr_multip)
  !
  !<TMP TEST
  ! slater_lgr_multip(1,:,:)=dreal(slater_lgr_multip(1,:,:))
  ! slater_lgr_multip(2,:,:)=zero
  ! proj_lgr_multip(:,:,:)=zero
  !TMP TEST>
  !
  call slater_minimization_fixed_lgr_superc(Rhop,Qhop,slater_lgr_multip,E_Hstar,n0_slater,slater_derivatives)
  !
  do is=1,Ns
     do js=1,Ns
        Rhop_lgr_multip(is,js) = -1.d0*slater_derivatives(1,is,js)/sqrt(n0_slater(1,js,js)*(1.d0-n0_slater(1,js,js)))
        Qhop_lgr_multip(is,js) = -1.d0*slater_derivatives(2,is,js)/sqrt(n0_slater(1,js,js)*(1.d0-n0_slater(1,js,js)))
     end do
  end do
  !
  do is=1,Ns
     tmp = 0.d0
     do js=1,Ns
        tmp = tmp + dreal(Rhop_lgr_multip(is,js)*Rhop(is,js)) + dreal(Qhop_lgr_multip(is,js)*Qhop(is,js))
     end do
     proj_lgr_multip(1,is,is) = -0.5*slater_lgr_multip(1,is,is) + & 
          0.5d0*(1.d0-2.d0*n0_slater(1,is,is))/sqrt(n0_slater(1,is,is)*(1.d0-n0_slater(1,is,is)))*tmp 
  end do
  !
  do is=1,Ns
     n0_diag(is) = n0_slater(1,is,is)
  end do
  !
  call gz_proj_minimization_fixed_lgr_hop_superc(n0_diag,proj_lgr_multip,Rhop_lgr_multip,Qhop_lgr_multip,E_Hloc,GZvect)
  !
  n0_GZproj = 0.d0
  do is=1,Ns
     do js=1,Ns
        n0_GZproj(1,is,js) = trace_phi_basis(GZvect,phi_traces_basis_dens(is,js,:,:))
        n0_GZproj(2,is,js) = trace_phi_basis(GZvect,phi_traces_basis_dens_anomalous(is,js,:,:))
        Rhop_GZProj(is,js) = trace_phi_basis(GZvect,phi_traces_basis_Rhop(is,js,:,:))
        Qhop_GZProj(is,js) = trace_phi_basis(GZvect,phi_traces_basis_Qhop(is,js,:,:))
     end do
  end do
  !  
  do is=1,Ns
     do js=1,Ns
        Rhop_out(is,js) = Rhop_GZproj(is,js)-Rhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js)))
        Qhop_out(is,js) = Qhop_GZproj(is,js)-Qhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js)))
        if(is.eq.js) then
           slater_out(1,is,js) = n0_slater(1,is,js) - n0_GZproj(1,is,js)
        else
           slater_out(1,is,js) = n0_slater(1,is,js)
        end if
        slater_out(2,is,js) = n0_slater(2,is,js)
        GZproj_out(1,is,js) = n0_GZproj(1,is,js)
        GZproj_out(2,is,js) = n0_GZproj(2,is,js)
     end do
  end do
  !
  call dump2vec_superc(tmp_x,Rhop_out,Qhop_out,slater_out,GZproj_out)
  Fout=tmp_x
  !
  delta_out=0.d0
  do is=1,Nopt
     delta_out = delta_out + Fout(is)*Fout(is)
  end do
  !
  if(GZmin_verbose) then
     write(root_unit,'(30F18.10)') delta_out,Fout(1:Nopt)
     write(xmin_unit,'(30F18.10)') x(1:Nopt)
  end if
  !
  if(optimization_flag) then
     !+- store final informations to global variables -+!                   
     GZ_vector      = GZvect
     GZ_opt_energy  = E_Hstar+E_Hloc
     GZ_opt_kinetic = E_Hstar
     GZ_opt_Eloc    = E_Hloc
     !+- here I should modify this
     GZ_opt_slater_lgr_superc = slater_lgr_multip
     GZ_opt_proj_lgr_superc = proj_lgr_multip
  end if
  !
end function R_Q_VDM_free_zeros_superc











function R_Q_VDM_free_zeros_superc_(x)  result(Fout)
  real(8),dimension(:)                :: x      
  real(8),dimension(size(x))          :: Fout
  !
  integer                             :: Nslater_lgr,NQhop,NRhop,Nproj_lgr,Nopt
  !
  real(8),dimension(:),allocatable                :: x_orig,Fout_orig      
  NRhop = 2*NRhop_opt
  NQhop = 2*NQhop_opt
  !
  Nslater_lgr = 2*(Nvdm_NC_opt+Nvdm_AC_opt)
  Nproj_lgr = 2*(Nvdm_NCoff_opt+Nvdm_AC_opt) 
  !
  Nopt = NRhop + NQhop + Nslater_lgr + Nproj_lgr
  !
  if(size(x).ne.Nopt_reduced) stop "R_Q_VDM_free_zeros_superc_: size input vector not equal to Nopt_reduced."
  allocate(x_orig(Nopt),Fout_orig(Nopt))
  call stride_zeros_red2orig(x,x_orig)
  Fout_orig = R_Q_VDM_free_zeros_superc(x_orig)
  call stride_zeros_orig2red(Fout_orig,Fout)
  !
end function R_Q_VDM_free_zeros_superc_






!
subroutine dump2vec(x,Rhop,slater_lgr,GZproj_lgr)
  real(8),dimension(:)                :: x
  complex(8),dimension(Ns,Ns)         :: Rhop
  complex(8),dimension(Ns,Ns)       :: slater_lgr
  complex(8),dimension(Ns,Ns)       :: GZproj_lgr
  integer                             :: NRhop,NQhop,Nslater_lgr,Nproj_lgr,Nopt  
  integer                             :: i,ix,i0,imap,is,js
  complex(8),dimension(:),allocatable :: dump_output
  !
  Nopt = 2*NRhop_opt + 2*Nvdm_NC_opt + 2*Nvdm_NCoff_opt 
  if(size(x).ne.Nopt) stop "wrong dimensions @ dump2vec"
  !
  i=0
  i0=i
  if(Nrhop_opt.gt.0) then
     allocate(dump_output(NRhop_opt))
     call Rhop_stride_m2v(Rhop,dump_output)
     do ix=1,NRhop_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+NRhop_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)
  end if
  !
  i0=i
  if(Nvdm_NC_opt.gt.0) then
     allocate(dump_output(Nvdm_NC_opt))
     call vdm_NC_stride_m2v(slater_lgr,dump_output)
     do ix=1,Nvdm_NC_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+Nvdm_NC_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)
  end if
  !
  i0=i
  if(Nvdm_NCoff_opt.gt.0) then
     allocate(dump_output(Nvdm_NCoff_opt))
     call vdm_NCoff_stride_m2v(GZproj_lgr,dump_output)
     do ix=1,Nvdm_NCoff_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+Nvdm_NCoff_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)
  end if
  !
end subroutine dump2vec
!
subroutine dump2mats(x,Rhop,slater_lgr_multip,proj_lgr_multip)
  real(8),dimension(:)                :: x
  complex(8),dimension(Ns,Ns)         :: Rhop
  complex(8),dimension(Ns,Ns)       :: slater_lgr_multip
  complex(8),dimension(Ns,Ns)       :: proj_lgr_multip
  integer                             :: NRhop,NQhop,Nslater_lgr,Nproj_lgr,Nopt  
  integer                             :: i,ix,i0,imap,is,js
  complex(8),dimension(:),allocatable :: dump_input
  Nopt = 2*NRhop_opt + 2*Nvdm_NC_opt + 2*Nvdm_NCoff_opt 
  if(size(x).ne.Nopt) stop "wrong dimensions @ dump2mats_superc"
  !
  Rhop = zero
  slater_lgr_multip = zero
  proj_lgr_multip = zero
  !
  i=0
  i0=i
  if(NRhop_opt.gt.0) then
     allocate(dump_input(NRhop_opt))
     do ix=1,NRhop_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+NRhop_opt); i=i+2
     end do
     call Rhop_stride_v2m(dump_input,Rhop)
     deallocate(dump_input)
  end if
  !
  i0=i
  if(Nvdm_NC_opt.gt.0) then
     allocate(dump_input(Nvdm_NC_opt))
     do ix=1,Nvdm_NC_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+Nvdm_NC_opt); i=i+2
     end do
     call vdm_NC_stride_v2m(dump_input,slater_lgr_multip)
     deallocate(dump_input)
  end if
  !
  i0=i
  if(Nvdm_NCoff_opt.gt.0) then
     allocate(dump_input(Nvdm_NCoff_opt))
     do ix=1,Nvdm_NCoff_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+Nvdm_NCoff_opt); i=i+2
     end do
     call vdm_NCoff_stride_v2m(dump_input,proj_lgr_multip)
     deallocate(dump_input)
  end if
  !
end subroutine dump2mats










subroutine dump2mats_superc(x,Rhop,Qhop,slater_lgr_multip,proj_lgr_multip)
  real(8),dimension(:)                :: x
  complex(8),dimension(Ns,Ns)         :: Rhop
  complex(8),dimension(Ns,Ns)         :: Qhop
  complex(8),dimension(2,Ns,Ns)       :: slater_lgr_multip
  complex(8),dimension(2,Ns,Ns)       :: proj_lgr_multip
  integer                             :: NRhop,NQhop,Nslater_lgr,Nproj_lgr,Nopt  
  integer                             :: i,ix,i0,imap,is,js
  complex(8),dimension(:),allocatable :: dump_input
  Nopt = 2*NRhop_opt + 2*NQhop_opt + 2*Nvdm_NC_opt + 2*Nvdm_AC_opt + 2*Nvdm_NCoff_opt + 2*Nvdm_AC_opt
  if(size(x).ne.Nopt) stop "wrong dimensions @ dump2mats_superc"
  !
  Rhop = zero
  Qhop = zero
  slater_lgr_multip = zero
  proj_lgr_multip = zero
  !
  i=0
  i0=i
  if(NRhop_opt.gt.0) then
     allocate(dump_input(NRhop_opt))
     do ix=1,NRhop_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+NRhop_opt); i=i+2
     end do
     call Rhop_stride_v2m(dump_input,Rhop)
     deallocate(dump_input)
  end if
  !
  i0=i
  if(NQhop_opt.gt.0) then
     allocate(dump_input(NQhop_opt))
     do ix=1,NQhop_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+NQhop_opt); i=i+2
     end do
     call Qhop_stride_v2m(dump_input,Qhop)
     deallocate(dump_input)
  end if
  !
  i0=i
  if(Nvdm_NC_opt.gt.0) then
     allocate(dump_input(Nvdm_NC_opt))
     do ix=1,Nvdm_NC_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+Nvdm_NC_opt); i=i+2
     end do
     call vdm_NC_stride_v2m(dump_input,slater_lgr_multip(1,:,:))
     deallocate(dump_input)
  end if
  i0=i
  if(Nvdm_AC_opt.gt.0) then
     allocate(dump_input(Nvdm_AC_opt))
     do ix=1,Nvdm_AC_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+Nvdm_AC_opt); i=i+2
     end do
     call vdm_AC_stride_v2m(dump_input,slater_lgr_multip(2,:,:))
     deallocate(dump_input)
  end if
  i0=i
  if(Nvdm_NCoff_opt.gt.0) then
     allocate(dump_input(Nvdm_NCoff_opt))
     do ix=1,Nvdm_NCoff_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+Nvdm_NCoff_opt); i=i+2
     end do
     call vdm_NCoff_stride_v2m(dump_input,proj_lgr_multip(1,:,:))
     deallocate(dump_input)
  end if
  i0=i
  if(Nvdm_AC_opt.gt.0) then
     allocate(dump_input(Nvdm_AC_opt))
     do ix=1,Nvdm_AC_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+Nvdm_AC_opt); i=i+2
     end do
     call vdm_AC_stride_v2m(dump_input,proj_lgr_multip(2,:,:))
     deallocate(dump_input)
  end if
  !
end subroutine dump2mats_superc
!
subroutine dump2vec_superc(x,Rhop,Qhop,slater_lgr,GZproj_lgr)
  real(8),dimension(:)                :: x
  complex(8),dimension(Ns,Ns)         :: Rhop
  complex(8),dimension(Ns,Ns)         :: Qhop
  complex(8),dimension(2,Ns,Ns)       :: slater_lgr
  complex(8),dimension(2,Ns,Ns)       :: GZproj_lgr
  integer                             :: NRhop,NQhop,Nslater_lgr,Nproj_lgr,Nopt  
  integer                             :: i,ix,i0,imap,is,js
  complex(8),dimension(:),allocatable :: dump_output
  !
  Nopt = 2*NRhop_opt + 2*NQhop_opt + 2*Nvdm_NC_opt + 2*Nvdm_AC_opt + 2*Nvdm_NCoff_opt + 2*Nvdm_AC_opt
  if(size(x).ne.Nopt) stop "wrong dimensions @ dump2vec_superc"
  !
  i=0
  i0=i
  if(Nrhop_opt.gt.0) then
     allocate(dump_output(NRhop_opt))
     call Rhop_stride_m2v(Rhop,dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map(is,js)
     !       if(imap.gt.0) dump_output(imap) = Rhop(is,js)
     !    end do
     ! end do
     do ix=1,NRhop_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+NRhop_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)
  end if
  !
  i0=i
  if(NQhop_opt.gt.0) then
     allocate(dump_output(NQhop_opt))
     call Qhop_stride_m2v(Qhop,dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map_anomalous(is,js)
     !       if(imap.gt.0) dump_output(imap) = Qhop(is,js)
     !    end do
     ! end do
     do ix=1,NQhop_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+NQhop_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)     
  end if
  !
  i0=i
  if(Nvdm_NC_opt.gt.0) then
     allocate(dump_output(Nvdm_NC_opt))
     call vdm_NC_stride_m2v(slater_lgr(1,:,:),dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map(is,js)
     !       if(imap.gt.0) dump_output(imap) = slater_lgr(1,is,js)
     !    end do
     ! end do
     do ix=1,Nvdm_NC_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+Nvdm_NC_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)
  end if
  i0=i
  if(Nvdm_AC_opt.gt.0) then
     allocate(dump_output(Nvdm_AC_opt))
     call vdm_AC_stride_m2v(slater_lgr(2,:,:),dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map_anomalous(is,js)
     !       if(imap.gt.0) dump_output(imap) = slater_lgr(2,is,js)
     !    end do
     ! end do
     do ix=1,Nvdm_AC_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+Nvdm_AC_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)     
  end if
  i0=i
  if(Nvdm_NCoff_opt.gt.0) then
     allocate(dump_output(Nvdm_NCoff_opt))
     call vdm_NCoff_stride_m2v(GZproj_lgr(1,:,:),dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map(is,js)
     !       if(imap.gt.0.and.is.ne.js) dump_output(imap) = GZproj_lgr(1,is,js)
     !    end do
     ! end do
     do ix=1,Nvdm_NCoff_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+Nvdm_NCoff_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)
  end if
  i0=i
  if(Nvdm_AC_opt.gt.0) then
     allocate(dump_output(Nvdm_AC_opt))
     call vdm_AC_stride_m2v(GZproj_lgr(2,:,:),dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map_anomalous(is,js)
     !       if(imap.gt.0) dump_output(imap) = GZproj_lgr(2,is,js)
     !    end do
     ! end do
     do ix=1,Nvdm_AC_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+Nvdm_AC_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)     
  end if
  !
end subroutine dump2vec_superc

















