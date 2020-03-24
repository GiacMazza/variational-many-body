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
  USE GZ_TWO_SITES_FOCK
  !
  implicit none
  real(8),dimension(:),allocatable :: epsik,hybik
  real(8) :: t
  integer :: Nx,out_unit,is,js,ik,it,itt,i,iorb,ispin,isite
  integer :: jorb,jspin,jsite
  integer :: ilm,jlm,iilm,jjlm,iw
  integer :: uio
  integer :: ix,jx
  integer :: i2,j2,k2
  integer :: nprint
  !
  character(len=200) :: store_dir,read_dir,read_optWF_dir,read_finSC_dir
  complex(8),dimension(:,:,:,:),allocatable :: slater_init
  complex(8),dimension(:),allocatable     :: gz_proj_init
  !
  complex(8),dimension(:,:),allocatable :: bcs_wf
  !
  complex(8),dimension(:),allocatable     :: psi_t,psi_bcs_t
  real(8),dimension(:,:),allocatable      :: Ut 
  real(8),dimension(:),allocatable      :: Jht
  real(8) :: r,s,tmpU,Ubcs,Ubcs0,Ubcsf
  !
  integer :: unit_neq_hloc
  integer :: unit_neq_local_dens
  integer :: unit_neq_local_dens_dens
  integer :: unit_neq_ene
  integer :: unit_neq_dens_constrSL
  integer :: unit_neq_dens_constrGZ
  integer :: unit_neq_dens_constrSLa
  integer :: unit_neq_dens_constrGZa
  integer :: unit_neq_constrU
  integer :: unit_neq_Rhop
  integer :: unit_neq_Qhop
  integer :: unit_neq_AngMom
  integer :: unit_neq_sc_order
  integer :: unit_neq_nqp
  integer :: unit_neq_bcs
  integer :: unit_proj
  !
  !+- observables -+!
  complex(8),dimension(:),allocatable   :: Rhop
  complex(8),dimension(:,:),allocatable   :: Rhop_matrix,Qhop_matrix
  complex(8),dimension(:,:),allocatable :: local_density_matrix
  real(8),dimension(:,:),allocatable    :: local_dens_dens
  complex(8),dimension(:,:,:),allocatable :: dens_constrSL
  complex(8),dimension(:,:,:),allocatable :: dens_constrGZ
  real(8)                               :: unitary_constr
  real(8),dimension(4)                  :: local_angular_momenta
  real(8),dimension(3)                  :: energies
  complex(8),dimension(:,:),allocatable             :: sc_order
  complex(8),dimension(:),allocatable :: neq_gzproj
  real(8),dimension(:,:),allocatable :: nqp 
  !
  real(8),dimension(:),allocatable      :: dump_vect
  real(8) :: fin_sc_dir
  real(8) :: Uneq,Uneq0,tStart_neqU,tRamp_neqU,tSin_neqU,dUneq
  real(8) :: Jhneq,Jhneq0,tStart_neqJ,tRamp_neqJ,tSin_neqJ,dJneq
  complex(8) :: bcs_sc_order,bcs_delta
  real(8) :: bcs_Kenergy,bcs_Uenergy,phiBCS,bcs_dens
  real(8) :: sc_phase
  logical :: bcs_neq
  logical :: linear_ramp,trpz
  !

  real(8),dimension(:,:,:,:),allocatable :: cc_ij,ccH_ij
  complex(8),dimension(:,:,:,:),allocatable :: ccHlm_ij
  complex(8) :: cc_tmp

  complex(8),dimension(:,:),allocatable :: aph_ij
  
  real(8),dimension(:,:,:,:),allocatable :: nn_ij
  complex(8),dimension(:,:),allocatable :: Hdw,Hlm,Hlm_cut

  complex(8),dimension(:,:),allocatable :: pnn
  complex(8),dimension(:,:),allocatable :: pnn_h2 
  
  real(8),dimension(:),allocatable :: Edw,eigv_lm,eigv_cut

  integer :: int_cut
  real(8) :: Uw,Vw,Utmp
  real(8),dimension(:,:),allocatable :: Uorb,Vorb

  real(8) :: dwell,Lwell
  real(8) :: xlmin,xlmax,xrmin,xrmax,xx,yy
  real(8),dimension(:),allocatable :: xr,wtr

  integer :: Lreal
  real(8) :: wmin,wmax
  
  integer :: Nph,Nh,Nh_cut

  real(8) :: wph
  real(8) :: Daa,Dpa
  
  integer :: Nxw,Nxd  
  real(8),dimension(:,:),allocatable :: ewells

  real(8),dimension(:,:,:),allocatable :: psi_w,dpsi_w
  real(8) :: kn
  
  !
  real(8) :: hbarc,mel
  real(8) :: chiLR,chi_tmp,chi_ij,zeta,chiLR_bare
  complex(8),dimension(:),allocatable :: gph
  complex(8) :: gtmp

  real(8) :: threshold

  real(8) :: leff !+- effective length of the cavity. i.e. (e A0/hbar) \sim 1/leff

  logical :: lm_coupling
  
  hbarc=Planck_constant_in_eV_s/2.d0/pi*speed_of_light_in_vacuum*1.d9 !+- speed of light in nm/s !!
  mel=electron_mass_energy_equivalent_in_Mev*1.d6
  
  call parse_input_variable(Uw,"Uw","inputIX.conf",default=1.d0)  
  call parse_input_variable(Vw,"Vw","inputIX.conf",default=1.d0)
  call parse_input_variable(int_cut,"int_cut","inputIX.conf",default=2)

  call parse_input_variable(dwell,"dwell","inputIX.conf",default=20.d0)
  call parse_input_variable(Lwell,"Lwell","inputIX.conf",default=20.d0)
  call parse_input_variable(leff,"leff","inputIX.conf",default=100.d0)

  call parse_input_variable(Nxw,"Nxw","inputIX.conf",default=1000)  
  call parse_input_variable(Nxd,"Nxd","inputIX.conf",default=100)

  call parse_input_variable(Lreal,"Lreal","inputIX.conf",default=500)

  call parse_input_variable(wmin,"WMIN","inputIX.conf",default=0.d0)
  call parse_input_variable(wmax,"WMAX","inputIX.conf",default=8.d0)

  call parse_input_variable(threshold,"threshold","inputIX.conf",default=1.d-9)
  !
  call parse_input_variable(Nph,"Nph","inputIX.conf",default=10)
  call parse_input_variable(wph,"wph","inputIX.conf",default=1.d0) !+- this is defined the ratio with the first inter-subband transition


  call parse_input_variable(lm_coupling,"lm_coupling","inputIX.conf",default=.true.)

  !
  call read_input("inputIX.conf")
  call save_input_file("inputIX.conf")

  write(*,*) 'hbarc^2/m [eV nm^2]',hbarc**2.0/mel

  if(Norb.gt.6) stop "Norb must be <= 6" 
  
  !+- initialize the fock space -+!
  call initialize_two_sites_fock_space
  !
  
  !+- build the single particle operators on the two-particle subspace -+!
  allocate(cc_ij(Ns,Ns,nh2,nh2))
  call build_cdgc_2p_states(cc_ij)
  !
  
  
  !
  !+- get all the double well wavefunctions and the dipole matrix elements -+!
  !
  xlmin=-Lwell-dwell*0.5d0
  xlmax=xlmin+Lwell  
  xrmin=dwell*0.5d0
  xrmax=xrmin+Lwell
  !
  allocate(xr(2*Nxw+Nxd),wtr(2*Nxw+Nxd)) 
  xr(1:Nxw)=linspace(xlmin,xlmax,Nxw,iend=.true.,istart=.true.)
  wtr(1:Nxw)=xr(2)-xr(1)
  !
  xr(Nxw+1:Nxw+Nxd)=linspace(xlmax,xrmin,Nxd,iend=.false.,istart=.false.)
  wtr(Nxw+1:Nxw+Nxd)=xr(Nxw+2)-xr(Nxw+1)
  !
  xr(Nxw+Nxd+1:2*Nxw+Nxd)=linspace(xrmin,xrmax,Nxw,iend=.true.,istart=.true.)
  wtr(Nxw+Nxd+1:2*Nxw+Nxd)=xr(Nxw+Nxd+2)-xr(Nxw+Nxd+1)
  !
  Nxw=Nxw*2+Nxd
  !
  !+- initialize the two-wells spectrum -+!
  out_unit=free_unit()  
  Lwell=Lwell
  allocate(eWells(Norb,2))
  open(unit=out_unit,file='dw_spectrum.data')
  do isite=1,2
     do iorb=1,Norb
        eWells(iorb,isite)=hbarc**2.d0/2.d0/mel/Lwell/Lwell*pi*pi*dble(iorb)**2.0
        write(out_unit,*) iorb,isite,eWells(iorb,isite)
     end do
  end do
  close(out_unit)
  !
  open(unit=out_unit,file='dw_wavefunctions.data')
  !
  allocate(psi_w(2,Norb,Nxw)); psi_w=0.d0
  allocate(dpsi_w(2,Norb,Nxw)); dpsi_w=0.d0
  do iorb=1,Norb
     psi_w(:,iorb,:)=0.d0
     do ix=1,Nxw
        if(xr(ix).ge.xlmin.and.xr(ix).le.xlmax) then
           kn=pi*dble(iorb)/Lwell
           psi_w(1,iorb,ix) = sqrt(2.d0/Lwell)*sin(kn*(xr(ix)-xlmin))
           dpsi_w(1,iorb,ix) = kn*sqrt(2.d0/Lwell)*cos(kn*(xr(ix)-xlmin))
        end if
        if(xr(ix).ge.xrmin.and.xr(ix).le.xrmax) then
           kn=pi*dble(iorb)/Lwell
           psi_w(2,iorb,ix) = sqrt(2.d0/Lwell)*sin(kn*(xr(ix)-xrmin))
        end if
        write(out_unit,'(10F18.10)') xr(ix),psi_w(1,iorb,ix),psi_w(2,iorb,ix)
     end do
     write(out_unit,*)
     write(out_unit,*)
  end do
  close(out_unit)
  !

  !+- get all the matrix elements -+!
  open(unit=out_unit,file='dw_pmatrix_elements.data')  
  allocate(pnn(Ns,Ns))
  do is=1,Ns
     do js=1,Ns
        pnn(is,js)=p_matrix_element(is,js)
        write(out_unit,'(2I4,10F18.10)') is,js,pnn(is,js) 
     end do
  end do
  !
  close(out_unit)



  !+- try to compute U and V -+!
  ! allocate(Uorb(Norb,Norb),Vorb(Norb,Norb))
  ! do iorb=1,Norb
  !    do jorb=1,Norb
  !       !
  !       Uorb(iorb,jorb) = 0.d0
  !       Vorb(iorb,jorb) = 0.d0
  !       !
  !       do ix=1,Nxw
  !          Utmp=0.d0
  !          do jx=1,Nxw
  !             xx=xr(ix)
  !             yy=xr(jx)          
  !             if(ix.ne.jx) then
  !                Utmp = Utmp + psi_w(1,iorb,jx)**2.d0/abs(xx-yy)*wtr(jx)*wtr(ix)
  !                !Vtmp = Vtmp + psi_w(1,iorb,jx)**2.d0/abs(xx-yy)*wtr(jx)*wtr(ix)
  !             end if
  !          end do
  !          !write(111,*) xx,Utmp,Utmp*psi_w(1,jorb,ix)**2.d0
  !          Uorb(iorb,jorb) = Uorb(iorb,jorb) + Utmp*psi_w(1,jorb,ix)**2.d0*wtr(ix)
  !          Vorb(iorb,jorb) = Vorb(iorb,jorb) + Utmp*psi_w(2,jorb,ix)**2.d0*wtr(ix)
  !       end do
  !       write(*,*) iorb,jorb,Uorb(iorb,jorb),Vorb(iorb,jorb)
  !    end do
  ! end do
  ! stop

  
  allocate(Hdw(nh2,nh2));  
  Hdw=0.d0
  do is=1,Ns
     ispin=ios_i(1,is)
     iorb=ios_i(2,is)
     isite=ios_i(3,is)
     !
     do i2=1,nh2
        Hdw(i2,i2)=Hdw(i2,i2)+eWells(iorb,isite)*cc_ij(is,is,i2,i2)
        
        do js=1,Ns
           jspin=ios_i(1,js)
           jorb=ios_i(2,js)
           jsite=ios_i(3,js)
           !
           if(iorb.le.int_cut.and.jorb.le.int_cut) then           
              if(isite.eq.jsite) then
                 !U
                 if(iorb.ne.jorb) then
                    Hdw(i2,i2) = Hdw(i2,i2) + 0.5d0*Uw*cc_ij(is,is,i2,i2)*cc_ij(js,js,i2,i2)
                 else
                    if(ispin.ne.jspin) then
                       Hdw(i2,i2) = Hdw(i2,i2) + 0.5d0*Uw*cc_ij(is,is,i2,i2)*cc_ij(js,js,i2,i2)
                    end if
                 end if
              else
                 Hdw(i2,i2) = Hdw(i2,i2) + 0.5d0*Vw*cc_ij(is,is,i2,i2)*cc_ij(js,js,i2,i2)
              end if
           end if
           !
        end do
     end do
  end do

  
  open(unit=out_unit,file='2p_spectrum.data')
  uio=free_unit()
  open(unit=uio,file='2p_eigenstates.data')  
  allocate(Edw(nh2))
  call eigh(Hdw,Edw)
  do i2=1,nh2
     write(out_unit,'(10F18.10)') Edw(i2)
     do j2=1,nh2
        write(uio,'(2I4,10F18.10)') i2,j2,Hdw(j2,i2)        
     end do
     write(uio,'(10F18.10)')
  end do
  close(out_unit)
  close(uio)
  !
  !
  allocate(ccH_ij(Ns,Ns,nh2,nh2))
  allocate(pnn_H2(nh2,nh2)) !+- momentum operator on the two-particle basis -+!

  pnn_H2 = 0.d0
  do is=1,Ns
     do js=1,Ns
        ccH_ij(is,js,:,:) = matmul(cc_ij(is,js,:,:),Hdw)
        ccH_ij(is,js,:,:) = matmul(transpose(conjg(Hdw)),ccH_ij(is,js,:,:))
        pnn_H2 = pnn_H2 + pnn(is,js)*ccH_ij(is,js,:,:)
     end do
  end do
  !
  

  
  !+- compute left/right susceptibility at zero frequency -+!
  zeta=0.d0
  do i2=1,nh2
     zeta=zeta+exp(-beta*(Edw(i2)-Edw(1)))
  end do   
  chiLR=0.d0; 
  isite=1
  jsite=2
  do ispin=1,2
     do iorb=1,Norb           
        do jorb=1,Norb
           is=i_ios(ispin,iorb,isite)
           js=i_ios(ispin,jorb,jsite)
           !+- here compute chi_{is,js} !+-> chiLR=\sum_{is,js} chi_{is,js}
           chi_ij=0.d0
           do i2=1,nh2
              do j2=1,nh2
                 !
                 if(abs(Edw(i2)-Edw(j2)).gt.1.d-10) then
                    !
                    chi_tmp = exp(-beta*(Edw(j2)-Edw(1)))-exp(-beta*(Edw(i2)-Edw(1)))
                    chi_tmp = chi_tmp/(Edw(i2)-Edw(j2))
                    chi_tmp = chi_tmp*ccH_ij(is,js,i2,j2)**2.d0                    
                    !
                 else
                    !
                    chi_tmp = beta*exp(-beta*(Edw(i2)-Edw(1)))*ccH_ij(is,js,i2,j2)**2.d0
                    !
                 end if                 
                 chi_ij = chi_ij - chi_tmp/zeta                 
              end do
           end do
           chiLR=chiLR+chi_ij           
        end do
     end do
  end do

  open(unit=out_unit,file='chiLR_bare.data')  
  write(out_unit,*) chiLR,zeta
  close(out_unit)
  chiLR_bare=chiLR

  if(Norb.eq.1.or.(.not.lm_coupling)) stop
  !
  wph=wph*(eWells(2,1)-eWells(1,1))
  Daa=0.5d0*hbarc**2.d0/mel/leff/leff!*0.d0
  Dpa=hbarc**2.d0/mel/leff!*0.d0
  open(unit=out_unit,file='cavity_info.data')
  write(out_unit,*) 'cavity effective lenght; defined as (e A_0/\hbar) = 1/l_eff',leff
  write(out_unit,*) ' (e A_0/\hbar) [nm^{-1}]',1.d0/leff
  write(out_unit,*) ' wph ',wph
  write(out_unit,*) ' Dpa',Dpa
  write(out_unit,*) ' Daa',Daa
  close(out_unit)


  open(unit=uio,file='2p_dipole_matrix_el.data')
  do i2=1,nh2
     do j2=1,nh2
        write(uio,'(2I4,30F8.4)') i2,j2,pnn_H2(i2,j2),Dpa*pnn_H2(i2,j2)/wph
     end do
  end do
  close(uio)
  
  !+-  here write down the light-matter hamiltonian -+!
  Nh=Nph*nh2
  !
  allocate(eigv_lm(Nh))
  allocate(Hlm(Nh,Nh)); 
  !
  Hlm=0.d0
  do i2=1,nh2
     do ik=1,Nph
        !
        ilm = (i2-1)*Nph+ik
        Hlm(ilm,ilm) = Edw(i2) + dble(ik-1)*(wph+2.d0*Daa)
        !
        do j2=1,nh2
           !
           if(ik.lt.Nph) then
              !
              ilm = (i2-1)*Nph + ik 
              jlm = (j2-1)*Nph + ik + 1                 
              Hlm(ilm,jlm) = Hlm(ilm,jlm) + Dpa*pnn_H2(i2,j2)*sqrt(dble(ik))
              !              
           end if
           !
           if(ik.gt.1) then
              !
              ilm = (i2-1)*Nph + ik 
              jlm = (j2-1)*Nph + ik - 1
              Hlm(ilm,jlm) = Hlm(ilm,jlm) + Dpa*pnn_H2(i2,j2)*sqrt(dble(ik-1))
              !
           end if
           !
           if(i2.eq.j2) then
              !
              if(ik.lt.Nph-1) then
                 ilm = (i2-1)*Nph + ik 
                 jlm = (j2-1)*Nph + ik + 2                 
                 Hlm(ilm,jlm) = Hlm(ilm,jlm) + Daa*sqrt(dble(ik))*sqrt(dble(ik+1))                    
              end if
              if(ik.gt.2) then
                 ilm = (i2-1)*Nph + ik 
                 jlm = (j2-1)*Nph + ik - 2                 
                 Hlm(ilm,jlm) = Hlm(ilm,jlm) + Daa*sqrt(dble(ik))*sqrt(dble(ik-1))                    
              end if
              !
           end if
        end do
     end do
  end do


  !
  call eigh(Hlm,eigv_lm)
  !
  open(unit=out_unit,file='lm2p_spectrum.data')
  uio=free_unit()
  open(unit=uio,file='lm2p_eigenstates.data')  
  do ilm=1,Nh
     write(out_unit,'(10F18.10)') eigv_lm(ilm)
     do jlm=1,Nh
        write(uio,'(2I4,10F18.10)') ilm,jlm,Hlm(jlm,ilm)        
     end do
     write(uio,'(10F18.10)')
  end do
  close(out_unit)
  close(uio)
  


  allocate(aph_ij(Nh,Nh)); aph_ij=0.d0
  do ik=1,Nph-1     
     do i2=1,nh2
        !jk=ik+1
        ilm = (i2-1)*Nph + ik
        jlm = (i2-1)*Nph + ik+1
        !
        aph_ij(ilm,jlm) = sqrt(dble(ik)) ! 
        !
     end do
  end do
  !
  aph_ij=matmul(aph_ij,Hlm)
  aph_ij=matmul(transpose(conjg(Hlm)),aph_ij)
  !
  !+- compute the photonic greens function -+!
  allocate(wr(Lreal)); wr=linspace(wmin,wmax,Lreal)
  allocate(gph(Lreal)); gph=0.d0


  zeta=0.d0
  do ilm=1,Nh
     zeta=zeta+exp(-beta*(eigv_lm(ilm)-eigv_lm(1)))
  end do
  
  open(unit=out_unit,file='ph_gf.data')    
  do iw=1,Lreal
     gph(iw)=0.d0
     do ilm=1,Nh
        do jlm=1,Nh
           !
           if(abs(eigv_lm(ilm)-eigv_lm(jlm)).gt.1.d-15) then
              !
              gtmp = exp(-beta*(eigv_lm(jlm)-eigv_lm(1)))-exp(-beta*(eigv_lm(ilm)-eigv_lm(1)))
              gtmp = gtmp/(wr(iw)+eigv_lm(ilm)-eigv_lm(jlm)+xi*0.02)
              gtmp = gtmp*abs(aph_ij(ilm,jlm))**2.d0                    
              !
           else
              !
              gtmp = beta*exp(-beta*(eigv_lm(ilm)-eigv_lm(1)))*abs(aph_ij(ilm,jlm))**2.d0
              !
           end if
           gph(iw) = gph(iw) - gtmp/zeta                 
        end do
     end do
     write(out_unit,'(10F18.10)') wr(iw),gph(iw)     
  end do
  close(out_unit)
  
  !+-> recompute the chi without storing the cc matrix elements and by cutting-off the spectrum -+!
  zeta=0.d0
  do ilm=1,Nh
     if(exp(-beta*(eigv_lm(ilm)-eigv_lm(1))).gt.threshold) then
        zeta = zeta + exp(-beta*(eigv_lm(ilm)-eigv_lm(1)))
     else
        Nh_cut=ilm-1
        exit
     end if
  end do
  
  write(*,*) 'spectrum cut-off for computing chiLR',Nh_cut,Nh
  !
  !
  chiLR=0.d0
  isite=1
  jsite=2
  do ispin=1,2
     do iorb=1,Norb           
        do jorb=1,Norb
           is=i_ios(ispin,iorb,isite)
           js=i_ios(ispin,jorb,jsite)
           !+- here compute chi_{is,js} !+-> chiLR=\sum_{is,js} chi_{is,js}
                     
           chi_ij=0.d0
           do ilm=1,Nh_cut
              do jlm=1,Nh
                 !
                 chi_tmp=0.d0
                 if(abs(eigv_lm(ilm)-eigv_lm(jlm)).gt.1.d-12) then
                    !
                    chi_tmp = -exp(-beta*(eigv_lm(ilm)-eigv_lm(1)))
                    chi_tmp = chi_tmp/(eigv_lm(ilm)-eigv_lm(jlm))
                    !
                    !+- here I should compute the ccHlm_ij matrix elements w/out storing it before -+!
                    cc_tmp=0.d0
                    do i2=1,nh2
                       do j2=1,nh2
                          do ik=1,Nph
                             iilm = (i2-1)*Nph + ik
                             jjlm = (j2-1)*Nph + ik
                             cc_tmp = cc_tmp + conjg(Hlm(iilm,ilm))*ccH_ij(is,js,i2,j2)*Hlm(jjlm,jlm)                          
                          end do
                       end do
                    end do
                    !
                    chi_tmp = chi_tmp*abs(cc_tmp)**2.d0
                 end if
                 chi_ij = chi_ij - chi_tmp/zeta
                 !
              end do
           end do
           chiLR=chiLR+chi_ij


           
           chi_ij=0.d0
           do ilm=1,Nh
              do jlm=1,Nh_cut
                 !
                 chi_tmp=0.d0
                 if(abs(eigv_lm(ilm)-eigv_lm(jlm)).gt.1.d-12) then
                    !
                    chi_tmp = exp(-beta*(eigv_lm(jlm)-eigv_lm(1)))
                    chi_tmp = chi_tmp/(eigv_lm(ilm)-eigv_lm(jlm))
                    !
                    !+- here I should compute the ccHlm_ij matrix elements w/out storing it before -+!
                    cc_tmp=0.d0
                    do i2=1,nh2
                       do j2=1,nh2
                          do ik=1,Nph
                             iilm = (i2-1)*Nph + ik
                             jjlm = (j2-1)*Nph + ik
                             cc_tmp = cc_tmp + conjg(Hlm(iilm,ilm))*ccH_ij(is,js,i2,j2)*Hlm(jjlm,jlm)                          
                          end do
                       end do
                    end do
                    chi_tmp = chi_tmp*abs(cc_tmp)**2.d0
                 end if
                 chi_ij = chi_ij - chi_tmp/zeta
                 !
              end do
           end do
           chiLR=chiLR+chi_ij
           !
           chi_ij=0.d0
           do ilm=1,Nh_cut
              do jlm=1,Nh_cut
                 !
                 chi_tmp=0.d0
                 if(abs(eigv_lm(ilm)-eigv_lm(jlm)).lt.1.d-12) then
                    !
                    chi_tmp = beta*exp(-beta*(eigv_lm(ilm)-eigv_lm(1)))
                    !
                    !+- here I should compute the ccHlm_ij matrix elements w/out storing it before -+!
                    cc_tmp=0.d0
                    do i2=1,nh2
                       do j2=1,nh2
                          do ik=1,Nph
                             iilm = (i2-1)*Nph + ik
                             jjlm = (j2-1)*Nph + ik
                             cc_tmp = cc_tmp + conjg(Hlm(iilm,ilm))*ccH_ij(is,js,i2,j2)*Hlm(jjlm,jlm)                          
                          end do
                       end do
                    end do
                    chi_tmp = chi_tmp*abs(cc_tmp)**2.d0
                    !
                 end if
                 chi_ij = chi_ij - chi_tmp/zeta
              end do
           end do
           chiLR=chiLR+chi_ij
        end do
     end do
  end do
  
  open(unit=out_unit,file='chiLR.data')  
  write(out_unit,*) chiLR,zeta,chiLR_bare
  close(out_unit)

  
  stop
  

  

  
  ! stop
  
  ! allocate(ccHlm_ij(Ns,Ns,Nh,Nh)) !+- cc-operator transformation -+!
  ! do is=1,Ns
  !    do js=1,Ns
  !       do ilm=1,Nh
  !          do jlm=1,Nh
  !             ccHlm_ij(is,js,ilm,jlm)=0.d0

  !             do ik=1,Nph
  !                do i2=1,nh2
  !                   do j2=1,nh2
  !                      iilm = (i2-1)*Nph + ik
  !                      jjlm = (j2-1)*Nph + ik
  !                      !
  !                      ccHlm_ij(is,js,ilm,jlm) = ccHlm_ij(is,js,ilm,jlm) + &
  !                           conjg(Hlm(iilm,ilm))*ccH_ij(is,js,i2,j2)*Hlm(jjlm,jlm)
  !                      !
  !                   end do
  !                end do
  !             end do
              
  !          end do
  !       end do                
  !    end do
  ! end do
  

  
  
  
  !
CONTAINS


  ! function x_matrix_element(ik,jk) result(Hpk)
  !   integer :: ik,jk
  !   complex(8) :: Hpk
  !   real(8) :: k,kk
  !   integer ::ix
  !   !
  !   Hpk=0.d0
  !   do ix=1,Nx
  !      !
  !      k=sqrt(2.d0/(xmax-xmin))*cos(dble(ik)*pi*xr(ix)/(xmax-xmin))
  !      if(mod(ik,2).eq.0) then
  !         k=sqrt(2.d0/(xmax-xmin))*sin(dble(ik)*pi*xr(ix)/(xmax-xmin))
  !      end if
  !      !
  !      kk=sqrt(2.d0/(xmax-xmin))*cos(dble(jk)*pi*xr(ix)/(xmax-xmin))
  !      if(mod(jk,2).eq.0) then
  !         kk=sqrt(2.d0/(xmax-xmin))*sin(dble(jk)*pi*xr(ix)/(xmax-xmin))
  !      end if
  !      Hpk=Hpk+k*kk*xr(ix)*wtr(ix)
  !   end do
  !   !
  !   !write(111,'(10F18.10)') dble(ik),dble(jk),Hpk
  ! end function x_matrix_element

  ! function x2_matrix_element(ik,jk) result(Hpk)
  !   integer :: ik,jk
  !   complex(8) :: Hpk
  !   real(8) :: k,kk
  !   integer ::ix
  !   !
  !   Hpk=0.d0
  !   do ix=1,Nx
  !      !
  !      k=sqrt(2.d0/(xmax-xmin))*cos(dble(ik)*pi*xr(ix)/(xmax-xmin))
  !      if(mod(ik,2).eq.0) then
  !         k=sqrt(2.d0/(xmax-xmin))*sin(dble(ik)*pi*xr(ix)/(xmax-xmin))
  !      end if

  !      kk=sqrt(2.d0/(xmax-xmin))*cos(dble(jk)*pi*xr(ix)/(xmax-xmin))
  !      if(mod(jk,2).eq.0) then
  !         kk=sqrt(2.d0/(xmax-xmin))*sin(dble(jk)*pi*xr(ix)/(xmax-xmin))
  !      end if
  !      Hpk=Hpk+k*kk*xr(ix)*xr(ix)*wtr(ix)
  !   end do
  !   !
  !   !write(111,'(10F18.10)') dble(ik),dble(jk),Hpk
  ! end function x2_matrix_element

  
  function p_matrix_element(is,js) result(Hpk)
    integer :: is,js
    integer ::ix
    complex(8) :: Hpk,tmp

    integer :: isite,jsite,iorb,jorb,ispin,jspin

    ispin=ios_i(1,is)
    jspin=ios_i(1,js)
    !
    iorb=ios_i(2,is)
    jorb=ios_i(2,js)
    !
    isite=ios_i(3,is)
    jsite=ios_i(3,js)

    
    Hpk=0.d0
    if(ispin.eq.jspin) then
       do ix=1,Nxw-1 !+-> minus one to compute the derivative numerically. (Last point is zero anyway)
          !tmp=psi_w(isite,iorb,ix)*(psi_w(jsite,jorb,ix+1)-psi_w(jsite,jorb,ix))/(xr(ix+1)-xr(ix))*wtr(ix)
          tmp=psi_w(isite,iorb,ix)*dpsi_w(jsite,jorb,ix)*wtr(ix)
          Hpk = Hpk + tmp
       end do
    end if
    Hpk=-xi*Hpk
  end function p_matrix_element

  !
  subroutine build_lattice_model  
    implicit none
    !
    integer                          :: ix,iy,iz,ik,Nk,iorb,jorb,ispin,istate,jstate
    real(8),allocatable,dimension(:) :: kx
    real(8)                          :: ts,test_k,kx_,ky_,kz_,wini,wfin,de,n1,n2
    !
    !
    Lk=Nx
    allocate(epsik(Lk),wtk(Lk),hybik(Lk))    
    wini=-Wband/2.d0
    wfin= Wband/2.d0
    ! wini=-5.d0
    ! wfin= 5.d0
    epsik=linspace(wini,wfin,Lk,mesh=de)
    !
    test_k=0.d0
    do ix=1,Lk
       if(epsik(ix).gt.-Wband/2.d0.and.epsik(ix).lt.Wband/2) then
          wtk(ix)=4.d0/Wband/pi*sqrt(1.d0-(2.d0*epsik(ix)/Wband)**2.d0)*de
       else
          wtk(ix) = 0.d0
       end if
       !wtk(ix) = 1.d0/Wband*de
       !if(ix==1.or.ix==Lk) wtk(ix)=0.d0
       test_k=test_k+wtk(ix)
       write(77,*) epsik(ix),wtk(ix)
    end do
    hybik=0.d0
    ! write(*,*) test_k,de
    !
    ! allocate(kx(Nx))
    ! kx = linspace(0.d0,pi,Nx,.true.,.true.)
    ! Lk=Nx*Nx*Nx
    ! allocate(epsik(Lk),wtk(Lk),hybik(Lk))
    ! ik=0
    ! test_k=0.d0;n1=0.d0;n2=0.d0
    ! do ix=1,Nx
    !    do iy=1,Nx
    !       do iz=1,Nx
    !          ik=ik+1
    !          !kx_=dble(ix)/dble(Nx)*pi
    !          epsik(ik) = -2.d0/6.d0*(cos(kx(ix))+cos(kx(iy))+cos(kx(iz))) 
    !          hybik(ik) = 1.d0/6.d0*(cos(kx(ix))-cos(kx(iy)))*cos(kx(iz)) 
    !          wtk(ik) = 1.d0/dble(Lk)
    !          n1=n1+fermi(epsik(ik)+cfield*0.5,beta)*wtk(ik)
    !          n2=n2+fermi(epsik(ik)-cfield*0.5,beta)*wtk(ik)
    !       end do
    !    end do
    ! end do
    ! !
    ! write(*,*) 'n1/n2'
    ! write(*,*) n1,n2,n1+n2
    call get_free_dos(epsik,wtk,file='DOS_free.kgrid')
    !
    allocate(Hk_tb(Ns,Ns,Lk))
    !
    Hk_tb=0.d0
    do ik=1,Lk
       do iorb=1,Norb
          do jorb=1,Norb
             do ispin=1,2
                istate=index(ispin,iorb)
                jstate=index(ispin,jorb)
                if(iorb.eq.jorb)  then
                   Hk_tb(istate,jstate,ik) = epsik(ik)
                else
                   Hk_tb(istate,jstate,ik) = hybik(ik)
                end if
             end do
          end do
       end do
    end do
    !
  end subroutine build_lattice_model


  subroutine getHk(Hk,ik,time)
    complex(8),dimension(:,:) :: Hk
    integer                   :: ik
    real(8)                   :: time
    if(size(Hk,1).ne.size(Hk,2)) stop "wrong dimenions in getHk"
    if(size(Hk,1).ne.Ns) stop "wrong dimenions in getHk"
    Hk = Hk_tb(:,:,ik)
  end subroutine getHk







  subroutine print_output(vdm_simplex,vdm_opt)
    real(8),dimension(Ns+1,Ns),optional :: vdm_simplex
    real(8),dimension(Nvdm_NC_opt-Nvdm_NCoff_opt),optional :: vdm_opt
    integer :: out_unit,istate,jstate,iorb,iphi,ifock,jfock,is,js,ik
    integer,dimension(Ns) :: fock_state
    complex(8),dimension(Ns) :: tmp
    real(8) :: deltani,delta_tmp,vdm_tmp

    real(8),dimension(nFock,nFock) :: test_full_phi

    !+- STORE GZ PROJECTORS -+!
    out_unit=free_unit()
    open(out_unit,file='optimized_projectors.data')
    test_full_phi=0.d0
    do iphi=1,Nphi
       !
       test_full_phi = test_full_phi + GZ_vector(iphi)*phi_basis(iphi,:,:)
       write(out_unit,'(2F18.10)') GZ_vector(iphi)
    end do
    write(out_unit,*) '!+-----------------------------+!'
    write(out_unit,*) '!+-----------------------------+!'
    write(out_unit,*) '!+-----------------------------+!'
    do ifock=1,nFock
       do jfock=1,nFock
          write(out_unit,*) test_full_phi(ifock,jfock),ifock,jfock
       end do
    end do
    close(out_unit)    

    !+- STORE SLATER DETERMINANT -+!
    do is=1,Ns
       do js=1,Ns
          out_unit=free_unit()
          open(out_unit,file='optimized_slater_normal_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data')
          do ik=1,Lk
             write(out_unit,'(2F18.10)') dreal(GZ_opt_slater(is,js,ik)),dimag(GZ_opt_slater(is,js,ik))
          end do
          close(out_unit)
       end do
    end do


    !
    out_unit=free_unit()
    open(out_unit,file='optimized_internal_energy.data')
    write(out_unit,'(5F18.10)') GZ_opt_energy,GZ_opt_kinetic,GZ_opt_Eloc
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_variational_density_matrix.data')
    write(out_unit,*) 'NORMAL VDM'
    do istate=1,Ns
       tmp(istate)=GZ_opt_VDM(istate,istate)
       write(out_unit,'(20F18.10)') GZ_opt_VDM(istate,:)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') tmp(1:Ns)
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_Rhop_matrix.data')
    do istate=1,Ns
       tmp(istate)=GZ_opt_Rhop(istate,istate)
       write(out_unit,'(20F18.10)') GZ_opt_Rhop(istate,1:Ns)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') tmp(1:Ns)
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_density.data')
    do istate=1,Ns
       write(out_unit,'(20F18.10)') gz_dens_matrix(istate,1:Ns)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') gz_dens(1:Ns)    
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='local_density_density.data')
    do is=1,Ns
       write(out_unit,'(20F18.10)') gz_dens_dens(is,:)
    end do
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='local_angular_momenta.data')
    write(out_unit,'(20F18.10)') gz_spin2,gz_spinZ,gz_isospin2,gz_isospinZ
    close(out_unit)
    !
    if(present(vdm_simplex)) then
       out_unit=free_unit()
       open(out_unit,file='vdm_simplex.restart')
       do jstate=1,Ns+1
          if(jstate.le.Ns) then
             do istate=1,Ns
                write(out_unit,'(20F18.10)') vdm_simplex(jstate,istate)
             end do
             if(jstate.le.Ns) write(out_unit,*)  'x'
          else
             do istate=1,Ns
                deltani=vdm_simplex(jstate,istate)-0.5
                if(deltani.gt.0.d0) then
                   delta_tmp=0.9999-vdm_simplex(jstate,istate)
                   vdm_tmp=vdm_simplex(jstate,istate)+delta_tmp*0.1
                   write(out_unit,'(20F18.10)') vdm_tmp
                else
                   delta_tmp=vdm_simplex(jstate,istate)-0.0001
                   vdm_tmp=vdm_simplex(jstate,istate)-delta_tmp*0.1
                   write(out_unit,'(20F18.10)') vdm_tmp
                end if
             end do
          end if
       end do
       close(out_unit)
    end if
    !
    if(present(vdm_opt)) then
       out_unit=free_unit()
       open(out_unit,file='vdm_seed.restart')
       do istate=1,Nvdm_NC_opt-Nvdm_NCoff_opt
          write(out_unit,'(10F18.10)')  vdm_opt(istate)
       end do
       close(out_unit)
    end if

  end subroutine print_output







  subroutine Rhop_vec2mat(Rhop_indep,Rhop_mat)
    complex(8),dimension(:)   :: Rhop_indep
    complex(8),dimension(:,:) :: Rhop_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(Rhop_mat,1).ne.size(Rhop_mat,2)) stop "wrong stride"
    if(size(Rhop_mat,1).ne.Ns) stop "wrong stride"
    if(size(Rhop_indep).ne.NRhop_opt) stop "wrong stride!"    
    Rhop_mat = zero
    do iorb=1,Norb
       do ispin=1,2
          is=index(ispin,iorb)
          Rhop_mat(is,is) = Rhop_indep(iorb)
       end do
    end do
    !
  end subroutine Rhop_vec2mat
  subroutine Rhop_mat2vec(Rhop_mat,Rhop_indep)
    complex(8),dimension(:,:) :: Rhop_mat
    complex(8),dimension(:)   :: Rhop_indep
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    complex(8) :: test_stride
    real(8) :: test
    if(size(Rhop_mat,1).ne.size(Rhop_mat,2)) stop "wrong stride"
    if(size(Rhop_mat,1).ne.Ns) stop "wrong stride"
    if(size(Rhop_indep).ne.NRhop_opt) stop "wrong stride!"    
    !
    ispin=1;iorb=1;is=index(ispin,iorb)
    Rhop_indep(1)=Rhop_mat(is,is)
    ispin=1;iorb=2;is=index(ispin,iorb)
    Rhop_indep(2)=Rhop_mat(is,is)
    !
  end subroutine Rhop_mat2vec



  subroutine vdm_NC_vec2mat(vdm_NC_indep,vdm_NC_mat)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NC_opt) stop "wrong stride!"    
    !
    vdm_NC_mat = zero
    do iorb=1,Norb
       do ispin=1,2
          is=index(ispin,iorb)
          vdm_NC_mat(is,is) = vdm_NC_indep(iorb)
       end do
    end do
    !
  end subroutine vdm_NC_vec2mat
  subroutine vdm_NC_mat2vec(vdm_NC_mat,vdm_NC_indep)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NC_opt) stop "wrong stride!"    
    !
    ispin=1;iorb=1;is=index(ispin,iorb)
    vdm_NC_indep(1)=vdm_NC_mat(is,is)
    ispin=1;iorb=2;is=index(ispin,iorb)
    vdm_NC_indep(2)=vdm_NC_mat(is,is)
    !
  end subroutine vdm_NC_mat2vec



  subroutine vdm_NCoff_vec2mat(vdm_NC_indep,vdm_NC_mat)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NCoff_opt) stop "wrong stride!"    
    !
    vdm_NC_mat = zero
    !
  end subroutine vdm_NCoff_vec2mat
  subroutine vdm_NCoff_mat2vec(vdm_NC_mat,vdm_NC_indep)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NCoff_opt) stop "wrong stride!"    
    !
    vdm_NC_indep = zero
    !
  end subroutine vdm_NCoff_mat2vec


  subroutine init_bcs_wf(bcs_wf,U,phase)
    complex(8),dimension(:,:) :: bcs_wf
    real(8) :: U,bcs_sc_order,phase
    real(8) :: Ek
    real(8) :: sintk,costk,sinph,cosph
    complex(8) :: bcs_phi
    if(size(bcs_wf,1).ne.3) stop "error init bcs \sigma"
    if(size(bcs_wf,2).ne.Lk) stop "error init bcs Lk"
    !
    Ubcs=U
    !write(*,*) bcs_self_cons(0.d0),bcs_self_cons(1.d0);stop
    bcs_sc_order=fzero_brentq(bcs_self_cons,0.d0,1.d0)    
    bcs_phi=bcs_sc_order*exp(xi*phase*pi)
    sinph=dimag(bcs_phi);cosph=dreal(bcs_phi)
    do ik=1,Lk
       !
       Ek = sqrt(epsik(ik)**2.d0 + (bcs_sc_order*Ubcs)**2.d0)
       sintk=bcs_sc_order*Ubcs/Ek
       costk=epsik(ik)/Ek
       !
       bcs_wf(1,ik) = -cosph*sintk*tanh(beta*Ek*0.5d0)
       bcs_wf(2,ik) = -sinph*sintk*tanh(beta*Ek*0.5d0)
       bcs_wf(3,ik) = -costk*tanh(beta*Ek*0.5d0)
       !
    end do
    !
  end subroutine init_bcs_wf

  function bcs_self_cons(phi) result(x)
    real(8),intent(in) :: phi
    real(8) :: x
    real(8) :: Ek
    integer :: ik
    !
    x=0.d0
    do ik=1,Lk
       Ek = sqrt(epsik(ik)**2.d0 + (phi*Ubcs)**2.d0)
       x = x + wtk(ik)/Ek*tanh(beta*Ek*0.5d0)
    end do
    x=x*abs(Ubcs)*0.5d0
    x=x-1.d0
    !
  end function bcs_self_cons

  function bcs_order_param(xU) result(phi)
    real(8),intent(in) :: xU
    real(8) :: phi
    !
    Ubcs=xU
    phi=fzero_brentq(bcs_self_cons,0.d0,1.d0)    
    !
  end function bcs_order_param

  subroutine getUbcs(phi,xu)
    real(8) :: phi,xu,tmp
    phiBCS=phi
    ! tmp=delta_bcs_order_param(-10.d0)
    ! write(*,*) phiBCS,tmp
    ! tmp=delta_bcs_order_param(-0.2d0)
    ! write(*,*) phiBCS,tmp
    xu = fzero_brentq(delta_bcs_order_param,-10.d0,-0.2d0)
    !write(*,*) xu
  end subroutine getUbcs


  function delta_bcs_order_param(xU) result(phi)
    real(8),intent(in) :: xU
    real(8) :: phi
    !
    Ubcs=xU
    phi=fzero_brentq(bcs_self_cons,0.d0,1.d0)    
    phi=phi-phiBCS
    !
  end function delta_bcs_order_param





  subroutine read_SC_order(read_dir,sc_order)
    character(len=200)             :: read_dir
    real(8)                        :: tmp,sc_order
    integer :: is,js,read_unit,flen,ios,i
    character(len=200) :: file_name
    logical :: check_file
    !
    file_name=reg(read_dir)//'local_sc_order.data'
    inquire(file=file_name,exist=check_file)
    if(check_file) then
       read_unit=free_unit()
       open(read_unit,file=file_name,status='old')
       flen=0
       do
          read (read_unit,*,iostat=ios) tmp
          if (ios/=0) exit     
          flen=flen+1
       end do
       close(read_unit)                    
       if(flen.ne.Ns) stop "READING SC ORDER PARAM: number of lines /= Ns"
       open(read_unit,file=file_name,status='old')
       do i=1,flen
          read(read_unit,'(2F18.10)') tmp
          if(i.eq.flen) sc_order = sqrt(tmp*tmp)
       end do
       close(read_unit)
    else
       write(*,*) 'FILE',file_name,'does not exist!!!'
    end if
    !
  end subroutine read_SC_order


end program GUTZ_mb



!AMOEBA TEST





  !+- compute again the left-right susceptibility -+!
  ! zeta=0.d0
  ! do ilm=1,Nh
  !    zeta=zeta+exp(-beta*(eigv_lm(ilm)-eigv_lm(1)))
  ! end do   
  ! chiLR=0.d0
  ! isite=1
  ! jsite=2
  ! do ispin=1,2
  !    do iorb=1,Norb           
  !       do jorb=1,Norb
  !          is=i_ios(ispin,iorb,isite)
  !          js=i_ios(ispin,jorb,jsite)
  !          !+- here compute chi_{is,js} !+-> chiLR=\sum_{is,js} chi_{is,js}
  !          chi_ij=0.d0
  !          do ilm=1,Nh
  !             do jlm=1,Nh
  !                !
  !                if(abs(eigv_lm(ilm)-eigv_lm(jlm)).gt.1.d-15) then
  !                   !
  !                   chi_tmp = exp(-beta*(eigv_lm(jlm)-eigv_lm(1)))-exp(-beta*(eigv_lm(ilm)-eigv_lm(1)))
  !                   chi_tmp = chi_tmp/(eigv_lm(ilm)-eigv_lm(jlm))
  !                   chi_tmp = chi_tmp*abs(ccHlm_ij(is,js,ilm,jlm))**2.d0                    
  !                   !
  !                else
  !                   !
  !                   chi_tmp = beta*exp(-beta*(eigv_lm(ilm)-eigv_lm(1)))*abs(ccHlm_ij(is,js,ilm,jlm))**2.d0
  !                   !
  !                end if                 
  !                chi_ij = chi_ij - chi_tmp/zeta                 
  !             end do
  !          end do
  !          chiLR=chiLR+chi_ij           
  !       end do
  !    end do
  ! end do
  
  ! open(unit=out_unit,file='chiLR.data')  
  ! write(out_unit,*) chiLR,zeta,chiLR_bare
  ! close(out_unit)



  ! allocate(ccHlm_ij(Ns,Ns,Nh,Nh)) !+- cc-operator transformation -+!
  ! ccHlm_ij=0.d0
  ! do ik=1,Nph
  !    do i2=1,nh2
  !       do j2=1,nh2
  !          iilm = (i2-1)*Nph + ik
  !          jjlm = (j2-1)*Nph + ik
  !          ccHlm_ij(:,:,iilm,jjlm) = ccH_ij(:,:,i2,j2)           
  !       end do
  !    end do
  ! end do
  ! !
  ! do is=1,Ns
  !    do js=1,Ns
  !       ccHlm_ij(is,js,:,:) = matmul(ccHlm_ij(is,js,:,:),Hlm)
  !       ccHlm_ij(is,js,:,:) = matmul(transpose(conjg(Hlm)),ccHlm_ij(is,js,:,:))
  !    end do
  ! end do



  !
  !
  !
  ! chiLR=0.d0
  ! isite=1
  ! jsite=2
  ! do ispin=1,2
  !    do iorb=1,Norb           
  !       do jorb=1,Norb
  !          is=i_ios(ispin,iorb,isite)
  !          js=i_ios(ispin,jorb,jsite)
  !          !+- here compute chi_{is,js} !+-> chiLR=\sum_{is,js} chi_{is,js}
                     
  !          chi_ij=0.d0
  !          do ilm=1,Nh
  !             do jlm=1,Nh
  !                !
  !                chi_tmp=0.d0
  !                if(abs(eigv_lm(ilm)-eigv_lm(jlm)).gt.1.d-12) then
  !                   !
  !                   chi_tmp = exp(-beta*(eigv_lm(jlm)-eigv_lm(1)))-exp(-beta*(eigv_lm(ilm)-eigv_lm(1)))
  !                   chi_tmp = chi_tmp/(eigv_lm(ilm)-eigv_lm(jlm))
  !                   !
  !                else
  !                   !                    
  !                   chi_tmp = beta*exp(-beta*(eigv_lm(ilm)-eigv_lm(1)))
  !                   !
  !                end if
                 
  !                !+- here I should compute the ccHlm_ij matrix elements w/out storing it before -+!
  !                cc_tmp=0.d0
  !                do i2=1,nh2
  !                   do j2=1,nh2
  !                      do ik=1,Nph
  !                         iilm = (i2-1)*Nph + ik
  !                         jjlm = (j2-1)*Nph + ik
  !                         cc_tmp = cc_tmp + conjg(Hlm(iilm,ilm))*ccH_ij(is,js,i2,j2)*Hlm(jjlm,jlm)                          
  !                      end do
  !                   end do
  !                end do
  !                !
  !                chi_tmp = chi_tmp*abs(cc_tmp)**2.d0
  !                chi_ij = chi_ij - chi_tmp/zeta
  !                !
  !             end do
  !          end do           
  !          chiLR=chiLR+chi_ij           
  !       end do
  !    end do
  ! end do
  
  ! open(unit=out_unit,file='chiLR.data')  
  ! write(out_unit,*) chiLR,zeta,chiLR_bare
  ! close(out_unit)

