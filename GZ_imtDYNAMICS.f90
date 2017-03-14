MODULE GZ_imtDYNAMICS  
  ! scifor
  USE SF_LINALG
  USE SF_OPTIMIZE
  USE SF_SPECIAL
  USE RK_IDE
  ! GZ rooutines  
  USE GZ_VARS_GLOBAL
  USE GZ_AUX_FUNX
  USE GZ_LOCAL_FOCK
  USE GZ_MATRIX_BASIS
  USE GZ_EFFECTIVE_HOPPINGS
  USE GZ_neqAUX_FUNX
  USE GZ_LOCAL_HAMILTONIAN
  USE MATRIX_SPARSE
  USE GZ_ENERGY_MINIMIZATION
  implicit none
  private
  !
  interface step_imt_dynamics_superc_d
     module procedure  step_imt_dynamics_superc_NC_AC_dens_d,step_imt_dynamics_superc_NC_dens_d,step_imt_dynamics_superc_NC_d
  end interface step_imt_dynamics_superc_d
  !
  public :: step_imt_dynamics
  public :: step_imt_dynamics_superc_d,step_imt_dynamics_superc_z
  !
  public :: gz_imt_equations_of_motion
  public :: gz_imt_equations_of_motion_superc
  !
  public :: setup_imt_hamiltonian
  public :: beta0_init_imt_qpH,beta0_init_imt_qpH_superc
  public :: beta0_entropy
  !
  public :: setup_imt_dynamics
  public :: setup_imt_dynamics_superc
  !
  public :: gz_imt_measure,gz_imt_measure_superc
  !
  public :: get_imt_local_dens
  public :: get_imt_local_dens_dens
  public :: get_imt_energies  
  public :: get_imt_dens_constr_slater
  public :: get_imt_dens_constr_gzproj
  public :: get_imt_dens_constrA_slater
  public :: get_imt_dens_constrA_gzproj
  public :: get_imt_lgrA_slater
  public :: get_imt_lgrA_gzproj
  public :: get_imt_unitary_constr
  public :: get_imt_Rhop
  public :: get_imt_Qhop
  public :: get_imt_local_angular_momenta
  public :: get_imt_local_sc
  public :: get_imt_nqp
  !  
  complex(8),dimension(:,:),allocatable :: gz_imt_local_density_matrix !
  real(8),dimension(:,:),allocatable    :: gz_imt_local_dens_dens
  real(8),dimension(4)                  :: gz_imt_local_angular_momenta
  real(8),dimension(3)                  :: gz_imt_energies
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_constr_slater
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_constr_gzproj
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_constrA_slater
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_constrA_gzproj
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_lgrA_slater
  complex(8),dimension(:,:),allocatable :: gz_imt_dens_lgrA_gzproj
  real(8)                               :: gz_imt_unitary_constr
  complex(8),dimension(:,:),allocatable :: gz_imt_Rhop         
  complex(8),dimension(:,:),allocatable :: gz_imt_Qhop         
  complex(8),dimension(:,:),allocatable :: gz_imt_local_sc_order 
  real(8),dimension(:,:),allocatable    :: gz_imt_nqp        
  !

  complex(8),dimension(:,:,:),allocatable :: gz_imt_qpH,Hksave

  complex(8),dimension(:,:),allocatable,public :: imt_lgrNC


  complex(8),dimension(:,:,:),allocatable,public :: imt_lgr_NC
  complex(8),dimension(:,:,:),allocatable,public :: imt_lgr_AC
  real(8) :: imt_lgrU
  real(8) :: imt_lgr_local_dens
  !

  real(8) :: imt_tstep
CONTAINS
  !
  include 'gz_IMTeom.f90'
  !
  subroutine setup_imt_dynamics
    allocate(gz_imt_local_density_matrix(Ns,Ns)); gz_imt_local_density_matrix = 0.d0
    allocate(gz_imt_local_dens_dens(Ns,Ns)); gz_imt_local_dens_dens = 0.d0
    gz_imt_local_angular_momenta = 0.d0
    gz_imt_energies = 0.d0
    allocate(gz_imt_dens_constr_slater(Ns,Ns)); gz_imt_dens_constr_slater = 0.d0
    allocate(gz_imt_dens_constr_gzproj(Ns,Ns)); gz_imt_dens_constr_gzproj = 0.d0
    gz_imt_unitary_constr = 0.d0
    allocate(gz_imt_Rhop(Ns,Ns)); gz_imt_Rhop = 0.d0
  end subroutine setup_imt_dynamics

  subroutine beta0_init_imt_qpH(gzproj,Hqp,lgrNC_)
    real(8) :: beta0
    complex(8),dimension(Nphi)  :: gzproj
    complex(8),dimension(Ns,Ns),optional :: lgrNC_
    complex(8),dimension(Ns,Ns,Lk) :: Hqp
    complex(8),dimension(Ns,Ns) :: lgrNC,Hk,Rhop,Rhop_hc
    real(8),dimension(Ns) :: vdm_diag,gz_dens
    real(8) :: ntarget,infT_mu
    integer :: is,js,ik
    !
    beta0=0.d0
    lgrNC = zero
    if(present(lgrNC_)) lgrNC = lgrNC_
    allocate(gz_imt_qpH(Lk,Ns,Ns)); gz_imt_qpH=zero
    !
    ntarget = 0.d0
    do is=1,Ns
       vdm_diag(is) = dreal(trace_phi_basis(gzproj,phi_traces_basis_dens(is,is,:,:)))
       gz_dens(is) = dreal(trace_phi_basis(gzproj,phi_traces_basis_local_dens(is,is,:,:)))
       ntarget = ntarget + gz_dens(is)
    end do
    ntarget = ntarget/dble(Ns)
    !
    Rhop = hopping_renormalization_normal_sp(gzproj,vdm_diag)
    Rhop_hc = conjg(transpose(Rhop))
    !    

    infT_mu=fzero_brentq(beta0_fix_mu,-Wband*0.5d0,Wband*0.5d0)
    write(*,*) infT_mu,beta0_fix_mu(infT_mu)

    do ik=1,Lk
       Hk=matmul(Hk_tb(:,:,ik),Rhop)
       Hk=matmul(Rhop_hc,Hk)
       !
       Hqp(:,:,ik) = beta0*Hk
       do is=1,Ns
          Hqp(is,is,ik) = Hqp(is,is,ik) - infT_mu
       end do
    end do
    !
  contains
    !
    function beta0_fix_mu(mu_infT) result(delta)
      real(8),intent(in) :: mu_infT
      real(8) :: delta
      complex(8),dimension(Ns,Ns) :: Hk,Hqp
      complex(8),dimension(Ns,Ns,Lk) :: slater
      real(8),dimension(Ns) :: eHk
      real(8) :: unitary_constrGZ
      integer :: i0,i,is,js,ks
      real(8) :: n_infT
      !
      n_infT=0.d0
      do ik=1,Lk
         Hk=matmul(Hk_tb(:,:,ik),Rhop)
         Hk=matmul(Rhop_hc,Hk)
         !
         Hqp = beta0*Hk
         do is=1,Ns
            Hqp(is,is) = Hqp(is,is) - mu_infT
         end do
         call matrix_diagonalize(Hqp,eHk)
         !
         !
         do is=1,Ns
            do js=1,Ns
               slater(is,js,ik) = zero
               do ks=1,Ns
                  slater(is,js,ik) = slater(is,js,ik) + conjg(Hqp(is,ks))*Hqp(js,ks)*fermi(eHk(ks),1.d0)
               end do
            end do
            n_infT=n_infT+slater(is,is,ik)*wtk(ik)
         end do
         !
         !
         !
      end do
      n_infT=n_infT/dble(Ns)
      delta = n_infT - ntarget
      !
    end function beta0_fix_mu
    !
  end subroutine beta0_init_imt_qpH

  !
  !
  subroutine setup_imt_dynamics_superc    
    allocate(gz_imt_local_density_matrix(Ns,Ns)); gz_imt_local_density_matrix = 0.d0
    allocate(gz_imt_local_dens_dens(Ns,Ns)); gz_imt_local_dens_dens = 0.d0
    gz_imt_local_angular_momenta = 0.d0
    gz_imt_energies = 0.d0
    allocate(gz_imt_dens_constr_slater(Ns,Ns)); gz_imt_dens_constr_slater = 0.d0
    allocate(gz_imt_dens_constr_gzproj(Ns,Ns)); gz_imt_dens_constr_gzproj = 0.d0
    !
    allocate(gz_imt_dens_constrA_slater(Ns,Ns)); gz_imt_dens_constrA_slater = 0.d0
    allocate(gz_imt_dens_constrA_gzproj(Ns,Ns)); gz_imt_dens_constrA_gzproj = 0.d0
    allocate(gz_imt_dens_lgrA_slater(Ns,Ns)); gz_imt_dens_lgrA_slater = 0.d0
    allocate(gz_imt_dens_lgrA_gzproj(Ns,Ns)); gz_imt_dens_lgrA_gzproj = 0.d0
    !
    gz_imt_unitary_constr = 0.d0
    allocate(gz_imt_Rhop(Ns,Ns)); gz_imt_Rhop = 0.d0
    allocate(gz_imt_Qhop(Ns,Ns)); gz_imt_Qhop = 0.d0
    allocate(gz_imt_local_sc_order(Ns,Ns)); gz_imt_local_sc_order = 0.d0
    allocate(gz_imt_nqp(Ns,Lk)); gz_imt_nqp = 0.d0
  end subroutine setup_imt_dynamics_superc
  !+-
  subroutine beta0_init_imt_qpH_superc(gzproj,Hqp,lgrNC_)
    real(8) :: beta0
    complex(8),dimension(Nphi)  :: gzproj
    complex(8),dimension(Ns,Ns),optional :: lgrNC_
    complex(8),dimension(3,Ns,Ns,Lk) :: Hqp
    complex(8),dimension(Ns,Ns) :: lgrNC,Hk_tmp
    complex(8),dimension(Ns,Ns) :: Rhop,Rhop_hc,Qhop,Qhop_hc
    complex(8),dimension(2*Ns,2*Ns) :: Hk
    real(8),dimension(Ns) :: vdm_diag,gz_dens
    real(8) :: ntarget,infT_mu
    integer :: is,js,ik,iter,i,i0,Nopt
    !
    real(8),dimension(:),allocatable            :: lgr
    complex(8),dimension(:),allocatable         :: lgr_cmplx
    real(8),dimension(:),allocatable            ::   delta_out     !+- real indeendent lgr_vector -+!
    complex(8),dimension(2,Ns,Ns) :: lgr_multip
    !
    beta0=0.d0
    lgrNC = zero
    !
    allocate(gz_imt_qpH(Lk,Ns,Ns)); gz_imt_qpH=zero
    !
    do is=1,Ns
       vdm_diag(is) = dreal(trace_phi_basis(gzproj,phi_traces_basis_dens(is,is,:,:))) !+- this is the target -+!
    end do
    !
    Rhop = hopping_renormalization_normal_sp(gzproj,vdm_diag)
    Rhop_hc = conjg(transpose(Rhop))
    Qhop = hopping_renormalization_anomalous_sp(gzproj,vdm_diag)
    Qhop_hc = conjg(transpose(Qhop))    
    !    
    Nopt=2*Nvdm_NC_opt + 2*Nvdm_AC_opt;
    allocate(lgr(Nopt));allocate(delta_out(Nopt))
    lgr=0.d0
    !
    call fsolve(beta0_fix_mu,lgr,tol=1.d-10,info=iter)
    delta_out=beta0_fix_mu(lgr)
    !
    write(*,*)
    write(*,*)
    write(*,*) "initializing infinte temperature Hqp"
    write(*,*)
    write(*,*)
    write(*,*) "lgr params",lgr
    write(*,*) "error",delta_out
    !
    allocate(lgr_cmplx(Nvdm_NC_opt))
    do i=1,Nvdm_NC_opt
       lgr_cmplx(i) = lgr(i)+xi*lgr(i+Nvdm_NC_opt)
    end do
    call vdm_NC_stride_v2m(lgr_cmplx,lgr_multip(1,:,:))
    deallocate(lgr_cmplx)
    i0=2*Nvdm_NC_opt
    allocate(lgr_cmplx(Nvdm_AC_opt))  
    do i=1,Nvdm_AC_opt
       lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_AC_opt)
    end do
    call vdm_AC_stride_v2m(lgr_cmplx,lgr_multip(2,:,:))  
    !
    if(present(lgrNC_)) lgrNC_=lgr_multip(1,:,:)
    !
    do ik=1,Lk
       !
       Hk_tmp=matmul(Hk_tb(:,:,ik),Rhop)
       Hk_tmp=matmul(Rhop_hc,Hk_tmp)
       HK(1:Ns,1:Ns) = Hk_tmp 
       !
       Hk_tmp=matmul(Hk_tb(:,:,ik),Qhop)
       Hk_tmp=matmul(Rhop_hc,Hk_tmp)
       Hk(1:Ns,1+Ns:2*Ns) = Hk_tmp 
       Hk(1+Ns:2*Ns,1:Ns) = conjg(transpose(Hk_tmp)) 
       !
       Hk_tmp=matmul(Hk_tb(:,:,ik),Qhop)
       Hk_tmp=matmul(Qhop_hc,Hk_tmp)
       Hk(1+Ns:2*Ns,1+Ns:2*Ns) = Hk_tmp
       !
       !
       Hqp(1,:,:,ik) = beta0*Hk(1:Ns,1:Ns)  + lgr_multip(1,:,:)
       Hqp(2,:,:,ik) = beta0*Hk(1:Ns,1+Ns:2*Ns)  + lgr_multip(2,:,:)
       Hqp(3,:,:,ik) = beta0*Hk(1+Ns:2*Ns,1+Ns:2*Ns) 
       !
       !
    end do
    !
  contains
    !
    function beta0_fix_mu(lm_) result(delta)
      real(8),dimension(:)                :: lm_
      real(8),dimension(size(lm_))        :: delta
      complex(8),dimension(:),allocatable :: lm_cmplx,delta_cmplx
      complex(8),dimension(2,Ns,Ns)       :: lm   !+- this may be also (2,Ns,Ns)
      complex(8),dimension(2,Ns,Ns)       :: delta_local_density_matrix,local_density_matrix
      complex(8),dimension(2*Ns,2*Ns)     :: Hk,Hqp
      complex(8),dimension(Ns,Ns)         :: Hk_tmp
      real(8),dimension(2*Ns)             :: ek
      integer                             :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik,imap,jmap,i,i0
      integer                             :: is,js,ks
      real(8)                             :: nqp
      !
      allocate(lm_cmplx(Nvdm_NC_opt))
      do i=1,Nvdm_NC_opt
         lm_cmplx(i) = lm_(i)+xi*lm_(i+Nvdm_NC_opt)
      end do
      call vdm_NC_stride_v2m(lm_cmplx,lm(1,:,:))    
      deallocate(lm_cmplx)
      i0=2*Nvdm_NC_opt
      allocate(lm_cmplx(Nvdm_AC_opt))
      do i=1,Nvdm_AC_opt
         lm_cmplx(i) = lm_(i0+i)+xi*lm_(i0+i+Nvdm_AC_opt)
      end do
      call vdm_AC_stride_v2m(lm_cmplx,lm(2,:,:))    
      deallocate(lm_cmplx)
      !
      local_density_matrix=0.d0
      !
      do ik=1,Lk     
         Hk=0.d0
         ek=0.d0
         !
         ! hopping renormalization !
         !
         Hk_tmp=matmul(Hk_tb(:,:,ik),Rhop)
         Hk_tmp=matmul(Rhop_hc,Hk_tmp)
         Hk(1:Ns,1:Ns) = Hk_tmp 
         !
         Hk_tmp=matmul(Hk_tb(:,:,ik),Qhop)
         Hk_tmp=matmul(Rhop_hc,Hk_tmp)
         Hk(1:Ns,Ns+1:2*Ns) = Hk_tmp 
         !
         Hk_tmp=matmul(Hk_tb(:,:,ik),Rhop)
         Hk_tmp=matmul(Qhop_hc,Hk_tmp)
         Hk(Ns+1:2*Ns,1:Ns) = Hk_tmp
         !
         Hk_tmp=matmul(Hk_tb(:,:,ik),Qhop)
         Hk_tmp=matmul(Qhop_hc,Hk_tmp)
         Hk(Ns+1:2*Ns,Ns+1:2*Ns) = Hk_tmp
         !
         Hqp = beta0*Hk
         !add Lagrange multipliers !
         Hqp(1:Ns,1:Ns)=Hqp(1:Ns,1:Ns)+lm(1,:,:)
         Hqp(1:Ns,Ns+1:2*Ns)=Hqp(1:Ns,1+Ns:2*Ns)+lm(2,:,:)  !+--> probs here???
         do is=1,Ns
            do js=1,Ns
               Hqp(is+Ns,js) = Hqp(is+Ns,js) + conjg(lm(2,js,is))
            end do
         end do
         !
         !
         Hk=Hqp
         call matrix_diagonalize(Hk,ek)
         !
         !
         do is=1,Ns
            do js=1,Ns
               !
               do ks=1,Ns
                  nqp = fermi(ek(ks)-ek(ks+Ns),1.d0)                
                  !
                  local_density_matrix(1,is,js) = local_density_matrix(1,is,js) + &
                       conjg(Hk(is,ks))*Hk(js,ks)*nqp*wtk(ik) + conjg(Hk(is,ks+Ns))*Hk(js,ks+Ns)*(1.d0-nqp)*wtk(ik)
                  !
                  local_density_matrix(2,is,js) = local_density_matrix(2,is,js) + &
                       conjg(Hk(is+Ns,ks))*Hk(js,ks)*nqp*wtk(ik) + conjg(Hk(is+Ns,ks+Ns))*Hk(js,ks+Ns)*(1.d0-nqp)*wtk(ik)
               end do
               !
            end do
         end do
         !
         !
      end do
      write(*,*) local_density_matrix(1,:,:)
      delta_local_density_matrix = local_density_matrix
      do istate=1,Ns
         delta_local_density_matrix(1,istate,istate) = delta_local_density_matrix(1,istate,istate) - vdm_diag(istate)      
      end do
      delta=0.d0
      allocate(delta_cmplx(Nvdm_NC_opt))
      call vdm_NC_stride_m2v(delta_local_density_matrix(1,:,:),delta_cmplx)
      do i=1,Nvdm_NC_opt
         delta(i) = dreal(delta_cmplx(i))
         delta(i+Nvdm_NC_opt) = dimag(delta_cmplx(i))
      end do
      deallocate(delta_cmplx)
      i0 = 2*Nvdm_NC_opt
      allocate(delta_cmplx(Nvdm_AC_opt))
      call vdm_AC_stride_m2v(delta_local_density_matrix(2,:,:),delta_cmplx)
      do i=1,Nvdm_AC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
         delta(i0+i+Nvdm_AC_opt) = dimag(delta_cmplx(i))       
      end do
      deallocate(delta_cmplx)
      !
    end function beta0_fix_mu
    !
  end subroutine beta0_init_imt_qpH_superc






  subroutine beta0_entropy(beta_,gzproj,Hqp,entropy)
    real(8) :: entropy,slater_entropy,gzproj_entropy,beta_
    complex(8),dimension(Nphi)  :: gzproj
    complex(8),dimension(3,Ns,Ns,Lk) :: Hqp
    complex(8),dimension(2,Ns,Ns,Lk) :: slater
    complex(8),dimension(Ns,Ns) :: lgrNC,Hk_tmp
    complex(8),dimension(Ns,Ns) :: Rhop,Rhop_hc,Qhop,Qhop_hc
    complex(8),dimension(2*Ns,2*Ns) :: Hk
    real(8),dimension(Ns) :: vdm_diag,gz_dens
    real(8) :: ntarget,infT_mu,nqp
    integer :: is,js,ik,iter,i,i0,Nopt,ks
    !
    real(8),dimension(2*Ns)            ::    ek
    complex(8),dimension(2,Ns,Ns) :: lgr_multip
    complex(8),dimension(Nfock,Nfock) :: full_phi
    real(8),dimension(Nfock,Nfock) :: rho_phi
    !
    integer :: ifock,jfock,nstate,iphi
    integer,dimension(Ns) :: ivec
    real(8),dimension(NFock) :: P0

    type(sparse_matrix_csr_z)  :: phik_tmp
    real(8),dimension(Nphi) :: traces




    !+- occupation of single particle states -+!
    do is=1,Ns
       vdm_diag(is) = dreal(trace_phi_basis(gzproj,phi_traces_basis_dens(is,is,:,:))) 
    end do
    !+- variational occupation of local many-body states -+!
    do ifock=1,NFock
       P0(ifock)=1.d0
       !Fock(ifock)=ifock
       call bdecomp(ifock,ivec)
       call get_state_number(ivec,jfock)
       write(*,'(10I3)') ivec(:)
       do is=1,Ns
          nstate = ivec(is)
          P0(ifock) = P0(ifock)*vdm_diag(is)**dble(nstate)*(1.d0-vdm_diag(is))**dble(1-nstate)
       end do
       write(*,*) P0(ifock)
    end do
    !
    do iphi=1,Nphi
       traces(iphi)=zero
       call sp_load_matrix(phi_basis(iphi,:,:),phik_tmp)
       do i=1,phik_tmp%Nnz
          traces(iphi) = traces(iphi) + conjg(phik_tmp%values(i))*phik_tmp%values(i)
       end do
       call sp_delete_matrix(phik_tmp)
    end do
    full_phi=0.d0
    do iphi=1,Nphi
       !                                                                                                                                                                            
       full_phi = full_phi + gzproj(iphi)*phi_basis(iphi,:,:)/sqrt(traces(iphi))
    end do
    !
    rho_phi=dreal(matmul(conjg(transpose(full_phi)),full_phi))    
    gzproj_entropy=0.d0
    do ifock=1,Nfock
       do jfock=1,Nfock
          if(rho_phi(ifock,jfock)/=0.d0) then
             gzproj_entropy = gzproj_entropy - rho_phi(ifock,jfock)*log(rho_phi(jfock,ifock)/P0(jfock))
          end if
       end do
    end do
    !    
    ! !
    slater_entropy=0.d0
    slater=zero
    do ik=1,Lk
       !
       HK(1:Ns,1:Ns) = Hqp(1,:,:,ik)
       Hk(1:Ns,1+Ns:2*Ns) = conjg(transpose(Hqp(2,:,:,ik))) 
       Hk(1+Ns:2*Ns,1:Ns) = Hqp(2,:,:,ik)
       Hk(1+Ns:2*Ns,1+Ns:2*Ns) = Hqp(3,:,:,ik)
       !
       call matrix_diagonalize(Hk,ek)
       !
       do is=1,Ns
          slater_entropy = slater_entropy + log(exp(-ek(is))+exp(-ek(is+Ns)))*wtk(ik)
       end do
       !       
       do is=1,Ns
          do js=1,Ns
             !
             do ks=1,Ns
                nqp = fermi(ek(ks)-ek(ks+Ns),1.d0)                
                !
                slater(1,is,js,ik) = slater(1,is,js,ik) + &
                     conjg(Hk(is,ks))*Hk(js,ks)*nqp + conjg(Hk(is,ks+Ns))*Hk(js,ks+Ns)*(1.d0-nqp)
                !
                slater(2,is,js,ik) = slater(2,is,js,ik) + &
                     conjg(Hk(is+Ns,ks))*Hk(js,ks)*nqp + conjg(Hk(is+Ns,ks+Ns))*Hk(js,ks+Ns)*(1.d0-nqp)
             end do
             !
          end do
       end do
    end do
    ! !    
    do ik=1,Lk
       do is=1,Ns
          do js=1,Ns
             !
             slater_entropy = slater_entropy + beta_*Hqp(1,is,js,ik)*slater(1,is,js,ik)*wtk(ik)
             slater_entropy = slater_entropy + beta_*conjg(Hqp(2,is,js,ik))*conjg(slater(2,js,is,ik))*wtk(ik)
             slater_entropy = slater_entropy + beta_*Hqp(2,is,js,ik)*slater(2,is,js,ik)*wtk(ik)
             slater_entropy = slater_entropy - beta_*Hqp(3,is,js,ik)*slater(1,js,is,ik)*wtk(ik)
             if(is.eq.js) then
                slater_entropy = slater_entropy + beta_*Hqp(3,is,js,ik)*wtk(ik)
             end if
             !
          end do
       end do
    end do   
    entropy = slater_entropy + gzproj_entropy    
    write(*,*) entropy,slater_entropy,gzproj_entropy,log(4.d0),beta_
    
    !beta=10000.d0
    slater_entropy=0.d0
    do ik=1,Lk       
       nqp = fermi(Hk_tb(1,1,ik),beta)

       if(nqp.gt.1.d-10) then
          slater_entropy = slater_entropy - ( nqp*log(nqp))*wtk(ik)*2.d0
       end if
       if(abs(1.d0-nqp).gt.1.d-10) then
          slater_entropy = slater_entropy - (1.d0-nqp)*log(1.d0-nqp )*wtk(ik)*2.d0
       end if
    end do
    write(*,*) slater_entropy
    !
  end subroutine beta0_entropy






  subroutine get_imt_local_dens(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_local_dens wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_local_dens wrong indeces'
    x = gz_imt_local_density_matrix(is,js)    
  end subroutine get_imt_local_dens
  !
  subroutine get_imt_local_dens_dens(is,js,x)
    integer :: is,js
    real(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_local_dens_dens wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_local_dens_dens wrong indeces'
    x = gz_imt_local_dens_dens(is,js)    
  end subroutine get_imt_local_dens_dens
  !
  subroutine get_imt_local_angular_momenta(x)
    real(8),dimension(4) :: x
    x = gz_imt_local_angular_momenta
  end subroutine get_imt_local_angular_momenta
  !
  subroutine get_imt_local_sc(is,js,x)
    integer :: is,js
    complex(8) :: x
    x=gz_imt_local_sc_order(is,js)
  end subroutine get_imt_local_sc
  !
  subroutine get_imt_energies(x)
    real(8),dimension(3) :: x
    x = gz_imt_energies
  end subroutine get_imt_energies
  !
  subroutine get_imt_dens_constr_slater(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_constr_slater wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_constr_slater wrong indeces'
    x = gz_imt_dens_constr_slater(is,js)
  end subroutine get_imt_dens_constr_slater
  !
  subroutine get_imt_dens_constr_gzproj(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_constr_slater wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_constr_slater wrong indeces'
    x = gz_imt_dens_constr_gzproj(is,js)    
  end subroutine get_imt_dens_constr_gzproj
  !
  subroutine get_imt_dens_constrA_slater(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_constr_slater wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_constr_slater wrong indeces'
    x = gz_imt_dens_constrA_slater(is,js)
  end subroutine get_imt_dens_constrA_slater
  !
  subroutine get_imt_dens_constrA_gzproj(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_constr_gzproj wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_constr_gzproj wrong indeces'
    x = gz_imt_dens_constrA_gzproj(is,js)    
  end subroutine get_imt_dens_constrA_gzproj
  !

  subroutine get_imt_lgrA_slater(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_lgrA_sl wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_lgrA_sl wrong indeces'
    x = gz_imt_dens_lgrA_slater(is,js)
  end subroutine get_imt_lgrA_slater
  !
  subroutine get_imt_lgrA_gzproj(is,js,x)
    integer :: is,js
    complex(8) :: x
    if(is.gt.Ns.or.js.gt.Ns) stop 'get_imt_dens_lgrA_GZ wrong indeces'
    if(is.lt.1.or.js.lt.1) stop 'get_imt_dens_lgrA_GZ wrong indeces'
    x = gz_imt_dens_lgrA_gzproj(is,js)    
  end subroutine get_imt_lgrA_gzproj
  !

  subroutine get_imt_unitary_constr(x)
    real(8) :: x
    x =gz_imt_unitary_constr
  end subroutine get_imt_unitary_constr
  !
  subroutine get_imt_Rhop(is,js,x)
    implicit none
    integer :: is,js
    complex(8),intent(inout) :: x
    x=gz_imt_Rhop(is,js)
  end subroutine get_imt_Rhop
  !
  subroutine get_imt_Qhop(is,js,x)
    integer :: is,js
    complex(8) :: x
    x=gz_imt_Qhop(is,js)
  end subroutine get_imt_Qhop
  !
  subroutine get_imt_nqp(is,ik,x)
    integer :: is,ik
    real(8) :: x
    x = gz_imt_nqp(is,ik)
  end subroutine get_imt_nqp
  !
  subroutine setup_imt_hamiltonian(Uloc_t_,Ust_t_,Jh_t_,Jph_t_,Jsf_t_,eLevels_t_)
    real(8),dimension(3,Nit_aux),optional           :: Uloc_t_
    real(8),dimension(Nit_aux),optional             :: Ust_t_
    real(8),dimension(Nit_aux),optional             :: Jh_t_
    real(8),dimension(Nit_aux),optional             :: Jph_t_
    real(8),dimension(Nit_aux),optional             :: Jsf_t_
    real(8),dimension(Ns,Nit_aux),optional          :: eLevels_t_
    !
    integer                                    :: it
    !
    allocate(Uloc_t(3,Nit_aux)); forall(it=1:Nit_aux) Uloc_t(:,it) = Uloc
    if(present(Uloc_t_)) Uloc_t = Uloc_t_
    !
    allocate(Ust_t(Nit_aux)); Ust_t = Ust
    if(present(Ust_t_)) Ust_t = Ust_t_
    !
    allocate(Jh_t(Nit_aux)); Jh_t = Jh
    if(present(Jh_t_)) Jh_t = Jh_t_
    !
    allocate(Jph_t(Nit_aux)); Jph_t = Jph
    if(present(Jph_t_)) Jph_t = Jph_t_
    !
    allocate(Jsf_t(Nit_aux)); Jsf_t = Jsf
    if(present(Jsf_t_)) Jsf_t = Jsf_t_
    !
    allocate(eLevels_t(Ns,Nit_aux)); forall(it=1:Nit_aux) eLevels_t(:,it) = eLevels
    if(present(eLevels_t_)) eLevels_t = eLevels_t_
    !
  end subroutine setup_imt_hamiltonian


  subroutine gz_imt_measure(psi_t,itime,slater_out,gzproj_out,Hqp_out)
    complex(8),dimension(nDynamics) :: psi_t
    real(8)                         :: itime
    complex(8),dimension(Ns,Ns,Lk),optional  :: slater_out
    complex(8),dimension(Ns,Ns,Lk),optional  :: Hqp_out
    complex(8),dimension(Nphi),optional :: gzproj_out
    complex(8),dimension(Ns,Ns,Lk)  :: slater,tmp_slater,Hqp !+- da finire da qui
    complex(8),dimension(Ns,Ns)     :: tmpHk,Hk
    complex(8),dimension(Nphi)      :: gzproj,gztmp
    real(8)                         :: Estar,Eloc,Egz
    real(8),dimension(Ns)           :: vdm_diag,tmp_eHk

    integer                         :: is,js,ik,iphi,iis,jjs,itt,ks
    !
    if(size(psi_t).ne.nDynamics) stop "gz_imt_measure_/ wrong dimensions"
    call dynamicalVector_2_wfMatrix(psi_t,Hqp,gzproj)
    !
    if(present(gzproj_out)) gzproj_out=gzproj
    if(present(Hqp_out)) Hqp_out=Hqp
    !
    !
    do is=1,Ns
       do js=1,Ns
          gz_imt_local_density_matrix(is,js) = &
               trace_phi_basis_sp(gzproj,phi_spTraces_basis_local_dens(is,js))          
          gz_imt_local_dens_dens(is,js) = &
               trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens_dens(is,js))
          !
          gz_imt_dens_constr_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens(is,js))
          !
       end do
       vdm_diag(is) = gz_imt_dens_constr_gzproj(is,is)
    end do
    !
    itt=imt_t2it(itime,0.5d0*itstep,beta_init)
    !
    eLevels=0.d0
    !
    gztmp = get_local_hamiltonian_HLOCphi(gzproj,eLevels)
    Eloc = 0.d0
    do iphi=1,Nphi
       Eloc = Eloc + conjg(gzproj(iphi))*gztmp(iphi)
    end do
    !
    gz_imt_local_angular_momenta(1) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_spin2)
    gz_imt_local_angular_momenta(2) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_spinZ)
    gz_imt_local_angular_momenta(3) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_isoSpin2)
    gz_imt_local_angular_momenta(4) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_isoSpinZ)
    !
    gz_imt_Rhop = hopping_renormalization_normal_sp(gzproj,vdm_diag)
    !
    gz_imt_unitary_constr = 0.d0
    do iphi=1,Nphi
       gz_imt_unitary_constr = gz_imt_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
    end do
    !+- SLATER -+!
    Estar=0.d0
    !+- first thing first I have to exctract slater(is,js,ik) from the immaginary time integral of the QP hamiltonian
    do ik=1,Lk
       !
       tmpHk= Hqp(:,:,ik) 
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
       !
    end do
    if(present(slater_out)) slater_out=slater
    !
    gz_imt_dens_constr_slater=0.d0
    do ik=1,Lk
       Hk = Hk_tb(:,:,ik)
       do is=1,Ns
          do js=1,Ns
             gz_imt_dens_constr_slater(is,js) = gz_imt_dens_constr_slater(is,js) + slater(is,js,ik)*wtk(ik)
             do iis=1,Ns
                do jjs=1,Ns
                   Estar = Estar + conjg(gz_imt_Rhop(iis,is))*Hk(iis,jjs)*gz_imt_Rhop(jjs,js)*slater(is,js,ik)*wtk(ik)
                end do
             end do
          end do
       end do
    end do
    !
    gz_imt_energies(1) = Estar+Eloc
    gz_imt_energies(2) = Estar
    gz_imt_energies(3) = Eloc
    !  
  end subroutine gz_imt_measure
  !
  subroutine gz_imt_measure_constr(psi_t,itime)
    complex(8),dimension(nDynamics) :: psi_t
    real(8)                         :: itime
    complex(8),dimension(Ns,Ns,Lk)  :: slater,Hqp
    complex(8),dimension(Ns,Ns)     :: tmpHk,Hk
    complex(8),dimension(Nphi)      :: gzproj
    real(8)                         :: Estar,Eloc,Egz,b0
    real(8),dimension(Ns)           :: vdm_diag,tmp_eHk
    integer                         :: is,js,ik,it,iphi,iis,jjs,itt,ks
    !
    it=imt_t2it(itime,itstep)
    !
    if(size(psi_t).ne.nDynamics) stop "gz_imt_measure_constr/ wrong dimensions"
    call dynamicalVector_2_wfMatrix(psi_t,Hqp,gzproj)
    !
    do is=1,Ns
       do js=1,Ns
          gz_imt_dens_constr_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens(is,js))
          !
          gz_imt_local_density_matrix(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_local_dens(is,js))
       end do
    end do
    !
    gz_imt_unitary_constr = 0.d0
    do iphi=1,Nphi
       gz_imt_unitary_constr = gz_imt_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
    end do
    !
    Estar=0.d0
    !+- first thing first I have to exctract slater(is,js,ik) from the immaginary time integral of the QP hamiltonian
    do ik=1,Lk
       !
       tmpHk= Hqp(:,:,ik)
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
       !       
    end do
    !
    gz_imt_dens_constr_slater=0.d0
    do ik=1,Lk
       Hk = Hk_tb(:,:,ik)
       do is=1,Ns
          do js=1,Ns
             gz_imt_dens_constr_slater(is,js) = gz_imt_dens_constr_slater(is,js) + slater(is,js,ik)*wtk(ik)
          end do
       end do
    end do
    ! !
  end subroutine gz_imt_measure_constr




  !+- SUPERC ROUTINES

  subroutine gz_imt_measure_superc(psi_t,itime,slater_out,gzproj_out,Hqp_out)
    complex(8),dimension(nDynamics) :: psi_t
    real(8)                         :: itime
    complex(8),dimension(2,Ns,Ns,Lk),optional  :: slater_out
    complex(8),dimension(3,Ns,Ns,Lk),optional  :: Hqp_out
    complex(8),dimension(Nphi),optional :: gzproj_out
    complex(8),dimension(2,Ns,Ns,Lk)  :: tmp_slater,slater
    complex(8),dimension(3,Ns,Ns,Lk) :: Hqp 
    complex(8),dimension(2*Ns,2*Ns)     :: tmpHk
    complex(8),dimension(Ns,Ns)     :: Hk
    complex(8),dimension(Nphi)      :: gzproj,gztmp
    real(8)                         :: Estar,Eloc,Egz,nqp
    real(8),dimension(Ns)           :: vdm_diag
    real(8),dimension(2*Ns) :: tmp_eHk

    integer                         :: is,js,ik,iphi,iis,jjs,itt,ks
    !
    if(size(psi_t).ne.nDynamics) stop "gz_imt_measure_/ wrong dimensions"
    call imt_dynamicalVector_2_wfMatrix_superc(psi_t,Hqp,gzproj)
    !
    if(present(gzproj_out)) gzproj_out=gzproj
    if(present(Hqp_out)) Hqp_out=Hqp
    !
    !
    do is=1,Ns
       do js=1,Ns
          gz_imt_local_density_matrix(is,js) = &
               trace_phi_basis_sp(gzproj,phi_spTraces_basis_local_dens(is,js))          
          gz_imt_local_dens_dens(is,js) = &
               trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens_dens(is,js))
          !
          gz_imt_dens_constr_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens(is,js))
          gz_imt_dens_constrA_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens_anomalous(is,js))
          !                                                                                                                                             
          gz_imt_local_sc_order(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_sc_order(is,js))
          !
       end do
       vdm_diag(is) = gz_imt_dens_constr_gzproj(is,is)
    end do
    !
    itt=imt_t2it(itime,0.5d0*itstep,beta_init)
    !
    eLevels=0.d0
    !
    gztmp = get_local_hamiltonian_HLOCphi(gzproj,eLevels)
    Eloc = 0.d0
    do iphi=1,Nphi
       Eloc = Eloc + conjg(gzproj(iphi))*gztmp(iphi)
    end do
    !
    gz_imt_local_angular_momenta(1) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_spin2)
    gz_imt_local_angular_momenta(2) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_spinZ)
    gz_imt_local_angular_momenta(3) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_isoSpin2)
    gz_imt_local_angular_momenta(4) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_isoSpinZ)
    !
    gz_imt_Rhop = hopping_renormalization_normal_sp(gzproj,vdm_diag)
    gz_imt_Qhop = hopping_renormalization_anomalous_sp(gzproj,vdm_diag)
    !
    gz_imt_unitary_constr = 0.d0
    do iphi=1,Nphi
       gz_imt_unitary_constr = gz_imt_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
    end do
    !+- SLATER -+!
    Estar=0.d0
    !+- first thing first I have to exctract slater(is,js,ik) from the immaginary time integral of the QP hamiltonian
    slater=zero
    do ik=1,Lk
       !
       tmpHk(1:Ns,1:Ns) = Hqp(1,:,:,ik) 
       tmpHk(1+Ns:2*Ns,1:Ns) = Hqp(2,:,:,ik) 
       tmpHk(1:Ns,1+Ns:2*Ns) = conjg(transpose(Hqp(2,:,:,ik)))
       tmpHk(1+Ns:2*Ns,1+Ns:2*Ns) = Hqp(3,:,:,ik) 
       !
       call matrix_diagonalize(tmpHk,tmp_eHk)       
       !
       do is=1,Ns
          do js=1,Ns
             do ks=1,Ns
                nqp = fermi(tmp_eHk(ks)-tmp_eHk(ks+Ns),1.d0)                
                slater(1,is,js,ik) = slater(1,is,js,ik) + &
                     conjg(tmpHk(is,ks))*tmpHk(js,ks)*nqp + conjg(tmpHk(is,ks+Ns))*tmpHk(js,ks+Ns)*(1.d0-nqp)
                slater(2,is,js,ik) = slater(2,is,js,ik) + &
                     conjg(tmpHk(is+Ns,ks))*tmpHk(js,ks)*nqp + conjg(tmpHk(is+Ns,ks+Ns))*tmpHk(js,ks+Ns)*(1.d0-nqp)
             end do
             !
          end do
       end do
       !
    end do
    if(present(slater_out)) slater_out=slater
    !
    gz_imt_dens_constr_slater=0.d0
    gz_imt_dens_constrA_slater=0.d0
    do ik=1,Lk
       Hk = Hk_tb(:,:,ik)
       do is=1,Ns
          do js=1,Ns
             gz_imt_dens_constr_slater(is,js) = gz_imt_dens_constr_slater(is,js) + slater(1,is,js,ik)*wtk(ik)
             gz_imt_dens_constrA_slater(is,js) = gz_imt_dens_constrA_slater(is,js) + slater(2,is,js,ik)*wtk(ik)
             do iis=1,Ns
                do jjs=1,Ns
                   !
                   !
                   Estar = Estar + conjg(gz_imt_Rhop(iis,is))*Hk(iis,jjs)*gz_imt_Rhop(jjs,js)*slater(1,is,js,ik)*wtk(ik)
                   Estar = Estar + conjg(gz_imt_Rhop(iis,is))*Hk(iis,jjs)*gz_imt_Qhop(jjs,js)*conjg(slater(2,js,is,ik))*wtk(ik)
                   Estar = Estar + conjg(gz_imt_Qhop(iis,is))*Hk(iis,jjs)*gz_imt_Rhop(jjs,js)*slater(2,is,js,ik)*wtk(ik)
                   Estar = Estar - conjg(gz_imt_Qhop(iis,is))*Hk(iis,jjs)*gz_imt_Qhop(jjs,js)*slater(1,js,is,ik)*wtk(ik)
                   !
                   !
                   if(is.eq.js) then
                      Estar = Estar + conjg(gz_imt_Qhop(iis,is))*Hk(iis,jjs)*gz_imt_Qhop(jjs,js)*wtk(ik)
                   end if
                end do
             end do
          end do
       end do
    end do
    !
    gz_imt_energies(1) = Estar+Eloc
    gz_imt_energies(2) = Estar
    gz_imt_energies(3) = Eloc
    !  
  end subroutine gz_imt_measure_superc



  subroutine gz_imt_measure_constr_superc(psi_t,itime)
    complex(8),dimension(nDynamics) :: psi_t
    real(8)                         :: itime
    complex(8),dimension(2,Ns,Ns,Lk)  :: tmp_slater,slater
    complex(8),dimension(3,Ns,Ns,Lk) :: Hqp 
    complex(8),dimension(2*Ns,2*Ns)     :: tmpHk
    complex(8),dimension(Ns,Ns)     :: Hk
    complex(8),dimension(Nphi)      :: gzproj,gztmp
    real(8)                         :: Estar,Eloc,Egz,nqp
    real(8),dimension(Ns)           :: vdm_diag
    real(8),dimension(2*Ns) :: tmp_eHk

    integer                         :: is,js,ik,iphi,iis,jjs,itt,ks
    !
    if(size(psi_t).ne.nDynamics) stop "gz_imt_measure_/ wrong dimensions"
    call imt_dynamicalVector_2_wfMatrix_superc(psi_t,Hqp,gzproj)
    !
    !
    !
    do is=1,Ns
       do js=1,Ns
          !
          gz_imt_dens_constr_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens(is,js))
          gz_imt_dens_constrA_gzproj(is,js) = trace_phi_basis_sp(gzproj,phi_spTraces_basis_dens_anomalous(is,js))
          !                                                                                                                                             
          gz_imt_local_density_matrix(is,js) = &
               trace_phi_basis_sp(gzproj,phi_spTraces_basis_local_dens(is,js))          
       end do
       vdm_diag(is) = gz_imt_dens_constr_gzproj(is,is)
    end do
    !
    itt=imt_t2it(itime,0.5d0*itstep,beta_init)
    !
    eLevels=0.d0
    !
    gztmp = get_local_hamiltonian_HLOCphi(gzproj,eLevels)
    Eloc = 0.d0
    do iphi=1,Nphi
       Eloc = Eloc + conjg(gzproj(iphi))*gztmp(iphi)
    end do
    !
    !
    gz_imt_unitary_constr = 0.d0
    do iphi=1,Nphi
       gz_imt_unitary_constr = gz_imt_unitary_constr + gzproj(iphi)*conjg(gzproj(iphi))
    end do
    !+- SLATER -+!
    Estar=0.d0
    !+- first thing first I have to exctract slater(is,js,ik) from the immaginary time integral of the QP hamiltonian
    slater=zero
    do ik=1,Lk
       !
       tmpHk(1:Ns,1:Ns) = Hqp(1,:,:,ik) 
       tmpHk(1+Ns:2*Ns,1:Ns) = conjg(transpose(Hqp(2,:,:,ik)))
       tmpHk(1:Ns,1+Ns:2*Ns) = Hqp(2,:,:,ik) 
       tmpHk(1+Ns:2*Ns,1+Ns:2*Ns) = Hqp(3,:,:,ik) 
       call matrix_diagonalize(tmpHk,tmp_eHk)       
       !
       do is=1,Ns
          do js=1,Ns
             do ks=1,Ns
                nqp = fermi(tmp_eHk(ks)-tmp_eHk(ks+Ns),1.d0)                
                slater(1,is,js,ik) = slater(1,is,js,ik) + &
                     conjg(tmpHk(is,ks))*tmpHk(js,ks)*nqp + conjg(tmpHk(is,ks+Ns))*tmpHk(js,ks+Ns)*(1.d0-nqp)
                slater(2,is,js,ik) = slater(2,is,js,ik) + &
                     conjg(tmpHk(is+Ns,ks))*tmpHk(js,ks)*nqp + conjg(tmpHk(is+Ns,ks+Ns))*tmpHk(js,ks+Ns)*(1.d0-nqp)
             end do
             !
          end do
       end do
       !
    end do
    !
    gz_imt_dens_constr_slater=0.d0
    gz_imt_dens_constrA_slater=0.d0
    do ik=1,Lk
       Hk = Hk_tb(:,:,ik)
       do is=1,Ns
          do js=1,Ns
             gz_imt_dens_constr_slater(is,js) = gz_imt_dens_constr_slater(is,js) + slater(1,is,js,ik)*wtk(ik)
             gz_imt_dens_constrA_slater(is,js) = gz_imt_dens_constrA_slater(is,js) + slater(2,is,js,ik)*wtk(ik)
          end do
       end do
    end do
    !  
  end subroutine gz_imt_measure_constr_superc
















  subroutine step_imt_dynamics(nsys,tstep,t,yt,itd_lgrNC,itd_lgr_local_dens,eom_funct,ndens,fix_constr_) 
    integer                                   :: nsys
    real(8)                                   :: tstep,t
    complex(8),dimension(nsys),intent(inout)  :: yt
    complex(8),dimension(nsys)                :: yt_old,yt_new
    complex(8),dimension(Ns,Ns)                ::lgrNC_old,lgrNC_new
    complex(8),dimension(Ns,Ns),intent(inout) :: itd_lgrNC
    real(8),intent(inout)                     :: itd_lgr_local_dens
    real(8),intent(in) :: ndens
    logical,optional :: fix_constr_
    logical :: fix_constr

    real(8),dimension(:),allocatable          :: lgr,delta_out
    complex(8),dimension(:),allocatable       :: lgr_cmplx
    integer                                   :: iter,Nopt
    integer                                   :: i,i0
    real(8)                                   :: delta
    interface
       function eom_funct(t,y,Nsys)
         implicit none
         integer                              :: Nsys
         real(8)                              :: t   
         complex(8),dimension(Nsys)           :: eom_funct
         complex(8),dimension(Nsys)           :: y
       end function eom_funct
    end interface
    !
    fix_constr=.false.
    if(present(fix_constr_)) fix_constr=fix_constr_
    !
    !
    !+- defines the initial time of the step (used for computing the increment integral)
    imt_tstep = t
    !+- number of lagrange params to be optimized -+!
    Nopt = Nvdm_NC_opt + 1
    allocate(lgr_cmplx(Nvdm_NC_opt))
    allocate(lgr(2*Nopt));allocate(delta_out(2*Nopt))
    !
    call vdm_NC_stride_m2v(itd_lgrNC,lgr_cmplx)    
    do i=1,Nvdm_NC_opt
       lgr(i) = dreal(lgr_cmplx(i))
       lgr(i+Nvdm_NC_opt) = dimag(lgr_cmplx(i))
    end do
    lgr(2*Nvdm_NC_opt+1) = itd_lgr_local_dens
    lgr(2*Nvdm_NC_opt+2) = 0.d0    
    !
    yt_old = yt
    !
    if(allocated(imt_lgrNC)) deallocate(imt_lgrNC)
    allocate(imt_lgrNC(Ns,Ns)); 
    imt_lgrNC=zero
    imt_lgr_local_dens=0.d0
    !
    if(fix_constr) then
       call fsolve(imt_fix_constr,lgr,tol=1.d-04,info=iter)    
       ! 
       delta_out = imt_fix_constr(lgr)
       delta=0.d0
       do i=1,2*Nopt
          delta = delta + delta_out(i)**2.d0
       end do
       !
       write(*,*) 'Immaginary Time Dynamics -> lagrange parameters'
       write(*,*) lgr
       write(*,*) delta_out
       write(*,*)
    else
       yt_new = RK4_step(nDynamics,4,tstep,t,yt_old,eom_funct)
    end if
    yt=yt_new
    ! 
    lgr_cmplx=zero
    do i=1,Nvdm_AC_opt
       lgr_cmplx(i) = lgr(i)+xi*lgr(i+Nvdm_NC_opt)
    end do
    call vdm_NC_stride_v2m(lgr_cmplx,itd_lgrNC)
    itd_lgr_local_dens = lgr(2*Nvdm_NC_opt+1)
    !
  contains

    function imt_fix_constr(lgr) result(delta)
      real(8),dimension(:) :: lgr
      real(8),dimension(size(lgr)) :: delta
      complex(8),dimension(:),allocatable  :: lgr_cmplx,delta_cmplx
      complex(8),dimension(Ns,Ns) :: normal_constrGZ,normal_constrSL,delta_constr
      real(8) :: unitary_constrGZ
      integer :: i0,i,is,js
      complex(8) :: tot_dens,tmp_dens
      !
      !
      !+- dump slater_lgr_multipliers
      allocate(lgr_cmplx(Nvdm_NC_opt))
      do i=1,Nvdm_NC_opt
         lgr_cmplx(i) = lgr(i) + xi*lgr(i+Nvdm_NC_opt)
      end do
      call vdm_NC_stride_v2m(lgr_cmplx,imt_lgrNC)   
      imt_lgr_local_dens = lgr(2*Nvdm_NC_opt+1)
      !
      yt_new = RK4_step(nDynamics,4,itstep,t,yt_old,eom_funct)
      !
      call gz_imt_measure_constr(yt_new,t) 
      !
      tot_dens=0.d0
      do is=1,Ns
         do js=1,Ns
            call get_imt_dens_constr_slater(is,js,normal_constrSL(is,js))
            call get_imt_dens_constr_gzproj(is,js,normal_constrGZ(is,js))
            delta_constr(is,js) = normal_constrGZ(is,js)-normal_constrSL(is,js)
         end do
         call get_imt_local_dens(is,is,tmp_dens)
         tot_dens = tot_dens + tmp_dens
      end do
      !
      !
      delta=0.d0
      allocate(delta_cmplx(Nvdm_NC_opt))
      call vdm_NC_stride_m2v(delta_constr,delta_cmplx)
      do i=1,Nvdm_NC_opt
         delta(i) = dreal(delta_cmplx(i))
         delta(i+Nvdm_NC_opt) = dimag(delta_cmplx(i))
      end do
      deallocate(delta_cmplx)
      delta(2*Nvdm_NC_opt+1) = dreal(tot_dens-ndens)
      delta(2*Nvdm_NC_opt+2) = dimag(tot_dens-ndens)
      !      
      if(GZneq_verbose) then
         write(*,*) 'imt-dep LAGRANGE'
         do is=1,Ns
            write(*,'(20F18.10)') imt_lgrNC(is,:)         
         end do
         write(*,'(20F18.10)') imt_lgrU         
         write(*,*)
         write(*,*) 'deviation from constraint conservation'
         write(*,'(20F18.10)') delta
      end if
      !
    end function imt_fix_constr
    !
  end subroutine step_imt_dynamics

  
  !+----------------------------+!
  !+- SUPERCONDUCTING ROUTINES -+!
  !+----------------------------+!
  subroutine step_imt_dynamics_superc_NC_d(nsys,tstep,t,yt,itd_lgrNC,imt_vdm,eom_funct,ndens,fix_constr_) 
    integer                                     :: nsys
    real(8)                                     :: tstep,t
    complex(8),dimension(nsys),intent(inout)    :: yt
    complex(8),dimension(nsys)                  :: yt_old,yt_new
    complex(8),dimension(2,Ns,Ns),intent(inout) :: itd_lgrNC !+- (1,:,:) -> VDM ; (2,:,:) -> gz=sl
    complex(8),dimension(Ns,Ns)                 :: imt_vdm
    real(8),intent(in)                          :: ndens

    real(8),dimension(:),allocatable            :: lgr,delta_out
    complex(8),dimension(:),allocatable         :: lgr_cmplx
    integer                                     :: iter,Nopt
    integer                                     :: i,i0,is
    real(8)                                     :: delta
    logical,optional :: fix_constr_
    logical :: fix_constr
    interface
       function eom_funct(t,y,Nsys)
         implicit none
         integer                                :: Nsys
         real(8)                                :: t   
         complex(8),dimension(Nsys)             :: eom_funct
         complex(8),dimension(Nsys)             :: y
       end function eom_funct
    end interface
    !
    fix_constr=.false.
    if(present(fix_constr_)) fix_constr=fix_constr_
    !+- defines the initial time of the step (used for computing the increment integral)
    imt_tstep = t
    !+- number of lagrange params to be optimized -+!
    Nopt = 2*Nvdm_NC_opt 
    allocate(lgr(Nopt));allocate(delta_out(Nopt))
    !
    allocate(lgr_cmplx(Nvdm_NC_opt))
    i0=0
    call vdm_NC_stride_m2v(itd_lgrNC(1,:,:),lgr_cmplx)    
    do i=1,Nvdm_NC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
    end do
    i0=i0+Nvdm_NC_opt
    call vdm_NC_stride_m2v(itd_lgrNC(2,:,:),lgr_cmplx)    
    do i=1,Nvdm_NC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
    end do
    deallocate(lgr_cmplx)
    !
    yt_old = yt
    !
    if(allocated(imt_lgr_NC)) deallocate(imt_lgr_NC)
    allocate(imt_lgr_NC(2,Ns,Ns));  imt_lgr_NC=zero
    !
    if(allocated(imt_lgr_AC)) deallocate(imt_lgr_AC)
    allocate(imt_lgr_AC(2,Ns,Ns)); imt_lgr_AC=zero
    imt_lgr_local_dens=0.d0
    !
    if(fix_constr) then
       call fsolve(imt_fix_constr,lgr,info=iter)    
       delta_out = imt_fix_constr(lgr)
       delta=0.d0
       do i=1,Nopt
          delta = delta + delta_out(i)**2.d0
       end do
       ! ! !
       write(*,*) 'Immaginary Time Dynamics -> lagrange parameters'
       write(*,*) lgr
       write(*,*) 'Immaginary Time Dynamics -> error'
       write(*,*) delta_out
       write(*,*)
       ! ! ! !       
    else
       yt_new = RK4_step(nDynamics,4,tstep,t,yt_old,eom_funct)
    end if
    !
    yt=yt_new
    ! 
    i0=0
    allocate(lgr_cmplx(Nvdm_NC_opt))    
    do i=1,Nvdm_NC_opt
       lgr_cmplx(i) = lgr(i0+i)
    end do
    call vdm_NC_stride_v2m(lgr_cmplx,itd_lgrNC(1,:,:))
    i0=i0+Nvdm_NC_opt
    do i=1,Nvdm_NC_opt
       lgr_cmplx(i) = lgr(i0+i)
    end do
    call vdm_NC_stride_v2m(lgr_cmplx,itd_lgrNC(2,:,:))
    deallocate(lgr_cmplx)
    !
  contains
    !
    function imt_fix_constr(lgr) result(delta)
      real(8),dimension(:)                :: lgr
      real(8),dimension(size(lgr))        :: delta
      complex(8),dimension(:),allocatable :: lgr_cmplx,delta_cmplx
      complex(8),dimension(Ns,Ns)         :: delta_VDM,delta_constr
      complex(8),dimension(2,Ns,Ns)       :: constrGZ,constrSL,delta_anomalous
      real(8)                             :: unitary_constrGZ
      integer                             :: i0,i,is,js
      complex(8)                          :: tot_dens,tmp_dens
      !
      !
      !+- dump slater_lgr_multipliers
      allocate(lgr_cmplx(Nvdm_NC_opt))
      i0=0
      do i=1,Nvdm_NC_opt
         lgr_cmplx(i) = lgr(i0+i) 
      end do
      call vdm_NC_stride_v2m(lgr_cmplx,imt_lgr_NC(1,:,:))   
      i0=i0+Nvdm_NC_opt
      do i=1,Nvdm_NC_opt
         lgr_cmplx(i) = lgr(i0+i) 
      end do
      call vdm_NC_stride_v2m(lgr_cmplx,imt_lgr_NC(2,:,:))   
      !
      yt_new = RK4_step(nDynamics,4,itstep,t,yt_old,eom_funct)
      !
      call gz_imt_measure_constr_superc(yt_new,t) 
      !
      tot_dens=0.d0
      do is=1,Ns
         do js=1,Ns
            call get_imt_dens_constr_slater(is,js,constrSL(1,is,js))
            call get_imt_dens_constrA_slater(is,js,constrSL(2,is,js))
            call get_imt_dens_constr_gzproj(is,js,constrGZ(1,is,js))
            call get_imt_dens_constrA_gzproj(is,js,constrGZ(2,is,js))
            delta_VDM(is,js) = constrSL(1,is,js) - imt_vdm(is,js)
            delta_constr(is,js) = constrGZ(1,is,js)-imt_vdm(is,js)
            delta_anomalous(1,is,js) = constrSL(2,is,js)
            delta_anomalous(2,is,js) = constrGZ(2,is,js)
         end do
         call get_imt_local_dens(is,is,tmp_dens)
         tot_dens = tot_dens + tmp_dens         
      end do
      !
      !
      delta=0.d0
      i0=0
      allocate(delta_cmplx(Nvdm_NC_opt))     
      call vdm_NC_stride_m2v(delta_VDM,delta_cmplx)
      do i=1,Nvdm_NC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
      end do
      i0=i0+Nvdm_NC_opt
      call vdm_NC_stride_m2v(delta_constr,delta_cmplx)
      do i=1,Nvdm_NC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
      end do
      deallocate(delta_cmplx)
      !      
      if(GZneq_verbose) then
         write(*,*) 'imt-dep LAGRANGE'
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_NC(1,is,:)         
         end do
         write(*,'(20F10.6)')
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_NC(2,is,:)         
         end do
         write(*,'(20F10.6)')
         write(*,*) 'deviation from constraint conservation'
         write(*,'(20F18.10)') delta
         write(*,*) '!+----------------+!'
      end if
      !
    end function imt_fix_constr
    !
  end subroutine step_imt_dynamics_superc_NC_d


  subroutine step_imt_dynamics_superc_NC_dens_d(nsys,tstep,t,yt,itd_lgrNC,itd_lgr_local_dens,imt_vdm,eom_funct,ndens,fix_constr_) 
    integer                                     :: nsys
    real(8)                                     :: tstep,t
    complex(8),dimension(nsys),intent(inout)    :: yt
    complex(8),dimension(nsys)                  :: yt_old,yt_new
    complex(8),dimension(2,Ns,Ns),intent(inout) :: itd_lgrNC !+- (1,:,:) -> VDM ; (2,:,:) -> gz=sl
    real(8),intent(inout)                       :: itd_lgr_local_dens
    complex(8),dimension(Ns,Ns)                 :: imt_vdm
    real(8),intent(in)                          :: ndens
    real(8),dimension(:),allocatable            :: lgr,delta_out
    complex(8),dimension(:),allocatable         :: lgr_cmplx
    integer                                     :: iter,Nopt
    integer                                     :: i,i0,is
    real(8)                                     :: delta
    logical,optional :: fix_constr_
    logical :: fix_constr
    interface
       function eom_funct(t,y,Nsys)
         implicit none
         integer                                :: Nsys
         real(8)                                :: t   
         complex(8),dimension(Nsys)             :: eom_funct
         complex(8),dimension(Nsys)             :: y
       end function eom_funct
    end interface
    !
    fix_constr=.false.
    if(present(fix_constr_)) fix_constr=fix_constr_
    !+- defines the initial time of the step (used for computing the increment integral)
    imt_tstep = t
    !+- number of lagrange params to be optimized -+!
    Nopt = 2*Nvdm_NC_opt + 1
    allocate(lgr(Nopt));allocate(delta_out(Nopt))
    !
    allocate(lgr_cmplx(Nvdm_NC_opt))
    i0=0
    call vdm_NC_stride_m2v(itd_lgrNC(1,:,:),lgr_cmplx)    
    do i=1,Nvdm_NC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
    end do
    i0=i0+Nvdm_NC_opt
    call vdm_NC_stride_m2v(itd_lgrNC(2,:,:),lgr_cmplx)    
    do i=1,Nvdm_NC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
    end do
    deallocate(lgr_cmplx)
    i0=i0+Nvdm_NC_opt
    lgr(i0+1) = itd_lgr_local_dens
    !
    yt_old = yt
    !
    if(allocated(imt_lgr_NC)) deallocate(imt_lgr_NC)
    allocate(imt_lgr_NC(2,Ns,Ns));  imt_lgr_NC=zero
    !
    if(allocated(imt_lgr_AC)) deallocate(imt_lgr_AC)
    allocate(imt_lgr_AC(2,Ns,Ns)); imt_lgr_AC=zero
    imt_lgr_local_dens=0.d0
    !
    if(fix_constr) then
       call fsolve(imt_fix_constr,lgr,info=iter)    
       delta_out = imt_fix_constr(lgr)
       delta=0.d0
       do i=1,Nopt
          delta = delta + delta_out(i)**2.d0
       end do
       ! ! !
       write(*,*) 'Immaginary Time Dynamics -> lagrange parameters'
       write(*,*) lgr
       write(*,*) 'Immaginary Time Dynamics -> error'
       write(*,*) delta_out
       write(*,*)
       ! ! ! !       
    else
       yt_new = RK4_step(nDynamics,4,tstep,t,yt_old,eom_funct)
    end if
    !
    yt=yt_new
    !
    i0=0
    allocate(lgr_cmplx(Nvdm_NC_opt))    
    do i=1,Nvdm_NC_opt
       lgr_cmplx(i) = lgr(i0+i)
    end do
    call vdm_NC_stride_v2m(lgr_cmplx,itd_lgrNC(1,:,:))
    i0=i0+Nvdm_NC_opt
    do i=1,Nvdm_NC_opt
       lgr_cmplx(i) = lgr(i0+i)
    end do
    call vdm_NC_stride_v2m(lgr_cmplx,itd_lgrNC(2,:,:))
    deallocate(lgr_cmplx)
    i0=i0+Nvdm_NC_opt
    itd_lgr_local_dens = lgr(i0+1)
    !
  contains
    !
    function imt_fix_constr(lgr) result(delta)
      real(8),dimension(:)                :: lgr
      real(8),dimension(size(lgr))        :: delta
      complex(8),dimension(:),allocatable :: lgr_cmplx,delta_cmplx
      complex(8),dimension(Ns,Ns)         :: delta_VDM,delta_constr
      complex(8),dimension(2,Ns,Ns)       :: constrGZ,constrSL,delta_anomalous
      real(8)                             :: unitary_constrGZ
      integer                             :: i0,i,is,js
      complex(8)                          :: tot_dens,tmp_dens
      !
      !
      !+- dump slater_lgr_multipliers
      allocate(lgr_cmplx(Nvdm_NC_opt))
      i0=0
      do i=1,Nvdm_NC_opt
         lgr_cmplx(i) = lgr(i0+i) 
      end do
      call vdm_NC_stride_v2m(lgr_cmplx,imt_lgr_NC(1,:,:))   
      i0=i0+Nvdm_NC_opt
      do i=1,Nvdm_NC_opt
         lgr_cmplx(i) = lgr(i0+i) 
      end do
      call vdm_NC_stride_v2m(lgr_cmplx,imt_lgr_NC(2,:,:))   
      i0=i0+Nvdm_NC_opt
      imt_lgr_local_dens = lgr(i0+1)
      !
      yt_new = RK4_step(nDynamics,4,itstep,t,yt_old,eom_funct)
      !
      call gz_imt_measure_constr_superc(yt_new,t) 
      !
      tot_dens=0.d0
      do is=1,Ns
         do js=1,Ns
            call get_imt_dens_constr_slater(is,js,constrSL(1,is,js))
            call get_imt_dens_constrA_slater(is,js,constrSL(2,is,js))
            call get_imt_dens_constr_gzproj(is,js,constrGZ(1,is,js))
            call get_imt_dens_constrA_gzproj(is,js,constrGZ(2,is,js))
            delta_VDM(is,js) = constrSL(1,is,js) - imt_vdm(is,js)
            delta_constr(is,js) = constrGZ(1,is,js)-imt_vdm(is,js)
            delta_anomalous(1,is,js) = constrSL(2,is,js)
            delta_anomalous(2,is,js) = constrGZ(2,is,js)
         end do
         call get_imt_local_dens(is,is,tmp_dens)
         tot_dens = tot_dens + tmp_dens         
      end do
      !
      !
      delta=0.d0
      i0=0
      allocate(delta_cmplx(Nvdm_NC_opt))     
      call vdm_NC_stride_m2v(delta_VDM,delta_cmplx)
      do i=1,Nvdm_NC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
      end do
      i0=i0+Nvdm_NC_opt
      call vdm_NC_stride_m2v(delta_constr,delta_cmplx)
      do i=1,Nvdm_NC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
      end do
      deallocate(delta_cmplx)
      i0=i0+Nvdm_NC_opt
      delta(i0+1) = dreal(tot_dens-ndens)
      !      
      if(GZneq_verbose) then
         write(*,*) 'imt-dep LAGRANGE'
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_NC(1,is,:)         
         end do
         write(*,'(20F10.6)')
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_NC(2,is,:)         
         end do
         write(*,'(20F10.6)')
         write(*,'(20F10.6)')
         write(*,'(20F10.6)') imt_lgr_local_dens
         write(*,*)
         write(*,*) 'deviation from constraint conservation'
         write(*,'(20F18.10)') delta
         write(*,*) '!+----------------+!'
      end if
      !
    end function imt_fix_constr
    !
  end subroutine step_imt_dynamics_superc_NC_dens_d




  subroutine step_imt_dynamics_superc_NC_AC_dens_d(nsys,tstep,t,yt,itd_lgrNC,itd_lgrAC,itd_lgr_local_dens,imt_vdm,eom_funct,ndens,fix_constr_) 
    integer                                     :: nsys
    real(8)                                     :: tstep,t
    complex(8),dimension(nsys),intent(inout)    :: yt
    complex(8),dimension(nsys)                  :: yt_old,yt_new
    complex(8),dimension(2,Ns,Ns),intent(inout) :: itd_lgrNC !+- (1,:,:) -> VDM ; (2,:,:) -> gz=sl
    complex(8),dimension(2,Ns,Ns),intent(inout)   :: itd_lgrAC !+- (1,:,:) -> SL ; (2,:,:) -> GZ
    real(8),intent(inout)                       :: itd_lgr_local_dens
    complex(8),dimension(Ns,Ns)                 :: imt_vdm
    real(8),intent(in)                          :: ndens

    real(8),dimension(:),allocatable            :: lgr,delta_out
    complex(8),dimension(:),allocatable         :: lgr_cmplx
    integer                                     :: iter,Nopt
    integer                                     :: i,i0,is
    real(8)                                     :: delta
    logical,optional :: fix_constr_
    logical :: fix_constr
    interface
       function eom_funct(t,y,Nsys)
         implicit none
         integer                                :: Nsys
         real(8)                                :: t   
         complex(8),dimension(Nsys)             :: eom_funct
         complex(8),dimension(Nsys)             :: y
       end function eom_funct
    end interface
    !
    fix_constr=.false.
    if(present(fix_constr_)) fix_constr=fix_constr_
    !+- defines the initial time of the step (used for computing the increment integral)
    imt_tstep = t
    !+- number of lagrange params to be optimized -+!
    Nopt = 2*Nvdm_NC_opt + 2*Nvdm_AC_opt + 1
    allocate(lgr(Nopt));allocate(delta_out(Nopt))
    !
    allocate(lgr_cmplx(Nvdm_NC_opt))
    i0=0
    call vdm_NC_stride_m2v(itd_lgrNC(1,:,:),lgr_cmplx)    
    do i=1,Nvdm_NC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
    end do
    i0=i0+Nvdm_NC_opt
    call vdm_NC_stride_m2v(itd_lgrNC(2,:,:),lgr_cmplx)    
    do i=1,Nvdm_NC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
    end do
    deallocate(lgr_cmplx)
    i0=i0+Nvdm_NC_opt
    allocate(lgr_cmplx(Nvdm_AC_opt))   
    call vdm_AC_stride_m2v(itd_lgrAC(1,:,:),lgr_cmplx)    
    do i=1,Nvdm_AC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
    end do
    i0=i0+Nvdm_AC_opt
    call vdm_AC_stride_m2v(itd_lgrAC(2,:,:),lgr_cmplx)    
    do i=1,Nvdm_AC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
    end do
    deallocate(lgr_cmplx)
    i0=i0+Nvdm_AC_opt
    lgr(i0+1) = itd_lgr_local_dens
    !
    yt_old = yt
    !
    if(allocated(imt_lgr_NC)) deallocate(imt_lgr_NC)
    allocate(imt_lgr_NC(2,Ns,Ns));  imt_lgr_NC=zero
    !
    if(allocated(imt_lgr_AC)) deallocate(imt_lgr_AC)
    allocate(imt_lgr_AC(2,Ns,Ns)); imt_lgr_AC=zero
    imt_lgr_local_dens=0.d0
    !
    if(fix_constr) then
       call fsolve(imt_fix_constr,lgr,info=iter)    
       delta_out = imt_fix_constr(lgr)
       delta=0.d0
       do i=1,Nopt
          delta = delta + delta_out(i)**2.d0
       end do
       ! ! !
       write(*,*) 'Immaginary Time Dynamics -> lagrange parameters'
       write(*,*) lgr
       write(*,*) 'Immaginary Time Dynamics -> error'
       write(*,*) delta_out
       write(*,*)
       ! ! ! !       
    else
       yt_new = RK4_step(nDynamics,4,tstep,t,yt_old,eom_funct)
    end if
    !
    yt=yt_new
    !
    i0=0
    allocate(lgr_cmplx(Nvdm_NC_opt))    
    do i=1,Nvdm_NC_opt
       lgr_cmplx(i) = lgr(i0+i)
    end do
    call vdm_NC_stride_v2m(lgr_cmplx,itd_lgrNC(1,:,:))
    i0=i0+Nvdm_NC_opt
    do i=1,Nvdm_NC_opt
       lgr_cmplx(i) = lgr(i0+i)
    end do
    call vdm_NC_stride_v2m(lgr_cmplx,itd_lgrNC(2,:,:))
    deallocate(lgr_cmplx)
    i0=i0+Nvdm_NC_opt
    allocate(lgr_cmplx(Nvdm_AC_opt))    
    do i=1,Nvdm_AC_opt
       lgr_cmplx(i) = lgr(i0+i)
    end do
    call vdm_AC_stride_v2m(lgr_cmplx,itd_lgrAC(1,:,:))
    i0=i0+Nvdm_AC_opt
    do i=1,Nvdm_AC_opt
       lgr_cmplx(i) = lgr(i0+i)
    end do
    call vdm_AC_stride_v2m(lgr_cmplx,itd_lgrAC(2,:,:))
    i0=i0+Nvdm_AC_opt
    itd_lgr_local_dens = lgr(i0+1)
    !
  contains
    !
    function imt_fix_constr(lgr) result(delta)
      real(8),dimension(:)                :: lgr
      real(8),dimension(size(lgr))        :: delta
      complex(8),dimension(:),allocatable :: lgr_cmplx,delta_cmplx
      complex(8),dimension(Ns,Ns)         :: delta_VDM,delta_constr
      complex(8),dimension(2,Ns,Ns)       :: constrGZ,constrSL,delta_anomalous
      real(8)                             :: unitary_constrGZ
      integer                             :: i0,i,is,js
      complex(8)                          :: tot_dens,tmp_dens
      !
      !
      !+- dump slater_lgr_multipliers
      allocate(lgr_cmplx(Nvdm_NC_opt))
      i0=0
      do i=1,Nvdm_NC_opt
         lgr_cmplx(i) = lgr(i0+i) !+ xi*lgr(i0+i+Nvdm_NC_opt)
      end do
      call vdm_NC_stride_v2m(lgr_cmplx,imt_lgr_NC(1,:,:))   
      i0=i0+Nvdm_NC_opt
      do i=1,Nvdm_NC_opt
         lgr_cmplx(i) = lgr(i0+i) !+ xi*lgr(i0+i+Nvdm_NC_opt)
      end do
      call vdm_NC_stride_v2m(lgr_cmplx,imt_lgr_NC(2,:,:))   
      deallocate(lgr_cmplx)
      i0=i0+Nvdm_NC_opt
      allocate(lgr_cmplx(Nvdm_AC_opt))
      do i=1,Nvdm_AC_opt
         lgr_cmplx(i) = lgr(i0+i) !+ xi*lgr(i0+i+Nvdm_NC_opt)
      end do
      call vdm_AC_stride_v2m(lgr_cmplx,imt_lgr_AC(1,:,:))   
      i0=i0+Nvdm_AC_opt
      do i=1,Nvdm_AC_opt
         lgr_cmplx(i) = lgr(i0+i) !+ xi*lgr(i0+i+Nvdm_NC_opt)
      end do
      call vdm_AC_stride_v2m(lgr_cmplx,imt_lgr_AC(2,:,:))   
      i0=i0+Nvdm_AC_opt
      deallocate(lgr_cmplx)
      imt_lgr_local_dens = lgr(i0+1)
      !
      !
      yt_new = RK4_step(nDynamics,4,itstep,t,yt_old,eom_funct)
      !
      !
      call gz_imt_measure_constr_superc(yt_new,t) 
      !
      tot_dens=0.d0
      do is=1,Ns
         do js=1,Ns
            call get_imt_dens_constr_slater(is,js,constrSL(1,is,js))
            call get_imt_dens_constrA_slater(is,js,constrSL(2,is,js))
            call get_imt_dens_constr_gzproj(is,js,constrGZ(1,is,js))
            call get_imt_dens_constrA_gzproj(is,js,constrGZ(2,is,js))
            delta_VDM(is,js) = constrSL(1,is,js) - imt_vdm(is,js)
            delta_constr(is,js) = constrGZ(1,is,js)-imt_vdm(is,js)
            delta_anomalous(1,is,js) = constrSL(2,is,js)
            delta_anomalous(2,is,js) = constrGZ(2,is,js)
         end do
         call get_imt_local_dens(is,is,tmp_dens)
         tot_dens = tot_dens + tmp_dens         
      end do
      !
      delta=0.d0
      i0=0
      allocate(delta_cmplx(Nvdm_NC_opt))     
      call vdm_NC_stride_m2v(delta_VDM,delta_cmplx)
      do i=1,Nvdm_NC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
      end do
      i0=i0+Nvdm_NC_opt
      call vdm_NC_stride_m2v(delta_constr,delta_cmplx)
      do i=1,Nvdm_NC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
      end do
      deallocate(delta_cmplx)
      i0=i0+Nvdm_NC_opt
      !
      allocate(delta_cmplx(Nvdm_AC_opt))     
      call vdm_AC_stride_m2v(delta_anomalous(1,:,:),delta_cmplx)
      do i=1,Nvdm_AC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
      end do
      i0 = i0 + Nvdm_AC_opt
      call vdm_AC_stride_m2v(delta_anomalous(2,:,:),delta_cmplx)
      do i=1,Nvdm_AC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
      end do
      i0 = i0 + Nvdm_AC_opt
      deallocate(delta_cmplx)
      delta(i0+1) = dreal(tot_dens-ndens)
      !            
      if(GZneq_verbose) then
         write(*,*) 'imt-dep LAGRANGE'
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_NC(1,is,:)         
         end do
         write(*,'(20F10.6)')
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_NC(2,is,:)         
         end do
         write(*,'(20F10.6)')
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_AC(1,is,:)         
         end do
         write(*,'(20F10.6)')
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_AC(2,is,:)         
         end do
         write(*,'(20F10.6)')
         write(*,'(20F10.6)') imt_lgr_local_dens
         write(*,*)
         write(*,*) 'deviation from constraint conservation'
         write(*,'(20F18.10)') delta
         write(*,*) '!+----------------+!'
      end if
      !
    end function imt_fix_constr
    !
  end subroutine step_imt_dynamics_superc_NC_AC_dens_d




  subroutine step_imt_dynamics_superc_z(nsys,tstep,t,yt,itd_lgrNC,itd_lgrAC,itd_lgr_local_dens,imt_vdm,eom_funct,ndens,fix_constr_) 
    integer                                     :: nsys
    real(8)                                     :: tstep,t
    complex(8),dimension(nsys),intent(inout)    :: yt
    complex(8),dimension(nsys)                  :: yt_old,yt_new
    complex(8),dimension(2,Ns,Ns),intent(inout) :: itd_lgrNC !+- (1,:,:) -> VDM ; (2,:,:) -> gz=sl
    complex(8),dimension(2,Ns,Ns),intent(inout)   :: itd_lgrAC !+- (1,:,:) -> SL ; (2,:,:) -> GZ
    real(8),intent(inout)                       :: itd_lgr_local_dens
    complex(8),dimension(Ns,Ns)                 :: imt_vdm
    real(8),intent(in)                          :: ndens

    real(8),dimension(:),allocatable            :: lgr,delta_out
    complex(8),dimension(:),allocatable         :: lgr_cmplx
    integer                                     :: iter,Nopt
    integer                                     :: i,i0,is
    real(8)                                     :: delta
    logical,optional :: fix_constr_
    logical :: fix_constr
    interface
       function eom_funct(t,y,Nsys)
         implicit none
         integer                                :: Nsys
         real(8)                                :: t   
         complex(8),dimension(Nsys)             :: eom_funct
         complex(8),dimension(Nsys)             :: y
       end function eom_funct
    end interface
    !
    fix_constr=.false.
    if(present(fix_constr_)) fix_constr=fix_constr_
    !+- defines the initial time of the step (used for computing the increment integral)
    imt_tstep = t
    !+- number of lagrange params to be optimized -+!
    Nopt = 2*Nvdm_NC_opt + 2*Nvdm_AC_opt + 1
    Nopt = Nopt*2
    allocate(lgr(Nopt));allocate(delta_out(Nopt))
    !
    allocate(lgr_cmplx(Nvdm_NC_opt))
    i0=0
    call vdm_NC_stride_m2v(itd_lgrNC(1,:,:),lgr_cmplx)    
    do i=1,Nvdm_NC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
       lgr(i0+i+Nvdm_NC_opt) = dimag(lgr_cmplx(i))
    end do
    i0=i0+2*Nvdm_NC_opt
    call vdm_NC_stride_m2v(itd_lgrNC(2,:,:),lgr_cmplx)    
    do i=1,Nvdm_NC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
       lgr(i0+i+Nvdm_NC_opt) = dimag(lgr_cmplx(i))
    end do
    deallocate(lgr_cmplx)
    i0=i0+2*Nvdm_NC_opt
    allocate(lgr_cmplx(Nvdm_AC_opt))   
    call vdm_AC_stride_m2v(itd_lgrAC(1,:,:),lgr_cmplx)    
    do i=1,Nvdm_AC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
       lgr(i0+i+Nvdm_AC_opt) = dimag(lgr_cmplx(i))
    end do
    i0=i0+2*Nvdm_AC_opt
    call vdm_AC_stride_m2v(itd_lgrAC(2,:,:),lgr_cmplx)    
    do i=1,Nvdm_AC_opt
       lgr(i0+i) = dreal(lgr_cmplx(i))
       lgr(i0+i+Nvdm_NC_opt) = dimag(lgr_cmplx(i))
    end do
    deallocate(lgr_cmplx)
    i0=i0+2*Nvdm_AC_opt
    lgr(i0+1) = itd_lgr_local_dens
    lgr(i0+2) = itd_lgr_local_dens
    !
    yt_old = yt
    !
    if(allocated(imt_lgr_NC)) deallocate(imt_lgr_NC)
    allocate(imt_lgr_NC(2,Ns,Ns));  imt_lgr_NC=zero
    !
    if(allocated(imt_lgr_AC)) deallocate(imt_lgr_AC)
    allocate(imt_lgr_AC(2,Ns,Ns)); imt_lgr_AC=zero
    imt_lgr_local_dens=0.d0
    !
    if(fix_constr) then
       call fsolve(imt_fix_constr,lgr,info=iter)    
       delta_out = imt_fix_constr(lgr)
       delta=0.d0
       do i=1,Nopt
          delta = delta + delta_out(i)**2.d0
       end do
       ! ! !
       write(*,*) 'Immaginary Time Dynamics -> lagrange parameters'
       write(*,*) lgr
       write(*,*) 'Immaginary Time Dynamics -> error'
       write(*,*) delta_out
       write(*,*)
       ! ! ! !       
    else
       yt_new = RK4_step(nDynamics,4,tstep,t,yt_old,eom_funct)
    end if
    !
    yt=yt_new
    !
    i0=0
    allocate(lgr_cmplx(Nvdm_NC_opt))    
    do i=1,Nvdm_NC_opt
       lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_NC_opt)
    end do
    call vdm_NC_stride_v2m(lgr_cmplx,itd_lgrNC(1,:,:))
    i0=i0+2*Nvdm_NC_opt
    do i=1,Nvdm_NC_opt
       lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_NC_opt)
    end do
    call vdm_NC_stride_v2m(lgr_cmplx,itd_lgrNC(2,:,:))
    deallocate(lgr_cmplx)
    i0=i0+2*Nvdm_NC_opt
    allocate(lgr_cmplx(Nvdm_AC_opt))    
    do i=1,Nvdm_AC_opt
       lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_AC_opt)
    end do
    call vdm_AC_stride_v2m(lgr_cmplx,itd_lgrAC(1,:,:))
    i0=i0+2*Nvdm_AC_opt
    do i=1,Nvdm_AC_opt
       lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_AC_opt)
    end do
    call vdm_AC_stride_v2m(lgr_cmplx,itd_lgrAC(2,:,:))
    i0=i0+2*Nvdm_AC_opt
    itd_lgr_local_dens = lgr(i0+1)+xi*lgr(i0+2)
    !
  contains
    !
    function imt_fix_constr(lgr) result(delta)
      real(8),dimension(:)                :: lgr
      real(8),dimension(size(lgr))        :: delta
      complex(8),dimension(:),allocatable :: lgr_cmplx,delta_cmplx
      complex(8),dimension(Ns,Ns)         :: delta_VDM,delta_constr
      complex(8),dimension(2,Ns,Ns)       :: constrGZ,constrSL,delta_anomalous
      real(8)                             :: unitary_constrGZ
      integer                             :: i0,i,is,js
      complex(8)                          :: tot_dens,tmp_dens
      !
      i0=0
      allocate(lgr_cmplx(Nvdm_NC_opt))    
      do i=1,Nvdm_NC_opt
         lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_NC_opt)
      end do
      call vdm_NC_stride_v2m(lgr_cmplx,imt_lgr_NC(1,:,:))
      i0=i0+2*Nvdm_NC_opt
      do i=1,Nvdm_NC_opt
         lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_NC_opt)
      end do
      call vdm_NC_stride_v2m(lgr_cmplx,imt_lgr_NC(2,:,:))
      deallocate(lgr_cmplx)
      i0=i0+2*Nvdm_NC_opt
      allocate(lgr_cmplx(Nvdm_AC_opt))    
      do i=1,Nvdm_AC_opt
         lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_AC_opt)
      end do
      call vdm_AC_stride_v2m(lgr_cmplx,imt_lgr_AC(1,:,:))
      i0=i0+2*Nvdm_AC_opt
      do i=1,Nvdm_AC_opt
         lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_AC_opt)
      end do
      call vdm_AC_stride_v2m(lgr_cmplx,imt_lgr_AC(2,:,:))
      i0=i0+2*Nvdm_AC_opt      
      imt_lgr_local_dens = lgr(i0+1)+xi*lgr(i0+2)
      !
      yt_new = RK4_step(nDynamics,4,itstep,t,yt_old,eom_funct)
      !
      call gz_imt_measure_constr_superc(yt_new,t) 
      !
      tot_dens=0.d0
      do is=1,Ns
         do js=1,Ns
            call get_imt_dens_constr_slater(is,js,constrSL(1,is,js))
            call get_imt_dens_constrA_slater(is,js,constrSL(2,is,js))
            call get_imt_dens_constr_gzproj(is,js,constrGZ(1,is,js))
            call get_imt_dens_constrA_gzproj(is,js,constrGZ(2,is,js))
            delta_VDM(is,js) = constrSL(1,is,js) - imt_vdm(is,js)
            delta_constr(is,js) = constrGZ(1,is,js)-imt_vdm(is,js)
            delta_anomalous(1,is,js) = constrSL(2,is,js)
            delta_anomalous(2,is,js) = constrGZ(2,is,js)
         end do
         call get_imt_local_dens(is,is,tmp_dens)
         tot_dens = tot_dens + tmp_dens         
      end do
      !
      delta=0.d0
      i0=0
      allocate(delta_cmplx(Nvdm_NC_opt))     
      call vdm_NC_stride_m2v(delta_VDM,delta_cmplx)
      do i=1,Nvdm_NC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
         delta(i0+i+Nvdm_NC_opt) = dimag(delta_cmplx(i))
      end do
      i0=i0+2*Nvdm_NC_opt
      call vdm_NC_stride_m2v(delta_constr,delta_cmplx)
      do i=1,Nvdm_NC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
         delta(i0+i+Nvdm_NC_opt) = dimag(delta_cmplx(i))
      end do
      deallocate(delta_cmplx)
      i0=i0+2*Nvdm_NC_opt
      !
      allocate(delta_cmplx(Nvdm_AC_opt))     
      call vdm_AC_stride_m2v(delta_anomalous(1,:,:),delta_cmplx)
      do i=1,Nvdm_AC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
         delta(i0+i+Nvdm_NC_opt) = dimag(delta_cmplx(i))
      end do
      i0 = i0 + 2*Nvdm_AC_opt
      call vdm_AC_stride_m2v(delta_anomalous(2,:,:),delta_cmplx)
      do i=1,Nvdm_AC_opt
         delta(i0+i) = dreal(delta_cmplx(i))
         delta(i0+i+Nvdm_NC_opt) = dimag(delta_cmplx(i))
      end do
      i0 = i0 + 2*Nvdm_AC_opt
      deallocate(delta_cmplx)
      delta(i0+1) = dreal(tot_dens-ndens)
      delta(i0+2) = dimag(tot_dens-ndens)
      !      
      if(GZneq_verbose) then
         write(*,*) 'imt-dep LAGRANGE'
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_NC(1,is,:)         
         end do
         write(*,'(20F10.6)')
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_NC(2,is,:)         
         end do
         write(*,'(20F10.6)')
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_AC(1,is,:)         
         end do
         write(*,'(20F10.6)')
         do is=1,Ns
            write(*,'(20F10.6)') imt_lgr_AC(2,is,:)         
         end do
         write(*,'(20F10.6)')
         write(*,'(20F10.6)') imt_lgr_local_dens
         write(*,*)
         write(*,*) 'deviation from constraint conservation'
         write(*,'(20F18.10)') delta
         write(*,*) '!+----------------+!'
      end if
      !
    end function imt_fix_constr
    !
  end subroutine step_imt_dynamics_superc_z
  





  function imt_t2it(time,delta_t,t0_) result(it)
    real(8) :: time,delta_t
    real(8),optional :: t0_
    real(8) :: t0
    integer :: it
    real(8) :: eps=1.d-8
    t0=beta_init
    if(present(t0_)) t0=t0_
    it = (time-t0+eps*delta_t)/delta_t + 1
  end function imt_t2it





END MODULE GZ_IMTDYNAMICS





  ! subroutine step_imt_dynamics(nsys,tstep,t,yt,itd_lgr,eom_funct) 
  !   integer :: nsys
  !   real(8) :: tstep,t
  !   complex(8),dimension(nsys),intent(inout) :: yt
  !   complex(8),dimension(nsys) :: yt_old,yt_new
  !   complex(8),dimension(Ns,Ns),intent(inout) :: itd_lgr
  !   real(8),dimension(:),allocatable            :: lgr,delta_out
  !   complex(8),dimension(:),allocatable :: lgr_cmplx
  !   integer :: iter,Nopt
  !   integer :: i,i0
  !   real(8) :: delta
  !   interface
  !      function eom_funct(t,y,Nsys)
  !        implicit none
  !        integer                    :: Nsys
  !        real(8)                    :: t   
  !        complex(8),dimension(Nsys) :: eom_funct
  !        complex(8),dimension(Nsys) :: y
  !      end function eom_funct
  !   end interface
  !   !+- number of lagrange params to be optimized -+!
  !   Nopt = 2*Nvdm_NC_opt
  !   allocate(lgr_cmplx(Nvdm_NC_opt))
  !   allocate(lgr(2*Nopt+1));allocate(delta_out(2*Nopt+1))
  !   !
  !   call vdm_NC_stride_m2v(itd_lgr(:,:),lgr_cmplx)    
  !   do i=1,Nvdm_NC_opt
  !      lgr(i) = dreal(lgr_cmplx(i))
  !      lgr(i+Nvdm_NC_opt) = dimag(lgr_cmplx(i))
  !   end do    
  !   !
  !   yt_old = yt
  !   !
  !   ! delta_out = fix_anomalous_vdm(lgr)
  !   ! delta=0.d0
  !   ! do i=1,2*Nopt
  !   !    delta = delta + delta_out(i)**2.d0
  !   ! end do   
  !   !
  !   select case(lgr_method)
  !   case('CG_min')
  !      call fmin_cgminimize(lgr,fix_anomalous_vdm_,iter,delta,ftol=1.d-04,itmax=20)    
  !   case('f_zero')
  !      call fsolve(fix_anomalous_vdm,lgr,tol=1.d-04,info=iter)    
  !   case default
  !      call fsolve(fix_anomalous_vdm,lgr,tol=1.d-04,info=iter)    
  !   end select
  !   !
  !   delta_out = fix_anomalous_vdm(lgr)
  !   delta=0.d0
  !   do i=1,2*Nopt
  !      delta = delta + delta_out(i)**2.d0
  !   end do
  !   !
  !   write(*,*) 'Time Dependent lagrange parameters: error'
  !   write(*,*) delta
  !   write(*,*)
  !   !
  !   yt=yt_new
  !   !
  !   lgr_cmplx=zero
  !   do i=1,Nvdm_NC_opt
  !      lgr_cmplx(i) = lgr(i)+xi*lgr(i+Nvdm_NC_opt)
  !   end do
  !   call vdm_NC_stride_v2m(lgr_cmplx,td_lgr(1,:,:))
  !   lgr_cmplx(+2*Nvdm_NC_opt)

  !   i0=2*Nvdm_AC_opt
  !   do i=1,Nvdm_AC_opt
  !      lgr_cmplx(i) = lgr(i0+i)+xi*lgr(i0+i+Nvdm_AC_opt)
  !   end do
  !   call vdm_AC_stride_v2m(lgr_cmplx,td_lgr(2,:,:))  
  !   !    
  ! contains

  !   function fix_anomalous_vdm(lgr) result(delta)
  !     real(8),dimension(:) :: lgr
  !     real(8),dimension(size(lgr)) :: delta
  !     complex(8),dimension(:),allocatable  :: lgr_cmplx,delta_cmplx
  !     complex(8),dimension(Ns,Ns) :: anomalous_constrGZ,anomalous_constrSL
  !     integer :: i0,i,is,js
  !     real(8) :: tmp_test
  !     !
  !     if(allocated(neq_lgr)) deallocate(neq_lgr)
  !     allocate(neq_lgr(2,Ns,Ns)); neq_lgr=zero
  !     !
  !     !+- dump slater_lgr_multipliers
  !     allocate(lgr_cmplx(Nvdm_AC_opt))
  !     do i=1,Nvdm_AC_opt
  !        lgr_cmplx(i) = lgr(i) + xi*lgr(i+Nvdm_AC_opt)
  !     end do
  !     call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(1,:,:))
  !     i0 = 2*Nvdm_AC_opt
  !     !+- dump gzproj_lgr_multipliers
  !     do i=1,Nvdm_AC_opt
  !        lgr_cmplx(i) = lgr(i+i0) + xi*lgr(i+i0+Nvdm_AC_opt)
  !     end do
  !     call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(2,:,:))
  !     !

  !     !+- major things to change are here -+!

      
  !     yt_new = RK_step(nDynamics,4,tstep,t,yt_old,eom_funct)
  !     !
  !     call gz_imt_measure_constr_superc_sp(yt_new,t) 
  !     !
  !     do is=1,Ns
  !        do js=1,Ns
  !           call get_imt_dens_constrA_slater(is,js,anomalous_constrSL(is,js))
  !           call get_imt_dens_constrA_gzproj(is,js,anomalous_constrGZ(is,js))
  !        end do
  !     end do
  !     !
  !     delta=0.d0
  !     allocate(delta_cmplx(Nvdm_AC_opt))
  !     call vdm_AC_stride_m2v(anomalous_constrSL,delta_cmplx)
  !     do i=1,Nvdm_AC_opt
  !        delta(i) = dreal(delta_cmplx(i))
  !        delta(i+Nvdm_AC_opt) = dimag(delta_cmplx(i))
  !     end do
  !     deallocate(delta_cmplx)
  !     i0 = 2*Nvdm_AC_opt
  !     allocate(delta_cmplx(Nvdm_AC_opt))
  !     call vdm_AC_stride_m2v(anomalous_constrGZ,delta_cmplx)
  !     do i=1,Nvdm_AC_opt
  !        delta(i0+i) = dreal(delta_cmplx(i))
  !        delta(i0+i+Nvdm_AC_opt) = dimag(delta_cmplx(i))       
  !     end do
  !     deallocate(delta_cmplx)
  !     tmp_test=0.d0
  !     do i=1,size(lgr)
  !        tmp_test=tmp_test+delta(i)**2.d0
  !     end do
  !     write(*,*) 't-dep LAGRANGE fsolve'
  !     do is=1,Ns
  !        write(*,'(20F18.10)') neq_lgr(1,is,:)         
  !     end do
  !     write(*,*) 
  !     do is=1,Ns
  !        write(*,'(20F18.10)') neq_lgr(2,is,:)         
  !     end do
  !     write(*,*) 'deviation from constraint conservation'
  !     write(*,'(20F18.10)') delta
  !     !
  !   end function fix_anomalous_vdm
  !   !
  !   function fix_anomalous_vdm_(lgr) result(delta)
  !     real(8),dimension(:) :: lgr
  !     real(8) :: delta
  !     complex(8),dimension(:),allocatable  :: lgr_cmplx,delta_cmplx
  !     complex(8),dimension(Ns,Ns) :: anomalous_constrGZ,anomalous_constrSL
  !     integer :: i0,i,is,js
  !     !
  !     if(allocated(neq_lgr)) deallocate(neq_lgr)
  !     allocate(neq_lgr(2,Ns,Ns)); neq_lgr=zero

  !     !+- dump slater_lgr_multipliers
  !     allocate(lgr_cmplx(Nvdm_AC_opt))
  !     do i=1,Nvdm_AC_opt
  !        lgr_cmplx(i) = lgr(i) + xi*lgr(i+Nvdm_AC_opt)
  !     end do
  !     call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(1,:,:))
  !     i0 = 2*Nvdm_AC_opt
  !     !+- dump gzproj_lgr_multipliers
  !     do i=1,Nvdm_AC_opt
  !        lgr_cmplx(i) = lgr(i+i0) + xi*lgr(i+i0+Nvdm_AC_opt)
  !     end do
  !     call vdm_AC_stride_v2m(lgr_cmplx,neq_lgr(2,:,:))
  !     !
  !     !yt_new = RK_step(nDynamics,4,tstep,t,yt_old,gz_equations_of_motion_superc_lgr_sp)
  !     yt_new = RK_step(nDynamics,4,tstep,t,yt_old,eom_funct)
  !     !
  !     call gz_imt_measure_constr_superc_sp(yt_new,t)
  !     !
  !     do is=1,Ns
  !        do js=1,Ns
  !           call get_imt_dens_constrA_slater(is,js,anomalous_constrSL(is,js))
  !           call get_imt_dens_constrA_gzproj(is,js,anomalous_constrGZ(is,js))
  !        end do
  !     end do
  !     !
  !     delta=0.d0
  !     do is=1,Ns
  !        do js=1,Ns
  !           delta = delta + anomalous_constrSL(is,js)*conjg(anomalous_constrSL(is,js))
  !           delta = delta + anomalous_constrGZ(is,js)*conjg(anomalous_constrGZ(is,js))
  !        end do
  !     end do

  !     write(*,*) 'deviation from constraint conservation',delta
  !     !
  !   end function fix_anomalous_vdm_

  ! end subroutine step_imt_dynamics
