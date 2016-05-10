program GUTZ_mb
  USE SCIFOR
  !
  USE DMFT_MISC
  USE SF_PARSE_INPUT
  !
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_LOCAL_HAMILTONIAN
  USE GZ_OPTIMIZED_ENERGY
  USE GZ_ENERGY_MINIMIZATION
  USE GZ_EFFECTIVE_HOPPINGS
  !
  USE GZ_MATRIX_BASIS
  !
  implicit none
  character(len=200) :: store_dir,read_dir
  !
  !+- hamiltonian details -+!
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")
  !
  call initialize_local_fock_space
  !
  store_dir='/mnt/home/giacomo.mazza/project_data/OFFICINA/test_gz_mb/TEST_STORE_PHI_TRACES/READ_PHI/'
  read_dir='/homepmc/giacomo.mazza/etc_local/GZ_basis/'
  !
  store_dir='/homepmc/giacomo.mazza/etc_local/GZ_basis/'
  ! read_dir='./STORE_DIR/'
  !
  call init_variational_matrices(wf_symmetry,store_dir_=store_dir)
  !call init_variational_matrices(wf_symmetry,store_dir_=store_dir,read_dir_=read_dir)  
  !call get_local_hamiltonian_trace
  !
  
  !
end program GUTZ_mb



!AMOEBA TEST


