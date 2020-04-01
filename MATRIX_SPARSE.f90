!!Some notes here:
!!this moule provides an interface to sparse matrix in ll format
!!the routines perform some generic action on the sparse_matrix object
!!but we still miss some elementary action such as +update_value, +check_value_exist, +delete_value
!!that I am not gonna use in ED code (for which this module is developed).
MODULE MATRIX_SPARSE
  USE SF_CONSTANTS, only:zero
#ifdef _MPI
  USE MPI
#endif
  implicit none
  private

  type sparse_element
     private
     real(8)                               :: val  !value of the entry: double precision
     complex(8)                            :: cval !value of the entry: double complex
     integer                               :: col  !col connected to this compress value
     type(sparse_element),pointer          :: next !link to next entry in the row
  end type sparse_element

  type sparse_row
     private
     integer                               :: size    !size of the list
     type(sparse_element),pointer          :: root    !head/root of the list\== list itself
  end type sparse_row

  type sparse_matrix
     integer                               :: Nrow
     logical                               :: status=.false.
     type(sparse_row),dimension(:),pointer :: row
  end type sparse_matrix

  type sparse_matrix_csr
     logical                          :: status=.false.
     integer                          :: Nrow
     integer                          :: Nnz
     real(8),dimension(:),allocatable :: values
     integer,dimension(:),allocatable :: columns
     integer,dimension(:),allocatable :: rowIndex
  end type sparse_matrix_csr

  type sparse_matrix_csr_z
     logical                          :: status=.false.
     integer                          :: Nrow
     integer                          :: Nnz
     complex(8),dimension(:),allocatable :: values
     integer,dimension(:),allocatable :: columns
     integer,dimension(:),allocatable :: rowIndex
  end type sparse_matrix_csr_z
  !
  public :: sparse_matrix
  public :: sparse_matrix_csr
  public :: sparse_matrix_csr_z
  !


  !INIT SPARSE MATRICES (LL,CSR)
  interface sp_init_matrix
     module procedure sp_init_matrix_ll,sp_init_matrix_csr,sp_init_matrix_csr_z
  end interface sp_init_matrix
  !
  public :: sp_init_matrix      !init the sparse matrix   !checked



  !DELETE SPARSE MATRIX (LL,CSR) OR ONE OF ITS ELEMENTS (LL)
  interface sp_delete_matrix
     module procedure sp_delete_matrix_ll,sp_delete_matrix_csr,sp_delete_matrix_csr_z
  end interface sp_delete_matrix
  !
  public :: sp_delete_matrix    !delete the sparse matrix !checked
  public :: sp_delete_element   !delete n-th/last element !checked



  !GET NUMBER OF NON-ZERO ELEMENTS
  interface sp_get_nnz
     module procedure sp_get_nnz_ll,sp_get_nnz_csr,sp_get_nnz_csr_z
  end interface sp_get_nnz
  !
  public :: sp_get_nnz



  !INSERT ELEMENTS (D,C) IN LL-SPARSE MATRIX
  interface sp_insert_element
     module procedure sp_insert_element_d,sp_insert_element_c
  end interface sp_insert_element
  !
  public :: sp_insert_element   !insert an element        !checked



  !INSERT DIAGONAL ENTRY IN LL-SPARSE MATRIX
  interface sp_insert_diag
     module procedure sp_insert_diag_d,sp_insert_diag_c
  end interface sp_insert_diag
  !
  public :: sp_insert_diag      !insert a vector at diag  !checked



  !GET ELEMENTS ALONG THE DIAGONAL
  interface sp_get_diagonal
     module procedure sp_get_diagonal_d,sp_get_diagonal_c,sp_get_diagonal_csr
  end interface sp_get_diagonal
  !
  public :: sp_get_diagonal     !get diagonal elements    !checked



  !LOAD STANDARD MATRIX INTO SPARSE MATRICES
  interface sp_load_matrix
     module procedure sp_load_matrix_d,sp_load_matrix_c,sp_load_matrix_csr,sp_load_matrix_csr_z
  end interface sp_load_matrix
  !
  public :: sp_load_matrix      !create sparse from array !checked



  !DUMP SPARSE MATRIX INTO STANDARD MATRIX
  interface sp_dump_matrix
     module procedure sp_dump_matrix_d,sp_dump_matrix_c
  end interface sp_dump_matrix
  !
  public :: sp_dump_matrix      !dump sparse into array   !checked



  !PRETTY PRINTING
  interface sp_print_matrix
     module procedure sp_print_matrix_ll,sp_print_matrix_csr,sp_print_matrix_csr_z
  end interface sp_print_matrix
  public :: sp_print_matrix     !print sparse             !checked
  
  !TEST
  public :: sp_test_symmetric

  !HOMEBREW MAT-VEC PRODUCT
  public :: sp_matrix_vector_product_dd !checked
  public :: sp_matrix_vector_product_dc !checked
  public :: sp_matrix_vector_product_cc !checked
  public :: sp_matrix_vector_product_csr!checked
  public :: sp_matrix_vector_product_csr_z ! checked

  interface sp_scalar_matrix_csr
     module procedure sp_scalar_matrix_csr_dd,sp_scalar_matrix_csr_dz,sp_scalar_matrix_csr_zd,sp_scalar_matrix_csr_zz
  end interface sp_scalar_matrix_csr
  public :: sp_scalar_matrix_csr !CHECKED

#ifdef _MPI
  public :: sp_matrix_vector_product_mpi_dd
  public :: sp_matrix_vector_product_mpi_dc
  public :: sp_matrix_vector_product_mpi_cc
#endif


  !GET ELEMENT FROM SPARSE MATRIX
  public :: sp_get_element_d    !pop n-th/last element    !checked
  public :: sp_get_element_c    !""                       !checked
  public :: sp_get_element_csr  !pop n-th/last element    !checked
  !INQUIRE IF ELEMENT EXISTS
  public :: sp_inquire_element  !inquire an element       !checked
  !TRANSFORMATION FROM LL TO CSR SPARSE MATRIX FORM
  public :: sp_copy_ll2csr
  public :: sp_move_ll2csr




contains       




  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  subroutine sp_init_matrix_ll(sparse,N)
    type(sparse_matrix),intent(inout) :: sparse
    integer                           :: i,N
    !put here a delete statement to avoid problems
    if(sparse%status)stop "sp_init_matrix: alreay allocate can not init"
    sparse%Nrow=N
    sparse%status=.true.
    allocate(sparse%row(N))
    do i=1,N
       allocate(sparse%row(i)%root)
       sparse%row(i)%root%next => null()
       sparse%row(i)%size=0
    end do
  end subroutine sp_init_matrix_ll
  !
  subroutine sp_init_matrix_csr(sparse,Nnz,Nrow)
    type(sparse_matrix_csr) :: sparse
    integer                 :: Nnz,Nrow
    if(sparse%status)stop "sp_init_matrix: alreay allocate can not init"
    allocate(sparse%values(Nnz))
    allocate(sparse%columns(Nnz))
    allocate(sparse%rowIndex(Nrow+1))
    sparse%nnz=nnz
    sparse%nrow=nrow
    sparse%values=0.d0
    sparse%columns=0
    sparse%rowIndex=0
    sparse%status=.true.
  end subroutine sp_init_matrix_csr
  !
  subroutine sp_init_matrix_csr_z(sparse,Nnz,Nrow)
    type(sparse_matrix_csr_z) :: sparse
    integer                 :: Nnz,Nrow
    if(sparse%status)stop "sp_init_matrix: alreay allocate can not init"
    allocate(sparse%values(Nnz))
    allocate(sparse%columns(Nnz))
    allocate(sparse%rowIndex(Nrow+1))
    sparse%nnz=nnz
    sparse%nrow=nrow
    sparse%values=zero
    sparse%columns=0
    sparse%rowIndex=0
    sparse%status=.true.
  end subroutine sp_init_matrix_csr_z








  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_matrix_ll(sparse)    
    type(sparse_matrix),intent(inout) :: sparse
    integer                           :: i
    if(.not.sparse%status)stop "Warning SPARSE/sp_delete_matrix: sparse not allocated already."
    do i=1,sparse%Nrow
       call delete_row(sparse%row(i))
       deallocate(sparse%row(i)%root)
    end do
    deallocate(sparse%row)
    sparse%Nrow=0
    sparse%status=.false.
  end subroutine sp_delete_matrix_ll
  !
  subroutine sp_delete_matrix_csr(sparse)    
    type(sparse_matrix_csr),intent(inout) :: sparse
    if(.not.sparse%status)stop "Warning SPARSE/sp_delete_matrix: sparse not allocated already."
    deallocate(sparse%values)
    deallocate(sparse%columns)
    deallocate(sparse%rowIndex)
    sparse%nnz=0
    sparse%nrow=0
    sparse%status=.false.
  end subroutine sp_delete_matrix_csr
  !
  subroutine sp_delete_matrix_csr_z(sparse)    
    type(sparse_matrix_csr_z),intent(inout) :: sparse
    if(.not.sparse%status)stop "Warning SPARSE/sp_delete_matrix: sparse not allocated already."
    deallocate(sparse%values)
    deallocate(sparse%columns)
    deallocate(sparse%rowIndex)
    sparse%nnz=0
    sparse%nrow=0
    sparse%status=.false.
  end subroutine sp_delete_matrix_csr_z
  !+------------------------------------------------------------------+
  !PURPOSE: delete a single element at (i,j) from the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_element(matrix,i,j)
    type(sparse_matrix),intent(inout) :: matrix
    integer,intent(in)                :: i,j
    logical :: delete
    delete = delete_element_from_row(matrix%row(i),col=j)
    if(.not.delete)write(*,"(A,I3,I3)")"sp_delete_element: can not delete element in",i,j
  end subroutine sp_delete_element

  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire row from the sparse matrix (private)
  !+------------------------------------------------------------------+
  subroutine delete_row(row)
    type(sparse_row),intent(inout) :: row
    type(sparse_element),pointer   :: p,c
    do
       p => row%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next=>null()
       deallocate(c)
    end do
  end subroutine delete_row

  !This shoud be better tested!
  !+------------------------------------------------------------------+
  !PURPOSE: delete a given element from a row of the sparse matrix (private)
  !+------------------------------------------------------------------+
  function delete_element_from_row(row,n,col) result(found)
    type(sparse_row),intent(inout)    :: row
    integer,optional                  :: n
    integer,optional                  :: col
    integer                           :: i,pos
    type(sparse_element),pointer      :: p,c
    logical                           :: found
    pos= row%size ; if(present(n))pos=n
    p => row%root
    c => p%next
    found = .false.
    if(present(col))then
       do 
          if(found .OR. .not.associated(c))return
          if(col == c%col)then
             found=.true.
             exit
          else
             p => c
             c => c%next
          endif
       end do
       if(found)then
          p%next => c%next !reallocate skipping the deleted link
          deallocate(c)           !free link
          row%size=row%size-1
       endif
    else
       do i=1,pos 
          if(.not.associated(c))return !empty list
          p => c
          c => c%next
       end do
       found=.true.
       p%next => c%next !reallocate skipping the deleted link
       deallocate(c)           !free link
       row%size=row%size-1
    endif
  end function delete_element_from_row










  !+------------------------------------------------------------------+
  !PURPOSE:  return total number of non-zero elements stored in sparse
  !+------------------------------------------------------------------+
  function sp_get_nnz_ll(sparse) result(Nnz)
    type(sparse_matrix) :: sparse
    integer             :: i
    integer             :: Nnz
    Nnz=0
    do i=1,sparse%Nrow
       Nnz=Nnz+sparse%row(i)%size
    enddo
  end function sp_get_nnz_ll

  function sp_get_nnz_csr(sparse) result(Nnz)
    type(sparse_matrix_csr) :: sparse
    integer                 :: Nnz
    Nnz=size(sparse%values)
  end function sp_get_nnz_csr


  function sp_get_nnz_csr_z(sparse) result(Nnz)
    type(sparse_matrix_csr_z) :: sparse
    integer                 :: Nnz
    Nnz=size(sparse%values)
  end function sp_get_nnz_csr_z










  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element_d(sparse,value,i,j)
    type(sparse_matrix),intent(inout) :: sparse
    real(8),intent(in)                :: value
    integer,intent(in)                :: i,j
    call insert_element_in_row_d(sparse%row(i),value,j)
  end subroutine sp_insert_element_d
  subroutine sp_insert_element_c(sparse,value,i,j)
    type(sparse_matrix),intent(inout) :: sparse
    complex(8),intent(in)             :: value
    integer,intent(in)                :: i,j
    call insert_element_in_row_c(sparse%row(i),value,j)
  end subroutine sp_insert_element_c








  !+------------------------------------------------------------------+
  !PURPOSE: insert a vector of elements at the diagonal of the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_diag_d(sparse,diag)
    type(sparse_matrix),intent(inout)  :: sparse
    real(8),intent(in),dimension(:)    :: diag
    integer                            :: i
    if(size(diag)/=sparse%Nrow)stop "sp_insert_diag: error in dimensions"
    do i=1,size(diag)
       call insert_element_in_row_d(sparse%row(i),diag(i),i)
    enddo
  end subroutine sp_insert_diag_d
  subroutine sp_insert_diag_c(sparse,diag)
    type(sparse_matrix),intent(inout)  :: sparse
    complex(8),intent(in),dimension(:) :: diag
    integer                            :: i
    if(size(diag)/=sparse%Nrow)stop "sp_insert_diag: error in dimensions"
    do i=1,size(diag)
       call insert_element_in_row_c(sparse%row(i),diag(i),i)
    enddo
  end subroutine sp_insert_diag_c








  !+------------------------------------------------------------------+
  !PURPOSE: insert an element in a given row (private) 
  !+------------------------------------------------------------------+
  subroutine insert_element_in_row_d(row,value,column)
    type(sparse_row),intent(inout)    :: row
    real(8) ,intent(in)               :: value
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: p,c
    logical                           :: iadd
    p => row%root
    c => p%next
    iadd = .false.                !check if column already exist    
    do                            !traverse the list
       if(.not.associated(c))exit !empty list or end of the list
       if(c%col == column)then
          iadd=.true.
          exit
       endif
       if(c%col > column)exit
       p => c
       c => c%next
    end do
    if(iadd)then
       c%val=c%val + value
    else
       allocate(p%next)                !Create a new element in the list
       p%next%val = value
       p%next%col = column
       row%size   = row%size+1
       if(.not.associated(c))then !end of the list special case (current=>current%next)
          p%next%next  => null()
       else
          p%next%next  => c      !the %next of the new node come to current
       end if
    endif
  end subroutine insert_element_in_row_d
  !
  subroutine insert_element_in_row_c(row,value,column)
    type(sparse_row),intent(inout)    :: row
    complex(8) ,intent(in)            :: value
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: p,c
    logical :: iadd
    !
    p => row%root
    c => p%next
    iadd = .false.                !check if column already exist
    do                            !traverse the list
       if(.not.associated(c))exit !empty list or end of the list
       if(c%col == column)then
          iadd=.true.
          exit
       endif
       if(column <= c%col)exit
       p => c
       c => c%next
    end do
    if(iadd)then
       c%cval=c%cval + value
    else
       allocate(p%next)           !create a new element in the list
       p%next%cval= value
       p%next%col = column
       row%size   = row%size+1
       if(.not.associated(c))then !end of the list special case (current=>current%next)
          p%next%next  => null()
       else
          p%next%next  => c       !the %next of the new node come to current
       end if
    endif
  end subroutine insert_element_in_row_c













  !+------------------------------------------------------------------+
  !PURPOSE: get the diagonal elements of the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_get_diagonal_d(sparse,diag)
    type(sparse_matrix),intent(inout) :: sparse
    real(8),dimension(:)              :: diag
    integer                           :: Ndim,i
    Ndim=size(diag);if(Ndim/=sparse%Nrow)stop "sp_get_diagonal: error in diag dimension." 
    do i=1,Ndim
       call get_element_from_row_d(sparse%row(i),diag(i),i)
    enddo
  end subroutine  sp_get_diagonal_d

  subroutine sp_get_diagonal_c(sparse,diag)
    type(sparse_matrix),intent(inout) :: sparse
    complex(8),dimension(:)           :: diag
    integer                           :: Ndim,i
    Ndim=size(diag);if(Ndim/=sparse%Nrow)stop "sp_get_diagonal: error in diag dimension." 
    do i=1,Ndim
       call get_element_from_row_c(sparse%row(i),diag(i),i)
    enddo
  end subroutine sp_get_diagonal_c

  subroutine sp_get_diagonal_csr(sparse,diag)
    type(sparse_matrix_csr),intent(inout) :: sparse
    real(8),dimension(:)               :: diag
    integer                               :: Nrow,i
    Nrow=size(diag);if(Nrow/=sparse%Nrow)stop "sp_get_diagonal: error in diag dimension." 
    do i=1,Nrow
       diag(i)=sp_get_element_csr(sparse,i,i)
    enddo
  end subroutine sp_get_diagonal_csr










  !+------------------------------------------------------------------+
  !PURPOSE: get an element from position (i,j) of the sparse matrix
  !+------------------------------------------------------------------+
  function sp_get_element_d(sparse,i,j) result(value)
    type(sparse_matrix),intent(inout) :: sparse    
    integer,intent(in)                :: i,j
    real(8)                           :: value
    call get_element_from_row_d(sparse%row(i),value,j)
  end function sp_get_element_d

  function sp_get_element_c(sparse,i,j) result(value)
    type(sparse_matrix),intent(inout) :: sparse    
    integer,intent(in)                :: i,j
    complex(8)                        :: value
    call get_element_from_row_c(sparse%row(i),value,j)
  end function sp_get_element_c

  function sp_get_element_csr(sparse,i,j) result(value)
    type(sparse_matrix_csr),intent(inout) :: sparse    
    integer,intent(in)                :: i,j
    real(8)                           :: value
    integer :: pos
    value=0.d0
    do pos=sparse%rowIndex(i),sparse%rowIndex(i+1)-1
       if(j==sparse%columns(pos))value=sparse%values(pos)
    enddo
  end function sp_get_element_csr


  function sp_get_element_csr_z(sparse,i,j) result(value)
    type(sparse_matrix_csr_z),intent(inout) :: sparse    
    integer,intent(in)                :: i,j
    complex(8)                           :: value
    integer :: pos
    value=0.d0
    do pos=sparse%rowIndex(i),sparse%rowIndex(i+1)-1
       if(j==sparse%columns(pos))value=sparse%values(pos)
    enddo
  end function sp_get_element_csr_z












  !+------------------------------------------------------------------+
  !PURPOSE: get an element from a given row of the matrix (private)
  !+------------------------------------------------------------------+
  subroutine get_element_from_row_d(row,value,column)
    type(sparse_row),intent(inout)    :: row
    real(8)                           :: value
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: c
    c => row%root%next
    value=0.d0
    do                            !traverse the list
       if(.not.associated(c))return !empty list or end of the list
       if(c%col == column)exit
       c => c%next
    end do
    value = c%val
  end subroutine get_element_from_row_d

  subroutine get_element_from_row_c(row,value,column)
    type(sparse_row),intent(inout)    :: row
    complex(8)                        :: value
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: c
    c => row%root%next
    value=cmplx(0.d0,0.d0,8)
    do                            !traverse the list
       if(.not.associated(c))return !empty list or end of the list
       if(c%col == column)exit
       c => c%next
    end do
    !
    value = c%cval
  end subroutine get_element_from_row_c












  !+------------------------------------------------------------------+
  !PURPOSE: check if a given element exists
  !+------------------------------------------------------------------+
  function sp_inquire_element(sparse,i,j) result(exist)
    type(sparse_matrix),intent(inout) :: sparse    
    integer,intent(in)                :: i,j
    logical                           :: exist
    exist = inquire_element_from_row(sparse%row(i),j)
  end function sp_inquire_element

  !+------------------------------------------------------------------+
  !PURPOSE: check if an element in a given row of the matrix exist (private)
  !+------------------------------------------------------------------+
  function inquire_element_from_row(row,column) result(exist)
    type(sparse_row),intent(inout)    :: row
    logical                           :: exist
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: c
    c => row%root%next
    exist=.false.
    do                            !traverse the list
       if(.not.associated(c))return !empty list or end of the list
       if(c%col == column)exit
       c => c%next
    end do
    exist=.true.
  end function inquire_element_from_row




















  !+------------------------------------------------------------------+
  !PURPOSE: load a regular matrix (2dim array) into a sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_load_matrix_d(matrix,sparse)
    real(8),dimension(:,:),intent(in)  :: matrix
    type(sparse_matrix),intent(inout)  :: sparse    
    integer                            :: i,j,Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    if(Ndim1/=Ndim2)print*,"Warning: SPARSE/load_matrix Ndim1.ne.Ndim2"
    if(sparse%Nrow /= Ndim1) stop "Warning SPARSE/load_matrix: dimensions error"
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=0.d0)call sp_insert_element_d(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine sp_load_matrix_d

  subroutine sp_load_matrix_c(matrix,sparse)
    complex(8),dimension(:,:),intent(in)  :: matrix
    type(sparse_matrix),intent(inout)     :: sparse    
    integer                               :: i,j,Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    if(Ndim1/=Ndim2)print*,"Warning: SPARSE/load_matrix Ndim1.ne.Ndim2"
    if(sparse%Nrow /= Ndim1) stop "Warning SPARSE/load_matrix: dimensions error"
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=cmplx(0.d0,0.d0,8))call sp_insert_element_c(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine sp_load_matrix_c

  subroutine sp_load_matrix_csr(matrix,sparse)
    real(8),dimension(:,:),intent(in)     :: matrix
    type(sparse_matrix_csr),intent(inout) :: sparse
    type(sparse_matrix)                   :: A
    integer                               :: i,j,Ndim1,Ndim2,Nnz
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    if(Ndim1/=Ndim2) stop  "SPARSE/load_matrix Ndim1.ne.Ndim2: modify the code."
    if(sparse%status) stop "SPARSE/load_matrix CSR matrix should not be init on call load. I'll take care of this."
    call sp_init_matrix(A,Ndim1)
    Nnz=0
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=0.d0) then
             Nnz=Nnz+1
             call sp_insert_element_d(A,matrix(i,j),i,j)
          endif
       enddo
    enddo
    call sp_init_matrix(sparse,nnz,Ndim1)
    call sp_move_ll2csr(A,sparse)
  end subroutine sp_load_matrix_csr




  
  subroutine sp_load_matrix_csr_z(matrix,sparse)
    complex(8),dimension(:,:),intent(in)     :: matrix
    type(sparse_matrix_csr_z),intent(inout) :: sparse
    type(sparse_matrix)                   :: A
    integer                               :: i,j,Ndim1,Ndim2,Nnz
    real(8) :: tmp
    complex(8) :: tmp_c
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    if(Ndim1/=Ndim2) stop  "SPARSE/load_matrix Ndim1.ne.Ndim2: modify the code."
    if(sparse%status) stop "SPARSE/load_matrix CSR matrix should not be init on call load. I'll take care of this."
    call sp_init_matrix(A,Ndim1)
    Nnz=0
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=cmplx(0.d0,0.d0,8)) then
             Nnz=Nnz+1
             tmp=dreal(matrix(i,j))
             tmp_c=matrix(i,j)
             call sp_insert_element_c(A,tmp_c,i,j)
             !call sp_insert_element_d(A,tmp,i,j)
          end if
       enddo
    enddo
    !
    call sp_init_matrix(sparse,nnz,Ndim1)
    call sp_move_ll2csr_z(A,sparse)
  end subroutine sp_load_matrix_csr_z






  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  test if a sparse matrix is symmetric 
  !+-----------------------------------------------------------------------------+
  subroutine sp_test_symmetric(sparse,type)
    type(sparse_matrix)                   :: sparse
    logical                               :: is_symmetric
    real(8),dimension(:,:),allocatable    :: rM,rMtra
    complex(8),dimension(:,:),allocatable :: cM,cMher
    integer                               :: Nrow,Ncol
    character(len=1),optional             :: type
    character(len=1)                      :: type_
    type_='d';if(present(type))type_=type
    Nrow=sparse%Nrow
    Ncol=Nrow
    is_symmetric=.false.
    select case(type_)
    case ("d")
       allocate(rM(Nrow,Ncol))
       call sp_dump_matrix(sparse,rM)
       if(abs(maxval(rM-transpose(rM))) < 1.d-12)is_symmetric=.true.
    case("c")
       allocate(cM(Nrow,Ncol))
       call sp_dump_matrix(sparse,cM)
       if( maxval(abs(cM-conjg(transpose(cM))) ) < 1.d-12)is_symmetric=.true.
    end select
    if(is_symmetric)then
       write(*,"(A)")"Matrix IS symmetric"
    else
       write(*,"(A)")"Matrix IS NOT symmetric"
    endif
  end subroutine sp_test_symmetric






  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix_d(sparse,matrix)
    type(sparse_matrix),intent(in)        :: sparse
    real(8),dimension(:,:),intent(inout)  :: matrix
    type(sparse_element),pointer          :: c
    integer                               :: i,Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !if(Ndim1/=Ndim2)print*,"Warning: SPARSE_MATRIX/sp_dump_matrix_d: Ndim1/=Ndim2"
    if(sparse%Nrow /= Ndim1) stop "Warning SPARSE_MATRIX/sp_dump_matrix_d: dimensions error"
    matrix=0.d0
    do i=1,Ndim1
       c => sparse%row(i)%root%next
       do 
          if(.not.associated(c))exit
          matrix(i,c%col) = c%val
          c => c%next  !traverse list
       enddo
    enddo
  end subroutine sp_dump_matrix_d

  subroutine sp_dump_matrix_c(sparse,matrix)
    type(sparse_matrix),intent(in)          :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    type(sparse_element),pointer            :: c
    integer                                 :: i,Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !if(Ndim1/=Ndim2)print*,"Warning: SPARSE_MATRIX/sp_dump_matrix_d: Ndim1/=Ndim2"
    if(sparse%Nrow /= Ndim1)stop "Warning SPARSE/load_matrix: dimensions error"
    matrix=0.d0
    do i=1,Ndim1
       c => sparse%row(i)%root%next
       do 
          if(.not.associated(c))exit
          matrix(i,c%col) = c%cval
          c => c%next  !traverse list
       enddo
    enddo
  end subroutine sp_dump_matrix_c









  !+------------------------------------------------------------------+
  !PURPOSE: copy a sparse LL matrix into a sparse CSR
  !+------------------------------------------------------------------+
  subroutine sp_copy_ll2csr(sparse,M)
    type(sparse_matrix)          :: sparse
    type(sparse_matrix_csr)      :: M
    type(sparse_element),pointer :: c
    integer                      :: i,count,Nrow,Nnz
    if(.not.sparse%status) stop "sp_dump_sparse: sparse not allocated"
    if(.not.M%status) stop "sp_dump_sparse: M not allocated"
    Nnz  = sp_get_nnz(sparse)
    Nrow = sparse%Nrow
    if(Nnz /= size(M%values)) stop "sp_dump_sparse: dimension mismatch"
    count=1
    do i=1,Nrow
       c => sparse%row(i)%root%next
       M%rowIndex(i)=count
       do 
          if(.not.associated(c))exit
          M%values(count)=c%val
          M%columns(count)=c%col
          c => c%next
          count=count+1
       enddo
    enddo
    M%rowIndex(Nrow+1)=Nnz+1
  end subroutine sp_copy_ll2csr
  !








  !+------------------------------------------------------------------+
  !PURPOSE: copy a sparse LL matrix into a sparse CSR
  !+------------------------------------------------------------------+
  subroutine sp_move_ll2csr(sparse,M)
    type(sparse_matrix)          :: sparse
    type(sparse_matrix_csr)      :: M
    type(sparse_element),pointer :: c,p
    integer                      :: i,count,Nrow,Nnz
    if(.not.sparse%status) stop "sp_dump_sparse: sparse not allocated"
    if(.not.M%status) stop "sp_dump_sparse: M not allocated"
    Nnz  = sp_get_nnz(sparse)
    if(Nnz /= size(M%values)) stop "sp_dump_sparse: dimension mismatch"
    Nrow = sparse%Nrow
    count= 1
    do i=1,Nrow
       p => sparse%row(i)%root
       M%rowIndex(i)=count
       do 
          c => p%next
          if(.not.associated(c))exit
          M%values(count)=c%val
          M%columns(count)=c%col
          count=count+1
          p%next => c%next !
          c%next=>null()
          deallocate(c)
       enddo
       nullify(sparse%row(i)%root)
    enddo
    M%rowIndex(Nrow+1)=Nnz+1
    deallocate(sparse%row)
    sparse%Nrow=0
    sparse%status=.false.
  end subroutine sp_move_ll2csr



  subroutine sp_move_ll2csr_z(sparse,M)
    type(sparse_matrix)          :: sparse
    type(sparse_matrix_csr_z)      :: M
    type(sparse_element),pointer :: c,p
    integer                      :: i,count,Nrow,Nnz
    if(.not.sparse%status) stop "sp_dump_sparse: sparse not allocated"
    if(.not.M%status) stop "sp_dump_sparse: M not allocated"
    Nnz  = sp_get_nnz(sparse)
    if(Nnz /= size(M%values)) stop "sp_dump_sparse: dimension mismatch"
    Nrow = sparse%Nrow
    count= 1
    do i=1,Nrow
       p => sparse%row(i)%root
       M%rowIndex(i)=count
       do 
          c => p%next
          if(.not.associated(c))exit
          M%values(count)=c%cval
          M%columns(count)=c%col
          count=count+1
          p%next => c%next !
          c%next=>null()
          deallocate(c)
       enddo
       nullify(sparse%row(i)%root)
    enddo
    M%rowIndex(Nrow+1)=Nnz+1
    deallocate(sparse%row)
    sparse%Nrow=0
    sparse%status=.false.
  end subroutine sp_move_ll2csr_z













  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+
  subroutine sp_print_matrix_ll(sparse,unit,fmt,type,full)
    type(sparse_matrix)            :: sparse
    integer,optional               :: unit
    integer                        :: i,j,unit_,Ns
    character(len=*),optional      :: fmt
    character(len=64)              :: fmt_
    character(len=1),optional      :: type
    character(len=1)               :: type_
    logical,optional               :: full
    logical                        :: full_
    unit_=6;if(present(unit))unit_=unit
    fmt_='F6.2';if(present(fmt))fmt_=fmt
    type_='d';if(present(type))type_=type
    full_=.false.;if(present(full))full_=full
    Ns=sparse%Nrow
    select case(type_)
    case('d')
       if(full_)then
          write(*,*)"Print sparse matrix (full mode < 100) ->",unit_
          do i=1,Ns
             write(unit_,"(100"//trim(fmt_)//",1X)")(sp_get_element_d(sparse,i,j),j=1,Ns)
          enddo
       else
          write(*,*)"Print sparse matrix (compact mode) ->",unit_
          do i=1,Ns
             call print_row_d(sparse%row(i),unit_,fmt_)
          enddo
       endif
    case('c')
       if(full_)then
          write(*,*)"Print sparse matrix (full mode < 100) ->",unit_
          do i=1,Ns
             write(unit_,"(100("//trim(fmt_)//",A1,"//trim(fmt_)//",2X))")(&
                  real(sp_get_element_c(sparse,i,j)),",",imag(sp_get_element_c(sparse,i,j)),j=1,Ns)
          enddo
       else
          write(*,*)"Print sparse matrix (compact mode) ->",unit_
          do i=1,Ns
             call print_row_c(sparse%row(i),unit_,fmt_)
          enddo
       endif
    end select
    write(unit_,*)
  end subroutine sp_print_matrix_ll

  subroutine sp_print_matrix_csr(sparse,unit,fmt,full)
    type(sparse_matrix_csr)        :: sparse
    integer,optional               :: unit
    integer                        :: i,j,unit_,Ns
    character(len=*),optional      :: fmt
    character(len=64)              :: fmt_
    logical,optional               :: full
    logical                        :: full_
    unit_=6;if(present(unit))unit_=unit
    fmt_='F8.3';if(present(fmt))fmt_=fmt
    full_=.false.;if(present(full))full_=full
    Ns=sparse%Nrow
    write(*,*)"Print sparse matrix (full mode < 100) ->",unit_
    do i=1,Ns
       write(*,"(100"//trim(fmt_)//",1X)")(sp_get_element_csr(sparse,i,j),j=1,Ns)
    enddo
    write(unit_,*)
  end subroutine sp_print_matrix_csr

  subroutine sp_print_matrix_csr_z(sparse,unit,fmt,full)
    type(sparse_matrix_csr_z)        :: sparse
    integer,optional               :: unit
    integer                        :: i,j,unit_,Ns
    character(len=*),optional      :: fmt
    character(len=64)              :: fmt_
    logical,optional               :: full
    logical                        :: full_
    unit_=6;if(present(unit))unit_=unit
    fmt_='F8.3';if(present(fmt))fmt_=fmt
    full_=.false.;if(present(full))full_=full
    Ns=sparse%Nrow
    write(*,*)"Print sparse matrix (full mode < 100) ->",unit_
    do i=1,Ns
       write(*,"(100"//trim(fmt_)//",1X)")(sp_get_element_csr_z(sparse,i,j),j=1,Ns)
    enddo
    write(unit_,*)
  end subroutine sp_print_matrix_csr_z








  !+------------------------------------------------------------------+
  !PURPOSE: print an entire row of the sparse matrix (private)
  !+------------------------------------------------------------------+
  subroutine print_row_d(row,unit,fmt)
    type(sparse_row),intent(in)  :: row
    type(sparse_element),pointer :: c
    integer                      :: count=0
    integer                      :: unit
    character(len=*)             :: fmt
    c => row%root%next   !assume is associated,ie list exists
    do
       if(.not.associated(c))exit
       count=count+1       
       write(unit,"("//trim(fmt)//",A1,I3,1X)",advance='no')c%val,',',c%col
       c => c%next  !traverse list
    end do
    write(unit,*)
  end subroutine print_row_d

  subroutine print_row_c(row,unit,fmt)
    type(sparse_row),intent(in)   :: row
    type(sparse_element),pointer  :: c
    integer                       :: count=0
    integer,optional :: unit
    integer          :: unit_
    character(len=*),optional :: fmt
    character(len=64)         :: fmt_
    unit_=6;if(present(unit))unit_=unit
    fmt_='F15.9';if(present(fmt))fmt_=fmt
    c => row%root%next   !assume is associated,ie list exists
    do
       if(.not.associated(c))exit
       count=count+1
       write(unit_,"(2"//trim(fmt_)//",A1,I3,3X)",advance='no')c%cval,',',c%col
       c => c%next  !traverse list
    end do
    write(unit_,*)
  end subroutine print_row_c












  !+------------------------------------------------------------------+
  !PURPOSE: given a vector vin, perform the matrix-vector multiplication
  ! H_sparse * vin and put the result in vout.
  !+------------------------------------------------------------------+
  subroutine sp_matrix_vector_product_dd(sparse,Ndim,vin,vout)
    integer                               :: Ndim
    type(sparse_matrix),intent(in)        :: sparse
    real(8),dimension(Ndim),intent(in)    :: vin
    real(8),dimension(Ndim),intent(inout) :: vout
    type(sparse_element),pointer          :: c
    integer                               :: i
    vout=0.d0
    !$omp parallel do shared (Ndim,sparse,vout,vin) private (i,c) schedule(static,1) if(Ndim>5000)
    do i=1,Ndim
       c => sparse%row(i)%root%next       
       matmul: do  
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next  !traverse list
       end do matmul
    end do
    !$omp end parallel do
  end subroutine sp_matrix_vector_product_dd
  !+------------------------------------------------------------------+
  subroutine sp_matrix_vector_product_dc(sparse,Ndim,vin,vout)
    integer                                  :: Ndim
    type(sparse_matrix),intent(in)           :: sparse
    complex(8),dimension(Ndim),intent(in)    :: vin
    complex(8),dimension(Ndim),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    vout=cmplx(0.d0,0.d0,8)
    !$omp parallel do shared (Ndim,sparse,vout,vin) private (i,c) schedule(static,1) if(Ndim>5000)
    do i=1,Ndim
       c => sparse%row(i)%root%next       
       matmul: do
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next  !traverse list
       end do matmul
    end do
    !$omp end parallel do
  end subroutine sp_matrix_vector_product_dc
  !+------------------------------------------------------------------+


  
  function sp_matrix_vector_product_csr(Nrow,sparse,vin) result(vout)
    integer                            :: Nrow
    real(8),dimension(Nrow)            :: vout
    type(sparse_matrix_csr),intent(in) :: sparse
    real(8),dimension(Nrow),intent(in) :: vin
    integer                            :: i,pos
    vout=0.d0
    do i=1,Nrow
       do pos=sparse%rowIndex(i),sparse%rowIndex(i+1)-1
          vout(i) = vout(i) + sparse%values(pos)*vin(sparse%columns(pos))
       end do
    end do
  end function sp_matrix_vector_product_csr
  !
  function sp_matrix_vector_product_csr_z(Nrow,sparse,vin) result(vout)
    integer                               :: Nrow
    complex(8),dimension(Nrow)            :: vout
    type(sparse_matrix_csr_z),intent(in)  :: sparse
    complex(8),dimension(Nrow),intent(in) :: vin
    integer                               :: i,pos
    vout=zero
    do i=1,Nrow
       do pos=sparse%rowIndex(i),sparse%rowIndex(i+1)-1
          vout(i) = vout(i) + sparse%values(pos)*vin(sparse%columns(pos))
       end do
    end do
  end function sp_matrix_vector_product_csr_z
  !
  !
  function sp_scalar_matrix_csr_dd(sparse_in,x) result(sparse_out)
    real(8) :: x
    type(sparse_matrix_csr) :: sparse_in,sparse_out
    integer :: i,nnz_,nrow_
    !
    nnz_=sparse_in%nnz
    nrow_=sparse_in%nrow
    call sp_init_matrix(sparse_out,nnz_,nrow_)
    !
    do i=1,nnz_
       sparse_out%values(i) = x*sparse_in%values(i)
       sparse_out%columns(i) = sparse_in%columns(i)
    end do
    do i=1,nrow_+1
       sparse_out%rowIndex(i) = sparse_in%rowIndex(i)
    end do
  end function sp_scalar_matrix_csr_dd
  !
  function sp_scalar_matrix_csr_dz(sparse_in,x) result(sparse_out)
    real(8) :: x
    type(sparse_matrix_csr_z) :: sparse_in,sparse_out
    integer :: i,nnz_,nrow_
    !
    nnz_=sparse_in%nnz
    nrow_=sparse_in%nrow
    call sp_init_matrix(sparse_out,nnz_,nrow_)
    !
    do i=1,nnz_
       sparse_out%values(i) = x*sparse_in%values(i)
       sparse_out%columns(i) = sparse_in%columns(i)
    end do
    do i=1,nrow_+1
       sparse_out%rowIndex(i) = sparse_in%rowIndex(i)
    end do
  end function sp_scalar_matrix_csr_dz
  !
  function sp_scalar_matrix_csr_zd(sparse_in,x) result(sparse_out)
    complex(8)              :: x
    type(sparse_matrix_csr) :: sparse_in
    type(sparse_matrix_csr_z) :: sparse_out
    integer :: i,nnz_,nrow_
    !
    nnz_=sparse_in%nnz
    nrow_=sparse_in%nrow
    call sp_init_matrix(sparse_out,nnz_,nrow_)
    !
    do i=1,nnz_
       sparse_out%values(i) = x*sparse_in%values(i)
       sparse_out%columns(i) = sparse_in%columns(i)
    end do
    do i=1,nrow_+1
       sparse_out%rowIndex(i) = sparse_in%rowIndex(i)
    end do
  end function sp_scalar_matrix_csr_zd
  !
  function sp_scalar_matrix_csr_zz(sparse_in,x) result(sparse_out)
    complex(8) :: x
    type(sparse_matrix_csr_z),intent(in) :: sparse_in
    type(sparse_matrix_csr_z) ::    sparse_out
    integer :: i,nnz_,nrow_
    !
    nnz_=sparse_in%nnz
    nrow_=sparse_in%nrow
    call sp_init_matrix(sparse_out,nnz_,nrow_)
    !
    do i=1,nnz_
       sparse_out%values(i) = x*sparse_in%values(i)
       sparse_out%columns(i) = sparse_in%columns(i)
    end do
    do i=1,nrow_+1
       sparse_out%rowIndex(i) = sparse_in%rowIndex(i)
    end do
  end function sp_scalar_matrix_csr_zz



  !

  subroutine sp_matrix_vector_product_cc(sparse,Ndim,vin,vout)
    integer                                  :: Ndim
    type(sparse_matrix),intent(in)           :: sparse
    complex(8),dimension(Ndim),intent(in)    :: vin
    complex(8),dimension(Ndim),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    vout=cmplx(0.d0,0.d0,8)
    !$omp parallel do shared (Ndim,sparse,vout,vin) private (i,c) schedule(static,1) if(Ndim>5000)
    do i=1,Ndim
       c => sparse%row(i)%root%next
       matmul: do  
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%cval*vin(c%col)
          c => c%next  !traverse list
       end do matmul
    end do
    !$omp end parallel do
  end subroutine sp_matrix_vector_product_cc



#ifdef _MPI
  subroutine sp_matrix_vector_product_mpi_dd(sparse,Ndim,vin,Nloc,vout)
    integer                               :: Ndim,Nloc
    type(sparse_matrix),intent(in)        :: sparse
    real(8),dimension(Ndim),intent(in)    :: vin
    real(8),dimension(Nloc),intent(inout) :: vout
    type(sparse_element),pointer          :: c
    integer                               :: i
    integer                               :: Nini,Nfin
    vout=0.d0
    do i=1,Nloc
       c => sparse%row(i)%root%next
       matmul: do
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next
       end do matmul
    end do
  end subroutine sp_matrix_vector_product_mpi_dd
  !+------------------------------------------------------------------+
  subroutine sp_matrix_vector_product_mpi_dc(sparse,Ndim,vin,Nloc,vout)
    integer                                  :: Ndim,Nloc
    type(sparse_matrix),intent(in)           :: sparse
    complex(8),dimension(Ndim),intent(in)    :: vin
    complex(8),dimension(Nloc),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    integer                                  :: Nini,Nfin
    vout=zero
    do i=1,Nloc
       c => sparse%row(i)%root%next       
       matmul: do
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next
       end do matmul
    end do
  end subroutine sp_matrix_vector_product_mpi_dc
  !+------------------------------------------------------------------+
  subroutine sp_matrix_vector_product_mpi_cc(sparse,Ndim,vin,Nloc,vout)
    integer                                  :: Ndim,Nloc
    type(sparse_matrix),intent(in)           :: sparse
    complex(8),dimension(Ndim),intent(in)    :: vin
    complex(8),dimension(Nloc),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    integer                                  :: Nini,Nfin
    vout=zero
    do i=1,Nloc
       c => sparse%row(i)%root%next       
       matmul: do
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%cval*vin(c%col)
          c => c%next
       end do matmul
    end do
  end subroutine sp_matrix_vector_product_mpi_cc
#endif

end module MATRIX_SPARSE
