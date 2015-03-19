module norm
  use, intrinsic :: iso_c_binding
  implicit none
  
  
  integer, public, parameter :: norm_infinity   = 1
  integer, public, parameter :: norm_one        = 2
  integer, public, parameter :: norm_frobenius  = 3
  
  integer, public, parameter :: norm_error_norm = -1
  integer, public, parameter :: norm_error_m    = -2
  integer, public, parameter :: norm_error_n    = -3
  integer, public, parameter :: norm_error_oom  = -2147483647
  contains
  
  function matnorm(norm, m, n, x, ret) &
  bind(c, name='matnorm')
    implicit none
    ! In/Out
    integer(kind=c_int), intent(in) :: norm
    integer(kind=c_int), intent(in) :: m, n
    real(kind=c_double), intent(in) :: x(m, n)
    real(kind=c_double), intent(out) :: ret
    integer(kind=c_int) :: matnorm
    ! Local
    character(len=1) :: normchar
    integer :: lwork, allocerr
    double precision :: ldx
    double precision, allocatable :: work(:)
    ! Functions
    double precision :: dlange
    
    
    ret = 0.0d0
    lwork = 0
    matnorm = 0
    
    if (m < 1) then
      matnorm = norm_error_m
      return
    else if (n < 1) then
      matnorm = norm_error_n
      return
    end if
    
    if (norm == norm_infinity) then
      normchar = "I"
    else if (norm == norm_one) then
      normchar = "O"
    else if (norm == norm_frobenius) then
      normchar = "F"
    else 
      matnorm = norm_error_norm
      return
    end if
    
    if (normchar == "I") then
      lwork = m
      allocerr = 0
      allocate(work(lwork), stat=allocerr)
      if (allocerr /= 0) then
        matnorm = norm_error_oom
        return
      end if
    end if
    
    ldx = max(1, m)
    
    ret = dlange(normchar, m, n, x, ldx, work)
    
    if (allocated(work)) deallocate(work)
    
    return
  end
  
  
end module

