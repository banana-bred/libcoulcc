! ================================================================================================================================ !
module libcoulcc_steed
  !!  Replaces the COMMON /STEED/ in the original COULCC
  !!   common blocks are for information & storage only.
  !!   (they are not essential to working of the code)

  use iso_fortran_env, only: dp => real64

  implicit none

  save

  private

  public :: NFP
  public :: N11
  public :: NPQ
  public :: N20
  public :: KAS
  public :: RERR

  integer :: NFP
  integer :: N11
  integer :: NPQ(2)
  integer :: N20
  integer :: KAS(2)
  real(dp) :: RERR

! ================================================================================================================================ !
end module libcoulcc_steed
! ================================================================================================================================ !
