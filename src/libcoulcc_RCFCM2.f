! ================================================================================================================================ !
module libcoulcc__rcfcm2
  !!  Replaces the COMMON /RCFCM2/ in the original COULCC

  use iso_fortran_env, only: dp => real64

  implicit none

  save

  private

  public :: EVEN
  public :: M
  public :: M2M1
  public :: MP12
  public :: X1

  logical :: EVEN
  integer :: M
  integer :: M2M1
  integer :: MP12
  complex(dp) :: X1

! ================================================================================================================================ !
end module libcoulcc__rcfcm2
! ================================================================================================================================ !
