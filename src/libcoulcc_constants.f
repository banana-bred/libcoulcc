! ================================================================================================================================ !
module libcoulcc_constants
  !! Contains constants used throughout coulcc

  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: ZERO
  public :: QUART
  public :: HALF
  public :: ONE
  public :: TWO
  public :: FOUR
  public :: CI

  real(dp),    parameter :: ZERO  = 0.0_dp
  real(dp),    parameter :: QUART = 0.25_dp
  real(dp),    parameter :: HALF  = 0.5_dp
  real(dp),    parameter :: ONE   = 1.0_dp
  real(dp),    parameter :: TWO   = 2.0_dp
  real(dp),    parameter :: FOUR  = 4.0_dp
  complex(dp), parameter :: CI = (zero, one)

! ================================================================================================================================ !
end module libcoulcc_constants
! ================================================================================================================================ !
