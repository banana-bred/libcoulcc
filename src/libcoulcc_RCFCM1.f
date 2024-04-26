! ================================================================================================================================ !
module libcoulcc_rcfcm1
  !!  Replaces the COMMON /RCFCM1/ in the original COULCC

  use iso_fortran_env, only: dp => real64

  implicit none

  save

  private

  public :: PK
  public :: EK
  public :: CLGAA
  public :: CLGAB
  public :: CLGBB
  public :: DSIG
  public :: TPK1
  public :: W
  public :: RL
  public :: FCL1
  public :: Q
  public :: GAM
  public :: HCL
  public :: HPL
  public :: FCM
  public :: HCL1
  public :: ALPHA
  public :: BETA
  public :: PL


  complex(dp) :: PK
  complex(dp) :: EK
  complex(dp) :: CLGAA
  complex(dp) :: CLGAB
  complex(dp) :: CLGBB
  complex(dp) :: DSIG
  complex(dp) :: TPK1
  complex(dp) :: W
  complex(dp) :: RL
  complex(dp) :: FCL1
  complex(dp) :: Q
  complex(dp) :: GAM
  complex(dp) :: HCL
  complex(dp) :: HPL
  complex(dp) :: FCM
  complex(dp) :: HCL1
  complex(dp) :: ALPHA
  complex(dp) :: BETA
  complex(dp) :: PL

! ================================================================================================================================ !
end module libcoulcc_rcfcm1
! ================================================================================================================================ !
