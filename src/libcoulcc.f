! ================================================================================================================================ !
module libcoulcc
  !! Provides interfaces and wrappers for calling COULCC that can
  !! be used as they are or as an example for implenting something
  !! similar. COULCC is a program (contained as a submodule) that
  !! calculates complex Coulomb functions using Steed's method.
  !! COULCC can also calculate other fuctions, such as Bessel
  !! and Hankel functions. See the implementation of the procedure
  !! "coulcc_wrapper()" in this file or the implementation of
  !! the procedure "COULCC()" in coulcc_submodule.f for more
  !! details.

  implicit none

  private

  public :: coulf
  public :: coulg
  public :: coulcc
  public :: coulcc_wrapper

  interface
    module subroutine COULCC(XX, ETA1, ZLMIN, NL, FC, GC, FCP, GCP, SIG, MODE1, KFN, IFAIL)
      use iso_fortran_env, only: dp => real64
      integer, intent(in)    :: NL
      integer, intent(in)    :: MODE1
      integer, intent(in)    :: KFN
      integer, intent(inout) :: IFAIL
      complex(dp), intent(in)   :: XX
      complex(dp), intent(in)   :: ETA1
      complex(dp), intent(in)   :: ZLMIN
      complex(dp), intent(out)  :: FC(NL)
      complex(dp), intent(out)  :: GC(NL)
      complex(dp), intent(out)  :: FCP(NL)
      complex(dp), intent(out)  :: GCP(NL)
      complex(dp), intent(out)  :: SIG(NL)
    end subroutine
  end interface

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  impure elemental function coulf(l, eta, x) result(f)
    !! Return the function F_l(η, x)
    use iso_fortran_env, only: dp => real64
    complex(dp), intent(in) :: l
    complex(dp), intent(in) :: x
    complex(dp), intent(in) :: eta
    complex(dp) :: f
    integer, parameter :: nl   = 1
    integer, parameter :: kfn  = 0
    integer, parameter :: mode = 4
    complex(dp) :: zlmin
    complex(dp) :: ceta
    complex(dp) :: cx
    complex(dp) :: fc(1:nl)
    complex(dp) :: gc(1:nl)
    complex(dp) :: fcp(1:nl)
    complex(dp) :: gcp(1:nl)
    complex(dp) :: sig(1:nl)
    ! cxx   = cmplx(xx,  kind = dp)
    ! ceta  = cmplx(eta, kind = dp)
    ! zlmin = cmplx(l,   kind = dp)
    ! call coulcc_wrapper(zlmin, nl, ceta, cxx, fc, fcp, gc, gcp, sig, kfn, mode)
    call coulcc_wrapper(l, nl, eta, x, fc, fcp, gc, gcp, sig, kfn, mode)
    f = fc(nl)
  end function coulf

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  impure elemental function coulg(l, eta, x) result(f)
    !! Return the function F_l(η, x)
    use iso_fortran_env, only: dp => real64
    complex(dp), intent(in) :: l
    complex(dp), intent(in) :: x
    complex(dp), intent(in) :: eta
    complex(dp) :: f
    integer, parameter :: nl   = 1
    integer, parameter :: kfn  = 0
    integer, parameter :: mode = 4
    complex(dp) :: zlmin
    complex(dp) :: ceta
    complex(dp) :: cx
    complex(dp) :: fc(1:nl)
    complex(dp) :: gc(1:nl)
    complex(dp) :: fcp(1:nl)
    complex(dp) :: gcp(1:nl)
    complex(dp) :: sig(1:nl)
    ! cxx   = cmplx(xx,  kind = dp)
    ! ceta  = cmplx(eta, kind = dp)
    ! zlmin = cmplx(l,   kind = dp)
    ! call coulcc_wrapper(zlmin, nl, ceta, cxx, fc, fcp, gc, gcp, sig, kfn, mode)
    call coulcc_wrapper(l, nl, eta, x, fc, fcp, gc, gcp, sig, kfn, mode)
    f = gc(nl)
  end function coulg

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  subroutine coulcc_wrapper(zlmin, nl, eta, xx, f, fp, g, gp, sig, kfn, mode)
    !! A wrapper routine to call COULCC. The integers MODE and KFN
    !! must be set as follows to obtain the desired output from
    !! COULCC :
    !!
    !! |MODE| means the absolute value of the MODE. MODE
    !! and KFN can be used independently. The value of the array
    !! F is always set.
    !! _________________________________
    !!| |MODE|  |  F  |  G  | FP  | GP  |
    !!| 1 11 21 | set | set | set | set |
    !!| 2 12 22 | set | set | --- | --- |
    !!|    3    | set | --- | set | --- |
    !!|    4    | set | --- | --- | --- |
    !! ---------------------------------
    !!
    !! The contents of F, FP, G, and GP (if set by MODE) are
    !! controlled by KFN. If KFN = 0, SIG holds the Coulomb phase
    !! shifts. Otherwise, it does not.
    !! _________________________________________________________
    !!| ARRAY | |MODE|  |                  KFN                  |
    !!|       |         | -1, 0 |   1       |   2       |   3   |
    !!|_________________|_______|___________|___________|_______|
    !!|   F   |   all   | F_l   |  j_l      |  J_l      |  I_l  |
    !!|   G   | 1 2 3 4 | G_l   |  y_l      |  Y_l      |  K_l  |
    !!|   G   |  11 12  | H_l^+ |  h_l^(1)  |  H_l^(1)  |  ---  |
    !!|   G   |  21 22  | H_l^- |  h_l^(2)  |  H_l^(2)  |  ---  |
    !! ---------------------------------------------------------
    !!
    !! If MODE1 < 0, the values returned are scaled by an
    !! exponential ! factor (dependent only on XX) to bring nearer
    !! unity the ! the functions for large |XX|, small ETA, and
    !! |ZL| < |XX|.
    !!   Define SCALE = (  0        if MODE1 > 0
    !!                  (  IMAG(XX) if MODE1 < 0  &  KFN < 3
    !!                  (  REAL(XX) if MODE1 < 0  &  KFN = 3
    !!   then FC = EXP(-ABS(SCALE)) * ( F, j, J, or I)
    !!    and GC = EXP(-ABS(SCALE)) * ( G, y, or Y )
    !!          or EXP(SCALE)       * ( H+, H(1), or K)
    !!          or EXP(-SCALE)      * ( H- or H(2) )
    !!
    use iso_fortran_env, only: dp => real64

    integer,  intent(in)    :: nl
      !! The number of l values for F_l(η ,r)
    integer,  intent(in)    :: mode
      !! Controls which arrays are set on output
    integer,  intent(in)    :: kfn
      !! Controls which functions are represented by the output
      !! arrays

    complex(dp), intent(in) :: eta
      !! The quantity η in F_l(η, r)
    complex(dp), intent(in) :: xx
      !! The quantity r in F_l(η, r), the quantity XX in COULCC
    complex(dp), intent(in) :: zlmin

    complex(dp), intent(inout) :: f(:)
      !! The array containing values of F. Expected bounds 1:nl
    complex(dp), intent(inout) :: g(:)
      !! The array containing values of G. Expected bounds 1:nl
    complex(dp), intent(inout) :: fp(:)
      !! The array containing values of FP. Expected bounds 1:nl
    complex(dp), intent(inout) :: gp(:)
      !! The array containing values of GP. Expected bounds 1:nl
    complex(dp), intent(inout) :: sig(:)
      !! The array containing the Coulomb phase shifts.
      !! Expected bounds 1:nl

    integer :: ifail

    call coulcc(xx, eta, zlmin, nl, f, g, fp, gp, sig, mode, kfn, &
    ifail)

    select case(ifail)
    case(0)
      continue
    case(1:)
      call die("An arithmetic error occurkred during the final recursion")
    case(-1)
      call die("One of the continued fraction calculations failed or there was an arithmetic error")
    case(-2)
      call die("Argument out of range")
    case(-3)
      call die("One or more subroutines corkresponding to lmin could not be calculated. Some values l > lmin may be correct.")
    case(-4)
      call die("Excessive internal cancellation probably renders the result meaningless")
    case default
      call die("unknown error code from COULCC in WCLBES")
    end select

  end subroutine coulcc_wrapper

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  subroutine die(message)
    !! Stop program execution with a message

    use iso_fortran_env, only: stderr => error_unit

    character(*), intent(in), optional :: message

    if( .NOT. present(message)) error stop

    write(stderr,'("STOP",X,"::",X,A)') message

    write(stderr,*)

    error stop

  end subroutine die

! ================================================================================================================================ !
end module libcoulcc
! ================================================================================================================================ !
