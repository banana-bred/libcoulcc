module libcoulcc__functions
  use iso_fortran_env, only: dp => real64, qp => real128

  implicit none

  private

  public :: log_gamma

  interface log_gamma
    module procedure :: l_gamma_cdp
  end interface log_gamma

contains
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function l_gamma_cdp(z) result (res)
  !
  ! log_gamma function for any complex number, excluding negative whole number
  ! "Computation of special functions", Shanjie Zhang & Jianmin Jin, 1996, p.48
  ! "Computing the principal branch of log-gamma", D.E.G. Hare,
  ! J. of Algorithms, 25(2), 1997 p. 221–236
  !
  ! Fortran 90 program by Jim-215-Fisher
  !
      complex(dp), intent(in) :: z
      complex(dp) :: res, z1, z2
      real(dp) :: d
      integer :: m, i
      complex(qp) :: zr, zr2, sum, s
      real(dp), parameter :: z_limit = 10.0_dp, zero_k1 = 0.0_dp
      integer, parameter :: n = 20
      real(qp), parameter :: zero = 0.0_qp, one = 1.0_qp,              &
                           pi = acos(-one), ln2pi = log(2 * pi)
      real(qp), parameter :: a(n) = [                                          &
                         .8333333333333333333333333333333333333333E-1_qp,&
                        -.2777777777777777777777777777777777777778E-2_qp,&
                         .7936507936507936507936507936507936507937E-3_qp,&
                        -.5952380952380952380952380952380952380952E-3_qp,&
                         .8417508417508417508417508417508417508418E-3_qp,&
                        -.1917526917526917526917526917526917526918E-2_qp,&
                         .6410256410256410256410256410256410256410E-2_qp,&
                        -.2955065359477124183006535947712418300654E-1_qp,&
                         .1796443723688305731649384900158893966944E+0_qp,&
                        -.1392432216905901116427432216905901116427E+1_qp,&
                         .1340286404416839199447895100069013112491E+2_qp,&
                        -.1568482846260020173063651324520889738281E+3_qp,&
                         .2193103333333333333333333333333333333333E+4_qp,&
                        -.3610877125372498935717326521924223073648E+5_qp,&
                         .6914722688513130671083952507756734675533E+6_qp,&
                        -.1523822153940741619228336495888678051866E+8_qp,&
                         .3829007513914141414141414141414141414141E+9_qp,&
                       -.1088226603578439108901514916552510537473E+11_qp,&
                        .3473202837650022522522522522522522522523E+12_qp,&
                       -.1236960214226927445425171034927132488108E+14_qp]
      ! parameters from above reference

      z2 = z

      if(z % re < zero_k1) then

          z2 = cmplx(abs(z % re), - z % im, kind = dp) + 1

      end if

      d = hypot(z2 % re, z2 % im)
      z1 = z2
      m = 0

      if(d <= z_limit) then                       !for small |z|

          m = ceiling(z_limit - d)
          z1 = z2 + m

      end if

      zr = one / z1
      zr2 = zr * zr

      sum = (((a(20) * zr2 + a(19)) * zr2 + a(18)) * zr2 + a(17)) * zr2
      sum = (((sum + a(16)) * zr2 + a(15)) * zr2 + a(14)) * zr2
      sum = (((sum + a(13)) * zr2 + a(12)) * zr2 + a(11)) * zr2
      sum = (((sum + a(10)) * zr2 + a(9)) * zr2 + a(8)) * zr2
      sum = (((sum + a(7)) * zr2 + a(6)) * zr2 + a(5)) * zr2
      sum = (((sum + a(4)) * zr2 + a(3)) * zr2 + a(2)) * zr2
      sum = (sum + a(1)) * zr + ln2pi / 2 - z1 + (z1 - 0.5_qp) * log(z1)

      if(m /= 0) then

          s = cmplx(zero, zero, kind = qp)

          do i = 1, m

              s = s + log(cmplx(z1, kind = qp) - i)

          end do

          sum = sum - s

      end if

      if(z % re < zero_k1) then

          sum = log(pi) - log(sin(pi * z)) - sum
          m = ceiling((2 * z % re - 3) / 4)
          sum % im = sum % im + 2 * pi * m * sign(1.0_dp, z % im)

      end if

      res = cmplx(sum, kind = dp)
  end function l_gamma_cdp

end module libcoulcc__functions
