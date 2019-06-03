module randomdistmod
    ! Random number generators and distributions module
    !
    ! This module contains subroutines and functions to compute the gamma and
    ! normal distribution, as well as a random number generator. Most of the
    ! subroutines have been translated from other languages to FORTRAN. The main
    ! references are
    !
    ! - http://www.johndcook.com/SimpleRNG.cpp for the random number generators
    ! - http://people.sc.fsu.edu/~jburkardt/f_src/prob/prob.html by John
    !   Burkardt for the gamma and normal distribution functions

    use parametersmod, only : sp, i4

    implicit none

    public  :: ran_seed
    public  :: ranur
    public  :: ran_normal
    public  :: ran_gamma
    public  :: gamma_cdf
    public  :: gamma_pdf

    private :: refill

    !----------------
    !module variables

    integer(i4), parameter :: qsiz  = 10   !41265_i4
    integer(i4), parameter :: cmul  = 69609_i4
    integer(i4), parameter :: coffs =   123_i4

    real(sp), parameter :: rng1 = 1. / (2. * real(huge(i4)))  !scales the random integer to -0.5,0.5
    real(sp), parameter :: rng2 = 1. / real(huge(i4))         !scales the random integer to -1,1

    real(sp), parameter :: one    = 1.
    real(sp), parameter :: half   = 0.5
    real(sp), parameter :: vsmall = tiny(1.)
    real(sp), parameter :: zero   = 0.

    type randomstate
        integer(i4), dimension(qsiz) :: q
        integer(i4) :: carry =       362_i4
        integer(i4) :: xcng  =   1236789_i4
        integer(i4) :: xs    = 521288629_i4  !default seed
        integer(i4) :: indx  = qsiz + 1
        logical     :: have  = .false.
    end type randomstate  !5+qsiz elements = 15 elements

    contains

    !-----------------------------------------------------------------------

    integer(i4) function ranu(state)
        ! Generates a uniformly distributed random 4 byte integer with the range (-huge(i4),+huge(i4))
        ! based on the 32-bit super KISS random number generator by George Marsaglia, published online
        ! and translated to Fortran 90 by user "mecej4" and Marsaglia, http://forums.silverfrost.com/viewtopic.php?t=1480
        ! Further modifications to pass the complete state of the generator as an argument by J.O. Kaplan, 2011

        implicit none

        type(randomstate), target, intent(inout) :: state ! state of the uniform random number generator

        integer(i4) :: supr

        integer(i4), pointer :: indx
        integer(i4), pointer :: xcng
        integer(i4), pointer :: xs
        integer(i4), pointer, dimension(:) :: q

        !---------------------

        indx => state%indx
        q    => state%q
        xcng => state%xcng
        xs   => state%xs

        !---

        if (indx <= qsiz) then
            supr = q(indx)
            indx = indx + 1
        else                     !reset the generator
            supr = refill(state)
        end if

        !---

        xcng = xcng * cmul + coffs
        xs   = ieor(xs,ishft(xs, 13))
        xs   = ieor(xs,ishft(xs,-17))
        xs   = ieor(xs,ishft(xs, -5))

        ranu = xcng + xs + supr

    end function ranu

    !-----------------------------------------------------------------------

    real(sp) function ranur(state)
        ! generate a random number in the range (0,1)

        implicit none

        type(randomstate), target, intent(inout) :: state ! state of the uniform random number generator

        !----

        ranur = real(ranu(state)) * rng1 + half

    end function ranur

    !-----------------------------------------------------------------------

    function refill(state) result(s)
        ! reset a random state

        implicit none

        type(randomstate), target, intent(inout) :: state ! state of the uniform random number generator

        integer(i4) :: s
        integer(i4) :: z
        integer(i4) :: h

        integer     :: i

        integer(i4), pointer :: indx
        integer(i4), pointer :: carry
        integer(i4), pointer, dimension(:) :: q

        !---------------------

        indx  => state%indx
        carry => state%carry
        q     => state%q

        !---

        do i = 1,qsiz

            h = iand(carry,1_i4)
            z = ishft(ishft(q(i),9),-1) + ishft(ishft(q(i),7),-1) + ishft(carry,-1)

            carry = ishft(q(i),-23) + ishft(q(i),-25) + ishft(z,-31)

            q(i) = not(ishft(z,1)+h)

        end do

        indx = 2
        s = q(1)

    end function refill

    !-----------------------------------------------------------------------

    subroutine ran_seed(sval,state)
        ! Set the seed of the random state

        implicit none

        integer(i4),               intent(in)    :: sval  ! state of the uniform random number generator
        type(randomstate), target, intent(inout) :: state ! the random state

        integer     :: i

        integer(i4), pointer :: xcng
        integer(i4), pointer :: xs
        integer(i4), pointer, dimension(:) :: q

        !---------------------

        q    => state%q
        xcng => state%xcng
        xs   => state%xs

        xs = sval

        !---

        do i = 1,qsiz

            xcng = xcng * cmul + coffs
            xs = ieor(xs,ishft(xs, 13))
            xs = ieor(xs,ishft(xs,-17))
            xs = ieor(xs,ishft(xs, -5))
            q(i) = xcng + xs

        end do

    end subroutine ran_seed

    !-----------------------------------------------------------------------

    subroutine ran_normal(state,nval)
        ! Sampler for the normal distribution centered at 0 with std. dev. of unity,
        ! based on Marsaglia polar method

        implicit none

        type(randomstate), intent(inout) :: state ! state of the uniform random number generator
        real(sp),          intent(out)   :: nval  ! output: The random number from the normal distribution

        !---

        real(sp), dimension(2) :: vals

        integer(i4), dimension(2) :: u
        real(sp),    dimension(2) :: v

        real(sp) :: s
        real(sp) :: a

        !---------------------

        do

            u(1) = ranu(state)
            u(2) = ranu(state)

            v = real(u) * rng2   !convert integer (-huge,+huge) to (-1,+1)

            s = sum(v**2)

            if (s < 1.) exit

        end do

            a = sqrt(-2. * log(s) / s)

            vals = v * a

        if (state%have) then

            state%have = .false.
            nval = vals(2)

        else

            nval = vals(1)
            state%have = .true.

        end if

    end subroutine ran_normal

    !----------------------------------------------------------------------------------------------------------

    real(sp) recursive function ran_gamma(state,first,shape,scale) result(ret)
        ! Select a random number from a Gamma distribution
        !
        ! adapted from the cpp adaptation of the Marsaglia & Tsang random gamma algorithm in:
        ! http://www.johndcook.com/SimpleRNG.cpp
        !
        ! Uses the algorithm in
        ! Marsaglia, G. and Tsang, W.W. (2000), *A simple method for generating
        ! gamma variables*, Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

        implicit none

        !arguments

        type(randomstate), intent(inout) :: state  ! state of the uniform random number generator
        logical,           intent(in)    :: first  ! flag if this is the first call to the distribution with this shape
        real(sp),          intent(in)    :: shape  ! shape parameter of the Gamma distribution (k or alpha, unitless) > 0
        real(sp),          intent(in)    :: scale  ! scale parameter of the Gamma distribution (theta = 1/beta)

        !local variables

        real(sp), save  :: c
        real(sp), save  :: d
        real(sp)        :: u
        real(sp)        :: v
        real(sp)        :: x

        !--------

        if (shape <= 0.) then

            write(0,*) 'shape parameter value must be positive'
            stop

        else if (shape >= 1.) then

            if (first) then
                d = shape - 1. / 3.
                c = 1. / sqrt(9. * d)
            end if

            do

                do  !generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

                  call ran_normal(state,x)
                  v = (1. + c * x)**3
                  if (v > 0.) exit

                end do

                u = ranur(state)  ! generate uniform variable u in the range (0,1)

                if (u < 1. - 0.0331 * x**4 .or. log(u) < half * x**2 + d*(1. - v + log(v))) then
                  ret = scale * d * v
                  exit
                end if
            end do

        else

            ret = scale * ran_gamma(state,first,shape + 1.,1.) * ranur(state)**(1. / shape)

        end if

    end function ran_gamma

    !----------------------------------------------------------------------------------------------------------

    real(sp) function ran_gamma_gp(state,first,shape,scale,thresh,shape_gp,scale_gp) result(ret)
        ! Select a random number from a hybrid Gamma-GP distribution

        implicit none

        !arguments

        type(randomstate), intent(inout) :: state     ! state of the uniform random number generator
        logical,           intent(in)    :: first     ! flag if this is the first call to the distribution with this shape
        real(sp),          intent(in)    :: shape     ! shape parameter of the Gamma distribution (k or alpha, unitless) > 0
        real(sp),          intent(in)    :: scale     ! scale parameter of the Gamma distribution (theta = 1/beta)
        real(sp),          intent(in)    :: thresh    ! the threshold above which to choose the GP distribution
        real(sp),          intent(in)    :: shape_gp  ! shape parameter of the GP distribution
        real(sp),          intent(in)    :: scale_gp  ! scale parameter of the GP distribution

        ret = ran_gamma(state,first,shape,scale)
        if (ret > thresh) then
            ret = ran_gp(state,shape_gp,scale_gp,thresh)
        endif

    end function ran_gamma_gp

    !-------------------------------

    real(sp) function ran_gp(state,shape,scale,loc)
        ! Select a random number from a generalized pareto (GP) distribution

        implicit none

        type(randomstate), intent(inout) :: state  ! state of the uniform random number generator
        real(sp),          intent(in)    :: shape  ! shape parameter of the GP distribution (k or alpha, unitless) > 0
        real(sp),          intent(in)    :: scale  ! scale parameter of the GP distribution (theta = 1/beta)
        real(sp),          intent(in)    :: loc    ! the location of the GP distribution

        real(sp)        :: u

        u = ranur(state)  ! generate uniform variable u in the range (0,1)

        if (shape == 0.0_sp) then
            ran_gp = loc - scale * log(u)
        else
            ran_gp = loc + scale * (U**(-shape) - 1) / shape
        endif

    end function ran_gp

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! -------------------------------------------------------------------------
    ! --------- Subroutines for gamma distribution and density functions ------
    ! -------------------------------------------------------------------------
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! The following functions and subroutines are taken from R, version 3.3.2 and
    ! converted from the original src/nmath/qgamma.c to fortran by Philipp Sommer,
    ! 2016, November 26th


    real(kind = 8) function qchisq_appr(p, nu, g, tol)
        ! chi-square approximation for the :f:func:`gamma_cdf_inv` function

        real(kind = 8), intent(in) :: p   ! the quantile
        real(kind = 8), intent(in) :: nu  ! twice the gamma shape
        real(kind = 8), intent(in) :: g   ! the logarithm of the gamma function at the gamma shape
        real(kind = 8), intent(in) :: tol ! the tolerance for the approximation

        real(kind = 8) :: alpha, a, c, ch, p1
        real(kind = 8) :: p2, q, t, x, lgam1pa

        alpha = 0.5 * nu
        c = alpha - 1.0

        p1 = log(p)

        if (nu < (-1.24) * p1) then
            ! for small chi-squared */
            !    log(alpha) + g = log(alpha) + log(gamma(alpha)) =
            !       = log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
            !    catastrophic cancellation when alpha << 1
            if (alpha < 0.5) then
                lgam1pa = gamma_log(alpha + 1.0)
            else
                lgam1pa = log(alpha) + g
            end if
            ch = exp((lgam1pa + p1)/alpha + log(2.0))
        else if (nu > 0.32) then !  using Wilson and Hilferty estimate
            call normal_cdf_inv(p, real(0.0, kind=8), real(1.0, kind=8), x)
            p1 = 2. / (9.0 * nu)
            ch = nu * (x * sqrt(p1) + 1.0 - p1) ** 3
            ! approximation for p tending to 1:
            if (ch > 2.2 * nu + 6) ch = -2.0 * (log(1 - p) - c * log(0.5 * ch) + g)
        else
            ch = 0.4
            a = log(1 - p) + g + c * log(2.0)
            do while (abs(q - ch) > tol * abs(ch))
                q = ch
                p1 = 1. / (1 + ch * (4.67 + ch))
                p2 = ch * (6.73 + ch * (6.66 + ch))
                t = -0.5 + (4.67 + 2 * ch) * p1 - (6.73 + ch*(13.32 + 3 * ch)) / p2
                ch = ch - (1 - exp(a + 0.5 * ch) * p2 * p1) / t
            end do
        end if
        qchisq_appr = ch

    end function qchisq_appr

    !-------------------------------

    real(kind = 8) function gamma_cdf_inv(p, alpha, scale)
        !	Compute the quantile function of the gamma distribution.
        !
        ! This function is based on the Applied Statistics Algorithm AS 91
        ! ("ppchi2") and via pgamma(.) AS 239.
        !
        ! References
        !	    Best, D. J. and D. E. Roberts (1975).
        !	    Percentage Points of the Chi-Squared Distribution.
        !	    Applied Statistics 24, page 385.
        !
        ! .. note::
        !
        !     Compared to the original R function, we do not use the final
        !     newton step which might lead to values going to infinity for
        !     quantiles close to 1

        real(kind = 8), intent(in) :: p ! the quantile between 0 and 1
        real(kind = 8), intent(in) :: alpha ! the shape of the gamma distribution
        real(kind = 8), intent(in) :: scale ! the scale of the gamma distribution

        real(kind = 8) :: a, b, c, g, ch, ch0, p1
        real(kind = 8) :: p2, q, s1, s2, s3, s4, s5, s6, t, x

        real(kind = 8), parameter :: EPS1 = 1.0e-2, EPS2 = 5.0e-7, &
                                     pMIN = 1.0e-25, pMAX = (1-1e-14), &
                                     pNEGinf = -10e34
        integer, parameter :: MAXIT = 1000
        integer :: i

        if (alpha == 0) then
            gamma_cdf_inv = 0
            goto 10002
        end if
        g = gamma_log(alpha)
        ! ----- Phase I : Starting Approximation
        ch = qchisq_appr(p, 2.0 * alpha, g, EPS1)

        if ((ch < EPS2) .or. (p > pMAX) .or. (p < pMIN)) then
            goto 10002
        end if

        ! ----- Phase II: Iteration
        ! Call pgamma() [AS 239]	and calculate seven term taylor series
        c = alpha - 1.0
        s6 = (120.0 + c * (346.0 + 127.0* c)) / 5040.0

        ch0 = ch  ! save initial approx.
        do i=1,MAXIT
            q = ch
            p1 = 0.5 * ch

            call gamma_cdf(p1 * scale, real(0.0, kind=8), scale, alpha, p2)
            p2 = p - p2

            if (ch <= 0.0) then
                ch = ch0
                goto 10001
            end if

            if ((p2 < pNEGinf) .or. (ch <= 0)) then
                ch = ch0
                goto 10001
            end if

            t = p2*exp(alpha * log(2.0) + g + p1 - c * log(ch))
            b = t / ch
            a = 0.5 * t - b * c
            s1 = (210.0 + a * (140.0 + a * (105.0 + a * (84.0 + a * ( &
                70.0 + 60.0 * a))))) / 420.0
            s2 = (420.0 +  a * (735.0 + a * (966.0 + a * (&
                1141.0 + 1278.0 * a)))) / 2520.
            s3 = (210.0 + a * (462.0 + a * (707.0 + 932.0 * a))) / 2520.0
            s4 = (252.0 + a * (672.0 + 1182.0 * a) + c * (294.0 +a * ( &
                889.0 + 1740.0 * a))) / 5040.0
            s5 = (84.0 + 2264.0 * a + c*(1175.0 + 606.0 * a)) / 2520.0

            ch = ch +  t * (1.0 + 0.5 * t * s1 - b * c * ( &
                s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))))

            if (abs(q - ch) < EPS2 * ch) goto 10001

            if (abs(q - ch) > 0.1 * ch) then
                if (ch < q) then
                    ch = 0.9 * q
                else
                    ch = 1.1 * q
                end if
            end if
        end do

    10001 continue

        gamma_cdf_inv = 0.5  * scale * ch

    10002 continue

    end function gamma_cdf_inv

    ! The following functions and subroutines to calculate the gamma distribution
    ! function (gamma_cdf) and the gamma density (gamma_pdf) are taken from
    ! John Burkardt, Department of Scientific Computing at Florida State University (FSU)
    ! http://people.sc.fsu.edu/~jburkardt/f_src/prob/prob.html
    ! on February 23rd, 2016. Last revised on 21 August 2013.

    subroutine gamma_cdf ( x, a, b, c, cdf )
        ! Evaluate the Gamma CDF.
        !
        !  Licensing:
        !    This code is distributed under the GNU LGPL license.
        !  Modified:
        !    02 January 2000
        !    Extracted: June, 2016
        !  Author:
        !    John Burkardt
        !    Extracted by Philipp Sommer
        implicit none

        real ( kind = 8 ), intent(in) :: x    ! the input value for which to compute the CDF
        real ( kind = 8 ), intent(in) :: a    ! the location (< `x`) of the gamma distribution (usually 0)
        real ( kind = 8 ), intent(in) :: b    ! the shape (> 0.0) of the distribution
        real ( kind = 8 ), intent(in) :: c    ! the scale (>0.0) of the distribution
        real ( kind = 8 ), intent(out) :: cdf ! the returned value of the CDF
        real ( kind = 8 ) p2
        real ( kind = 8 ) x2

        x2 = ( x - a ) / b
        p2 = c

        cdf = gamma_inc ( p2, x2 )

        return
    end subroutine gamma_cdf

    !------------------------------------------------------------------------------

    real ( kind = 8 ) function gamma_inc ( p, x )
        ! Compute the incomplete Gamma function.
        !
        ! Formulas:
        !
        !     .. math::
        !
        !         \Gamma_{inc}(P, 0) = 0
        !
        !     .. math::
        !
        !         \Gamma_{inc}(P, \infty) = 1.
        !
        !     .. math::
        !
        !         \Gamma_{inc}(P,X) = \int_0^x{T^{P-1} \exp{(-T)} \mathrm{d}t} / \Gamma(P)
        !
        ! Licensing:
        !     This code is distributed under the GNU LGPL license.
        ! Modified:
        !     - 01 May 2001
        !     - Extracted: June, 2016
        ! Author:
        !     - Original FORTRAN77 version by B L Shea.
        !     - FORTRAN90 version by John Burkardt
        !     - Extracted by Philipp Sommer
        ! Reference:
        !    BL Shea,
        !    Chi-squared and Incomplete Gamma Integral,
        !    Algorithm AS239,
        !    Applied Statistics,
        !    Volume 37, Number 3, 1988, pages 466-473.
        implicit none

        real ( kind = 8 ), intent(in) :: p  ! the exponent parameter (0.0 < P)
        real ( kind = 8 ), intent(in) :: x  ! the integral limit parameter. If X is less than or equal to 0, GAMMA_INC is returned as 0.
        real ( kind = 8 ) a
        real ( kind = 8 ) arg
        real ( kind = 8 ) b
        real ( kind = 8 ) c
        real ( kind = 8 ) cdf
        real ( kind = 8 ), parameter :: exp_arg_min = -88.0D+00
        real ( kind = 8 ), parameter :: overflow = 1.0D+37
        real ( kind = 8 ), parameter :: plimit = 1000.0D+00
        real ( kind = 8 ) pn1
        real ( kind = 8 ) pn2
        real ( kind = 8 ) pn3
        real ( kind = 8 ) pn4
        real ( kind = 8 ) pn5
        real ( kind = 8 ) pn6
        real ( kind = 8 ) rn
        real ( kind = 8 ), parameter :: tol = 1.0D-07
        real ( kind = 8 ), parameter :: xbig = 1.0D+08

        gamma_inc = 0.0D+00

        if ( p <= 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'GAMMA_INC - Fatal error!'
            write ( *, '(a)' ) '  Parameter P <= 0.'
            stop 1
        end if

        if ( x <= 0.0D+00 ) then
            gamma_inc = 0.0D+00
            return
        end if
        !
        !  Use a normal approximation if PLIMIT < P.
        !
        if ( plimit < p ) then
            pn1 = 3.0D+00 * sqrt ( p ) * ( ( x / p ) ** ( 1.0D+00 / 3.0D+00 ) &
              + 1.0D+00 / ( 9.0D+00 * p ) - 1.0D+00 )
            call normal_01_cdf ( pn1, cdf )
            gamma_inc = cdf
            return
        end if
        !
        !  Is X extremely large compared to P?
        !
        if ( xbig < x ) then
            gamma_inc = 1.0D+00
            return
        end if
        !
        !  Use Pearson's series expansion.
        !  (P is not large enough to force overflow in the log of Gamma.
        !
        if ( x <= 1.0D+00 .or. x < p ) then

            arg = p * log ( x ) - x - gamma_log ( p + 1.0D+00 )
            c = 1.0D+00
            gamma_inc = 1.0D+00
            a = p

            do

                a = a + 1.0D+00
                c = c * x / a
                gamma_inc = gamma_inc + c

                if ( c <= tol ) then
                    exit
                end if

            end do

            arg = arg + log ( gamma_inc )

            if ( exp_arg_min <= arg ) then
                gamma_inc = exp ( arg )
            else
                gamma_inc = 0.0D+00
            end if

        else
            !
            !  Use a continued fraction expansion.
            !
            arg = p * log ( x ) - x - gamma_log ( p )
            a = 1.0D+00 - p
            b = a + x + 1.0D+00
            c = 0.0D+00
            pn1 = 1.0D+00
            pn2 = x
            pn3 = x + 1.0D+00
            pn4 = x * b
            gamma_inc = pn3 / pn4

            do

                a = a + 1.0D+00
                b = b + 2.0D+00
                c = c + 1.0D+00
                pn5 = b * pn3 - a * c * pn1
                pn6 = b * pn4 - a * c * pn2

                if ( 0.0D+00 < abs ( pn6 ) ) then

                    rn = pn5 / pn6

                    if ( abs ( gamma_inc - rn ) <= min ( tol, tol * rn ) ) then

                        arg = arg + log ( gamma_inc )

                        if ( exp_arg_min <= arg ) then
                            gamma_inc = 1.0D+00 - exp ( arg )
                        else
                            gamma_inc = 1.0D+00
                        end if

                    return

                    end if

                    gamma_inc = rn

                end if

                pn1 = pn3
                pn2 = pn4
                pn3 = pn5
                pn4 = pn6
                !
                !  Rescale terms in continued fraction if terms are large.
                !
                if ( overflow <= abs ( pn5 ) ) then
                    pn1 = pn1 / overflow
                    pn2 = pn2 / overflow
                    pn3 = pn3 / overflow
                    pn4 = pn4 / overflow
                end if

            end do

        end if

        return
    end function gamma_inc

    !------------------------------------------------------------------------------

    real ( kind = 8 ) function gamma_log ( x )
        ! Calculate the natural logarithm of GAMMA ( X ).
        !
        ! Computation is based on an algorithm outlined in references 1 and 2.
        ! The program uses rational functions that theoretically approximate
        ! :math:`\log(\Gamma(X))` to at least 18 significant decimal digits.  The
        ! approximation for 12 < X is from Hart et al, while approximations
        ! for X < 12.0D+00 are similar to those in Cody and Hillstrom,
        ! but are unpublished.
        !
        ! The accuracy achieved depends on the arithmetic system, the compiler,
        ! intrinsic functions, and proper selection of the machine dependent
        ! constants.
        !
        !  Licensing:
        !    This code is distributed under the GNU LGPL license.
        !  Modified:
        !    - 16 June 1999
        !    - Extracted June, 2016
        !  Author:
        !    - Original FORTRAN77 version by William Cody, Laura Stoltz.
        !    - FORTRAN90 version by John Burkardt.
        !    - Extracted by Philipp Sommer
        !  Reference:
        !      - William Cody, Kenneth Hillstrom,
        !        Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
        !        Mathematics of Computation,
        !        Volume 21, 1967, pages 198-203.
        !      - Kenneth Hillstrom,
        !        ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
        !        May 1969.
        !      - John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
        !        Charles Mesztenyi, John Rice, Henry Thacher, Christoph Witzgall,
        !        Computer Approximations, Wiley, 1968.

        !  Local Parameters:
        !
        !    Local, real ( kind = 8 ) BETA, the radix for the floating-point
        !    representation.
        !
        !    Local, integer MAXEXP, the smallest positive power of BETA that overflows.
        !
        !    Local, real ( kind = 8 ) XBIG, the largest argument for which
        !    LN(GAMMA(X)) is representable in the machine, the solution to the equation
        !      LN(GAMMA(XBIG)) = BETA**MAXEXP.
        !
        !    Local, real ( kind = 8 ) FRTBIG, a rough estimate of the fourth root
        !    of XBIG.
        !
        !  Approximate values for some important machines are:
        !
        !                            BETA      MAXEXP         XBIG     FRTBIG
        !
        !  CRAY-1        (S.P.)        2        8191       9.62D+2461  3.13D+615
        !  Cyber 180/855 (S.P.)        2        1070       1.72D+319   6.44D+79
        !  IEEE (IBM/XT) (S.P.)        2         128       4.08D+36    1.42D+9
        !  IEEE (IBM/XT) (D.P.)        2        1024       2.55D+305   2.25D+76
        !  IBM 3033      (D.P.)       16          63       4.29D+73    2.56D+18
        !  VAX D-Format  (D.P.)        2         127       2.05D+36    1.20D+9
        !  VAX G-Format  (D.P.)        2        1023       1.28D+305   1.89D+76
        !
        implicit none

        real ( kind = 8 ), intent(in) :: x ! the argument of the Gamma function (> 0.0)

        real ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
            -1.910444077728D-03, &
             8.4171387781295D-04, &
            -5.952379913043012D-04, &
             7.93650793500350248D-04, &
            -2.777777777777681622553D-03, &
             8.333333333333333331554247D-02, &
             5.7083835261D-03 /)
        real ( kind = 8 ) corr
        real ( kind = 8 ), parameter :: d1 = -5.772156649015328605195174D-01
        real ( kind = 8 ), parameter :: d2 =  4.227843350984671393993777D-01
        real ( kind = 8 ), parameter :: d4 =  1.791759469228055000094023D+00
        integer ( kind = 4 ) i
        real ( kind = 8 ), parameter :: frtbig = 1.42D+09
        real ( kind = 8 ), parameter, dimension ( 8 ) :: p1 = (/ &
            4.945235359296727046734888D+00, &
            2.018112620856775083915565D+02, &
            2.290838373831346393026739D+03, &
            1.131967205903380828685045D+04, &
            2.855724635671635335736389D+04, &
            3.848496228443793359990269D+04, &
            2.637748787624195437963534D+04, &
            7.225813979700288197698961D+03 /)
        real ( kind = 8 ), parameter, dimension ( 8 ) :: p2 = (/ &
            4.974607845568932035012064D+00, &
            5.424138599891070494101986D+02, &
            1.550693864978364947665077D+04, &
            1.847932904445632425417223D+05, &
            1.088204769468828767498470D+06, &
            3.338152967987029735917223D+06, &
            5.106661678927352456275255D+06, &
            3.074109054850539556250927D+06 /)
        real ( kind = 8 ), parameter, dimension ( 8 ) :: p4 = (/ &
            1.474502166059939948905062D+04, &
            2.426813369486704502836312D+06, &
            1.214755574045093227939592D+08, &
            2.663432449630976949898078D+09, &
            2.940378956634553899906876D+10, &
            1.702665737765398868392998D+11, &
            4.926125793377430887588120D+11, &
            5.606251856223951465078242D+11 /)
        real ( kind = 8 ), parameter :: pnt68 = 0.6796875D+00
        real ( kind = 8 ), parameter, dimension ( 8 ) :: q1 = (/ &
            6.748212550303777196073036D+01, &
            1.113332393857199323513008D+03, &
            7.738757056935398733233834D+03, &
            2.763987074403340708898585D+04, &
            5.499310206226157329794414D+04, &
            6.161122180066002127833352D+04, &
            3.635127591501940507276287D+04, &
            8.785536302431013170870835D+03 /)
        real ( kind = 8 ), parameter, dimension ( 8 ) :: q2 = (/ &
            1.830328399370592604055942D+02, &
            7.765049321445005871323047D+03, &
            1.331903827966074194402448D+05, &
            1.136705821321969608938755D+06, &
            5.267964117437946917577538D+06, &
            1.346701454311101692290052D+07, &
            1.782736530353274213975932D+07, &
            9.533095591844353613395747D+06 /)
        real ( kind = 8 ), parameter, dimension ( 8 ) :: q4 = (/ &
            2.690530175870899333379843D+03, &
            6.393885654300092398984238D+05, &
            4.135599930241388052042842D+07, &
            1.120872109616147941376570D+09, &
            1.488613728678813811542398D+10, &
            1.016803586272438228077304D+11, &
            3.417476345507377132798597D+11, &
            4.463158187419713286462081D+11 /)
        real ( kind = 8 ) res
        real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
        real ( kind = 8 ), parameter :: xbig = 4.08D+36
        real ( kind = 8 ) xden
        real ( kind = 8 ) xm1
        real ( kind = 8 ) xm2
        real ( kind = 8 ) xm4
        real ( kind = 8 ) xnum
        real ( kind = 8 ) xsq
        !
        !  Return immediately if the argument is out of range.
        !
        if ( x <= 0.0D+00 .or. xbig < x ) then
            gamma_log = huge ( gamma_log )
            return
        end if

        if ( x <= epsilon ( x ) ) then

            res = -log ( x )

        else if ( x <= 1.5D+00 ) then

            if ( x < pnt68 ) then
                corr = - log ( x )
                xm1 = x
            else
                corr = 0.0D+00
                xm1 = ( x - 0.5D+00 ) - 0.5D+00
            end if

            if ( x <= 0.5D+00 .or. pnt68 <= x ) then

                xden = 1.0D+00
                xnum = 0.0D+00

                do i = 1, 8
                    xnum = xnum * xm1 + p1(i)
                    xden = xden * xm1 + q1(i)
                end do

                res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

            else

                xm2 = ( x - 0.5D+00 ) - 0.5D+00
                xden = 1.0D+00
                xnum = 0.0D+00
                do i = 1, 8
                    xnum = xnum * xm2 + p2(i)
                    xden = xden * xm2 + q2(i)
                end do

                res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

            end if

        else if ( x <= 4.0D+00 ) then

            xm2 = x - 2.0D+00
            xden = 1.0D+00
            xnum = 0.0D+00
            do i = 1, 8
                xnum = xnum * xm2 + p2(i)
                xden = xden * xm2 + q2(i)
            end do

            res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

        else if ( x <= 12.0D+00 ) then

            xm4 = x - 4.0D+00
            xden = - 1.0D+00
            xnum = 0.0D+00
            do i = 1, 8
                xnum = xnum * xm4 + p4(i)
                xden = xden * xm4 + q4(i)
            end do

            res = d4 + xm4 * ( xnum / xden )

        else

            res = 0.0D+00

            if ( x <= frtbig ) then

                res = c(7)
                xsq = x * x

                do i = 1, 6
                    res = res / xsq + c(i)
                end do

            end if

            res = res / x
            corr = log ( x )
            res = res + sqrtpi - 0.5D+00 * corr
            res = res + x * ( corr - 1.0D+00 )

        end if

        gamma_log = res

        return
    end function gamma_log

    !------------------------------------------------------------------------------

    subroutine gamma_pdf ( x, a, b, c, pdf )
        ! Evaluate the Gamma PDF.
        !
        ! .. math::
        !
        !     PDF(a,b,c;x) = \exp({-(x-a)/b}) \cdot ((x-a)/b)^{c-1} / (b \cdot \Gamma(c))
        !
        ! - GAMMA_PDF(A,B,C;X), where C is an integer, is the Erlang PDF.
        ! - GAMMA_PDF(A,B,1;X) is the Exponential PDF.
        ! - GAMMA_PDF(0,2,C/2;X) is the Chi Squared PDF with C degrees of freedom.
        !
        ! Licensing:
        !     This code is distributed under the GNU LGPL license.
        ! Modified:
        !     - 02 January 2000
        !     - Extracted: June, 2016
        ! Author:
        !     - John Burkardt
        !     - Extracted by Philipp Sommer
        implicit none

        real ( kind = 8 ), intent(in) :: x    ! the argument of the PDF. A <= X
        real ( kind = 8 ), intent(in) :: a    ! the location of the peak;  A is often chosen to be 0.0.
        real ( kind = 8 ), intent(in) :: b    ! the "scale" parameter; 0.0 < B, and is often 1.0.
        real ( kind = 8 ), intent(in) :: c    ! the "shape" parameter; 0.0 < C, and is often 1.0.
        real ( kind = 8 ), intent(out) :: pdf ! the returned value of the PDF.
        real ( kind = 8 ) y

        if ( x <= a ) then

            pdf = 0.0D+00

        else

            y = ( x - a ) / b

            pdf = y**( c - 1.0D+00 ) / ( b * r8_gamma ( c ) * exp ( y ) )

        end if

        return
    end subroutine gamma_pdf

    !------------------------------------------------------------------------------

    real ( kind = 8 ) function r8_gamma ( x )
        ! Evaluate Gamma(X) for a real argument.
        !
        ! This routine calculates the gamma function for a real argument X.
        !
        ! Computation is based on an algorithm outlined in reference 1.
        ! The program uses rational functions that approximate the gamma
        ! function to at least 20 significant decimal digits.  Coefficients
        ! for the approximation over the interval (1,2) are unpublished.
        ! Those for the approximation for 12 <= X are from reference 2.
        !
        ! Modified:
        !     - 11 February 2008
        !     - Extracted: June, 2016
        !
        ! Author:
        !     - Original FORTRAN77 version by William Cody, Laura Stoltz.
        !     - FORTRAN90 version by John Burkardt.
        !     - Extracted by Philipp Sommer
        ! Reference:
        !     - William Cody,
        !       An Overview of Software Development for Special Functions,
        !       in Numerical Analysis Dundee, 1975,
        !       edited by GA Watson, Lecture Notes in Mathematics 506,
        !       Springer, 1976.
        !     - John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
        !       Charles Mesztenyi, John Rice, Henry Thatcher,
        !       Christoph Witzgall, Computer Approximations, Wiley, 1968,
        !       LC: QA297.C64.

        implicit none

        real ( kind = 8 ), intent(in) :: x ! the argument of the function.
        !
        !  Coefficients for minimax approximation over (12, INF).
        !
        real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
           -1.910444077728D-03, &
            8.4171387781295D-04, &
           -5.952379913043012D-04, &
            7.93650793500350248D-04, &
           -2.777777777777681622553D-03, &
            8.333333333333333331554247D-02, &
            5.7083835261D-03 /)
        real ( kind = 8 ), parameter :: eps = 2.22D-16
        real ( kind = 8 ) fact
        real ( kind = 8 ), parameter :: half = 0.5D+00
        integer ( kind = 4 ) i
        integer ( kind = 4 ) n
        real ( kind = 8 ), parameter :: one = 1.0D+00
        real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
            -1.71618513886549492533811D+00, &
             2.47656508055759199108314D+01, &
            -3.79804256470945635097577D+02, &
             6.29331155312818442661052D+02, &
             8.66966202790413211295064D+02, &
            -3.14512729688483675254357D+04, &
            -3.61444134186911729807069D+04, &
             6.64561438202405440627855D+04 /)
        logical parity
        real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
        real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
            -3.08402300119738975254353D+01, &
             3.15350626979604161529144D+02, &
            -1.01515636749021914166146D+03, &
            -3.10777167157231109440444D+03, &
             2.25381184209801510330112D+04, &
             4.75584627752788110767815D+03, &
            -1.34659959864969306392456D+05, &
            -1.15132259675553483497211D+05 /)
        real ( kind = 8 ) res
        real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
        real ( kind = 8 ) sum
        real ( kind = 8 ), parameter :: twelve = 12.0D+00
        real ( kind = 8 ), parameter :: two = 2.0D+00
        real ( kind = 8 ), parameter :: xbig = 171.624D+00
        real ( kind = 8 ) xden
        real ( kind = 8 ), parameter :: xinf = 1.0D+30
        real ( kind = 8 ), parameter :: xminin = 2.23D-308
        real ( kind = 8 ) xnum
        real ( kind = 8 ) y
        real ( kind = 8 ) y1
        real ( kind = 8 ) ysq
        real ( kind = 8 ) z
        real ( kind = 8 ), parameter :: zero = 0.0D+00

        parity = .false.
        fact = one
        n = 0
        y = x
        !
        !  Argument is negative.
        !
        if ( y <= zero ) then

            y = - x
            y1 = aint ( y )
            res = y - y1

            if ( res /= zero ) then

                if ( y1 /= aint ( y1 * half ) * two ) then
                    parity = .true.
                end if

                fact = - pi / sin ( pi * res )
                y = y + one

            else

                res = xinf
                r8_gamma = res
                return

            end if

        end if
        !
        !  Argument is positive.
        !
        if ( y < eps ) then
            !
            !  Argument < EPS.
            !
            if ( xminin <= y ) then
                res = one / y
            else
                res = xinf
                r8_gamma = res
                return
            end if

        else if ( y < twelve ) then

            y1 = y
            !
            !  0.0 < argument < 1.0.
            !
            if ( y < one ) then

                z = y
                y = y + one
            !
            !  1.0 < argument < 12.0.
            !  Reduce argument if necessary.
            !
            else

                n = int ( y ) - 1
                y = y - real ( n, kind = 8 )
                z = y - one

            end if
            !
            !  Evaluate approximation for 1.0 < argument < 2.0.
            !
            xnum = zero
            xden = one
            do i = 1, 8
                xnum = ( xnum + p(i) ) * z
                xden = xden * z + q(i)
            end do

            res = xnum / xden + one
            !
            !  Adjust result for case  0.0 < argument < 1.0.
            !
            if ( y1 < y ) then

                res = res / y1
            !
            !  Adjust result for case 2.0 < argument < 12.0.
            !
            else if ( y < y1 ) then

                do i = 1, n
                    res = res * y
                    y = y + one
                end do

            end if

        else
            !
            !  Evaluate for 12.0 <= argument.
            !
            if ( y <= xbig ) then

                ysq = y * y
                sum = c(7)
                do i = 1, 6
                    sum = sum / ysq + c(i)
                end do
                sum = sum / y - y + sqrtpi
                sum = sum + ( y - half ) * log ( y )
                res = exp ( sum )

            else

                res = xinf
                r8_gamma = res
                return

            end if

        end if
        !
        !  Final adjustments and return.
        !
        if ( parity ) then
            res = - res
        end if

        if ( fact /= one ) then
            res = fact / res
        end if

        r8_gamma = res

        return
    end function r8_gamma

    !------------------------------------------------------------------------------

    subroutine normal_cdf_inv ( cdf, a, b, x )
        ! Invert the Normal CDF.
        !
        ! Licensing:
        !     This code is distributed under the GNU LGPL license.
        ! Modified:
        !     - 23 February 1999
        !     - Extracted: November, 2016
        ! Author:
        !     - John Burkardt
        !     - Extracted by Philipp Sommer
        !
        implicit none

        real ( kind = 8 ), intent(in) :: cdf ! the value of the CDF. 0.0 <= CDF <= 1.0.
        real ( kind = 8 ), intent(in) :: a   ! the mean of the pdf
        real ( kind = 8 ), intent(in) :: b   ! the standard deviation of the pdf
        real ( kind = 8 ), intent(out) :: x  ! the corresponding argument
        real ( kind = 8 ) x2

        if ( cdf < 0.0D+00 .or. 1.0D+00 < cdf ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'NORMAL_CDF_INV - Fatal error!'
            write ( *, '(a)' ) '  CDF < 0 or 1 < CDF.'
            stop 1
        end if

        call normal_01_cdf_inv ( cdf, x2 )

        x = a + b * x2

        return
    end subroutine normal_cdf_inv

    !------------------------------------------------------------------------------

    subroutine normal_01_cdf_inv ( p, x )
        ! Invert the standard normal CDF.
        !
        ! Licensing:
        !     This code is distributed under the GNU LGPL license.
        ! Modified:
        !    - 05 June 2007
        !    - Extracted: November, 2016
        ! Author:
        !     - Original FORTRAN77 version by Michael Wichura.
        !     - FORTRAN90 version by John Burkardt.
        !     - Extracted by Philipp Sommer
        ! Reference:
        !     Michael Wichura,
        !     Algorithm AS241:
        !     The Percentage Points of the Normal Distribution,
        !     Applied Statistics,
        !     Volume 37, Number 3, pages 477-484, 1988.
        !
        ! .. note::
        !
        !     The result is accurate to about 1 part in 10^16.
        !
        implicit none

        real ( kind = 8 ), intent(in) :: p ! the value of the cumulative probability densitity function.  0 < P < 1.
                                           ! If P is outside this range, an "infinite" value will be returned.
        real ( kind = 8 ), intent(out) :: x ! the normal deviate value with the property that the probability of a
                                            ! standard normal deviate being less than or equal to the value is P.

        real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
            3.3871328727963666080D+00, &
            1.3314166789178437745D+02, &
            1.9715909503065514427D+03, &
            1.3731693765509461125D+04, &
            4.5921953931549871457D+04, &
            6.7265770927008700853D+04, &
            3.3430575583588128105D+04, &
            2.5090809287301226727D+03 /)
        real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ &
            1.0D+00, &
            4.2313330701600911252D+01, &
            6.8718700749205790830D+02, &
            5.3941960214247511077D+03, &
            2.1213794301586595867D+04, &
            3.9307895800092710610D+04, &
            2.8729085735721942674D+04, &
            5.2264952788528545610D+03 /)
        real ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ &
            1.42343711074968357734D+00, &
            4.63033784615654529590D+00, &
            5.76949722146069140550D+00, &
            3.64784832476320460504D+00, &
            1.27045825245236838258D+00, &
            2.41780725177450611770D-01, &
            2.27238449892691845833D-02, &
            7.74545014278341407640D-04 /)
        real ( kind = 8 ), parameter :: const1 = 0.180625D+00
        real ( kind = 8 ), parameter :: const2 = 1.6D+00
        real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
            1.0D+00, &
            2.05319162663775882187D+00, &
            1.67638483018380384940D+00, &
            6.89767334985100004550D-01, &
            1.48103976427480074590D-01, &
            1.51986665636164571966D-02, &
            5.47593808499534494600D-04, &
            1.05075007164441684324D-09 /)
        real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ &
            6.65790464350110377720D+00, &
            5.46378491116411436990D+00, &
            1.78482653991729133580D+00, &
            2.96560571828504891230D-01, &
            2.65321895265761230930D-02, &
            1.24266094738807843860D-03, &
            2.71155556874348757815D-05, &
            2.01033439929228813265D-07 /)
        real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ &
            1.0D+00, &
            5.99832206555887937690D-01, &
            1.36929880922735805310D-01, &
            1.48753612908506148525D-02, &
            7.86869131145613259100D-04, &
            1.84631831751005468180D-05, &
            1.42151175831644588870D-07, &
            2.04426310338993978564D-15 /)
        real ( kind = 8 ) q
        real ( kind = 8 ) r
        real ( kind = 8 ), parameter :: split1 = 0.425D+00
        real ( kind = 8 ), parameter :: split2 = 5.0D+00

        if ( p <= 0.0D+00 ) then
            x = - huge ( x )
            return
        end if

        if ( 1.0D+00 <= p ) then
            x = huge ( x )
            return
        end if

        q = p - 0.5D+00

        if ( abs ( q ) <= split1 ) then

            r = const1 - q * q
            x = q * r8poly_value_horner ( 7, a, r ) / r8poly_value_horner ( 7, b, r )

        else

            if ( q < 0.0D+00 ) then
                r = p
            else
                r = 1.0D+00 - p
            end if

            if ( r <= 0.0D+00 ) then

                x = huge ( x )

            else

                r = sqrt ( - log ( r ) )

                if ( r <= split2 ) then

                    r = r - const2
                    x = r8poly_value_horner ( 7, c, r ) / r8poly_value_horner ( 7, d, r )

                else

                    r = r - split2
                    x = r8poly_value_horner ( 7, e, r ) / r8poly_value_horner ( 7, f, r )

                end if

            end if

            if ( q < 0.0D+00 ) then
                x = -x
            end if

        end if

        return
    end subroutine normal_01_cdf_inv

    !------------------------------------------------------------------------------

    real ( kind = 8 ) function r8poly_value_horner ( m, c, x )
        ! Evaluate a polynomial using Horner's method.
        !
        ! The polynomial
        !
        ! .. math::
        !
        !     p(x) = c_0 + c_1 * x + c_2 * x^2 + ... + c_m * x^m
        !
        ! is to be evaluated at the value X.
        !
        ! Licensing:
        !     This code is distributed under the GNU LGPL license.
        ! Modified:
        !     - 02 January 2014
        !     - Extracted: November, 2016
        ! Author:
        !     - John Burkardt
        !     - Extracted by Philipp Sommer
        implicit none

        integer ( kind = 4 ), intent(in) :: m   ! the degree
        real ( kind = 8 ), intent(in) :: c(0:m) ! the polynomial coefficients. C(I) is the coefficient of  :math:`X^I`
        real ( kind = 8 ), intent(in) :: x      ! the polynomial value

        integer ( kind = 4 ) i
        real ( kind = 8 ) value

        value = c(m)
        do i = m - 1, 0, -1
            value = value * x + c(i)
        end do

        r8poly_value_horner = value

        return
    end function r8poly_value_horner

    !------------------------------------------------------------------------------

    subroutine normal_01_cdf ( x, cdf )
        ! evaluate the Normal 01 CDF.
        !
        ! Licensing:
        !     This code is distributed under the GNU LGPL license.
        ! Modified:
        !    - 10 February 1999
        !    - Extracted: June, 2016
        ! Author:
        !     - John Burkardt
        !     - Extracted by Philipp Sommer
        ! Reference:
        !     AG Adams,
        !     Algorithm 39,
        !     Areas Under the Normal Curve,
        !     Computer Journal,
        !     Volume 12, pages 197-198, 1969.
        implicit none

        real ( kind = 8 ), intent(in) :: x ! the argument of the CDF.
        real ( kind = 8 ), intent(out) :: cdf ! the value of the CDF.

        real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
        real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
        real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
        real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
        real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
        real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
        real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
        real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
        real ( kind = 8 ), parameter :: b1 = 3.8052D-08
        real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
        real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
        real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
        real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
        real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
        real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
        real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
        real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
        real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
        real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
        real ( kind = 8 ) q
        real ( kind = 8 ) y
        !
        !  |X| <= 1.28.
        !
        if ( abs ( x ) <= 1.28D+00 ) then

            y = 0.5D+00 * x * x

            q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
              + a6 / ( y + a7 ) ) ) )
        !
        !  1.28 < |X| <= 12.7
        !
        else if ( abs ( x ) <= 12.7D+00 ) then

            y = 0.5D+00 * x * x

            q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
                + b2 / ( abs ( x ) + b3 &
                + b4 / ( abs ( x ) - b5 &
                + b6 / ( abs ( x ) + b7 &
                - b8 / ( abs ( x ) + b9 &
                + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
        !
        !  12.7 < |X|
        !
        else

            q = 0.0D+00

        end if
        !
        !  Take account of negative X.
        !
        if ( x < 0.0D+00 ) then
            cdf = q
        else
            cdf = 1.0D+00 - q
        end if

        return
    end subroutine normal_01_cdf

end module randomdistmod
