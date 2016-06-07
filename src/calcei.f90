subroutine calcei ( arg, result, int )

!*****************************************************************************80
!
!! CALCEI computes exponential integrals.
!
!  Discussion:
!
!    This routine computes the exponential integrals Ei(x),
!    E1(x), and  exp(-x)*Ei(x)  for real arguments  x  
!    where, if x > 0,
!      Ei(x) = integral (from t=-infinity to t=x) (exp(t)/t),  x > 0,
!    while, if x < 0,
!      Ei(x) = -integral (from t=-x to t=infinity) (exp(t)/t),  x < 0,
!    and where the first integral is a principal value integral.
!
!    The packet contains three function type subprograms: EI, EONE,
!    and EXPEI;  and one subroutine type subprogram: CALCEI.  The
!    calling statements for the primary entries are:
!      Y = EI(X),    where  X /= 0,
!      Y = EONE(X),      where  X > 0,
!      Y = EXPEI(X),     where  X /= 0,
!
!    and where the entry points correspond to the functions Ei(x),
!    E1(x), and exp(-x)*Ei(x), respectively.  The routine CALCEI
!    is intended for internal packet use only, all computations within
!    the packet being concentrated in this routine.  The function
!    subprograms invoke CALCEI with the statement
!      CALL CALCEI(ARG,RESULT,INT)
!    where the parameter usage is as follows
!
!      Function      Parameters for CALCEI
!      Call           ARG         RESULT         INT
!
!      EI(X)          X /= 0      Ei(X)            1
!      EONE(X)        X > 0      -Ei(-X)           2
!      EXPEI(X)       X /= 0      exp(-X)*Ei(X)    3
!
!    The main computation involves evaluation of rational Chebyshev
!    approximations published in Math. Comp. 22, 641-649 (1968), and
!    Math. Comp. 23, 289-303 (1969) by Cody and Thacher.  This
!    transportable program is patterned after the machine-dependent
!    FUNPACK packet  NATSEI,  but cannot match that version for
!    efficiency or accuracy.  This version uses rational functions
!    that theoretically approximate the exponential integrals to
!    at least 18 significant decimal digits.  The accuracy achieved
!    depends on the arithmetic system, the compiler, the intrinsic
!    functions, and proper selection of the machine-dependent
!    constants.
!
!  Modified:
!
!    10 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
! Explanation of machine-dependent constants.  Let
!
!   beta = radix for the floating-point system.
!   minexp = smallest representable power of beta.
!   maxexp = smallest power of beta that overflows.
!
! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.
!
!   XBIG = largest argument acceptable to EONE; solution to
!      equation:
!         exp(-x)/x * (1 + 1/x) = beta ** minexp.
!   XINF = largest positive machine number; approximately
!         beta ** maxexp
!   XMAX = largest argument acceptable to EI; solution to
!      equation:  exp(x)/x * (1 + 1/x) = beta ** maxexp.
!
! Error returns
!
!  The following table shows the types of error that may be
!  encountered in this routine and the function value supplied
!  in each case.
!
!   Error   Argument   function values for
!        Range     EI  EXPEI     EONE
!
!     UNDERFLOW  (-)X > XBIG     0    -     0
!     OVERFLOW  X >= XMAX    XINF  -     -
!     ILLEGAL X   X = 0   -XINF    -XINF     XINF
!     ILLEGAL X  X < 0   -    -     USE ABS(X)
!
  implicit none

  integer ( kind = 4 ) i,int
  real ( kind = 8 ) &
    a,arg,b,c,d,exp40,e,ei,f,four,fourty,frac,half,one,p, &
    plg,px,p037,p1,p2,q,qlg,qx,q1,q2,r,result,s,six,sump, &
    sumq,t,three,twelve,two,two4,w,x,xbig,xinf,xmax,xmx0, &
    x0,x01,x02,x11,y,ysq,zero
  dimension  a(7),b(6),c(9),d(9),e(10),f(10),p(10),q(10),r(10), &
    s(9),p1(10),q1(9),p2(10),q2(9),plg(4),qlg(4),px(10),qx(10)
!
!  Mathematical constants
!  EXP40 = exp(40)
!  X0 = zero of Ei
!  X01/X11 + X02 = zero of Ei to extra precision
!
  data zero,p037,half,one,two/0.0d0,0.037d0,0.5d0,1.0d0,2.0d0/, &
       three,four,six,twelve,two4/3.0d0,4.0d0,6.0d0,12.d0,24.0d0/, &
       fourty,exp40/40.0d0,2.3538526683701998541d17/, &
       x01,x11,x02/381.5d0,1024.0d0,-5.1182968633365538008d-5/, &
       x0/3.7250741078136663466d-1/
!
!  Machine-dependent constants
!
  data xinf/1.79d+308/,xmax/716.351d0/,xbig/701.84d0/
!
! Coefficients  for -1.0 <= X < 0.0
!
  data a/1.1669552669734461083368d2, 2.1500672908092918123209d3, &
     1.5924175980637303639884d4, 8.9904972007457256553251d4, &
     1.5026059476436982420737d5,-1.4815102102575750838086d5, &
     5.0196785185439843791020d0/
  data b/4.0205465640027706061433d1, 7.5043163907103936624165d2, &
     8.1258035174768735759855d3, 5.2440529172056355429883d4, &
     1.8434070063353677359298d5, 2.5666493484897117319268d5/
!
! Coefficients for -4.0 <= X < -1.0
!
  data c/3.828573121022477169108d-1, 1.107326627786831743809d+1, &
     7.246689782858597021199d+1, 1.700632978311516129328d+2, &
     1.698106763764238382705d+2, 7.633628843705946890896d+1, &
     1.487967702840464066613d+1, 9.999989642347613068437d-1, &
     1.737331760720576030932d-8/
  data d/8.258160008564488034698d-2, 4.344836335509282083360d+0, &
     4.662179610356861756812d+1, 1.775728186717289799677d+2, &
     2.953136335677908517423d+2, 2.342573504717625153053d+2, &
     9.021658450529372642314d+1, 1.587964570758947927903d+1, &
     1.000000000000000000000d+0/
!
! Coefficients for X < -4.0
!
  data e/1.3276881505637444622987d+2,3.5846198743996904308695d+4, &
     1.7283375773777593926828d+5,2.6181454937205639647381d+5, &
     1.7503273087497081314708d+5,5.9346841538837119172356d+4, &
     1.0816852399095915622498d+4,1.0611777263550331766871d03, &
     5.2199632588522572481039d+1,9.9999999999999999087819d-1/
  data f/3.9147856245556345627078d+4,2.5989762083608489777411d+5, &
     5.5903756210022864003380d+5,5.4616842050691155735758d+5, &
     2.7858134710520842139357d+5,7.9231787945279043698718d+4, &
     1.2842808586627297365998d+4,1.1635769915320848035459d+3, &
     5.4199632588522559414924d+1,1.0d0/
!
!  Coefficients for rational approximation to ln(x/a), |1-x/a| < .1
!
  data plg/-2.4562334077563243311d+01,2.3642701335621505212d+02, &
       -5.4989956895857911039d+02,3.5687548468071500413d+02/
  data qlg/-3.5553900764052419184d+01,1.9400230218539473193d+02, &
       -3.3442903192607538956d+02,1.7843774234035750207d+02/
!
! Coefficients for  0.0 < X < 6.0,
!  ratio of Chebyshev polynomials
!
  data p/-1.2963702602474830028590d01,-1.2831220659262000678155d03, &
     -1.4287072500197005777376d04,-1.4299841572091610380064d06, &
     -3.1398660864247265862050d05,-3.5377809694431133484800d08, &
      3.1984354235237738511048d08,-2.5301823984599019348858d10, &
      1.2177698136199594677580d10,-2.0829040666802497120940d11/
  data q/ 7.6886718750000000000000d01,-5.5648470543369082846819d03, &
      1.9418469440759880361415d05,-4.2648434812177161405483d06, &
      6.4698830956576428587653d07,-7.0108568774215954065376d08, &
      5.4229617984472955011862d09,-2.8986272696554495342658d10, &
      9.8900934262481749439886d10,-8.9673749185755048616855d10/
!
! J-fraction coefficients for 6.0 <= X < 12.0
!
  data r/-2.645677793077147237806d00,-2.378372882815725244124d00, &
     -2.421106956980653511550d01, 1.052976392459015155422d01, &
      1.945603779539281810439d01,-3.015761863840593359165d01, &
      1.120011024227297451523d01,-3.988850730390541057912d00, &
      9.565134591978630774217d00, 9.981193787537396413219d-1/
  data s/ 1.598517957704779356479d-4, 4.644185932583286942650d00, &
      3.697412299772985940785d02,-8.791401054875438925029d00, &
      7.608194509086645763123d02, 2.852397548119248700147d01, &
      4.731097187816050252967d02,-2.369210235636181001661d02, &
      1.249884822712447891440d00/
!
! J-fraction coefficients for 12.0 <= X < 24.0
!
  data p1/-1.647721172463463140042d00,-1.860092121726437582253d01, &
      -1.000641913989284829961d01,-2.105740799548040450394d01, &
      -9.134835699998742552432d-1,-3.323612579343962284333d01, &
       2.495487730402059440626d01, 2.652575818452799819855d01, &
      -1.845086232391278674524d00, 9.999933106160568739091d-1/
  data q1/ 9.792403599217290296840d01, 6.403800405352415551324d01, &
       5.994932325667407355255d01, 2.538819315630708031713d02, &
       4.429413178337928401161d01, 1.192832423968601006985d03, &
       1.991004470817742470726d02,-1.093556195391091143924d01, &
       1.001533852045342697818d00/
!
! J-fraction coefficients for  X >= 24.0
!
  data p2/ 1.75338801265465972390d02,-2.23127670777632409550d02, &
      -1.81949664929868906455d01,-2.79798528624305389340d01, &
      -7.63147701620253630855d00,-1.52856623636929636839d01, &
      -7.06810977895029358836d00,-5.00006640413131002475d00, &
      -3.00000000320981265753d00, 1.00000000000000485503d00/
  data q2/ 3.97845977167414720840d04, 3.97277109100414518365d00, &
       1.37790390235747998793d02, 1.17179220502086455287d02, &
       7.04831847180424675988d01,-1.20187763547154743238d01, &
      -7.99243595776339741065d00,-2.99999894040324959612d00, &
       1.99999999999048104167d00/

  x = arg

  if (x == zero) then

    ei = -xinf
    if (int == 2) ei = -ei

  else if ((x < zero) .or. (int == 2)) then
!
! Calculate EI for negative argument or for E1.
!
    y = abs(x)

    if (y <= one) then

      sump = a(7) * y + a(1)
      sumq = y + b(1)
      do i = 2, 6
        sump = sump * y + a(i)
        sumq = sumq * y + b(i)
      end do
      ei = log(y) - sump / sumq
      if (int == 3) ei = ei * exp(y)

    else if (y <= four) then

      w = one / y
      sump = c(1)
      sumq = d(1)
      do i = 2, 9
        sump = sump * w + c(i)
        sumq = sumq * w + d(i)
      end do
      ei = - sump / sumq
      if (int /= 3) ei = ei * exp(-y)

    else

      if ((y > xbig) .and. (int < 3)) then
        ei = zero
      else
        w = one / y
        sump = e(1)
        sumq = f(1)
        do i = 2, 10
          sump = sump * w + e(i)
          sumq = sumq * w + f(i)
        end do
        ei = -w * (one - w * sump / sumq )
        if (int /= 3) ei = ei * exp(-y)
      end if

    end if

    if (int == 2) ei = -ei

  else if (x < six) then
!
!  To improve conditioning, rational approximations are expressed
!  in terms of Chebyshev polynomials for 0 <= X < 6, and in
!  continued fraction form for larger X.
!
    t = x + x
    t = t / three - two
    px(1) = zero
    qx(1) = zero
    px(2) = p(1)
    qx(2) = q(1)
    do i = 2, 9
      px(i+1) = t * px(i) - px(i-1) + p(i)
      qx(i+1) = t * qx(i) - qx(i-1) + q(i)
    end do
    sump = half * t * px(10) - px(9) + p(10)
    sumq = half * t * qx(10) - qx(9) + q(10)
    frac = sump / sumq
    xmx0 = (x - x01/x11) - x02

    if (abs(xmx0) >= p037) then

      ei = log(x/x0) + xmx0 * frac

      if (int == 3) ei = exp(-x) * ei

    else
!
!  Special approximation to  ln(X/X0)  for X close to X0
!
      y = xmx0 / (x + x0)
      ysq = y*y
      sump = plg(1)
      sumq = ysq + qlg(1)
      do i = 2, 4
        sump = sump*ysq + plg(i)
        sumq = sumq*ysq + qlg(i)
      end do
      ei = (sump / (sumq*(x+x0)) + frac) * xmx0
      if (int == 3) ei = exp(-x) * ei

    end if

  else if (x < twelve) then

    frac = zero
    do i = 1, 9
      frac = s(i) / (r(i) + x + frac)
    end do
    ei = (r(10) + frac) / x
    if (int /= 3) ei = ei * exp(x)

  else if (x <= two4) then

    frac = zero
    do i = 1, 9
      frac = q1(i) / (p1(i) + x + frac)
    end do
    ei = (p1(10) + frac) / x
    if (int /= 3) ei = ei * exp(x)

  else

    if ((x >= xmax) .and. (int < 3)) then

      ei = xinf

    else

      y = one / x
      frac = zero
      do i = 1, 9
        frac = q2(i) / (p2(i) + x + frac)
      end do
      frac = p2(10) + frac
      ei = y + y * y * frac

      if (int /= 3) then
        if (x <= xmax-two4) then
          ei = ei * exp(x)
        else
!
!  Calculation reformulated to avoid premature overflow
!
          ei = (ei * exp(x-fourty)) * exp40
        end if

      end if

    end if

  end if

  result = ei

  return
end
