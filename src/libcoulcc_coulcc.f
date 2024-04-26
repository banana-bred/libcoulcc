! ================================================================================================================================ !
submodule (libcoulcc) libcoulcc_coulcc
  !! Submodule containing the implementation of the COULCC procedure

  implicit none

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine COULCC(XX, ETA1, ZLMIN, NL, FC, GC, FCP, GCP, SIG, MODE1, KFN, IFAIL)
    !! COMPLEX COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD
    !!
    !! A. R. Barnett           Manchester March 1981
    !! modified I.J. Thompson  Daresbury, Sept. 1983 for Complex Functions
    !!
    !! original program  RCWFN       in    CPC  8 (1974) 377-395
    !!                +  RCWFF       in    CPC 11 (1976) 141-142
    !!                +  COULFG      in    CPC 27 (1982) 147-166
    !! description of real algorithm in    CPC 21 (1981) 297-314
    !! description of complex algorithm    JCP XX (1985) YYY-ZZZ
    !! this version written up       in    CPC XX (1985) YYY-ZZZ
    !!
    !! COULCC returns F,G,G',G',SIG for complex XX, ETA1, and ZLMIN,
    !!  for NL integer-spaced lambda values ZLMIN to ZLMIN+NL-1 inclusive,
    !!  thus giving  complex-energy solutions to the Coulomb Schrodinger
    !!  equation,to the Klein-Gordon equation and to suitable forms of
    !!  the Dirac equation ,also spherical & cylindrical Bessel equations
    !!
    !! if /MODE1/= 1  get F,G,F',G'   for integer-spaced lambda values
    !!           = 2      F,G      unused arrays must be dimensioned in
    !!           = 3      F,  F'          call to at least length (1)
    !!           = 4      F
    !!           = 11 get F,H+,F',H+' ) if KFN=0, H+ = G + i.F        )
    !!           = 12     F,H+        )       >0, H+ = J + i.Y = H(1) ) in
    !!           = 21 get F,H-,F',H-' ) if KFN=0, H- = G - i.F        ) GC
    !!           = 22     F,H-        )       >0, H- = J - i.Y = H(2) )
    !!
    !!    if MODE1<0 then the values returned are scaled by an exponential
    !!               factor (dependent only on XX) to bring nearer unity
    !!               the functions for large /XX/, small ETA & /ZL/ < /XX/
    !!       Define SCALE = (  0        if MODE1 > 0
    !!                      (  IMAG(XX) if MODE1 < 0  &  KFN < 3
    !!                      (  REAL(XX) if MODE1 < 0  &  KFN = 3
    !!       then FC = EXP(-ABS(SCALE)) * ( F, j, J, or I)
    !!        and GC = EXP(-ABS(SCALE)) * ( G, y, or Y )
    !!              or EXP(SCALE)       * ( H+, H(1), or K)
    !!              or EXP(-SCALE)      * ( H- or H(2) )
    !!
    !! if  KFN  =  0,-1  complex Coulomb functions are returned   F & G
    !!          =  1   spherical Bessel      "      "     "       j & y
    !!          =  2 cylindrical Bessel      "      "     "       J & Y
    !!          =  3 modified cyl. Bessel    "      "     "       I & K
    !!
    !!         and where Coulomb phase shifts put in SIG if KFN=0 (not -1)
    !!
    !! The use of MODE and KFN is independent
    !!   (except that for KFN=3,  H(1) & H(2) are not given)
    !!
    !! With negative orders lambda, COULCC can still be used but with
    !!   reduced accuracy as CF1 becomes unstable. The user is thus
    !!   strongly advised to use reflection formulae based on
    !!   H+-(ZL,,) = H+-(-ZL-1,,) * exp +-i(sig(ZL)-sig(-ZL-1)-(ZL+1/2)pi)
    !!
    !! Precision:  results to within 2-3 decimals of 'machine accuracy',
    !!              but if CF1A fails because X too small or ETA too large
    !!              the F solution  is less accurate if it decreases with
    !!              decreasing lambda (e.g. for lambda.LE.-1 & ETA.NE.0)
    !!             RERR in COMMON/STEED/ traces the main roundoff errors.
    !!
    !!  COULCC is coded for real*8 on IBM or equivalent  ACCUR >= 10**-14
    !!         with a section of doubled REAL*16 for less roundoff errors.
    !!         (If no doubled precision available, increase JMAX to eg 100)
    !!  Use IMPLICIT COMPLEX*32 & REAL*16 on VS compiler ACCUR >= 10**-32
    !!  For single precision CDC (48 bits) reassign REAL(double)=REAL etc.
    !!
    !!  IFAIL  on input   = 0 : no printing of error messages
    !!                   ne 0 : print error messages on file 6
    !!  IFAIL  in output = -2 : argument out of range
    !!                   = -1 : one of the continued fractions failed,
    !!                          or arithmetic check before final recursion
    !!                   =  0 : All Calculations satisfactory
    !!                   ge 0 : results available for orders up to & at
    !!                            position NL-IFAIL in the output arrays.
    !!                   = -3 : values at ZLMIN not found as over/underflow
    !!                   = -4 : roundoff errors make results meaningless
    !! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !!     Machine dependent constants :
    !!
    !!     ACCUR    target bound on relative error (except near 0 crossings)
    !!               (ACCUR should be at least 100 * ACC8)
    !!     ACC8     smallest number with 1+ACC8 .ne.1 in REAL(double)  arithmetic
    !!     ACC16    smallest number with 1+ACC16.ne.1 in REAL*16 arithmetic
    !!     FPMAX    magnitude of largest floating point number * ACC8
    !!     FPMIN    magnitude of smallest floating point number / ACC8
    !!     FPLMAX   LOG(FPMAX)
    !!     FPLMIN   LOG(FPMIN)
    !!
    !!     ROUTINES CALLED :       log_gamma/digamma
    !!                             F20, CF1A, RCF, CF1C, CF2, F11, CF1R
    !!     Intrinsic functions :   MIN, MAX, SQRT, REAL, IMAG, ABS, LOG, EXP
    !!      (Generic names)        NINT, MOD, ATAN, ATAN2, COS, SIN, DCMPLX,
    !!                             SIGN, CONJG, INT, TANH
    !!     Note: Statement fntn.   NINTC = integer nearest to a complex no.
    !!
    !!     Parameters determining region of calculations :
    !!
    !!        R20      estimate of (2F0 iterations)/(CF2 iterations)
    !!        ASYM     minimum X/(ETA**2+L) for CF1A to converge easily
    !!        XNEAR    minimum ABS(X) for CF2 to converge accurately
    !!        LIMIT    maximum no. iterations for CF1, CF2, and 1F1 series
    !!        JMAX     size of work arrays for Pade accelerations
    !!        NDROP    number of successive decrements to define instability
    !!
    !! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    use libcoulcc_steed,               only: RERR, NFP, N11, NPQ, N20, KAS
    use libcoulcc_rcfcm1,              only: PK, EK, CLGAA, CLGAB, CLGBB, DSIG, TPK1, W, RL, FCL1, Q, GAM, HCL, HPL, FCM, HCL1, &
                                             ALPHA, BETA, PL
    use libcoulcc_constants,           only: ZERO, HALF, ONE, TWO, CI
    use iso_fortran_env,               only: dp => real64, qp => real128, stderr => error_unit
    use stdlib_specialfunctions_gamma, only: log_gamma

    ! IMPLICIT COMPLEX(dp) (A-H,O-Z)
    integer, intent(in)    :: NL
    integer, intent(in)    :: MODE1
    integer, intent(in)    :: KFN
    integer, intent(inout) :: IFAIL
    complex(dp), intent(in)   :: XX
    complex(dp), intent(in)   :: ETA1
    complex(dp), intent(in)   :: ZLMIN
    complex(dp), intent(out)  :: FC(nl)
    complex(dp), intent(out)  :: GC(nl)
    complex(dp), intent(out)  :: FCP(nl)
    complex(dp), intent(out)  :: GCP(nl)
    complex(dp), intent(out)  :: SIG(nl)

    integer, parameter :: jmax = 50
    integer, parameter :: LIMIT = 2000
    integer, parameter :: NDROP = 5

    real(dp), parameter :: FMAX = 1e60_dp
    real(dp), parameter :: FPMIN = 1e-60_dp
    real(dp), parameter :: FPLMAX = 140.0_dp
    real(dp), parameter :: FPLMIN = -140.0_dp
    real(dp), parameter :: R20 = 3.0_dp
    real(dp), parameter :: ASYM = 3.0_dp
    real(dp), parameter :: XNEAR = 0.5_dp
    ! real(dp), parameter :: ACCUR = 1e-14_dp
    ! real(dp), parameter :: ACC8 = 1e-14_dp
    ! real(dp), parameter :: ACC16 = 1e-33_qp
    real(dp), parameter :: ACC8 = epsilon(1.0_dp)
    real(dp), parameter :: ACC16 = epsilon(1.0_qp)
    real(dp), parameter :: ACCUR = 100 * ACC8

    complex(dp) :: XRCF(jmax,4)

    LOGICAL :: PR, ETANE0, IFCP, RLEL, DONEM, UNSTAB, ZLNEG, AXIAL, NOCF2, NPINT
    REAL(dp) :: ERR,RERR,ABSC,ACCUR,ACCT,ACC8,ACCH,ACC16,ACCB, XNEAR, &
                HPI,TLOG,FPMAX,FPMIN,FPLMIN,FPLMAX,                   &
                PACCQ,EPS,OFF,SCALE,SF,SFSH,TA,RK,OMEGA,ABSX

    integer :: I, ID, IH, KASE, L, L1, LAST, LF, LH, M1, MODE, MONO, N
    integer :: NINTC
    complex(dp) :: AA, AB, BB, CHI, CIK, CLL, DELL, DF, ETA, ETAI, ETAP, F, FIRST, FPL, F11V, F20V, FCL, FESL, FEST
    complex(dp) :: P, P11, PQ1, PQ2, PM, SIGMA, SL, THETA, THETAM, X, XI, XLOG, Z11, ZID, ZL, ZLL, ZLM, ZLOG, ZM1
    complex(dp) :: WW

    !COMMON /STEED/ RERR,NFP,N11,NPQ(2),N20,KAS(2)
    !**  common blocks are for information & storage only.
    !     (they are not essential to working of the code)
    ! COMMON /RCFCM1/ PK,EK,CLGAA,CLGAB,CLGBB,DSIG,TPK1,W,RL,FCL1,Q,GAM,HCL,HPL,FCM,HCL1,ALPHA,BETA,PL

    NINTC(WW) = NINT(REAL(REAL(WW)))
    ABSC(WW) = ABS(REAL(WW)) + ABS(IMAG(WW))
    NPINT(WW,ACCB) = ABSC(NINTC(WW)-WW) < ACCB .AND. REAL(WW) < HALF

    ! EQUIVALENCE (PK, XRCF(1, 1))
    XRCF(1, 1) = PK ! -- instead of EQUIVALENCE because PK is used as a submodule variable

    MODE = MOD(ABS(MODE1),10)
    IFCP = MOD(MODE,2) == 1
    PR = IFAIL /= 0
    IFAIL = -2
    N11   = 0
    NFP   = 0
    KAS(1)   = 0
    KAS(2)   = 0
    NPQ(1)   = 0
    NPQ(2)   = 0
    N20 = 0
    HPI = TWO*ATAN(ONE)
    TLOG = LOG(TWO)
    ! ACCUR = MAX(ACCUR, 50*ACC8) ! -- ACCUR is already set to 100 * ACC8
    ACCT = ACCUR * .5
    ! -- initialise the log_gamma function :
    ! CALL LOGAM(ACC8)
    ACCH  = SQRT(ACCUR)
    ACCB  = SQRT(ACCH)
    RERR = ACCT

    CIK = ONE
    IF(KFN >= 3) CIK = CI * SIGN(ONE,ACC8-IMAG(XX))
    X     = XX * CIK
    ETA   = ETA1
    IF(KFN > 0) ETA = ZERO
    ETANE0  = ABSC(ETA) > ACC8
    ETAI = ETA*CI
    DELL  = ZERO
    IF(KFN >= 2)  DELL = HALF
    ZM1   = ZLMIN - DELL
    SCALE = ZERO
    IF(MODE1 < 0) SCALE = IMAG(X)

    M1 = 1
    L1  = M1 + NL - 1
    RLEL = ABS(IMAG(ETA)) + ABS(IMAG(ZM1)) < ACC8
    ABSX = ABS(X)
    AXIAL = RLEL .AND. ABS(IMAG(X)) < ACC8 * ABSX
    IF(MODE <= 2 .AND. ABSX < FPMIN) GO TO 310
    XI  = ONE/X
    XLOG = LOG(X)
    ! -- log with cut along the negative real axis] see also OMEGA
    ID = 1
    DONEM = .FALSE.
    UNSTAB = .FALSE.
    LF = M1
    IFAIL = -1
10  ZLM = ZM1 + LF - M1
    ZLL = ZM1 + L1 - M1

    ! -- ZLL  is final lambda value, or 0.5 smaller for J,Y Bessels

    Z11 = ZLL
    IF(ID < 0) Z11 = ZLM
    P11 = CI*SIGN(ONE,ACC8-IMAG(ETA))
    LAST = L1

    ! -- Find phase shifts and Gamow factor at lambda = ZLL

    PK = ZLL + ONE
    AA = PK - ETAI
    AB = PK + ETAI
    BB = TWO*PK
    ZLNEG = NPINT(BB,ACCB)
    CLGAA = log_gamma(AA)
    ! CLGAA = CLOGAM(AA)
    CLGAB = CLGAA
    IF(ETANE0 .AND. .NOT. RLEL)  CLGAB = log_gamma(AB)
    ! IF(ETANE0 .AND. .NOT. RLEL)  CLGAB = CLOGAM(AB)
    IF(ETANE0 .AND.     RLEL)  CLGAB = CONJG(CLGAA)
    SIGMA = (CLGAA - CLGAB) * CI*HALF
    IF(KFN == 0) SIG(L1) = SIGMA
    ! IF( .NOT. ZLNEG) CLL = ZLL*TLOG- HPI*ETA - CLOGAM(BB) + (CLGAA + CLGAB) * HALF
    IF( .NOT. ZLNEG) CLL = ZLL*TLOG- HPI*ETA - log_gamma(BB) + (CLGAA+CLGAB)*HALF
    THETA  = X - ETA*(XLOG+TLOG) - ZLL*HPI + SIGMA

    TA = (IMAG(AA)**2+IMAG(AB)**2+ABS(REAL(AA))+ABS(REAL(AB)))*HALF
    IF(ID > 0 .AND. ABSX < TA*ASYM .AND. .NOT. ZLNEG) GO TO 20

    ! -- use CF1 instead of CF1A, if predicted to converge faster,
    !    (otherwise using CF1A as it treats negative lambda &
    !    recurrence-unstable cases properly)

    RK = SIGN(ONE, REAL(X) + ACC8)
    P =  THETA
    IF(RK < 0) P = -X + ETA*(LOG(-X)+TLOG)-ZLL*HPI-SIGMA
    F = RK * CF1A(X*RK, ETA*RK, ZLL, P, ACCT, JMAX, NFP, FEST, ERR, FPMAX, XRCF, XRCF(1, 3), XRCF(1, 4))
    PK = XRCF(1, 1) ! -- instead of EQUIVALENCE because PK is used as a submodule variable

    FESL = LOG(FEST) + ABS(IMAG(X))
    NFP = - NFP
    IF(NFP < 0   .OR. (UNSTAB .AND. ERR < ACCB)) GO TO 40
    IF( .NOT. ZLNEG .OR. UNSTAB .AND. ERR > ACCB)  GO TO 20
    IF(PR) WRITE(stderr,1060) '-L',ERR
    IF(ERR > ACCB) GO TO 280
    GO TO 40

    ! --  evaluate CF1  =  f   =  F'(ZLL,ETA,X)/F(ZLL,ETA,X)

20  IF(AXIAL) THEN
    ! -- REAL VERSION
        F = complex( &
        CF1R(real(X),real(ETA),real(ZLL),ACC8,SF,RK,ETANE0,LIMIT,ERR, &
        NFP, &
        ACCH,FPMIN,FPMAX,PR,'COULCC'),0._dp)
        FCL = SF
        TPK1= RK
    ELSE
    ! -- COMPLEX VERSION
        F = CF1C(X,ETA,ZLL,ACC8,FCL,TPK1,ETANE0,LIMIT,ERR,NFP, &
        ACCH,FPMIN,FPMAX,PR,'COULCC')
    ENDIF
    IF(ERR > ONE) GO TO 390

    ! --  Make a simple check for CF1 being badly unstable:

    IF(ID < 0) GO TO 30
    UNSTAB = REAL((ONE-ETA*XI)*CI*IMAG(THETA)/F) > ZERO &
     .AND. .NOT. AXIAL .AND. ABS(IMAG(THETA)) > -LOG(ACC8)*.5 &
     .AND. ABSC(ETA)+ABSC(ZLL) < ABSC(X)
    IF(UNSTAB) GO TO 60

    ! -- compare accumulated phase FCL with asymptotic phase for G(k+1) :
    !    to determine estimate of F(ZLL) (with correct sign) to start recu

30  W   =  X*X  *(HALF/TPK1 + ONE/TPK1**2) + ETA*(ETA-TWO*X)/TPK1
    FESL   = (ZLL+ONE) * XLOG + CLL - W - LOG(FCL)
    40 FESL = FESL - ABS(SCALE)
    RK   =        MAX(REAL(FESL), FPLMIN*HALF)
    FESL = DCMPLX(MIN(RK,   FPLMAX*HALF ) , IMAG(FESL))
    FEST= EXP(FESL)

    RERR = MAX(RERR, ERR, ACC8 * ABS(REAL(THETA)) )

    FCL = FEST
    FPL = FCL*F
    IF(IFCP) FCP(L1) = FPL
    FC (L1) = FCL

    ! -- downward recurrence to lambda = ZLM. array GC,if present,stores R

    I  = MAX(-ID, 0)
    ZL  = ZLL + I
    MONO = 0
    OFF = ABS(FCL)
    TA = ABSC(SIGMA)
    DO  L  = L1-ID,LF,-ID
        IF(ETANE0) THEN
            IF(RLEL) THEN
                DSIG = ATAN2(REAL(ETA),REAL(ZL))
                RL = SQRT(REAL(ZL)**2 + REAL(ETA)**2)
            ELSE
                AA = ZL - ETAI
                BB = ZL + ETAI
                IF(ABSC(AA) < ACCH .OR. ABSC(BB) < ACCH) GOTO 50
                DSIG = (LOG(AA) - LOG(BB)) * CI*HALF
                RL = AA * EXP(CI*DSIG)
            ENDIF
            IF(ABSC(SIGMA) < TA*HALF) THEN
            ! -- re-calculate SIGMA because of accumulating roundoffs:
                SL =(log_gamma(ZL+I-ETAI)-log_gamma(ZL+I+ETAI))*CI*HALF
                ! SL =(CLOGAM(ZL+I-ETAI)-CLOGAM(ZL+I+ETAI))*CI*HALF
                RL = (ZL - ETAI) * EXP(CI*ID*(SIGMA - SL))
                SIGMA = SL
                TA = ZERO
            ELSE
                SIGMA = SIGMA - DSIG * ID
            ENDIF
            TA = MAX(TA, ABSC(SIGMA))
            SL    =  ETA  + ZL*ZL*XI
            PL = ZERO
            IF(ABSC(ZL) > ACCH) PL = (SL*SL - RL*RL)/ZL
            FCL1  = (FCL *SL + ID*ZL*FPL)/RL
            SF = ABS(FCL1)
            IF(SF > FPMAX) GO TO 350
            FPL   = (FPL *SL + ID*PL*FCL)/RL
            IF(MODE <= 1) GCP(L+ID)= PL * ID
        ELSE
            ! -- ETA = 0, including Bessels.  NB RL==SL
            RL = ZL* XI
            FCL1 = FCL * RL + FPL*ID
            SF = ABS(FCL1)
            IF(SF > FPMAX) GO TO 350
            FPL  =(FCL1* RL - FCL) * ID
        ENDIF
            ! -- IF(ABSC(FCL1).LT.ABSC(FCL)) THEN
        IF(SF < OFF) THEN
            MONO = MONO + 1
        ELSE
            MONO = 0
        ENDIF
        FCL   =  FCL1
        OFF = SF
        FC(L) =  FCL
        IF(IFCP) FCP(L)  = FPL
        IF(KFN == 0) SIG(L) = SIGMA
        IF(MODE <= 2) GC(L+ID) = RL
        ZL = ZL - ID
        ! IF(MONO < NDROP) GO TO 70
        IF(MONO < NDROP) cycle
        ! IF(AXIAL .OR. REAL(ZLM)*ID > -NDROP .AND. .NOT. ETANE0) GO TO 70
        IF(AXIAL .OR. REAL(ZLM)*ID > -NDROP .AND. .NOT. ETANE0) cycle
        UNSTAB = .TRUE.

        ! -- take action if cannot or should not recur below this ZL:
        50 ZLM = ZL
        LF = L
        IF(ID < 0) GO TO 380
        IF( .NOT. UNSTAB) LF = L + 1
        IF(L+MONO < L1-2 .OR. ID < 0 .OR. .NOT. UNSTAB) GO TO 80
        ! otherwise, all L values (for stability) should be done
        ! in the reverse direction:
        GO TO 60
 ! 70 END DO
    END DO
    GO TO 80
    60 ID = -1
    LF = L1
    L1 = M1
    RERR = ACCT
    GO TO 10
    80 IF(FCL == ZERO) FCL = + ACC8
    F  = FPL/FCL

    ! -- Check, if second time around, that the 'f' values agree]

    IF(ID > 0) FIRST = F
    IF(DONEM) RERR = MAX(RERR, ABSC(F-FIRST)/ABSC(F))
    IF(DONEM) GO TO 90

    NOCF2 = .FALSE.
    THETAM  = X - ETA*(XLOG+TLOG) - ZLM*HPI + SIGMA

    ! --  on left x-plane, determine OMEGA by requiring cut on -x axis
    !     on right x-plane, choose OMEGA (using estimate based on THETAM)
    !     so H(omega) is smaller and recurs upwards accurately.
    !     (x-plane boundary is shifted to give CF2(LH) a chance to converge

    OMEGA = SIGN(ONE,IMAG(X)+ACC8)
    IF(REAL(X) >= XNEAR) OMEGA = SIGN(ONE,IMAG(THETAM)+ACC8)

    SFSH = EXP(OMEGA*SCALE - ABS(SCALE))
    OFF=EXP(MIN(TWO * MAX(ABS(IMAG(X)),ABS(IMAG(THETAM)), &
    ABS(IMAG(ZLM))*3 ) , FPLMAX) )
    EPS = MAX(ACC8 , ACCT * HALF / OFF)

    ! -- Try first estimated omega, then its opposite,
    !    to find the H(omega) linearly independent of F
    !    i.e. maximise  CF1-CF2 = 1/(F H(omega)) , to minimise H(omega)

    90 DO L=1,2
        LH = 1
        IF(OMEGA < ZERO) LH = 2
        PM = CI*OMEGA
        ETAP = ETA * PM
        IF(DONEM) GO TO 130
        PQ1 = ZERO
        PACCQ = ONE
        KASE = 0

    ! -- Check for small X, i.e. whether to avoid CF2 :

        IF(MODE >= 3 .AND. ABSX < ONE ) GO TO 190
        IF(MODE < 3 .AND. (NOCF2 .OR. ABSX < XNEAR .AND. &
        ABSC(ETA)*ABSX < 5 .AND. ABSC(ZLM) < 4)) THEN
            KASE = 5
            GO TO 120
        ENDIF

    ! -- Evaluate   CF2 : PQ1 = p + i.omega.q  at lambda = ZLM

        PQ1 = CF2(X,ETA,ZLM,PM,EPS,LIMIT,ERR,NPQ(LH),ACC8,ACCH, &
        PR,ACCUR,DELL,'COULCC')

        ERR = ERR * MAX(ONE,ABSC(PQ1)/MAX(ABSC(F-PQ1),ACC8) )
        IF(ERR < ACCH)       GO TO 110

    ! -- check if impossible to get F-PQ accurately because of cancellatio
        NOCF2 = REAL(X) < XNEAR .AND. ABS(IMAG(X)) < -LOG(ACC8)
    !    original guess for OMEGA (based on THETAM) was wrong
    !    Use KASE 5 or 6 if necessary if Re(X) < XNEAR
        OMEGA = - OMEGA
    END DO
    IF(UNSTAB) GO TO 360
    IF(REAL(X) < -XNEAR .AND. PR) WRITE(stderr,1060) '-X',ERR
    110 RERR = MAX(RERR,ERR)

    ! -- establish case of calculation required for irregular solution

    120 IF(KASE >= 5) GO TO 130
    IF(REAL(X) > XNEAR) THEN
    ! estimate errors if KASE 2 or 3 were to be used:
        PACCQ = EPS * OFF * ABSC(PQ1) / MAX(ABS(IMAG(PQ1)),ACC8)
    ENDIF
    IF(PACCQ < ACCUR) THEN
        KASE = 2
        IF(AXIAL) KASE = 3
    ELSE
        KASE = 1
        IF(NPQ(1) * R20 < JMAX)     KASE = 4
    ! i.e. change to kase=4 if the 2F0 predicted to converge
    ENDIF
    ! 130 GO TO (190,140,150,170,190,190),  ABS(KASE)
    130 select case(abs(KASE))
      case(1)      ; go to 190
      case(2)      ; go to 140
      case(3)      ; go to 150
      case(4)      ; go to 170
      case(5:6)    ; go to 190
      case default ; continue
    end select
    ! -- Evaluate   CF2 : PQ2 = p - i.omega.q  at lambda = ZLM   (Kase 2)
    140 IF( .NOT. DONEM)  PQ2 = CF2(X,ETA,ZLM,-PM,EPS,LIMIT,ERR,NPQ(3-LH),ACC8,ACCH, PR,ACCUR,DELL,'COULCC')

    P     = (PQ2 + PQ1) * HALF
    Q     = (PQ2 - PQ1) * HALF*PM
    GO TO 160
    150 P     = REAL(PQ1)
    Q     = IMAG(PQ1)

    ! --   With Kase = 3 on the real axes, P and Q are real & PQ2 = PQ1*
    PQ2 = CONJG(PQ1)

    ! -- solve for FCM = F at lambda = ZLM,then find norm factor W=FCM/FCL
    160 W   = (PQ1 - F) * (PQ2 - F)
    SF = EXP(-ABS(SCALE))
    FCM = SQRT(Q / W) * SF
    ! any SQRT given here is corrected by
    ! using sign for FCM nearest to phase of FCL
    IF(REAL(FCM/FCL) < ZERO) FCM  = - FCM
    GAM = (F - P)/Q
    TA = ABSC(GAM + PM)
    PACCQ= EPS * MAX(TA,ONE/TA)
    HCL = FCM * (GAM + PM) * (SFSH/(SF*SF))

    IF(PACCQ > ACCUR .AND. KASE > 0) THEN
        ! -- Consider a KASE = 1 Calculation
        F11V= F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)
        IF(ERR < PACCQ) GO TO 200
    ENDIF
    RERR=MAX(RERR,PACCQ)
    GO TO 230

    ! -- Arrive here if KASE = 4
    !    to evaluate the exponentially decreasing H(LH) directly.

    170 IF(DONEM) GO TO 180
    AA = ETAP - ZLM
    BB = ETAP + ZLM + ONE
    F20V = F20(AA,BB,-HALF*PM*XI, ACCT,JMAX,ERR,FPMAX,N20,XRCF)
    PK = XRCF(1, 1) ! -- instead of EQUIVALENCE because PK is used as a submodule variable
    IF(N20 <= 0) GO TO 190
    RERR = MAX(RERR,ERR)
    HCL = FPMIN
    IF(ABS(REAL(PM*THETAM)+OMEGA*SCALE) > FPLMAX) GO TO 330
    180 HCL = F20V * EXP(PM * THETAM + OMEGA*SCALE)
    FCM = SFSH / ((F - PQ1) * HCL )
    GO TO 230

    ! -- Arrive here if KASE=1   (or if 2F0 tried mistakenly & failed)

    !   for small values of X, calculate F(X,SL) directly from 1F1
    !   using REAL*16 arithmetic if possible.
    !   where Z11 = ZLL if ID>0, or = ZLM if ID<0

190 F11V = F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)

200 IF(N11 < 0) THEN
        ! -- F11 failed from BB = negative integer
        WRITE(stderr,1060) '-L',ONE
        GO TO 390
    ENDIF
    IF(ERR > PACCQ .AND. PACCQ < ACCB) THEN
        ! -- Consider a KASE 2 or 3 calculation :
        KASE = -2
        IF(AXIAL) KASE = -3
        GO TO 130
    ENDIF
    RERR = MAX(RERR, ERR)
    IF(ERR > FPMAX) GO TO 370
    IF(ID < 0) CLL = Z11*TLOG- HPI*ETA - log_gamma(BB) + log_gamma(Z11 + ONE + P11*ETA) - P11*SIGMA
    ! IF(ID < 0) CLL = Z11*TLOG- HPI*ETA - CLOGAM(BB) + CLOGAM(Z11 + ONE + P11*ETA) - P11*SIGMA
    EK   = (Z11+ONE)*XLOG - P11*X + CLL  - ABS(SCALE)
    IF(ID > 0) EK = EK - FESL + LOG(FCL)
    IF(REAL(EK) > FPLMAX) GO TO 350
    IF(REAL(EK) < FPLMIN) GO TO 340
    FCM = F11V * EXP(EK)

    IF(KASE >= 5) THEN
        IF(ABSC(ZLM+ZLM-NINTC(ZLM+ZLM)) < ACCH) KASE = 6

        ! --  For abs(X) < XNEAR, then CF2 may not converge accurately, so
        ! --      use an expansion for irregular soln from origin :

        SL = ZLM
        ZLNEG = REAL(ZLM) < -ONE + ACCB
        IF(KASE == 5 .OR. ZLNEG) SL = - ZLM - ONE
        PK = SL + ONE
        AA = PK - ETAP
        AB = PK + ETAP
        BB = TWO*PK
        CLGAA = log_gamma(AA)
        ! CLGAA = CLOGAM(AA)
        CLGAB = CLGAA
        IF(ETANE0) CLGAB = log_gamma(AB)
        ! IF(ETANE0) CLGAB = CLOGAM(AB)
        CLGBB = log_gamma(BB)
        ! CLGBB = CLOGAM(BB)
        IF(KASE == 6 .AND. .NOT. ZLNEG) THEN
            IF(NPINT(AA,ACCUR)) CLGAA = CLGAB - TWO*PM*SIGMA
            IF(NPINT(AB,ACCUR)) CLGAB = CLGAA + TWO*PM*SIGMA
        ENDIF
        CLL = SL*TLOG- HPI*ETA - CLGBB + (CLGAA + CLGAB) * HALF
        DSIG = (CLGAA - CLGAB) * PM*HALF
        IF(KASE == 6) P11 = - PM
        EK  = PK * XLOG - P11*X + CLL  - ABS(SCALE)
        SF = EXP(-ABS(SCALE))
        CHI = ZERO
        IF( .NOT. ( KASE == 5 .OR. ZLNEG ) ) GO TO 210

        ! -- Use  G(l)  =  (cos(CHI) * F(l) - F(-l-1)) /  sin(CHI)
        !     where CHI = sig(l) - sig(-l-1) - (2l+1)*pi/2

        CHI = SIGMA - DSIG - (ZLM-SL) * HPI
        F11V=F11(X,ETA,SL,P11,ACCT,LIMIT,0,ERR,NPQ(1),FPMAX,ACC8,ACC16)
        RERR = MAX(RERR,ERR)
        IF(KASE == 6) GO TO 210
        FESL = F11V * EXP( EK )
        FCL1 = EXP(PM*CHI) * FCM
        HCL = FCL1 - FESL
        RERR=MAX(RERR,ACCT*MAX(ABSC(FCL1),ABSC(FESL))/ABSC(HCL))
        HCL = HCL / SIN(CHI) * (SFSH/(SF*SF))
        GO TO 220

        ! -- Use the logarithmic expansion for the irregular solution (KASE 6)
        !       for the case that BB is integral so sin(CHI) would be zero.

        210 RL = BB - ONE
        N  = NINTC(RL)
        ZLOG = XLOG + TLOG - PM*HPI
        CHI = CHI + PM * THETAM + OMEGA * SCALE + AB * ZLOG
        AA  = ONE - AA
        IF(NPINT(AA,ACCUR)) THEN
            HCL = ZERO
        ELSE
            IF(ID > 0 .AND. .NOT. ZLNEG) F11V = FCM * EXP(-EK)
            HCL = EXP(CHI - CLGBB - log_gamma(AA)) * (-1)**(N+1) &
            ! HCL = EXP(CHI - CLGBB - CLOGAM(AA)) * (-1)**(N+1) &
            * ( F11V * ZLOG + &
            F11(X,ETA,SL,-PM,ACCT,LIMIT,2,ERR,NPQ(2),FPMAX,ACC8,ACC16))
            RERR = MAX(RERR,ERR)
        ENDIF
        IF(N > 0) THEN
            EK = CHI + log_gamma(RL) - CLGAB - RL*ZLOG
            ! EK = CHI + CLOGAM(RL) - CLGAB - RL*ZLOG
            DF =F11(X,ETA,-SL-ONE,-PM,ZERO,N,0,ERR,L,FPMAX,ACC8,ACC16)
            HCL = HCL + EXP(EK) * DF
        ENDIF

        220 PQ1 = F - SFSH/(FCM * HCL)
    ELSE
        IF(MODE <= 2) HCL = SFSH/((F - PQ1) * FCM)
        KASE = 1
    ENDIF

    ! --  Now have absolute normalisations for Coulomb Functions
    !         FCM & HCL  at lambda = ZLM
    !     so determine linear transformations for Functions required :

    230 IH = ABS(MODE1) / 10
    IF(KFN == 3) IH = (3-IMAG(CIK))/2  + HALF
    P11 = ONE
    IF(IH == 1) P11 = CI
    IF(IH == 2) P11 = -CI
    DF = - PM
    IF(IH >= 1) DF = - PM + P11
    IF(ABSC(DF) < ACCH) DF = ZERO

    ! -- Normalisations for spherical or cylindrical Bessel functions

    ALPHA = ZERO
    IF(KFN  == 1) ALPHA = XI
    IF(KFN  >= 2) ALPHA = XI*HALF
    BETA  = ONE
    IF(KFN  == 1) BETA  = XI
    IF(KFN  >= 2) BETA  = SQRT(XI/HPI)
    IF(KFN  >= 2 .AND. REAL(BETA) < ZERO) BETA  = - BETA

    AA = ONE
    IF(KFN > 0) AA = -P11 * BETA
    IF(KFN >= 3) THEN
        ! Calculate rescaling factors for I & K output
        P = EXP((ZLM+DELL) * HPI * CIK)
        AA= BETA * HPI * P
        BETA = BETA / P
        Q = CIK * ID
    ENDIF
    ! Calculate rescaling factors for GC output
    IF(IH == 0) THEN
        TA = ABS(SCALE) + IMAG(PM)*SCALE
        RK = ZERO
        IF(TA < FPLMAX) RK = EXP(-TA)
    ELSE
        TA = ABS(SCALE) + IMAG(P11)*SCALE

        IF(ABSC(DF) > ACCH .AND. TA > FPLMAX) GO TO 320
        IF(ABSC(DF) > ACCH) DF = DF * EXP(TA)
        SF = TWO * (LH-IH) * SCALE
        RK = ZERO
        IF(SF > FPLMAX) GO TO 320
        IF(SF > FPLMIN) RK = EXP(SF)
    ENDIF

    KAS((3-ID)/2) = KASE
    W = FCM / FCL
    IF(LOG(ABSC(W))+LOG(ABSC(FC(LF))) < FPLMIN) GO TO 340
    IF(MODE >= 3) GO TO 240
    IF(ABSC(F-PQ1) < ACCH*ABSC(F) .AND. PR) &
    WRITE(stderr,1020) LH,ZLM+DELL
    HPL = HCL * PQ1
    IF(ABSC(HPL) < FPMIN .OR. ABSC(HCL) < FPMIN) GO TO 330

    ! -- IDward recurrence from HCL,HPL(LF) (stored GC(L) is RL if reqd)
    ! -- renormalise FC,FCP at each lambda
    ! --    ZL   = ZLM - MIN(ID,0) here

  240 DO L = LF,L1,ID
        FCL = W* FC(L)
        IF(ABSC(FCL) < FPMIN) GO TO 340
        IF(IFCP) FPL = W*FCP(L)
        FC(L)  = BETA * FCL
        IF(IFCP) FCP(L) = BETA * (FPL - ALPHA * FCL) * CIK
        FC(L)  = TIDY(FC(L),ACCUR)
        IF(IFCP) FCP(L) = TIDY(FCP(L),ACCUR)
        IF(MODE >= 3) GO TO 260
        IF(L == LF)  GO TO 250
        ZL = ZL + ID
        ZID= ZL * ID
        RL = GC(L)
        IF(ETANE0)   THEN
            SL = ETA + ZL*ZL*XI
            IF(MODE == 1) THEN
                PL = GCP(L)
            ELSE
                PL = ZERO
                IF(ABSC(ZL) > ACCH) PL = (SL*SL - RL*RL)/ZID
            ENDIF
            HCL1     = (SL*HCL - ZID*HPL) / RL
            HPL      = (SL*HPL - PL *HCL) / RL
        ELSE
            HCL1 = RL * HCL - HPL * ID
            HPL  = (HCL - RL * HCL1) * ID
        ENDIF
        HCL      = HCL1
        IF(ABSC(HCL) > FPMAX) GO TO 320
        250 GC(L) = AA * (RK * HCL + DF * FCL)
        IF(MODE == 1) GCP(L) = (AA *(RK*HPL +DF*FPL) - ALPHA * GC(L)) *CI
        GC(L) = TIDY(GC(L),ACCUR)
        IF(MODE == 1) GCP(L) = TIDY(GCP(L),ACCUR)
        IF(KFN >= 3) AA = AA * Q
        260 IF(KFN >= 3) BETA = - BETA * Q
        LAST = MIN(LAST,(L1 - L)*ID)
    END DO

    ! -- Come here after all soft errors to determine how many L values ok

280 IF(ID > 0 .OR.  LAST == 0) IFAIL = LAST
    IF(ID < 0 .AND. LAST /= 0) IFAIL = -3
    ! IF(ID.GT.0 .OR.  LAST.EQ.0) write(stderr,'("r = ",E30.22)') XX
    ! IF(ID.LT.0 .AND. LAST.NE.0) write(stderr,'("r = ",E30.22)') XX


    ! -- Come here after ALL errors for this L range (ZLM,ZLL)

290 IF(ID > 0 .AND. LF /= M1) GO TO 300
    IF(IFAIL < 0) RETURN
    IF(RERR > ACCB) WRITE(stderr,1070) RERR
    IF(RERR > 0.1) IFAIL = -4
    RETURN

    ! -- so on first block, 'F' started decreasing monotonically,
    !                       or hit bound states for low ZL.
    !    thus redo M1 to LF-1 in reverse direction
    !     i.e. do CF1A at ZLMIN & CF2 at ZLM (midway between ZLMIN & ZLMAX

300 ID = -1
    IF( .NOT. UNSTAB) LF = LF - 1
    DONEM = UNSTAB
    LF = MIN(LF,L1)
    L1 = M1
    GO TO 10

    ! --    error messages
310 IF(PR) WRITE (stderr,1000) XX
1000 FORMAT(/' COULCC: CANNOT CALCULATE IRREGULAR SOLUTIONS FOR X =', &
    &  1P,2D10.2,', AS ABS(X) IS TOO SMALL'/)
    RETURN
320 IF(PR) WRITE(stderr,1010) ZL+DELL,'IR',HCL,'MORE',FPMAX
1010 FORMAT(' COULCC: AT ZL =',2F8.3,' ',A2,'REGULAR SOLUTION (',1P, &
    &  2E10.1,') WILL BE ',A4,' THAN',E10.1)
    GO TO 280
330 IF(PR) WRITE(stderr,1010) ZL+DELL,'IR',HCL,'LESS',FPMIN
    GO TO 280
340 IF(PR) WRITE(stderr,1010) ZL+DELL,'  ',FCL,'LESS',FPMIN
    GO TO 280
350 IF(PR) WRITE(stderr,1010) ZL+DELL,'  ',FCL,'MORE',FPMAX
    GO TO 280
1020 FORMAT('0COULCC WARNING: LINEAR INDEPENDENCE BETWEEN ''F'' AND '' &
    &(',I1,')'' IS LOST AT ZL =',2F7.2,' (EG. COULOMB EIGENSTATE, OR C &
    & 1 UNSTABLE)'/)
360 IF(PR) WRITE(stderr,1030) ZLL+DELL
1030 FORMAT(' COULCC: (ETA&L)/X TOO LARGE FOR CF1A, AND CF1 UNSTABLE AT  L =',2F8.2)
    GO TO 280
370 IF(PR) WRITE(stderr,1040) Z11,I
1040 FORMAT(' COULCC: OVERFLOW IN 1F1 SERIES AT ZL =',2F8.3,' AT TERM', I5)
    GO TO 390
380 IF(PR) WRITE(stderr,1050) ZLMIN,ZLM,ZLM+ONE,ZLMIN+NL-ONE
1050 FORMAT(' COULCC: BOTH BOUND-STATE POLES AND F-INSTABILITIES' &
    'OCCUR, OR MULTIPLE INSTABILITIES PRESENT.' &
    ,' TRY CALLING TWICE,  FIRST FOR ZL FROM',2F8.3,' TO',2F8.3, &
    ' (INCL.)',/,20X,     'SECOND FOR ZL FROM',2F8.3,' TO',2F8.3)
    ! GO TO 390
390 IFAIL = -1
    GO TO 290
1060 FORMAT('0COULCC WARNING: AS ''',A2,''' REFLECTION RULES NOT USED,' &
    ' ERRORS CAN BE UP TO',1P,D12.2/)
1070 FORMAT('0COULCC WARNING: OVERALL ROUNDOFF ERROR APPROX.',1P,E11.1)
  END subroutine coulcc

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  FUNCTION TIDY(Z,ACC)
    !! TIDY A COMPLEX NUMBER
    use iso_fortran_env, only: dp => real64
    REAL(dp) X,Y,ACC,AZ
    COMPLEX(dp) Z,TIDY

    X = REAL(Z)
    Y = IMAG(Z)
    AZ= (ABS(X) + ABS(Y)) * ACC * 5
    IF(ABS(X) < AZ) X = 0D+0
    IF(ABS(Y) < AZ) Y = 0D+0
    TIDY = DCMPLX(X,Y)
    RETURN
  END FUNCTION TIDY

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  FUNCTION CF1C(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP, ACCH,FPMIN,FPMAX,PR,CALLER)
    !! Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)
    !! using complex arithmetic

    use libcoulcc_constants, only: ONE, TWO
    use iso_fortran_env,     only: dp => real64, stderr => error_unit

    ! IMPLICIT COMPLEX(dp)(A-H,O-Z)
    complex(dp) :: CF1C
    LOGICAL :: PR,ETANE0
    integer, intent(in) :: LIMIT
    integer, intent(inout) :: NFP
    complex(dp), intent(in) :: X
    complex(dp), intent(in) :: ETA
    complex(dp), intent(in) :: ZL
    complex(dp), intent(inout) :: FCL
    complex(dp), intent(inout) :: TPK1
    real(dp) :: EPS
    real(dp) :: ERR
    real(dp) :: ACCH
    real(dp) :: FPMIN
    real(dp) :: FPMAX
    real(dp) :: ABSC
    real(dp) :: SMALL
    real(dp) :: RK
    real(dp) :: PX
    complex(dp) :: W
    complex(dp) :: D
    complex(dp) :: DF
    complex(dp) :: EK
    complex(dp) :: F
    complex(dp) :: PK
    complex(dp) :: PK1
    complex(dp) :: PK2
    complex(dp) :: SL
    complex(dp) :: TK
    complex(dp) :: XI
    complex(dp) :: RK2
    CHARACTER(6) :: CALLER
    ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))

    FCL = ONE
    XI = ONE/X
    PK  = ZL + ONE
    PX  = PK  + LIMIT
    EK  = ETA / PK
    10 continue
    RK2 =          ONE + EK*EK
    F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
    PK1 =  PK + ONE
    TPK1 = PK + PK1
    TK  = TPK1*(XI + EK/PK1)
    IF(ETANE0) THEN
    ! --   test ensures b1 .ne. zero for negative ETA etc.; fixup is exact
        IF(ABSC(TK) > ACCH)  GO TO 20
        FCL  = RK2/(ONE + (ETA/PK1)**2)
        SL   = TPK1*XI * (TPK1+TWO)*XI
        PK   =  TWO + PK
        GO TO 10
    ENDIF
    20 D   =  ONE/TK
    DF  = -FCL*RK2*D
    IF(REAL(PK) > REAL(ZL)+TWO) FCL = - RK2 * SL
    FCL = FCL * D * TPK1 * XI
    F   =  F  + DF

    ! --  begin CF1 loop on PK = k = lambda + 1

    RK    = ONE
    SMALL    = SQRT(FPMIN)
    30 PK    = PK1
    PK1 = PK1 + ONE
    TPK1 = PK + PK1
    IF(ETANE0) THEN
        EK  = ETA / PK
        RK2 =          ONE + EK*EK
    ENDIF
    TK  = TPK1*(XI + EK/PK1)
    D   =  TK - D*RK2
    IF(ABSC(D) > ACCH)             GO TO 40
    IF(PR) WRITE (stderr,1000) CALLER,D,DF,ACCH,PK,EK,ETA,X
    RK= RK +   ONE
    IF( RK > TWO )                  GO TO 50
    40 D     = ONE/D
    FCL = FCL * D * TPK1*XI
    IF(ABSC(FCL) < SMALL) FCL = FCL / SMALL
    IF(ABSC(FCL) > FPMAX) FCL = FCL / FPMAX
    DF  = DF*(D*TK - ONE)
    F   = F  + DF
    IF( REAL(PK) > PX ) GO TO 50
    IF(ABSC(DF) >= ABSC(F)*EPS)             GO TO 30
    NFP = PK - ZL - 1
    ERR = EPS * SQRT(REAL(NFP))
    CF1C = F
    RETURN
1000 FORMAT(/' ',A6,': CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ', /1X,1P,13D9.2/)
50  IF(PR) WRITE (stderr,1010) CALLER,LIMIT,ABS(X)
1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERATIONS AS ABS(X) =',F15.0)
    ERR = TWO
    RETURN
  END FUNCTION CF1C

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  FUNCTION CF2(X, ETA, ZL, PM, EPS, LIMIT, ERR, NPQ, ACC8, ACCH, PR, ACCUR, DELL, CALLER)
    !!                                    (omega)        (omega)
    !! -- Evaluate  CF2  = p + PM.q  =  H   (ETA,X)' / H   (ETA,X)
    !!                                    ZL             ZL
    !!     where PM = omega.i

    use libcoulcc_constants, only: ZERO, HALF, ONE, TWO
    use iso_fortran_env,     only: dp => real64, stderr => error_unit

    ! IMPLICIT COMPLEX(dp)(A-H,O-Z)
    LOGICAL :: PR
    complex(dp) :: CF2
    integer, intent(in) :: LIMIT
    integer, intent(inout) :: NPQ
    complex(dp), intent(in) :: X
    complex(dp), intent(in) :: ETA
    complex(dp), intent(in) :: ZL
    complex(dp), intent(in) :: PM
    complex(dp), intent(in) :: DELL
    real(dp) :: EPS
    real(dp) :: ERR
    real(dp) :: ACC8
    real(dp) :: ACCH
    real(dp) :: ACCUR
    real(dp) :: TA
    real(dp) :: RK
    real(dp) :: ABSC
    complex(dp) :: W
    complex(dp) :: AA
    complex(dp) :: BB
    complex(dp) :: DD
    complex(dp) :: DL
    complex(dp) :: E2MM1
    complex(dp) :: ETAP
    complex(dp) :: PQ
    complex(dp) :: RL
    complex(dp) :: WI
    complex(dp) :: XI
    CHARACTER(6) :: CALLER
    ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))

    TA = TWO*LIMIT
    E2MM1 = ETA*ETA + ZL*ZL + ZL
    ETAP = ETA * PM
    XI = ONE/X
    WI = TWO*ETAP
    RK = ZERO
    PQ = (ONE - ETA*XI) * PM
    AA = -E2MM1 + ETAP
    BB = TWO*(X - ETA + PM)
    RL = XI * PM
    IF(ABSC(BB) < ACCH) THEN
        RL = RL * AA / (AA + RK + WI)
        PQ = PQ + RL * (BB + TWO*PM)
        AA = AA + TWO*(RK+ONE+WI)
        BB = BB + (TWO+TWO)*PM
        RK = RK + (TWO+TWO)
    ENDIF
    DD = ONE/BB
    DL = AA*DD* RL
    10 PQ    = PQ + DL
    RK = RK + TWO
    AA = AA + RK + WI
    BB = BB + TWO*PM
    DD = ONE/(AA*DD + BB)
    DL = DL*(BB*DD - ONE)
    ERR = ABSC(DL)/ABSC(PQ)
    IF(ERR >= MAX(EPS,ACC8*RK*HALF) .AND. RK <= TA) GO TO 10

    NPQ   = RK/TWO
    PQ    = PQ + DL
    IF(PR .AND. NPQ >= LIMIT-1 .AND. ERR > ACCUR) &
    WRITE(stderr,1000) CALLER,INT(IMAG(PM)),NPQ,ERR,ZL+DELL
1000 FORMAT(' ',A6,': CF2(',I2,') NOT CONVERGED FULLY IN ',I7, &
    ' ITERATIONS, SO ERROR IN IRREGULAR SOLUTION =',1P,D11.2,' AT ZL =', 0P,2F8.3)
    CF2 = PQ
    RETURN
  END FUNCTION CF2

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  FUNCTION F11(X, ETA, ZL, P, EPS, LIMIT, KIND, ERR, NITS, FPMAX, ACC8, ACC16)
    !! evaluate the HYPERGEOMETRIC FUNCTION 1F1
    !!                                       i
    !!           F (AA;BB; Z) = SUM  (AA)   Z / ( (BB)  i] )
    !!          1 1              i       i            i
    !!
    !!    to accuracy EPS with at most LIMIT terms.
    !! If KIND = 0 : using extended precision but real arithmetic only,
    !!           1 : using normal precision in complex arithmetic,
    !!  or       2 : using normal complex arithmetic, but with CDIGAM factor
    !!
    !! where
    !!   AA = ZL+ONE - ETA*P
    !!   BB = TWO*(ZL+ONE)
    !! and
    !!   Z  = TWO*P*X

    use libcoulcc_constants, only: ZERO, ONE, TWO, CI
    use iso_fortran_env,     only: dp => real64, qp => real128

    ! IMPLICIT REAL(dp)(A-H,O-Z)

    integer, intent(in) :: LIMIT
    integer, intent(inout) :: NITS
    integer, intent(in) :: KIND
    real(dp), intent(in) :: EPS
    real(dp), intent(inout) :: ERR
    real(dp), intent(in) :: FPMAX
    real(dp), intent(in) :: ACC8
    real(dp), intent(in) :: ACC16
    complex(dp), intent(in) :: X
    complex(dp), intent(in) :: ETA
    complex(dp), intent(in) :: ZL
    complex(dp), intent(in) :: P
    integer :: I
    integer :: NINTC
    real(dp) :: ABSC
    real(dp) :: R
    real(dp) :: RK
    real(dp) :: TA
    complex(dp) :: AA
    complex(dp) :: BB
    complex(dp) :: Z
    complex(dp) :: F11
    ! complex(dp) :: CI
    COMPLEX(dp) :: DD
    complex(dp) :: G
    complex(dp) :: F
    complex(dp) :: AI
    complex(dp) :: BI
    complex(dp) :: T
    LOGICAL :: ZLLIN
    real(qp) :: AR
    real(qp) :: BR
    real(qp) :: GR
    real(qp) :: GI
    real(qp) :: DR
    real(qp) :: DI
    real(qp) :: TR
    real(qp) :: TI
    real(qp) :: UR
    real(qp) :: UI
    real(qp) :: FI
    real(qp) :: FI1
    real(qp) :: DEN
    ABSC(AA) = ABS(REAL(AA)) + ABS(IMAG(AA))
    NINTC(AA) = NINT(REAL(REAL(AA)))

    AA = ZL+ONE - ETA*P
    BB = TWO*(ZL+ONE)
    Z  = TWO*P*X

    ZLLIN = REAL(BB) <= ZERO .AND. ABS(BB-NINTC(BB)) < ACC8**0.25
    IF( .NOT. ZLLIN .OR. REAL(BB)+LIMIT < 1.5) GO TO 10
    NITS = -1
    RETURN
    10 IF(LIMIT <= 0) THEN
        F11 = ZERO
        ERR = ZERO
        NITS= 1
        RETURN
    ENDIF
    TA = ONE
    RK = ONE
    IF(KIND <= 0 .AND. ABSC(Z)*ABSC(AA) > ABSC(BB) * 1.0) THEN
        DR = ONE
        DI = ZERO
        GR = ONE
        GI = ZERO
        AR = REAL(AA)
        BR = REAL(BB)
        FI = ZERO
        DO I=2,LIMIT
            FI1 = FI + ONE
            TR = BR * FI1
            TI = IMAG(BB) * FI1
            DEN= ONE / (TR*TR + TI*TI)
            UR = (AR*TR + IMAG(AA)*TI) * DEN
            UI = (IMAG(AA)*TR - AR*TI) * DEN
            TR = UR*GR - UI*GI
            TI = UR*GI + UI*GR
            GR = REAL(Z) * TR - IMAG(Z)*TI
            GI = REAL(Z) * TI + IMAG(Z)*TR
            DR = DR + GR
            DI = DI + GI
            ERR = ABS(GR) + ABS(GI)
            IF(ERR > FPMAX) GO TO 60
            RK  = ABS(DR) + ABS(DI)
            TA = MAX(TA,RK)
            IF(ERR < RK*EPS .OR. I >= 4 .AND. ERR < ACC16) GO TO 30
            FI = FI1
            AR = AR + ONE
            BR = BR + ONE
        END DO

        30 F11 = DR + CI * DI
        ERR = ACC16 * TA / RK

    ELSE
    !* ---------------------------------- alternative code
    !*    If REAL*16 arithmetic is not available, (or already using it]),
    !*    then use KIND > 0
        G = ONE
        F = ONE
        ! IF(KIND >= 2) F = CDIGAM(AA) - CDIGAM(BB) - CDIGAM(G)
        IF(KIND >= 2) F = digamma(AA) - digamma(BB) - digamma(G)
        DD = F
        DO I=2,LIMIT
            AI = AA + (I-2)
            BI = BB + (I-2)
            R  = I-ONE
            G = G * Z * AI / (BI * R)
            IF(KIND >= 2) &
        !                              multiply by (psi(a+r)-psi(b+r)-psi(1+r))
            F = F + ONE/AI - ONE/BI - ONE/R
            T  = G * F
            DD = DD + T
            ERR = ABSC(T)
            IF(ERR > FPMAX) GO TO 60
            RK = ABSC(DD)
            TA = MAX(TA,RK)
            IF(ERR < RK*EPS .OR. ERR < ACC8 .AND. I >= 4) GO TO 50
        END DO

        50 ERR = ACC8 * TA / RK
        F11 = DD
    !* ------------------------------------------- end of alternative code
    ENDIF
    60 NITS = I
    RETURN
  END FUNCTION F11

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  FUNCTION CF1R(X, ETA, ZL, EPS, FCL, TPK1, ETANE0, LIMIT, ERR, NFP, ACCH, FPMIN, FPMAX, PR, CALLER)
    !! --    Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)
    !!        using real arithmetic

    ! IMPLICIT REAL(dp)(A-H,O-Z)
    use libcoulcc_constants, only: ONE, TWO
    use iso_fortran_env,     only: dp => real64

    integer, intent(in) :: LIMIT
    integer, intent(inout) :: NFP
    real(dp), intent(in) :: X
    real(dp), intent(in) :: ETA
    real(dp), intent(in) :: ZL
    real(dp), intent(in) :: EPS
    real(dp), intent(inout) :: TPK1
    real(dp), intent(inout) :: ERR
    real(dp), intent(in) :: ACCH
    real(dp), intent(in) :: FPMIN
    real(dp), intent(in) :: FPMAX
    real(dp), intent(out) :: FCL

    LOGICAL :: PR, ETANE0
    real(dp) :: CF1R, D, DF, EK, F, PK, PK1, PX, RK, RK2, SL, SMALL, TK, XI
    CHARACTER(6) :: CALLER

    FCL = ONE
    XI = ONE/X
    PK  = ZL + ONE
    PX  = PK  + LIMIT
    10 EK  = ETA / PK
    RK2 =          ONE + EK*EK
    F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
    PK1 =  PK + ONE
    TPK1 = PK + PK1
    TK  = TPK1*(XI + EK/PK1)
    IF(ETANE0) THEN
    ! --   test ensures b1 .ne. zero for negative ETA etc.; fixup is exact
        IF(ABS(TK) > ACCH)  GO TO 20
        FCL  = RK2/(ONE + (ETA/PK1)**2)
        SL   = TPK1*XI * (TPK1+TWO)*XI
        PK   =  TWO + PK
        GO TO 10
    ENDIF
    20 D   =  ONE/TK
    DF  = -FCL*RK2*D
    IF(PK > ZL+TWO) FCL = - RK2 * SL
    FCL = FCL * D * TPK1 * XI
    F   =  F  + DF

! --   begin CF1 loop on PK = k = lambda + 1

    RK    = ONE
    SMALL    = SQRT(FPMIN)
    30 PK    = PK1
    PK1 = PK1 + ONE
    TPK1 = PK + PK1
    IF(ETANE0) THEN
        EK  = ETA / PK
        RK2 =          ONE + EK*EK
    ENDIF
    TK  = TPK1*(XI + EK/PK1)
    D   =  TK - D*RK2
    IF(ABS(D) > ACCH)             GO TO 40
    IF(PR) WRITE (6,1000) CALLER,D,DF,ACCH,PK,EK,ETA,X
    RK= RK +   ONE
    IF( RK > TWO )                  GO TO 50
    40 D     = ONE/D
    FCL = FCL * D * TPK1*XI
    IF(ABS(FCL) < SMALL) FCL = FCL / SMALL
    IF(ABS(FCL) > FPMAX) FCL = FCL / FPMAX
    DF  = DF*(D*TK - ONE)
    F   = F  + DF
    IF( PK > PX ) GO TO 50
    IF(ABS(DF) >= ABS(F)*EPS)             GO TO 30
    NFP = PK - ZL - 1
    ERR = EPS * SQRT(REAL(NFP))
    CF1R = F
    RETURN
    1000 FORMAT(/' ',A6,': CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ', &
    /1X,1P,7D9.2/)
    50 IF(PR) WRITE (6,1010) CALLER,LIMIT,ABS(X)
    1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERATIONS AS ABS(X) =',F15.0)
    ERR = TWO
    RETURN
  END FUNCTION CF1R

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  FUNCTION F20(AA, BB, Z, EPS, JMAX, RE, FPMAX, N, X)
    !!     evaluate the HYPERGEOMETRIC FUNCTION 2F0
    !!                                             i
    !!            F (AA,BB;;Z) = SUM  (AA)  (BB)  Z / i]
    !!           2 0              i       i     i
    !!
    !!     to accuracy EPS with at most JMAX terms.
    !!
    !!     if the terms start diverging,
    !!     the corresponding continued fraction is found by RCF
    !!    & evaluated progressively by Steed's method to obtain convergence
    !!
    !!      useful number also input:  FPMAX = near-largest f.p. number

    ! IMPLICIT COMPLEX(dp)(A-H,O-Z)
    use libcoulcc_constants, only: ZERO, ONE
    use iso_fortran_env,     only: dp => real64

    complex(dp) :: F20
    integer, intent(in) :: JMAX
    integer, intent(out) :: N
    real(dp), intent(in) :: EPS
    real(dp), intent(in) :: FPMAX
    real(dp), intent(out) :: RE
    complex(dp), intent(in) :: Z
    complex(dp), intent(in) :: AA
    complex(dp), intent(in) :: BB
    complex(dp), intent(inout) :: X(JMAX,4)
    LOGICAL :: FINITE
    integer :: I
    integer :: IRCF
    integer :: IMAX
    integer :: K
    integer :: J
    integer :: MA
    integer :: MB
    integer :: NINTC
    real(dp) :: EP
    real(dp) :: AT
    real(dp) :: ATL
    real(dp) :: ABSC
    ! DATA ONE,ZERO / (1D+0,0D+0), (0D+0,0D+0) /
    complex(dp) :: D
    complex(dp) :: DF
    complex(dp) :: F
    complex(dp) :: CSUM
    complex(dp) :: WW
    ABSC(WW) = ABS(REAL(WW)) + ABS(IMAG(WW))
    NINTC(WW) = NINT(REAL(REAL(WW)))

    RE = 0.0
    X(1,1) = ONE
    CSUM = X(1,1)
    ATL = ABSC(X(1,1))
    F    = CSUM
    D = ONE
    DF   = CSUM
    J = 0
    EP = EPS * JMAX *10.
    MA = - NINTC(AA)
    MB = - NINTC(BB)
    FINITE = ABS(ABS(REAL(AA))-MA) < EP .AND. ABS(IMAG(AA)) < EP .OR. ABS(ABS(REAL(BB))-MB) < EP .AND. ABS(IMAG(BB)) < EP
    IMAX = JMAX
    IF(FINITE .AND. MA >= 0) IMAX = MIN(MA+1,IMAX)
    IF(FINITE .AND. MB >= 0) IMAX = MIN(MB+1,IMAX)
    DO I=2,IMAX
        X(I,1) = X(I-1,1) * Z * (AA+I-2) * (BB+I-2) / (I-1)
        IF(ABSC(X(I,1)) > FPMAX) GO TO 40
        AT = ABSC(X(I,1))
        IF(J == 0) THEN
            CSUM = CSUM + X(I,1)
            IF(AT < ABSC(CSUM)*EPS) GO TO 20
        ENDIF
        IF(FINITE) GO TO 10
        IF(J > 0 .OR. AT > ATL .OR. I >= JMAX-2) J = J + 1
        IF(J == 0) GO TO 10
        ! -- use IRCF instead of I so that we don't change the loop iterator
        ! CALL RCF(X(1,1),X(1,2),J,I,X(1,3),EPS)
        ! IF(I < 0) GO TO 40
        IRCF = I
        CALL RCF(X(1,1),X(1,2),J,IRCF,X(1,3),EPS)
        if(IRCF < 0) GO TO 40
        DO K=MAX(J,2),I
            D = ONE/(D*X(K,2) + ONE)
            DF = DF*(D - ONE)
            F = F + DF
            IF(ABSC(DF) < ABSC(F)*EPS) GO TO 30
            IF(DF == ZERO .AND. F == ZERO .AND. I >= 4) GO TO 30
        END DO
        J = I
     10 ATL = AT
    END DO
    IF( .NOT. FINITE) I = -JMAX
    20 N = I
    F20 = CSUM
    IF( .NOT. FINITE) RE  = AT / ABSC(CSUM)
    RETURN
    30 F20 = F
    RE = ABSC(DF) / ABSC(F)
    N = K
    RETURN
    40 I = 0
    GO TO 20
  END FUNCTION F20

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  FUNCTION CF1A(RHO,ETA,XL,PSI,EPS,NMAX,NUSED,FCL,RE,FPMAX,XX,G,C)
    !!    evaluate the ASYMPTOTIC EXPANSION for the
    !!           LOGARITHMIC DERIVATIVE OF THE REGULAR SOLUTION
    !!
    !! --        CF1A  =  f   =  F'(XL,ETA,RHO)/F(XL,ETA,RHO)
    !!
    !!     that is valid for REAL(RHO)>0, and best for RHO >> ETA**2, XL,
    !!     and is derived from the 2F0 expansions for H+ and H-
    !!     e.g. by Froeberg (Rev. Mod. Physics Vol 27, p399 , 1955)
    !!     Some lines of this subprogram are for convenience copied from
    !!          Takemasa, Tamura & Wolter CPC 17 (1979) 351.
    !!
    !!    Evaluate to accuracy EPS with at most NMAX terms.
    !!
    !!    If the terms start diverging,
    !!    the corresponding continued fraction is found by RCF
    !!    & evaluated progressively by Steed's method to obtain convergence
    !!
    !!     useful number also input:  FPMAX = near-largest f.p. number

    ! IMPLICIT COMPLEX(dp)(A-H,O-Z)

    use libcoulcc_constants, only: ZERO, ONE, TWO, CI
    use iso_fortran_env,     only: dp => real64

    complex(dp) :: CF1A
    integer, intent(in) :: NMAX
    integer, intent(out) :: NUSED
    complex(dp), intent(out) :: FCL
    complex(dp), intent(in) :: RHO
    complex(dp), intent(in) :: ETA
    complex(dp), intent(in) :: XL
    complex(dp), intent(in) :: PSI
    complex(dp), intent(inout) :: XX(2, NMAX)
    complex(dp), intent(inout) :: G(NMAX)
    complex(dp), intent(inout) :: C(NMAX)
    real(dp), intent(inout) :: RE
    real(dp), intent(in) :: EPS
    integer :: J
    integer :: K
    integer :: N
    integer :: NRCF
    real(dp) :: T1
    real(dp) :: T2
    real(dp) :: T3
    real(dp) :: AT
    real(dp) :: ATL
    real(dp) :: ABSC
    real(dp), intent(in) :: FPMAX
    complex(dp) :: C1
    complex(dp) :: C2
    complex(dp) :: COSL
    complex(dp) :: W
    complex(dp) :: F
    complex(dp) :: D
    complex(dp) :: DF
    complex(dp) :: DENOM
    complex(dp) :: ETASQ
    complex(dp) :: GLAST
    complex(dp) :: GSUM
    complex(dp) :: HPI
    complex(dp) :: SC
    complex(dp) :: SC1
    complex(dp) :: SC2
    complex(dp) :: SL
    complex(dp) :: SL1
    complex(dp) :: SL2
    complex(dp) :: TANL
    complex(dp) :: TC
    complex(dp) :: TC1
    complex(dp) :: TC2
    complex(dp) :: TL
    complex(dp) :: TL1
    complex(dp) :: TL2
    complex(dp) :: TL12
    complex(dp) :: XLL1
    ! DATA ZERO,ONE,TWO,CI / 0D+0, 1D+0, 2D+0, (0D+0,1D+0) /
    ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))

    HPI = TWO*ATAN(ONE)
    T1 = SIN(REAL(PSI))
    T2 = COS(REAL(PSI))
    ATL= TANH(IMAG(PSI))
!             GIVE COS(PSI)/COSH(IM(PSI)), WHICH ALWAYS HAS CORRECT SIG
    COSL = DCMPLX( T2 , -T1 * ATL )
    TANL = DCMPLX(T1,T2*ATL) / COSL
    RE = ZERO
    XLL1= XL*(XL+ONE)
    ETASQ = ETA*ETA
    SL1=ONE
    SL=SL1
    SC1=ZERO
    SC=SC1
    TL1=SC
    TL=TL1
    TC1=ONE-ETA/RHO
    TC=TC1
    FCL  = TL + SL*TANL
    G(1) = (TC + SC*TANL) / FCL
    GLAST = G(1)
    ATL = ABSC(GLAST)
    F    = GLAST
    D = ONE
    DF   = GLAST
    J = 0
    DO N=2,NMAX
        T1=N-1
        T2=TWO*T1-ONE
        T3=T1*(T1-ONE)
        DENOM=TWO*RHO*T1
        C1=(ETA*T2)/DENOM
        C2=(ETASQ+XLL1-T3)/DENOM
        SL2=C1*SL1-C2*TL1
        TL2=C1*TL1+C2*SL1
        SC2=C1*SC1-C2*TC1-SL2/RHO
        TC2=C1*TC1+C2*SC1-TL2/RHO
        SL=SL+SL2
        TL=TL+TL2
        SC=SC+SC2
        TC=TC+TC2
        SL1=SL2
        TL1=TL2
        SC1=SC2
        TC1=TC2
        FCL  =  TL + SL*TANL
        IF(ABSC(FCL) > FPMAX .OR. ABSC(FCL) < 1./FPMAX) GO TO 40
        GSUM = (TC + SC*TANL) / FCL
        G(N) = GSUM - GLAST
        GLAST = GSUM
        AT = ABSC(G(N))
        IF(AT < ABSC(GSUM)*EPS) GO TO 20
        IF(J > 0 .OR. AT > ATL .OR. N >= NMAX-2) J = J + 1
        IF(J == 0) GO TO 10
        ! -- use NRCF instead of I so that we don't change the loop iterator
        ! CALL RCF(G,C,J,N,XX,EPS)
        ! IF(N < 0) GO TO 40
        NRCF = N
        CALL RCF(G,C,J,NRCF,XX,EPS)
        IF(NRCF < 0) GO TO 40
        DO K=MAX(J,2),N
            D = ONE/(D*C(K) + ONE)
            DF = DF*(D - ONE)
            F = F + DF
            IF(ABSC(DF) < ABSC(F)*EPS) GO TO 30
            IF(DF == ZERO .AND. F == ZERO .AND. N >= 4) GO TO 30
        END DO
        J = N
       10 ATL = AT
    END DO
    K = -NMAX
    GO TO 30
    20 FCL = FCL * COSL
    CF1A = GSUM
    RE = AT / ABSC(GSUM)
    NUSED = N
    RETURN
    30 CF1A = F
    FCL = FCL * COSL
    RE = ABSC(DF) / ABSC(F)
    NUSED = K
    RETURN
    40 CF1A = G(1)
    FCL = 1.0
    RE = 1.0
    NUSED = 0
    RETURN
  END FUNCTION CF1A

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  SUBROUTINE RCF(A,B,IBEG,INUM,XX,EPS)
    !! converts polynomial A to the corresponding continued
    !!        fraction, in 'normal'  form with coefficients B
    !!        by the 'P algorithmn' of Patry & Gupta
    !!
    !!  A(z) = A1/z + A2/z**3 + A3/z**5 + ... + An/z**(2n-1)
    !!
    !!  B(z) = B1/z+ B2/z+ B3/z+ .../(z+ Bn/z)
    !!
    !! data:
    !!  A     vector A(k), k=1,INUM         input
    !!  B     vector B(k), k=IBEG,INUM      output
    !!  IBEG  order of first coef. calc.    input
    !!  INUM  order of A, even or odd       input
    !!  XX    auxiliary vector of length .ge. length of vector B
    !!        caller provides space for A,B,XX
    !!    Note that neither of the first two terms A(1) A(2) should be zero
    !!            & the user can start the calculation with any value of
    !!               IBEG provided the c.f. coefs have been already
    !!               calculated up to INUM = IBEG-1
    !!            & the method breaks down as soon as the absolute value
    !!               of a c.f. coef. is less than EPS.    At the time of th
    !!               break up XX(1) has been replaced by 1E-50, and INUM ha
    !!               been replaced by minus times the number of this coef.
    !!  algorithm: J.Patry & S.Gupta,
    !!             EIR-bericht nr. 247,
    !!             Eidg. Institut fur Reaktorforschung Wuerenlingen
    !!             Wueringlingen, Schweiz.
    !!             November 1973
    !!  see also:  Haenggi,Roesel & Trautmann,
    !!             Jnl. Computational Physics, vol 137, pp242-258 (1980)
    !!  note:      restart procedure modified by I.J.Thompson
    !!
    !! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    use libcoulcc_rcfcm2, only: MP12, X1, EVEN, M, M2M1
    use iso_fortran_env,  only: dp => real64

    ! IMPLICIT COMPLEX(dp)(A-H,O-Z)
    integer, intent(in) :: IBEG
    integer, intent(inout) :: INUM
    complex(dp), intent(in) :: A(100)
    complex(dp), intent(inout) :: B(100)
    complex(dp), intent(inout) :: XX(2, 100)
    integer :: K
    integer :: IBN
    complex(dp) :: X0
    REAL(dp) EPS
    ! COMMON /RCFCM2/ X1,M2M1,MP12,EVEN,M
!    ibn = ibeg + inum - 1
    IBN = INUM
!                            B(IBN) is last value set on this call
    IF(IBEG > 4 .AND. M /= IBEG-1) GO TO 90
    ! -- B(M) is last value set in previous call
    IF(IBEG > 4) GO TO 50
    IF(IBEG == 4) GO TO 20
    B(1) = A(1)
    IF(IBN >= 2) B(2) = - A(2)/A(1)
    IF(IBN < 3) GO TO 10
    X0 = A(3) / A(2)
    XX(2,1) = B(2)
    XX(1,1) = - X0
    XX(1,2) = 0.
    B(3) = -X0 - B(2)
    X0 = -B(3) * A(2)
    M = 3
    MP12 = 2
    EVEN = .TRUE.
    IF(IBN > 3) GO TO 20
    10 RETURN
    20 IF(ABS(B(3)) < EPS*ABS(X0)) GOTO 80
    M = 4
    30 X1 = A(M)
    M2M1 = MP12
    MP12 = M2M1 + 1
    IF(EVEN) MP12 = M2M1
    DO K=2,MP12
        X1 = X1 + A(M-K+1) * XX(1,K-1)
    END DO
    B(M) = - X1/X0
    IF(M >= IBN) RETURN
    50 IF(ABS(B(M)) < EPS*ABS(X0)) GO TO 80
    K = M2M1
    60 XX(2,K) = XX(1,K) + B(M) * XX(2,K-1)
    K = K-1
    IF(K > 1) GO TO 60
    XX(2,1) = XX(1,1) + B(M)
    DO K=1,M2M1
        X0 = XX(2,K)
        XX(2,K) = XX(1,K)
        XX(1,K) = X0
    END DO
    X0 = X1
    XX(1,M2M1+1) = 0.
    M = M+1
    EVEN = .NOT. EVEN
    GO TO 30
    80 INUM = -M
    ! XX(1,1) = 1.E-50
    ! PRINT 1000,M
    ! 1000 FORMAT('0RCF: ZERO CF COEFFICIENT AT POSITION ',I4/)
    RETURN
    90 PRINT 1000,M,IBEG-1
    1000 FORMAT('0RCF: LAST CALL SET M =',I4,', BUT RESTART REQUIRES',I4)
    STOP
  END SUBROUTINE RCF

    ! ------------------------------------------------------------------------------------------------------------------------------ !
  recursive impure elemental function digamma(z) result (res)
    !!
    !! Returns the digamma function for any complex number, excluding negative
    !! whole numbers, by reflection (z < 0), upward recurence (z % re < z_limit), or a
    !! truncated Stirling / de Moivre series
    !!
    use iso_fortran_env, only: dp => real64, qp => real128
    complex(dp), intent(in) :: z
    complex(dp) :: res
    integer, parameter :: n = 7
    integer :: k
    complex(qp) :: z2, zr, zr2, series
    complex(qp) :: res2
    real(dp), parameter :: z_limit = 10.0_dp
    real(dp), parameter :: zero_k1 = 0.0_dp
    real(qp), parameter :: one = 1.0_qp, two = 2.0_qp, pi = acos(-one)
    real(qp), parameter :: a(n) = [                     &
                         -one / 12.0_qp,             &
                          one / 120.0_qp,            &
                         -one / 252.0_qp,            &
                          one / 240.0_qp,            &
                         -one / 132.0_qp,            &
                          691.0_qp / 32760.0_qp, &
                         -one / 12.0_qp]

    z2 = z * one

    if( z % re <= zero_k1 ) then

      ! -- reflection ψ(1 - x) - ψ(x) = π cot(πx), including the imaginary axis
      res = digamma(1.0_dp - z)
      res2 = - pi * cotan(pi * z2)
      res = res + cmplx(res2 % re, res2 % im, kind = dp)

      return

    elseif( z % re > z_limit ) then

      zr  = one / z2
      zr2 = zr * zr
      series = log(z2) - zr / two  + sum( [(  a(k) * (zr2 ** k), k = 1, n )] )
      res = cmplx(series % re, series % im, kind = dp)

      return

    endif

    !-- recurrence relation ψ(z + 1) = ψ(z) + 1 / z
    res = digamma(z + 1) - 1 / z

  end function digamma

! ================================================================================================================================ !
end submodule libcoulcc_coulcc
! ================================================================================================================================ !
