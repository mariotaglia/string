!*==ZSPOW.spg  processed by SPAG 6.72Dc at 15:10 on 27 Jul 2016
!   IMSL ROUTINE NAME   - ZSPOW
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDCFT5/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - SOLVE A SYSTEM OF NONLINEAR EQUATIONS
!
!   USAGE               - CALL ZSPOW (FCN,NSIG,N,ITMAX,PAR,X,FNORM,
!                           WK,IER)
!
!   ARGUMENTS    FCN    - THE NAME OF A USER-SUPPLIED SUBROUTINE WHICH
!                           EVALUATES THE SYSTEM OF EQUATIONS TO BE
!                           SOLVED. FCN MUST BE DECLARED EXTERNAL IN
!                           THE CALLING PROGRAM AND MUST HAVE THE
!                           FOLLOWING FORM,
!                             SUBROUTINE FCN(X,F,N,PAR)
!                             REAL X(N),F(N),PAR(1)
!                             F(1)=
!                              .
!                             F(N)=
!                             RETURN
!                             END
!                           GIVEN X(1)...X(N), FCN MUST EVALUATE THE
!                           FUNCTIONS F(1)...F(N) WHICH ARE TO BE MADE
!                           ZERO. X SHOULD NOT BE ALTERED BY FCN. THE
!                           PARAMETERS IN VECTOR PAR (SEE ARGUMENT
!                           PAR BELOW) MAY ALSO BE USED IN THE
!                           CALCULATION OF F(1)...F(N).
!                NSIG   - THE NUMBER OF DIGITS OF ACCURACY DESIRED
!                           IN THE COMPUTED ROOT. (INPUT)
!                N      - THE NUMBER OF EQUATIONS TO BE SOLVED AND
!                           THE NUMBER OF UNKNOWNS. (INPUT)
!                ITMAX  - THE MAXIMUM ALLOWABLE NUMBER OF ITERATIONS.
!                           (INPUT) THE MAXIMUM NUMBER OF CALLS TO FCN
!                           IS ITMAX*(N+1). SUGGESTED VALUE = 200.
!                PAR    - PAR CONTAINS A PARAMETER SET WHICH IS
!                           PASSED TO THE USER-SUPPLIED FUNCTION FCN.
!                           PAR MAY BE USED TO PASS ANY AUXILIARY
!                           PARAMETERS NECESSARY FOR COMPUTATION OF
!                           THE FUNCTION FCN. (INPUT)
!                X      - A VECTOR OF LENGTH N. (INPUT/OUTPUT) ON INPUT,
!                           X IS THE INITIAL APPROXIMATION TO THE ROOT.
!                           ON OUTPUT, X IS THE BEST APPROXIMATION TO
!                           THE ROOT FOUND BY ZSPOW.
!                FNORM  - ON OUTPUT, FNORM IS EQUAL TO
!                           F(1)**2+...F(N)**2 AT THE POINT X.
!                WK     - WORK VECTOR OF LENGTH N*(3*N+15)/2
!                IER    - ERROR PARAMETER. (OUTPUT)
!                         TERMINAL ERROR
!                           IER = 129 INDICATES THAT THE NUMBER OF
!                             CALLS TO FCN HAS EXCEEDED ITMAX*(N+1).
!                             THE USER MAY TRY A NEW INITIAL GUESS.
!                           IER = 130 INDICATES THAT NSIG IS TOO
!                             LARGE.  NO FURTHER IMPROVEMENT IN THE
!                             APPROXIMATE SOLUTION X IS POSSIBLE.
!                             THE USER SHOULD DECREASE NSIG.
!                           IER = 131 INDICATES THAT THE ITERATION
!                             HAS NOT MADE GOOD PROGRESS.  THE USER
!                             MAY TRY A NEW INITIAL GUESS.
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO,VBLA=SNRM2,ZSPWA,
!                           ZSPWB,ZSPWC,ZSPWD,ZSPWE,ZSPWF,ZSPWG
!                       - DOUBLE/UERTST,UGETIO,VBLA=DNRM2,ZSPWA,
!                           ZSPWB,ZSPWC,ZSPWD,ZSPWE,ZSPWF,ZSPWG
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ZSPOW(Nsig,N,Itmax,Par,X,Fnorm,Wk,Ier)
      IMPLICIT NONE
!*--ZSPOW87
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER Nsig , N , Itmax , Ier
      real*8 Par(*) , X(*) , Fnorm , Wk(*)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER index2 , index , info , i , j , lr , maxfev , ml , mode , &
            & mu , nfev , nprint
      REAL*8 epsfcn , factor , one , xtol , zero
      DATA factor , one , zero/1.0D2 , 1.0D0 , 0.0D0/
!                                  FIRST EXECUTABLE STATEMENT

      info = 0
!                                  CALL ZSPWA
      maxfev = Itmax*(N+1)
      xtol = 0.1**Nsig
      ml = N - 1
      mu = N - 1
      epsfcn = zero
      mode = 2
      DO j = 1 , N
         Wk(j) = one
      ENDDO
      nprint = 0
      lr = (N*(N+1))/2
      index = 7*N + lr
      CALL ZSPWA(N,X,Wk(6*N+1),xtol,maxfev,ml,mu,epsfcn,Wk(1),mode, &
               & factor,nprint,info,nfev,Wk(index+1),N,Wk(7*N+1),lr,    &
               & Wk(N+1),Wk(2*N+1),Wk(3*N+1),Wk(4*N+1),Wk(5*N+1),Par)
      IF ( info.EQ.5 ) info = 4
      Fnorm = 0.0
      DO i = 1 , N
         index2 = 6*N + i
         Fnorm = Fnorm + Wk(index2)*Wk(index2)
      ENDDO
      Ier = 0
      IF ( info.EQ.2 ) Ier = 129
      IF ( info.EQ.3 ) Ier = 130
      IF ( info.EQ.4 ) Ier = 131
      IF ( Ier.GT.0 ) CALL UERTST(Ier,'ZSPOW ')
      CONTINUE
      END
!*==ZSPWA.spg  processed by SPAG 6.72Dc at 15:10 on 27 Jul 2016
!DECK,ZSPWA
!   IMSL ROUTINE NAME   - ZSPWA
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDCFT5/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - SINGLE/VBLA=SNRM2,ZSPWB,ZSPWC,ZSPWD,ZSPWE,
!                           ZSPWF,ZSPWG
!                       - DOUBLE/VBLA=DNRM2,ZSPWB,ZSPWC,ZSPWD,ZSPWE,
!                           ZSPWF,ZSPWG
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ZSPWA(N,X,Fvec,Xtol,Maxfev,Ml,Mu,Epsfcn,Diag,Mode, &
                     & Factor,Nprint,Info,Nfev,Fjac,Ldfjac,R,Lr,Qtf,Wa1,&
                     & Wa2,Wa3,Wa4,Par)
      IMPLICIT NONE
!*--ZSPWA164
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER N , Maxfev , Ml , Mu , Mode , Nprint , Info , Nfev ,      &
            & Ldfjac , Lr
      REAL*8 X(N) , Fvec(N) , Xtol , Epsfcn , Diag(N) , Factor ,        &
           & Fjac(Ldfjac,N) , R(Lr) , Qtf(N) , Wa1(N) , Wa2(N) , Wa3(N) &
           & , Wa4(N) , Par(*)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER iflag , iter , iwa(1) , i , jm1 , j , l , msum , ncfail , &
            & ncsuc , nslow1 , nslow2
      REAL*8 actred , delta , epsmch , fnorm1 , fnorm , one , p0001 ,   &
           & p001 , p1 , p5 , pnorm , prered , ratio , spmpar , sum ,   &
           & temp , xnorm , zero
      REAL*8 DNRM2
      LOGICAL jeval , sing
      DATA spmpar/.710543D-14/
      DATA one , p1 , p5 , p001 , p0001 , zero/1.0D0 , 1.0D-1 , 5.0D-1 ,&
         & 1.0D-3 , 1.0D-4 , 0.0D0/
!                                  EPSMCH IS THE MACHINE PRECISION.
!                                  FIRST EXECUTABLE STATEMENT
      epsmch = spmpar
      Info = 0
      iflag = 0
      Nfev = 0
!                                  CHECK THE INPUT PARAMETERS FOR
!                                  ERRORS.
      IF ( N.LE.0 .OR. Xtol.LT.zero .OR. Maxfev.LE.0 .OR. Ml.LT.0 .OR.  &
         & Mu.LT.0 .OR. Factor.LE.zero .OR. Ldfjac.LT.N .OR.            &
         & Lr.LT.(N*(N+1))/2 ) GOTO 1300
      IF ( Mode.NE.2 ) GOTO 100
      DO j = 1 , N
         IF ( Diag(j).LE.zero ) GOTO 1300
      ENDDO
 100  CONTINUE
!                                  EVALUATE THE FUNCTION AT THE STARTING
!                                  POINT AND CALCULATE ITS NORM.
      iflag = 1
      CALL FCN(X,Fvec,N,Par)
      Nfev = 1
      IF ( iflag.LT.0 ) GOTO 1300
      fnorm = DNRM2(N,Fvec,1)
!                                  DETERMINE THE NUMBER OF CALLS TO FCN
!                                  NEEDED TO COMPUTE THE JACOBIAN
!                                  MATRIX.
!
      msum = MIN0(Ml+Mu+1,N)
!
!                                  INITIALIZE ITERATION COUNTER AND
!                                  MONITORS.
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
!                                  BEGINNING OF THE OUTER LOOP.
 200  CONTINUE
      jeval = .TRUE.
!                                  CALCULATE THE JACOBIAN MATRIX.
      iflag = 2
      CALL ZSPWB(N,X,Fvec,Fjac,Ldfjac,iflag,Ml,Mu,Epsfcn,Wa1,Wa2,   &
               & Par)
      Nfev = Nfev + msum
      IF ( iflag.LT.0 ) GOTO 1300
!                                  COMPUTE THE QR FACTORIZATION OF THE
!                                  JACOBIAN.
      CALL ZSPWG(N,N,Fjac,Ldfjac,.FALSE.,iwa,1,Wa1,Wa2,Wa3)
!                                  ON THE FIRST ITERATION AND IF MODE IS
!                                  1, SCALE ACCORDING TO THE NORMS OF   ZSPR
!                                  THE COLUMNS OF THE INITIAL JACOBIAN.
      IF ( iter.NE.1 ) GOTO 400
      IF ( Mode.EQ.2 ) GOTO 300
      DO j = 1 , N
         Diag(j) = Wa2(j)
         IF ( Wa2(j).EQ.zero ) Diag(j) = one
      ENDDO
 300  CONTINUE
!                                  ON THE FIRST ITERATION, CALCULATE THE
!                                  NORM OF THE SCALED X AND INITIALIZE
!                                  THE STEP BOUND DELTA.
      DO j = 1 , N
         Wa3(j) = Diag(j)*X(j)
      ENDDO
      xnorm = DNRM2(N,Wa3,1)
      delta = Factor*xnorm
      IF ( delta.EQ.zero ) delta = Factor
 400  CONTINUE
!                                  FORM (Q TRANSPOSE)*FVEC AND STORE IN
!                                  QTF.
      DO i = 1 , N
         Qtf(i) = Fvec(i)
      ENDDO
      DO j = 1 , N
         IF ( Fjac(j,j).EQ.zero ) GOTO 450
         sum = zero
         DO i = j , N
            sum = sum + Fjac(i,j)*Qtf(i)
         ENDDO
         temp = -sum/Fjac(j,j)
         DO i = j , N
            Qtf(i) = Qtf(i) + Fjac(i,j)*temp
         ENDDO
 450     CONTINUE
      ENDDO
!                                  COPY THE TRIANGULAR FACTOR OF THE QR
!                                  FACTORIZATION INTO R.
      sing = .FALSE.
      DO j = 1 , N
         l = j
         jm1 = j - 1
         IF ( jm1.LT.1 ) GOTO 500
         DO i = 1 , jm1
            R(l) = Fjac(i,j)
            l = l + N - i
         ENDDO
 500     CONTINUE
         R(l) = Wa1(j)
         IF ( Wa1(j).EQ.zero ) sing = .TRUE.
      ENDDO
!                                  ACCUMULATE THE ORTHOGONAL FACTOR IN
!                                  FJAC.
      CALL ZSPWF(N,N,Fjac,Ldfjac,Wa1)
!                                  RESCALE IF NECESSARY.
      IF ( Mode.EQ.2 ) GOTO 600
      DO j = 1 , N
         Diag(j) = DMAX1(Diag(j),Wa2(j))
      ENDDO
 600  CONTINUE
!                                  BEGINNING OF THE INNER LOOP.
 700  CONTINUE
!                                  IF REQUESTED, CALL FCN TO ENABLE
!                                  PRINTING OF ITERATES.
      IF ( Nprint.LE.0 ) GOTO 800
      iflag = 0
      IF ( iflag.LT.0 ) GOTO 1300
 800  CONTINUE
!                                  DETERMINE THE DIRECTION P.
      CALL ZSPWC(N,R,Lr,Diag,Qtf,delta,Wa1,Wa2,Wa3)
!                                  STORE THE DIRECTION P AND X + P.
!                                  CALCULATE THE NORM OF P.
      DO j = 1 , N
         Wa1(j) = -Wa1(j)
         Wa2(j) = X(j) + Wa1(j)
         Wa3(j) = Diag(j)*Wa1(j)
      ENDDO
      pnorm = DNRM2(N,Wa3,1)
!                                  ON THE FIRST ITERATION, ADJUST THE
!                                  INITIAL STEP BOUND.
      IF ( iter.EQ.1 ) delta = DMIN1(delta,pnorm)
!                                  EVALUATE THE FUNCTION AT X + P AND
!                                  CALCULATE ITS NORM.
      iflag = 1
      CALL FCN(Wa2,Wa4,N,Par)
      Nfev = Nfev + 1
      IF ( iflag.LT.0 ) GOTO 1300
      fnorm1 = DNRM2(N,Wa4,1)
!                                  COMPUTE THE SCALED ACTUAL REDUCTION.
      actred = -one
      IF ( fnorm1.LT.fnorm ) actred = one - (fnorm1/fnorm)**2
!                                  COMPUTE THE SCALED PREDICTED
!                                  REDUCTION.
      l = 1
      DO i = 1 , N
         sum = zero
         DO j = i , N
            sum = sum + R(l)*Wa1(j)
            l = l + 1
         ENDDO
         Wa3(i) = Qtf(i) + sum
      ENDDO
      temp = DNRM2(N,Wa3,1)
      prered = one
      IF ( temp.LT.fnorm ) prered = one - (temp/fnorm)**2
!                                  COMPUTE THE RATIO OF THE ACTUAL TO
!                                  THE PREDICTED REDUCTION.
      ratio = zero
      IF ( prered.GT.zero ) ratio = actred/prered
!                                  UPDATE THE STEP BOUND.
      IF ( ratio.GE.p1 ) GOTO 900
      ncsuc = 0
      ncfail = ncfail + 1
      delta = p5*delta
      GOTO 1000
 900  CONTINUE
      ncfail = 0
      ncsuc = ncsuc + 1
      IF ( ratio.GE.p5 .OR. ncsuc.GT.1 ) delta = DMAX1(delta,pnorm/p5)
      IF ( ABS(ratio-one).LE.p1 ) delta = pnorm/p5
 1000 CONTINUE
!                                  TEST FOR SUCCESSFUL ITERATION.
      IF ( ratio.LT.p0001 ) GOTO 1100
!                                  SUCCESSFUL ITERATION. UPDATE X, FVEC,
!                                  AND THEIR NORMS.
      DO j = 1 , N
         X(j) = Wa2(j)
         Wa2(j) = Diag(j)*X(j)
         Fvec(j) = Wa4(j)
      ENDDO
      xnorm = DNRM2(N,Wa2,1)
      fnorm = fnorm1
      iter = iter + 1
 1100 CONTINUE
!                                  DETERMINE THE PROGRESS OF THE
!                                  ITERATION.
      nslow1 = nslow1 + 1
      IF ( actred.GE.p001 ) nslow1 = 0
      IF ( jeval ) nslow2 = nslow2 + 1
      IF ( actred.GE.p1 ) nslow2 = 0
!                                  TEST FOR CONVERGENCE.
      IF ( delta.LE.Xtol*xnorm .OR. fnorm.EQ.zero ) Info = 1
      IF ( Info.NE.0 ) GOTO 1300
!                                  TESTS FOR TERMINATION AND STRINGENT
!                                  TOLERANCES.
      IF ( Nfev.GE.Maxfev ) Info = 2
      IF ( p1*DMAX1(p1*delta,pnorm).LE.epsmch*xnorm ) Info = 3
      IF ( nslow2.EQ.5 ) Info = 4
      IF ( nslow1.EQ.10 ) Info = 5
      IF ( Info.NE.0 ) GOTO 1300
!                                  CRITERION FOR RECALCULATING JACOBIAN
!                                  APPROXIMATION BY FORWARD DIFFERENCES.
      IF ( ncfail.EQ.2 ) GOTO 1200
!                                  CALCULATE THE RANK ONE MODIFICATION
!                                  TO THE JACOBIAN AND UPDATE QTF IF
!                                  NECESSARY.
      DO j = 1 , N
         sum = zero
         DO i = 1 , N
            sum = sum + Fjac(i,j)*Wa4(i)
         ENDDO
         Wa2(j) = (sum-Wa3(j))/pnorm
         Wa1(j) = Diag(j)*((Diag(j)*Wa1(j))/pnorm)
         IF ( ratio.GE.p0001 ) Qtf(j) = sum
      ENDDO
!                                  COMPUTE THE QR FACTORIZATION OF THE
!                                  UPDATED JACOBIAN.
      CALL ZSPWE(N,N,R,Lr,Wa1,Wa2,Wa3,sing)
      CALL ZSPWD(N,N,Fjac,Ldfjac,Wa2,Wa3)
      CALL ZSPWD(1,N,Qtf,1,Wa2,Wa3)
!                                  END OF THE INNER LOOP.
      jeval = .FALSE.
      GOTO 700
 1200 CONTINUE
!                                  END OF THE OUTER LOOP.
      GOTO 200
 1300 CONTINUE
!                                  TERMINATION, EITHER NORMAL OR USER
!                                  IMPOSED.
      IF ( iflag.LT.0 ) Info = iflag
      iflag = 0
      CONTINUE
      END
!*==ZSPWB.spg  processed by SPAG 6.72Dc at 15:10 on 27 Jul 2016
!DECK,ZSPWB
!   IMSL ROUTINE NAME   - ZSPWB
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDCFT5/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ZSPWB(N,X,Fvec,Fjac,Ldfjac,Iflag,Ml,Mu,Epsfcn,Wa1, &
                     & Wa2,Par)
      IMPLICIT NONE
!*--ZSPWB447
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER N , Ldfjac , Iflag , Ml , Mu
      REAL*8 X(N) , Fvec(N) , Fjac(Ldfjac,N) , Epsfcn , Wa1(N) , Wa2(N) &
           & , Par(*)
      EXTERNAL FCN
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER i , j , k , msum
      REAL*8 epsmch , eps , h , spmpar , temp , zero
      DATA spmpar/.710543D-14/
      DATA zero/0.0D0/
!                                  EPSMCH IS THE MACHINE PRECISION.
!                                  FIRST EXECUTABLE STATEMENT
      epsmch = spmpar
      eps = SQRT(DMAX1(Epsfcn,epsmch))
      msum = Ml + Mu + 1
      IF ( msum.LT.N ) GOTO 200
!                                  COMPUTATION OF DENSE APPROXIMATE
!                                  JACOBIAN.
      DO j = 1 , N
         temp = X(j)
         h = eps*ABS(temp)
         IF ( h.EQ.zero ) h = eps
         X(j) = temp + h
         CALL FCN(X,Wa1,N,Par)
         IF ( Iflag.LT.0 ) GOTO 100
         X(j) = temp
         DO i = 1 , N
            Fjac(i,j) = (Wa1(i)-Fvec(i))/h
         ENDDO
      ENDDO
 100  CONTINUE
      GOTO 400
 200  CONTINUE
!                                  COMPUTATION OF BANDED APPROXIMATE
!                                  JACOBIAN.
      DO k = 1 , msum
         DO j = k , N , msum
            Wa2(j) = X(j)
            h = eps*ABS(Wa2(j))
            IF ( h.EQ.zero ) h = eps
            X(j) = Wa2(j) + h
         ENDDO
         CALL FCN(X,Wa1,N,Par)
         IF ( Iflag.LT.0 ) GOTO 300
         DO j = k , N , msum
            X(j) = Wa2(j)
            h = eps*ABS(Wa2(j))
            IF ( h.EQ.zero ) h = eps
            DO i = 1 , N
               Fjac(i,j) = zero
               IF ( i.GE.j-Mu .AND. i.LE.j+Ml ) Fjac(i,j)               &
                  & = (Wa1(i)-Fvec(i))/h
            ENDDO
         ENDDO
      ENDDO
 300  CONTINUE
 400  CONTINUE
      CONTINUE
      END
!*==ZSPWC.spg  processed by SPAG 6.72Dc at 15:10 on 27 Jul 2016
!DECK,ZSPWC
!   IMSL ROUTINE NAME   - ZSPWC
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDCFT5/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - SINGLE/VBLA=SNRM2
!                       - DOUBLE/VBLA=DNRM2
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ZSPWC(N,R,Lr,Diag,Qtb,Delta,X,Wa1,Wa2)
      IMPLICIT NONE
!*--ZSPWC539
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER N , Lr
      REAL*8 R(Lr) , Diag(N) , Qtb(N) , Delta , X(N) , Wa1(N) , Wa2(N)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER i , jj , jp1 , j , k , l
      REAL*8 alpha , bnorm , epsmch , gnorm , one , qnorm , sgnorm ,    &
           & spmpar , sum , temp , zero
      REAL*8 DNRM2
      DATA spmpar/.710543D-14/
      DATA one , zero/1.0D0 , 0.0D0/
!                                  EPSMCH IS THE MACHINE PRECISION.
!                                  FIRST EXECUTABLE STATEMENT
      epsmch = spmpar
!                                  FIRST, CALCULATE THE GAUSS-NEWTON
!                                  DIRECTION.
      jj = (N*(N+1))/2 + 1
      DO k = 1 , N
         j = N - k + 1
         jp1 = j + 1
         jj = jj - k
         l = jj + 1
         sum = zero
         IF ( N.LT.jp1 ) GOTO 50
         DO i = jp1 , N
            sum = sum + R(l)*X(i)
            l = l + 1
         ENDDO
 50      CONTINUE
         temp = R(jj)
         IF ( temp.NE.zero ) GOTO 100
         l = j
         DO i = 1 , j
            temp = DMAX1(temp,ABS(R(l)))
            l = l + N - i
         ENDDO
         temp = epsmch*temp
         IF ( temp.EQ.zero ) temp = epsmch
 100     CONTINUE
         X(j) = (Qtb(j)-sum)/temp
      ENDDO
!                                  TEST WHETHER THE GAUSS-NEWTON
!                                  DIRECTION IS ACCEPTABLE.
      DO j = 1 , N
         Wa1(j) = zero
         Wa2(j) = Diag(j)*X(j)
      ENDDO
      qnorm = DNRM2(N,Wa2,1)
      IF ( qnorm.LE.Delta ) GOTO 300
!                                  THE GAUSS-NEWTON DIRECTION IS NOT
!                                  ACCEPTABLE. NEXT, CALCULATE THE
!                                  SCALED GRADIENT DIRECTION.
      l = 1
      DO j = 1 , N
         temp = Qtb(j)
         DO i = j , N
            Wa1(i) = Wa1(i) + R(l)*temp
            l = l + 1
         ENDDO
         Wa1(j) = Wa1(j)/Diag(j)
      ENDDO
!                                  CALCULATE THE NORM OF THE SCALED
!                                  GRADIENT AND TEST FOR THE SPECIAL
!                                  CASE IN WHICH THE SCALED GRADIENT IS
!                                  ZERO.
      gnorm = DNRM2(N,Wa1,1)
      sgnorm = zero
      alpha = Delta/qnorm
      IF ( gnorm.EQ.zero ) GOTO 200
!                                  CALCULATE THE POINT ALONG THE SCALED
!                                  GRADIENT AT WHICH THE QUADRATIC IS
!                                  MINIMIZED.                           ZSPT
      DO j = 1 , N
         Wa1(j) = (Wa1(j)/gnorm)/Diag(j)
      ENDDO
      l = 1
      DO j = 1 , N
         sum = zero
         DO i = j , N
            sum = sum + R(l)*Wa1(i)
            l = l + 1
         ENDDO
         Wa2(j) = sum
      ENDDO
      temp = DNRM2(N,Wa2,1)
      sgnorm = (gnorm/temp)/temp
!                                  TEST WHETHER THE SCALED GRADIENT
!                                  DIRECTION IS ACCEPTABLE.
      alpha = zero
      IF ( sgnorm.GE.Delta ) GOTO 200
!                                  THE SCALED GRADIENT DIRECTION IS NOT
!                                  ACCEPTABLE. FINALLY, CALCULATE THE
!                                  POINT ALONG THE DOGLEG AT WHICH THE
!                                  QUADRATIC IS MINIMIZED.
      bnorm = DNRM2(N,Qtb,1)
      temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/Delta)
      temp = temp - (Delta/qnorm)*(sgnorm/Delta)                        &
           & **2 + SQRT((temp-(Delta/qnorm))**2+(one-(Delta/qnorm)**2)  &
           & *(one-(sgnorm/Delta)**2))
      alpha = ((Delta/qnorm)*(one-(sgnorm/Delta)**2))/temp
 200  CONTINUE
!                                  FORM APPROPRIATE CONVEX COMBINATION
!                                  OF THE GAUSS-NEWTON DIRECTION AND THE
!                                  SCALED GRADIENT DIRECTION.
      temp = (one-alpha)*DMIN1(sgnorm,Delta)
      DO j = 1 , N
         X(j) = temp*Wa1(j) + alpha*X(j)
      ENDDO
 300  CONTINUE
      CONTINUE
      END
!*==ZSPWD.spg  processed by SPAG 6.72Dc at 15:10 on 27 Jul 2016
!DECK,ZSPWD
!   IMSL ROUTINE NAME   - ZSPWD
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDCFT5/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ZSPWD(M,N,A,Lda,V,W)
      IMPLICIT NONE
!*--ZSPWD681
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER M , N , Lda
      REAL*8 A(Lda,N) , V(N) , W(N)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER i , j , nm1 , nmj
      REAL*8 temp1 , one , temp2 , temp
      DATA one/1.0D0/
!                                  APPLY THE FIRST SET OF GIVENS
!                                  ROTATIONS TO A.
!                                  FIRST EXECUTABLE STATEMENT
      nm1 = N - 1
      IF ( nm1.LT.1 ) GOTO 100
      DO nmj = 1 , nm1
         j = N - nmj
         IF ( ABS(V(j)).GT.one ) temp1 = one/V(j)
         IF ( ABS(V(j)).GT.one ) temp2 = SQRT(one-temp1**2)
         IF ( ABS(V(j)).LE.one ) temp2 = V(j)
         IF ( ABS(V(j)).LE.one ) temp1 = SQRT(one-temp2**2)
         DO i = 1 , M
            temp = temp1*A(i,j) - temp2*A(i,N)
            A(i,N) = temp2*A(i,j) + temp1*A(i,N)
            A(i,j) = temp
         ENDDO
      ENDDO
!                                  APPLY THE SECOND SET OF GIVENS
!                                  ROTATIONS TO A.
      DO j = 1 , nm1
         IF ( ABS(W(j)).GT.one ) temp1 = one/W(j)
         IF ( ABS(W(j)).GT.one ) temp2 = SQRT(one-temp1**2)
         IF ( ABS(W(j)).LE.one ) temp2 = W(j)
         IF ( ABS(W(j)).LE.one ) temp1 = SQRT(one-temp2**2)
         DO i = 1 , M
            temp = temp1*A(i,j) + temp2*A(i,N)
            A(i,N) = -temp2*A(i,j) + temp1*A(i,N)
            A(i,j) = temp
         ENDDO
      ENDDO
 100  CONTINUE
      CONTINUE
      END
!*==ZSPWE.spg  processed by SPAG 6.72Dc at 15:10 on 27 Jul 2016
!DECK,ZSPWE
!   IMSL ROUTINE NAME   - ZSPWE
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDCFT5/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ZSPWE(M,N,S,Ls,U,V,W,Sing)
      IMPLICIT NONE
!*--ZSPWE753
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER M , N , Ls
      REAL*8 S(Ls) , U(M) , V(N) , W(M)
      LOGICAL Sing
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER i , jj , j , l , nm1 , nmj
      REAL*8 temp1 , temp2 , giant , one , p25 , p5 , temp3 , temp4 ,   &
           & tau , temp , zero
!     DATA               GIANT /.12650140831E+323/
      DATA giant/.1D+38/
      DATA one , p5 , p25 , zero/1.0D0 , 5.0D-1 , 2.5D-1 , 0.0D0/
!                                  INITIALIZE THE DIAGONAL ELEMENT
!                                  POINTER.
!                                  FIRST EXECUTABLE STATEMENT
      jj = (N*(2*M-N+1))/2 - (M-N)
!                                  MOVE THE NONTRIVIAL PART OF THE LAST
!                                  COLUMN OF S INTO W.
      l = jj
      DO i = N , M
         W(i) = S(l)
         l = l + 1
      ENDDO
!                                  ROTATE THE VECTOR V INTO A MULTIPLE
!                                  OF THE N-TH UNIT VECTOR IN SUCH A WAY
!                                  THAT A SPIKE IS INTRODUCED INTO W.
      nm1 = N - 1
      IF ( nm1.LT.1 ) GOTO 200
      DO nmj = 1 , nm1
         j = N - nmj
         jj = jj - (M-j+1)
         W(j) = zero
         IF ( V(j).EQ.zero ) GOTO 150
!                                  DETERMINE A GIVENS ROTATION WHICH
!                                  ELIMINATES THE J-TH ELEMENT OF V.
         IF ( ABS(V(N)).GE.ABS(V(j)) ) GOTO 50
         temp2 = V(N)/V(j)
         temp3 = p5/SQRT(p25+p25*temp2**2)
         temp1 = temp3*temp2
         tau = one
         IF ( ABS(temp1)*giant.GT.one ) tau = one/temp1
         GOTO 100
 50      CONTINUE
         temp4 = V(j)/V(N)
         temp1 = p5/SQRT(p25+p25*temp4**2)
         temp3 = temp1*temp4
         tau = temp3
 100     CONTINUE
!                                  APPLY THE TRANSFORMATION TO V AND
!                                  STORE THE INFORMATION NECESSARY TO
!                                  RECOVER THE GIVENS ROTATION.
         V(N) = temp3*V(j) + temp1*V(N)
         V(j) = tau
!                                  APPLY THE TRANSFORMATION TO S AND
!                                  EXTEND THE SPIKE IN W.
         l = jj
         DO i = j , M
            temp = temp1*S(l) - temp3*W(i)
            W(i) = temp3*S(l) + temp1*W(i)
            S(l) = temp
            l = l + 1
         ENDDO
 150     CONTINUE
      ENDDO
 200  CONTINUE
!                                  ADD THE SPIKE FROM THE RANK 1 UPDATE
!                                  TO W.
      DO i = 1 , M
         W(i) = W(i) + V(N)*U(i)
      ENDDO
!                                  ELIMINATE THE SPIKE.
      Sing = .FALSE.
      IF ( nm1.LT.1 ) GOTO 400
      DO j = 1 , nm1
         IF ( W(j).EQ.zero ) GOTO 350
!                                  DETERMINE A GIVENS ROTATION WHICH
!                                  ELIMINATES THE J-TH ELEMENT OF THE
!                                  SPIKE.
         IF ( ABS(S(jj)).GE.ABS(W(j)) ) GOTO 250
         temp2 = S(jj)/W(j)
         temp3 = p5/SQRT(p25+p25*temp2**2)
         temp1 = temp3*temp2
         tau = one
         IF ( ABS(temp1)*giant.GT.one ) tau = one/temp1
         GOTO 300
 250     CONTINUE
         temp4 = W(j)/S(jj)
         temp1 = p5/SQRT(p25+p25*temp4**2)
         temp3 = temp1*temp4
         tau = temp3
 300     CONTINUE
!                                  APPLY THE TRANSFORMATION TO S AND
!                                  REDUCE THE SPIKE IN W.
         l = jj
         DO i = j , M
            temp = temp1*S(l) + temp3*W(i)
            W(i) = -temp3*S(l) + temp1*W(i)
            S(l) = temp
            l = l + 1
         ENDDO
!                                  STORE THE INFORMATION NECESSARY TO
!                                  RECOVER THE GIVENS ROTATION.
         W(j) = tau
 350     CONTINUE
!                                  TEST FOR ZERO DIAGONAL ELEMENTS IN
!                                  THE OUTPUT S.
         IF ( S(jj).EQ.zero ) Sing = .TRUE.
         jj = jj + (M-j+1)
      ENDDO
 400  CONTINUE
!                                  MOVE W BACK INTO THE LAST COLUMN OF
!                                  THE OUTPUT S.
      l = jj
      DO i = N , M
         S(l) = W(i)
         l = l + 1
      ENDDO
      IF ( S(jj).EQ.zero ) Sing = .TRUE.
      CONTINUE
      END
!*==ZSPWF.spg  processed by SPAG 6.72Dc at 15:10 on 27 Jul 2016
!DECK,ZSPWF
!   IMSL ROUTINE NAME   - ZSPWF
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDCFT5/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ZSPWF(M,N,Q,Ldq,Wa)
      IMPLICIT NONE
!*--ZSPWF904
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER M , N , Ldq
      REAL*8 Q(Ldq,M) , Wa(M)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER i , jm1 , j , k , l , minmn , np1
      REAL*8 one , sum , temp , zero
      DATA one , zero/1.0D0 , 0.0D0/
!                                  ZERO OUT UPPER TRIANGLE OF Q IN THE
!                                  FIRST MIN(M,N) COLUMNS.
!                                  FIRST EXECUTABLE STATEMENT
      minmn = MIN0(M,N)
      IF ( minmn.LT.2 ) GOTO 100
      DO j = 2 , minmn
         jm1 = j - 1
         DO i = 1 , jm1
            Q(i,j) = zero
         ENDDO
      ENDDO
 100  CONTINUE
!                                  INITIALIZE REMAINING COLUMNS TO THOSE
!                                  OF THE IDENTITY MATRIX.
      np1 = N + 1
      IF ( M.LT.np1 ) GOTO 200
      DO j = np1 , M
         DO i = 1 , M
            Q(i,j) = zero
         ENDDO
         Q(j,j) = one
      ENDDO
 200  CONTINUE
!                                  ACCUMULATE Q FROM ITS FACTORED FORM.
      DO l = 1 , minmn
         k = minmn - l + 1
         DO i = k , M
            Wa(i) = Q(i,k)
            Q(i,k) = zero
         ENDDO
         Q(k,k) = one
         IF ( Wa(k).EQ.zero ) GOTO 250
         DO j = k , M
            sum = zero
            DO i = k , M
               sum = sum + Q(i,j)*Wa(i)
            ENDDO
            temp = sum/Wa(k)
            DO i = k , M
               Q(i,j) = Q(i,j) - temp*Wa(i)
            ENDDO
         ENDDO
 250     CONTINUE
      ENDDO
      CONTINUE
      END
!*==ZSPWG.spg  processed by SPAG 6.72Dc at 15:10 on 27 Jul 2016
!DECK,ZSPWG
!   IMSL ROUTINE NAME   - ZSPWG
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDCFT5/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1982
!
!   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ZSPOW
!
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!                       - SINGLE/H36,H48,H60
!
!   REQD. IMSL ROUTINES - SINGLE/VBLA=SNRM2
!                       - DOUBLE/VBLA=DNRM2
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
!
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ZSPWG(M,N,A,Lda,Pivot,Ipvt,Lipvt,Rdiag,Acnorm,Wa)
      IMPLICIT NONE
!*--ZSPWG990
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER M , N , Lda , Lipvt , Ipvt(Lipvt)
      REAL*8 A(Lda,N) , Rdiag(N) , Acnorm(N) , Wa(N)
      LOGICAL Pivot
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER i , jp1 , j , kmax , k , minmn
      REAL*8 ajnorm , epsmch , one , p05 , spmpar , sum , temp , zero
      REAL*8 DNRM2
      DATA spmpar/.710543D-14/
      DATA one , p05 , zero/1.0D0 , 5.0D-2 , 0.0D0/
!                                  EPSMCH IS THE MACHINE PRECISION.
!                                  FIRST EXECUTABLE STATEMENT
      epsmch = spmpar
!                                  COMPUTE THE INITIAL COLUMN NORMS AND
!                                  INITIALIZE SEVERAL ARRAYS.
      DO j = 1 , N
         Acnorm(j) = DNRM2(M,A(1,j),1)
         Rdiag(j) = Acnorm(j)
         Wa(j) = Rdiag(j)
         IF ( Pivot ) Ipvt(j) = j
      ENDDO
!                                  REDUCE A TO R WITH HOUSEHOLDER
!                                  TRANSFORMATIONS.
      minmn = MIN0(M,N)
      DO j = 1 , minmn
         IF ( .NOT.Pivot ) GOTO 50
!                                  BRING THE COLUMN OF LARGEST NORM INTO
!                                  THE PIVOT POSITION.
         kmax = j
         DO k = j , N
            IF ( Rdiag(k).GT.Rdiag(kmax) ) kmax = k
         ENDDO
         IF ( kmax.EQ.j ) GOTO 50
         DO i = 1 , M
            temp = A(i,j)
            A(i,j) = A(i,kmax)
            A(i,kmax) = temp
         ENDDO
         Rdiag(kmax) = Rdiag(j)
         Wa(kmax) = Wa(j)
         k = Ipvt(j)
         Ipvt(j) = Ipvt(kmax)
         Ipvt(kmax) = k
 50      CONTINUE
!                                  COMPUTE THE HOUSEHOLDER
!                                  TRANSFORMATION TO REDUCE THE J-TH
!                                  COLUMN OF A TO A MULTIPLE OF THE J-TH
!                                  UNIT VECTOR.
         ajnorm = DNRM2(M-j+1,A(j,j),1)
         IF ( ajnorm.EQ.zero ) GOTO 100
         IF ( A(j,j).LT.zero ) ajnorm = -ajnorm
         DO i = j , M
            A(i,j) = A(i,j)/ajnorm
         ENDDO
         A(j,j) = A(j,j) + one
!                                  APPLY THE TRANSFORMATION TO THE
!                                  REMAINING COLUMNS AND UPDATE THE
!                                  NORMS.
         jp1 = j + 1
         IF ( N.LT.jp1 ) GOTO 100
         DO k = jp1 , N
            sum = zero
            DO i = j , M
               sum = sum + A(i,j)*A(i,k)
            ENDDO
            temp = sum/A(j,j)
            DO i = j , M
               A(i,k) = A(i,k) - temp*A(i,j)
            ENDDO
            IF ( .NOT.Pivot .OR. Rdiag(k).EQ.zero ) GOTO 60
            temp = A(j,k)/Rdiag(k)
            Rdiag(k) = Rdiag(k)*SQRT(DMAX1(zero,one-temp**2))
            IF ( p05*(Rdiag(k)/Wa(k))**2.GT.epsmch ) GOTO 60
            Rdiag(k) = DNRM2(M-j,A(jp1,k),1)
            Wa(k) = Rdiag(k)
 60         CONTINUE
         ENDDO
 100     CONTINUE
         Rdiag(j) = -ajnorm
      ENDDO
      CONTINUE
      END
!*==UERTST.spg  processed by SPAG 6.72Dc at 15:10 on 27 Jul 2016
 
 
 
 
!   IMSL ROUTINE NAME   - UERTST
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDC/SINGLE
!
!   LATEST REVISION     - NOVEMBER 1, 1979
!
!   PURPOSE             - PRINT A MESSAGE REFLECTING AN ERROR CONDITION
!
!   USAGE               - CALL UERTST (IER,NAME)
!
!   ARGUMENTS    IER    - ERROR PARAMETER. (INPUT)
!                           IER = I+J WHERE
!                             I = 128 IMPLIES TERMINAL ERROR,
!                             I =  64 IMPLIES WARNING WITH FIX, AND
!                             I =  32 IMPLIES WARNING.
!                             J = ERROR CODE RELEVANT TO CALLING
!                                 ROUTINE.
!                NAME   - A SIX CHARACTER LITERAL STRING GIVING THE
!                           NAME OF THE CALLING ROUTINE. (INPUT)
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - UGETIO
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS      THE ERROR MESSAGE PRODUCED BY UERTST IS WRITTEN
!                ONTO THE STANDARD OUTPUT UNIT. THE OUTPUT UNIT
!                NUMBER CAN BE DETERMINED BY CALLING UGETIO AS
!                FOLLOWS..   CALL UGETIO(1,NIN,NOUT).
!                THE OUTPUT UNIT NUMBER CAN BE CHANGED BY CALLING
!                UGETIO AS FOLLOWS..
!                                NIN = 0
!                                NOUT = NEW OUTPUT UNIT NUMBER
!                                CALL UGETIO(3,NIN,NOUT)
!                SEE THE UGETIO DOCUMENT FOR MORE DETAILS.
!
!
!
!-----------------------------------------------------------------------
!
      SUBROUTINE UERTST(Ier,Name)
      IMPLICIT NONE
!*--UERTST1125
!*** Start of declarations inserted by SPAG
      INTEGER i , ieqdf , Ier , iounit , level , levold , nin
!*** End of declarations inserted by SPAG
!                                  SPECIFICATIONS FOR ARGUMENTS
      CHARACTER*1 ieq
      CHARACTER*2 Name(3)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      CHARACTER*2 namset(3) , nameq(3)
!     DATA               NAMSET/2HUE,2HRS,2HET/
!     DATA               NAMEQ/2H  ,2H  ,2H  /
!                                  FIRST EXECUTABLE STATEMENT
      DATA level/4/ , ieqdf/0/
!     DATA		 IEQ/1H=/
      ieq(1:1) = '='
      namset(1) = 'UE'
      namset(2) = 'RS'
      namset(3) = 'ET'
      nameq(1) = '  '
      nameq(2) = '  '
      nameq(3) = '  '
      IF ( Ier.GT.999 ) GOTO 400
      IF ( Ier.LT.-32 ) GOTO 600
      IF ( Ier.LE.128 ) GOTO 100
      IF ( level.LT.1 ) GOTO 500
!                                  PRINT TERMINAL MESSAGE
      CALL UGETIO(1,nin,iounit)
      IF ( ieqdf.EQ.1 ) WRITE (iounit,99001) Ier , nameq , ieq , Name
      IF ( ieqdf.EQ.0 ) WRITE (iounit,99001) Ier , Name
      GOTO 500
 100  IF ( Ier.LE.64 ) GOTO 200
      IF ( level.LT.2 ) GOTO 500
!                                  PRINT WARNING WITH FIX MESSAGE
      CALL UGETIO(1,nin,iounit)
      IF ( ieqdf.EQ.1 ) WRITE (iounit,99002) Ier , nameq , ieq , Name
      IF ( ieqdf.EQ.0 ) WRITE (iounit,99002) Ier , Name
      GOTO 500
 200  IF ( Ier.LE.32 ) GOTO 300
!                                  PRINT WARNING MESSAGE
      IF ( level.LT.3 ) GOTO 500
      CALL UGETIO(1,nin,iounit)
      IF ( ieqdf.EQ.1 ) WRITE (iounit,99003) Ier , nameq , ieq , Name
      IF ( ieqdf.EQ.0 ) WRITE (iounit,99003) Ier , Name
      GOTO 500
 300  CONTINUE
!                                  CHECK FOR UERSET CALL
      DO i = 1 , 3
         IF ( Name(i).NE.namset(i) ) GOTO 400
      ENDDO
      levold = level
      level = Ier
      Ier = levold
      IF ( level.LT.0 ) level = 4
      IF ( level.GT.4 ) level = 4
      GOTO 500
 400  CONTINUE
      IF ( level.LT.4 ) GOTO 500
!                                  PRINT NON-DEFINED MESSAGE
      CALL UGETIO(1,nin,iounit)
      IF ( ieqdf.EQ.1 ) WRITE (iounit,99004) Ier , nameq , ieq , Name
      IF ( ieqdf.EQ.0 ) WRITE (iounit,99004) Ier , Name
 500  ieqdf = 0
      RETURN
!                                  SAVE P FOR P = R CASE
!                                    P IS THE PAGE NAME
!                                    R IS THE ROUTINE NAME
 600  ieqdf = 1
      DO i = 1 , 3
         nameq(i) = Name(i)
      ENDDO
      CONTINUE
99001 FORMAT (' *** TERMINAL ERROR',10X,'(IER = ',I3,                   &
             &') FROM IMSL ROUTINE ',3A2,A1,3A2)
99002 FORMAT (' *** WARNING WITH FIX ERROR  (IER = ',I3,                &
             &') FROM IMSL ROUTINE ',3A2,A1,3A2)
99003 FORMAT (' *** WARNING ERROR',11X,'(IER = ',I3,                    &
             &') FROM IMSL ROUTINE ',3A2,A1,3A2)
99004 FORMAT (' *** UNDEFINED ERROR',9X,'(IER = ',I5,                   &
             &') FROM IMSL ROUTINE ',3A2,A1,3A2)
      END
!*==UGETIO.spg  processed by SPAG 6.72Dc at 15:10 on 27 Jul 2016
 
 
!   IMSL ROUTINE NAME   - UGETIO
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDC/SINGLE
!
!   LATEST REVISION     - JUNE 1, 1981
!
!   PURPOSE             - TO RETRIEVE CURRENT VALUES AND TO SET NEW
!                           VALUES FOR INPUT AND OUTPUT UNIT
!                           IDENTIFIERS.
!
!   USAGE               - CALL UGETIO(IOPT,NIN,NOUT)
!
!   ARGUMENTS    IOPT   - OPTION PARAMETER. (INPUT)
!                           IF IOPT=1, THE CURRENT INPUT AND OUTPUT
!                           UNIT IDENTIFIER VALUES ARE RETURNED IN NIN
!                           AND NOUT, RESPECTIVELY.
!                           IF IOPT=2, THE INTERNAL VALUE OF NIN IS
!                           RESET FOR SUBSEQUENT USE.
!                           IF IOPT=3, THE INTERNAL VALUE OF NOUT IS
!                           RESET FOR SUBSEQUENT USE.
!                NIN    - INPUT UNIT IDENTIFIER.
!                           OUTPUT IF IOPT=1, INPUT IF IOPT=2.
!                NOUT   - OUTPUT UNIT IDENTIFIER.
!                           OUTPUT IF IOPT=1, INPUT IF IOPT=3.
!
!   PRECISION/HARDWARE  - SINGLE/ALL
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!   REMARKS      EACH IMSL ROUTINE THAT PERFORMS INPUT AND/OR OUTPUT
!                OPERATIONS CALLS UGETIO TO OBTAIN THE CURRENT UNIT
!                IDENTIFIER VALUES. IF UGETIO IS CALLED WITH IOPT=2 OR
!                IOPT=3, NEW UNIT IDENTIFIER VALUES ARE ESTABLISHED.
!                SUBSEQUENT INPUT/OUTPUT IS PERFORMED ON THE NEW UNITS.
!
!
!
!-----------------------------------------------------------------------
!
      SUBROUTINE UGETIO(Iopt,Nin,Nout)
      IMPLICIT NONE
!*--UGETIO1255
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER Iopt , Nin , Nout
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER nind , noutd
!cc   the following line is replaced by the next
!     DATA               NIND/1/,NOUTD/2/
      DATA nind/5/ , noutd/6/
!                                  FIRST EXECUTABLE STATEMENT
      IF ( Iopt.EQ.3 ) GOTO 200
      IF ( Iopt.EQ.2 ) GOTO 100
      IF ( Iopt.NE.1 ) GOTO 300
      Nin = nind
      Nout = noutd
      GOTO 300
 100  nind = Nin
      GOTO 300
 200  noutd = Nout
 300  CONTINUE
      END
!*==DNRM2.spg  processed by SPAG 6.72Dc at 15:10 on 27 Jul 2016
 
 
 
!   IMSL ROUTINE NAME   - VBLA=DNRM2
!
!-----------------------------------------------------------------------
!
!   COMPUTER            - CDC/DOUBLE
!
!   LATEST REVISION     - JANUARY 1, 1978
!
!   PURPOSE             - COMPUTE THE EUCLIDEAN LENGTH OR L2 NORM
!                           OF A DOUBLE PRECISION VECTOR
!
!   USAGE               - FUNCTION DNRM2 (N,DX,INCX)
!
!   ARGUMENTS    DNRM2  - DOUBLE PRECISION SQUARE ROOT OF THE SUM FROM
!                           I=1 TO N OF X(I)**2. (OUTPUT)
!                           X(I) REFERS TO A SPECIFIC ELEMENT OF DX.
!                           SEE INCX ARGUMENT DESCRIPTION.
!                N      - LENGTH OF VECTOR X. (INPUT)
!                DX     - DOUBLE PRECISION VECTOR OF LENGTH N*INCX.
!                           (INPUT)
!                INCX   - DISPLACEMENT BETWEEN ELEMENTS OF DX. (INPUT)
!                           X(I) IS DEFINED TO BE DX(1+(I-1)*INCX).
!                           INCX MUST BE GREATER THAN ZERO.
!
!   PRECISION/HARDWARE  - DOUBLE/ALL
!
!   REQD. IMSL ROUTINES - NONE REQUIRED
!
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!
!
!
!-----------------------------------------------------------------------
!
      DOUBLE PRECISION FUNCTION DNRM2(N,Dx,Incx)
      IMPLICIT NONE
!*--DNRM21317
!
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER N , Incx
      DOUBLE PRECISION Dx(*)
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER i , j , next , nn
      DOUBLE PRECISION cutlo , cuthi , sum , xmax , zero , one , hitest
      DATA zero , one/0.0D0 , 1.0D0/
      DATA cutlo , cuthi/8.232D-11 , 1.304D19/
!                                  FIRST EXECUTABLE STATEMENT
      IF ( N.GT.0 ) GOTO 100
      DNRM2 = zero
      GOTO 1300
!
 100  ASSIGN 300 TO next
      sum = zero
      nn = N*Incx
!                                  BEGIN MAIN LOOP
      i = 1
 200  GOTO next
 300  IF ( DABS(Dx(i)).GT.cutlo ) GOTO 1100
      ASSIGN 400 TO next
      xmax = zero
!                                  PHASE 1. SUM IS ZERO
 400  IF ( Dx(i).EQ.zero ) GOTO 1200
      IF ( DABS(Dx(i)).GT.cutlo ) GOTO 1100
!                                  PREPARE FOR PHASE 2.
      ASSIGN 700 TO next
      GOTO 600
!                                  PREPARE FOR PHASE 4.
 500  i = j
      ASSIGN 800 TO next
      sum = (sum/Dx(i))/Dx(i)
 600  xmax = DABS(Dx(i))
      GOTO 900
!                                  PHASE 2. SUM IS SMALL. SCALE TO
!                                    AVOID DESTRUCTIVE UNDERFLOW.
 700  IF ( DABS(Dx(i)).GT.cutlo ) GOTO 1000
!                                  COMMON CODE FOR PHASES 2 AND 4. IN
!                                    PHASE 4 SUM IS LARGE. SCALE TO
!                                    AVOID OVERFLOW.
 800  IF ( DABS(Dx(i)).LE.xmax ) GOTO 900
      sum = one + sum*(xmax/Dx(i))**2
      xmax = DABS(Dx(i))
      GOTO 1200
!
 900  sum = sum + (Dx(i)/xmax)**2
      GOTO 1200
!                                  PREPARE FOR PHASE 3.
 1000 sum = (sum*xmax)*xmax
!                                  FOR REAL OR D.P. SET HITEST =
!                                    CUTHI/N FOR COMPLEX SET HITEST =
!                                    CUTHI/(2*N)
 1100 hitest = cuthi/FLOAT(N)
!                                  PHASE 3. SUM IS MID-RANGE. NO
!                                    SCALING.
      DO j = i , nn , Incx
         IF ( DABS(Dx(j)).GE.hitest ) GOTO 500
         sum = sum + Dx(j)**2
      ENDDO
      DNRM2 = DSQRT(sum)
      GOTO 1300
!
 1200 CONTINUE
      i = i + Incx
      IF ( i.LE.nn ) GOTO 200
!                                  END OF MAIN LOOP. COMPUTE SQUARE
!                                    ROOT AND ADJUST FOR SCALING.
      DNRM2 = xmax*DSQRT(sum)
 1300 CONTINUE
      CONTINUE
      END
 
