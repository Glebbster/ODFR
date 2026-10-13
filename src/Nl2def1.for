c-begin-text nl2defu
*-module d1mach
!                            * * * * * * * * *
!                            *    D1MACH     *
!                            * * * * * * * * *


cAB      real*8 function D1MACH( No )
cAB 
cAB      real*8        EPS, RLEAST, RLARGE
cAB      common / MACHIN / EPS, RLEAST, RLARGE, NSD
cAB
cAB      SAVE
cAB
c     Machine-dependent constants. See also X2MDC1 subroutine.
cAB      call X2MDC1( EPS, RLEAST, RLARGE, NSD )
cAB
cAB      D1MACH  =  0.0
cAB      if ( No .EQ. 1 ) then
cAB         D1MACH  =  RLEAST              ! The least representable pos. num.
cAB      else if ( No .EQ. 2 ) then
cAB         D1MACH  =  RLARGE              ! The largest one.
cAB      else if ( No .EQ. 4 ) then
cAB         D1MACH  =  EPS                 ! Machine epsilon.
cAB      end if
cAB      return
cAB      end
*-end-module d1mach
*-module i1mach
!                            * * * * * * * * *
!                            *    I1MACH     *
!                            * * * * * * * * *

      integer function I1MACH( No )

c     I/O channel numbers.

      if ( No .EQ. 1 ) then
         I1MACH  =  5                  ! Input.
      else if ( No .EQ. 2 ) then
         I1MACH  =  6                  ! Output.
      else if ( No .EQ. 3) then
         I1MACH  =  6
      else if ( No .EQ. 4 ) then
         I1MACH  =  6
      end if
      return
      end
*-end-module i1mach
C*BEGIN-TEXT x2mdc1
C                          * * * * * * * * * *
C                          *     X2MDC1      *
C                          * * * * * * * * * *
C
C     USER'S GUIDE TO X2MDC1: COMPUTATION OF MACHINE DEPENDENT
C                             CONSTANTS.
C
C     EPS ------ THE SMALLEST (IN MAGNITUDE) FLOATING POINT NUMBER
C                WHICH, WHEN ADDED TO THE FLOATING POINT NUMBER 1.0,
C                PRODUCES A FLOATING POINT RESULT STILL DIFFERENT FROM
C                1.0 (MACHINE ACCURACY).
C     RLARGE --- LARGEST POSITIVE REAL NUMBER SUCH THAT BOTH "RLARGE"
C                AND "-RLARGE" CAN BE HELD IN THE MACHINE.
C     RLEAST --- SMALLEST POSITIVE REAL FLOATING POINT NUMBER REPRESEN-
C                TABLE ON THE COMPUTER.
C     NSD ------ NUMBER OF SIGNIFICANT DIGITS (DECIMAL MANTISSA LENGTH).
C
C     WRITTEN BY V.V.VOLKOV                           REVISED 09.08.89
C     ----------------------------------------------------------------
C
      SUBROUTINE X2MDC1( EPS, RLEAST, RLARGE, NSD )
C
      real*8            Eps, rLeast, rLarge, ONE, TEN, rMult, Temp, Sum

      DATA              ONE /1.0d+000/,  TEN /1.0d+001/
C
      EPS    =  ONE
      SUM    =  EPS + ONE
      RMULT  =  0.1d+000
      NSD    =  0
C
C   PERFORM DOUBLE PROTECTED SCHEME FOR "EPS":
C
      DO WHILE (RMULT .LT. ONE)
         DO WHILE (SUM .GT. ONE)
            TEMP = EPS*RMULT
            IF (TEMP .GE. EPS) GO TO 10
            EPS  = TEMP
            SUM  = EPS + ONE
         END DO
         EPS   = EPS/RMULT
         SUM   = EPS + ONE
         NSD   = NSD + 1
         RMULT = ONE - TEN**(-NSD)
      END DO
C
C     COMPUTE THE DECIMAL MANTISSA LENGTH:
 10   NSD = INT( DABS( DLOG10( EPS ) ) + 0.5 )
C
C   FIND THE APPROXIMATE TO THE LEAST POSITIVE REAL NUMBER:
C
CAI   MSDOS FORTRAN (1 LINE):
      RLEAST = 1.0d-306
CAI   SINTRANIII FORTRAN (15 LINES):
CS      RLEAST = ONE
CS      IEXP = 0
CSC     DIRECT MOVE:
CS      DO WHILE (RLEAST .GT. 0.0)
CS          IEXP = IEXP + 1
CS          RLEAST = TEN**(-IEXP)
CS      END DO
CSC     REVERSE MOVE:
CS      DO WHILE (RLEAST .EQ. 0.0)
CS          IEXP = IEXP - 1
CS          RLEAST = TEN**(-IEXP)
CS      END DO
CSC     INCREASE THE PRODUCED NUMBER JUST IN CASE:
CS      IEXP = IEXP - 2
CS      RLEAST = TEN**(-IEXP)
C
C   FIND THE APPROXIMATE TO THE LARGEST REAL NUMBER:
C
CAI   MSDOS FORTRAN (1 LINE):
      RLARGE = 1.0d+306
CAI   SINTRANIII FORTRAN (1 LINE):
CS      RLARGE = TEN**(IEXP)
C
      RETURN
      END
C*END-TEXT x2mdc1
*-module dfault
!                            * * * * * * * * *
!                            *    DFAULT     *
!                            * * * * * * * * *

      subroutine DFAULT( iv, v )

c     Default values used by NL2SOL.

C     ARRAY ARGUMENTS
      include 'dimnls.inc'
      include 'constNL2SOL.inc'
C
C     LOCAL SCALARS
      real*8            machep, mepcrt, one, sqteps, three
     *,                 stepJac, stepBound
      integer
     +   AFCTOL,COSMIN,COVPRT,COVREQ,D0INIT,DECFAC,DELTA0,DFAC,
     +   DINIT,DLTFDC,DLTFDJ,DTYPE,EPSLON,FUZZ,INCFAC,INITS,JTINIT,
     +   LMAX0,MXFCAL,MXITER,OUTLEV,PARPRT,PHMNFC,PHMXFC,PRUNIT,
     +   RDFCMN,RDFCMX,RFCTOL,RLIMIT,SOLPRT,STATPR,TUNER1,TUNER2,
     +   TUNER3,TUNER4,TUNER5,X0PRT,XCTOL,XFTOL
C

C     EXTERNAL FUNCTIONS
      real*8
     +   RMDCON
      INTEGER
     +   IMDCON
      EXTERNAL RMDCON,IMDCON
C
C     INTRINSIC FUNCTIONS
      INTRINSIC MAX
C
C  ***  SUPPLY NL2SOL (VERSION 2.2) DEFAULT VALUES TO IV AND V  ***
C
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER AFCTOL, COSMIN, COVPRT, COVREQ, DECFAC, DELTA0, DFAC,
C    1        DINIT,  DLTFDC, DLTFDJ, DTYPE,  D0INIT, EPSLON, FUZZ,
C    2        INCFAC, INITS,  JTINIT, LMAX0,  MXFCAL, MXITER, OUTLEV,
C    3        PARPRT, PHMNFC, PHMXFC, PRUNIT, RDFCMN, RDFCMX,
C    4        RFCTOL, RLIMIT, SOLPRT, STATPR, TUNER1, TUNER2, TUNER3,
C    5        TUNER4, TUNER5, XCTOL,  XFTOL,  X0PRT
C
      common / STEPJ  / stepJac, stepBound
	 DOUBLE PRECISION PARFCTOL, PARXCTOL
	 COMMON /CONVERG/ PARFCTOL, PARXCTOL !ax!
       integer keyCov                      !ax!
       common/KeyCov/keyCov                !ax!		          
c
      SAVE
c
      DATA ONE/1.0D0/, THREE/3.0D0/
C
C  ***  IV SUBSCRIPT VALUES  ***
C
      DATA COVPRT/14/, COVREQ/15/, DTYPE/16/, INITS/25/,
     +     MXFCAL/17/, MXITER/18/, OUTLEV/19/,
     +     PARPRT/20/, PRUNIT/21/, SOLPRT/22/,
     +     STATPR/23/, X0PRT/24/
C
C  ***  V SUBSCRIPT VALUES  ***
C
      DATA AFCTOL/31/, COSMIN/43/, DECFAC/22/,
     +     DELTA0/44/, DFAC/41/,   DINIT/38/,  DLTFDC/40/,
     +     DLTFDJ/36/, D0INIT/37/, EPSLON/19/, FUZZ/45/,
     +     INCFAC/23/, JTINIT/39/, LMAX0/35/,  PHMNFC/20/,
     +     PHMXFC/21/, RDFCMN/24/, RDFCMX/25/,
     +     RFCTOL/32/, RLIMIT/42/, TUNER1/26/,
     +     TUNER2/27/, TUNER3/28/, TUNER4/29/,
     +     TUNER5/30/, XCTOL/33/,  XFTOL/34/
C
C-----------------------------------------------------------------------
C
c       StepJAC=1.0D-8
cAXV        STEPBOUND=1.0
c       PARFCTOL=1.0D-6
c       PARXCTOL=1.0D-6
c       keyCov=2
cAXV      write(*,*) "Entering module DFAULT for NL2SOL options"
      IV(1)         =  12
      IV(COVPRT)    =  0                  ! Supress covar. matr. output.
cAXV:
       if	(keyCov.gt.3) then
	    IV(COVREQ)    =  2
       else
		  IV(COVREQ)    =  keyCov
	 endif
      IV(DTYPE)     =  1
      IV(INITS)     =  0
      IV( MXFCAL )  =  1000
      IV( MXITER )  =  MAX( 1000, 1 )
      IV(OUTLEV)    =  -1                 ! Iteration summary info.
      IV(PARPRT)    =  1
      IV(PRUNIT)    =  IMDCON(1)
!!!      IV(PRUNIT)    =  0                         ! Supress ANY output.
      IV(SOLPRT)    =  0                  ! Supress solution printout.
      IV(STATPR)    =  0                  ! Supress statistics printout.
      IV(X0PRT)     =  0                  ! Supress init. set. printout.

      MACHEP        = RMDCON(3)
      V(AFCTOL)     = MACHEP
      IF (MACHEP .GT. 1.0D-10) V(AFCTOL) = MACHEP**2
      V(COSMIN)     = MAX(1.0D-7, 1.0D2 * MACHEP)
      V(DECFAC)     = 0.5D0
      SQTEPS        = RMDCON(4)
      V(DELTA0)     = SQTEPS
      V(DFAC)       = 0.6D0
      V(DINIT)      = 0.0D0
      MEPCRT        = MACHEP ** (ONE/THREE)
      V(DLTFDC)     = MEPCRT
!!!      V(DLTFDJ)     = SQTEPS
      V(DLTFDJ)     = stepJac          ! Finite-difference REL/ABS step for
                                       ! the Jacobian estimation.
      V(D0INIT)     = 1.0D0
      V(EPSLON)     = 0.1D0
      V(FUZZ)       = 1.5D0
      V(INCFAC)     = 2.0D0
      V(JTINIT)     = 1.0D-6

      V(LMAX0)      = MIN( 100.0D0, stepBound )
!!!      V(LMAX0)      = 0.1D+000      ! Upper bound on the length of the
                                       ! initial search step.
      V(PHMNFC)     = -0.1D0
      V(PHMXFC)     = 0.1D0
      V(RDFCMN)     = 0.1D0
      V(RDFCMX)     = 4.0D0
cAXV      V(RFCTOL)     = MAX(1.0D-10, MEPCRT**2)
cAXV	      V(RFCTOL)     = 1.0D-6 !My probe
		V(RFCTOL)     =	PARFCTOL
      V(RLIMIT)     = RMDCON(5)
      V(TUNER1)     = 0.1D0
      V(TUNER2)     = 1.0D-4
      V(TUNER3)     = 0.75D0
      V(TUNER4)     = 0.5D0
      V(TUNER5)     = 0.75D0
cAXV      V(XCTOL)      = SQTEPS
cAXV           V(XCTOL)      = 1.0D-6 
             V(XCTOL)      = PARXCTOL
      V(XFTOL)      = 10.0D0 * MACHEP
C      V(XFTOL)      = 1.1*PARXCTOL

      RETURN
C  ***  LAST CARD OF DFAULT FOLLOWS  ***
      END
c-end-text nl2defu
