c Version 1.5.1 beta 03/07/96
c----------------------------------------------------------------------
c                    =========================
c                         common ERRMSG
c                    =========================
c
c Contains text for error messages and definitions of symbols used
c to define error message numbers. Errors after ZEROB0 are defined
c as fatal, and will cause the calculation to halt.
c
c Error message texts are set in a block data section in nlstxt.f
c Definitions of the symbols are as follows:
c
c ZEROIN2     Zero nuclear spin assumed
c POSB0       Positive B0 assumed
c BADHIFLD    High-field approx. may not apply
c BADJMP      Nonzero potential with jump/free: ipdf set to 0
c BADAV       Zero potential with aniso. viscos.: ipdf set to 0
c BADDJ       Discrete jumps with aniso. viscos.: ipdf set 0
c LEMXHI      lemx must be 48 or less with potenti
c BSSADJ      Basis set adjusted to (iii,iii,iii,iii,iii,iii,ii)
c NSTEPHI     NSTEP too large: adjusted to iiiii
c CGTOLHI     CGTOL too small: adjusted to xxxxxxxxx
c SHIFTHI     SHIFTR too small: adjusted to xxxxxxxxx
c NOCONVRG    CG did not converge after iiiii steps
c NOALPHAM    Angle alm not needed with axial g-matrix
c NOALPHAD    Angle ald not needed with axial diffusion tensor
c ZEROB0      Zero B0
c ZEROG       Zero g-tensor element
c TDGERR      Error in tridiagonal matrix storage
c MTXBIG      Too many matrix elements
c DIMBIG      Matrix dimension too large
c MTXHLT      User halt during matrix calculation
c CGHLT       User halt during CG tridiagonalization
c TDGBIG      Too many tridiagonal matrix elements
c NOMIN       Bracketing routine (search command) did not find a minimum
c ZERONFLD    NFIELD is zero (no points in spectrum)
c ZRORANGE    RANGE is zero 
c----------------------------------------------------------------------

      integer FATAL,NEPERR
      parameter (FATAL=16,NEPERR=26)
      character*50 eprerr(NEPERR)
      character*32 minerr(0:11)
c
       integer ZEROIN2,POSB0,BADHIFLD,BADJMP,BADAV,BADDJ,LEMXHI,BSSADJ,
     #         NSTEPHI,CGTOLHI,SHIFTHI,NOCONVRG,NOALPHAM,NOALPHAD,
     #         BADBESS,ZEROB0,ZEROG,TDGERR,MTXBIG,DIMBIG,MTXHLT,CGHLT,
     #         TDGBIG,NOMIN,ZERONFLD,ZRORANGE
c
       parameter (ZEROIN2=1,POSB0=2,BADHIFLD=3,BADJMP=4,BADAV=5,BADDJ=6,
     #            LEMXHI=7,BSSADJ=8,NSTEPHI=9,CGTOLHI=10,SHIFTHI=11,
     #            NOCONVRG=12,NOALPHAM=13,NOALPHAD=14,BADBESS=15,
     #            ZEROB0=16,ZEROG=17,TDGERR=18,MTXBIG=19,
     #            DIMBIG=20,MTXHLT=21,CGHLT=22,TDGBIG=23,NOMIN=24,
     #            ZERONFLD=25,ZRORANGE=26)

       data minerr / 'Does not complete',
     1               'Chi-squared convergence',
     2               'X vector convergence',
     3               'Chi-squared/X vector convergence',
     4               'Gradient convergence',
     5               'Max function evaluations reached',
     5               'Max iterations reached',
     6               'FTOL too small',
     7               'XTOL too small',
     8               'Zero Jacobian or GTOL too small',
     9               'Terminated internally',
     A               'Single calculation completed'/

      data eprerr /
     1     'Zero nuclear spin assumed',
     2     'Positive B0 assumed',
     3     'High-field approx. may not apply',
     4     'Nonzero potential with jump/free: ipdf set to 0',
     5     'Zero potential with aniso. viscos.: ipdf set to 0',
     6     'Discrete jumps with aniso. viscos.: ipdf set 0',
     7     'lemx must be 48 or less with potential',
     8     'Basis set adjusted to (iii,iii,iii,iii,iii,iii,ii)',
     9     'NSTEP too large: adjusted to iiiii',
     A     'CGTOL too small: adjusted to xxxxxxxxx',
     1     'SHIFTR too small: adjusted to xxxxxxxxx',
     2     'CG did not converge after iiiii steps',
     3     'Angle alpham not needed w/ axial g-matrix',
     4     'Angle alphad not needed w/ axial diffusion tensor',
     5     'Taylor series for Bessel function did not converge',
     6     'Zero B0',
     7     'Zero g-tensor element',
     8     'Error in tridiagonal matrix storage',
     9     'Too many matrix elements',
     A     'Matrix dimension too large',
     1     'User halt during matrix calculation',
     2     'User halt during CG tridiagonalization',
     3     'Too many tridiagonal matrix elements',
     4     'Bracketing did not find minimum within step bound',
     5     'Number of points in spectrum is zero',
     6     'Field range of spectrum is zero' /
     
      
c      common /errmsg/ eprerr,minerr



