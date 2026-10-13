c based on NLS Version 1.5.1 beta 2/3/96
c*********************************************************************
c                    =========================
c                      include file DIMNLS
c                    =========================
c
c  Defines array dimensioning parameters for NLS programs
c
c MXSTEP  Maximum number of CG/Lanczos steps
c MXDIM   Maximum dimension of matrix
c MXEL    Maximum number of off-diagonal elements in matrix
c MXINP   Maximum number of points in an input file
c MXSPC   Maximum number of spectra that can be fit at once
c MXSITE  Maximum number of sites that can be fit for a given spectrum
c MXSPT   Maximum number of points in an individual spectrum
c         NOTE: THIS SHOULD BE A POWER OF 2
c MXPT    Maximum total number of data points
c MXVAR   Maximum number of parameters that may be varied
c MXJCOL  Maximum number of columns in Jacobian matrix for L-M least
c         squares (should be MXVAR+MXSITE) 
c      integer NFPRM	 length fparm
c			   NIPRM	 lenght iparm
cAXV additionally:
c MAXIT   Maximum number of iteration in Levenberg-Marquardt fitting
c MAXEV	  Maximum number of function evaluation in course of 
c         Levenberg-Marquardt fitting
c MAXNORT Maximum number of orientations in MOMD/ODF calculation
c 
c*********************************************************************
c
      integer MXSTEP,MXDIM,MXDIM1,MXEL,MXSPC,MXSPT,MXPT,MXVAR,MXJCOL,
     #        MXSITE,MXINP,
     #        MAXIT,MAXEV,MaxAddNucl,
     #   	  MAXNORT,MAXStorMat
c
      parameter (MXDIM=50000,MXDIM1=MXDIM+1,
     #           MXSTEP=MXDIM, ! MXSTEP=2000,
     #           MXEL=15000000,
     #           MXSPC=20,
     #           MXSPT=4096,
     #           MXPT=MXSPT*MXSPC,
     #           MXVAR=120,
     #           MXINP=4096,
     #           MXSITE=3,
     #           MXJCOL=MXVAR+MXSITE,
     #           MAXIT=1000,
     #           MAXEV=1000000,
     #           MaxAddNucl=15,
     #           MAXNORT=40,
cdb     #           MAXNORT=40,
c     #           MAXStorMat=4*MAXNORT )     
     #           MAXStorMat=MXSPC*MAXNORT )
c
      integer NFPRM,NIPRM
      parameter( NFPRM=60+2*MaxAddNucl, NIPRM=38+MaxAddNucl)
cAXV  ! fparms: 46-s - AMPFAC, 47-s -SPSHIFT, 48-s - ANGFI,49-s - ANGKSI,
cAXV            53 - Anglib,60 - IAadd,61 - Aadd(MaxAddNucl)
cAXV  ! iparms:
       INTEGER, PARAMETER :: MaxCALC=10    !max number of consecutive simulations 
       INTEGER, PARAMETER :: MaxParODF=150 !max number of parameters in ODF description (199)
       INTEGER, parameter :: Jmax=30       !max order of ODF Wigner decomposition (24)
       INTEGER, PARAMETER :: mxlval=190	   !LEMX<mxval, used in w3j.f
c
c*********************************************************************
c

