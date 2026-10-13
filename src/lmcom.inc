c Version 1.4  12/29/31
c----------------------------------------------------------------------
c                       =========================
c                          common block LMCOM
c                       =========================
c     Work arrays for Levenberg-Marquardt least squares fitting
c     procedure and parameters controlling termination criteria.
c
c----------------------------------------------------------------------
      double precision corr,diag,fjac,fvec,fvecPure,gnvec,gradf,qtf,
     #                work1,work2,work3,work4,xfit,ftol,xtol,gtol,factor
cAXV    tcov,
      integer ipvt,nprint
c
      common /lmcom/ fjac(MXPT,MXJCOL),fvec(MXPT),fvecPure(MXPT),
     #               xfit(MXJCOL), diag(MXJCOL),qtf(MXJCOL),
     #               corr(MXJCOL,MXJCOL), work1(MXJCOL),work2(MXJCOL),
     #               work3(MXJCOL),work4(MXPT),gnvec(MXJCOL),
     #               gradf(MXJCOL),ftol,gtol,xtol,factor,
     #               ipvt(MXJCOL),nprint
cAXV     #               tcov(MXTV,MXTV)
c
