c based on iterat.inc from NLSL Version 1.4 1/2/95
c----------------------------------------------------------------------
c                    =========================
c                      common block ITERAT
c                    =========================
c
c   Iteration counter, residual norm, and related statistical quantities
c   for nonlinear least-squares procedure
c----------------------------------------------------------------------
c
cAXV: rdsqCom  - sum of squire deviations
cAXV: rdsqSp  -  array, sum of squire deviations	for spectrum
cAXV: rdchsqW  - Residual weighted according to noise for each spectrum
cAXV: rdchsqUn - Unweighted residuals
cAXV: lossSseW/lossSpW   - reporting SSE from final optimizer residual
cAXV: lossSseUn/lossSpUn - reporting SSE from pure unweighted residual
  
      integer iter,iwflag
      double precision fnorm,chisqr,rdchsqW,rdchsqUn,ch2bnd,chnbnd,
     #       confid,delchi1,qfit,f2bnd,fnbnd,tbound,
     #       rdsqCom,rdsqSp(MXSPC),lossSseW,lossSseUn,
     #       lossRedW,lossRedUn,lossSpW(MXSPC),lossSpUn(MXSPC)
      common /iterat/ fnorm,chisqr,rdchsqW,rdchsqUn,ch2bnd,chnbnd,
     #      confid,delchi1,qfit,f2bnd,fnbnd,tbound,
     #      rdsqCom,rdsqSp,lossSseW,lossSseUn,lossRedW,lossRedUn,
     #      lossSpW,lossSpUn,iter,iwflag







