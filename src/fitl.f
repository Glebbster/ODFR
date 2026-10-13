cAXV based on FITL.f from NLSL Version 1.5 beta 11/23/95
c----------------------------------------------------------------------
c                    =========================
c                         program FITL
c                    =========================
c
c Carry out nonlinear least-squares analysis of slow-motional spectra
c the EPRLL calculation. Uses a modification of the Levenberg-Marquardt
c trust-region algorithm available in MINPACK. 
c
c----------------------------------------------------------------------
c
      subroutine fitl(info)
      implicit none
c
      include 'dimnls.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'lmcom.inc'
      include 'errmsg.inc'
      include 'iterat.inc'
c
      integer nfev,njev,info
      double precision enorm,gammq
      double precision chisqrPure,rdchsqPure     
      common /rdchsqPure/rdchsqPure  
      double precision   DndatNprm     
      common /DndatNprm/ DndatNprm         

      external lfun, xpack, enorm
	logical convrg
c #####################################################################
c    
        call lmnls( lfun,ndatot,nprm,xfit,fvec,fjac,MXPT,ftol,xtol,gtol,
     1          MAXEV,MAXIT,diag,prscl,factor,1,info,nfev,njev,
     2          ipvt,qtf,gnvec,gradf,work1,work2,work3,work4,fvecPure )

c                                         *** MINPACK error return
           convrg = info.ne.0 .and. info .le.3
           if (.not. convrg) then
              write (2,2034) minerr(info)
c              write (luout,2037) nfev,cpu 
              write (*,2034) minerr(info)
 2034 format(/2x,'MINPACK could not find a solution: ', a32 )
c                                         *** Normal return
           else
              write(2,2035) minerr(info)
              write(*,2035) minerr(info)
 2035 format(/1x,'MINPACK completed: ', a32) 
              fnorm=enorm(ndatot,fvec)
              chisqr=fnorm*fnorm
              if (iwflag.ne.0) then
c                                         --- weighted residual fit ---
                 rdchsqW=chisqr/DndatNprm 
            qfit=gammq( 0.5d0*DndatNprm, 0.5*chisqr )
            write (2,2038) qfit
 2038 format(12x,'Goodness of fit from chi-squared (Q)=',g13.6)
              else
c                                         --- unweighted residual fit ---
                 rdchsqUn=chisqr/DndatNprm
              end if
              fnorm=enorm(ndatot,fvecPure)
              chisqrPure=fnorm*fnorm         
              rdchsqPure=chisqrPure/DndatNprm              

               call covrpt( 2 )
           endif
c
       return
       end
c

