c NLSL Version 1.5 beta 11/23/95
c----------------------------------------------------------------------
c                    =========================
c                       subroutine COVRPT  
c                    =========================
c
c    Produces a report of the covariance matrix and parameter 
c    uncertainties on the specified logical unit.
c    Some reducted by AXV
c----------------------------------------------------------------------
c
      subroutine covrpt( lu )
      implicit none
      integer lu
c      
      include 'dimnls.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'lmcom.inc'
      include 'iterat.inc'
      include 'rndoff.inc'
cAXV      include 'stdio.inc'
c
      integer i,j,jj	!,k,isi,isp,mrgn

      double precision denom   !
cAXV      character*132 line
c
      integer itrim
      double precision chi,ts,fs
      external itrim,chi,ts,fs
c
c-----------------------------------------------------------------------
c  Calculate residuals and covariance matrix
c  COVAR calculates covariance matrix from the R matrix of QR 
c  factorization, which lmder stored in upper triangular form in fjac 
c  ipvt is permutation array returned by lmder
c------------------------------------------------------------------------
c
         call covar( nprm,fjac,MXPT,ipvt,xtol*xtol,work1 )
c                                    --- weighted residual fit
         if (iwflag.ne.0) then
            delchi1=chi(confid,1)
            ch2bnd=chi(confid,2)
            chnbnd=chi(confid,nprm)
c                                    --- unweighted residual fit
         else
            tbound=ts( (1.0D0-confid)/2.0D0, ndatot-nprm )
            f2bnd=2.0D0*fs( 1.0D0-confid, 2, ndatot-2 )*fnorm*fnorm/
     #            float(ndatot-2)
            fnbnd=dfloat(nprm)*fs( 1.0D0-confid, nprm, ndatot-nprm )*
     #            fnorm*fnorm/dfloat(ndatot-nprm)
         end if
c
      do i=1,nprm
CAXV         line=' '
         do j=i,nprm
c
c----------------------------------------------------------------------
c Estimate of error bounds:
c     For weighted residuals, use chi-squared statistics
c     For unweighted residuals use t-statistics
c----------------------------------------------------------------------
            if (i.eq.j) then
               if (iwflag.ne.0) then
                  xerr(i)=sqrt(delchi1*fjac(i,i))
               else
                  xerr(i)=tbound*fnorm*sqrt(fjac(i,i))
     #                    /sqrt(dfloat(ndatot-nprm))
               end if
            end if

            denom = dsqrt( fjac(i,i) * fjac(j,j) )
            if (denom .ne. 0.0D0) corr(i,j) = fjac(i,j)/denom
cAXV            k=10*(j-1)+mrgn
         end do
         write(lu,2003) (corr(i,jj),jj=1,nprm)
cAXV         write (lu,2001) line(:10*(nprm+1)+mrgn)
      end do
c     
c    --------------------------------
c     output final fit of parameters 
c    --------------------------------
cAXV      write (lu,1000)
      if (iwflag.ne.0) then
         write (lu,2004) 'Chi'
      else
         write(lu,2004) 'T'
      end if
c
cAXV
      write (lu,2007)(i,xfit(i),xerr(i),i=1,nprm)
c
c                              --- weighted residual fit ---
      if (iwflag.ne.0) then
         write (lu,2008) confid,delchi1,ch2bnd,nprm,chnbnd
c
c                              --- weighted residual fit ---
      else
         write (lu,2009) confid,tbound,f2bnd,nprm,fnbnd
      end if
c

      return
c
c##### Formats #########################################################
c
 1000 format(/,2x,70('-'),/)
 1010 format('SITE',i1)
 1011 format('SPCTRM',i1)
 2000 format(/24x,'*** Correlation Matrix ***'/)
 2001 format(a)
 2002 format(12(2x,a8))
 2003 format(20(1x,f8.4))
cAXV 2004 format(/9x,'*** Final Parameters ***',/,
cAXV     1       7x,' Parameter      Value',13x,'Uncertainty (',a,
cAXV     2       ' statistics)'/7x,46('-')/)

 2004 format(7x,' Parameter      Value',13x,'Uncertainty (',a,
     2       ' statistics)')
cAXV:
 2007 format(2x,i3,5x,g14.7,' +/- ',g14.7)
 2008 format(/7x,'Confidence= ',f6.3,3x,'Del-Chi-Sqr= ',g12.5/
     #        7x,'Del-Chi-Sqr for confidence regions:'/
     #        7x,'2 parameters: ',g12.5,4x,i2,' parameters: ',g12.5)
 2009 format(/7x,'Confidence = ',f6.3,3x,'T-bound=',g12.5/
     #        7x,'F-bound for confidence region:'/
     #        7x,'2 parameters: ',g12.5,4x,i2,' parameters: ',g12.5)
 1060 format(/10x,'Order parameters for site',i2,':')
c
      end


