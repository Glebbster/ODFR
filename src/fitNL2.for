cAXV  18/01/16
c----------------------------------------------------------------------
c                    =========================
c                         program fitNL2
c                    =========================
c
c Carry out nonlinear least-squares fitting using NL2sol
c
      subroutine fitNL2
      use optimization_steps, only: steps_begin,steps_end
c      implicit none
      
      include 'dimnls.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'lmcom.inc'      
c      
      include 'constNL2SOL.inc'
      include 'mpi_state.inc'
      integer itersDone      
      common / NL2CON /itersDone
      DOUBLE PRECISION PARFCTOL, PARXCTOL
      COMMON /CONVERG/ PARFCTOL, PARXCTOL 
      DOUBLE PRECISION stepJac, stepBound      
      common / STEPJ  / stepJac, stepBound
      double precision xbest(MXVAR), bestsum
      common /BESTFIT/ xbest, bestsum
c      
      character*14 AA0
      character*14 CorTmp(MXVAR)            
c            
      external      lfun
      external      MADR1
      
c
            stepJac=gtol
            stepBound      =  1.0d0    
c
            itersDone  =  0            ! NL2SOL iteration counter.
            PARFCTOL = ftol
            PARXCTOL = xtol
            
        do i = 1, MIVD
         iv(i)  =  0 
        end do
        do i = 1, MVD
         v(i)  =  0
        end do
        call DFAULT(IV,V)
            iv(1)=12
            iv(21)=2        
      if (mpi_rank.eq.0) OPEN ( UNIT=2,FILE='history',
     #                          POSITION='APPEND')

      bestsum = 1.0d300
      do i=1,nprm
         xbest(i)=xfit(i)
      enddo
      
      call steps_begin(ndatot,nprm)
      call NL2SNO(ndatot,nprm,xfit,MADR1,iv,v,uIparm,uRparm,lfun)
      call steps_end()
      do i=1,nprm
         xfit(i)=xbest(i)
      enddo

c           
              ncovfirst=iv(26)
          if ((ncovfirst.lt.0).and.(mpi_rank.eq.0))
     #       write(2,*)'  !!! Hessian is indefined  !!!'
          if ((ncovfirst.eq.0).and.(mpi_rank.eq.0))
     #       write(2,*)'  No attempt to compute covariance was made' 
      if (ncovfirst.gt.0)then
          do i=1,nprm
              do j=1,i
              corr(i,j)=v(ncovfirst-1+j)
              enddo
          ncovfirst=ncovfirst+i
          enddo
      endif 

      if ((iv(1).gt.0).and.(ncovfirst.gt.0)) then
       do i=1,nprm
          if (mpi_rank.eq.0) write(2,111)xfit(i),sqrt(corr(i,i))
  111       FORMAT(1x,'var=',G11.4,' Sigma=',G11.4)
          enddo
       AA0=' undefined   ,' 
       if ((iv(1).gt.2).and.(iv(1).lt.7)) then
         do i=1,nprm
           do j=1,i
              if ((corr(i,i)*corr(j,j)).eq. 0.0D0) then
               CorTmp(j)=AA0  
            else 
               write(CorTmp(j),505) corr(i,j)/sqrt(corr(i,i)*corr(j,j))
            endif
           enddo
               if (mpi_rank.eq.0) write(2,506)(CorTmp(jj),jj=1,i)
         enddo
       endif
      endif
  505        FORMAT (G12.4,';')
  506        FORMAT (1X,'C',1X,40(A14))
       
  507    return
      end       

