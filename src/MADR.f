cAXV, for ODFR
      subroutine MADR1( n, p, x, nf, r, uiParm, urParm, lfun, iv )
!c
!c     Intermediate program calling lfun.
!c
      include 'dimnls.inc'
      include 'constNL2SOL.inc'
      include 'expdat.inc'
      include 'iterat.inc'
      include 'mpi_state.inc'
      external      lfun      
      integer       n,nf,p
      integer       iflag,ldfjac,itersDone,i,iter_now,isp
!c
      DOUBLE PRECISION  x(p), r(MXPT), fvecPure(MXPT)
      double precision xbest(MXVAR), bestsum, sse, loss_unweighted
      common /BESTFIT/ xbest, bestsum
      common / NL2CON /itersDone
!c
      ldfjac = n    !it does not matter

      IF (ITERSDONE.LT.IV(31))THEN
          iter_now = iv(31)

          iflag=0
          call lfun(n,p,x,r,fvecPure,fjac,ldfjac,iflag)

          iflag=1
          call lfun(n,p,x,r,fvecPure,fjac,ldfjac,iflag)

          sse = 0.0d0
          do i=1,n
             sse = sse + r(i)*r(i)
          enddo

          loss_unweighted = 0.0d0
          do i=1,n
             loss_unweighted = loss_unweighted + fvecPure(i)*fvecPure(i)
          enddo

          if ((iter_now.eq.1).and.(mpi_rank.eq.0)) then
             write(*,'(A)') 
     #      '  it        LOSS_W        LOSS_UNW   (spectrum contrib below)'
          endif
          if (mpi_rank.eq.0) write(*,'(I4,2X,1PE12.5,2X,1PE12.5)')
     #      iter_now, sse, loss_unweighted
          if (mpi_rank.eq.0) then
             write(*,'(A)') '  Spectrum SSE contributions:'
             do isp=1,nspc
                write(*,'(4X,I4,8X,1PE14.6,8X,1PE14.6)')
     #             isp,lossSpW(isp),lossSpUn(isp)
             enddo
          endif

          itersDone = iter_now

          if (sse .lt. bestsum) then
             bestsum = sse
             do i=1,p
                xbest(i) = x(i)
             enddo
          endif

          return
      ENDIF

      iflag=1
      call lfun(n,p,x,r,fvecPure,fjac,ldfjac,iflag)

      sse = 0.0d0
      do i=1,n
         sse = sse + r(i)*r(i)
      enddo

      if (sse .lt. bestsum) then
         bestsum = sse
         do i=1,p
            xbest(i) = x(i)
         enddo
      endif
      
!c 
      RETURN
      end
