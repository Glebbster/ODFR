cAXV based on LFUN.f fromNLSL Version 1.5 beta 11/23/95
c----------------------------------------------------------------------
c                         =======================
c                            subroutine LFUN
c                         =======================
c
c Subroutine for interfacing EPRLL spectral calculations with the MINPACK 		 
c version of the Levenberg-Marquardt nonlinear least squares algorithm 
c for fitting experimental spectra. The function performed
c function determined by the IFLAG argument. 
c IFLAG=0 (Printing of iterates)
c IFLAG=1 (Evaluation of residuals function)
c IFLAG=2 (Jacobian evaluation)
c IFLAG=3 (Fit is terminating. )
c---------------------------------------------------------------------- 
c
      subroutine lfun(mdum,ndum,xfit,fvec,fvecPure,fjac,ldfjac,iflag)
      implicit none
      integer iflag,ldfjac,mdum,ndum	!,lthtd
c
      include 'dimnls.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'errmsg.inc'
      include 'pidef.inc'
      include 'prmeqv.inc'
      include 'iterat.inc'
      include 'rndoff.inc'
      include 'dimodfl.inc'	 ! Jmax is the maximum rank of decomposition
      include 'parODF.inc'
      include 'mpi_state.inc'
c      include 'temp.inc'      
cAXV:
      double precision xfit(MXJCOL),resSp(MXSPT)
      double precision fvec(MXPT),fjac(LDFJAC,MXJCOL)
      double precision fvecPure(MXPT),chisqrPure,rdchsqPure 
      common /rdchsqPure/rdchsqPure           
cdb            double precision fjac1(LDFJAC)

      integer icalc,ierr,isi,isp,ix,ixs,iort,Kcontrol,numODFs	!,ixt
      double precision xtemp
      double precision  shtraftmp,delshtraf
      double precision  enorm,fnormtmp
c
      integer FULL,CFONLY
      double precision ZERO
      parameter (FULL=111,CFONLY=0,ZERO=0.0D0)
c
c      double precision ??? spectmp(MXSPT,MXSITE)
     	integer jpt
      logical flfirst,DiscTrid     	
c   For preservation the values in course of
c   consecutive calling the subroutine lfunR.f:     	
      double precision xfitOld, shtraf
      common /templfun/ xfitOld(MXJCOL),shtraf(MXSPC),DiscTrid
     	double precision PromSpecDat(MXSPT)
c   Add weight flag and list of weight for weighted calc of residuals      
      double precision wgt(MXPT)
      double precision wtmp
      logical haveWgt
      common /wgtpt/ wgt, haveWgt

      integer iparms(NIPRM,MXSITE,MXSPC)
      DOUBLE PRECISION fparms(NFPRM,MXSITE,MXSPC)
        COMMON /COMPARMS/ iparms,fparms
      integer jsim,itemp,ii,jtmp,jnext,numstor(MXSITE,MXSPC)
      save numstor
        common/flfst/flfirst,jsim
      logical depOK,useMOMD
      integer modtd,iflvar,filenprm
        common/depOK/ depOK(MXVAR,MXSITE,MXSPC),modtd(MXSITE,MXSPC),
     #                useMOMD(MXSITE,MXSPC) 
      common /iflvar/iflvar(NFPRM,MXSITE,MXSPC),filenprm(MXSPC)
            
c
      integer jODF,ioerr
        common /A_jmn/A_jmn,jODF,numODFs
        common /control/Kcontrol
      double precision   DndatNprm     
      common /DndatNprm/ DndatNprm         
        
c
        external  enorm	
c
       character*33  TempName   
       character*15  tmpNams(MXSITE*MXSPC*MAXNORT)
       character*2  NamTmp,NamTmp1,NamTmp2
       character*4, PARAMETER :: temp='temp' 
       character*30, PARAMETER :: DummyName=' '        
c for debugging
c       double precision TempField(MXINP), TempSpec(MXINP)
c######################################################################
c
      ierr=0
c      
      OPEN ( UNIT=1,FILE='temp 1.par',status='old',IOSTAT=ioerr)
      if (ioerr.eq.0) then
     	       READ(1,*,IOSTAT=ioerr) Kcontrol
             if (ioerr.ne.0) Kcontrol=3
         CLOSE ( UNIT=1 )
      else
         Kcontrol=3
      endif
      if (Kcontrol.eq.1) iflag=3
c**********************************************************************
c     IFLAG=0: Output iteration information
c**********************************************************************
      if (iflag.eq.0) then
cAXV         write(*,*) 'N'
        if (mpi_rank.eq.0) then
         if (Kcontrol.eq.3)then
1021        format(1x,20(1x,G15.8))
         else
102         format(1x,i3,20(1x,G15.8))
         endif
         do isp=1,nspc
           write(NamTmp,101)isp
101        FORMAT (I2)
           CLOSE ( UNIT=2)
           OPEN ( UNIT=2,FILE='history', POSITION='APPEND')

          TempName=temp//NamTmp//'.dat'	 	  
          call SpecOUT(DummyName,TempName,isp,TempName,
     #                 DummyName,DummyName,0)
          TempName=temp//NamTmp//'.par'
          call PARWrite(TempName,isp)
         enddo ! isp
          TempName=temp//'ODF.par'	 	  			 	  
          call ODFwrite(TempName)
        endif
c**********************************************************************
c     IFLAG=1: Function evaluation
c**********************************************************************
      else if (iflag.eq.1) then
c              
ctemp        if (Kcontrol.ne.2) call xpack(2,itemp,xfit) 
        call xpack(2,itemp,xfit) 
c
        if (flfirst) then
c  using the disc has removed        
           disctrid=.false. 
c   Number of independent orientation functions:            
           numODFs=0


c  ----------------------------------------------------------
c   Optional point-weights (one coefficient per data point).
c   Enabled when iwflag == 4 (weights only) or iwflag == 5
c   (weights + /rmsn). If file is missing/short -> defaults to 1.
c  ----------------------------------------------------------
c --- default: no weights
        haveWgt = .false.
        do jpt = 1, ndatot
          wgt(jpt) = 1.d0
        enddo

        if (iwflag.eq.4 .or. iwflag.eq.5) then
          open(unit=97,file='weight',status='old',err=909)

c ---   read exactly ndatot values
915       do jpt = 1, ndatot
              read(97,*,end=911,err=911) wgt(jpt)
          enddo

c --- check for extra values
          read(97,*,end=913,err=912) wtmp

c --- too many values
912       continue
          close(97)
          haveWgt = .false.
          do jpt = 1, ndatot
              wgt(jpt) = 1.d0
          enddo
          if (mpi_rank.eq.0) write(*,*)
     #       'Weight file size mismatch: too many values, expected',
     #       ndatot
          goto 914

c --- exact match
913       continue
          close(97)
          haveWgt = .true.
          if (mpi_rank.eq.0) write(*,*) 'Weights OK, ndatot =',ndatot
          goto 914

c --- too few values
911       continue
          close(97)
          haveWgt = .false.
          do jpt = 1, ndatot
              wgt(jpt) = 1.d0
          enddo
          if (mpi_rank.eq.0) write(*,*)
     #       'Weight file size mismatch: too few values, expected',
     #       ndatot
          goto 914

c --- fallback file name
909       continue
          open(unit=97,file='weight.dat',status='old',err=910)
          goto 915

c --- file not found
910       continue
          haveWgt = .false.
          if (mpi_rank.eq.0) write(*,*)
     #       'Weight file not found, using wgt=1'

914       continue
        endif
        do isi=1,nnsite(1)
           if ((ODFtype(1).ge.2).and.(keyODF(isi).eq.0))
     #                 numODFs=numODFs+1  
        enddo                   
c allocate room for ort spectra storage (alpha,betta) 
          jnext=1
          do isp=1,nspc 
            do isi=1,nnsite(isp)
              if (useMOMD(isi,isp)) then 
c spectra/sites for calculation using iort spectra for first spectrum:              
                modtd(isi,isp)=0	
                numstor(isi,isp)=numstor(isi,1)
              else
                if (iparms(IItypeSite,isi,isp).gt.0) then
                    modtd(isi,isp)=0
                    numstor(isi,isp)=
     #                    numstor(iparms(IItypeSite,isi,isp),isp)
                else      
                    modtd(isi,isp)=1 ! mark spectra for first full calculation
                    numstor(isi,isp)=jnext
                    jnext=jnext+iparms(INORT,isi,isp)
                endif    
              endif !useMOMD(isi,isp)                
            enddo
          enddo
c
          if ((jnext-1).gt.MAXStorMat) then
            if (mpi_rank.eq.0) write(*,*)'Insufficient room
     #            for tridiagonal matrices storage'
            if (mpi_rank.eq.0) write(2,*)'Insufficient room
     #            for tridiagonal matrices storage'
             stop            
c
          endif
c
          fldShftSite(1)=0.d0
          flfirst=.false.
        else  !not first
          do ix=1, nprm+nprm3
            if (ixx(ix).gt.0) then
              SELECT CASE (ixpr(ix))
              CASE (ISPSHFT,IAMPFAC,IGIB0,IGIB2,ILBL,ILBG)
                continue	 
              CASE (IAadd:NFPRM)
                continue	 
              CASE DEFAULT    
                if (xfit(ix).ne.xfitOld(ix))	then
                  modtd(ixst2(ix),ixst1(ix))=1
                  if (ixst1(ix).eq.1) then 
                    do ii=2,nspc
                      if (.not.useMOMD(ixst2(ix),ii))                    
     #            	        modtd(ixst2(ix),ii)=1
                    enddo
                  endif                
                endif  !xfit(ix).ne.xfitOld(ix)
              END SELECT
            endif !(ixx>0)
            if (ixx(ix).lt.0) then
              if (xfit(ix).ne.xfitOld(ix)) then
                do isp=1,nspc
                  do isi=1,nnsite(isp)
                    if ((ODFtype(isp).ge.2).and.
     #                  (keyODF(isi).ge.0)) then
                      modtd(isi,isp)=1
                    endif
                  enddo
                enddo
              endif
            endif !(ixx<0)
          enddo
        endif  !first
c       -----------------------------------
c        Loop over spectra 
c       -----------------------------------
        rdsqCom=0.0d0
        do isp=1,nspc        
c
          ixs=ixsp(isp)
          do jpt=1,npts(isp)
               ClSpecDat(jpt,isp,1)=0.0d0
          enddo
          shtraf(isp)=0.d0
c          --------------------------------------
c           Loop over all sites for the spectrum
c          --------------------------------------
          do isi=1,nnsite(isp)
c            fill ODF array:
             if (ODFtype(isp).ge.2) then
                call ODFLInit(isi,jODF,A_jmn)
                 if (keyODF(isi).le.1) call ODFint(knet,anegODF(isi))
                 if (PenFLG.GT.0) shtraf(isp)=
     #                      shtraf(isp)+anegODF(isi)*PenLev 
             endif  
c                        
             icalc=FULL
c Positive line width:
     	       fparms(IGIB0,isi,isp)=dabs(fparms(IGIB0,isi,isp)) 
     	       fparms(IGIBX,isi,isp)=dabs(fparms(IGIBX,isi,isp)) 
     	       fparms(IWXX,isi,isp)=dabs(fparms(IWXX,isi,isp)) 
     	       fparms(IWYY,isi,isp)=dabs(fparms(IWYY,isi,isp)) 
     	       fparms(IWZZ,isi,isp)=dabs(fparms(IWZZ,isi,isp)) 
c              ---------------------------------
c              Calculate an individual spectrum
c              ---------------------------------
cdb For debug:
cdb             icalc=FULL 
c             
             call momdls( fparms(1,isi,isp),iparms(1,isi,isp),icalc,
     #              resSp,nft(isp),isi,isp,numstor,tmpNams,
     #              DiscTrid,ierr )
c
             do jpt=1,npts(isp)
                  ClSpecDat(jpt,isp,isi+1)=resSp(jpt)
     &                        *fparms(IAMPFAC,isi,isp)*keySp(isp,2)
                  ClSpecDat(jpt,isp,1)=ClSpecDat(jpt,isp,1)+
     &	                             ClSpecDat(jpt,isp,isi+1)
             enddo
             modtd(isi,isp)=0
          enddo !isi
c                ----------------------------------------------------------
c                 For the current experimental spectrum calculate residuals
c                ----------------------------------------------------------
          rdsqSp(isp)=0.d0
          if (numODFs.gt.0) shtraf(isp)=shtraf(isp)/numODFs
          do jpt=1,npts(isp)
               jtmp=ixs+jpt-1
c               fvec(jtmp)=fitdata(jtmp)
               fvec(jtmp)=fitdata(jtmp)-ClSpecDat(jpt,isp,1)
               fvecPure(jtmp)=fvec(jtmp)               
c      Penalty for negative part of ODF:              
               fvec(jtmp)=fvec(jtmp)*(1.d0+shtraf(isp))
               if (haveWgt) fvec(jtmp)=fvec(jtmp)*wgt(jtmp)
               rdsqSp(isp)=rdsqSp(isp)+fvec(jtmp)**2
          end do

          rdsqCom=rdsqCom+rdsqSp(isp)
c     
c for debugging:
c          do jpt=1,npts(isp)
c             jtmp=ixs+jpt-1
c	       TempSpec(jpt)=fitdata(jtmp)
cc            TempSpec(jpt)=ClSpecDat(jpt,isp,1)
c	       TempField(jpt)=jpt*1.d0
c          enddo
c         call promout(npts(isp),TempField,TempSpec)
c	write(*,*)isp
c          do jpt=1,npts(isp)
c             jtmp=ixs+jpt-1
cc	       TempSpec(jpt)=fitdata(jtmp)
c            TempSpec(jpt)=ClSpecDat(jpt,isp,1)
c	       TempField(jpt)=jpt*1.d0
c          enddo
c         call promout(npts(isp),TempField,TempSpec)
c	write(*,*)isp
        end do  !*** End of loop over spectra 
        do ii=1, nprm+nprm3
	    xfitOld(ii)=xfit(ii) 
        enddo
c     ----------------------------------------------------------
c     For weighted residual fit, normalize differences in fvec to the spectral noise 
c     (then |fvec|^2 will approximate  the true chi^2)
c     ----------------------------------------------------------
         if ((iwflag.ge.1 .and. iwflag.le.3) .or. (iwflag.eq.5)) then
            do isp=1,nspc
               ixs=ixsp(isp)           
               do jpt=1,npts(isp)
                  jtmp=ixs+jpt-1
                  if (dabs(rmsn(isp)).gt.RNDOFF) then
                     fvec(jtmp)=fvec(jtmp)/rmsn(isp)
                  endif
               end do
            end do
         end if
c
         lossSseW=0.d0
         lossSseUn=0.d0
         rdsqCom=0.d0
         do isp=1,nspc
            lossSpW(isp)=0.d0
            lossSpUn(isp)=0.d0
            ixs=ixsp(isp)
            do jpt=1,npts(isp)
               jtmp=ixs+jpt-1
               lossSpW(isp)=lossSpW(isp)+fvec(jtmp)**2
               lossSpUn(isp)=lossSpUn(isp)+fvecPure(jtmp)**2
            enddo
            rdsqSp(isp)=lossSpW(isp)
            lossSseW=lossSseW+lossSpW(isp)
            lossSseUn=lossSseUn+lossSpUn(isp)
            rdsqCom=rdsqCom+rdsqSp(isp)
         enddo

         chisqr=lossSseW
         chisqrPure=lossSseUn
         if (dabs(DndatNprm).gt.RNDOFF) then
            lossRedW=lossSseW/DndatNprm
            lossRedUn=lossSseUn/DndatNprm
         else
            lossRedW=lossSseW
            lossRedUn=lossSseUn
         endif
         rdchsqW=lossRedW
         rdchsqUn=lossRedUn
         rdchsqPure=lossRedUn

c	write(*,*)rdchsqUn
      continue
c**********************************************************************
c IFLAG = 2: Jacobian evaluation by forward-difference approximation:
c             fjac = (f(x0+dx)-f(x0))/dx
c Note that f(x) = fitdata-calc(x), so that f(x0+dx)-f(x0) = calc(x0)-calc(x0+dx)
c**********************************************************************
      else if (iflag.eq.2) then
c        -------------------------------------------------------------
c         Loop over all parameters, introducing the forward-difference
c         step into each parameter
c        -------------------------------------------------------------
        do ix=1,nprm
          xtemp=xfit(ix)
          xfit(ix)=xfit(ix)+xfdstp(ix)
          call xpack(2,itemp,xfit) 
c             ---------------------
c              Initialize Jacobian
c             ---------------------
          do jpt=1,ndatot
               fjac(jpt,ix)=ZERO
cdb               fjac1(jpt)=ZERO
          end do
cAXV: Select the cases when the full recalculation of sprectrum
c    is not required.
          if (ixx(ix).gt.0) then
            SELECT CASE (ixpr(ix))
            CASE (ISPSHFT,IAMPFAC,IGIB0,IGIB2,ILBL,ILBG)
              continue
            CASE (IAadd:NFPRM)
              continue	 				 
            CASE DEFAULT    
              modtd(ixst2(ix),ixst1(ix))=1
              if (ixst1(ix).eq.1) then 
                do isi=1,nnsite(1)
                   if (iflvar(ixpr(ix),isi,1).gt.90) 
     #                       modtd(isi,1)=modtd(1,1)
                enddo     
                do ii=2,nspc
                  do isi=1,nnsite(1)                
                   if (.not.useMOMD(isi,ii))                    
     #                       modtd(isi,ii)=modtd(1,1)
                  enddo     
                enddo
              endif                
            END SELECT
		endif !(ixx>0)
          if (ixx(ix).lt.0) then
            do isp=1,nspc
              do isi=1,nnsite(isp)
                if ((ODFtype(isp).ge.2).and.
     #              (keyODF(isi).ge.0)) then
                  modtd(isi,isp)=1
                endif
              enddo
            enddo
          endif !(ixx<0)
c         ------------------------------------
c          Loop over all spectra 
c         ------------------------------------
          do isp=1,nspc
            ixs=ixsp(isp)
            shtraftmp=0.d0
            do jpt=1,npts(isp)
              PromSpecDat(jpt)=0.d0 !is necessary for straf calculation!
            enddo                        
c             ----------------------------------------
c              Loop over all sites in the spectrum
c             ----------------------------------------
            do isi=1,nnsite(isp)
c                fill ODF array:          
                 if (ODFtype(isp).ge.2) then
                    call ODFLInit(isi,jODF,A_jmn)
                    if (keyODF(isi).le.1) call ODFint(knet,anegODF(isi))
                    if (PenFLG.GT.0) shtraftmp= 
     #                         shtraftmp + anegODF(isi)*PenLev 
                 endif  
c  calculating the spectra/sites marked as dependent on the parameter
              if (depOK(ix,isi,isp).or.
     #            ((ixx(ix).lt.0).and.(ODFtype(isp).ge.2).and.
     #             (keyODF(isi).ge.0))) then
                 icalc=CFONLY
                 if (modtd(isi,isp).ne.0)  icalc=FULL
c

                 call momdls( fparms(1,isi,isp),iparms(1,isi,isp),icalc,
     #                    resSp,nft(isp),isi,isp,numstor,tmpNams,
     #                    DiscTrid,ierr )
                  modtd(isi,isp)=0
                  
                  do jpt=1,npts(isp)
                    PromSpecDat(jpt)=PromSpecDat(jpt)+resSp(jpt)
     &                        *fparms(IAMPFAC,isi,isp)*keySp(isp,2)
                  enddo
cdb                   do jpt=1,npts(isp)
cdb                           jtmp=ixs+jpt-1
cdb                 fjac1(jtmp)=fjac1(jtmp)+(ClSpecDat(jpt,isp,isi+1)-
cdb     #                  resSp(jpt)*fparms(IAMPFAC,isi,isp)*keySp(isp,2))
cdb     #                             *(1.d0+shtraftmp)     
cdb                   end do
cdb        continue          
              else  ! not depOK
                  do jpt=1,npts(isp)
                    PromSpecDat(jpt)=PromSpecDat(jpt)+
     &	                             ClSpecDat(jpt,isp,isi+1)
                  enddo              
              endif !depOK
            end do	!*** end of loop over sites in spectrum
c           ------------------------------------------
c            Add difference spectrum
c           ------------------------------------------
                 if (numODFs.gt.0) shtraftmp=shtraftmp/numODFs
                 if ((iwflag.ge.1 .and. iwflag.le.3) .or. (iwflag.eq.5)) then
c                                               --- weighted ---
                   do jpt=1,npts(isp)
                           jtmp=ixs+jpt-1
                 fjac(jtmp,ix)=fjac(jtmp,ix)+((ClSpecDat(jpt,isp,1)
     #                -PromSpecDat(jpt))*(1.d0+shtraftmp)
     #        +(fitdata(jtmp)-PromSpecDat(jpt))*(shtraftmp-shtraf(isp)))
     #                 /max(dabs(rmsn(isp)),RNDOFF)
                 if (haveWgt) fjac(jtmp,ix)=fjac(jtmp,ix)*wgt(jtmp)
                   end do
c                                               --- unweighted ---
                 else
                   do jpt=1,npts(isp)
                      jtmp=ixs+jpt-1
                 fjac(jtmp,ix)=fjac(jtmp,ix)+(ClSpecDat(jpt,isp,1)
     #                -PromSpecDat(jpt))*(1.d0+shtraftmp)
     #        +(fitdata(jtmp)-PromSpecDat(jpt))*(shtraftmp-shtraf(isp))
                 if (haveWgt) fjac(jtmp,ix)=fjac(jtmp,ix)*wgt(jtmp)
                   end do
                 end if	!iwflag  
cdb       continue                                                   
          end do    !*** end of loop over spectra
c           -------------------------------------------------------  
c            Use the step size to calculate the partial derivative 
c            w.r.t. xfit(ix) for the current spectrum
c           -------------------------------------------------------  
          do jpt=1,ndatot
               fjac(jpt,ix)=fjac(jpt,ix)/xfdstp(ix)
          end do
          xfit(ix)=xtemp
        end do !*** end of loop over parameters
        continue         
c**********************************************************************
c     IFLAG=3
c     Fit has terminated: print result 
c**********************************************************************
      else if (abs(iflag).eq.3) then
         call xpack(2,itemp,xfit)
c              Final reporting uses metrics populated by the latest
c              iflag=1 pass, after all objective residual normalization.
         if (mpi_rank.eq.0) then
            call LossReport(6)
            call LossReport(2)
         endif
c
         if (Kcontrol.eq.1) then
                if (mpi_rank.eq.0) write(*,*)
     #          ' Fitting is terminated by user'
                if (mpi_rank.eq.0) write(2,*)
     #          ' Fitting is terminated by user'
                iflag=-1
         endif  
      end if
c
      return
      end

c**********************************************************************
c     Print loss metrics collected from the latest LFUN iflag=1 pass.
c**********************************************************************
      subroutine LossReport(lu)
      implicit none
      integer lu,isp
      double precision frac,redw,redun,sumw,sumun,tol
      include 'dimnls.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'iterat.inc'
      include 'rndoff.inc'
      double precision DndatNprm
      common /DndatNprm/ DndatNprm

      write(lu,900) lossSseW,lossSseUn,lossRedW,lossRedUn,
     #              ndatot,nprm,DndatNprm
      write(lu,901)

      sumw=0.d0
      sumun=0.d0
      do isp=1,nspc
         frac=0.d0
         redw=0.d0
         redun=0.d0
         if (dabs(lossSseW).gt.RNDOFF) frac=lossSpW(isp)/lossSseW
         if (npts(isp).gt.0) then
            redw=lossSpW(isp)/dfloat(npts(isp))
            redun=lossSpUn(isp)/dfloat(npts(isp))
         endif
         write(lu,902) isp,npts(isp),lossSpW(isp),lossSpUn(isp),
     #                 redw,redun,frac
         sumw=sumw+lossSpW(isp)
         sumun=sumun+lossSpUn(isp)
      enddo

      tol=1.d-8*max(1.d0,dabs(lossSseW),dabs(lossSseUn))
      if ((dabs(sumw-lossSseW).gt.tol).or.
     #    (dabs(sumun-lossSseUn).gt.tol)) then
         write(lu,903)
      endif

 900  format(/1x,'Final loss reporting:',/,
     #       3x,'SSE_W=',1PE13.6,2x,'SSE_UNW=',1PE13.6,/,
     #       3x,'RED_W=',1PE13.6,2x,'RED_UNW=',1PE13.6,/,
     #       3x,'N_TOTAL=',I8,2x,'NPRM=',I8,2x,'DOF=',1PE13.6)
 901  format(1x,'Spectrum contributions:',/,
     #       3x,'isp',3x,'npts',7x,'SSE_W',9x,'SSE_UNW',
     #       9x,'RED_W',10x,'RED_UNW',7x,'FRACTION_W')
 902  format(3x,I3,2x,I6,2x,1PE13.6,2x,1PE13.6,2x,
     #       1PE13.6,2x,1PE13.6,2x,1PE11.4)
 903  format(1x,'WARNING: spectrum contributions do not sum to total loss')

      return
      end

c#######################################################################
