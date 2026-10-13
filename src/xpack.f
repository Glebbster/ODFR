cAXV for ODFR 01.12.21014
c This program pack the varied parameters in the vector xfit and unpack
c the vector to set of parameters.
c   inout=1 - pack parameters from fparm array and parameter array
c             into xfit vector
c   inout=2 - back transfer
c  nprm is number of parameters for fitting  (in common block)
c  nprm=nprm1+nprm2
c  nprm1 is number of parameters which are varied for all spectra together
c  nprm2 is number of parameters which are varied independently for distinct spectrum
c  nprm3 is number of parameters which are not varied but are changed
c  nprm4 is the number ODF parameters being varied
c order of parameters in xfit-array is nprm1, nprm2, nprm4;
c nprm3 indirectly changing parameters are not keeping in xfit-array 
c but their properties are stored in arrays ixpr, ixst1,ixst2,ixx with indexes 
c after the same properties for varied parameters 
c----------------------------------------------------------------------
      subroutine xpack(inout,JVARout,xfit)
      implicit none
c
      include 'dimnls.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'rndoff.inc'
      include 'parODF.inc'

      double precision xfit(MXJCOL)
	double precision MagPar0
      common /magpar0/magpar0(6,MXSITE)
c   For preservation the values in course of
c   consecutive calling the subroutine lfunR.f:     	
      double precision angfi0,angpsi0,angksi0
      double precision gav, Aav,Dav,ParRatio      
      common /tempxpack/ angfi0(MXSPC),angpsi0(MXSPC),angksi0(MXSPC),
     #      gav(MXSPC,MXSITE), Aav(MXSPC,MXSITE),Dav(MXSPC,MXSITE),
     #      ParRatio(NFPRM)
c
      integer inout,JVAR,JVARout,jsiter,ii,jspec,js,ipar
      integer iflvar,filenprm,iparms(NIPRM,MXSITE,MXSPC),iRatio
      DOUBLE PRECISION fparms(NFPRM,MXSITE,MXSPC),prom
      double precision SAFELOGMIN,SAFEEXPMAX
      parameter (SAFELOGMIN=1.0D0,SAFEEXPMAX=13.0D0)
      COMMON /COMPARMS/ iparms,fparms
      common /iflvar/iflvar(NFPRM,MXSITE,MXSPC),filenprm(MXSPC)
      integer modtd
      logical depOK,useMOMD
        common/depOK/ depOK(MXVAR,MXSITE,MXSPC),modtd(MXSITE,MXSPC),
     #                useMOMD(MXSITE,MXSPC)        
      integer Naver,keyiparm,nInSite 
       common /Naver/Naver,keyiparm(MXSPC)     
        integer IEQs,Kcontrol
        COMMON /IEQ/ IEQs(MaxAddNucl,MXSITE,MXSPC)
        common /control/Kcontrol
c        
c ########################################################################
c
      SELECT CASE (inout)
	CASE (1)  ! pack varied parameters in xfit
        JVAR=0
        nprm1=0
        nprm2=0
        nprm3=0
c
c  for librations
        do js=1,nnsite(1)
          if (iparms(IKEYLIB,js,1).gt.0) then
            do ipar=0,5
              MagPar0(ipar+1,js)=fparms(IGXX+ipar,js,1)
            enddo	  !ipar
            call librat(iparms(IKEYLIB,js,1),fparms(IANGLIB,js,1),
     #                  MagPar0(1,js),fparms(IGXX,js,1))
            do jspec=2,nspc
              do ipar=0,5
                fparms(IGXX+ipar,js,jspec)=fparms(IGXX+ipar,js,1)
              enddo !ipar
            enddo   !jspec
          endif     !iparms(IKEYLIB,js,1)
        enddo !js
c treatment of flags 6, 7 and 8                 
        do js=1,nnsite(1)		!loop over sites
          do ipar=1,NFPRM	!loop over parameters
            if (iflvar(ipar,js,1).eq.6) then
               iflvar(ipar,js,1)=1
               do jspec=2,nspc
                 if ((iflvar(ipar,js,jspec).ne.11).and.
     #              (iflvar(ipar,js,jspec).ne.12)) then
                     filenprm(jspec)=filenprm(jspec)+1
                     nprm = nprm+1
                     iflvar(ipar,js,jspec)=11                     
                 endif                         
               enddo
            endif 
            if (iflvar(ipar,js,1).eq.7) then
               iflvar(ipar,js,1)=1
               do jspec=2,nspc
                if ((iflvar(ipar,js,jspec).eq.11).or.
     #              (iflvar(ipar,js,jspec).eq.12)) then
                    filenprm(jspec)=filenprm(jspec)-1
                    nprm = nprm-1
                 endif 
                 iflvar(ipar,js,jspec)=0                   
               enddo
            endif 
            if (iflvar(ipar,js,1).eq.8) then
               iflvar(ipar,js,1)=0
               do jspec=2,nspc
                if ((iflvar(ipar,js,jspec).eq.11).or.
     #              (iflvar(ipar,js,jspec).eq.12)) then
                    filenprm(jspec)=filenprm(jspec)-1
                    nprm = nprm-1
                 endif 
                 iflvar(ipar,js,jspec)=10                   
               enddo
            endif             
          enddo
        enddo                
        if (Kcontrol.eq.2) go to 100
c          
c    search independent variables for first spectrum 
        do js=1,nnsite(1)		!loop over sites
          do ipar=1,NFPRM	!loop over parameters
c                 
		    SELECT CASE (iflvar(ipar,js,1))
            CASE (1) ! varied directly
             JVAR=JVAR+1
             ixpr(JVAR)=ipar
             ixst1(JVAR)=1
             ixst2(JVAR)=js
             ixx(jvar)=iflvar(ipar,js,1)
		      SELECT CASE (ipar)
c               CASE (IANGFI) 
c                do jspec=1,nspc
c			     angfi0(jspec)=fparms(ipar,js,jspec)
c	          enddo
c                xfit(JVAR)=0.d0    
              CASE (IPSI) 
                do jspec=1,nspc
			     angpsi0(jspec)=fparms(ipar,js,jspec)
                 iparms(IIPSI0,js,jspec)=1			
	          enddo
                xfit(JVAR)=0.d0 
              CASE (IANGKSI) 
c                do jspec=1,nspc
c			     angksi0(jspec)=fparms(ipar,js,jspec)
c	          enddo
c                xfit(JVAR)=0.d0 
              CASE DEFAULT    			   			   			   
                xfit(JVAR)=fparms(ipar,js,1) 
              END SELECT			      
            CASE (2)  !	varied as transformed
             JVAR=JVAR+1	
             ixpr(JVAR)=ipar
             ixst1(JVAR)=1
             ixst2(JVAR)=js
             ixx(jvar)=iflvar(ipar,js,1)
             SELECT CASE (ipar)
             CASE(IGXX:IGZZ)  !g-values
               xfit(JVAR)=(fparms(ipar,js,1) -2.D0)*1000.D0
c          other transformations is inserted here 
             CASE(IDX:IDZ,IDJF:IOSS)  !diffusion coefficients and Heisenberg exchange
	               xfit(JVAR)=dlog10(dmax1(fparms(ipar,js,1),SAFELOGMIN))
c          other transformations is inserted here 
c            CASE  (?)
c		     CASE DEFAULT 
c
             END SELECT !ipar
            END SELECT	!iflvar(ipar,js,1) for first spectrum
          enddo !ipar
        enddo !js
c    search independent variables for other spectra 
c    flag=11,12
        nprm1=JVAR
        do jspec=1,nspc				!loop over spectra
          do js=1,nnsite(jspec)		!loop over sites
            do ipar=1,NFPRM			!loop over parameters
		      if (iflvar(ipar,js,jspec).EQ.11) then
                nprm2=nprm2+1				
                ixpr(nprm1+nprm2)=ipar
                ixst1(nprm1+nprm2)=jspec
                ixst2(nprm1+nprm2)=js
                ixx(nprm1+nprm2)=iflvar(ipar,js,jspec)		    
                xfit(nprm1+nprm2)=fparms(ipar,js,jspec)
              endif
		      if (iflvar(ipar,js,jspec).EQ.12) then
                nprm2=nprm2+1				
                ixpr(nprm1+nprm2)=ipar
                ixst1(nprm1+nprm2)=jspec
                ixst2(nprm1+nprm2)=js
                ixx(nprm1+nprm2)=iflvar(ipar,js,jspec)
                SELECT CASE (ipar)
                CASE(IGXX:IGZZ)  !g-values
                  xfit(nprm1+nprm2)=(fparms(ipar,js,jspec) -2.D0)
     #                             *1000.D0
c          other transformations is inserted here 
                CASE(IDX:IDZ,IDJF:IOSS)  
c !diffusion coefficients and Heisenberg exchange                
                  xfit(nprm1+nprm2)=
     #                 dlog10(dmax1(fparms(ipar,js,jspec),SAFELOGMIN))
			    END SELECT					    
              endif
            enddo !ipar
          enddo !js
        enddo !jspec
        JVAR=nprm1+nprm2
c  search the parameters Of ODF being varied for the sites	(they are the same for all spectra!)
        do js=1,nnsite(1)		!loop over sites
          if ((keyODF(js).eq.0).and.(ODFcontrol.lt.2)) then     
            do ipar=1,nprmODF(js)
		      if (iflagODF(ipar,js).gt.0)then
                JVAR=JVAR+1	
                ixpr(JVAR)=ipar
                ixst1(JVAR)=1	 !! for all cpectra!!
                ixst2(JVAR)=js
c value ixx(jvar) is negative for varied ODF parameters!              
                ixx(jvar)=-iflagODF(ipar,js)		    
                xfit(JVAR)=paramODF(ipar,js) 
              endif
            enddo !ipar
          endif !(keyODF(js).eq.0)          
        enddo  !js				           		  
        JVARout=JVAR
c for all spectrum, search the parameters which are changed in course of variation of others 
        iRatio=0
        do jspec=1,nspc				!loop over spectra
          do js=1,nnsite(jspec)		!loop over sites
            do ipar=1,NFPRM			!loop over parameters
		      SELECT CASE (iflvar(ipar,js,jspec))
              CASE (3,13)  !	averaged tenzors
               nprm3=nprm3+1
               ixpr(JVAR+nprm3)=ipar
               ixst1(JVAR+nprm3)=jspec
               ixst2(JVAR+nprm3)=js
               ixx(JVAR+nprm3)=iflvar(ipar,js,jspec)
               SELECT CASE (ipar)
               CASE(IGXX:IGZZ)  !averaged g-values.  
                 gav(jspec,js)=(fparms(IGXX,js,jspec)+ 
     *               fparms(IGYY,js,jspec)+fparms(IGZZ,js,jspec))    
               CASE(IAXX:IAZZ)  !averaged A-values
                 Aav(jspec,js)=(fparms(IAXX,js,jspec)+ 
     *               fparms(IAYY,js,jspec)+fparms(IAZZ,js,jspec))
               CASE(IDX:IDZ)  !averaged D-values	 
                 Dav(jspec,js)=(fparms(IDX,js,jspec)+ 
     *               fparms(IDY,js,jspec)+fparms(IDZ,js,jspec))
               CASE(IANGLIB+3:IANGLIB+5)!Libration frame=Rotatiion frame
                     fparms(ipar,js,jspec)=fparms(ipar-32,js,jspec)
c               CASE (IANGFI) 
c			     angfi0(jspec)=fparms(ipar,js,jspec)
               CASE (IPSI) 
			     angpsi0(jspec)=fparms(ipar,js,jspec)
c               CASE (IANGKSI) 
c			     angksi0(jspec)=fparms(ipar,js,jspec)
               END SELECT !ipar
c               
              CASE (4,14)  !	axial tenzors
               nprm3=nprm3+1
               ixpr(JVAR+nprm3)=ipar
               ixst1(JVAR+nprm3)=jspec
               ixst2(JVAR+nprm3)=js
               ixx(JVAR+nprm3)=iflvar(ipar,js,jspec)
               if (((ipar.eq.IGYY).or.(ipar.eq.IAYY)).or.(ipar.eq.IDY))
     #            fparms(ipar,js,jspec)=fparms(ipar-1,js,jspec)
              CASE (5,15)  !	isotropic tenzors
               nprm3=nprm3+1
               ixpr(JVAR+nprm3)=ipar
               ixst1(JVAR+nprm3)=jspec
               ixst2(JVAR+nprm3)=js
               ixx(JVAR+nprm3)=iflvar(ipar,js,jspec)
               if (ipar.eq.IGYY.or.ipar.eq.IAYY
     #            .or.ipar.eq.IDY.or.ipar.eq.IWYY)
     #            fparms(ipar,js,jspec)=fparms(ipar-1,js,jspec)
               if (ipar.eq.IGZZ.or.ipar.eq.IAZZ
     #            .or.ipar.eq.IDZ.or.ipar.eq.IWZZ)
     #            fparms(ipar,js,jspec)=fparms(ipar-2,js,jspec)
c
              CASE (91:95)  ! preparation of ratious for sites 
               nprm3=nprm3+1
               ixpr(JVAR+nprm3)=ipar
               ixst1(JVAR+nprm3)=jspec
               ixst2(JVAR+nprm3)=js
               ixx(JVAR+nprm3)=iflvar(ipar,js,jspec)
               jsiter=iflvar(ipar,js,jspec)-90
               if (dabs(fparms(ipar,jsiter,jspec)).lt.rndoff) then
                 WRITE (2,*)  ' inputted flag=',iflvar(ipar,js,jspec),
     #  ' is impossible, as base parameter is equal zero'
                 WRITE (*,*)  ' inputted flag=',iflvar(ipar,js,jspec),
     #  ' is impossible, as base parameter is equal zero'
                 stop
               endif               
               iRatio=iRatio+1
               ParRatio(iRatio)=fparms(ipar,js,jspec)
     #                     /fparms(ipar,jsiter,jspec)
              END SELECT	!iflvar(ipar,js,jspec) for first spectrum
            enddo !ipar
          enddo !js
        enddo !jspec
c copy the parameters of first spectrum to all other
c excepting all independent parameters c flag=10,11,12,13,90...
c and parameters ODF
        do JVAR=1,nprm1+nprm2
          do jspec=2,nspc
            if   ((ixx(jvar).ge.10) .or.(ixpr(JVAR).eq.IPSI).or.
c     #          (ixpr(JVAR).eq.IANGFI).or.(ixpr(JVAR).eq.IANGKSI).or.   
     #          (iflvar(ixpr(JVAR),ixst2(JVAR),jspec).ge.10)) then	 
                fparms(ixpr(JVAR),ixst2(JVAR),jspec)=
     #                            fparms(ixpr(JVAR),ixst2(JVAR),1)
            endif
          enddo !jspec
        enddo !JVAR
c copy indirectly changed values from first to all other spectra        
                do JVAR=JVARout+1,JVARout+nprm3
          do jspec=2,nspc
            if   ((ixx(jvar).ge.10) .or.(ixpr(JVAR).eq.IPSI).or.
c     #          (ixpr(JVAR).eq.IANGFI).or.(ixpr(JVAR).eq.IANGKSI).or.   
     #          (iflvar(ixpr(JVAR),ixst2(JVAR),jspec).ge.10)) then	 
                fparms(ixpr(JVAR),ixst2(JVAR),jspec)=
     #                            fparms(ixpr(JVAR),ixst2(JVAR),1)
            endif
          enddo !jspec
        enddo !JVAR
c mark the spectrum/site dependent on the parameter
        do JVAR=1,JVARout          
          depOK(JVAR,ixst2(JVAR),ixst1(JVAR))=.true.
c mark all spectra when first (master) spectrum is dependent or if it is ODF parameter:                             
          if (ixst1(JVAR).eq.1) then
	       do jspec=2,nspc
                depOK(JVAR,ixst2(JVAR),jspec)=
     #             (iflvar(abs(ixx(JVAR)),ixst2(JVAR),jspec).lt.10)
	       enddo
          endif
          if (ixx(JVAR).lt.0) then
	       do jspec=2,nspc
                depOK(JVAR,ixst2(JVAR),jspec)=.true.
 	       enddo
          endif          
        enddo 
c  mark as dependent spectra with indirect dependence      
        do JVAR=JVARout+1,JVARout+nprm3
          if (iflvar(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR)).gt.90) then
             jsiter=ixx(JVAR)-90
             do ipar=1,nprm1+nprm2 !the same dependences as master parameter
               if (ixpr(ipar).eq.ixpr(JVAR)) then
                   depOK(ipar,ixst2(JVAR),ixst1(JVAR))=
     #                      depOK(ipar,jsiter,ixst1(JVAR))
                  if (ixst1(JVAR).eq.1) then         
	              do jspec=2,nspc
                     depOK(ipar,ixst2(JVAR),jspec)=
     #                      depOK(ipar,jsiter,1)
                    enddo !jspec
                 endif   !ixst1(JVAR).eq.1
               endif     !(ixpr(ipar).eq.ixpr(JVAR))                    
             enddo
          endif
        enddo
c detect is the spectrum/site can be calculated using iort spectra for first file:  
        do jspec=2,nspc 
          do js=1,nnsite(jspec)      
            if ((keyiparm(jspec).eq.0).and.
     #                (iparms(INORT,js,jspec).gt.1))then 
              useMOMD(js,jspec)=.true.     
              do ipar=IWXX,IOSS   !5,22
                 useMOMD(js,jspec)=(useMOMD(js,jspec).and.
     #                (iflvar(ipar,js,jspec).lt.10).and.
     #                 (iflvar(ipar,js,1).lt.10)) 
              enddo              
              do ipar=IALD,IC44   !23,33
                 useMOMD(js,jspec)=(useMOMD(js,jspec).and.
     #                (iflvar(ipar,js,jspec).lt.10).and.
     #                 (iflvar(ipar,js,1).lt.10)) 
              enddo                         
              do ipar=IEXPL,Ireserved  ! 49,57
                 useMOMD(js,jspec)=(useMOMD(js,jspec).and.
     #                (iflvar(ipar,js,jspec).lt.10).and.
     #                 (iflvar(ipar,js,1).lt.10)) 
              enddo              
            endif !keyiparm & INORT
          enddo  !js
        enddo    !jspec
        if (JVARout+nprm3.gt.MXVAR) then
          OPEN ( UNIT=2,FILE='history', POSITION='APPEND')
            write(2,*)' Number of direct and undirect variable',
     #            JVARout+nprm3,' more then MXVAR= ',MXVAR
            write(*,*)' Number of direct and undirect variable',
     #            JVARout+nprm3,' more then MXVAR= ',MXVAR       
          CLOSE ( UNIT=2)
 	    stop
 	  endif          
 100    continue
c              
c__________________________________________________________________
c  
	CASE (2) !inout=2, transfer the parametrrs from xfit to fpars
	!loop over parameters varied directly for all spectra together
        do JVAR=1,nprm1
          SELECT CASE (ixx(JVAR))	! 
          CASE (1) ! varied directly
            SELECT CASE (ixpr(JVAR))
c            CASE (IANGFI) 
c              do jspec=1,nspc
c                 if (iflvar(IANGFI,1,jspec).NE.10)              
c     #               fparms(IANGFI,1,jspec)=angfi0(jspec)+xfit(JVAR)
c	        enddo
            CASE (IPSI) 
              do jspec=1,nspc
                 if (iflvar(IPSI,1,jspec).NE.10)   
     #               fparms(IPSI,1,jspec)=angpsi0(jspec)+xfit(JVAR)
	        enddo
c            CASE (IANGKSI)
c              do jspec=1,nspc
c                 if (iflvar(IANGKSI,1,jspec).NE.10)   
c     #               fparms(IANGKSI,1,jspec)=angksi0(jspec)+xfit(JVAR)
c	        enddo
            CASE (IANGLIB:IANGLIB+5)
              fparms(ixpr(JVAR),ixst2(JVAR),1)=xfit(JVAR)
              call librat(iparms(IKEYLIB,ixst2(JVAR),1),
     #          fparms(IANGLIB,ixst2(JVAR),1),MagPar0(1,ixst2(JVAR)),
     #          fparms(IGXX,ixst2(JVAR),1))
              do jspec=2,nspc
                do ipar=0,5
                 if (iflvar(IGXX+ipar,ixst2(JVAR),jspec).LT.10)  
     #              fparms(IGXX+ipar,ixst2(JVAR),jspec)=
     #              fparms(IGXX+ipar,ixst2(JVAR),1)
                enddo !ipar
              enddo   !jspec                  
            CASE (IGXX:IAZZ)
c libration case
              if (iparms(IKEYLIB,ixst2(JVAR),1).gt.0) then
                 magpar0(ixpr(JVAR)-IGXX+1,ixst2(JVAR))=xfit(JVAR)
                 call librat(iparms(IKEYLIB,ixst2(JVAR),1),
     #                fparms(IANGLIB,ixst2(JVAR),1),
     #                MagPar0(1,ixst2(JVAR)),fparms(IGXX,ixst2(JVAR),1))
                 do jspec=2,nspc
                   do ipar=0,5
                     fparms(IGXX+ipar,ixst2(JVAR),jspec)=
     #                    fparms(IGXX+ipar,ixst2(JVAR),1)
                   enddo !ipar
                 enddo   !jspec
              else   
                 fparms(ixpr(JVAR),ixst2(JVAR),1)=xfit(JVAR) 
              endif     !iparms(IKEYLIB,ixst2(JVAR),1)  
            CASE DEFAULT		    
                 fparms(ixpr(JVAR),ixst2(JVAR),1)=xfit(JVAR)
            END SELECT !ixpr(JVAR)		          
          CASE (2)  !	variation of transformed pararameter
            SELECT CASE (ixpr(JVAR))
            CASE(IGXX:IGZZ)  !g-values
c libration case
              if (iparms(IKEYLIB,ixst2(JVAR),1).gt.0) then
                 magpar0(ixpr(JVAR)-IGXX+1,ixst2(JVAR))=
     #                          2.D0+ xfit(JVAR)/1000.D0
                 call librat(iparms(IKEYLIB,ixst2(JVAR),1),
     #                fparms(IANGLIB,ixst2(JVAR),1),
     #                MagPar0(1,ixst2(JVAR)),fparms(IGXX,ixst2(JVAR),1))
                 do jspec=2,nspc
                   do ipar=0,2
                     fparms(IGXX+ipar,ixst2(JVAR),jspec)=
     #                     fparms(IGXX+ipar,ixst2(JVAR),1)
                   enddo !ipar
                 enddo   !jspec
              else   
               fparms(ixpr(JVAR),ixst2(JVAR),1)=2.D0+ xfit(JVAR)/1000.D0
              endif     !iparms(IKEYLIB,ixst2(JVAR),1)  
            CASE(IDX:IDZ,IOSS)  !diffusion coefficients, Heisenberg exchange
                prom=xfit(JVAR)
	                if (prom.gt.SAFEEXPMAX) then
                   fparms(ixpr(JVAR),ixst2(JVAR),1)=1.D13
                else
                   if	(prom.lt.1.D0) then
                   fparms(ixpr(JVAR),ixst2(JVAR),1)=10.D0
                   else
                   fparms(ixpr(JVAR),ixst2(JVAR),1)=10.D0**(prom)
                   endif 
                endif
c
c          other transformations is inserted here 
c            CASE  (?)
c            fparms(??) = xfit(JVAR)??
c		  CASE DEFAULT 
c           fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR)= xfit(JVAR)
            END SELECT !ixpr(JVAR)
          END SELECT ! (ixx(JVAR)
        enddo !JVAR
c  transfer the parametrrs from xfit to paramODF
        do JVAR=nprm1+nprm2+1,nprm
              paramODF(ixpr(JVAR),ixst2(JVAR))=xfit(JVAR) 
        enddo 
        JVARout=JVAR
c
        iRatio=0
        do JVAR=nprm+1,nprm+nprm3		!loop over calculated parameters    
          SELECT CASE (ixx(JVAR))	! 
	    CASE (3,13)  !from averaged tenzors
            SELECT CASE (ixpr(JVAR))
            CASE(IGXX:IGZZ)  !calculation g-value using gav 
            fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=
     #                               gav(ixst1(JVAR),ixst2(JVAR)) 
            do ii=IGXX,IGZZ
              if (ii.NE.ixpr(JVAR))
     #        fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))= 
     #                       fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))-
     #                  	   fparms(ii,ixst2(JVAR),ixst1(JVAR))
	      enddo
c
            CASE(IAXX:IAZZ)  !calculation A-value from Aav
            fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=
     #                               Aav(ixst1(JVAR),ixst2(JVAR)) 
            do ii=IAXX,IAZZ
              if (ii.NE.ixpr(JVAR)) 
     #           fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))= 
     #                       fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))-
     #                       fparms(ii,ixst2(JVAR),ixst1(JVAR))
	      enddo
c
            CASE(IDX:IDZ)  !averaged D-values
            fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))= 
     #                               Dav(ixst1(JVAR),ixst2(JVAR)) 
            do ii=IDX,IDZ
              if (ii.NE.ixpr(JVAR))
     #           fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))= 
     #                       fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))-
     #                  	   fparms(ii,ixst2(JVAR),ixst1(JVAR))
	      enddo
            CASE(IANGLIB+3:IANGLIB+5)  !Libration frame=Rotatiion frame	 
            fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=
     #         fparms(ixpr(JVAR)-32,ixst2(JVAR),ixst1(JVAR))
c            CASE (IANGFI) 
c			fparms(IANGFI,ixst2(JVAR),ixst1(JVAR))=angfi0(ixst1(JVAR))+
c     #			fparms(IANGFI,ixst2(JVAR),1)-angfi0(1)
            CASE (IPSI) 
			fparms(IPSI,ixst2(JVAR),ixst1(JVAR))=angpsi0(ixst1(JVAR))+
     #			fparms(IPSI,ixst2(JVAR),1)-angpsi0(1)  
c            CASE (IANGKSI) 
c			fparms(IANGKSI,ixst2(JVAR),ixst1(JVAR))=angksi0(ixst1(JVAR))+
c     #			fparms(IANGKSI,ixst2(JVAR),1)-angksi0(1)         
            END SELECT !ixpr(JVAR)
          CASE (4,14)  !	axial tenzors
            fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=
     #        fparms(ixpr(JVAR)-1,ixst2(JVAR),ixst1(JVAR))
c
          CASE (5,15)  !	isotropic tenzors
            if (ixpr(JVAR).eq.IGYY.or.ixpr(JVAR).eq.IAYY
     #	    .or.ixpr(JVAR).eq.IDY.or.ixpr(JVAR).eq.IWYY) 
     #        fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=
     #           fparms(ixpr(JVAR)-1,ixst2(JVAR),ixst1(JVAR))
            if (ixpr(JVAR).eq.IGZZ.or.ixpr(JVAR).eq.IAZZ
     #	    .or.ixpr(JVAR).eq.IDZ.or.ixpr(JVAR).eq.IWZZ)
     #        fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=
     #           fparms(ixpr(JVAR)-2,ixst2(JVAR),ixst1(JVAR))
c
          CASE (91:95)  ! calculation values using ratious for sites
            jsiter=ixx(JVAR)-90
            iRatio=iRatio+1
            fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=
     #         fparms(ixpr(JVAR),jsiter,1)*ParRatio(iRatio)
c
          END SELECT !ixx(JVAR)
        enddo	!JVAR=nprm+1,nprm+nprm3, loop over calculated parameters
c copy the parameters of first spectrum to all other
        do JVAR=1,nprm+nprm3
          do jspec=2,nspc
            do js=1,nnsite(jspec)
              if  (.not.((ixx(jvar).lt.0) .or.
     #            (ixx(jvar).ge.10) .or.(ixpr(JVAR).eq.IPSI).or.
c     #            (ixpr(JVAR).eq.IANGFI).or.(ixpr(JVAR).eq.IANGKSI).or.
     #            (iflvar(ixpr(JVAR),js,jspec).ge.10))) then
                    fparms(ixpr(JVAR),js,jspec)=
     #                              fparms(ixpr(JVAR),js,1)
              endif
            enddo !js              
          enddo !jspec
        enddo !JVAR
c   additional nuclea
        do jspec=1,nspc
          do js=1,nnsite(jspec)
            do ii=0,iparms(iNaddN,js,jspec)-1
             if (IEQs(ii+1,js,jspec).gt.0)then
              fparms(IAadd+2*ii,js,jspec)=
     #            fparms(IAadd+2*(IEQs(ii+1,js,1)-1),js,jspec)
              fparms(IAadd+2*ii+1,js,jspec)=
     #            fparms(IAadd+2*IEQs(ii+1,js,1)-1,js,jspec)
             endif
            enddo
          enddo !js
        enddo !jspec
c assign the independent variable for distinct spectrum
        do JVAR=nprm1+1,nprm1+nprm2
          if (ixx(JVAR).EQ.11) then
              fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=xfit(JVAR)
          endif
          if (ixx(JVAR).EQ.12) then
                SELECT CASE (ixpr(JVAR))
                CASE(IGXX:IGZZ)  !g-values
                 fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=
     #                               2.D0+ xfit(JVAR)/1000.D0
                CASE(IDX:IDZ,IOSS) !diffusion coefficients,Heisenberg exchange
                prom=xfit(JVAR)
	                if (prom.gt.SAFEEXPMAX) then
                   fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=10.D12
                else
                   if	(prom.lt.1.D0) then
                   fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=10.D0
                   else
                fparms(ixpr(JVAR),ixst2(JVAR),ixst1(JVAR))=10.D0**(prom)
                   endif 
                endif
                END SELECT
          endif
        enddo !JVAR
c  the common parameters for sites:
        do jspec=1,nspc
          do js=2,nnsite(jspec)
c            fparms(IANGFI,js,jspec)=fparms(IANGFI,1,jspec)
            fparms(IPSI,js,jspec)=fparms(IPSI,1,jspec)
c            fparms(IANGKSI,js,jspec)=fparms(IANGKSI,1,jspec)
            fparms(ISPSHFT,js,jspec)=fparms(ISPSHFT,1,jspec)
c the same content for additional sites 
            if (ixx(JVAR).EQ.0)          
     #       fparms(IAMPFAC,js,jspec)=fparms(IAMPFAC,1,jspec) 
     #           *fparms(IAMPFAC,js,1)/fparms(IAMPFAC,1,1)
c  the case of sites which differs by additional HFI only:
            if (iparms(IItypeSite,js,jspec).gt.0) then
              nInSite=iparms(IItypeSite,js,jspec)          
              fparms(1,js,jspec)=fparms(1,nInSite,jspec)
              do ipar=IWXX,IDFLD    !5,45 
                 fparms(ipar,js,jspec)=fparms(ipar,nInSite,jspec)
              enddo  !ipar  
              do ipar=ISPSHFT,Ireserved     !47,59 
                 fparms(ipar,js,jspec)=fparms(ipar,nInSite,jspec)
              enddo  !ipar                            
            endif !iparms(IItypeSite,js,jspec)
          enddo  !js
        enddo !jspec
c  identical ODF for sites:
        do js=2,nnsite(1)		!loop over sites
          if (keyODF(js).gt.0) then     
            do ipar=1,nprmODF(js)
              paramODF(ipar,js)=paramODF(ipar,1) 
            enddo !ipar
          endif !(keyODF(js).eq.0)          
        enddo  !js	      
c	  		   
	END SELECT	!INOUT
      return
      end
