C for ODFR 01/12/2014
	SUBROUTINE  PARREADN(FileErr,NAME,jfile)
c
      implicit none
      include 'dimnls.inc'
      include 'rndoff.inc'
      include 'eprprm.inc'
      include 'parcom.inc'
      include 'prmeqv.inc'
      include 'expdat.inc'
      include 'lmcom.inc'
      include 'iterat.inc'
      include 'errmsg.inc'
      include 'mpi_state.inc'
c
      double precision ZERO,HALF,ONE,TEN
      integer ierr
      parameter (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TEN=1.0D1)
c 
        integer iparms(NIPRM,MXSITE,MXSPC)
        double precision fparms(NFPRM,MXSITE,MXSPC)
        double precision gmax,gmin,hmax
        COMMON /COMPARMS/ iparms,fparms
        integer IEQs
        COMMON /IEQ/ IEQs(MaxAddNucl,MXSITE,MXSPC)
c
        integer Kcontrol,js,i,iprom
	  integer iflvar, filenprm
        common /iflvar/iflvar(NFPRM,MXSITE,MXSPC),filenprm(MXSPC)
c
      	 Logical FileErr,lpr
         logical axiala,axialg,axialw,axialr
c      	 
         common /control/Kcontrol
c
	character*30 str1,NAME
      integer jfile, jpar,ICON,ISTRING
      double precision parReserved,temp1,temp2,temp3,fldf
c
      integer ntemp,ntemp1,Naver,keyiparm,nInSite 
       common /Naver/Naver,keyiparm(MXSPC)
      integer ipar
       integer keyCov 
       common/KeyCov/keyCov		          
c      
      external ipar
cdb for information only
cdb  IPHASE=1,IGIB0=2,IGIB2=3,IGIBX=4,IWXX=5,IWYY=6,IWZZ=7,IGXX=8,
cdb  IGYY=9,IGZZ=10,IAXX=11,IAYY=12,IAZZ=13,IDX=14,IDY=15,
cdb  IDZ=16,IPML=17,IPMXY=18,IPMZZ=19,IDJF=20,IDJFPRP=21,
cdb  IOSS=22,IPSI=23,IALD=24,IBED=25,IGAD=26,IALM=27,IBEM=28,
cdb  IGAM=29,IC20=30,IC44=34,ILBL=35,ILBG=36,IDC20=37,IB0=38,
cdb  IGAMAN=39,ICGTOL=40,ISHIFTR=41,IRANGE=43,IFLDI=44,
cdb  IDFLD=45,IAMPFAC=46,ISPSHFT=47,IANGFI=48,IANGKSI=49,
cdb  IEXPL=50,IEXPKXY=51,IEXPKZZ=52,IANGLIB=53,Ireserved=59,
cdb  IAadd=60)
cbg     parameters for iparms(*,js,jfile)    
cdb     IIN2=1,IIPDF=2,IIST=3,IML=4,ILEMX=7,INORT=14,INSTEP=15,
cbd     INFLD=16,IIDERV=17,IItypeSite=18,IIFLGGIB=19,
cdb     IIAFLG=20,IIRFLG=21, !reserved, are not used!
cdb     INDIM=22,IIPT=23,IJKMN=24,IJMMN=25,IITM=26,
cdb     IITD=27,IITW=28,IIPSI0=29, Ilband=30,Ikband=31,
cdb     Ildelta=32,Ikdelta=33,Ilptmx=34,Ikptmx=35,
cdb     IKEYLIB=36,iNAddN=37,iin2Add=38)      
c#########################################################################       
c  
 1    ICON=0 
      FileErr=.FALSE.
      ierr=0
      ISTRING=0
      OPEN ( UNIT=1,FILE=NAME,IOSTAT=ICON,status='old')
      if (ICON.NE.0) then
         if (mpi_rank.eq.0) OPEN ( UNIT=2,FILE='history',
     #                              POSITION='APPEND')
         if (mpi_rank.eq.0) write(2,*)
     #      ' File - ',NAME,' is missing or cannot be opened'
         if (mpi_rank.eq.0) write(*,*)
     #      ' File - ',NAME,' is missing or cannot be opened'
         if (mpi_rank.eq.0) CLOSE ( UNIT=2)
         goto 37
      endif
             ISTRING=ISTRING+1
       if (jfile.eq.1) then            
     	   READ(1,*,ERR=38,IOSTAT=ICON) ntemp,keyiparm(jfile), ntemp1
             Kcontrol=ntemp
             keyCov=ntemp1     	   
     	 else  
     	   READ(1,*,ERR=38,IOSTAT=ICON) ntemp,keyiparm(jfile)
	 end if !jfile.eq.1
     	   
   38    if (ICON.NE.0) then
               if (mpi_rank.eq.0) OPEN ( UNIT=2,FILE='history',
     #                                    POSITION='APPEND')
          if (mpi_rank.eq.0) write(2,*)
     #       ' File - ',NAME,' is missing or wrong in first line'
          if (mpi_rank.eq.0) write(*,*)
     #       ' File - ',NAME,' is missing or wrong in first line'
               if (mpi_rank.eq.0) CLOSE ( UNIT=2)
  	         CLOSE (UNIT=1)
               goto 37
         endif 
c	             	       
         ISTRING=ISTRING+1
         READ(1,*,ERR=36,IOSTAT=ICON) filenprm(jfile),nnsite(jfile)
                 if ((jfile.eq.1).and.(nnsite(1).gt.MXSITE)) then
                  if (mpi_rank.eq.0) OPEN ( UNIT=2,FILE='history',
     #                                      POSITION='APPEND')
                  if (mpi_rank.eq.0) WRITE (2,*)
     #                   ' Number of sites ',nnsite(1),
     #                   'exceeds MXSITE= ', MXSITE
                  if (mpi_rank.eq.0) WRITE (*,*)
     #                   ' Number of sites ',nnsite(1),
     #                   'exceeds MXSITE= ', MXSITE
                  if (mpi_rank.eq.0) CLOSE ( UNIT=2)
                  goto 37
                endif          	       
         nprm=nprm+filenprm(jfile)
         ISTRING=ISTRING+1
     	   READ(1,*,ERR=36,IOSTAT=ICON) iparms(INFLD,1,jfile),ntemp
                if (jfile.eq.1) Naver=ntemp
                if (iparms(INFLD,1,jfile).gt.MXSPT) then
                  if (mpi_rank.eq.0) OPEN ( UNIT=2,FILE='history',
     #                                      POSITION='APPEND')
                  if (mpi_rank.eq.0) WRITE (2,*)
     #                   ' Number of points ',iparms(INFLD,1,jfile),
     #                   'is diminished to ', MXSPT
                  if (mpi_rank.eq.0) WRITE (*,*)
     #                   ' Number of points ',iparms(INFLD,1,jfile),
     #                   'is diminished to ', MXSPT
                  iparms(INFLD,1,jfile)= MXSPT     
                  if (mpi_rank.eq.0) CLOSE ( UNIT=2)
                endif                                                
         ISTRING=ISTRING+1
         READ(1,*,ERR=36,IOSTAT=ICON) keySp(jfile,1),keySp(jfile,2)
         ISTRING=ISTRING+1
         READ(1,*,ERR=36,IOSTAT=ICON) ShiftExp(jfile),rmsn(jfile),ntemp
         if (rmsn(jfile).le.RNDOFF) rmsn(jfile)=1.0d0
             if (jfile.eq.1) iwflag=ntemp
         ISTRING=ISTRING+1
         READ(1,*,ERR=36,IOSTAT=ICON) fparms(IB0,1,jfile),fldi,fldf
         fparms(IRANGE,1,jfile)=fldf-fldi
         if (iparms(INFLD,1,jfile).gt.1) then
           fparms(IDFLD,1,jfile)=fparms(IRANGE,1,jfile)/
     #                         (iparms(INFLD,1,jfile)-1)
         else
           fparms(IDFLD,1,jfile)=ZERO
         endif
         fparms(IFLDI,1,jfile)=fparms(IB0,1,jfile)+fldi
         iparms(IIDERV,1,jfile)=1	!IDERIV  (0th/1st derivative flag for spectrum)
	   ISTRING=ISTRING+1 
         READ(1,*,ERR=36,IOSTAT=ICON) temp1, temp2, temp3
             if (jfile.eq.1) then
               xtol=temp1
               ftol=temp2
               gtol=temp3
             endif
c  loop over sites
      do	js=1,nnsite(jfile)
         if (js.eq.1) then
           iparms(IItypeSite,js,jfile)=0         
         else
           ISTRING=ISTRING+1
     	     READ(1,*,ERR=36,IOSTAT=ICON) iparms(IItypeSite,js,jfile)
c  the parameters common for all sites:
           iparms(INFLD,js,jfile)=iparms(INFLD,1,jfile)
           fparms(IB0,js,jfile)=fparms(IB0,1,jfile)
           fparms(IFLDI,js,jfile)= fparms(IFLDI,1,jfile)
           fparms(IRANGE,js,jfile)= fparms(IRANGE,1,jfile)
           fparms(IDFLD,js,jfile)=fparms(IDFLD,1,jfile)
           iparms(IIDERV,js,jfile)=1
         endif 
         ISTRING=ISTRING+1
     	   READ(1,*,ERR=36,IOSTAT=ICON) 
     &      iparms(IIN2,js,jfile),fparms(IGAMAN,js,jfile)
         ISTRING=ISTRING+1
    	   READ(1,*,ERR=36,IOSTAT=ICON) 
     &      iparms(iNaddN,js,jfile),
     &      (iparms(iin2add-1+i,js,jfile),i=1,iparms(iNaddN,js,jfile))  
         ISTRING=ISTRING+1
         READ(1,*,ERR=36,IOSTAT=ICON) 
     &      iparms(INORT,js,jfile),iparms(IKEYLIB,js,jfile)
           if (iparms(INORT,js,jfile).eq.0) iparms(INORT,js,jfile)=1
                if (iparms(INORT,js,jfile).gt.MAXNORT) then
                  if (mpi_rank.eq.0) OPEN ( UNIT=2,FILE='history',
     #                                      POSITION='APPEND')
                  if (mpi_rank.eq.0) WRITE (2,*)
     #                   ' Indicated  NORT ',iparms(INORT,js,jfile),
     #                   'is diminished to ', MAXNORT
                  if (mpi_rank.eq.0) WRITE (*,*)
     #                   ' Indicated  NORT ',iparms(INORT,js,jfile),
     #                   'is diminished to ', MAXNORT                
                  iparms(INORT,js,jfile)= MAXNORT       
                  if (mpi_rank.eq.0) CLOSE ( UNIT=2)
                endif          	                 
         ISTRING=ISTRING+1
     	   READ(1,*,ERR=36,IOSTAT=ICON) 
     &      iparms(IIPDF,js,jfile),iparms(IIPT,js,jfile)
            if (iparms(IIPDF,js,jfile).gt.2) iparms(IIPDF,js,jfile)=2
            if (iparms(IIPDF,js,jfile).lt.0) iparms(IIPDF,js,jfile)=0

         ISTRING=ISTRING+1
     	   READ(1,*,ERR=36,IOSTAT=ICON)
     #      iparms(IITW,js,jfile),iparms(IIST,js,jfile)
            if (iparms(IIST,js,jfile).lt.0) iparms(IIST,js,jfile)=0
         ISTRING=ISTRING+1       
         READ(1,*,ERR=36,IOSTAT=ICON) 
     &      fparms(ICGTOL,js,jfile),iparms(INSTEP,js,jfile)
cAXV           fparms(ISHIFTR,js,jfile)=fparms(ICGTOL,js,jfile)
             fparms(ISHIFTR,js,jfile)=1.0D-3
         ISTRING=ISTRING+1
         READ(1,*,ERR=36,IOSTAT=ICON)iparms(ILEMX,js,jfile),
     &      iparms(ILEMX+1,js,jfile),iparms(ILEMX+3,js,jfile),
     &      iparms(ILEMX+5,js,jfile),iparms(ILEMX+6,js,jfile) 
cAXV     ILEMX+1  - lomx
cAXV     ILEMX+2  - kmn
cAXV     ILEMX+3  - kmx
cAXV     ILEMX+4  - mmn 
cAXV     ILEMX+5  - mmx
cAXV     ILEMX+6  - ipnmx
         ISTRING=ISTRING+1
         READ(1,*,ERR=36,IOSTAT=ICON) 
     &     iparms(IML,js,jfile),iparms(IML+1,js,jfile),
     &     iparms(IML+2,js,jfile)
         ISTRING=ISTRING+1
         Read(1,*) str1    !word 'flag'
c varied parameters     
         ISTRING=ISTRING+1
     	   READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(ISPSHFT,js,jfile),iflvar(ISPSHFT,js,jfile),
     &     fparms(IPHASE,js,jfile),iflvar(IPHASE,js,jfile), 	   
     &     fparms(IAMPFAC,js,jfile),iflvar(IAMPFAC,js,jfile)	  
         ISTRING=ISTRING+1
     	   READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(IGIB0,js,jfile),iflvar(IGIB0,js,jfile),
     &     fparms(IGIB2,js,jfile),iflvar(IGIB2,js,jfile), 	   
     &     fparms(IGIBX,js,jfile),iflvar(IGIBX,js,jfile),
     &     fparms(ILBG,js,jfile),iflvar(ILBG,js,jfile)
     	   fparms(IGIB0,js,jfile)=dabs(fparms(IGIB0,js,jfile))
         fparms(IGIBX,js,jfile)=dabs(fparms(IGIBX,js,jfile)) 
     	   if (fparms(IGIBX,js,jfile) .gt. 0.d0) then
     	      iparms(IIFLGGIB,js,jfile)=1
     	      fparms(IGIB2,js,jfile)=0.d0
     	   else
     	      iparms(IIFLGGIB,js,jfile)=0
     	   endif      
         ISTRING=ISTRING+1
     	   READ(1,*,ERR=36,IOSTAT=ICON)   
     &     fparms(IWXX,js,jfile),iflvar(IWXX,js,jfile),
     &     fparms(IWYY,js,jfile),iflvar(IWYY,js,jfile), 	   
     &     fparms(IWZZ,js,jfile),iflvar(IWZZ,js,jfile),
     &     fparms(ILBL,js,jfile),iflvar(ILBL,js,jfile) 
c    
     	  fparms(IWXX,js,jfile)=dabs(fparms(IWXX,js,jfile)) 
     	  fparms(IWYY,js,jfile)=dabs(fparms(IWYY,js,jfile)) 
     	  fparms(IWZZ,js,jfile)=dabs(fparms(IWZZ,js,jfile)) 
        axialw = (((fparms(IWXX,js,jfile).eq.fparms(IWYY,js,jfile)).and.
     #         ((iflvar(IWXX,js,jfile)+iflvar(IWYY,js,jfile)).eq.0)).or.
     #  ((iflvar(IWYY,js,jfile).eq.4).or.(iflvar(IWYY,js,jfile).eq.5))) 
c
        ISTRING=ISTRING+1
     	  READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(IGXX,js,jfile),iflvar(IGXX,js,jfile),
     &     fparms(IGYY,js,jfile),iflvar(IGYY,js,jfile), 
     &     fparms(IGZZ,js,jfile),iflvar(IGZZ,js,jfile) 
        axialg = (((fparms(IGXX,js,jfile).eq.fparms(IGYY,js,jfile)).and.
     #         ((iflvar(IGXX,js,jfile)+iflvar(IGYY,js,jfile)).eq.0)).or.
     #  ((iflvar(IGYY,js,jfile).eq.4).or.(iflvar(IGYY,js,jfile).eq.5))) 
c
        ISTRING=ISTRING+1
     	  READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(IAXX,js,jfile),iflvar(IAXX,js,jfile),
     &     fparms(IAYY,js,jfile),iflvar(IAYY,js,jfile), 	   
     &     fparms(IAZZ,js,jfile),iflvar(IAZZ,js,jfile) 
        axiala = (((fparms(IAXX,js,jfile).eq.fparms(IAYY,js,jfile)).and.
     #         ((iflvar(IAXX,js,jfile)+iflvar(IAYY,js,jfile)).eq.0)).or.
     #  ((iflvar(IAYY,js,jfile).eq.4).or.(iflvar(IAYY,js,jfile).eq.5))) 
c
        ISTRING=ISTRING+1
     	  READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(IALM,js,jfile),iflvar(IALM,js,jfile),
     &     fparms(IBEM,js,jfile),iflvar(IBEM,js,jfile), 	   
     &     fparms(IGAM,js,jfile),iflvar(IGAM,js,jfile) 
        ISTRING=ISTRING+1
     	  READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(IDX,js,jfile),iflvar(IDX,js,jfile),
     &     fparms(IDY,js,jfile),iflvar(IDY,js,jfile), 	   
     &     fparms(IDZ,js,jfile),iflvar(IDZ,js,jfile)
        axialr = (((fparms(IDX,js,jfile).eq.fparms(IDY,js,jfile)).and.
     #         ((iflvar(IDX,js,jfile)+iflvar(IDY,js,jfile)).eq.0)).or.
     #  ((iflvar(IDY,js,jfile).eq.4).or.(iflvar(IDY,js,jfile).eq.5))) 
      	       	  
        ISTRING=ISTRING+1
     	  READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(IALD,js,jfile),iflvar(IALD,js,jfile),
     &     fparms(IBED,js,jfile),iflvar(IBED,js,jfile), 	   
     &     fparms(IGAD,js,jfile),iflvar(IGAD,js,jfile) 	       	  
        ISTRING=ISTRING+1
     	  READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(IPML,js,jfile),iflvar(IPML,js,jfile),
     &     fparms(IPMXY,js,jfile),iflvar(IPMXY,js,jfile), 	   
     &     fparms(IPMZZ,js,jfile),iflvar(IPMZZ,js,jfile) 	  
        ISTRING=ISTRING+1
     	  READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(IDJF,js,jfile),iflvar(IDJF,js,jfile),
     &     fparms(IDJFPRP,js,jfile),iflvar(IDJFPRP,js,jfile), 	   
     &     fparms(IOSS,js,jfile),iflvar(IOSS,js,jfile) 	  
        ISTRING=ISTRING+1
     	  READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(IC20,js,jfile),iflvar(IC20,js,jfile),
     &     fparms(IC20+1,js,jfile),iflvar(IC20+1,js,jfile), 	   				
     &     fparms(IC20+2,js,jfile),iflvar(IC20+2,js,jfile)	  
        ISTRING=ISTRING+1
    	  READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(IC20+3,js,jfile),iflvar(IC20+3,js,jfile),
     &     fparms(IC44,js,jfile),iflvar(IC44,js,jfile), 	   
     &     fparms(IDC20,js,jfile),iflvar(IDC20,js,jfile)
        iprom=0
        do i=0,5
          if ((dabs(fparms(IC20+i,js,jfile)).ne.0.0d0).or.
     #        (iflvar(IC20+i,js,jfile).ne.0)) iprom=i+1 
        enddo
        if (iprom.ne.iparms(IIPT,js,jfile)) then
           if (mpi_rank.eq.0) write(*,*)' Number of potential members',
     #       iparms(IIPT,js,jfile), 'is changed to', iprom
           iparms(IIPT,js,jfile)=iprom
        endif            
        ISTRING=ISTRING+1
    	  READ(1,*,ERR=36,IOSTAT=ICON)
     &     fparms(IANGLIB,js,jfile),iflvar(IANGLIB,js,jfile),
     &     fparms(IANGLIB+1,js,jfile),iflvar(IANGLIB+1,js,jfile),
     &     fparms(IANGLIB+2,js,jfile),iflvar(IANGLIB+2,js,jfile) 
        ISTRING=ISTRING+1
    	  READ(1,*,ERR=36,IOSTAT=ICON)
     &     fparms(IANGLIB+3,js,jfile),iflvar(IANGLIB+3,js,jfile),
     &     fparms(IANGLIB+4,js,jfile),iflvar(IANGLIB+4,js,jfile),
     &     fparms(IANGLIB+5,js,jfile),iflvar(IANGLIB+5,js,jfile) 
        ISTRING=ISTRING+1
    	  READ(1,*,ERR=36,IOSTAT=ICON) 
     &     fparms(IANGFI,js,jfile),iflvar(IANGFI,js,jfile),
     &     fparms(IPSI,js,jfile),iflvar(IPSI,js,jfile),
     &     fparms(IANGKSI,js,jfile),iflvar(IANGKSI,js,jfile) 
c
        do i=0,iparms(iNaddN,js,jfile)-1
          ISTRING=ISTRING+1
    	    READ(1,*,ERR=36,IOSTAT=ICON) 
     &      fparms(IAadd+2*i,js,jfile),iflvar(IAadd+2*i,js,jfile),
     &      fparms(IAadd+2*i+1,js,jfile),iflvar(IAadd+2*i+1,js,jfile),
     &      IEQs(i+1,js,jfile) 
        enddo
        do i=0,iparms(iNaddN,js,jfile)-1
           if (iparms(iin2add+i,js,jfile).gt.2)
     #  write(*,*)' Only additional nuclea with s=1/2 and 1 are carried'
           if (IEQs(i+1,js,jfile).gt.0)then
              iparms(iin2add+i,js,jfile)=
     #            iparms(iin2add+IEQs(i+1,js,jfile)-1,js,jfile)
              fparms(IAadd+2*i,js,jfile)=
     #            fparms(IAadd+2*(IEQs(i+1,js,jfile)-1),js,jfile)
              iflvar(IAadd+2*i,js,jfile)=0
              fparms(IAadd+2*i+1,js,jfile)=
     #            fparms(IAadd+2*IEQs(i+1,js,jfile)-1,js,jfile)
              iflvar(IAadd+2*i+1,js,jfile)=0
           endif
        enddo    
      enddo !js
c----------------------------------------------------------------------
c     Check for zero field
c----------------------------------------------------------------------
      if (fparms(IB0,1,jfile) .lt. ZERO) then
         ierr=POSB0
      elseif (fparms(IB0,1,jfile) .eq. ZERO) then
         ierr=ZEROB0
      endif
	if (ierr.gt.0)write (*,*) eprerr(ierr)
c----------------------------------------------------------------------
c     Check for nonzero number of points and field range
c----------------------------------------------------------------------
      if (iparms(INFLD,1,jfile).le.0) then
         ierr=ZERONFLD
           OPEN ( UNIT=2,FILE='history', POSITION='APPEND')
           write (*,*) eprerr(ierr)
           write (2,*) eprerr(ierr)
           CLOSE ( UNIT=2)           
           goto 37         
      end if

  36   if (ICON.NE.0) then
    	 write(*,*)'Error in line No',ISTRING,' in file - ',NAME
       goto 37       
	 end if
  	         CLOSE (UNIT=1)
c
c   the parameters for sites: 
      if ((jfile.gt.1).and.(keyiparm(jfile).eq.0))
     #       nnsite(jfile)=nnsite(1)
 
      do js=2,nnsite(jfile)
c  the same parameters for all sites:      
           fparms(IANGFI,js,jfile)=fparms(IANGFI,1,jfile)
           fparms(IPSI,js,jfile)=fparms(IPSI,1,jfile)
           fparms(IANGKSI,js,jfile)=fparms(IANGKSI,1,jfile)
           fparms(ISPSHFT,js,jfile)=fparms(ISPSHFT,1,jfile)
c  the case of sites which differs by additional HFI only:
           if (iparms(IItypeSite,js,jfile).gt.0) then
              nInSite=iparms(IItypeSite,js,jfile)          
              do jpar=1,IIDERV   !1,17
                iparms(jpar,js,jfile)=iparms(jpar,nInSite,jfile)
              enddo  !jpar
              do jpar=IIFLGGIB,IKEYLIB   !19,36
                iparms(jpar,js,jfile)=iparms(jpar,nInSite,jfile)
              enddo  !jpar 
              fparms(1,js,jfile)=fparms(1,nInSite,jfile)
              do jpar=IWXX,IDFLD     !5,45
                 fparms(jpar,js,jfile)=fparms(jpar,nInSite,jfile)
                 if ((iflvar(jpar,js,jfile).eq.1).or.
     #               (iflvar(jpar,js,jfile).eq.2).or.
     #               (iflvar(jpar,js,jfile).eq.11).or.
     #               (iflvar(jpar,js,jfile).eq.12)) then
                     nprm=nprm-1         
                     filenprm(jfile)=filenprm(jfile)-1
                 endif                
                 iflvar(jpar,js,jfile)=0                 
              enddo  !jpar  
              do jpar=ISPSHFT,Ireserved     !47,59 
                 fparms(jpar,js,jfile)=fparms(jpar,nInSite,jfile)
                 if ((iflvar(jpar,js,jfile).eq.1).or.
     #               (iflvar(jpar,js,jfile).eq.2).or.
     #               (iflvar(jpar,js,jfile).eq.11).or.
     #               (iflvar(jpar,js,jfile).eq.12)) then
                     nprm=nprm-1         
                     filenprm(jfile)=filenprm(jfile)-1
                 endif                
                 iflvar(jpar,js,jfile)=0                               
              enddo  !jpar                            
           endif !iparms(IItypeSite,js,jfile)
      enddo  !js
      
c
c Some check of the parameters
c AXV The checking had moved here from file lcheck.f   
      if ((jfile.eq.1).or.(keyiparm(jfile).ne.0)) then
       do js=1,nnsite(jfile)
c----------------------------------------------------------------------
c     Check for zero values in g-tensor and zero A tensor
c----------------------------------------------------------------------
         if ( dabs(fparms(IGXX,js,jfile)) .lt. RNDOFF  .or.
     #     dabs(fparms(IGYY,js,jfile)) .lt. RNDOFF  .or.
     #     dabs(fparms(IGZZ,js,jfile)) .lt. RNDOFF)   then
           ierr=ZEROG
           OPEN ( UNIT=2,FILE='history', POSITION='APPEND')
           write (*,*) eprerr(ierr)
           write (2,*) eprerr(ierr)
           CLOSE ( UNIT=2)           
           goto 37
         endif
c
         if ( (iparms(IIN2,js,jfile) .le. 0) .or.
     #    (dabs(fparms(IAXX,js,jfile)) .lt. RNDOFF .and. 
     #     dabs(fparms(IAYY,js,jfile)) .lt. RNDOFF .and.
     #     dabs(fparms(IAZZ,js,jfile)) .lt. RNDOFF)     ) then
           iparms(IIN2,js,jfile)=0
           fparms(IAXX,js,jfile)=ZERO
           fparms(IAYY,js,jfile)=ZERO
           fparms(IAZZ,js,jfile)=ZERO
           ierr=ZEROIN2
           write (*,*) eprerr(ierr)
         endif
c----------------------------------------------------------------------
c Issue warning if high-field approximation has been violated
c (this is a very rough criterion for the high-field approx.!)
c----------------------------------------------------------------------
         g0=(fparms(IGXX,js,jfile)+fparms(IGYY,js,jfile)+
     #       fparms(IGZZ,js,jfile))/3.0D0
         gmin=dmin1(fparms(IGXX,js,jfile),fparms(IGYY,js,jfile),
     #           fparms(IGZZ,js,jfile) )
         gmax=dmax1( fparms(IGXX,js,jfile),fparms(IGYY,js,jfile),
     #           fparms(IGZZ,js,jfile) )
         hmax=dmax1(dabs(fparms(IAXX,js,jfile)),
     #           dabs(fparms(IAYY,js,jfile)),
     #           dabs(fparms(IAXX,js,jfile)))
         if (fparms(IB0,1,jfile) .lt.
     # 10.0D0*dmax1( hmax, (gmax-gmin)*fparms(IB0,1,jfile)/g0)) then
           ierr=BADHIFLD
           write (*,*) eprerr(ierr)
         endif
         if ( dabs(fparms(IALM,js,jfile)).gt.RNDOFF .or.
     #        dabs(fparms(IGAM,js,jfile)).gt.RNDOFF) then
	         iparms(IITM,js,jfile)=2
         else if	(dabs(fparms(IBEM,js,jfile)).gt.RNDOFF) then 
               iparms(IITM,js,jfile)=1
              else 
               iparms(IITM,js,jfile)=0
         endif
         if ( dabs(fparms(IALD,js,jfile)).gt.RNDOFF .or.
     #        dabs(fparms(IGAD,js,jfile)).gt.RNDOFF) then
               iparms(IITD,js,jfile)=2
         else if	(fparms(IBED,js,jfile).gt.RNDOFF) then 
               iparms(IITD,js,jfile)=1
              else 
               iparms(IITD,js,jfile)=0
         endif
c  Set exponents according to non-Brownian model
         fparms(IEXPL,js,jfile)=ONE
         fparms(IEXPKXY,js,jfile)=ONE
         fparms(IEXPKZZ,js,jfile)=ONE
         if (iparms(IML,js,jfile).eq.1)  fparms(IEXPL,js,jfile)=HALF
         if (iparms(IML+1,js,jfile).eq.1) fparms(IEXPKXY,js,jfile)=HALF
         if (iparms(IML+2,js,jfile).eq.1) fparms(IEXPKZZ,js,jfile)=HALF
c    Set residence times to zero for Brownian model
        if (iparms(IIPDF,js,jfile) .eq. 0) then
           fparms(IPML,js,jfile) = ZERO
           fparms(IPMXY,js,jfile) = ZERO
           fparms(IPMZZ,js,jfile) = ZERO
        else
c     *** Non-Brownian motion: must have at least one nonzero residence
c         time, and no potential specified
          if (iparms(IIPT,js,jfile) .gt. 0) then
            ierr=BADJMP
            iparms(IIPDF,js,jfile)=0
            write (*,*) eprerr(ierr)
          endif
        endif 
c     *** Anisotropic viscosity model: must have potential and
c         no discrete jump motion specified
        if (iparms(IIPDF,js,jfile) .eq. 2) then
         if (iparms(IIPT,js,jfile) .eq. 0) then
            ierr=BADAV
            iparms(IIPDF,js,jfile)=0
            write (*,*) eprerr(ierr)
         endif
c
         if (iparms(IIPT,js,jfile) .gt. 0) then
            ierr=BADDJ
            iparms(IIPDF,js,jfile)=0
            write (*,*) eprerr(ierr)
         endif
        endif
         if((iparms(ILEMX,js,jfile).gt.MXLVAL).and.
     #         (iparms(IIPT,js,jfile).ne.0)) then
            ierr=LEMXHI
            iparms(ILEMX,js,jfile)=MXLVAL
            write (*,*) eprerr(ierr)
         endif
c----------------------------------------------------------------------
c     Check specified tolerances and shifts for CG calculation
c----------------------------------------------------------------------
        if (dabs(fparms(ICGTOL,js,jfile)) .lt. TEN*RNDOFF) then
         fparms(ICGTOL,js,jfile)=TEN*RNDOFF
         ierr=CGTOLHI
         write (*,*) eprerr(ierr)
        end if
         if (dabs(fparms(ISHIFTR,js,jfile)) .lt. TEN*RNDOFF) then
         fparms(ISHIFTR,js,jfile)=10.0D0*RNDOFF
         ierr=SHIFTHI
         write (*,*) eprerr(ierr)
        end if 
        SELECT CASE (iparms(IIPT,js,jfile)) 
        CASE (0,1)
           iparms(Ilptmx,js,jfile)=2
        CASE (2)
           iparms(Ilptmx,js,jfile)=2
           iparms(Ikptmx,js,jfile)=2
        CASE (3)
           iparms(Ilptmx,js,jfile)=4
        CASE (4)
           iparms(Ilptmx,js,jfile)=4
           iparms(Ikptmx,js,jfile)=2
        CASE (5)
           iparms(Ilptmx,js,jfile)=4
           iparms(Ikptmx,js,jfile)=4
        CASE DEFAULT
            WRITE(*,*) '  Improper number of potential members!'
        END SELECT     
        if (iparms(Ilptmx,js,jfile).ge.2) then
           iparms(Ilband,js,jfile)=iparms(Ilptmx,js,jfile)*2
        else
           iparms(Ilband,js,jfile)=2
        end if
c
        iparms(Ikband,js,jfile)=iparms(Ikptmx,js,jfile)*2
cAXV nonaxial:
cAXV       if (axialr) 
        if (.not.axialr) 
     #    iparms(Ikband,js,jfile)=iparms(Ikband,js,jfile)+2
c----------------------------------------------------------------------
c  Determine basis set using rules in M,I,I,M & F
c  (1)   If there is no diffusion or magnetic tilt, only even K values
c        are needed
c----------------------------------------------------------------------
        if(iparms(IITM,js,jfile).eq.0 .and.
     #     iparms(IITD,js,jfile).eq.0) then
           iparms(Ikdelta,js,jfile)=2
        else
           iparms(Ikdelta,js,jfile)=1
        end if
c----------------------------------------------------------------------
c Use only even L values (ldelta=2) and no K values (kmx=0) in the case
c   (1) there is no magnetic or diffusion tilt (kdelta.eq.2)
c   (2) magnetic tensors are axial
c   (3) potential is axial (no K terms in potential)
c----------------------------------------------------------------------
        if (axiala .and. axialg .and. axialw .and.
     #      (iparms(Ikdelta,js,jfile).eq.2).and.  
     #      (iparms(Ikptmx,js,jfile).eq.0)) then
           iparms(Ildelta,js,jfile)=2
           iparms(ILEMX+3,js,jfile)=0
        else
           iparms(Ildelta,js,jfile)=1
        end if 
c    
         iparms(IIPSI0,js,jfile)=0
         if( ((dabs(fparms(IPSI,1,jfile)).gt.RNDOFF).and.
     #      (dabs(fparms(IPSI,1,jfile)-180.0d0).gt.RNDOFF))
     #      .or. iparms(INORT,js,jfile).gt.1) iparms(IIPSI0,js,jfile)=1
c    --- Allow antisymmetric K combinations if any tensor has
c        imaginary elements (nonzero alm, gam, ald, or gad)
         iparms(IJKMN,js,jfile)=1
         if (iparms(IITM,js,jfile).eq.2.or.
     #       iparms(IITD,js,jfile).eq.2) iparms(IJKMN,js,jfile)=-1
c    --- Allow antisymmetric M combinations if there is a nonzero
c        nuclear Zeeman interaction
         if (dabs(fparms(IGAMAN,js,jfile)).ge.RNDOFF .and.
     #       iparms(IIN2,js,jfile).gt.0) then
            iparms(IJMMN,js,jfile)=-1
         else
            iparms(IJMMN,js,jfile)=1
         end if
c----------------------------------------------------------------------
c     Check basis set parameters
c----------------------------------------------------------------------
         if(ipar(iparms(ILEMX,js,jfile)).ne.1)
     #       iparms(ILEMX,js,jfile)=iparms(ILEMX,js,jfile)-1
         if(ipar(iparms(ILEMX+1,js,jfile)).ne.-1)
     #       iparms(ILEMX+1,js,jfile)=iparms(ILEMX+1,js,jfile)-1
         if(iparms(ILEMX+1,js,jfile).gt.iparms(ILEMX,js,jfile))
     #        iparms(ILEMX+1,js,jfile)=iparms(ILEMX,js,jfile)-1
         if(iparms(ILEMX+3,js,jfile).gt.iparms(ILEMX,js,jfile))
     #        iparms(ILEMX+3,js,jfile)=iparms(ILEMX,js,jfile)
         if(iparms(ILEMX+5,js,jfile).gt.iparms(ILEMX,js,jfile))
     #        iparms(ILEMX+5,js,jfile)=iparms(ILEMX,js,jfile)
         if(ipar(iparms(ILEMX+3,js,jfile)).ne.1)
     #        iparms(ILEMX+3,js,jfile)=iparms(ILEMX+3,js,jfile)-1
         if(iparms(ILEMX+6,js,jfile).gt.iparms(IIN2,js,jfile))
     #        iparms(ILEMX+6,js,jfile)=iparms(IIN2,js,jfile)

         if(-iparms(ILEMX+4,js,jfile).gt.iparms(ILEMX+5,js,jfile))
     #        iparms(ILEMX+4,js,jfile)=-iparms(ILEMX+5,js,jfile)
         if(iparms(IJKMN,js,jfile).eq.1) 
     #        iparms(ILEMX+2,js,jfile)=max(iparms(ILEMX+2,js,jfile),0)
c     --- Get truncation value for antisymmetric K-space if needed
         if (iparms(IJKMN,js,jfile).eq.-1) 
     #        iparms(ILEMX+2,js,jfile)=-iparms(ILEMX+3,js,jfile)
cAXV It's redundant:
cAXV         if(-iparms(ILEMX+2,js,jfile).gt.iparms(ILEMX+3,js,jfile))
cAXV     #        iparms(ILEMX+2,js,jfile)=-iparms(ILEMX+3,js,jfile)

         if(iparms(IJMMN,js,jfile).eq.1) 
     #        iparms(ILEMX+4,js,jfile)=max(iparms(ILEMX+2,js,jfile),0)
c     Get truncation value for antisymmetric M-space if needed
         if (iparms(IJMMN,js,jfile).eq.-1)
     #        iparms(ILEMX+4,js,jfile)=-iparms(ILEMX+5,js,jfile)
c----------------------------------------------------------------------
c     check number of Lanczos/CG steps to execute 
c----------------------------------------------------------------------
         if (iparms(INSTEP,js,jfile).gt.MXSTEP) then
           ierr=NSTEPHI
           write (*,*) eprerr(ierr)
           iparms(INSTEP,js,jfile)=MXSTEP
         endif
c
c flags for g-values:
         lpr=(((iflvar(IGXX,js,jfile).ne.0).or.
     #            (iflvar(IGYY,js,jfile).ne.0))
     #  .and.((iflvar(IGZZ,js,jfile).ne.0).or.
     #            (iflvar(IAXX,js,jfile).ne.0))
     #  .and.((iflvar(IAYY,js,jfile).ne.0).or.
     #            (iflvar(IAZZ,js,jfile).ne.0)))
          SELECT CASE (iparms(IKEYLIB,js,jfile))
          CASE (5:)
             write(*,*)'Unknown libration model'
             goto 37
ctemp          CASE (1:4)
ctemp             iprom=0
ctemp             do i=0,5
ctemp               iprom=iprom+iflvar(IANGLIB+i,js,jfile)
ctemp             enddo
ctemp             if ((iprom.gt.0).and.lpr) then
ctemp               write(*,*)'  Libration angles and magnetic parameter'
ctemp               write(*,*)'  should not be varied simultaneously!!!'
ctemp             endif
          END SELECT
        lpr=(((iflvar(IGXX,js,jfile).eq.1).or.
     #            (iflvar(IGXX,js,jfile).eq.2))
     #  .and.((iflvar(IGYY,js,jfile).eq.1).or.
     #            (iflvar(IGYY,js,jfile).eq.2))
     #  .and.((iflvar(IGZZ,js,jfile).eq.1).or.
     #            (iflvar(IGZZ,js,jfile).eq.2)))
          if ((iflvar(ISPSHFT,js,jfile).GT.0).and.lpr) then
             write(*,*)'Sdvig and all g-values should not'
             write(*,*)'   be varied simultaneously'
             goto 37            
          endif
       enddo !js
      endif !jfile=1 or keyiparm <>0
cc  the parameters for second and other files is taken from the first one.
      if (jfile.gt.1) then
         do js=1,nnsite(jfile)
            do jpar=1,NFPRM 
              SELECT CASE (jpar)
		      CASE (IRANGE:IDFLD,IB0)
                if (keyiparm(jfile).eq.0) 
     #              fparms(jpar,js,jfile)=fparms(jpar,js,1)
              CASE (IANGFI:IANGFI)
                    continue 
              CASE (IPSI)
                    continue                         
              CASE default
                if (iflvar(jpar,js,jfile).lt.10)then	 
                  fparms(jpar,js,jfile)=fparms(jpar,js,1)
                  iflvar(jpar,js,jfile)=0
                endif
              END SELECT
            enddo  !jpar
            if (keyiparm(jfile).eq.0) then
              do jpar=1,NIPRM
                iparms(jpar,js,jfile)=iparms(jpar,js,1)
              enddo  !jpar
              fparms(IRANGE,js,jfile)=fparms(IRANGE,js,1)
              fparms(IDFLD,js,jfile)=fparms(IDFLD,js,1)
              fparms(IFLDI,js,jfile)=fparms(IFLDI,js,1)
              fparms(IEXPL,js,jfile)=fparms(IEXPL,js,1)
              fparms(IEXPKXY,js,jfile)=fparms(IEXPKXY,js,1)
              fparms(IEXPKZZ,js,jfile) =fparms(IEXPKZZ,js,1)
            endif
            iparms(IIPSI0,js,jfile)=0
            if( ((dabs(fparms(IPSI,1,jfile)).gt.RNDOFF).and.
     #        (dabs(fparms(IPSI,1,jfile)-180.0d0).gt.RNDOFF))
     #       .or. iparms(INORT,js,jfile).gt.1) iparms(IIPSI0,js,jfile)=1

         enddo !js
         keySp(jfile,1)=keySp(1,1)
         keySp(jfile,2)=keySp(1,2)
         ShiftExp(jfile)=ShiftExp(1)
         rmsn(jfile)=rmsn(1)
      endif  
      IF (FileErr) then
	    goto 37
	ELSE
          return
	ENDIF
37	FileErr=.true.
      if (mpi_rank.eq.0) OPEN ( UNIT=2,FILE='history',
     #                           POSITION='APPEND')
      if (mpi_rank.eq.0) WRITE (2,*)  'FILE - ',NAME,
     #                                'SHOULD BE IMPROVED'
      if (mpi_rank.eq.0) WRITE (*,*)  'FILE - ',NAME,
     #                                'SHOULD BE IMPROVED'
      if (mpi_rank.eq.0) CLOSE ( UNIT=2)
 	 return
  	 end
