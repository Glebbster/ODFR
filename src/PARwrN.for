C for ODFR 01/12/2014
	SUBROUTINE  PARWrite(AME,jfile)
c
          IMPLICIT DOUBLE PRECISION (A-H),DOUBLE PRECISION (O-Z)
	character*30 ame
c    
      include 'dimnls.inc'
      include 'eprprm.inc'
      include 'parcom.inc'
      include 'expdat.inc'
      include 'lmcom.inc'
      include 'iterat.inc'
      include 'mpi_state.inc'
c
        integer iparms(NIPRM,MXSITE,MXSPC)
        DOUBLE PRECISION fparms(NFPRM,MXSITE,MXSPC)
        COMMON /COMPARMS/ iparms,fparms
        integer IEQs
        COMMON /IEQ/ IEQs(MaxAddNucl,MXSITE,MXSPC)
c
	  integer iflvar, filenprm
        common /iflvar/iflvar(NFPRM,MXSITE,MXSPC),filenprm(MXSPC)
c
        integer Kcontrol
        common /control/Kcontrol
        integer Naver,js,ntemp,keyiparm 
       common /Naver/Naver,keyiparm(MXSPC)
	double precision MagPar0
      common /magpar0/magpar0(6,MXSITE)
       DOUBLE PRECISION  parReserved
       double precision function angfun
       integer keyCov 
       common/KeyCov/keyCov		          
c             
c ###########################################################################
c
      if (mpi_rank.ne.0) return
c   integer parameters common for all sites
             ntemp=0
	         OPEN ( UNIT=1,FILE=AME)
     		   WRITE(1,101,IOSTAT=ICON) Kcontrol,keyiparm(jfile),keyCov
 101  format(i2,4x,i2,4x,i2,'   Kcontrol, keyiparm, keyCov' )
     		   WRITE(1,102,IOSTAT=ICON) filenprm(jfile),nnsite(jfile)
 102  format(i2,4x,i2,'   number of var, number of sites')
     		   WRITE(1,103,IOSTAT=ICON) iparms(INFLD,1,jfile),Naver
 103  format(i4,2x,i2,'   number of points, NumAve')
     		   WRITE(1,104,IOSTAT=ICON) keySp(jfile,1),keySp(jfile,2)
 104  format(i2,4x,i2,'   key for exp.spectr,key for calc.spectr')
     		   WRITE(1,105,IOSTAT=ICON) ShiftExp(jfile),rmsn(jfile),iwflag
 105  format(F7.1,1x,F8.5,i2,
     # '   Shift for exp.spectr, exp.noise, weighted/unweighted')
     		   WRITE(1,106,IOSTAT=ICON) fparms(IB0,1,jfile),
     #                   fparms(IFLDI,1,jfile)-fparms(IB0,1,jfile),
     #                   fparms(IRANGE,1,jfile)+fparms(IFLDI,1,jfile)
     #                   -fparms(IB0,1,jfile)
 106  format(F7.1,2(2x,F5.1),'   B0, fieldi, fieldf')
     		   WRITE(1,107,IOSTAT=ICON) xtol, ftol, gtol
 107  format(G8.1,2(2x,G8.1),'   (XTOL),(FTOL),gradient(GTOL)')	
c  loop over sites
      do	js=1,nnsite(jfile)
         if (js.gt.1) 
     &     WRITE(1,*,IOSTAT=ICON) iparms(IItypeSite,js,jfile)
c     
     		   WRITE(1,108,IOSTAT=ICON)
     &      iparms(IIN2,js,jfile),fparms(IGAMAN,js,jfile)
 108  format(i2,4x,F9.6,'   twice the nuclear spin [in2], gamman')
    	       WRITE(1,110,IOSTAT=ICON)iparms(iNAddN,js,jfile),
     &      (iparms(iin2Add-1+i,js,jfile),i=1,iparms(iNAddN,js,jfile))
 110  format(i2,2x,15(2x,i2),'   number of additinal nuclea')
     		   WRITE(1,109,IOSTAT=ICON)
     &      iparms(INORT,js,jfile),iparms(IKEYLIB,js,jfile)
 109  format(i2,3x,i2,'   NORT, keyLib')
     		   WRITE(1,112,IOSTAT=ICON)
     &      iparms(IIPDF,js,jfile),iparms(IIPT,js,jfile)
 112  format(i2,4x,i2,'   diffusion model,number of potential terms')
      	   WRITE(1,113,IOSTAT=ICON) 
     #      iparms(IITW,js,jfile),iparms(IIST,js,jfile)
 113  format(i2,4x,i2,'   Lorentz frame,number of sites for jumps')
     		   WRITE(1,114,IOSTAT=ICON)
     &      fparms(ICGTOL,js,jfile),iparms(INSTEP,js,jfile)
 114  format(G8.1,4x,i5,
     #'   error tolerance [cgtol],number of Lanczos/CG steps')
       WRITE(1,*,IOSTAT=ICON)iparms(ILEMX,js,jfile),
     &  iparms(ILEMX+1,js,jfile),iparms(ILEMX+3,js,jfile),
     &  iparms(ILEMX+5,js,jfile),iparms(ILEMX+6,js,jfile)
	 WRITE(1,115,IOSTAT=ICON)  
     & iparms(IML,js,jfile),iparms(IML+1,js,jfile),
     & iparms(IML+2,js,jfile)
 115  format(i3,2(2x,i3),'   non-Brownian flags ml,mxy,mzz')  
       WRITE(1,*) '       flag' 	      
c parameters which can be varied:
       parReserved =0.d0
     	  WRITE(1,116,IOSTAT=ICON) 
     & fparms(ISPSHFT,js,jfile),iflvar(ISPSHFT,js,jfile),
     & fparms(IPHASE,js,jfile),iflvar(IPHASE,js,jfile), 	   
     & fparms(IAMPFAC,js,jfile),iflvar(IAMPFAC,js,jfile)
 116  format(3(F16.12,1x,i2,1x),' Shift, phase, content')
     	  WRITE(1,117,IOSTAT=ICON) 
     & fparms(IGIB0,js,jfile),iflvar(IGIB0,js,jfile),
     & fparms(IGIB2,js,jfile),iflvar(IGIB2,js,jfile), 	   
     & fparms(IGIBX,js,jfile),iflvar(IGIBX,js,jfile),
     & fparms(ILBG,js,jfile),iflvar(ILBG,js,jfile) 
 117  format(4(F16.12,1x,i2,1x),' gib0,gib2,gibX')      	  
     	  WRITE(1,118,IOSTAT=ICON)   
     & fparms(IWXX,js,jfile),iflvar(IWXX,js,jfile),
     & fparms(IWYY,js,jfile),iflvar(IWYY,js,jfile), 	   
     & fparms(IWZZ,js,jfile),iflvar(IWZZ,js,jfile),
     & fparms(ILBL,js,jfile),iflvar(ILBL,js,jfile) 
 118  format(4(F16.12,1x,i2,1x),' Lorentz width,rsv')   	  
      if (iparms(IKEYLIB,js,jfile).eq.0) then
     	  WRITE(1,119,IOSTAT=ICON) 
     & fparms(IGXX,js,jfile),iflvar(IGXX,js,jfile),
     & fparms(IGYY,js,jfile),iflvar(IGYY,js,jfile), 
     & fparms(IGZZ,js,jfile),iflvar(IGZZ,js,jfile)
 119  format(3(F16.12,1x,i2,1x),'  g-values')      
     	  WRITE(1,120,IOSTAT=ICON) 
     & fparms(IAXX,js,jfile),iflvar(IAXX,js,jfile),
     & fparms(IAYY,js,jfile),iflvar(IAYY,js,jfile),
     & fparms(IAZZ,js,jfile),iflvar(IAZZ,js,jfile)
 120  format(3(F16.12,1x,i2,1x),'  A-values')       
      else
     	  WRITE(1,119,IOSTAT=ICON) 
     & magpar0(1,js),iflvar(IGXX,js,jfile),
     & magpar0(2,js),iflvar(IGYY,js,jfile), 
     & magpar0(3,js),iflvar(IGZZ,js,jfile)
     	  WRITE(1,120,IOSTAT=ICON) 
     & magpar0(4,js),iflvar(IAXX,js,jfile),
     & magpar0(5,js),iflvar(IAYY,js,jfile), 	   
     & magpar0(6,js),iflvar(IAZZ,js,jfile)
      endif 
     	  WRITE(1,121,IOSTAT=ICON) 
     & angfun(fparms(IALM,js,jfile)),iflvar(IALM,js,jfile),
     & angfun(fparms(IBEM,js,jfile)),iflvar(IBEM,js,jfile), 	   
     & angfun(fparms(IGAM,js,jfile)),iflvar(IGAM,js,jfile)
 121  format(3(F16.12,1x,i2,1x),'  magnetic tilt [alm,bem,gam]')       
     	  WRITE(1,122,IOSTAT=ICON) 
     & fparms(IDX,js,jfile),iflvar(IDX,js,jfile),
     & fparms(IDY,js,jfile),iflvar(IDY,js,jfile), 	   
     & fparms(IDZ,js,jfile),iflvar(IDZ,js,jfile) 
 122  format(3(G16.6,1x,i2,1x),'  diffusion tensor [dx,dy,dz]')     	   
     	  WRITE(1,123,IOSTAT=ICON) 
     & angfun(fparms(IALD,js,jfile)),iflvar(IALD,js,jfile),
     & angfun(fparms(IBED,js,jfile)),iflvar(IBED,js,jfile), 	   
     & angfun(fparms(IGAD,js,jfile)),iflvar(IGAD,js,jfile) 
 123  format(3(F16.12,1x,i2,1x),'  diffusion tilt [ald,bed,gad]')     	 
     	  WRITE(1,124,IOSTAT=ICON) 
     & fparms(IPML,js,jfile),iflvar(IPML,js,jfile),
     & fparms(IPMXY,js,jfile),iflvar(IPMXY,js,jfile), 	   
     & fparms(IPMZZ,js,jfile),iflvar(IPMZZ,js,jfile) 
 124  format(3(F16.12,1x,i2,1x),'  [pml,pmxy,pmzz]')     	  
     	  WRITE(1,125,IOSTAT=ICON) 
     & fparms(IDJF,js,jfile),iflvar(IDJF,js,jfile),
     & fparms(IDJFPRP,js,jfile),iflvar(IDJFPRP,js,jfile), 	   
     & fparms(IOSS,js,jfile),iflvar(IOSS,js,jfile) 
 125  format(3(G16.6,1x,i2,1x),'  jumps [djf]),(djfprp),(oss)')      	  
     	  WRITE(1,126,IOSTAT=ICON) 
     & fparms(IC20,js,jfile),iflvar(IC20,js,jfile),
     & fparms(IC20+1,js,jfile),iflvar(IC20+1,js,jfile), 	   				
     & fparms(IC20+2,js,jfile),iflvar(IC20+2,js,jfile)
 126  format(3(F16.12,1x,i2,1x),'  potentials (2,0),(2,2),(4,0)')      	
    	  WRITE(1,127,IOSTAT=ICON) 
     & fparms(IC20+3,js,jfile),iflvar(IC20+3,js,jfile),
     & fparms(IC44,js,jfile),iflvar(IC44,js,jfile), 	   
     & fparms(IDC20,js,jfile),iflvar(IDC20,js,jfile) 
 127  format(3(F16.12,1x,i2,1x),'  potentials (4,2),(4,4),dlam')      	 
c
        do i=0,5
          fparms(IANGLIB+i,js,jfile)=angfun(fparms(IANGLIB+i,js,jfile))
        enddo  	  
	  WRITE(1,128,IOSTAT=ICON)
     & fparms(IANGLIB,js,jfile),iflvar(IANGLIB,js,jfile),
     & fparms(IANGLIB+1,js,jfile),iflvar(IANGLIB+1,js,jfile),
     & fparms(IANGLIB+2,js,jfile),iflvar(IANGLIB+2,js,jfile) 
    	  WRITE(1,128,IOSTAT=ICON)
     & fparms(IANGLIB+3,js,jfile),iflvar(IANGLIB+3,js,jfile),
     & fparms(IANGLIB+4,js,jfile),iflvar(IANGLIB+4,js,jfile),
     & fparms(IANGLIB+5,js,jfile),iflvar(IANGLIB+5,js,jfile) 
 128  format(3(F16.12,1x,i2,1x),'  librations')
        WRITE(1,129,IOSTAT=ICON) 
     &  fparms(IANGFI,js,jfile),iflvar(IANGFI,js,jfile),
     &  fparms(IPSI,js,jfile),iflvar(IPSI,js,jfile),
     &  fparms(IANGKSI,js,jfile),iflvar(IANGKSI,js,jfile)
 129  format(3(F16.12,1x,i2,1x),'  rsrv,psi,rsrv')      	  
c
         do i=0,iparms(iNaddN,js,jfile)-1
    	    WRITE(1,130,IOSTAT=ICON) 
     & fparms(IAadd+2*i,js,jfile),iflvar(IAadd+2*i,js,jfile),
     & fparms(IAadd+2*i+1,js,jfile),iflvar(IAadd+2*i+1,js,jfile),
     & IEQs(i+1,js,jfile)      
 130  format(2(F16.12,1x,i2,1x),2x,i2,'  AAdd')      	  
	  enddo
      enddo !js
       if (ICON.NE.0) then
    	 write(*,*)'Error in writing the parameter file'
	 end if
  	         CLOSE (UNIT=1)
	return
	end
