c based on eprprm.inc from NLSL Version 1.5.1 beta 2/3/96
c----------------------------------------------------------------------
c                    =========================
c                       common block EPRPRM
c                    =========================
c
c       Variable and common block definition for NLSL fitting of
c  experimental spectra using LBLL-family slow-motional calculations.
c
c
c IMPORTANT NOTES: 
c
c   The order in which the parameters appear in /eprprm/ is critical to 
c   the proper functioning of the NLSL program. If this common block is
c   to be changed, the following rules must be observed:
c
c   (1) To permit alias names to be used for the g,A, and diffusion
c       tensors, the parameters gxx,gyy,gzz,axx,ayy,azz,dx,dy,dz,wxx,wyy,wzz
c       should appear continuously in that order in the common block.
c
c   (4) The residence times tl,tkxy,tkzz have been replaced with the
c       parameters pml, pmxy, pmzz. This is to facilitate fitting using
c       non-Brownian models where it is desirable to vary the residence
c       time and the diffusion rate constant together (i.e. keeping the
c       product, which is related to the rms jump angle, constant).
c
c   (5) The non-Brownian model flags have been changed to ml, mzz, and
c       mxy for perpendicular, parallel, and internal motions, respectively.
c       The flags may be set to 0, 1, and 2 for Brownian, free, and
c       jump diffusion, respectively.
c
c   (7) A new spectral linewidth broadening parameter, "lb", has been
c       added
cAXV:
c      Line broadening is descibed  as Lorentzian broadening lbl 
c      and Gaussian broadening lbg.
c The additional hyperfine splitting is described by Aadd,NAddN,in2Add.
c The libration amlitudes are described by anglib,keylib.
c The orientation of sample in magnetic field is described by ANGFI,PSI,ANGKSI.
c Amplitude of spectrum for every site is  AmpFac.
c Field shift for every spectrum is SpcShift.
c  
c
c***********************************************************************
c
       integer keySp(MXSPC,2)
       common /keySp/ keySp

        double precision gxx,gyy,gzz,axx,ayy,azz,wxx,wyy,wzz,dx,dy,dz,
     #         pml,pmxy,pmzz,djf,djfprp,oss,psi,ald,bed,gad,alm,bem,gam,
     #         c20,c22,c40,c42,c44,gib0,gib2,gibx,lbl,lbg,b0,gamman,
     #         cgtol,shiftr,shifti,range,fldi,dfld,a0,g0,w0,expl,expkxy,
     #         expkzz,phase,dc20,
     #		   AmpFac,SpcShift,ANGFI,ANGKSI,Aadd,anglib,reserved
c
        double precision faa,fgm,fam,fww,fwm,fgd,fad,fwd,cpot,xlk,
     #         g0first, fldShftSite        
c
        integer in2,ipdf,ist,ml,mxy,mzz,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,
     #          nort,nstep,ideriv,ndim,
     #          nfld,jkmn,jmmn,ipt,itm,itd,itw,ipsi0,lband,kband,ldelta,
     #          kdelta,lptmx,kptmx,
     #          NAddN,in2Add,keylib,ItypeSite,IFLGGIB,nresrv1,nresrv2
cAXV	 ,neltot,nelv,nelre,nelim
c
      common /eprprm/ phase,gib0,gib2,gibx,
     #         wxx,wyy,wzz,
     #         gxx,gyy,gzz,
     #         axx,ayy,azz,
     #         dx,dy,dz,
     #         pml,pmxy,pmzz,
     #         djf,djfprp,oss,
     #         psi,ald,bed,gad,alm,bem,gam,				 
     #         c20,c22,c40,c42,c44,lbl,lbg,dc20,b0,gamman,	 
     #         cgtol,shiftr,shifti,range,fldi,dfld,
     #         AmpFac,SpcShift,ANGFI,ANGKSI,
     #         expl,expkxy,expkzz,anglib(6),reserved,Aadd(2*MaxAddNucl),
     #         in2,ipdf,ist,ml,mxy,mzz,lemx,lomx,kmn,kmx,mmn,mmx,ipnmx,
     #         nort,nstep,nfld,ideriv,ItypeSite,IFLGGIB,nresrv1,nresrv2,
     #	       ndim,ipt,jkmn,jmmn,itm,itd,itw,ipsi0,lband,kband,ldelta,
     #         kdelta,lptmx,kptmx,
     #         keylib,NAddN,in2Add(MaxAddNucl)
cAXV    The following block concerns to all spectra/sites, so momdls 
cAXV    recalculates these values many times!!!
      common /eprwork/ a0,g0,w0,faa(5),fam(2,5),fad(2,5),fgm(5),
     #         fgd(2,5),fww(5),fwm(2,5),fwd(2,5),cpot(5,5),xlk(5,5),
     #         g0first, fldShftSite(MXSITE)	
c
c *** The following constants identify the position of several
c     important parameters within the fepr (and fparm) arrays.
c     THESE CONSTANTS *MUST* BE REDEFINED IF THE PARAMETER ORDER 
c     IN EPRPRM IS CHANGED!!!
c
      integer IPHASE,IGIB0,IGIB2,IGIBX,IGXX,IGYY,IGZZ,IAXX,IAYY,IAZZ,
     #        IWXX,IWYY,IWZZ,IDX,IDY,IDZ,IPML,IPMXY,IPMZZ,IDJF,IDJFPRP,
     #        IOSS,IPSI,IALD,IBED,IGAD,IALM,IBEM,IGAM,IC20,IC44,ILBL,
     #        ILBG,IB0,IGAMAN,IML,ILEMX,IIPDF,IIST,IIN2,INDIM,
     #        INORT,INSTEP,ISHIFTR,ICGTOL,IFLDI,IDFLD,INFLD,IIDERV,
     #        IItypeSite,IIFLGGIB,IDC20,IRANGE,IIPT,IAMPFAC,ISPSHFT,
     #        IANGKSI,iNAddN,iin2Add,IAadd,IKEYLIB,IANGLIB,IEXPL,
     #        IANGFI,IEXPKXY,IEXPKZZ,IJKMN,IJMMN,IITM,IITD,IITW,IIPSI0,
     #        Ireserved,Ilband,Ikband,Ildelta,Ikdelta,Ilptmx,Ikptmx
c
cAXV names of parameters for fparms(*,js,jfile) and for iflvar(*,js,jfile):
      parameter(IPHASE=1,IGIB0=2,IGIB2=3,IGIBX=4,IWXX=5,IWYY=6,IWZZ=7,
     #         IGXX=8,IGYY=9,IGZZ=10,IAXX=11,IAYY=12,IAZZ=13,IDX=14,
     #         IDY=15,IDZ=16,IPML=17,IPMXY=18,IPMZZ=19,IDJF=20,
     #         IDJFPRP=21,IOSS=22,IPSI=23,IALD=24,IBED=25,IGAD=26,
     #         IALM=27,IBEM=28,IGAM=29,IC20=30,IC44=34,ILBL=35,
     #         ILBG=36,IDC20=37,IB0=38,IGAMAN=39,ICGTOL=40,ISHIFTR=41,
     #         IRANGE=43,IFLDI=44,IDFLD=45,IAMPFAC=46,ISPSHFT=47,
     #         IANGFI=48,IANGKSI=49,IEXPL=50,IEXPKXY=51,IEXPKZZ=52,
     #         IANGLIB=53,Ireserved=59,IAadd=60)
cAXV ( IGIBX=4) and  fparms(IGIBX,isi,isp) are used now
cAXV ( ISPSHFT=47) and  fparms(ISPSHFT,isi,isp) are used now
cAXV       number of real parameters NFPRM=60+MaxAddNucl*2
c
cAXV names of parameters for iparms(*,js,jfile)    
      parameter(IIN2=1,IIPDF=2,IIST=3,IML=4,ILEMX=7,INORT=14,INSTEP=15,
     #         INFLD=16,IIDERV=17,IItypeSite=18,IIFLGGIB=19,
cAXV    IIAFLG=20,IIRFLG=21, ! do not used!
     #         INDIM=22,IIPT=23,IJKMN=24,IJMMN=25,IITM=26,
     #         IITD=27,IITW=28,IIPSI0=29, Ilband=30,Ikband=31,
     #         Ildelta=32,Ikdelta=33,Ilptmx=34,Ikptmx=35,
     #         IKEYLIB=36,iNAddN=37,iin2Add=38)
cAXV       number of integer parameters NIPRM=38+MaxAddNucl,they are defined
cAXV       in dimnls.inc
c
cAXV     ILEMX+1  - lomx
cAXV     ILEMX+2  - kmn 
cAXV     ILEMX+3  - kmx
cAXV     ILEMX+4  - mmn 
cAXV     ILEMX+5  - mmx
cAXV     ILEMX+6  - ipnmx

