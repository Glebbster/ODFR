c based on parcom.inc from NLSL Version 1.5.1b 11/6/95
c----------------------------------------------------------------------
c                    =========================
c                         common PARCOM
c                    =========================
c
c     Common block containing a list of quantities associated with
c     each of the fitting parameters being varied (i.e. each parameter
c     in the x array for LMDER-family programs).
cc
c prscl  : Vector containing desired absolute accuracy for each parameter
c          in x array
c xfdstp : Vector containing size of forward-difference step for each
c          parameter in x array
c xerr   : Vector containing uncertainty estimates for each parameter in x
c          0=none, 1=minimum, 2=maximum, 3=both
c ixpr   : Index of each variable in the xfit array which indicates its number
c          into the fepr array in /eprprm/ or in array paramODF
c ixst1  : index of each parameter in the xfit array 
c          identifying which spectrum the parameter	is associated with
c ixst2  : index of each parameter appearing in the xfit array 
c          identifying which site in the spectrum the parameter
c          is associated with
c
c ixx    : =iflvar for this parameter, negative value 
c          means ODF parameter 
c nnsite(isp)  : Number of sites defined for a given spectrum
c nprm   : Number of parameters which are varied for all spectra together
c nprm1  : Number of parameters being varied directly
c nprm2  : Number of parameters which are varied independently for distinct spectrum
c nprm3	 : Number of parameters which are not varied but are changed
c nprm4	 : Number of parameters for ODF being varied
c nprmODF: Total number of ODF parameters
c----------------------------------------------------------------------
c
      double precision prscl,xfdstp,xerr
      integer ixpr,ixst1,ixst2,ixx, nnsite,
     #        nprm,nprm1,nprm2,nprm3,nprm4,nprmODF
c
      common /parcom/prscl(MXVAR),xfdstp(MXVAR),xerr(MXJCOL),
     #                ixx(MXVAR),nnsite(MXSPC),
     #    ixpr(MXVAR*MXSITE),ixst1(MXVAR*MXSITE),ixst2(MXVAR*MXSITE),
     #    nprm,nprm1,nprm2,nprm3,nprm4,nprmODF(MXSITE)
c
