c based on expdat.inc from NLSL Version 1.5.1 beta 11/25/95
c----------------------------------------------------------------------
c                         ====================
c                          common block EXPDAT
c                         ====================
c     Experimental data used for least-squares fitting and associated
c     parameters. 
c

cAXV    fitdata  : Experimental data points after spline
c     rmsn  : Noise level for each spectrum
c     npts  : Number of data points in each spectrum 
c     nft   : Smallest power of 2 larger than npts
c     idrv  : Flag=1 if data are in 1st derivative form				  !?
c     ixsp  : Starting index of each spectrum in data array
c     nspc  : Number of spectra
c     ndatot: Total number of points (all spectra)
c
c----------------------------------------------------------------------
      double precision fitdata,rmsn,ShiftExp			
      integer npts,nft,idrv,ixsp,nspc, ndatot
c
      common /expdat/ fitdata(MXPT),ShiftExp(MXSPC),rmsn(MXSPC),
     #      nft(MXSPC),npts(MXSPC),idrv(MXSPC),ixsp(MXSPC),
     #                nspc,ndatot
c
      double precision ClSpecDat
      common /mspctr/ ClSpecDat(MXSPT,MXSPC,MXSITE+1)

