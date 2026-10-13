C for ODFNLSL 19/11/2012
c*********************************************************************
c        Names of files
c      NamSpecInp - experimental spectra
c      NamSpecRes - results of simulation
c      ParInpEPR - initial parameters
c      ParResEPR - result parameters 
c      ParInpODF - initial parameters of ODF
c      ParResODF - result parameters of ODF
*********************************************************************

        character*30 NamSpecInp(MXSPC,MaxCALC),NamSpecRes(MXSPC,MaxCALC)
        character*30 ParInpEPR(MXSPC,MaxCALC),ParResEPR(MXSPC,MaxCALC)
        character*30 ParInpODF(MaxCALC),ParResODF(MaxCALC) 

       common /filenames/NamSpecInp,NamSpecRes,ParInpEPR,ParResEPR,
     # ParInpODF,ParResODF 
c
c*********************************************************************


