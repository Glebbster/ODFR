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

        integer, parameter :: MaxPath = 4096
        character(MaxPath) NamSpecInp(MXSPC,MaxCALC),NamSpecRes(MXSPC,MaxCALC)
        character(MaxPath) ParInpEPR(MXSPC,MaxCALC),ParResEPR(MXSPC,MaxCALC)
        character(MaxPath) ParInpODF(MaxCALC),ParResODF(MaxCALC)

       common /filenames/NamSpecInp,NamSpecRes,ParInpEPR,ParResEPR,
     # ParInpODF,ParResODF 
c
c*********************************************************************


