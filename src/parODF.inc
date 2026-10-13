c----------------------------------------------------------------------
c                    =========================
c                       common block parODF
c                    =========================
c
c      integer JTMP,ICON,js,ntemp,jord,jpar,ISTRING
      integer ODFtype,Jmx,keyODF,iflagODF,knet,PenFLG,ODFcontrol
      common/iparODF/ODFtype(MXSPC),Jmx(MXSITE),keyODF(MXSITE),
     # iflagODF(MaxParODF,MXSITE),knet,PenFLG,ODFcontrol
      integer ODFtypeAll
      common/iparODFAll/ODFtypeAll(MXSPC,MaxCALC)

      double precision  paramODF,PenLev,anegODF
      common/parODF/paramODF(MaxParODF,MXSITE),PenLev,anegODF(MXSITE)
