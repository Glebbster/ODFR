c----------------------------------------------------------------------
c                    =========================
c                       subroutine NLSINIT
c                    =========================
c
c----------------------------------------------------------------------
      subroutine nlsinit
      implicit none
c
      integer jpar,isi,isp,j,m,n,jODF,numODFs
c
      include 'dimnls.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'parODF.inc'
      include 'dimodfl.inc'

      include 'lmcom.inc'
      include 'iterat.inc'
        integer iparms(NIPRM,MXSITE,MXSPC)
        DOUBLE PRECISION fparms(NFPRM,MXSITE,MXSPC)
        COMMON /COMPARMS/ iparms,fparms
        integer IEQs
        COMMON /IEQ/ IEQs(MaxAddNucl,MXSITE,MXSPC)
      integer jsim,jpts,NptsExp       
      logical flfirst
      common/flfst/flfirst,jsim
       double precision FieldExp,SpecExp
	 common/iexp/NptsExp(MXSPC)
 	 common/exp/FieldExp(MXINP,MXSPC),SpecExp(MXINP,MXSPC)
       double precision Field(MXSPT,MXSPC)
       common /Field/ Field
      integer modtd
      logical depOK,useMOMD
        common/depOK/ depOK(MXVAR,MXSITE,MXSPC),modtd(MXSITE,MXSPC),
     #                useMOMD(MXSITE,MXSPC)        
      logical A_null
      dimension A_null(0:Jmax,-Jmax:Jmax,-Jmax:Jmax)
      common /A_jmn/A_jmn,jODF,numODFs
      common /WigODFIsNull/A_null
c
c----------------------------------------------------------------------
c    Initializations
c----------------------------------------------------------------------
      nspc=0
      ndatot=0
      nprm=0
      do isp=1,MXSPC
        do jpts=1,MXINP
           FieldExp(jpts,isp)=0.d0
         SpecExp(jpts,isp)=0.d0
        enddo
        do jpts=1,MXSPT
           Field(jpts,isp)=0.d0
        enddo
      enddo
c     Initialize parameter arrays
c----------------------------------------
      flfirst=.true.
      do isp=1,MXSPC
       do isi=1,MXSITE
         do jpar=1,NFPRM
           fparms(jpar,isi,isp)=0.0d0
         enddo !jpar
         do jpar=1,NIPRM
           iparms(jpar,isi,isp)=0
         enddo !jpar
         do jpar=1,MaxAddNucl
           IEQs(jpar,isi,isp)=0
         enddo !jpar
       enddo !isi
      enddo !isp
c
      do isp=1,MXSPC
         do isi=1,MXSITE
            modtd(isi,isp)=1
         end do
      end do
      do jpar=1,MXVAR
        do isp=1,MXSPC
           do isi=1,MXSITE
            depOK(jpar,isi,isp)=.false.
           end do
        end do
      enddo
        do isp=1,MXSPC
           do isi=1,MXSITE
            useMOMD(isi,isp)=.false.
           end do
        end do
c------------------------------------------------------------
c    Initialize data array parameters
c------------------------------------------------------------
      do isp=1,MXSPC
         nnsite(isp)=0
         npts(isp)=0
         ixsp(isp)=1
         idrv(isp)=1
         ODFtype(isp)=0
      end do
      knet=0
      PenFLG=0
      ODFcontrol=0
      PenLev=0.d0
      jODF=0
      numODFs=0
      do isi=1,MXSITE
         Jmx(isi)=0
         keyODF(isi)=0
         nprmODF(isi)=0
         anegODF(isi)=0.d0
         do jpar=1,MaxParODF
            paramODF(jpar,isi)=0.d0
            iflagODF(jpar,isi)=0
         enddo
      enddo
      do j=0,Jmax
         do m=-Jmax,Jmax
            do n=-Jmax,Jmax
               A_jmn(j,m,n)=(0.0d0,0.0d0)
               A_null(j,m,n)=.true.
            enddo
         enddo
      enddo
cAXV factor is used by minimizator (see comments in lmnls) 
      factor=1.0d2
      confid=0.683
      rdsqCom=0.d0
      lossSseW=0.d0
      lossSseUn=0.d0
      lossRedW=0.d0
      lossRedUn=0.d0
      rdchsqW=0.d0
      rdchsqUn=0.d0
      do isp=1,MXSPC
         rdsqSp(isp)=0.d0
         lossSpW(isp)=0.d0
         lossSpUn(isp)=0.d0
      enddo
c
      return
      end

c----------------------------------------------------------------------
c                       =========================
c                          function ITRIM
c                       =========================
c
c   Returns position of first blank in a string
c----------------------------------------------------------------------

      function itrim( string )
      implicit none
      integer itrim,j,lth
      character*(*) string
      lth=len(string)
      do 1 j=1,lth
         if (string(j:j).eq.' ') go to 2
 1    continue
      itrim=lth
      return
c
 2    itrim=j-1
      return
      end

c----------------------------------------------------------------------
