C for ODFR 01.12.21014
C Spectra OUTPUT
       SUBROUTINE SpecOUT(programName,FileName,jspec,NamSpec,
     *             ParEPR,NamParODF,Info)
c	 
       IMPLICIT DOUBLE PRECISION (A-H),DOUBLE PRECISION (O-Z)
C
       include 'dimnls.inc'
       include 'expdat.inc'
       include 'eprprm.inc'
       include 'parcom.inc'
       include 'iterat.inc'
      include 'errmsg.inc'
      include 'mpi_state.inc'
c
        integer iparms(NIPRM,MXSITE,MXSPC)
        DOUBLE PRECISION fparms(NFPRM,MXSITE,MXSPC)
        COMMON /COMPARMS/ iparms,fparms
c
       character*30 programName,FileName,NamSpec,ParEPR,NamParODF
        integer NptsExp
        double precision FieldExp,SpecExp
	 common/iexp/NptsExp(MXSPC)
 	 common/exp/FieldExp(MXINP,MXSPC),SpecExp(MXINP,MXSPC)
       integer jspec,isite,js,info
c
       double precision Field(MXSPT,MXSPC)
       common /Field/ Field
c       
      double precision fvecPure(MXPT),chisqrPure,rdchsqPure     
      common /rdchsqPure/rdchsqPure           
c
       character*16 AA0, AAexp,BBexp
       character*16 AAres,BBres(nnsite(jspec)+1)
c   ######################################################################
c   
       if (mpi_rank.ne.0) return
       OPEN ( UNIT=1,FILE=FileName)
      WRITE(1,115)programName,NamSpec,ParEPR,NamParODF
 115      Format(1x,A15,', ',A15,', ',A15,', ',A15)
              write(1,117) minerr(info)
 117	 format('MINPACK returned: ', a32)
      WRITE(1,116) lossSseW,lossSseUn,lossRedW,lossRedUn
 116      Format('C',1x,'SSE_W=',G15.8,1x,'SSE_UNW=',G15.8,1x,
     #'RED_W=',G15.8,1x,'RED_UNW=',G15.8)
      WRITE(1,118) jspec,npts(jspec),lossSpW(jspec),lossSpUn(jspec)
 118      Format('C',1x,'Spectrum=',I3,1x,'NPTS=',I7,1x,
     #'SSE_W_SP=',G15.8,1x,'SSE_UNW_SP=',G15.8)
      do js=1,nnsite(jspec)
      if (iparms(IKEYLIB,js,jspec).ne.0) then
     	  WRITE(1,119,IOSTAT=ICON)  
     & fparms(IGXX,js,jspec),fparms(IGYY,js,jspec), 
     & fparms(IGZZ,js,jspec),fparms(IAXX,js,jspec),
     & fparms(IAYY,js,jspec),fparms(IAZZ,js,jspec)
 119  format('  Eff. g:',3(F16.12,1x),/,'  Eff. A:',3(F16.12,1x))  
      endif
      enddo   
c
       AA0='             ,'     
            WRITE(1,*)'  Hexp       ,     Sexp     ,',
     *'     Hmod    ,      Smod     '
c 	 
	 DO 10 i=1,MAX(NptsExp(jspec),npts(jspec))
         if (i.gt.NptsExp(jspec)) then
            AAexp=AA0
            BBexp=AA0
         else 
            write(AAexp,501)FieldExp(i,jspec)
            write(BBexp,502)SpecExp(i,jspec)
         endif
         if (i.gt.npts(jspec)) then
            AAres=AA0 
            BBres(1)=AA0
               DO js=1,nnsite(jspec)
                  BBres(js+1)=AA0
               ENDDO
         else 
            write(AAres,501)real(Field(i,jspec))		!+fparms(ISPSHFT,1,jspec))
            write(BBres(1),502)real(ClSpecDat(i,jspec,1)) 
            if (nnsite(jspec).gt.1) then
              DO isite=1,nnsite(jspec)
               write(BBres(isite+1),502)real(ClSpecDat(i,jspec,isite+1))
              END DO
            else
               BBres(2)=''
            endif	!nnsite(jspec).gt.0
         endif ! Field(i,jspec).eq.0
 
           write(1,503)AAexp,BBexp,
     *  AAres,(BBres(js),js=1,nnsite(jspec)+1)
  10   CONTINUE
c       ******************************************************************
       CLOSE (UNIT=1)

501        FORMAT (1X,F10.4,' ')
502        FORMAT (1X,G14.5,' ')
503        FORMAT (1X,20(8(A16)))
c
        return
      end


