c cAXV based on LCHECK.f from NLSL Version 1.5 beta 11/23/95
cAXV Major part of LCHECK had shifted to file PARreadN
c----------------------------------------------------------------------
c                    =========================
c                        subroutine TCHANGE
c                    =========================
c       Uses:
c               cd2km
c----------------------------------------------------------------------
      subroutine tchange(jsite)
      implicit none
c
      include 'dimnls.inc'
      include 'eprprm.inc'
      integer i,j,jsite 
      double precision d2km(2,5,5)
c
      double precision ZERO,HALF,DSQ23
      parameter (ZERO=0.0D0,HALF=0.5D0)
      parameter (DSQ23=0.816496580927726d0)
      
cdb       double precision fgd1(5,5),fad1(5,5)
c
c######################################################################
c----------------------------------------------------------------------
c     initialize spherical tensors
c----------------------------------------------------------------------
      do i=1,5
         faa(i)=ZERO
         fgm(i)=ZERO
         fww(i)=ZERO
         do j=1,2
            fam(j,i)=ZERO
            fgd(j,i)=ZERO
            fad(j,i)=ZERO
            fwm(j,i)=ZERO
            fwd(j,i)=ZERO
         end do
      end do
c
c----------------------------------------------------------------------
c     Calculate spherical components of F_g and F_A
c----------------------------------------------------------------------
c                                        *** Electronic Zeeman ***
      g0=(gxx+gyy+gzz)/3.0D0
      fgm(1)=HALF*(gxx-gyy)*(b0/g0)
      fgm(3)=DSQ23*(gzz-HALF*(gxx+gyy))*(b0/g0)
      fgm(5)=fgm(1)
c      
c    For shifting the spectra of 2-nd ... sites       
      if (jsite.eq.1) then 
         g0first=g0
      else 
         fldShftSite(jsite)=(g0first-g0)*b0/g0      
      endif                  
c                                        *** Hyperfine ***
      a0=(axx+ayy+azz)/3.0D0
      faa(1)=HALF*(axx-ayy)
      faa(3)=DSQ23*(azz-HALF*(axx+ayy))
      faa(5)=faa(1)
c                                        *** Line-broadening ***
      w0=(wxx+wyy+wzz)/3.0D0
      fww(1)=HALF*(wxx-wyy)
      fww(3)=DSQ23*(wzz-HALF*(wxx+wyy))
      fww(5)=fww(1)
c----------------------------------------------------------------------
c  Set ipt and cpot array according to the potential coefficients given
c----------------------------------------------------------------------
         cpot(2,1)=c20
         cpot(2,2)=c22
         cpot(3,1)=c40
         cpot(3,2)=c42
         cpot(3,3)=c44
c----------------------------------------------------------------------
c     Set  A tensor in g frame
c----------------------------------------------------------------------
      if (itm .eq. 0) then
c                           *** no tilt: copy A tensor directly
         do j=1,5
            fam(1,j)=faa(j)
         end do
cAXV:
         if (itw.eq.1) then  !if Lorentz frame is a-frame
           do j=1,5
              fwm(1,j)=fww(j)
           end do
         endif
      else
c                          *** transform A tensor into g axis system
         call cd2km(d2km,alm,bem,gam)
         do i=1,5
            do j=1,5
                  fam(1,i)=fam(1,i)+d2km(1,i,j)*faa(j)
                  fam(2,i)=fam(2,i)+d2km(2,i,j)*faa(j)
            end do
         end do
c
         if (itw.eq.1) then  !if Lorentz frame is a-frame
            do i=1,5
               do j=1,5
                  fwm(1,i)=fwm(1,i)+d2km(1,i,j)*fww(j)
                  fwm(2,i)=fwm(2,i)+d2km(2,i,j)*fww(j)
               end do
            end do
         endif  !(itw.eq.1)
      end if
c----------------------------------------------------------------------
c     Set  g,A,W tensors in diffusion frame
c----------------------------------------------------------------------
      if (itd .eq. 0) then
c                           *** no tilt: copy tensors directly
cAXV:
         do i=1,5
            fgd(1,i)=fgm(i)
            fad(1,i)=fam(1,i)
            fad(2,i)=fam(2,i)
         end do
c
         if (itw.eq.0) then  !if Lorentz frame is g-frame =D-frame with no tilt
           do i=1,5
            fwd(1,i)=fww(i)
           end do
         endif !itw.eq.0
         if (itw.eq.1) then  !if Lorentz frame is a-frame
           do j=1,5
              fwd(1,j)=fwm(1,j)
              fwd(2,j)=fwm(2,j)
           end do
         endif !itw.eq.1
      else
c                    *** Transform A, g, W tensors into the diffusion frame
         call cd2km(d2km,ald,bed,gad)
           do i=1,5
              do j=1,5
                 fgd(1,i)=fgd(1,i)+d2km(1,i,j)*fgm(j)
                 fgd(2,i)=fgd(2,i)+d2km(2,i,j)*fgm(j)
                 fad(1,i)=fad(1,i)+d2km(1,i,j)*fam(1,j)
     #                            -d2km(2,i,j)*fam(2,j)
                 fad(2,i)=fad(2,i)+d2km(1,i,j)*fam(2,j)
     #                            +d2km(2,i,j)*fam(1,j)
              end do
           end do
cdb           
cdb          call cd2km(d2km,-gad,-bed,-ald)
cdb      do i=1,5
cdb         do j=1,2
cdb            fgd1(j,i)=ZERO
cdb            fad1(j,i)=ZERO
cdb         end do
cdb      end do
cdb           do i=1,5
cdb              do j=1,5
cdb                 fgd1(1,i)=fgd1(1,i)+d2km(1,i,j)*fgd(1,j)
cdb     #                               -d2km(2,i,j)*fgd(2,j)
cdb                 fgd1(2,i)=fgd1(2,i)+d2km(1,i,j)*fgd(2,j)
cdb     #                              +d2km(2,i,j)*fgd(1,j)
cdb                 fad1(1,i)=fad1(1,i)+d2km(1,i,j)*fam(1,j)
cdb     #                            -d2km(2,i,j)*fam(2,j)
cdb                 fad1(2,i)=fad1(2,i)+d2km(1,i,j)*fam(2,j)
cdb     #                            +d2km(2,i,j)*fam(1,j)
cdb              end do
cdb           end do           
           
         if (itw.eq.0) then	!!if Lorentz frame is g-frame
           do i=1,5
              do j=1,5
                 fwd(1,i)=fwd(1,i)+d2km(1,i,j)*fww(j)
                 fwd(2,i)=fwd(2,i)+d2km(2,i,j)*fww(j)
              end do
           end do
         endif  !	(itw.eq.0)
         if (itw.eq.1) then  !if Lorentz frame is A-frame
            do i=1,5
               do j=1,5
                  fwd(1,i)=fwd(1,i)+d2km(1,i,j)*fwm(1,j)
     #                             -d2km(2,i,j)*fwm(2,j)
                  fwd(2,i)=fwd(2,i)+d2km(1,i,j)*fwm(2,j)
     #                             +d2km(2,i,j)*fwm(1,j)
               end do
            end do
         endif !(itw.eq.1) 
      end if
c
      if (itw.eq.2) then  !Lorentz frame is D-frame - copy tensor w directly
            do j=1,5
              fwd(1,j)=fww(j)
            end do
      endif
c
      return
      end

