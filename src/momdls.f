cAXV based on MOMDLS.f from NLSL Version 1.5 beta 11/23/95
c    The orientation distribution explicitly presented as 
c     an expansion in a series of generalized spherical harmonics
c     is added.
c*********************************************************************
c
c                       ==================
c                        subroutine:MOMDLS
c                       ==================
c
c       Subroutine version of EPRLL family of programs.
c       This routine is intended for use with nonlinear least-squares
c       applications. The routine calculates a MOMD spectrum by
c       averaging spectra for a specified number of orientation using
c       the parameters passed in the fparm and iparm arrays. Spectra
c       are calculated using the Conjugate Gradients version of the Lanczos
c       algorithm.
c*********************************************************************
c
      subroutine momdls(fparm,iparm,icalc,spectr,nft,jsite,jsp,
     #                  numstor,tmpNams,DiscTrid,ierr,cache_mode,
     #                  reuseOneShot,first_site,last_site)
      use, intrinsic :: iso_fortran_env, only: int64
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      include 'mpif.h'
      include 'mpi_state.inc'
      include 'dimnls.inc'
      include 'eprprm.inc'
      include 'prmeqv.inc'
      include 'errmsg.inc'
      include 'pidef.inc'
      include 'rndoff.inc'
      include 'parODF.inc'      
c
      integer iparm(NIPRM),icalc,ierr,nft,jsite,jsp,ndone,	
     # numstor(MXSITE,MXSPC),zndone(MXSITE*MXSPC*MAXNORT)      
      integer ipnt,ioerr
      double precision fparm(NFPRM),spectr(MXSPT),
     #                 work(MXSPT),local_spectr(MXSPT)
      double complex store1(MXDIM,MAXStorMat),store2(MXDIM,MAXStorMat)
c   For preservation the values in course of
c   consecutive calling the subroutine lfunR.f:     	
      common /tempmomdls/store1,store2,zndone
c Euler angles SampleFrame -> Resonator frame:     
	 double precision alphaSF, betaSF, gammaSF
c
      double precision dlam,xi,cc,ss
      common /dfunc/ dlam,xi,cc,ss
c
      integer i,iort,id,j,nstp,nevl,jstore,iort_start,iort_step
      double precision cospsi,cosxi,delcos,dnorm,dwt,onorm,sinxi
      logical Full,DiscTrid,reuseOneShot,first_site,last_site
      integer cache_mode
c Each entry holds unweighted spectra for one calculator on this rank.
c CG indices refer to the existing slots during fitting and to saved
c coefficient prefixes when one-shot resets those slots each spectrum.
      type orientation_cache_entry
         integer(int64) real_key(NFPRM),field_shift_key
         integer int_key(NIPRM),fft_size
         double precision, allocatable :: spectra(:,:)
         integer cg_index(MAXNORT),cg_length(MAXNORT)
      end type
      type(orientation_cache_entry), save :: cache(MXSITE*MXSPC)
      integer, save :: entry_count,owner(MAXStorMat)
      integer, save :: full_count,hit_count
c At most 32 MiB per rank in addition to the saved spectra.
      integer, parameter :: max_saved_coefficients=1024*1024
      double complex, allocatable, save :: saved_coefficients(:,:)
      integer, save :: coefficients_used
      integer current_entry,nmatch,matches(MXSITE*MXSPC),local_count
c
      double precision ZERO,ONE,SMALL,D180,acc
      parameter(ZERO=0.0D0,ONE=1.0D0,SMALL=1.0D-16,D180=180.0D0,
     #  acc=1.0d-4)
c
      double precision fu20,fu20phi
      external fu20,fu20phi
c Function from WignerL.for for calculation of ODF value:	 
      external ODFL_RoReFast_Integrated
      double precision ODFL_RoReFast_Integrated
      double precision psiFr 
c
       character*15  tmpNams(MXSITE*MXSPC*MAXNORT)
c
cFor debugging:
cdb       double precision TempField(MXINP), TempSpec(MXINP)
cdb       integer flagname                  
c#####################################################################
c
      Full = mod(icalc,1000)/100 .ne. 0
      current_entry=0
      if (first_site) then
         entry_count=0
         owner=0
         coefficients_used=0
         full_count=0
         hit_count=0
      endif
c     ----------------------------------------------
c     Load values from parameter array into /eprprm/
c     ----------------------------------------------
      do i=1,NFPRM
         fepr(i)=fparm(i)
      end do
c
      do i=1,NIPRM
         iepr(i)=iparm(i)
      end do
c
      if (Full) then
         call tchange(jsite)
      else         
c  *** Nevetheless w0 which is used in eprls ***
         w0=(wxx+wyy+wzz)/3.0D0
      endif         
c
      do j=1,nfld
         spectr(j)=ZERO
         local_spectr(j)=ZERO
      end do
c
      nstp=nstep	 !user defined max number of steps for CG tridiagonalization
      if (nstep.le.0) nstp=MXDIM
c
      jstore=numstor(jsite,jsp)
      if (jstore.lt.1) jstore=1   
c
      if (nort.le.1 .or. ipt.le.0) then
c       --------------------------
c       --- Single orientation ---
c       --------------------------
         if ((.not.mpi_enabled).or.(mpi_rank.eq.0)) then
            call evaluate_orientation(1,local_spectr)
            if (ierr.ge.FATAL) then
               if (mpi_rank.eq.0) write(*,*) ' eprls gives FATAL error'
               if (mpi_enabled) call MPI_ABORT(MPI_COMM_WORLD,ierr,
     #              mpi_ierr)
               stop 1
            endif
         endif
         if (mpi_enabled) then
            call MPI_ALLREDUCE(local_spectr,spectr,nfld,
     #           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
         else
            do j=1,nfld
               spectr(j)=local_spectr(j)
            end do
         endif
c         
      else
c      --------------------------
c      ---- MOMD calculation ----
c      --------------------------
         onorm=ONE
         if (nort.gt.1) onorm=ONE/dfloat(nort-1)
c        ------------------------------------------------------
c        Partial MOMD
c        If there is ordering of the director, calculate 
c        normalization factor for the pseudopotential function
c        In the case of partial director ordering, the parameter
c        psi actually represents the angle between the director 
c        ordering axis and the field, whereas the calculation
c        subroutine always interprets psi as the director tilt.
c        Save the given psi value as the angle xi. 
c        ------------------------------------------------------
            dlam=ZERO
            dnorm=ONE
            xi=dmod(psi,D180)
         if (ODFtype(jsp).eq.1) then
            cosxi=dcos(RADIAN*xi)
            sinxi=dsin(RADIAN*xi)
            dlam=dc20
            call ccrint(ZERO,ONE,acc,SMALL,dnorm,nevl,fu20,id)
            if (dabs(dnorm).le.SMALL) dnorm=ONE
         end if
c        -----------------------
c        Loop over orientations
c        -----------------------
         cospsi=ZERO
         delcos=ZERO
         if (nort.gt.1) delcos=ONE/(nort-1)
c         
         iort_start=1
         iort_step=1
         if (mpi_enabled) then
            iort_start=1+mpi_rank
            iort_step=mpi_size
         endif
c
         if ((cache_mode.eq.1.or.cache_mode.eq.2).and.Full.and.
     #       iort_start.le.nort) call find_matching_calculators()

         do iort=iort_start,nort,iort_step
           if (mpi_enabled) then
              cospsi=min(dfloat(iort-1)*delcos,ONE)
           else
              cospsi=min(cospsi,ONE)
           endif
           psiFr=dacos(cospsi)
           fepr(IPSI)=psiFr/RADIAN
c
           call evaluate_orientation(iort,work)
c
            if (ierr.ge.FATAL) then 
               if (mpi_rank.eq.0) write(*,*) ' eprls gives FATAL error'
               if (mpi_enabled) call MPI_ABORT(MPI_COMM_WORLD,ierr,
     #              mpi_ierr)
               stop 1
            endif
c           --------------------------------------------------
c           If there is ordering of the director, calculate the
c           weighting coefficient for this tilt angle and
c           weight the spectrum
c           --------------------------------------------------
            if (ODFtype(jsp).eq.1) then
               if (xi.gt.RNDOFF) then
                  ss=sinxi*dsqrt(ONE-cospsi*cospsi)
                  cc=cosxi*cospsi
                  call ccrint(ZERO,PI,acc,SMALL,dwt,nevl,fu20phi,id)
                  dwt=dwt/(dnorm*PI)
               else
                  dwt=fu20(cospsi)/dnorm
               end if  
c
               do j=1,nfld
                  work(j)=work(j)*dwt
               end do
c			               
            end if   !ODFtype=1
c ODF
            if ((ODFtype(jsp).ge.2).and.(keyODF(jsite).ge.0)) then
c  
               alphaSF=0
               betaSF=RADIAN*xi
               gammaSF=0
c
               dwt=ODFL_RoReFast_Integrated
     *             (alphaSF,betaSF,gammaSF,psiFr)
               do j=1,nfld
                  work(j)=work(j)*dwt
               end do
cdb           write(*,*) jsp,iort,dwt              
cdb           call promout(jsp,iort,nfld,TempField,TempSpec)
             end if !	 ODFtype=2 or 3
c          --------------------------------------------------------
c          Trapezoidal rule: first and last pts. get factor of 1/2
c          --------------------------------------------------------
            if (iort.eq.1 .or. iort.eq.nort) then
               do j=1,nfld
                  work(j)=work(j)*0.5d0
               end do
            end if
c        ------------------------------------------------
c        Add spectrum for current orientation into total 
c        ------------------------------------------------
            do j=1,nfld
               local_spectr(j)=local_spectr(j)+work(j)
            end do
c
            if (.not.mpi_enabled) cospsi=cospsi+delcos
         end do !iort, loop over orientations 
c        -------------------
c        Normalize spectrum
c        -------------------
         do j=1,nfld
            local_spectr(j)=local_spectr(j)*onorm
         end do
         if (mpi_enabled) then
            call MPI_ALLREDUCE(local_spectr,spectr,nfld,
     #           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
         else
            do j=1,nfld
               spectr(j)=local_spectr(j)
            end do
         endif
c        -------------------------------------------
c        Return total number of steps actually taken
c        -------------------------------------------
      end if !(nort.le.1
c
      if (last_site.and.(cache_mode.eq.2.or.cache_mode.eq.3)) then
         write(*,'(A,I0,A,I0,A,I0,A,I0)') 'ORIENTATION_CACHE rank=',
     #        mpi_rank,' full=',full_count,' cfonly=',0,
     #        ' spectra=',hit_count
      endif
      return

      contains

      subroutine find_matching_calculators()
      double precision local_inputs(NFPRM)
      integer local_options(NIPRM),entry
      integer(int64) real_key(NFPRM),field_shift_key

c PSI becomes the local angle before EPRLS. DC20 and AMPFAC act on
c its result. ITYPESITE controls coefficient storage.
c All other prepared parameters must agree, including field and basis.
      local_inputs=fparm
      local_inputs(IPSI)=ZERO
      local_inputs(IDC20)=ZERO
      local_inputs(IAMPFAC)=ZERO
      if (.not.all(ieee_is_finite(local_inputs)).or.
     #    .not.ieee_is_finite(fldShftSite(jsite))) return
      local_options=iparm
      local_options(IItypeSite)=0
      real_key=transfer(local_inputs,0_int64,NFPRM)
      field_shift_key=transfer(fldShftSite(jsite),0_int64)
      nmatch=0
      do entry=1,entry_count
         if (all(cache(entry)%int_key.eq.local_options).and.
     #       all(cache(entry)%real_key.eq.real_key).and.
     #       cache(entry)%field_shift_key.eq.field_shift_key.and.
     #       cache(entry)%fft_size.eq.nft) then
            nmatch=nmatch+1
            matches(nmatch)=entry
         endif
      enddo
      entry_count=entry_count+1
      current_entry=entry_count
      local_count=1+(nort-iort_start)/iort_step
      cache(current_entry)%real_key=real_key
      cache(current_entry)%int_key=local_options
      cache(current_entry)%field_shift_key=field_shift_key
      cache(current_entry)%fft_size=nft
      cache(current_entry)%cg_index=0
      end subroutine find_matching_calculators

      subroutine evaluate_orientation(orientation,values)
      integer orientation,slot
      double precision values(MXSPT)

      slot=jstore+orientation-1
      call calculate_orientation(orientation,slot,store1(1,slot),
     #     store2(1,slot),zndone(slot),values)
      end subroutine evaluate_orientation

      subroutine calculate_orientation(orientation,slot,al,be,steps,values)
      integer orientation,slot,steps
      double complex al(MXDIM),be(MXDIM)
      double precision values(MXSPT)
      integer entry,src,src_index,match
      integer local_index,alloc_status,cg_index

      src=0
      if (current_entry.ne.0) then
         do match=1,nmatch
            entry=matches(match)
            src_index=cache(entry)%cg_index(orientation)
            if (src_index.eq.0) cycle
            if (.not.reuseOneShot) then
               if (owner(src_index).ne.cache_stamp(entry,orientation)) cycle
            endif
            src=entry
            exit
         enddo
      endif

      if (src.ne.0) then
         local_index=1+(orientation-iort_start)/iort_step
         values(1:nfld)=cache(src)%spectra(:,local_index)
         ierr=0
         call restore_coefficients(src,orientation,slot,al,be,steps)
         if (.not.reuseOneShot) owner(slot)=cache_stamp(src,orientation)
         hit_count=hit_count+1
         return
      endif

      call eprls(icalc,jsite,al,be,values,nft,nstp,steps,ierr)
      full_count=full_count+1

      if (Full.and..not.reuseOneShot) owner(slot)=0
      if (current_entry.eq.0) return
      if (ierr.ne.0) return
      if (.not.all(ieee_is_finite(values(1:nfld)))) return
      cg_index=slot
      if (reuseOneShot) then
         if (coefficients_used+steps.gt.max_saved_coefficients) return
         if (.not.allocated(saved_coefficients)) then
            allocate(saved_coefficients(max_saved_coefficients,2),
     #               stat=alloc_status)
            if (alloc_status.ne.0) return
         endif
         cg_index=coefficients_used+1
      endif
      if (allocated(cache(current_entry)%spectra)) then
         if (size(cache(current_entry)%spectra,1).ne.nfld.or.
     #       size(cache(current_entry)%spectra,2).ne.local_count)
     #       deallocate(cache(current_entry)%spectra)
      endif
      if (.not.allocated(cache(current_entry)%spectra)) then
         allocate(cache(current_entry)%spectra(nfld,local_count),
     #            stat=alloc_status)
         if (alloc_status.ne.0) return
      endif
      local_index=1+(orientation-iort_start)/iort_step
      cache(current_entry)%spectra(:,local_index)=values(1:nfld)
      cache(current_entry)%cg_index(orientation)=cg_index
      cache(current_entry)%cg_length(orientation)=steps
      if (reuseOneShot) then
         coefficients_used=coefficients_used+steps
         saved_coefficients(cg_index:coefficients_used,1)=al(1:steps)
         saved_coefficients(cg_index:coefficients_used,2)=be(1:steps)
      else
         owner(slot)=cache_stamp(current_entry,orientation)
      endif
      end subroutine calculate_orientation

      subroutine restore_coefficients(entry,orientation,slot,al,be,steps)
      integer entry,orientation,slot,steps
      double complex al(MXDIM),be(MXDIM)
      integer cg_index,cg_end

      cg_index=cache(entry)%cg_index(orientation)
      steps=cache(entry)%cg_length(orientation)
      if (reuseOneShot) then
c Preserve the writes FULL would make for linked sites and later jobs.
         cg_end=cg_index+steps-1
         al(1:steps)=saved_coefficients(cg_index:cg_end,1)
         be(1:steps)=saved_coefficients(cg_index:cg_end,2)
      else if (slot.ne.cg_index) then
         al(1:steps)=store1(1:steps,cg_index)
         be(1:steps)=store2(1:steps,cg_index)
      endif
      end subroutine restore_coefficients

      integer function cache_stamp(entry,orientation)
      integer entry,orientation
      cache_stamp=(entry-1)*MAXNORT+orientation
      end function cache_stamp
      end

c                         function FU20
c                    =========================
c  Function defining Boltzmann population at a given orientation
c  resulting from an orienting pseudopotential using only the c20
c  term. Used for calculating weighting as a function of angle
c  between the average director and the spectrometer field, for
c  the "partial MOMD" model. The potential is intended to describe
c  partial macroscopic orientation of the directors of local domains.
c
c  This is required to be a single-argument function so that it 
c  may be called by ccrint to calculate the normalization factor
c  for orientational weighting coefficients. The second parameter, the 
c  director c20 coefficient, is passed in common /dfunc/
c----------------------------------------------------------------------
      function fu20(cost)
      implicit none
      double precision fu20,cost
c
      double precision dlam,xi,cc,ss
      common /dfunc/ dlam,xi,cc,ss
c
      fu20 = exp( dlam*(1.5d0*cost*cost-0.5d0 ) )
      return
      end
c
c  phi is expected in radians
      function fu20phi(phi)
      implicit none
      double precision costht,fu20phi,phi
c
      include 'pidef.inc'
      double precision dlam,xi,cc,ss
      common /dfunc/ dlam,xi,cc,ss
c
      costht = ss*cos(phi)+cc
      fu20phi = exp( dlam*(1.5d0*costht*costht-0.5d0 ) )
      return
      end
