      program steps_probe
      use optimization_steps, only: steps_begin,steps_end
      include 'dimnls.inc'
      include 'expdat.inc'
      include 'mpi_state.inc'
      integer iv(62),ui(1),i,NptsExp,keyCov
      double precision v(500),x(2),ur(1),Field,FieldExp,SpecExp
      double precision stepJac,stepBound,PARFCTOL,PARXCTOL
      common /Field/ Field(MXSPT,MXSPC)
      common /iexp/ NptsExp(MXSPC)
      common /exp/ FieldExp(MXINP,MXSPC),SpecExp(MXINP,MXSPC)
      common /STEPJ/ stepJac,stepBound
      common /CONVERG/ PARFCTOL,PARXCTOL
      common /KeyCov/ keyCov
      external probe_residual,probe_unused
      mpi_rank=0
      mpi_size=1
      mpi_enabled=.false.
      nspc=2
      npts(1)=2
      npts(2)=1
      NptsExp(1)=2
      NptsExp(2)=1
      ndatot=3
      Field(1:2,1)=(/1.d0,2.d0/)
      Field(1,2)=3.d0
      FieldExp=Field
      SpecExp=0.d0
      stepJac=1.d-6
      stepBound=1.d0
      PARFCTOL=1.d-10
      PARXCTOL=1.d-10
      keyCov=1
      iv=0
      v=0.d0
      call dfault(iv,v)
      iv(21)=0
      x=(/-1.2d0,1.d0/)
      open(unit=71,file='evaluations.dat',status='replace')
      call steps_begin(3,2)
      call nl2sno(3,2,x,probe_residual,iv,v,ui,ur,probe_unused)
      call steps_end()
      close(71)
      open(unit=72,file='solution.dat',status='replace')
      write(72,*) iv(1),iv(6),iv(30),iv(31),iv(40),iv(41)
      write(72,'(2ES25.16E3)') x
      close(72)
      if (maxval(abs(x-1.d0)).gt.1.d-5) stop 1
      end

      subroutine probe_residual(n,p,x,nf,r,ui,ur,unused,iv)
      include 'dimnls.inc'
      include 'expdat.inc'
      include 'iterat.inc'
      integer n,p,nf,ui(*),iv(*)
      double precision x(p),r(n),ur(*)
      external unused
      r(1)=10.d0*(x(2)-x(1)**2)
      r(2)=1.d0-x(1)
c     Third observation has zero objective weight but a nonzero model.
      r(3)=0.d0
      ClSpecDat(1:2,1,1)=-r(1:2)
      ClSpecDat(1,2,1)=x(1)+x(2)
      lossSseW=sum(r*r)
      lossSseUn=lossSseW+(x(1)+x(2))**2
      write(71,'(4I8,5ES25.16E3)') nf,iv(1),iv(31),iv(38),x,r
      end

      subroutine probe_unused()
      end
