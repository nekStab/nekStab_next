c-----------------------------------------------------------------------
      subroutine nekStab_setDefault
!     specifying default values for nekStab

      implicit none
      include 'SIZE'
      include 'TOTAL'

      k_dim = 100               ! standard value, increas  in .usr 
      schur_tgt = 2             ! schur target for schur step factorizaiton
      eigen_tol = 1.0e-6        ! tolerance for eigenmodes convergence
      schur_del = 0.10d0        !
      maxmodes = 20             ! max number of converged modes to disk
      glob_skip = 10            ! global energy computation skip frequency

      bst_skp = 10              ! boostconv skip iterations
      bst_snp = 10              ! bootsconv residual subspace matrix size

      ifres  = .false.          ! outpost restart files (KRY*, HES*)
      ifvor  = .false.          ! outpost vorticity (vor* omega_x,omega_y,omega_z components)
      ifvox  = .false.          ! outpost vortex (vox*: q,lambda2,omega criterions)
      ifldbf = .true.           ! load base flow for stability computations
      ifbf2D = .false.          ! force 2D base flow solution
      ifstorebase = .true.      ! store base flow for Floquet analysis (dynamic allocated)
      ifdyntol = .false.        ! dynamical tolerances for SFD and Newton (potential speed-up)

      ifseed_nois = .true.      ! noise as initial seed
      ifseed_symm = .false.     ! symmetry initial seed
      ifseed_load = .false.     ! loading initial seed (e.g. Re_ )
!     if all false the 'useric' subroutine prescribes the initial seed

!     position for zero-crossing vertical velocity check ! 
      xck = 2.0D0  ; call bcast(xck, wdsize)
      yck = 0.0D0  ; call bcast(yck, wdsize)
      zck = 0.0D0  ; call bcast(zck, wdsize)

      xLspg   = 0.0d0; call bcast(xLspg, wdsize) ! x left
      xRspg   = 0.0d0; call bcast(xRspg, wdsize) ! x right
      yLspg   = 0.0d0; call bcast(yLspg, wdsize)
      yRspg   = 0.0d0; call bcast(yRspg, wdsize)
      zLspg   = 0.0d0; call bcast(zLspg, wdsize)
      zRspg   = 0.0d0; call bcast(zRspg, wdsize)
      acc_spg = 0.333d0; call bcast(acc_spg, wdsize) !percentage for the acceleration phase in the sponge (e.g. 1/3)
      spng_str = 0.0d0;  call bcast(spng_str, wdsize)

      evop = '_'

!     !Broadcast all defaults !
      call bcast(schur_tgt, isize) ! isize for integer
      call bcast(eigen_tol, wdsize) ! wdsize for real
      call bcast(schur_del, wdsize)
      call bcast(maxmodes, isize)
      call bcast(k_dim, isize)
      call bcast(bst_skp, isize)
      call bcast(bst_snp, isize)
      call bcast(glob_skip, isize)

      call bcast(ifres   , lsize) !lsize for boolean
      call bcast(ifvor   , lsize)    
      call bcast(ifvox   , lsize)
      call bcast(ifseed_nois  , lsize)
      call bcast(ifseed_symm  , lsize)
      call bcast(ifseed_load  , lsize)
      call bcast(ifldbf  , lsize)
      call bcast(ifbf2D  , lsize)
      call bcast(ifstorebase  , lsize)
      call bcast(ifdyntol  , lsize)

      return
      end subroutine nekStab_setDefault
c-----------------------------------------------------------------------
      subroutine nekStab_init
      use krylov_subspace
!     initialize arrays and variables defaults
      implicit none
      include 'SIZE'
      include 'TOTAL'
      logical ifto_sav, ifpo_sav
      real glmin,glmax
      integer i
      n = nx1*ny1*nz1*nelv

      call nekStab_setDefault
      call nekStab_usrchk       ! where user change defaults
      call nekStab_printNEKParams

      xmn = glmin(xm1,n); xmx = glmax(xm1,n)
      ymn = glmin(ym1,n); ymx = glmax(ym1,n)
      zmn = glmin(zm1,n); zmx = glmax(zm1,n)

      if (nid==0) then
         print *,'                 __   _____  __          __  '
         print *,'   ____   ___   / /__/ ___/ / /_ ____ _ / /_ '
         print *,'  / __ \ / _ \ / //_/\__ \ / __// __ `// __ \'
         print *,' / / / //  __// ,<  ___/ // /_ / /_/ // /_/ /'
         print *,'/_/ /_/ \___//_/|_|/____/ \__/ \__,_//_.___/ '  
         print *,'COPYRIGHT (c) 2020-2023 DynFluid Laboratoire Paris ',NSVERSION
         print *,'Nek5000 ', NVERSION
         print *,''
      endif

      call copy(bm1s, bm1, n)   ! never comment this !

      if(spng_str.ne.0)then !sponge on

         if(nid.eq.0)write(6,*)
         if(nid.eq.0)write(6,*)' Initializing sponge...'
         if(nid.eq.0)write(6,*)' Sponge strenght:',spng_str
         if(spng_str.lt.0)then
           spng_str=abs(spng_str) 
           if(nid.eq.0)write(6,*)' Ensure positive sponge strenght:',spng_str
         endif
         call spng_init

!     applying sponge to the BM1 matrix to remove the sponge zone from eigensolver
         do i=1,n
            if( spng_fun( i ) .ne. 0 ) bm1s( i,1,1,1 )=0.0d0
         enddo

!     outposting BM1s to disk for check
!     ifto_sav = ifto; ifpo_sav = ifpo
!     ifvo=.false.; ifpo = .false.; ifto = .true.
!     call outpost(vx,vy,vz,pr,bm1s,'BMS')
!     ifvo=.true.; ifpo = ifpo_sav; ifto = ifto_sav

         if(nid.eq.0)write(6,*)'Sponge activated.'
         if(nid.eq.0)write(6,*)
      endif
      ifbfcv = .false.

      return
      end subroutine nekStab_init
c-----------------------------------------------------------------------
      subroutine nekStab
!     nekStab main driver
      implicit none
      include 'SIZE'
      include 'TOTAL'

      if(istep.eq.0)call nekStab_init

      call oprzero(fcx,fcy,fcz) ! never comment this!
      call rzero(fct,nx1*ny1*nz1*nelv)

      select case (floor(uparam(1)))

      case(3)                   ! eigenvalue problem

!     sanity check
         if    (uparam(1).eq.3.1)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for direct LNSE...'
         elseif(uparam(1).eq.3.11)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for direct LNSE in Floquet...'
         elseif(uparam(1).eq.3.2)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for adjoint LNSE...'
         elseif(uparam(1).eq.3.21)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for adjoint LNSE in Floquet...'           
         elseif(uparam(1).eq.3.3)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for transient growth...'
         elseif(uparam(1).eq.3.31)then
            if(nid.eq.0)write(6,*)'Krylov-Schur for transient growth in Floquet...'
         else
            if(nid.eq.0)write(6,*)'Krylov-Schur MODE NOT CORRECTLY SPECIFIED!'
            call nek_end
         endif

         call krylov_schur
         call nek_end


      end select

      return
      end subroutine nekStab
c-----------------------------------------------------------------------
      subroutine nekStab_outpost
!     nekStab custom outpost routine
      use krylov_subspace
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real vort(lv,3),wo1(lv),wo2(lv)
      common /ugrad/ vort,wo1,wo2

      logical ifto_sav, ifpo_sav

      if((istep.eq.0).OR.ifoutfld.AND.(ifvor.or.ifvox))then

         ifto_sav = ifto; ifpo_sav = ifpo
         ifto = .false.; ifpo = .false.

!---  > Compute and oupost vorticity.
         if(ifvor)then
            call oprzero(vort(1,1),vort(1,2),vort(1,3))
            call comp_vort3(vort,wo1,wo2,vx,vy,vz)
            call outpost(vort(1,1),vort(1,2),vort(1,3),pr,t, 'vor')
         endif

         ifto = ifto_sav ; ifpo = ifpo_sav

      endif
!---  > Outpost initial condition.
      if(istep.eq.0)call outpost(vx,vy,vz,pr,t,'   ')

      return
      end subroutine nekStab_outpost
c-----------------------------------------------------------------------
      subroutine nekStab_comment
!     Comment the evoltuion of the simulation
      implicit none
      include 'SIZE'
      include 'TOTAL'
      real*8,save :: eetime0,eetime1,eetime2
      data           eetime0,eetime1,eetime2 /0.0d0, 0.0d0, 0.0d0/
      real, save :: deltatime
      real telapsed,tpernondt,tmiss,dnekclock,ttime
      integer ttime_stp

!     if extrapolation is not OIFS -> ifchar = false
!     if OIFS active -> ifchar = .true. and CFL 2-5
!     some cases can have CFL>1 in initial time steps
      if (courno.gt.10) then
         if (nio.eq.0)then
            write(6,*)
            write(6,*)'    CFL > 10 stopping code'
            write(6,*)
         endif
         call nek_end
      endif

      if (nio.ne.0) return

      if (eetime0.eq.0.0 .and. istep.eq.1)then
         eetime0=dnekclock()
         deltatime=time
      endif
      eetime1=eetime2
      eetime2=dnekclock()

      if (istep.gt.0 .and. lastep.eq.0 .and. iftran) then

         ttime_stp = eetime2-eetime1 ! time per timestep
         ttime     = eetime2-eetime0 ! sum of all timesteps

         if(istep.eq.1)then
            ttime_stp = 0.0d0; ttime = 0.0d0
         endif

         if (mod(istep,5).eq.0) then

            telapsed = ttime/3600.0d0
            tpernondt = (ttime/(time-deltatime))
            tmiss = (param(10)-time)*tpernondt/3600.0d0

            print *,''
            write(6,"('      Mean time per timestep: ',F8.4,'  dev:',I8,'ms')") ttime/istep,ceiling(((ttime/istep)-ttime_stp)*1000) !to ms
            write(6,"('      Remaining time: ',I8,' h ',I2,' min')") int(tmiss),ceiling((tmiss-int(tmiss))*60.)
            if(tpernondt.gt.60.)then
               write(6,"('      Time per nondimensional time: ',F8.2,' sec')") tpernondt
            else
               write(6,"('      Time per nondimensional time: ',F8.2,' min ')") tpernondt/60.0d0
            endif
            write(6,"('      Local time: ',F10.4,'  File:',I8)") time-deltatime, int((time-deltatime)/param(14))+1
            print *,''

         endif
      endif

      return
      end subroutine nekStab_comment
c-----------------------------------------------------------------------
      subroutine nekStab_printNEKParams
!     print parameters at initialization for sanity check
      implicit none
      include 'SIZE'
      include 'TOTAL'
      if(nid.eq.0)then
         write(6,*)'P01=',param(1),'density'
         write(6,*)'P02=',param(2),'viscosity (1/Re)'
         write(6,*)'P07=',param(7),'rhoCp'
         write(6,*)'P08=',param(8),'conductivity (1/(Re*Pr))'
         write(6,*)'P10=',param(10),'stop at endTime'
         write(6,*)'P10=',param(11),'stop at numSteps'
         write(6,*)'P14=',param(14),'io step'
         write(6,*)'P15=',param(15),'io time'
         write(6,*)'P21=',param(21),'pressure sol tol'
         write(6,*)'P22=',param(22),'velocity sol tol'
         write(6,*)'P26=',param(26),'target CFL number'
         write(6,*)'P27=',param(27),'order in time'
         write(6,*)'P28=',param(28),'use same torder for mesh solver'
         write(6,*)'P31=',param(31),'numberOfPerturbations'
         write(6,*)'P41=',param(41),'1 for multiplicative SEMG'
         write(6,*)'P42=',param(42),'lin solv for the pres equation 0:GMRES,1:CG'
         write(6,*)'P43=',param(43),'0:additive multilevel scheme 1:orig 2lvl sch'
         write(6,*)'P44=',param(44),'0=E-based addit Schwarz PnPn-2;1=A-based'
         write(6,*)'P93=',param(93),'num vectors for projection'
         write(6,*)'P94 =',param(94),'num projection for helmholz solves'
         write(6,*)'P95=',param(95),'projection for pressure solver on/off'
         write(6,*)'P101=',param(101),'no additional modes'
         write(6,*)'P103=',param(103),'filter weight'
         write(6,*)
         write(6,*)'uparam01=',uparam(1)
         write(6,*)'uparam02=',uparam(02)
         write(6,*)'uparam03=',uparam(03)
         write(6,*)'uparam04=',uparam(04)
         write(6,*)'uparam05=',uparam(05)
         write(6,*)'uparam06=',uparam(06)
         write(6,*)'uparam07=',uparam(07)
         write(6,*)'uparam08=',uparam(08)
         write(6,*)'uparam09=',uparam(09)
         write(6,*)'uparam10=',uparam(10)
         write(6,*)
         write(6,*)'x min,max,tot=',xmn,xmx,xmx-xmn
         write(6,*)'y min,max,tot=',ymn,ymx,ymx-xmn
         write(6,*)'z min,max,tot=',zmn,zmx,zmx-zmn
         write(6,*)
      endif
      end subroutine nekStab_printNEKParams