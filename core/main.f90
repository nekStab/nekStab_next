c---------------------------------------------------------------------
      subroutine nekStab_setDefault
      !     specifying default values for nekStab
      
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         k_dim = 100               ! standard value, increas  in .usr
         schur_tgt = int(2)        ! schur target for schur step factorizaiton
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
         spng_st = 0.0d0;  call bcast(spng_st, wdsize)
      
         evop = '_'
      
      !     !Broadcast all defaults !
         call bcast(eigen_tol, wdsize) ! wdsize for real
         call bcast(schur_del, wdsize)
      
         call bcast(schur_tgt, isize8) ! isize8 for integer
         call bcast(maxmodes, isize8)
         call bcast(k_dim, isize8)
         call bcast(bst_skp, isize8)
         call bcast(bst_snp, isize8)
         call bcast(glob_skip, isize8)
      
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
c---------------------------------------------------------------------
      subroutine nekStab_init
         use krylov_subspace
      !     initialize arrays and variables defaults
         implicit none
         include 'SIZE'
         include 'TOTAL'
         logical scal
         real glmin,glmax
         integer i
         nv = nx1*ny1*nz1*nelv
      
         call nekStab_setDefault
         call nekStab_usrchk       ! where user change defaults
         call nekStab_printNEKParams
      
         xmn = glmin(xm1,nv); xmx = glmax(xm1,nv)
         ymn = glmin(ym1,nv); ymx = glmax(ym1,nv)
         zmn = glmin(zm1,nv); zmx = glmax(zm1,nv)
      
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
      
         call copy(bm1s, bm1, nv)   ! never comment this !
      
         if(spng_st.ne.0)then !sponge on
      
            if(nid.eq.0)write(6,*)
            if(nid.eq.0)write(6,*)' Initializing sponge...'
            if(nid.eq.0)write(6,*)' Sponge strenght:',spng_st
            if(spng_st.lt.0)then
               spng_st=abs(spng_st)
               if(nid.eq.0)write(6,*)' Ensure positive sponge strenght:',spng_st
            endif
            call spng_init
      
      !     applying sponge to the BM1 matrix to remove the sponge zone from eigensolver
            do i=1,nv
               if( spng_fn( i ) .ne. 0 ) bm1s( i,1,1,1 )=0.0d0
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
      
         nof = 0
         scal = .false.
         do i = 1, size(ifpsco)
            if (ifpsco(i) .eqv. .true.) then
               scal = .true.
               nof = nof + 1
            endif
         enddo
         if(ifto.eqv..true..or.scal.eqv..true.)then
            if(nid.eq.0)write(6,*)'Scalars found:'
            if(nid.eq.0)write(6,*)' ifto=',ifto
            if(nid.eq.0)write(6,*)' ifpsco=',ifpsco
            if(ifto) nof = nof + 1
            if(nid.eq.0)write(6,*)'number of possible scalars (ldimt)=',ldimt
            if(nid.eq.0)write(6,*)'number of scalars (nof)=',nof, npscal
         endif
         return
      end subroutine nekStab_init
c---------------------------------------------------------------------
      subroutine nekStab
      !     nekStab main driver
         use LinearStab, only : linear_stability_analysis
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         if(istep.eq.0) call nekStab_init
      
         call oprzero(fcx,fcy,fcz) ! never comment this!
         call rzero(fct,nx1*ny1*nz1*nelv)
      
         select case (floor(uparam(1)))
      
          case(0)                   ! DNS
      
            call nekStab_outpost   ! outpost vorticity
            call nekStab_comment   ! print comments
            call nekStab_energy   (vx,vy,vz,t,'total_energy.dat',glob_skip)
            call nekStab_enstrophy(vx,vy,vz,t,'total_enstrophy.dat',glob_skip)
            if (lastep .eq. 1) call nek_end
      
          case(1)                   ! fixed points computation
      
            call nekStab_outpost   ! outpost vorticity
            call nekStab_comment   ! print comments
      
            if( uparam(1) .eq. 1.1)then
               call SFD
               if(uparam(5).eq.0)call nekStab_energy(vx,vy,vz,t,'total_energy.dat',glob_skip)
            elseif( uparam(1) .eq. 1.2)then
               if(nid.eq.0)write(6,*)'BOOSTCONV'
               call BoostConv
            elseif( uparam(1) .eq. 1.3)then
               if(nid.eq.0)write(6,*)'DMT'
               if(nid.eq.0)write(6,*)'stopping ! not yet ported to this version'; call nek_end
      ! call DMT
            elseif( uparam(1) .eq. 1.4)then
               if(nid.eq.0)write(6,*)'TDF'
               call TDF
            endif
      
            if(ifbfcv)call nek_end
      
          case(2) ! Newton-Krylov solver
      
      ! Initialize flags based on the value of uparam(1)
            isNewtonFP  = (uparam(1) .eq. 2.0)
            isNewtonPO  = (uparam(1) .eq. 2.1)
            isNewtonPO_T = (uparam(1) .eq. 2.2)
      
      ! Conditional statements for each Newton-Krylov case
            if (nid .eq. 0) then
               if (isNewtonFP) then
                  write(6,*) 'Newton-Krylov for fixed points...'
               elseif (isNewtonPO) then
                  write(6,*) 'Newton-Krylov for UPOs...'
               elseif (isNewtonPO_T) then
                  write(6,*) 'Newton-Krylov for forced UPOs...'
               else
                  write(6,*) 'Unrecognized option...'
                  call nek_end
               endif
            endif
      
      ! Proceed with Newton-Krylov computation
            call newton_krylov
            call nek_end
      
          case(3)                   ! eigenvalue problem
      
            isDirect = (uparam(1) .eq. 3.1)
            isAdjoint = (uparam(1) .eq. 3.2)
            isDirectAdjoint = (uparam(1) .eq. 3.3)
      
            isFloquetDirect = (uparam(1) .eq. 3.11)
            isFloquetAdjoint = (uparam(1) .eq. 3.21)
            isFloquetDirectAdjoint = (uparam(1) .eq. 3.31)
      
      ! Conditional statements for each case
            if (nid .eq. 0) then
               if (isDirect) then
                  write(6,*) 'Krylov-Schur for Direct LNSE...'
               elseif (isAdjoint) then
                  write(6,*) 'Krylov-Schur for Adjoint LNSE...'
               elseif (isDirectAdjoint) then
                  write(6,*) 'Krylov-Schur for Direct-Adjoint LNSE...'
               elseif (isFloquetDirect) then
                  write(6,*) 'Krylov-Schur for Direct LNSE in Floquet...'
               elseif (isFloquetAdjoint) then
                  write(6,*) 'Krylov-Schur for Adjoint LNSE in Floquet...'
               elseif (isFloquetDirectAdjoint) then
                  write(6,*) 'Krylov-Schur for Direct-Adjoint LNSE in Floquet...'
               else
                  write(6,*) 'Unrecognized option...'
                  call nek_end
               endif ! isDirect
            endif ! nid == 0
      
            call linear_stability_analysis
      !call krylov_schur
            call nek_end
      
          case(4)                   ! in postprocessing.f
      
            if(uparam(01) .eq. 4.0) then ! all
               call stability_energy_budget
               call wave_maker
               call bf_sensitivity
            endif
      
      !     -----> Direct mode kinetic energy budget.
            if(uparam(01) .eq. 4.1) call stability_energy_budget
      
      !     -----> Wavemaker computation.
            if(uparam(01) .eq. 4.2) call wave_maker
      
      !     -----> Baseflow sensitivity.
            if(uparam(01) .eq. 4.3) call bf_sensitivity
      
      !     -----> Sensitivity to steady force.
            if(uparam(01) .eq. 4.41 .or. uparam(01) .eq. 4.42) call ts_steady_force_sensitivity
            if(uparam(01) .eq. 4.43) call delta_forcing
      
            call nek_end
      
         end select
      
         return
      end subroutine nekStab
c---------------------------------------------------------------------