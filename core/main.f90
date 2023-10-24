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
      
      
      
      
      subroutine linear_stability_analysis
         use NekVectors
         use LinearOperators
         use LightKrylov, krylov_atol => atol
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         class(exponential_prop), allocatable :: A
         class(real_nek_vector), allocatable :: X(:)
         complex :: eigvecs(k_dim, k_dim), eigvals(k_dim)
         real :: residuals(k_dim)
         integer :: info
         logical :: transpose
            
         allocate(X(1:k_dim + 1))

         call prepare_base_flow
         call prepare_linearized_solver      
         call prepare_perturbation_seed(X(1))
      
         if      (isDirect .or. isFloquetDirect) then
            evop = 'd'
            transpose = .false.
         elseif (isAdjoint .or. isFloquetAdjoint) then
            evop = 'a'
            transpose = .true.
         end if
      
         A = exponential_prop(fintim)

         call eigs(A, X, eigvecs, eigvals, residuals, info)!, nev, tolerance, verbosity, transpose)
      
         call outpost_eigenspectrum(eigvals, residuals, 'Spectre_H'//trim(evop)//'.dat')
      
         eigvals = log(eigvals)/A%t
      
         call outpost_eigenspectrum(eigvals, residuals, 'Spectre_NS'//trim(evop)//'.dat')
         call outpost_eigenvectors()
      
      end subroutine linear_stability_analysis
      
      subroutine prepare_base_flow
         use NekVectors
         use LinearOperators
         use LightKrylov, krylov_atol => atol
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         character(len=30) :: filename
      
         time = 0.0d0 ! might be overwritten by load_fld
      
         if (ifldbf) then !skip loading if single run
      
            if (nid .eq. 0) write (*, *) ' Loading base flow from disk:'
            write (filename, '(a,a,a)') 'BF_', trim(SESSION), '0.f00001'
            call load_fld(filename)
      
            if (nid .eq. 0) write (*, *) ' Number os scalars found (npscal): ', npscal
            if (nid .eq. 0) write (*, *) ' ifldbf done.'
      
         else
      
            if (nid .eq. 0) write (*, *) 'Baseflow prescribed by the useric function in the .usr'
      
         end if

         call nopcopy(ubase, vbase, wbase, pbase, tbase, vx, vy, vz, pr, t)
      
         return
      end subroutine prepare_base_flow
      
      subroutine prepare_linearized_solver
      
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         call nekgsync
      
      ! here we are taking the time from the base flow file
         if (istep .eq. 0 .and. (isFloquetDirect .or. isFloquetAdjoint .or. isFloquetDirectAdjoint)) then
            param(10) = time  ! Update UPO period in the field
            if (nid .eq. 0) then
               write (6, *) 'Floquet mode activated.'
               write (6, *) 'Getting endTime from file: endTime = ', param(10)
            end if
         end if
      ! Broadcast the UPO period to all processors
         call bcast(param(10), wdsize)
      
      !     --> Force only single perturbation mode.
         if (param(31) .gt. 1) then
            write(6, *) 'nekStab not ready for npert > 1 -- jp loops not yet implemented. Stoppiing.'
            call nek_end
         endif
         param(31) = 1 ; npert = param(31)
      
      !     --> Force deactivate OIFS.
         if (ifchar) write(6, *) 'OIFS not working with linearized solver. Turning it off.'
         ifchar = .false. ; call bcast(ifchar, lsize)
      
      !     --> Enforce CFL target for EXTk
         if (param(26) .gt. 1.0) then
            write(6, *) "Forcing target CFL to 0.5!"
            param(26) = 0.5d0
         endif
      
      !     --> Set nsteps/endTime accordingly.
         if (param(10) .gt. 0) then
            if(nid.eq.0)write(6,*)'param(10),time=',param(10),time
            if(nid.eq.0)write(6,*)'endTime specified! Recomputing dt and nsteps to match endTime'
            call compute_cfl(ctarg, vx, vy, vz, 1.0d0) ! ctarg contains the sum ( ux_i / dx_i )
            if(nid.eq.0)write(6,*)'max spatial restriction:',ctarg
            dt = param(26) / ctarg ! dt given target CFL
            nsteps = ceiling(param(10) / dt) ! computing a safe value of nsteps
            dt = param(10) / nsteps ! reducing dt to match param(10)
            if(nid.eq.0)write(6,*)' new timeStep dt=',dt
            if(nid.eq.0)write(6,*)' new numSteps nsteps=',nsteps
            if(nid.eq.0)write(6,*)' resulting sampling period =',nsteps*dt
            param(12) = dt
            call compute_cfl(ctarg, vx, vy, vz, dt) ! C=sum(ux_i/dx_i)*dt
            if(nid.eq.0)write(6,*)' current CFL and target=',ctarg,param(26)
            lastep = 0             ! subs1.f:279
            fintim = nsteps*dt
         endif
      
      !     --> Force constant time step.
         param(12) = -abs(param(12))
      
      !     --> Broadcast parameters.
         call bcast(param, 200*wdsize)
      
         return
      end subroutine prepare_linearized_solver
      
      subroutine prepare_perturbation_seed(seed)
         use NekVectors
         use LinearOperators
         use LightKrylov, krylov_atol => atol
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
         class(real_nek_vector), intent(inout) :: seed
         character(len=30) :: filename
         integer :: m
      
         if (ifseed_nois) then ! noise as initial seed
            if(nid.eq.0)write(6,*)'Filling fields with noise...'
      
            call op_add_noise(seed%vx,seed%vy,seed%vz)
            if (ifto) call add_noise_scal(seed%t(:,1),9.0e4, 3.0e3, 4.0e5)
            if (ldimt.gt.1) then
               do m = 2, ldimt
                  if(ifpsco(m-1)) call add_noise_scal(seed%t(:,1),9.0e1*m, 3.0e2*m, 4.0e1*m)
               enddo
            endif
      
         elseif (ifseed_load) then ! loading initial seed (e.g. dRe )
      
            if (isDirect.or.isFloquetDirect) then
               write (filename, '(a,a,a)') 'dRe', trim(SESSION), '0.f00001'
            elseif (isAdjoint.or.isFloquetAdjoint) then
               write (filename, '(a,a,a)') 'aRe', trim(SESSION), '0.f00001'
            end if
      
            if (nid .eq. 0) write (*, *) 'Load real part of mode 1 as seed: ', filename
            call load_fld(filename)
            call nopcopy(seed%vx,seed%vy,seed%vz,seed%pr,seed%t, vx,vy,vz,pr,t)
      
         else ! base flow is prescribed in usric
            call nopcopy(seed%vx,seed%vy,seed%vz,seed%pr,seed%t, ubase,vbase,wbase,pbase,tbase)
         end if
      
      !call nopcopy(X(1)%vx,X(1)%vy,X(1)%vz,X(1)%pr,X(1)%t, wrk%vx,wrk%vy,wrk%vz,wrk%pr,wrk%t)
      
         return ! seed returns to X(1)
      end subroutine prepare_perturbation_seed
      
      
      subroutine outpost_eigenspectrum(eigvals, residuals, filename)
         implicit none
         include 'SIZE'
         complex, intent(in) :: eigvals(k_dim)
         real, intent(in) :: residuals(k_dim)
         character(len=*), intent(in) :: filename
         integer :: i

         open(unit=10, file=trim(filename), status='replace', form='formatted')
         do i = 1, size(eigvals)
            write(10, "(3E15.7)") real(eigvals(i)), aimag(eigvals(i)), residuals(i)
         end do
         close(10)

      end subroutine outpost_eigenspectrum

      subroutine outpost_eigenvectors(X, eigvecs, eigvals, residuals)
         use NekVectors
         use LinearOperators
         use LightKrylov, krylov_atol => atol
         implicit none
         include 'SIZE'
         include 'TOTAL'

         class(real_nek_vector), intent(in) :: X(1:k_dim + 1)
         complex, intent(in) :: eigvecs(k_dim, k_dim), eigvals(k_dim)
         real, intent(in) :: residuals(k_dim)
         class(abstract_vector), allocatable :: nek_vector

         integer :: i
         character(len=3)  :: nRe,nIm,nRv

         nRe = trim(evop) // 'Re'
         nIm = trim(evop) // 'Im'
         nRv = trim(evop) // 'Rv'
      
         do i = 1, maxmodes
      
            if(nid.eq.0)write(6,*)'Outposting eigenvector:',i!,'/',maxmodes
            if(nid.eq.0)write(6,*)'  sigma=',real(eigvals(i))
            if(nid.eq.0)write(6,*)'  omega=',aimag(eigvals(i))
      
      !     ----- Output the real part -----
            call get_vec(nek_vector, X(1:k_dim), real(eigvecs(:, i)))
            select type(nek_vector)
             type is(real_nek_vector)
               ! alpha = nek_vector%norm()
               ! call nek_vector%scal(1.0 / alpha)
               call nopcopy(vx,vy,vz,pr,t, nek_vector%vx,nek_vector%vy,nek_vector%vz,nek_vector%pr,nek_vector%t)
      !call nopcmult(vx,vy,vz,pr,t, beta)
               call outpost2(vx,vy,vz,pr,t, nof, nRe)
               call outpost_vort(vx,vy,vz,nRv)
            end select
      
      !     ----- Output the real part -----
            call get_vec(nek_vector, X(1:k_dim), imag(eigvecs(:, i)))
            select type(nek_vector)
             type is(real_nek_vector)
               call nopcopy(vx,vy,vz,pr,t, nek_vector%vx,nek_vector%vy,nek_vector%vz,nek_vector%pr,nek_vector%t)
      !call nopcmult(vx,vy,vz,pr,t, beta)
               call outpost2(vx,vy,vz,pr,t, nof, nIm)
            end select
      
         enddo ! i = 1, maxmodes

      end subroutine outpost_eigenvectors
