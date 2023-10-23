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
      spng_st = 0.0d0;  call bcast(spng_st, wdsize)

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
c-----------------------------------------------------------------------
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
                 write(6,*) 'NEWTON MODE NOT CORRECTLY SPECIFIED!'
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
         if(uparam(01) .eq. 4.41 .or. 
     $      uparam(01) .eq. 4.42) call ts_steady_force_sensitivity
         if(uparam(01) .eq. 4.43) call delta_forcing

         call nek_end

      end select

      return
      end subroutine nekStab
c-----------------------------------------------------------------------
      ! this is the new krylov_schur routine    
      subroutine linear_stability_analysis
      use NekVectors
      use LinearOperators
      use LightKrylov, krylov_atol => atol
      implicit none
      include 'SIZE'
      include 'TOTAL'
      !> Exponential Propagator.
      class(exponential_prop), allocatable :: A
      !> Krylov subspace.
      class(real_nek_vector), allocatable :: X(:)
      !> Coordinates of the eigenvectors in the Krylov basis.
      !complex*16, allocatable, dimension(:,:) :: vecs
      complex :: eigvecs(k_dim, k_dim)
      !> Eigenvalues.
      complex :: eigvals(k_dim)
      !> Residual.
      real :: residuals(k_dim)
      !> Information flag.
      integer :: info
      !> Miscellaneous.
      integer :: i, j, k
      real :: alpha
      class(abstract_vector), allocatable :: wrk
      
      character(len=20) :: fich1, fich2, fich3, fmt2, fmt3, fmt4, fmt5, fmt6
      character(len=30) :: filename
      !>
      !A = exponential_prop(1.0_wp)
      ! --> Initialize Krylov subspace.
      allocate (X(1:k_dim + 1))
      ! call random_number(X(1)%x)
      ! alpha = X(1)%norm()
      ! call X(1)%scal(1.0_wp / alpha)
      !     ----- Loading baseflow from disk (optional) -----
      time = 0.0d0 ! might be overwritten by load_fld
      if (ifldbf) then !skip loading if single run
      if (nid .eq. 0) write (*, *) 'Loading base flow from disk:'
      write (filename, '(a,a,a)') 'BF_', trim(SESSION), '0.f00001'
      call load_fld(filename)
      if (nid .eq. 0) write (*, *) ' Number os scalars found (npscal): ', npscal
      if (nid .eq. 0) write (*, *) ' ifldbf done.'
      else
      if (nid .eq. 0) write (*, *) 'Baseflow prescribed by the useric function in the .usr'
      end if
      !     ----- Save baseflow to disk (recommended) -----
      call nopcopy(ubase, vbase, wbase, pbase, tbase, vx, vy, vz, pr, t)
      !     ----- Prepare stability parameters -----
      if (istep .eq. 0 .and. (isFloquetDirect .or. isFloquetAdjoint .or. isFloquetDirectAdjoint)) then
          param(10) = time  ! Update UPO period in the field
          if (nid .eq. 0) then
              write (6, *) 'Floquet mode activated.'
              write (6, *) 'Getting endTime from file: endTime = ', param(10)
          end if
      end if
      ! Broadcast the UPO period to all processors
      call bcast(param(10), wdsize)
      !     ----- First vector (new from noise or restart) -----
      ! if (uparam(2) .eq. 0) then
      ! if (nid .eq. 0) write (6, *) 'Starting first Arnoldi decomposition...'
      ! !     ----- Creates seed vector for the Krylov subspace -----
      ! if (ifseed_nois) then    ! noise as initial seed
      !    if(nid.eq.0)write(6,*)'Filling fields with noise...'
      !    call op_add_noise(wrk2%vx,wrk2%vy,wrk2%vz)
      !    if (ifto) call add_noise_scal(wrk2%t(:,1),9.0e4, 3.0e3, 4.0e5)
      !    if (ldimt.gt.1) then
      !     do m = 2, ldimt
      !      if(ifpsco(m-1)) call add_noise_scal(wrk2%t(:,1),9.0e1*m, 3.0e2*m, 4.0e1*m)
      !     enddo
      !    endif
      !    call k_normalize(wrk2, alpha)
      !    call matvec(wrk, wrk2)
      ! elseif(ifseed_symm)then ! symmetry initial seed
      ! This should be loaded in useric instead for specific cases !
      !    if(nid.eq.0)write(6,*)'Enforcing symmetric seed perturb...'
      !    call add_symmetric_seed(wrk%vx,wrk%vy,wrk%vz,wrk%t(:,1))
      ! elseif (ifseed_load) then ! loading initial seed (e.g. dRe )
      ! if (uparam(01) .ge. 3.0 .and. uparam(01) .lt. 3.2) then
      ! write (filename, '(a,a,a)') 'dRe', trim(SESSION), '0.f00001'
      ! elseif (uparam(01) .ge. 3.2 .and. uparam(01) .lt. 3.3) then
      ! write (filename, '(a,a,a)') 'aRe', trim(SESSION), '0.f00001'
      ! end if
      ! if (nid .eq. 0) write (*, *) 'Load real part of mode 1 as seed: ', filename
      ! !    call load_fld(filename)
      ! !    call nopcopy(wrk2%vx,wrk2%vy,wrk2%vz,wrk2%pr,wrk2%t, vx,vy,vz,pr,t)
      ! !    call k_normalize(wrk2, alpha)
      ! !    call matvec(wrk, wrk2)
      ! else ! base flow is prescribed in usric
      ! !    call nopcopy(wrk%vx,wrk%vy,wrk%vz,wrk%pr,wrk%t, ubase,vbase,wbase,pbase,tbase)
      ! !    call k_normalize(wrk, alpha)
      ! end if
      ! --> Eigenvalue analysis.
      if (isDirect .or. isFloquetDirect) then
          evop = 'd'
          call eigs(A, X, eigvecs, eigvals, residuals, info, nev=schur_tgt, tol=eigen_tol, verbosity=.true., transpose=.false.)
      
      elseif (isAdjoint .or. isFloquetAdjoint) then
          evop = 'a'
          call eigs(A, X, eigvecs, eigvals, residuals, info, nev=schur_tgt, tol=eigen_tol, verbosity=.true. ,transpose=.true.)
      end if
      fich1 = 'Spectre_H'//trim(evop)//'.dat'
      fich2 = 'Spectre_NS'//trim(evop)//'.dat'
      fich3 = 'Spectre_NS'//trim(evop)//'_conv.dat'
      if (nid == 0) then
          open (unit=10, file=fich1, form='formatted')
          open (unit=20, file=fich2, form='formatted')
          open (unit=30, file=fich3, form='formatted')
          !outpost the eigenspectrum of the Hessenberg matrix.
          do i = 1, size(eigvals)
          write (10, "(3E15.7)") real(eigvals(i)), aimag(eigvals(i)), residuals(i)
          end do
          close (10)
      end if
      ! --> Transform eigenspectrum. (old outpost_ks)
      eigvals = log(eigvals)/A%t ! A%t is the sampling period dt*nsteps
      ! --> Save the eigenspectrum.
      if (nid .eq. 0) then
          do i = 1, size(eigvals)
          write (20, "(3E15.7)") real(eigvals(i)), aimag(eigvals(i)), residuals(i)
          end do
          close (20)
      end if
      ! --> Save leading eigenvector.
      call get_vec(wrk, X(1:k_dim), real(eigvecs(:, 1)))
      ! select type(wrk)
      ! type is(real_nek_vector)
      ! !call save_npy("Leading_eigenvector.npy", wrk%x)
      ! end select
      return
      end subroutine linear_stability_analysis