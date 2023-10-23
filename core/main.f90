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
