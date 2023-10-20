
    ! this is the new krylov_schur routine
subroutine linear_stability_analysis(adjoint)

    use LinearOperators

    !> Exponential Propagator.
    class(exponential_prop), allocatable :: A
    
    !> Krylov subspace.
    integer, parameter :: kdim = int(uparam(7))

    class(real_nek_vector), allocatable :: X(:)
    !> Coordinates of the eigenvectors in the Krylov basis.
    complex :: v(kdim, kdim)

    !> Eigenvalues.
    complex :: lambda(kdim)
    
    !> Residual.
    real :: residuals(kdim)

    !> Information flag.
    integer :: info

    !> Miscellaneous.
    integer :: i, j, k
    real :: alpha
    class(abstract_vector), allocatable :: wrk
    complex :: eigenvector(2*nx)
    logical :: adjoint
    character(len=20) :: fich1,fich2,fich3,fmt2,fmt3,fmt4,fmt5,fmt6

    !>
    A = exponential_prop(1.0_wp)

    ! --> Initialize Krylov subspace.
    allocate(X(1:kdim+1)) ; call random_number(X(1)%x)
    alpha = X(1)%norm() ; call X(1)%scal(1.0_wp / alpha)

    !     ----- Allocate arrays -----
    allocate(Q(k_dim+1))
    allocate(H(k_dim+1,k_dim),b_vec(1,k_dim),vals(k_dim),vecs(k_dim,k_dim),residual(k_dim))

    time   = 0.0d0
    H(:,:)  = 0.0d0; b_vec  = 0.0d0 ; residual = 0.0d0
    call k_zero(Q(1:k_dim+1))

!     ----- Loading baseflow from disk (optional) -----

    if(ifldbf)then            !skip loading if single run
       if(nid.eq.0)write(*,*)'Loading base flow from disk:'
       write(filename,'(a,a,a)')'BF_',trim(SESSION),'0.f00001'
       call load_fld(filename)
       if(nid.eq.0)write(*,*)' Number os scalars found (npscal): ',npscal
       if(nid.eq.0)write(*,*)' ifldbf done.'
    else
       if(nid.eq.0)write(*,*)'Baseflow prescribed by the useric function in the .usr'
    endif

!     ----- Save baseflow to disk (recommended) -----
    call nopcopy(ubase,vbase,wbase,pbase,tbase, vx,vy,vz,pr,t)

    ! --> Eigenvalue analysis.
    if (adjoint.eqv..true.) then
       evop = 'a'
       call eigs(A, X, v, lambda, residuals, info, nev=schur_tgt, transpose=.true.)
    else
       evop = 'd'
       call eigs(A, X, v, lambda, residuals, info, nev=schur_tgt, transpose=.false.)
    endif

    fich1 = 'Spectre_H' // trim(evop) // '.dat'
    fich2 = 'Spectre_NS' // trim(evop) // '.dat'
    fich3 = 'Spectre_NS' // trim(evop) // '_conv.dat'
    
    if (nid == 0) then
      open(unit=10, file=fich1, form='formatted')
      open(unit=20, file=fich2, form='formatted')
      open(unit=30, file=fich3, form='formatted')
    endif


    !outpost the eigenspectrum of the Hessenberg matrix.
    if (nid.eq.0) then
        do i = 1, size(lambda)
            write(10,"(3E15.7)") real(lambda(i)), aimag(lambda(i)), residuals(i)
        enddo
        close(10)
    endif

    ! --> Transform eigenspectrum. (old outpost_ks)
    lambda = log(lambda) / A%t ! A%t is the sampling period dt*nsteps

    ! --> Save the eigenspectrum.
    if (nid.eq.0) then
        do i = 1, size(lambda)
            write(20, "(3E15.7)") real(lambda(i)), aimag(lambda(i)), residuals(i)
        enddo
        close(20)
    endif

    ! --> Save leading eigenvector.
    call get_vec(wrk, X(1:kdim), real(v(:, 1)))
    select type(wrk)
    type is(real_nek_vector)
       !call save_npy("Leading_eigenvector.npy", wrk%x)
    end select

    return
  end subroutine linear_stability_analysis

!   subroutine transient_growth(times)
!     !> Time instants at which to compute the optimal growth.
!     real, intent(in) :: times(:)
!     !> Exponential propagator.
!     type(exponential_prop), allocatable :: A
!     !> Krylov subspaces.
!     integer, parameter :: kdim = 2*nx
!     type(vector), allocatable :: U(:), V(:)
!     !> Singular triplets.
!     real :: sigma(kdim), uvecs(kdim, kdim), vvecs(kdim, kdim)
!     real :: residuals(kdim)
!     !> Information flag.
!     integer :: info
!     !> Gains.
!     real :: gains(size(times), 5)
!     !> Miscellaneous.
!     integer :: i, j, k
!     real :: alpha

!     ! --> Initialize variables.
!     allocate(U(kdim+1)) ; allocate(V(kdim+1))
!     A = exponential_prop(0.0_wp)
!     gains = 0.0_wp ; sigma = 0.0_wp ; uvecs = 0.0_wp ; vvecs = 0.0_wp ; residuals = 0.0_wp

!     do i = 1, size(times)

!        write(*, *) "--: Integration time : ", times(i)

!        !> Set integration time for the exponential propagator.
!        A%t = times(i)

!        !> Initialize Krylov subspaces.
!        do j = 1, size(U)
!           call U(j)%zero() ; call V(j)%zero()
!        enddo
!        call random_number(U(1)%x) ; alpha = U(1)%norm() ; call U(1)%scal(1.0_wp / alpha)

!        !> Singular value computation.
!        call svds(A, U, V, uvecs, vvecs, sigma, residuals, info, nev=10, tolerance=rtol)

!        !> Store computed gains.
!        do k = 1, size(gains, 2)
!           gains(i, k) = sigma(2*k-1)**2
!        enddo
!        write(*, *) "    Optimal gains :", gains(i, :), sigma(1:2*size(gains, 2))**2
!        write(*, *)
!     enddo

!     call save_npy("Optimal_gains.npy", gains)

!     return
!   end subroutine transient_growth

!   subroutine resolvent_analysis(omegas)
!     !> Time instants at which to compute the optimal growth.
!     real, intent(in) :: omegas(:)
!     !> Exponential propagator.
!     type(resolvent_op), allocatable :: R
!     !> Krylov subspaces.
!     integer, parameter :: kdim = 2*nx
!     type(vector), allocatable :: U(:), V(:)
!     !> Singular triplets.
!     real :: sigma(kdim), uvecs(kdim, kdim), vvecs(kdim, kdim)
!     real :: residuals(kdim)
!     !> Information flag.
!     integer :: info
!     !> Gains.
!     real :: gains(size(omegas), 5)
!     !> Miscellaneous.
!     integer :: i, j, k
!     real :: alpha

!     ! --> Initialize variables.
!     allocate(U(kdim+1)) ; allocate(V(kdim+1))
!     R = resolvent_op(0.0_wp)
!     gains = 0.0_wp ; sigma = 0.0_wp ; uvecs = 0.0_wp ; vvecs = 0.0_wp ; residuals = 0.0_wp

!     do i = 1, size(omegas)

!        write(*, *) "--: Circular frequency : ", omegas(i)

!        !> Set the forcing frequency for the Resolvent.
!        R%omega = omegas(i)

!        !> Initialize Krylov subspaces.
!        do j = 1, size(U)
!           call U(j)%zero() ; call V(j)%zero()
!        enddo
!        call random_number(U(1)%x) ; alpha = U(1)%norm() ; call U(1)%scal(1.0_wp / alpha)

!        !> Singular value computation.
!        call svds(R, U, V, uvecs, vvecs, sigma, residuals, info, nev=10, tolerance=rtol)

!        !> Store computed gains.
!        do k = 1, size(gains, 2)
!           gains(i, k) = sigma(2*k-1)
!        enddo
!        write(*, *) "    Optimal gains :", gains(i, :)
!        write(*, *)
!     enddo

!     call save_npy("Resolvent_gains.npy", gains)

!     return
!   end subroutine resolvent_analysis

!   subroutine resolvent_single_freq(omega, elapsed_time)
!     !> Time instants at which to compute the optimal growth.
!     real, intent(in) :: omega
!     !> Elapsed time.
!     real, intent(out) :: elapsed_time
!     real              :: start_time = 0.0_wp, end_time = 0.0_wp
!     !> Exponential propagator.
!     type(resolvent_op), allocatable :: R
!     !> Krylov subspaces.
!     integer, parameter :: kdim = 2*nx
!     type(vector), allocatable :: U(:), V(:)
!     !> Singular triplets.
!     real :: sigma(kdim), uvecs(kdim, kdim), vvecs(kdim, kdim)
!     real :: residuals(kdim)
!     !> Information flag.
!     integer :: info
!     !> Miscellaneous.
!     integer :: i, j, k
!     real :: alpha

!     ! --> Initialize variables.
!     allocate(U(1:kdim+1)) ; allocate(V(1:kdim+1))
!     R = resolvent_op(omega)
!     sigma = 0.0_wp ; uvecs = 0.0_wp ; vvecs = 0.0_wp ; residuals = 0.0_wp

!     !> Initialize Krylov subspaces.
!     do j = 1, size(U)
!        call U(j)%zero() ; call V(j)%zero()
!     enddo
!     call random_number(U(1)%x) ; alpha = U(1)%norm() ; call U(1)%scal(1.0_wp / alpha)

!     !> Singular value computation.
!     call cpu_time(start_time)
!     call svds(R, U, V, uvecs, vvecs, sigma, residuals, info, nev=10, tolerance=rtol)
!     call cpu_time(end_time)

!     !> Elapsed time.
!     elapsed_time = end_time - start_time

!     return
!   end subroutine resolvent_single_freq
