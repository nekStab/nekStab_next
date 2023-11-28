
subroutine harmonic_resolvent_precond

    use krylov_subspace
    implicit none
    include 'SIZE'
    include 'TOTAL'

!     ----- Krylov basis V and W for the projection M*V = V*H and M*W = W*H-----
    type(krylov_vector), allocatable, dimension(:) :: Q, W

!     ----- Upper Hessenberg matrix -----
    complex*16, allocatable, dimension(:,:) :: H
    complex*16, allocatable, dimension(:,:) :: P

!     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----
    complex*16, allocatable, dimension(:) :: vals
    complex*16, allocatable, dimension(:,:) :: vecs

    integer :: k_dim, i

!     ----- Allocate arrays -----
    allocate(Q(k_dim+1))
    allocate(W(k_dim+1))
    allocate(H(k_dim+1,k_dim), vals(k_dim), vecs(k_dim,k_dim))
    allocate(P(k_dim+1,k_dim))

    ! load base flow

    ! create initial seed




    ! rank-k_dim approximation creates a preconditioner P(mT) for M(mT)
    ! selecting a subset of r right eigenvectors (corresponding to the r least stable Floquet multipliers) and truncating the matrices V, W, and Λ to these r eigenvectors and corresponding eigenvalues. 
    
    ! preconditioner P(mT) = I - V_r e^Λr mT W*_r.

    ! inverse of the preconditioner, P(mT)^(-1)
    ! serves as an approximate inverse for M(mT), 
    ! is then derived using the matrix inversion lemma. 
    ! This is expressed as P(mT)^(-1) = I + V_r (e^-Λr mT - I)^(-1) W*_r. 
    ! This operator can be used to solve (3) for any periodic g(t)

    ! the inversion of the factor (e^-Λr mT - I) is an O(r) operation since 
    ! Λr is diagonal, which makes the process computationally efficient when k_dim is small


    ! spectral decomposition of Φ(T,0) = V e^ΛT W*
    ! Λ is a diagonal matrix of the Floquet exponents, 
    ! V is a matrix of right eigenvectors, 
    ! W* is a matrix of left eigenvectors (conjugate transpose of W)
    ! W*V = I - bi-orthogonal condition

    ! compute Arnoldi factorization for direct problem
    uparam(1) = 3.11; call bcast(uparam(1), wdsize)
    call arnoldi_factorization(Q, H, 1, k_dim, k_dim)
    
    ! compute eigenvalues and eigenvectors
    call eig(H(1:k_dim, 1:k_dim), vecs, vals, k_dim)

    ! get the NS eigenvalues 

    ! bi ortoghonalize the convergend ones (check ordering)

    ! compute Arnoldi factorization for adjoint problem
    uparam(1) = 3.21; call bcast(uparam(1), wdsize)
    call arnoldi_factorization(W, H, 1, k_dim, k_dim)

    
    call biorthogonalize(vx_dRe, vy_dRe, vz_dRe, pr, t, 
    $     vx_dIm, vy_dIm, vz_dIm, pr, t, 
    $     vx_aRe, vy_aRe, vz_aRe, pr, t, 
    $     vx_aIm, vy_aIm, vz_aIm, pr, t)
    
    ! Save P for use in GMRES computation, for example, by writing it to a file.
    ! Writing to a file is not included in this example and should be tailored to your specific needs.

    ! Deallocation
    deallocate(Q, W, H, vals, vecs, P)
    
end subroutine harmonic_resolvent_precond
