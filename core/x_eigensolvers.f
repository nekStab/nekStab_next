!-----------------------------------------------------------------------




      
      subroutine inner_product(alpha, px, py, pz, pp, qx, qy, qz, qp)

!     This function provides the user-defined inner product to be used throughout
!     the computation.
!
!     INPUTS
!     ------
!
!     px, py, pz : nek arrays of size lt = lx1*ly1*lz1*lelv.
!     Arrays containing the velocity fields of the first vector.
!
!     pp : nek array of size lt2 = lx2*ly2*lz2*lelv
!     Array containing the pressure field of the first vector.
!     
!     qx, qy, qz : nek arrays of size lt = lx1*ly1*lz1*lelv.
!     Arrays containing the velocity fields of the second vector.
!
!     qp : nek array of size lt2 = lx2*ly2*lz2*lelv
!     Array containing the pressure field of the second vector.
!
!     RETURN
!     ------
!
!     alpha : real
!     Value of the inner-product alpha = <p, q>.
!
!     Last edit : April 2nd 2020 by JC Loiseau.
      
      implicit none
      include "SIZE"
      include "TOTAL"

      integer, parameter :: lt = lx1*ly1*lz1*lelt
      integer, parameter :: lt2 = lx2*ly2*lz2*lelt

      real, dimension(lt), intent(in) :: px, py, pz
      real, dimension(lt), intent(in) :: qx, qy, qz
      real, dimension(lt2), intent(in) :: pp, qp

      real, intent(out) :: alpha
      real :: glsc3
      integer :: n

      n = nx1 * ny1 * nz1 * nelv

      alpha = glsc3(px, bm1, qx, n) + glsc3(py, bm1, qy, n)
      if (if3d) alpha = alpha + glsc3(pz, bm1, qz, n)
      
      return
      end




      
!--------------------------------------------------------------------------




      
      subroutine norm(qx, qy, qz, qp, alpha) ! Compute vector norm
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter               :: lt  = lx1*ly1*lz1*lelt
      integer, parameter               :: lt2 = lx2*ly2*lz2*lelt
      real, intent(in), dimension(lt)  :: qx, qy, qz
      real, intent(in), dimension(lt2) :: qp
      real, intent(out)                :: alpha

      call inner_product(alpha, qx, qy, qz, qp, qx, qy, qz, qp)
      alpha = dsqrt(alpha)

      return
      end




      
!----------------------------------------------------------------------




      
      subroutine normalize(qx, qy, qz, qp, alpha)

!     This function normalizes the state vector [qx, qy, qz, qp]^T where
!     qx, qy and qz are the streamwise, cross-stream and spanwise velocity
!     components while qp is the corresponding pressure field.
!
!     INPUTS / OUTPUTS
!     ----------------
!
!     qx, qy, qz : nek arrays of size lt = lx1*ly1*lz1*lelv.
!     Arrays storing the velocity components.
!
!     qp : nek array of size lt2 = lx2*ly2*lz2*lelv
!     Array storing the corresponding pressure field.
!
!     alpha : real
!     Norm of the vector.
!
!     Last edit : April 2nd 2020 by JC Loiseau.
      
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                  :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                  :: lt2 = lx2*ly2*lz2*lelt
      real, dimension(lt), intent(inout)  :: qx, qy, qz
      real, dimension(lt2), intent(inout) :: qp
      real, intent(out)                   :: alpha
      real                                :: beta
      integer n2

      n2 = nx2*ny2*nz2*nelt

!     --> Compute the user-defined norm.
      call norm(qx, qy, qz, qp, alpha)
      beta = 1.0d0/alpha

!     --> Normalize the vector.
      call opcmult(qx, qy, qz, beta)
      if(ifpo) call cmult(qp, beta, n2)

      return
      end
c-----------------------------------------------------------------------




      
      subroutine krylov_schur()

!
!
!     
      
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

!     ----- Krylov basis V for the projection M*V = V*H -----
      real, dimension(lt,k_dim+1)        :: qx, qy, qz
      real, dimension(lt2,k_dim+1)       :: qp

!     ----- Upper Hessenberg matrix -----
      real, dimension(k_dim+1,k_dim)     :: H
      real, dimension(1,k_dim)           :: b_vec

!     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----
      complex*16, dimension(k_dim)       :: vals
      complex*16, dimension(k_dim,k_dim) :: vecs

!     ----- Miscellaneous -----
      real wo1(lt),wo2(lt),wo3(lt)
      common /ugrad2/ wo1,wo2,wo3

      integer :: mstart, cnt
      real, dimension(k_dim)             :: residual
      real                               :: alpha, beta, dumm, glsc3, bcast_size
      logical                            :: converged
      integer                            :: have_converged, n, n2, i, j
      integer                            :: ifich3 = 30
      character(len=30)                  :: filename
      
      n      = nx1*ny1*nz1*nelt
      n2     = nx2*ny2*nz2*nelt
      time   = 0.0d0
      H  = 0.0d0
      b_vec  = 0.0d0
      residual = 0.0d0

      call oprzero(wo1, wo2, wo3)
      do i = 1, k_dim+1
         call oprzero(qx(:, i), qy(:, i), qz(:, i))
         if (ifpo) call rzero(qp(:, i), n2)
      enddo

!     ----- Loading baseflow from disk (optional) -----

      if(.not. ifbfcv)then !skip loading if single run
       if(nid.eq.0)write(6,*)'Loading base flow: ',bf_handle
       call load_fld(bf_handle)
      endif

!     ----- Save baseflow to disk (recommended) -----

      ifto = .true. ;ifpo = .true.
      !call comp_vort3(t(1,1,1,1,2),wo1,wo2,vx,vy,vz)
      call lambda2(t(1,1,1,1,2))
      call outpost(vx,vy,vz,pr,t,'BF_') !outpost for sanity check
      call opcopy(ubase,vbase,wbase,vx,vy,vz)

!     ----- First vector (new from noise or restart) -----

      if (uparam(2).eq.0) then

        if(nid.eq.0)write(6,*)'Starting new Arnoldi decomposition...'

!     ----- Creates seed vector for the Krylov subspace -----

        if(ifnois)then

         call add_noise(vxp,vyp,vzp)
         call normalize(vxp,vyp,vzp,prp,alpha)
         call matrix_vector_product(vxp,vyp,vzp,prp, vxp,vyp,vzp)

         else !ifnois.eq..false.

          call opcopy(vxp,vyp,vzp,ubase,vbase,wbase)

        endif

!     ----- Normalized to unit-norm -----

         call normalize(vxp,vyp,vzp,prp,alpha)

         mstart = 1; istep = 1; time = 0.0d0

         call opcopy(qx(:,1), qy(:,1), qz(:,1), vxp, vyp, vzp)
         if(ifpo) call copy(qp(:,1), prp(:,1), n2)
         
         call whereyouwant('KRY',1)
         ifto = .false.; ifpo= .true.; time = 0.0d0
         call outpost(qx(:,1), qy(:,1), qz(:,1), qp(:,1), t, 'KRY')

      elseif(uparam(2).gt.0)then

      mstart = uparam(2)

      if(nid.eq.0)then
            write(6,*)'Restarting from:',mstart
            write(6,'(a,a,i4.4)')' Loading Hessenberg matrix: HES',trim(SESSION),mstart
            write(filename,'(a,a,i4.4)')'HES',trim(SESSION),mstart

            open(67, file=trim(filename), status='unknown', form='formatted')
             if (k_dim.lt.mstart) then !subsampling
              do i = 1,k_dim+1
               do j = 1,mstart
                if (j.le.k_dim)read(67,"(1E15.7)") H(i,j)
                enddo
              enddo
             else
              read(67, *)  (( H(i, j), j=1, mstart ), i = 1, mstart+1 )
             endif
            close(67)
            write(6,*)'Broadcast H matrix to all procs...'
      endif !nid.eq.0
      call bcast(H,(k_dim+1)*k_dim*wdsize) !broadcast H matrix to all procs

      mstart=mstart+1  !careful here!
      ifpo = .true.; ifto = .false.
      call load_files(qx, qy, qz, qp, mstart, k_dim+1, 'KRY')
      !k_dim+1 is the dim of V_x

      !mstart = mstart+1

      endif !(uparam(02))

!     ======================================
!     =====                            =====
!     ===== Krylov-Schur decomposition =====
!     =====                            =====
!     ======================================

      converged = .false.
      do while ( .not. converged )
!     --> Arnoldi factorization.
         call arnoldi_factorization(qx, qy, qz, qp, H, mstart)
         
!     --> Compute the eigenspectrum of the Hessenberg matrix.
         call eig(H(1:k_dim, 1:k_dim), vecs, vals, k_dim)
         
!     --> Check the residual of the eigenvalues.
         residual = abs(H(k_dim+1, k_dim) * vecs(k_dim, :))
         cnt = count(residual .lt. eigen_tol)
         
!     --> Select whether to stop or apply Schur condensation depending on schur_tgt.
         select case (schur_tgt)
            
!     --> k-step Arnoldi factorization completed.
         case (:0)
            converged = .true.
            if (nid .eq. 0) write(6, *) 'Arnoldi factorization completed.'
            
!     --> Krylov-Schur factorization.
         case (1:)
            if (cnt .ge. schur_tgt) then ! Krylov-Schur factorization completed.
               converged = .true.
            else                ! Apply Schur condensation before restarting the factorization.
               if (nid .eq. 0) write(6, *) 'Starting Schur condensation phase.'
               call schur_condensation(mstart, H, qx, qy, qz, qp)
            endif
            
         end select
         
      enddo
      
!     ----- Output all the spectrums and converged eigenmodes -----
      if(nid.eq.0)write(6,*)'Exporting modes...'
      call outpost_ks( vals , vecs , qx , qy , qz , qp , residual )
      
      if(nid.eq.0)write(6,*)'Stopping code...'
      call exitt

      return
      end




      
!-----------------------------------------------------------------------




      
      subroutine arnoldi_factorization(qx, qy, qz, qp, H, mstart)

!     This function implements the k-step Arnoldi factorization of the linearized
!     Navier-Stokes operator. The rank k of the Arnoldi factorization is set as a user
!     parameter in x_SIZE.f (see parameter k_dim).
!
!     INPUT
!     -----
!
!     mstart : integer
!     Index at which to start the Arnoldi factorization. By default, it should be set to 1.
!     Note that it changes when the Arnoldi factorization is used as part of the Krylov-Schur
!     algorithm.
!
!     RETURNS
!     -------
!
!     qx, qy, qz : nek arrays of size (lx1*ly1*lz1*lelv, k_dim).
!     Arrays containing the various Krylov vectors associated to each velocity component.
!
!     qp : nek arrays of size (lx2*ly2*lz2*lelv, k_dim)
!     Arrays containing the various Krylov vectors associated to the pressure field.
!
!     H : k x k real matrix.
!     Upper Hessenberg matrix resulting from the Arnoldi factorization of the linearized
!     Navier-Stokes operator.
!
!     Last edit : April 3rd 2020 by JC Loiseau.
      
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelv
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelv

!     ----- Miscellaneous -----
      real                               :: beta
      integer                            :: mstep, mstart
      integer                            :: n, n2

!     ----- Krylov basis V for the projection MQ = QH -----
      real, dimension(lt,k_dim+1)        :: qx, qy, qz
      real, dimension(lt2,k_dim+1)       :: qp

!     ----- Orthogonal residual f = w - (Q,w)*Q -----
      real, dimension(lt)                :: f_xr, f_yr, f_zr
      real, dimension(lt2)               :: f_pr

!     ----- Upper Hessenberg matrix -----
      real, dimension(k_dim+1, k_dim)       :: H

      n  = nx1*ny1*nz1*nelv ; n2 = nx2*ny2*nz2*nelv

!     --> Initialize arrays.
      call oprzero(f_xr, f_yr, f_zr)
      if (ifpo) call rzero(f_pr, n2)

!     --> Arnoldi factorization.
      do mstep = mstart, k_dim
         if (nid.EQ.0) write(6,*) 'Begin iteration :', mstep , '/' , k_dim

!     --> Matrix-vector product f = M * v (e.g. calling the linearized Navier-Stokes solver).
         call matrix_vector_product(f_xr, f_yr, f_zr, f_pr, qx(:,mstep), qy(:,mstep), qz(:,mstep))
         
!     --> Update Hessenberg matrix and compute the orthogonal residual f.
         call update_hessenberg_matrix(
     $        H(1:mstep, 1:mstep),
     $        f_xr, f_yr, f_zr, f_pr,
     $        qx(:, 1:mstep), qy(:, 1:mstep), qz(:, 1:mstep), qp(:, 1:mstep),
     $        mstep)
         
!     --> Normalise the residual vector.
         call normalize(f_xr, f_yr, f_zr, f_pr, beta)
         
!     --> Update the Hessenberg matrix.
         H(mstep+1, mstep) = beta
         
!     --> Add the residual vector as the new Krylov vector.
         call opcopy(qx(:,mstep+1), qy(:,mstep+1), qz(:,mstep+1), f_xr, f_yr, f_zr)
         if(ifpo) call copy(qp(:,mstep+1), f_pr, n2)
         
!     --> Save checkpoint for restarting/run-time analysis.
         call arnoldi_checkpoint(f_xr, f_yr, f_zr, f_pr, H(1:mstep+1, 1:mstep), mstep)
         
         if(nid.EQ.0) write(6,*) "End of iteration:", mstep
      enddo
      
      return
      end




      




      
c-----------------------------------------------------------------------




      
      subroutine schur_condensation(mstart, H, qx, qy, qz, qp)

      implicit none
      include 'SIZE'
      include 'TOTAL'
      
c     =================================================
c     =====                                       =====
c     ===== Declaration of the required variables =====
c     =====                                       =====
c     =================================================

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

c     ----- Krylov basis V for the projection M*V = V*H -----

      real, dimension(lt,k_dim+1)        :: qx, qy, qz
      real, dimension(lt2,k_dim+1)       :: qp

c     ----- Upper Hessenberg matrix -----

      real, dimension(k_dim+1, k_dim)     :: H
      real, dimension(k_dim)           :: b_vec

c     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----

      complex*16, dimension(k_dim)       :: vals

c     ----- Miscellaneous -----

      integer :: mstart,n,n2,i,j
      real beta
      logical, dimension(k_dim)          :: selected

c     ----- Schur and Hessenberg decomposition -----

      real, dimension(k_dim,k_dim)       :: vecs

      n  = nx1*ny1*nz1*nelv ; n2 = nx2*ny2*nz2*nelv

!     --> Initialize arrays.
      b_vec = 0.0D0 ; b_vec(k_dim) = H(k_dim+1, k_dim)
      vals = ( 0.0D0 , 0.0D0 ) ; vecs = 0.0D0

!     --> Perform the Schur decomposition.
      call schur(H(1:k_dim, 1:k_dim), vecs, vals, k_dim)

!     --> Partition the eigenvalues in wanted / unwanted.
      call select_eigenvalues(selected, mstart, vals, schur_del, schur_tgt, k_dim)
      if ( nid.EQ.0 ) write(*,*) mstart, 'Ritz eigenpairs have been selected.'
      
!     --> Re-order the Schur decomposition based on the partition.
      call ordschur(H(1:k_dim, 1:k_dim), vecs, selected, k_dim)

!     --> Zero-out the unwanted blocks of the Schur matrix.
      H(1:mstart,mstart+1:k_dim) = 0.0D0
      H(mstart+1:k_dim+1, :)     = 0.0D0

!     --> Re-order the Krylov basis accordingly.
      qx(:, 1:k_dim) = matmul(qx(:, 1:k_dim), vecs)
      qy(:, 1:k_dim) = matmul(qy(:, 1:k_dim), vecs)
      if (if3d) qz(:, 1:k_dim) = matmul(qz(:, 1:k_dim), vecs)
      if (ifpo) qp(:, 1:k_dim) = matmul(qp(:, 1:k_dim), vecs)
      
!     --> Update the Schur matrix with b.T @ Q corresponding to
!     the residual beta in the new basis.
      b_vec = matmul(b_vec, vecs)      
      H(mstart+1, :) = b_vec

!     --> Add the last generated Krylov vector as the new starting one.
      mstart = mstart + 1
      call opcopy(qx(:,mstart),  qy(:,mstart),  qz(:,mstart),
     $     qx(:,k_dim+1), qy(:,k_dim+1), qz(:,k_dim+1))
      
      return
      end




      
c----------------------------------------------------------------------




      
      subroutine matrix_vector_product(fx, fy, fz, fp, qx, qy, qz)
      implicit none
      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'
      include 'ADJOINT'
      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt
      real, dimension(lt)                :: fx, fy, fz, qx, qy, qz
      real, dimension(lt2)               :: fp
      integer n, n2
      n = nx1*ny1*nz1*nelt; n2 = nx2*ny2*nz2*nelt

!     ----- Initial condition -----
      call opcopy(vxp, vyp, vzp, qx, qy, qz)

!     ----- Time-stepper matrix-vector product -----
      time = 0.0d0
      do istep = 1, nsteps

!     ----- Prepare the base flow to vx,vy,vz

        if(ifldbf)then

          if(nid.eq.0)write(6,*)'Copying base flow!'
          call opcopy(vx,vy,vz,ubase,vbase,wbase)

        else

          if(nid.eq.0)then
            if(ifbase.eqv..true.)then
              write(6,*)'Running DNS alongside stability!'
            else
             write(6,*)'Base flow prescribed in useric!'
            endif
          endif

        endif

!      ----- Integrate in time vxp,vyp,vzp on top of vx,vy,vz

         call nek_advance

      enddo

      call opcopy(fx, fy, fz, vxp, vyp, vzp)
      if(ifpo) call copy(fp, prp, n2)

      return
      end




      
c----------------------------------------------------------------------




      
      subroutine outpost_ks(vals, vecs, qx, qy, qz, qp, residual) !outposting vectors
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter                 :: lt  = lx1*ly1*lz1*lelt
      integer, parameter                 :: lt2 = lx2*ly2*lz2*lelt

c     ----- Krylov basis V for the projection M*V = V*H -----

      real wo1(lt),wo2(lt),wo3(lt),vort(lt,3)
      common /ugrad/ wo1,wo2,wo3,vort

      real, dimension(lt,k_dim+1)        :: qx, qy, qz
      real, dimension(lt2,k_dim+1)       :: qp

c     ----- Eigenvalues (VP) and eigenvectors (FP) of the Hessenberg matrix -----

      complex*16, dimension(k_dim)       :: vals
      complex*16, dimension(k_dim,k_dim) :: vecs

c     ----- Arrays for the storage/output of a given eigenmode of the NS operator -----

      complex*16, dimension(lt)          :: fp_cx, fp_cy, fp_cz
      complex*16, dimension(lt2)         :: fp_cp

c     ----- Miscellaneous -----

      integer :: ifich1 = 10, ifich2 = 20, ifich3 = 30
      integer :: n, n2, i, j

      real                               :: sampling_period
      real, dimension(k_dim)             :: residual
      real                               :: alpha, alpha_r, alpha_i
      complex :: log_transform

      sampling_period = dt*nsteps

      n = nx1*ny1*nz1*nelt ; n2 = nx2*ny2*nz2*nelt

c     ----- Output all the spectrums and converged eigenmodes -----

      if(nid.EQ.0) then
         open(unit=ifich1, file='Spectre_H.dat'      , form='formatted')
         open(unit=ifich2, file='Spectre_NS.dat'     , form='formatted')
         open(unit=ifich3, file='Spectre_NS_conv.dat', form='formatted')
      endif
      
      do i = 1, k_dim
         
         time = real(i)         !here for the outposted files have different times
         
         if(nid.EQ.0) then
            
!     --> Outpost the eigenspectrum of the Hessenberg matrix.
            write(ifich1,"(3E15.7)") real(vals(i)), aimag(vals(i)), residual(i)
            
!     --> Outpost the log-transform spectrum (i.e. eigenspectrum of the linearized Navier-Stokes operator).
            write(ifich2, "(3E15.7)") real(log_transform(vals(i))) / sampling_period,
     $           aimag(log_transform(vals(i))) / sampling_period,
     $           residual(i)
            
         endif
         
         if (residual(i) .lt. eigen_tol) then !just the converged oned
            
!     --> Outpost only the converged part of the log-transformed spectrum.
            if (nid.EQ.0) then
               write(ifich3, "(2E15.7)") real(log_transform(vals(i))) / sampling_period,
     $              aimag(log_transform(vals(i))) / sampling_period
               
            endif
            
c     ----- Computation of the corresponding eigenmode -----
            fp_cx = matmul(qx(:, 1:k_dim), vecs(:, i))
            fp_cy = matmul(qy(:, 1:k_dim), vecs(:, i))
            if (if3d) fp_cz = matmul(qz(:, 1:k_dim), vecs(:, i))
            if (ifpo) fp_cp = matmul(qp(:, 1:k_dim), vecs(:, i))
            
c     ----- Normalization to be unit-norm -----
c     Note: volume integral of FP*conj(FP) = 1.
            call norm(real(fp_cx), real(fp_cy), real(fp_cz), real(fp_cp), alpha_r)
            call norm(aimag(fp_cx), aimag(fp_cy), aimag(fp_cz), aimag(fp_cp), alpha_i)
            alpha = alpha_r**2 + alpha_i**2
            alpha = 1.0d0/sqrt(alpha)
            
c     ----- Output the imaginary part -----
            call opcopy(vx, vy, vz, aimag(fp_cx), aimag(fp_cy), aimag(fp_cz))
            if(ifpo) call copy(pr, aimag(fp_cp), n2)
            
            call opcmult(vx, vy, vz, alpha)
            if(ifpo) call cmult(pr, alpha, n2)
            ifto = .false. ; ifpo = .true.
            call outpost(vx, vy, vz, pr, t, 'Im_')
            
c     ----- Output the real part -----
            call opcopy(vx, vy, vz, real(fp_cx), real(fp_cy), real(fp_cz))
            if (ifpo) call copy(pr, real(fp_cp), n2)
            
            call opcmult(vx, vy, vz, alpha)
            if (ifpo) call cmult(pr, alpha, n2)
            ifto = .false. ; ifpo = .true.
            call outpost(vx, vy, vz, pr, t, 'Re_')
            
c     ----- Output vorticity from real part -----
            call oprzero(wo1, wo2, wo3)
            call comp_vort3(vort, wo1, wo2, vx, vy, vz)
            ifto = .false. ; ifpo = .false.
            call outpost(vort(1,1), vort(1,2), vort(1,3), pr, t, 'Rev')
            ifto = .true. ; ifpo = .true.
            
         endif
         
      enddo
      
      if (nid .eq. 0) then
         close(ifich1) ; close(ifich2) ;  close(ifich3)
      endif
      
      return
      end




      
c-----------------------------------------------------------------------




      
      subroutine select_eigenvalues(selected, cnt, vals, delta, nev, n)

!     This function selects the eigenvalues to be placed in the upper left corner
!     during the Schur condensation phase.
!
!     INPUTS
!     ------
!
!     vals : n-dimensional complex array.
!     Array containing the eigenvalues.

!     delta : real
!     All eigenvalues outside the circle of radius 1-delta will be selected.
!
!     nev : integer
!     Number of desired eigenvalues. At least nev+4 eigenvalues will be selected
!     to ensure "smooth" convergence of the Krylov-Schur iterations.
!
!     n : integer
!     Total number of eigenvalues.
!
!     RETURNS
!     -------
!
!     selected : n-dimensional logical array.
!     Array indicating which eigenvalue has been selected (.true.).
!
!     cnt : integer
!     Number of selected eigenvalues. cnt >= nev + 4.
!
!     Last edit : April 2nd 2020 by JC Loiseau.
      
      implicit none

!     ----- Input arguments -----
      integer :: nev
      integer :: n
      integer, dimension(n) :: idx
      complex*16, dimension(n) :: vals
      real :: delta

!     ----- Output argument -----
      logical, dimension(n) :: selected
      integer :: cnt

!     ----- Miscellaneous -----
      integer :: i

!     --> Sort eigenvalues based on their magnitudes (Note : the array vals itself is not sorted).
      do i = 1, n
         idx(i) = i
      enddo
      call quicksort2(n, abs(vals), idx)

!     --> Select eigenvalues closer to the unit circle.
      selected = abs(vals) .ge. (1.0d0-delta)

!     --> Select at least the nev+4 largest eigenvalues.
      selected(idx(n-(nev+3):n)) = .true.
      if (aimag(vals(idx(n-(nev+3)))) .eq. -aimag(vals(idx(n-(nev+4))))) then
         selected(idx(n-(nev+4))) = .true.
      endif

      cnt = count(selected)

      return
      end





!     -------------------------------------------------------------------------





      subroutine update_hessenberg_matrix(
     $     H,
     $     f_xr, f_yr, f_zr, f_pr,
     $     qx, qy, qz, qp,
     $     k)

!     This function orthonormalizes the latest Krylov vector f w.r.t. all of the
!     previous ones and updates the entries of the Hessenberg matrix accordingly.
!
!     INPUTS
!     ------
!
!     k : int
!     Current step of the Arnoldi factorization.
!
!     f_xr, f_yr, f_zr : nek arrays of size lt = lx1*ly1*lz1*lelt.
!     Velocity components of the latest Krylov vector.
!     When returned, it has been orthonormalized w.r.t. to all previous
!     Krylov vectors.
!
!     f_pr : nek array of size lt2 = lx2*ly2*lz2*lelt.
!     Pressure component of the latest Krylov vector.
!
!     qx, qy, qz : nek arrays of size (lt, k)
!     Velocity components of the Krylov basis.
!     
!     qp : nek array of size (lt2, k).
!     Pressure component of the Krylov basis.
!
!     H : k x k real matrix.
!     Upper Hessenberg matrix.
!
!     Last edit : April 3rd 2020 by JC Loiseau.
      
      implicit none
      include "SIZE"
      include "TOTAL"

      integer, parameter :: lt = lx1*ly1*lz1*lelv
      integer, parameter :: lt2 = lx2*ly2*lz2*lelv
      
      integer, intent(in) :: k
      real, dimension(k, k), intent(inout) :: H

      real, dimension(lt, k), intent(in) :: qx, qy, qz
      real, dimension(lt2, k), intent(in) :: qp

      real, dimension(lt), intent(inout) :: f_xr, f_yr, f_zr
      real, dimension(lt2), intent(inout) :: f_pr

      integer i, j, n, n2
      real alpha

      real, dimension(lt) :: work1, work2, work3
      real, dimension(lt2) :: workp
      real, dimension(k) :: h_vec

      n = nx1*ny1*nz1*nelv;      n2 = nx2*ny2*nz2*nelv

!     --> Initialize array.
      call rzero(h_vec, k)

!     --> Orthonormalize f w.r.t the Krylov basis.
      do i = 1, k

!     --> Copy the i-th Krylov vector to the working arrays.
         call opcopy(work1, work2, work3, qx(:, i), qy(:, i), qz(:, i))
         if (ifpo) call copy(workp, qp(:, i), n2)

!     --> Orthogonalize f w.r.t. to q_i.
         call inner_product(alpha, f_xr, f_yr, f_zr, f_pr, work1, work2, work3, workp)

         call opcmult(work1, work2, work3, alpha)
         if (ifpo) call cmult(workp, alpha, n2)

         call opsub2(f_xr, f_yr, f_zr, work1, work2, work3)
         if (ifpo) call sub2(f_pr, workp, n2)

!     --> Update the corresponding entry in the Hessenberg matrix.
         H(i, k) = alpha

      enddo
         
      return
      end





!     ------------------------------------------------------------------------------------





      subroutine arnoldi_checkpoint(f_xr, f_yr, f_zr, f_pr, H, k)

!     This function implements a fairly simple checkpointing procedure in case one
!     would need to restart the computation (e.g. in case of cluster shutdown).
!
!     INPUTS
!     ------
!
!     f_xr, f_yr, f_zr : nek arrays of size lt = lx1*ly1*lz1*lelv
!     Velocity components of the latest Krylov vector.
!
!     f_pr : nek array of size lt2 = lx2*ly2*lz2*lelv
!     Pressure field of the latest Krylov vector.
!
!     H : k+1 x k real matrix.
!     Current upper Hessenberg matrix resulting from the k-step Arnoldi factorization.
!
!     k : int
!     Current iteration of the Arnoldi factorization.
!
!     Last edit : April 3rd 2020 by JC Loiseau.

      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer, parameter :: lt = lx1*ly1*lz1*lelv
      integer, parameter :: lt2 = lx2*ly2*lz2*lelv

      real, dimension(lt), intent(in) :: f_xr, f_yr, f_zr
      real, dimension(lt2), intent(in) :: f_pr

      integer, intent(in) :: k
      real, dimension(k+1,k), intent(in) :: H

      integer :: i, j, cnt
      real :: beta
      complex*16, dimension(k) :: vals
      complex*16, dimension(k,k) :: vecs
      real, dimension(k) :: residual
      character(len=80) filename

!     -->
      if(nid .eq. 0) then
         write(6, *)
         write(6, *) 'Outposting Krylov subspace to'
      endif

!     --> Outpost the latest Krylov vector.
      call whereyouwant("KRY", k+1)
      ifto = .false.; ifpo = .true.; time = time * k !order in ParaView
      call outpost(f_xr, f_yr, f_zr, f_pr, t, "KRY")

!     --> Compute the eigenvalues and eigenvectors of the current Hessenberg matrix.
      call eig(H(1:k, 1:k), vecs, vals, k)

!     --> Compute the residual (assume H results from classical Arnoldi factorization).
      residual = abs(H(k+1, k) * vecs(k, :))
      cnt = count(residual .lt. eigen_tol)

      if (nid.EQ.0) then
!     --> Outpost the eigenspectrum and residuals of the current Hessenberg matrix.      
         write(filename, '(a, i4.4)') "H", k
         write(6, *) 'Outposting the eigenspectrum of the Hessenberg matrix to : ', filename

         open(67, file=trim(filename), status='unknown', form='formatted')
         write(67, '(3E15.7)') (real(vals(i)), aimag(vals(i)), residual(i), i=1, k)
         close(67)

!     --> Outpost the Hessenberg matrix for restarting purposes (if needed).
         write(filename, '(a, a, i4.4)') 'HES', trim(SESSION), k
         write(6, *) 'Outposting the current Hessenberg matrix to :', filename

         open(67, file=trim(filename), status='unknown', form='formatted')
         write(67, *) ((H(i, j), j=1, k), i=1, k+1)
         close(67)

!     --> Write to logfile the current number of converged eigenvalues.
         write(6, *) 'So far', cnt, ' eigenvalues have converged.'
      endif

      return
      end

      function log_transform(x)
      implicit none

!     -->
      complex :: x, y
      complex :: log_transform

!     -->
      log_transform = log(x)
      if (aimag(x) .eq. 0) log_transform = real(log_transform)
      
      return
      end
