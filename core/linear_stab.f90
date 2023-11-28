      module LinearStab
         use LightKrylov!, krylov_atol => atol
         use LinearOperators
         use NekVectors
         implicit none

         private
         public :: linear_stability_analysis, transient_growth_analysis, resolvent_analysis

      contains

         subroutine linear_stability_analysis
            implicit none
            include 'SIZE'
            include 'TOTAL'

            class(exponential_prop), allocatable :: A
            class(real_nek_vector), allocatable :: X(:)
            complex :: eigvecs(k_dim, k_dim), eigvals(k_dim)
            real :: residuals(k_dim)
            integer :: info, i
            logical :: transpose

            if(.not.isNekStabinit) call nekStab_init
            call prepare_base_flow

            call set_linear_solver
            A = exponential_prop(fintim)

            allocate(X(k_dim + 1))
            eigvecs(:,:) = 0.0D0 ; eigvals(:) = 0.0D0 ; residuals(:) = 0.0D0
            do i = 1, size(X)
               call X(i)%zero()
            end do
            call prepare_seed(X)

            if (isDirect .or. isFloquetDirect) then
               evop = 'd'; transpose = .false.
               if (findiff_order>1) evop = 'f'
            elseif (isAdjoint .or. isFloquetAdjoint) then
               evop = 'a'; transpose = .true.
            end if

            call eigs(A, X, eigvecs, eigvals, residuals, info, verbosity=.false., nev=schur_tgt, tolerance=eigen_tol, transpose=transpose)

            call outpost_eigenspectrum(eigvals, residuals, 'Spectrum_H'//trim(evop)//'.dat')

            eigvals = log(eigvals)/A%t

            call outpost_eigenspectrum(eigvals, residuals, 'Spectrum_NS'//trim(evop)//'.dat')
            call outpost_eigenvectors(X, eigvals, eigvecs, residuals)

            if(nid.eq.0)write(*,*)'Linear stability finished.'
         end subroutine linear_stability_analysis

         subroutine transient_growth_analysis
            implicit none
            include 'SIZE'
            include 'TOTAL'

            class(exponential_prop), allocatable :: A
            type(real_nek_vector), allocatable :: U(:), V(:)
            real :: sigma(k_dim), uvecs(k_dim, k_dim), vvecs(k_dim, k_dim)
            real :: residuals(k_dim), alpha
            integer :: info, i

            if(.not.isNekStabinit) call nekStab_init
            call prepare_base_flow
            call set_linear_solver
            A = exponential_prop(fintim)

            allocate(U(k_dim+1)) ; allocate(V(k_dim+1))
            sigma(:) = 0.0D0 ; uvecs(:,:) = 0.0D0 ; vvecs(:,:) = 0.0D0 ; residuals(:) = 0.0D0
            do i = 1, size(U)
               call U(i)%zero() ; call V(i)%zero()
            end do
            call prepare_seed(U)

            evop = 'p'
            call svds(A, U, V, uvecs, vvecs, sigma, residuals, info, nev=schur_tgt, tolerance=eigen_tol)
            sigma = sigma ** 2
            call outpost_singvals(sigma(:), residuals(:), 'Spectrum_S'//trim(evop)//'.dat')
            call outpost_singvectors(U, V, uvecs, vvecs, sigma, residuals)
            if(nid.eq.0)write(*,*)'Transient growth finished.'

         end subroutine transient_growth_analysis

         subroutine resolvent_analysis
            implicit none
            include 'SIZE'
            include 'TOTAL'

            type(resolvent_op), allocatable :: R
            type(cmplx_nek_vector), allocatable :: U(:), V(:)
            real :: sigma(k_dim), uvecs(k_dim, k_dim), vvecs(k_dim, k_dim)
            real :: residuals(k_dim), alpha
            integer :: info, i

            if(.not.isNekStabinit) call nekStab_init
            call prepare_base_flow
            call set_linear_solver
            R = resolvent_op(fintim)

            allocate(U(k_dim+1)) ; allocate(V(k_dim+1))
            sigma(:) = 0.0D0 ; uvecs(:,:) = 0.0D0 ; vvecs(:,:) = 0.0D0 ; residuals(:) = 0.0D0
            do i = 1, size(U)
               call U(i)%zero() ; call V(i)%zero()
            end do
            call prepare_seed(U%real)
            call prepare_seed(U%imag)
            
            evop = 'r'
            if(nid.eq.0)write(*,*)'Resolvent started: k_dim, schur_target, eigen_tol = ', k_dim, schur_tgt, eigen_tol
            call svds(R, U, V, uvecs, vvecs, sigma, residuals, info, nev=schur_tgt, tolerance=eigen_tol)
            sigma = sigma ** 2
            
            call outpost_singvals(sigma(:), residuals(:), 'Spectrum_S'//trim(evop)//'.dat')
            !call outpost_singvectors(U%real, V%real, uvecs, vvecs, sigma, residuals)          
            if(nid.eq.0)write(*,*)'Resolvent finished.'
         end subroutine resolvent_analysis

         subroutine prepare_base_flow
            implicit none
            include 'SIZE'
            include 'TOTAL'
            character(len=30) :: filename

            if (ifldbf) then
               if (nid .eq. 0) write (*,*) ' Loading base flow from disk:'
               write (filename, '(a,a,a)') 'BF_', trim(SESSION), '0.f00001'
               call load_fld(filename)
               !time might containt the exact period of the PO (if Floquet)
            else
               if (nid .eq. 0) write (*,*) 'Baseflow prescribed by the useric function in the .usr'
            end if

            call nopcopy(ubase, vbase, wbase, pbase, tbase, vx, vy, vz, pr, t)

         end subroutine prepare_base_flow

         subroutine set_linear_solver
            implicit none
            include 'SIZE'
            include 'TOTAL'

            call nekgsync

            ! Update UPO period from base flow file if Floquet mode is activated
            if (istep == 0 .and. (isFloquetDirect .or. isFloquetAdjoint .or. isFloquetTransientGrowth .or. isFloquetResolvent)) then
               param(10) = time
               if (nid == 0) then
                  write(*,*) 'Floquet mode activated. Getting endTime from file: ', param(10)
               end if
            end if

            ! Limit number of perturbation modes to 1
            if (param(31) > 1) then
               write(*,*) 'Multiple perturbation (npert > 1) not supported. Exiting.'
               call nek_end
            end if
            param(31) = 1 ; npert = param(31)

            ! Deactivate OIFS for linearized solver
            if (ifchar) then
               write(*,*) 'OIFS not compatible with linearized solver. Deactivating.'
               ifchar = .false.
               call bcast(ifchar, lsize)
            end if

            ! Limit target CFL for EXTk to 0.5
            if (param(26) > 1.0) then
               param(26) = 0.5d0
               write(*,*) 'Target CFL limited to 0.5'
            end if

            ! Recompute time steps and endTime to match UPO period
            if (param(10) > 0.0d0) then
               call compute_cfl(ctarg, vx, vy, vz, 1.0d0)
               dt = param(26) / ctarg
               nsteps = ceiling(param(10) / dt)
               dt = param(10) / nsteps
               if (nid == 0) then
                  write(*,*) 'New timeStep dt=', dt
                  write(*,*) 'New numSteps nsteps=', nsteps
               end if
               param(12) = dt
               call compute_cfl(ctarg, vx, vy, vz, dt) ! C=sum(ux_i/dx_i)*dt
               if(nid.eq.0)write(*,*)' current CFL and target=',ctarg,param(26)
               lastep = 0 ! subs1.f:279
               fintim = nsteps*dt
            end if

            param(12) = -abs(param(12)) ! force constant time step
            call bcast(param, 200*wdsize) ! broadcast all parameters

         end subroutine set_linear_solver

         subroutine prepare_seed(X)
            implicit none
            include 'SIZE'
            include 'TOTAL'

            class(real_nek_vector), intent(inout) :: X(k_dim + 1)
            type(real_nek_vector) :: seed
            character(len=30) :: filename
            real :: alpha
            integer :: m

            if (ifseed_nois) then ! noise as initial seed

               if(nid.eq.0)write(*,*)'Filling fields with noise...'

               call op_add_noise(seed%vx,seed%vy,seed%vz)
               if (ifto) call add_noise_scal(seed%t(:,1),9.0e4, 3.0e3, 4.0e5)
               if (ldimt.gt.1) then
                  do m = 2, ldimt
                     if(ifpsco(m-1)) call add_noise_scal(seed%t(:,1),9.0e1*m, 3.0e2*m, 4.0e1*m)
                  end do
               end if

            elseif (ifseed_load) then ! loading initial seed (e.g. dRe )

               if (isDirect.or.isFloquetDirect) then
                  write (filename, '(a,a,a)') 'dRe', trim(SESSION), '0.f00001'
               elseif (isAdjoint.or.isFloquetAdjoint) then
                  write (filename, '(a,a,a)') 'aRe', trim(SESSION), '0.f00001'
               end if

               if (nid .eq. 0) write (*,*) 'Load real part of mode 1 as seed: ', filename
               call load_fld(filename)
               call nopcopy(seed%vx,seed%vy,seed%vz,seed%pr,seed%t, vx,vy,vz,pr,t)

            else ! seed is prescribed in the base flow

               ! might need to call this first so we can load the BF after and subscribe the seed
               ! not sure if this is needed
               !call nopcopy(seed%vx,seed%vy,seed%vz,seed%pr,seed%t, ubase,vbase,wbase,pbase,tbase)
               call nopcopy(seed%vx,seed%vy,seed%vz,seed%pr,seed%t, vx,vy,vz,pr,t) ! seed comes from useric 

            end if

            call nopcopy(X(1)%vx,X(1)%vy,X(1)%vz,X(1)%pr,X(1)%t, seed%vx,seed%vy,seed%vz,seed%pr,seed%t)
            alpha = X(1)%norm()
            call X(1)%scal(1.0D0 / alpha)

         end subroutine prepare_seed

         subroutine outpost_eigenspectrum(eigvals, residuals, filename)
            implicit none
            include 'SIZE'

            complex, intent(in) :: eigvals(k_dim)
            real, intent(in) :: residuals(k_dim)
            character(len=*), intent(in) :: filename
            integer :: i

            if (nid.eq.0) then
               open(unit=10, file=trim(filename), status='replace', form='formatted')
               do i = 1, size(eigvals)
                  write(10, "(3E15.7)") real(eigvals(i)), aimag(eigvals(i)), residuals(i)
               end do
               close(10)
            end if

         end subroutine outpost_eigenspectrum

         subroutine outpost_singvals(sigma, residuals, filename)
            implicit none
            include 'SIZE'

            real, intent(in) :: sigma(k_dim)
            real, intent(in) :: residuals(k_dim)
            character(len=*), intent(in) :: filename
            integer :: i

            if (nid.eq.0) then
               open(unit=10, file=trim(filename), status='replace', form='formatted')
               do i = 1, size(sigma)
                  write(10, "(2E15.7)") sigma(i), residuals(i)
               end do
               close(10)
            end if

         end subroutine outpost_singvals

         subroutine outpost_eigenvectors(X, eigvals, eigvecs, residuals)
            implicit none
            include 'SIZE'
            include 'TOTAL'

            class(real_nek_vector), intent(in) :: X(k_dim + 1)
            complex, intent(in) :: eigvals(k_dim)
            complex, intent(in) :: eigvecs(k_dim,k_dim)
            real, intent(in) :: residuals(k_dim)
            class(abstract_vector), allocatable :: nek_vector

            integer :: i
            character(len=3)  :: nRe,nIm,nRv,nIv

            nRe = trim(evop) // 'Re' ! real part
            nIm = trim(evop) // 'Im' ! imaginary part
            nRv = trim(evop) // 'Rv' ! real part of vorticity
            nIv = trim(evop) // 'Iv' ! imaginary part of vorticity

            do i = 1, maxmodes

               if (nid == 0) then
                  write(*,*) "Outposting eigenvector: ", i, "/", maxmodes
                  write(*,*) "  Sigma: ", real(eigvals(i))
                  write(*,*) "  Omega: ", aimag(eigvals(i))
                  write(*,*) "  Residual: ", residuals(i)
              endif
               !     ----- Output the real part -----
               call get_vec(nek_vector, X(1:k_dim), real(eigvecs(:, i)))
               select type(nek_vector)
                type is(real_nek_vector)
                  call nopcopy(vx,vy,vz,pr,t, nek_vector%vx,nek_vector%vy,nek_vector%vz,nek_vector%pr,nek_vector%t)
                  ! need to normalize to unit norm?
                  call outpost2(vx,vy,vz,pr,t, nof, nRe) 
                  call outpost_vort(vx,vy,vz,nRv)
               end select

               !     ----- Output the imaginary part -----
               call get_vec(nek_vector, X(1:k_dim), aimag(eigvecs(:, i)))
               select type(nek_vector)
                type is(real_nek_vector)
                  call nopcopy(vx,vy,vz,pr,t, nek_vector%vx,nek_vector%vy,nek_vector%vz,nek_vector%pr,nek_vector%t)
                  ! need to normalize to unit norm?
                  call outpost2(vx,vy,vz,pr,t, nof, nIm)
                  !call outpost_vort(vx,vy,vz,nIv)
               end select

            end do ! i = 1, maxmodes

         end subroutine outpost_eigenvectors

         subroutine outpost_singvectors(U, V, uvec, vvec, sigma, residuals)
             implicit none
             include 'SIZE'
             include 'TOTAL'

             class(real_nek_vector), intent(in) :: U(k_dim + 1)
             class(real_nek_vector), intent(in) :: V(k_dim + 1)
             real, intent(in) :: uvec(k_dim, k_dim), vvec(k_dim, k_dim), sigma(k_dim)
             real, intent(in) :: residuals(k_dim)
             class(abstract_vector), allocatable :: nek_vector

             integer :: i
             character(len=3)  :: nU,nV

             nU = trim(evop) // 'U_' ! U matrix
             nV = trim(evop) // 'V_' ! V matrix

             do i = 1, maxmodes

                if (nid == 0) then
                   write(*,*) "Outposting singular vector: ", i, "/", maxmodes
                   write(*,*) "  Sigma: ", sigma(i)
                   write(*,*) "  Residual: ", residuals(i)
               endif
                !     ----- Output the U matrix -----
                call get_vec(nek_vector, U(1:k_dim), uvec(:, i))
                select type(nek_vector)
                 type is(real_nek_vector)
                   call nopcopy(vx,vy,vz,pr,t, nek_vector%vx,nek_vector%vy,nek_vector%vz,nek_vector%pr,nek_vector%t)
                   ! need to normalize to unit norm?
                   call outpost2(vx,vy,vz,pr,t, nof, nU) 
                end select

                !     ----- Output the V matrix -----
                call get_vec(nek_vector, V(1:k_dim), vvec(:, i))
                select type(nek_vector)
                 type is(real_nek_vector)
                   call nopcopy(vx,vy,vz,pr,t, nek_vector%vx,nek_vector%vy,nek_vector%vz,nek_vector%pr,nek_vector%t)
                   ! need to normalize to unit norm?
                   call outpost2(vx,vy,vz,pr,t, nof, nV)
                end select

             end do ! i = 1, maxmodes

         end subroutine outpost_singvectors

      end module LinearStab