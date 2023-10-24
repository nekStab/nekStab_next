      subroutine linear_stability_analysis
         use NekVectors
         use LinearOperators
         use LightKrylov, krylov_atol => atol
         implicit none
         include 'SIZE'
         include 'TOTAL'
         
         class(exponential_prop), allocatable :: A
         class(real_nek_vector), allocatable :: X(:),seed
         complex :: eigvecs(k_dim, k_dim), eigvals(k_dim)
         real :: residuals(k_dim)
         integer :: info
         logical :: transpose

         allocate(X(1:k_dim + 1),seed)

         call prepare_base_flow
         call prepare_linearized_solver
         A = exponential_prop(fintim)
         call prepare_perturbation_seed(seed)
         call nopcopy(X(1)%vx,X(1)%vy,X(1)%vz,X(1)%pr,X(1)%t, seed%vx,seed%vy,seed%vz,seed%pr,seed%t)
         
         if      (isDirect .or. isFloquetDirect) then
            evop = 'd'
            transpose = .false.
         elseif (isAdjoint .or. isFloquetAdjoint) then
            evop = 'a'
            transpose = .true.
         end if
         call eigs(A, X, eigvecs, eigvals, residuals, info)!, nev, tolerance, verbosity, transpose)
        
         call outpost_eigenspectrum(eigvals, residuals, 'Spectre_H'//trim(evop)//'.dat')
        
         eigvals = log(eigvals)/A%t
        
         call outpost_eigenspectrum(eigvals, residuals, 'Spectre_NS'//trim(evop)//'.dat')
         !call outpost_eigenvectors()
        
      end subroutine linear_stability_analysis

      subroutine prepare_base_flow
         implicit none
         include 'SIZE'
         include 'TOTAL'

         character(len=30) :: filename      
         time = 0.0d0 ! might be overwritten by load_fld
      
         if (ifldbf) then !skip loading if single run
         
            if (nid .eq. 0) write (*, *) ' Loading base flow from disk:'
            write (filename, '(a,a,a)') 'BF_', trim(SESSION), '0.f00001'
            call load_fld(filename)

            !time might containt the exact period of the PO (if Floquet)
         
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
      
         type(real_nek_vector), intent(inout) :: seed
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
           
         else ! seed is prescribed in the base flow >????

            ! might need to call this first so we can load the BF after and subscribe the seed
            ! not sure if this is needed
            call nopcopy(seed%vx,seed%vy,seed%vz,seed%pr,seed%t, ubase,vbase,wbase,pbase,tbase)
         
         end if
        
         !call nopcopy(X(1)%vx,X(1)%vy,X(1)%vz,X(1)%pr,X(1)%t, wrk%vx,wrk%vy,wrk%vz,wrk%pr,wrk%t)
        
         return ! seed returns to X(1)
      end subroutine prepare_perturbation_seed

      subroutine outpost_eigenspectrum(eigvals, residuals, filename)
         implicit none
         include 'SIZE'
         include 'TOTAL'

         complex, intent(in) :: eigvals(k_dim)
         real, intent(in) :: residuals(k_dim)
         character(len=*), intent(in) :: filename
         integer :: i

         if (nid.eq.0)then
         open(unit=10, file=trim(filename), status='replace', form='formatted')
         do i = 1, size(eigvals)
            write(10, "(3E15.7)") real(eigvals(i)), aimag(eigvals(i)), residuals(i)
         end do
         close(10)
         endif

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