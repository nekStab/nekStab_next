      module LinearOperators

         use lightkrylov, krylov_atol => atol
         use NekVectors

         implicit none
         include 'SIZE'
         include 'TOTAL'
         include 'ADJOINT'

         private

         !------------------------------------------
         !-----     EXPONENTIAL PROPAGATOR     -----
         !------------------------------------------

         type, extends(abstract_linop), public :: exponential_prop
            real :: t
         contains
            private
            procedure, pass(self), public :: matvec => direct_solver
            procedure, pass(self), public :: rmatvec => adjoint_solver
         end type exponential_prop

         !--------------------------------------
         !-----     RESOLVENT OPERATOR     -----
         !--------------------------------------

         type, extends(abstract_linop), public :: resolvent_op
            real :: t
         contains
            private
            procedure, pass(self), public :: matvec  => direct_map
            procedure, pass(self), public :: rmatvec => adjoint_map
         end type resolvent_op

      contains

         subroutine direct_solver(self, vec_in, vec_out)

            implicit none
            include 'SIZE'
            include 'TOTAL'
            include 'ADJOINT'

            class(exponential_prop), intent(in)  :: self
            class(abstract_vector), intent(in)  :: vec_in
            class(abstract_vector), intent(out) :: vec_out

            logical, save :: init
            data             init /.false./

            select type (vec_in)
             type is (real_nek_vector)
               select type (vec_out)
                type is (real_nek_vector)

                  if (findiff_order>1) then
                     call setup_FD
                     call setupLinearSolver('DIRECT FD', init)
                     call nekstab_solver('DIRECT FD', init, vec_in, vec_out)
                  else
                     call setupLinearSolver('DIRECT', init)
                     call nekstab_solver('DIRECT', init, vec_in, vec_out)
                  end if

               end select
            end select

         contains

            subroutine setup_FD
               implicit none
               include 'SIZE'
               include 'TOTAL'

               ! For 2nd-order finite difference
               if (findiff_order == 2) then ! coefs and amps allocated with 4 elements
                  ampls = [1.0D0, -1.0D0, 0.0D0, 0.0D0]
                  coefs = [1.0D0, -1.0D0, 0.0D0, 0.0D0] / 2.0D0
                  ! For 4th-order finite difference
               else if (findiff_order == 4) then
                  ampls = [1.0D0, -1.0D0, 2.0D0, -2.0D0]
                  coefs = [8.0D0, -8.0D0, -1.0D0, 1.0D0] / 12.0D0
               endif

            end subroutine setup_FD

         end subroutine direct_solver

         subroutine adjoint_solver(self, vec_in, vec_out)

            implicit none
            include 'SIZE'
            include 'TOTAL'
            include 'ADJOINT'

            class(exponential_prop), intent(in)  :: self
            class(abstract_vector), intent(in)  :: vec_in
            class(abstract_vector), intent(out) :: vec_out

            logical, save :: init
            data             init /.false./

            select type (vec_in)
             type is (real_nek_vector)
               select type (vec_out)
                type is (real_nek_vector)

                  call setupLinearSolver('ADJOINT', init)
                  call nekstab_solver('ADJOINT', init, vec_in, vec_out)

               end select
            end select

         end subroutine adjoint_solver

         subroutine setupLinearSolver(solver_type,init)

            implicit none
            include 'SIZE'
            include 'TOTAL'
            include 'ADJOINT'

            character(len=*), intent(in) :: solver_type
            logical, intent(inout) :: init

            ! turn on linearized solver (default)
            ! integrates vxp, vyp, vzp, prp, tp
            ifpert = .true.
            if (solver_type == 'DIRECT FD')  ifpert = .false. ! no need for linearized solver

            ! turn off adjoint solver (default)
            ifadj = .false.
            if (solver_type == 'ADJOINT') ifadj = .true.

            ! turn off nonlinear solver (default)
            ! integrates vx, vy, vz, pr, t
            ifbase = .false.

            if (isFloquetDirect .or. isFloquetAdjoint) then
               ifbase = .true.

            else if (isFloquetDirectAdjoint) then

               if (solver_type == 'ADJOINT') then
                  init = .true.
               end if

               ifbase = .true.

            else if (isNewtonPO .or. isNewtonPO_T) then

               init = .true.
               ifbase = .true.

            end if

            if (isFloquetDirect.or.isFloquetAdjoint.or.isFloquetDirectAdjoint) then
               if(nid.eq.0)write(*,*)'Floquet mode'

               ! turn off nonlinear solver if criteria met
               if (ifstorebase .and. init) ifbase = .false.

               ! allocate memory for periodic base flow
               if (ifstorebase .and. ifbase .and. .not. init) then

                  if(nid .eq. 0) write(*,*) 'ALLOCATING ORBIT WITH NSTEPS:', nsteps

                  allocate(uor(lv, nsteps), vor(lv, nsteps))

                  if (if3d) then
                     allocate(wor(lv, nsteps))
                  else
                     allocate(wor(1, 1))
                  end if

                  if (ifto .or. ldimt > 1) allocate(tor(lt, nsteps, ldimt))

               end if

               ! activate Floquet for intracycle transient growth
               ! in the direct/adjoint mode at this point
               ! the nonlinear solution is already stored
               if (solver_type == 'ADJOINT' .and. isFloquetDirectAdjoint) init=.true.

            end if ! Floquet mode

            ! just broadcast the flags !
            call bcast(ifpert, lsize)
            call bcast(ifadj, lsize)
            call bcast(ifbase, lsize)
            call bcast(ifstorebase, lsize)

         end subroutine

         subroutine nekstab_solver(solver_type, init, in, out)
            use NekVectors
            implicit none
            include 'SIZE'
            include 'TOTAL'
            include 'ADJOINT'

            character(len=*), intent(in) :: solver_type
            logical, intent(inout) :: init
            type(real_nek_vector), intent(in) :: in
            type(real_nek_vector), intent(out) :: out

            type(real_nek_vector) :: work, pert
            real :: epsilon0, work_norm
            integer :: i, m
            integer, save :: krylov_counter, matvec_counter
            data             krylov_counter /1/
            data             matvec_counter /0/

            nt = nx1*ny1*nz1*nelt

            ! --> Setup the Nek parameters for the finite-differences approximation.
            if (solver_type == 'DIRECT FD') then

               call nopcopy(work%vx, work%vy, work%vz, work%pr, work%t, ubase, vbase, wbase, pbase, tbase)

               work_norm = work%norm()
               epsilon0 = epsilon * work_norm
               ampls = ampls * epsilon0

               if(nid.eq.0) then
                  write(*,*) 'epsilon, work_norm=',epsilon,work_norm
                  write(*,*) 'epsilon0=',epsilon0
                  write(*,*) 'ampls(:)=',ampls(:)
                  write(*,*) 'coefs(:)=',coefs(:)
               endif

            else

               ! --> Pass the initial condition for the perturbation.
               call nopcopy(vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,:,1), in%vx, in%vy, in%vz, in%pr, in%t)

            end if

            do i = 1, findiff_order

               if (solver_type == 'DIRECT FD') then

                  ! --> Scale the perturbation.
                  call nopcopy(pert%vx, pert%vy, pert%vz, pert%pr, pert%t, in%vx, in%vy, in%vz, in%pr, in%t)
                  call nopcmult(pert%vx, pert%vy, pert%vz, pert%pr, pert%t, ampls(i))

                  ! --> Initial condition for the each evaluation.
                  call nopcopy(vx, vy, vz, pr, t, ubase, vbase, wbase, pbase, tbase)
                  call nopadd2(vx, vy, vz, pr, t, pert%vx, pert%vy, pert%vz, pert%pr, pert%t)
                  if(ifbf2d .and. if3d) call rzero(vz, nv)

               endif

               ! --> Integrate in time.
               time = 0.0D0
               do istep = 1, nsteps

                  ! --> Output current info to logfile.
                  if(nid .eq. 0) then
                     write(*,*)''
                     if(findiff_order>1)then
                        write(6,"(' ', A17,' (',I1,'/',I1,')' ':', I6, '/', I6, ' from', I5, '/', I5, ' [', I1, ']')") solver_type,i,findiff_order,istep, nsteps, krylov_counter, k_dim, schur_cnt
                     else
                        write(6,"(' ', A17, ':', I6, '/', I6, ' from', I5, '/', I5, ' [', I1, ']')") solver_type, istep, nsteps, krylov_counter, k_dim, schur_cnt
                     end if
                  end if

                  ! --> Integrate in time.
                  omega_t = 0.0D0
                  if (INDEX(solver_type, 'RESOLVENT') > 0) then
                     omega_t = (2.0D0 * pi / fintim) * time
                     !write(*,*) 'omega_t=',omega_t
                  endif
                  call nekstab_usrchk()
                  call nek_advance()
                  matvec_counter = matvec_counter + 1

                  if (ifstorebase.and.ifbase.and..not.init)then !storing first time

                     if(nid.eq.0)write(*,*)'storing first series:',istep,'/',nsteps
                     call opcopy(uor(:,istep),vor(:,istep),wor(:,istep),vx,vy,vz)
                     if (ifto) call copy(tor(:,istep,1),t(:,:,:,:,1),nt)
                     if (ldimt.gt.1) then
                        do m = 2,ldimt
                           if(ifpsco(m-1)) call copy(tor(:,istep,m),t(:,:,:,:,m),nt)
                        end do
                     end if

                  elseif(ifstorebase.and.init.and..not.ifbase)then !just moving in memory

                     if(nid.eq.0)write(*,*)'using stored baseflow'
                     call opcopy(vx,vy,vz,uor(:,istep),vor(:,istep),wor(:,istep))
                     if (ifto) call copy(t(:,:,:,:,1),tor(:,istep,1),nt)
                     if (ldimt.gt.1) then
                        do m = 2,ldimt
                           if(ifpsco(m-1)) call copy(t(:,:,:,:,m),tor(:,istep,m),nt)
                        end do
                     end if

                  end if

               end do ! istep

               if( ifstorebase .and. .not. init .and. ifbase )then
                  ifbase=.false. ; init=.true.
               end if

               if (solver_type == 'DIRECT FD') then
                  ! --> Copy the solution and compute the approximation of the Fréchet dérivative.
                  call nopcopy(work%vx, work%vy, work%vz, work%pr, work%t, vx, vy, vz, pr, t)
                  call nopcmult(work%vx, work%vy, work%vz, work%pr, work%t, coefs(i))
                  call nopadd2(out%vx, out%vy, out%vz, out%pr, out%t, work%vx, work%vy, work%vz, work%pr, work%t)
               endif

            enddo

            if (solver_type == 'DIRECT FD') then
               ! Rescale the approximate Fréchet derivative with the step size.
               call nopcmult(out%vx, out%vy, out%vz, out%pr, out%t, 1.0D0/epsilon0)
            else
               ! --> Copy the solution.
               call nopcopy(out%vx,out%vy,out%vz,out%pr,out%t, vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,:,1))
            endif

            if (nid.eq.0)then
               write(*,*) 'matvec_counter=',matvec_counter
               !if (krylov_counter==1) then
               !   open(unit=10, file=trim('matvec_counter.dat'), status='replace', position='append', form='formatted')
               !else
               open(unit=10, file=trim('matvec_counter.dat'), status='old', position='append', form='formatted')
               !endif
               write(10, "(2I10)") krylov_counter, matvec_counter
               close(10)
            endif

            krylov_counter = krylov_counter + 1

         end subroutine nekstab_solver

         !-------------------------------------------------------------------
         !-----     TYPE-BOUND PROCEDURE FOR THE RESOLVENT OPERATOR     -----
         !-------------------------------------------------------------------

         subroutine direct_map(self, vec_in, vec_out)

            implicit none
            class(resolvent_op)   , intent(in)  :: self
            class(abstract_vector), intent(in)  :: vec_in
            class(abstract_vector), intent(out) :: vec_out

            logical, save :: init
            data             init /.false./

            call resolvent_solver('DIRECT', init, vec_in, vec_out)

         end subroutine direct_map

         subroutine adjoint_map(self, vec_in, vec_out)

            implicit none
            class(resolvent_op)   , intent(in)  :: self
            class(abstract_vector), intent(in)  :: vec_in
            class(abstract_vector), intent(out) :: vec_out

            logical, save :: init
            data             init /.false./

            call resolvent_solver('ADJOINT', init, vec_in, vec_out)

         end subroutine adjoint_map

         subroutine resolvent_solver(solver_type, init, vec_in, vec_out)
            use NekVectors
            implicit none
            include 'SIZE'
            include 'TOTAL'
            include 'ADJOINT'

            character(len=*), intent(in) :: solver_type
            logical, intent(inout) :: init
            class(abstract_vector), intent(in) :: vec_in
            class(abstract_vector), intent(out) :: vec_out

            !> Exponential propagator.
            class(exponential_prop), allocatable :: A
            class(identity_linop)  , allocatable :: Id
            class(axpby_linop), allocatable :: S
            class(real_nek_vector), allocatable :: b
            type(gmres_opts)          :: opts
            integer                   :: info
            real :: tol
            character(len=17) :: case

            logical, save :: init2
            data             init2 /.false./

            tol = max(param(21),param(22))
            Id = identity_linop()
            A = exponential_prop(fintim)

            select type (vec_in)
             type is (real_nek_vector)

               !> Compute the right-hand side vector for gmres.
               allocate(b, source=vec_in) ; call b%zero()

               !> Sets the forcing spacial support.
               call nopcopy(fReu,fRev,fRew,fRep,fRet, real(vec_in%vx),real(vec_in%vy),real(vec_in%vz),real(vec_in%pr),real(vec_in%t))
               call nopcopy(fImu,fImv,fImw,fImp,fImt, vec_in%vx,vec_in%vy,vec_in%vz,vec_in%pr,vec_in%t)

            end select

            !> Compute int_{0}^time exp( (t-s)*A ) * f(s) ds
            call setupLinearSolver(solver_type, init2)
            case = 'RESOLVENT ' // trim(solver_type)

            select type (vec_in_casted => vec_in)
             type is (real_nek_vector)
               call nekstab_solver(case, init2, vec_in_casted, b)
            end select

            !> GMRES solve to compute the post-transient response.
            opts = gmres_opts(verbose=.true., atol=tol, rtol=tol)
            if (solver_type == 'DIRECT') then
               S = axpby_linop(Id, A, 1.0D0, -1.0D0)
               call gmres(S, b, vec_out, info, options=opts)

            elseif(solver_type == 'ADJOINT') then

               S = axpby_linop(Id, A, 1.0D0, -1.0D0, .false., .true.)
               call gmres(S, b, vec_out, info, options=opts, transpose=.false.)

            endif

         end subroutine resolvent_solver

      end module LinearOperators