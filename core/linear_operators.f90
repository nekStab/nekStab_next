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

        !write(*,*)'direct solver, init=',init
        call setupLinearSolver('DIRECT', init)
        call nekstab_solver('DIRECT', init, vec_in, vec_out)

        end select
        end select
        
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

        character(len=8), intent(in) :: solver_type
        logical, intent(inout) :: init

        write(*,*)'setup solver, init=',init

        if(init)then
            if(nid.eq.0)write(*,*)'initializing linear solver'
            if(nid.eq.0)write(*,*)'solver type:',solver_type
        end if

        ! turn on linearized solver (default)
        ! integrates vxp, vyp, vzp, prp, tp
        ifpert = .true.
        call bcast(ifpert, lsize)

        ! turn on adjoint solver if needed
        ifadj = .false.
        if (solver_type == 'ADJOINT') ifadj = .true.
        call bcast(ifadj, lsize)

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

        end subroutine

        subroutine nekstab_solver(solver_type, init, in, out)
        use NekVectors
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'ADJOINT'

        character(len=7), intent(in) :: solver_type
        logical, intent(inout) :: init
        type(real_nek_vector), intent(in) :: in
        type(real_nek_vector), intent(out) :: out
        
        integer, save :: mstep_counter 
        data             mstep_counter /1/

        integer m
        nt = nx1*ny1*nz1*nelt

        ! --> Pass the initial condition for the perturbation.
        call nopcopy(vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,:,1), in%vx, in%vy, in%vz, in%pr, in%t)

        time = 0.0D0
        do istep = 1, nsteps

        ! --> Output current info to logfile.
        if(nid .eq. 0) then
            write(*,*)''
            write(6,"(' ', A7, ':', I6, '/', I6, ' from', I6, '/', I6, ' (', I3, ')')") trim(solver_type), istep, nsteps, mstep_counter, k_dim, schur_cnt
        end if

        ! --> Integrate in time.
        call nekstab_usrchk()
        call nek_advance()

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

        ! --> Copy the solution.
        call nopcopy(out%vx,out%vy,out%vz,out%pr,out%t, vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,:,1))
        mstep_counter = mstep_counter + 1

        end subroutine nekstab_solver

        !-------------------------------------------------------------------
        !-----     TYPE-BOUND PROCEDURE FOR THE RESOLVENT OPERATOR     -----
        !-------------------------------------------------------------------

        subroutine direct_map(self, vec_in, vec_out)
            class(resolvent_op)   , intent(in)  :: self
            class(abstract_vector), intent(in)  :: vec_in
            class(abstract_vector), intent(out) :: vec_out

            !> Exponential propagator.
            class(exponential_prop), allocatable :: A
            class(identity_linop)  , allocatable :: Id
            class(axpby_linop), allocatable :: S
            class(real_nek_vector), allocatable :: b
            type(real_nek_vector) :: forcing
            type(gmres_opts)          :: opts
            integer                   :: info
            real, parameter :: pi = 4.0D0 * atan(1.0D0)
            real :: omega, tol

            logical, save :: init
            data             init /.false./

            tol = max(param(21),param(22)) 
            Id = identity_linop() ; A = exponential_prop(fintim)
            omega = 2.0D0 * pi / A%t

            select type (vec_in)
            type is (real_nek_vector)
            select type (vec_out)
            type is (real_nek_vector)
        
                !> Compute int_{0}^t exp( (t-s)*A ) * f(s) ds
                call setupLinearSolver('DIRECT', init)
                call nekstab_solver('DIRECT', init, vec_in, vec_out)
        
                !> Compute the right-hand side vector for gmres.
                allocate(b, source=vec_in) ; call b%zero()
                !> Sets the forcing spatial support.
                !forcing(:) = vec_in(:)
            
            end select
            end select
            
            !> GMRES solve to compute the post-transient response.
            opts = gmres_opts(verbose=.false., atol=tol, rtol=tol)
            S = axpby_linop(Id, A, 1.0D0, -1.0D0)
            call gmres(S, b, vec_out, info, options=opts)
            
        contains

            subroutine forced_rhs!(me, t, x, f)
            ! !> Time-integrator.
            ! class(rk_class), intent(inout) :: me
            ! !> Current time.
            ! real  , intent(in)                :: t
            ! !> State vector.
            ! real  , dimension(:), intent(in)  :: x
            ! !> Time-derivative.
            ! real  , dimension(:), intent(out) :: f
            
            ! integer :: i, j, k
            ! real, dimension(nx) :: u, du
            ! real, dimension(nx) :: v, dv
            ! real                :: d2u, d2v, cu, cv

            ! f = 0.0D0
            ! u = x(1:nx)      ; du = f(1:nx)
            ! v = x(nx+1:2*nx) ; dv = f(nx+1:2*nx)
            
            ! !> Add the forcing term.
            ! du(:) = du(:) + forcing(1:nx)*cos(self%omega*t) - forcing(nx+1:2*nx)*sin(self%omega*t)
            ! dv(:) = dv(:) + forcing(1:nx)*sin(self%omega*t) + forcing(nx+1:2*nx)*cos(self%omega*t)
            
            ! !> Copy results to the output array.
            ! f(1:nx) = du ; f(nx+1:2*nx) = dv

            end subroutine forced_rhs
            
        end subroutine direct_map

        subroutine adjoint_map(self, vec_in, vec_out)
            class(resolvent_op)   , intent(in)  :: self
            class(abstract_vector), intent(in)  :: vec_in
            class(abstract_vector), intent(out) :: vec_out

            !> Exponential propagator.
            class(exponential_prop), allocatable :: A
            class(identity_linop)  , allocatable :: Id

            ! call gmres(S, b, vec_out, info, options=opts, transpose=.false.)

        contains

            subroutine forced_adjoint_rhs
                            
            end subroutine forced_adjoint_rhs
            
        end subroutine adjoint_map

        end module LinearOperators