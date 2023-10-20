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
        procedure, pass(self), public :: matvec => direct_linearized_map
        procedure, pass(self), public :: rmatvec => adjoint_linearized_map
        end type exponential_prop

        contains

        ! Subroutine for setting up the linear solver
        subroutine setupLinearSolver(solver_type,init)
        
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'ADJOINT'

        character(len=7), intent(in) :: solver_type
        logical, intent(inout) :: init

        ! turn on linearized solver
        ifpert = .true.
        call bcast(ifpert, lsize)

        ! turn on adjoint solver if needed
        ifadj = .false.
        if (solver_type == 'adjoint') ifadj = .true.
        call bcast(ifadj, lsize)

        ! turn off nonlinear solver
        ifbase = .false.

        if (uparam(01) == 3.11 .or. uparam(01) == 3.21) then
            ifbase = .true.
            
        else if (uparam(01) == 3.31) then

            if (solver_type == 'adjoint') then
                init = .true.
            endif

            ifbase = .true.
            
        else if (uparam(01) == 2.1 .or. uparam(01) == 2.2) then
            init = .true.
            ifbase = .true.
            
        endif
        
        ! turn off nonlinear solver if criteria met
        if (ifstorebase .and. init) ifbase = .false.

        ! allocate memory for periodic base flow
        if (ifstorebase .and. ifbase .and. .not. init) then

            if(nid .eq. 0) write(6,*) 'ALLOCATING ORBIT WITH NSTEPS:', nsteps
            
            allocate(uor(lv, nsteps), vor(lv, nsteps))
            
            if (if3d) then
                allocate(wor(lv, nsteps))
            else
                allocate(wor(1, 1))
            endif
            
            if (ifto .or. ldimt > 1) allocate(tor(lt, nsteps, ldimt))
        
        endif

        ! activate Floquet for intracycle transient growth
        ! in the direct/adjoint mode at this point
        ! the nonlinear solution is already stored
        if (solver_type == 'adjoint' .and. uparam(01) .eq. 3.31) init=.true. 

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
        
        integer m
        nt = nx1*ny1*nz1*nelt

        ! --> Pass the initial condition for the perturbation.
        call nopcopy(vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,:,1), in%vx, in%vy, in%vz, in%pr, in%t)

        time = 0.0D0
        do istep = 1, nsteps
        ! --> Output current info to logfile.
        !if(nid .eq. 0) write(6,"(' ", A7, ":',I6,'/',I6,' from',I6,'/',I6,' (',I3,')')") trim(solver_type), istep, nsteps, mstep, k_dim, schur_cnt

        ! --> Integrate in time.
        call nekstab_usrchk()
        call nek_advance()

        if (ifstorebase.and.ifbase.and..not.init)then !storing first time

            if(nid.eq.0)write(6,*)'storing first series:',istep,'/',nsteps
            call opcopy(uor(:,istep),vor(:,istep),wor(:,istep),vx,vy,vz)
            if (ifto) call copy(tor(:,istep,1),t(:,:,:,:,1),nt)
            if (ldimt.gt.1) then
                do m = 2,ldimt
                    if(ifpsco(m-1)) call copy(tor(:,istep,m),t(:,:,:,:,m),nt)
                enddo
            endif

        elseif(ifstorebase.and.init.and..not.ifbase)then !just moving in memory

            if(nid.eq.0)write(6,*)'using stored baseflow'
            call opcopy(vx,vy,vz,uor(:,istep),vor(:,istep),wor(:,istep))
            if (ifto) call copy(t(:,:,:,:,1),tor(:,istep,1),nt)
            if (ldimt.gt.1) then
                do m = 2,ldimt
                    if(ifpsco(m-1)) call copy(t(:,:,:,:,m),tor(:,istep,m),nt)
                enddo
            endif
        
        endif
        
        end do ! istep

        if( ifstorebase .and. .not. init .and. ifbase )then
            ifbase=.false. ; init=.true.
        endif

        ! --> Copy the solution.
        call nopcopy(out%vx,out%vy,out%vz,out%pr,out%t, vxp(:,1),vyp(:,1),vzp(:,1),prp(:,1),tp(:,:,1))

        end subroutine nekstab_solver

        subroutine forward_linearized_map(self, vec_in, vec_out)

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

        call setupLinearSolver('forward', init)
        call nekstab_solver('forward', init, vec_in, vec_out)

        end select
        end select
        return
        end subroutine forward_linearized_map

        subroutine adjoint_linearized_map(self, vec_in, vec_out)

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

        call setupLinearSolver('adjoint', init)
        call nekstab_solver('adjoint', init, vec_in, vec_out)

        end select
        end select
        return
        end subroutine adjoint_linearized_map

        end module LinearOperators