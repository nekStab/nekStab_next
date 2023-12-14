      !-----------------------------------------------------------------------
      subroutine whereyouwant(resnam, posfil) !file numbering suffix counter
         implicit none
         character(len=3) :: resnam
         integer :: posfil, i_find_prefix, nopen
         common/RES_WANT/nopen(99, 2)
      
      ! Update the corresponding entry in nopen
         nopen(i_find_prefix(resnam, 99), 1) = posfil - 1
      end subroutine whereyouwant
      !-----------------------------------------------------------------------
      subroutine load_files(Q, mstart, kd, fname)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
      
      ! ----- Krylov basis V for the projection M*V = V*H -----
         type(krylov_vector), dimension(kd) :: Q
      
         integer :: kd, mstart, i, j, m, ndgts
         character(len=3) :: fname
         character(len=7) :: tl
         character(len=20) :: fmt
         character(len=60) :: filename
      
         n2 = nx2*ny2*nz2*nelv
         nt = nx1*ny1*nz1*nelt
      
      ! ----- Upload the snapshots -----
      
         do i = 1, mstart
            ndgts = 1  ! Determine the number of digits in i
            if (i >= 10) ndgts = 2
            if (i >= 100) ndgts = 3
            if (i >= 1000) ndgts = 4
            if (i >= 10000) ndgts = 5
      
      ! Set the format for writing i
            write (tl, '(i5)') ndgts
            fmt = '('//'i'//trim(tl)//'.'//trim(tl)//')'
      
      ! Write i using the determined format
            j = i
            write (tl, fmt) j
      
      ! Set the filename based on the number of digits in i
            select case (ndgts)
            case (1)
               filename = trim(fname)//trim(SESSION)//'0.f0000'//trim(tl)
            case (2)
               filename = trim(fname)//trim(SESSION)//'0.f000'//trim(tl)
            case (3)
               filename = trim(fname)//trim(SESSION)//'0.f00'//trim(tl)
            case (4)
               filename = trim(fname)//trim(SESSION)//'0.f0'//trim(tl)
            case (5)
               filename = trim(fname)//trim(SESSION)//'0.f'//trim(tl)
            end select
      
            call load_fld(filename)
            call opcopy(Q(i)%vx, Q(i)%vy, Q(i)%vz, vx, vy, vz)
            if (ifpo) call copy(Q(i)%pr, pr, n2)
            if (ifto) call copy(Q(i)%t(:, 1), t(:, :, :, :, 1), nt)
            if (ldimt > 1) then
            do m = 2, ldimt
               if (ifpsco(m - 1)) call copy(Q(i)%t(:, m), t(:, :, :, :, m), nt)
            end do
            end if
         end do
      
         return
      end subroutine load_files
