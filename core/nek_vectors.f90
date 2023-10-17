      module nek_vectors

      use LightKrylov
      implicit none
      include 'SIZE'

      interface 
      real function glsc3(a,b,mult,n)
      real, intent(in) :: a(*),b(*),mult(*)
      integer, intent(in) :: n
      end function glsc3
      end interface

      private

      integer, public, parameter :: lv = lx1*ly1*lz1*lelv
      integer, public, parameter :: lt = lx1*ly1*lz1*lelt
      integer, public, parameter :: lp = lx2*ly2*lz2*lelv ! lp is used for norm -> lp
      integer, save, public :: nv,nt,n2 ! np conflits with mpi variable -> n2

      type, extends(abstract_vector), public :: real_nek_vector
            real, dimension(lv) :: vx, vy, vz
            real, dimension(lp) :: pr
            real, dimension(lv,ldimt) :: t
            real :: time
      contains
            private
            procedure, pass(self), public :: zero => real_zero
            procedure, pass(self), public :: dot => real_dot
            procedure, pass(self), public :: scal => real_scal
            procedure, pass(self), public :: axpby => real_axpby
      end type real_nek_vector

      type, extends(abstract_vector), public :: cplx_nek_vector
            complex, dimension(lv) :: vx, vy, vz
            complex, dimension(lp) :: pr
            complex, dimension(lv,ldimt) :: t
            complex :: time
      contains
            private
            procedure, pass(self), public :: zero => cplx_zero
            procedure, pass(self), public :: dot => cplx_dot
            procedure, pass(self), public :: scal => cplx_scal
            procedure, pass(self), public :: axpby => cplx_axpby
      end type cplx_nek_vector

      type(real_nek_vector), save, public :: ic_nwt, fc_nwt
      real,save,allocatable,dimension(:, :),public ::uor,vor,wor
      real,save,allocatable,dimension(:, :, :),public :: tor

      contains
       
            ! -----------------------------------------------------------
            ! -----                                                 -----
            ! -----     DEFINITION OF THE TYPE-BOUND PROCEDURES     -----
            ! -----                                                 -----
            ! -----------------------------------------------------------
            
            ! Real prodcedures definitions.

            subroutine real_zero(self)
            implicit none
            include 'SIZE'      
            include 'TOTAL'
            class(real_nek_vector), intent(inout) :: self
            call noprzero(self%vx,self%vy,self%vz,self%pr,self%t)
            self%time = 0.0D0
            return
            end subroutine real_zero
            
            real function real_dot(self, vec) result(alpha)
            implicit none
            include 'SIZE'      
            include 'TOTAL'

            class(real_nek_vector), intent(in)      :: self
            class(abstract_vector), intent(in) :: vec
            integer m
            nv = nx1*ny1*nz1*nelv
            nt = nx1*ny1*nz1*nelt
            select type(vec)
            type is(real_nek_vector)

            ! --> Kinetic energy.
            alpha = glsc3(self%vx, bm1s, vec%vx, nv)
            alpha = alpha + glsc3(self%vy, bm1s, vec%vy, nv)
            if (if3d) alpha = alpha + glsc3(self%vz, bm1s, vec%vz, nv)
            
            ! --> Potential energy.
            if (ifto) alpha = alpha + glsc3(self%t(:,1), bm1s, vec%t(:,1), nt)
            if (ldimt.gt.1) then
               do m = 2,ldimt
                 if(ifpsco(m-1)) alpha = alpha + glsc3(self%t(:,m), bm1s, vec%t(:,m), nt)
               enddo
            endif
            if ( uparam(1) .eq. 2.1 ) then
                  alpha = alpha + self%time * vec%time
            end if
            if ( isnan(alpha) ) then 
                  if (nid.eq.0) write(6,*) 'NaN detected in dot product'
                  call nek_end
            end if
            end select
            return
            end function real_dot
           
            subroutine real_scal(self, alpha)
            implicit none
            include 'SIZE'      
            include 'TOTAL'
            class(real_nek_vector), intent(inout) :: self
            real, intent(in) :: alpha
            call nopcmult(self%vx,self%vy,self%vz,self%pr,self%t, alpha)
            self%time = self%time * alpha
            return
            end subroutine real_scal
            
            subroutine real_axpby(self, alpha, vec, beta)
            implicit none
            include 'SIZE'      
            include 'TOTAL'
            class(real_nek_vector), intent(inout) :: self
            class(abstract_vector), intent(in) :: vec
            real , intent(in) :: alpha, beta
            select type(vec)
            type is(real_nek_vector)
            call nopaxpby(self%vx,self%vy,self%vz,self%pr,self%t,alpha,vec%vx,vec%vy,vec%vz,vec%pr,vec%t,beta)
            end select
            return
            end subroutine real_axpby
      
            ! Complex prodcedures definitions.

            subroutine cplx_zero(self)
            implicit none
            include 'SIZE'      
            include 'TOTAL'
            class(cplx_nek_vector), intent(inout) :: self
            call noprzero(self%vx,self%vy,self%vz,self%pr,self%t)
            self%time = 0.0D0
            return
            end subroutine cplx_zero
            
            real function cplx_dot(self, vec) result(alpha)
            implicit none
            include 'SIZE'      
            include 'TOTAL'
            class(cplx_nek_vector), intent(in) :: self
            class(abstract_vector), intent(in) :: vec
            integer m
            nv = nx1*ny1*nz1*nelv
            nt = nx1*ny1*nz1*nelt
            select type(vec)
            type is(cplx_nek_vector)
            ! bm1s is real, glsc3 is real
            end select
            return
            end function cplx_dot
            
            subroutine cplx_scal(self, alpha)
            implicit none
            include 'SIZE'      
            include 'TOTAL'
            class(cplx_nek_vector), intent(inout) :: self
            real, intent(in) :: alpha
            call nopcmult(self%vx,self%vy,self%vz,self%pr,self%t, alpha)
            self%time = self%time * alpha
            return
            end subroutine cplx_scal
            
            subroutine cplx_axpby(self, alpha, vec, beta)
            implicit none
            include 'SIZE'      
            include 'TOTAL'
            class(cplx_nek_vector), intent(inout) :: self
            class(abstract_vector), intent(in) :: vec
            real , intent(in) :: alpha, beta
            select type(vec)
            type is(cplx_nek_vector)
            call nopaxpby(self%vx,self%vy,self%vz,self%pr,self%t,alpha,vec%vx,vec%vy,vec%vz,vec%pr,vec%t,beta)
            end select
            return
            end subroutine cplx_axpby

      end module nek_vectors

c-----------------------------------------------------------------------

      subroutine noprzero(a1,a2,a3,a4,a5)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer n,k
      real, intent(inout) :: a1(1),a2(1),a3(1),a4(1),a5(lx1*ly1*lz1*lelt,1)
      n=nx1*ny1*nz1*nelv
      call rzero(a1,n)
      call rzero(a2,n)
      if (if3D) call rzero(a3,n)
      if (ifpo) call rzero(a4,nx2*ny2*nz2*nelv)
      if (ifto) call rzero(a5(1,1),lx1*ly1*lz1*nelfld(2))
      if (ldimt.gt.1) then 
         do k=1,npscal
               if(ifpsco(k)) call rzero(a5(1,k+1),lx1*ly1*lz1*nelfld(k+2))
         enddo
      endif
      return
      end subroutine noprzero
c-----------------------------------------------------------------------
      subroutine nopcmult(a1,a2,a3,a4,a5,c)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer n,k
      real, intent(inout) :: a1(1),a2(1),a3(1),a4(1),a5(lx1*ly1*lz1*lelt,1)
      real, intent(in) :: c
      n=nx1*ny1*nz1*nelv
      call cmult(a1,c,n)
      call cmult(a2,c,n)
      if (if3D) call cmult(a3,c,n)
      if (ifpo) call cmult(a4,c,nx2*ny2*nz2*nelv)
      if (ifto) call cmult(a5(1,1),c,lx1*ly1*lz1*nelfld(2))
      if (ldimt.gt.1) then 
         do k=1,npscal
               if(ifpsco(k)) call cmult(a5(1,k+1),c,lx1*ly1*lz1*nelfld(k+2))
         enddo
      endif
      return
      end subroutine nopcmult
c-----------------------------------------------------------------------
      subroutine axpby(x, alpha, y, beta, n)
            real x(1),y(1),alpha,beta
            do i=1,n
               x(i) = x(i)*alpha + y(i)*beta
            enddo
            return
      end subroutine axpby
      subroutine nopaxpby(a1,a2,a3,a4,a5,alpha, b1,b2,b3,b4,b5,beta)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer n,k
      real, intent(inout) :: a1(1),a2(1),a3(1),a4(1),a5(lx1*ly1*lz1*lelt,1)
      real, intent(in) :: b1(1),b2(1),b3(1),b4(1),b5(lx1*ly1*lz1*lelt,1)
      real, intent(in) :: alpha, beta
      n=nx1*ny1*nz1*nelv
      call axpby(a1,alpha,b1,beta,n)
      call axpby(a2,alpha,b2,beta,n)
      if (if3D) call axpby(a3,alpha,b3,beta,n)
      if (ifpo) call axpby(a4,alpha,b4,beta,nx2*ny2*nz2*nelv)
      if (ifto) call axpby(a5(1,1),alpha,b5(1,1),beta,lx1*ly1*lz1*nelfld(2))
      if (ldimt.gt.1) then 
         do k=1,npscal
               if(ifpsco(k)) call axpby(a5(1,k+1),alpha,b5(1,k+1),beta,lx1*ly1*lz1*nelfld(k+2))
         enddo
      endif
      return
      end subroutine nopaxpby
c-----------------------------------------------------------------------


      ! EXTRA DEFINITONS USED IN NKESTAB ! 



c-----------------------------------------------------------------------
      subroutine nopcopy(a1,a2,a3,a4,a5, b1,b2,b3,b4,b5)
         implicit none
         include 'SIZE'
         include 'TOTAL'
         integer :: n,k
         real, intent(inout) :: a1(1),a2(1),a3(1),a4(1),a5(lx1*ly1*lz1*lelt,1)
         real, intent(in) :: b1(1),b2(1),b3(1),b4(1),b5(lx1*ly1*lz1*lelt,1)
         n = nx1*ny1*nz1*nelv
         call copy(a1,b1,n)
         call copy(a2,b2,n)
         if (if3D) call copy(a3,b3,n)
         if (ifpo) call copy(a4,b4,nx2*ny2*nz2*nelv)
         if (ifto) call copy(a5(1,1),b5(1,1),lx1*ly1*lz1*nelfld(2))
         if (ldimt.gt.1) then 
            do k=1,npscal
               if(ifpsco(k)) call copy(a5(1,k+1),b5(1,k+1),lx1*ly1*lz1*nelfld(k+2))
            enddo
         endif
         return
      end subroutine nopcopy       
c-----------------------------------------------------------------------
      subroutine nopsub2(a1,a2,a3,a4,a5, b1,b2,b3,b4,b5)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer n,k
      real, intent(inout) :: a1(1),a2(1),a3(1),a4(1),a5(lx1*ly1*lz1*lelt,1)
      real, intent(in) :: b1(1),b2(1),b3(1),b4(1),b5(lx1*ly1*lz1*lelt,1)
      n=nx1*ny1*nz1*nelv
      call sub2(a1,b1,n)
      call sub2(a2,b2,n)
      if (if3D) call sub2(a3,b3,n)
      if (ifpo) call sub2(a4,b4,nx2*ny2*nz2*nelv)
      if (ifto) call sub2(a5(1,1),b5(1,1),lx1*ly1*lz1*nelfld(2))
      if (ldimt.gt.1) then 
         do k=1,npscal
               if(ifpsco(k)) call sub2(a5(1,k+1),b5(1,k+1),lx1*ly1*lz1*nelfld(k+2))
         enddo
      endif
      return
      end subroutine nopsub2
c-----------------------------------------------------------------------
      subroutine nopsub3(c1,c2,c3,c4,c5, a1,a2,a3,a4,a5, b1,b2,b3,b4,b5)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer n,k
      real, intent(inout) :: c1(1),c2(1),c3(1),c4(1),c5(lx1*ly1*lz1*lelt,1)
      real, intent(in) :: a1(1),a2(1),a3(1),a4(1),a5(lx1*ly1*lz1*lelt,1)
      real, intent(in) :: b1(1),b2(1),b3(1),b4(1),b5(lx1*ly1*lz1*lelt,1)
      n=nx1*ny1*nz1*nelv
      call sub3(c1,a1,b1,n)
      call sub3(c2,a2,b2,n)
      if (if3D) call sub3(c3,a3,b3,n)
      if (ifpo) call sub3(c4,a4,b4,nx2*ny2*nz2*nelv)
      if (ifto) call sub3(c5(1,1),a5(1,1),b5(1,1),lx1*ly1*lz1*nelfld(2))
      if (ldimt.gt.1) then 
         do k=1,npscal
               if(ifpsco(k)) call sub3(c5(1,k+1),a5(1,k+1),b5(1,k+1),lx1*ly1*lz1*nelfld(k+2))
         enddo
      endif
      return
      end subroutine nopsub3
c-----------------------------------------------------------------------
      subroutine nopadd2(a1,a2,a3,a4,a5, b1,b2,b3,b4,b5)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer n,k
      real, intent(inout) :: a1(1),a2(1),a3(1),a4(1),a5(lx1*ly1*lz1*lelt,1)
      real, intent(in) :: b1(1),b2(1),b3(1),b4(1),b5(lx1*ly1*lz1*lelt,1)
      n=nx1*ny1*nz1*nelv
      call add2(a1,b1,n)
      call add2(a2,b2,n)
      if (if3D) call add2(a3,b3,n)
      if (ifpo) call add2(a4,b4,nx2*ny2*nz2*nelv)
      if (ifto) call add2(a5(1,1),b5(1,1),lx1*ly1*lz1*nelfld(2))
      if (ldimt.gt.1) then 
         do k=1,npscal
               if(ifpsco(k)) call add2(a5(1,k+1),b5(1,k+1),lx1*ly1*lz1*nelfld(k+2))
         enddo
      endif
      return
      end subroutine nopadd2
c-----------------------------------------------------------------------
      subroutine opadd3(a1,a2,a3,b1,b2,b3,c1,c2,c3)
      implicit none
      include 'SIZE'
      integer n
      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1),c1(1),c2(1),c3(1)
      n=nx1*ny1*nz1*nelv
      call add3(a1,b1,c1,n)
      call add3(a2,b2,c2,n)
      if (ndim.eq.3) call add3(a3,b3,c3,n)
      return
      end subroutine opadd3
c-----------------------------------------------------------------------
      subroutine opaddcol3(a1,a2,a3,b1,b2,b3,c1,c2,c3)
      implicit none
      include 'SIZE'
      integer n
      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1),c1(1),c2(1),c3(1)
      n=nx1*ny1*nz1*nelv
      call addcol3(a1,b1,c1,n)
      call addcol3(a2,b2,c2,n)
      if (ndim.eq.3) call addcol3(a3,b3,c3,n)
      return
      end subroutine opaddcol3