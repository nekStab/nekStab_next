c-----------------------------------------------------------------------
!     function to check the velocity value in a node
!     posiz   = position in the velocity vector of the good grid point
!     procmin = processor that contains the good grid point
      subroutine pointcheck(posiz,procmin) 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer m,n,procmin,posiz
      real chk(lx1*ly1*lz1*lelt),chkmin,glmin
      n = nx1*ny1*nz1*nelt
      procmin = 0
      if(nid.eq.0)write(6,*)'Evaluating probe at:',xck,yck,zck
      do m = 1,n
         chk(m)=(xm1(m,1,1,1)-xck)**2+  (ym1(m,1,1,1)-yck)**2
         if(if3D)      chk(m)= chk(m) + (zm1(m,1,1,1)-zck)**2
      enddo
      chkmin = glmin(chk,n)
      do m = 1,n
         if (chkmin.eq.chk(m)) then
            procmin = 1
            posiz = m
            print *, 'Point found: ',m
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine quicksort2(n, arr, idx)
      implicit none
      integer n, m, nstack
      real arr(n)
      integer idx(n)
      parameter (m=7, nstack=50)
      integer i, ir, j, jstack, k, l, istack(nstack)
      real a, temp
      integer b, temp_int

      jstack = 0
      l = 1
      ir = n

 1    if (ir-l .lt. m) then
         do 12 j = l+1, ir
            a = arr(j)
            b = idx(j)
            do 11 i = j-1, l, -1
               if (arr(i) .le. a) goto 2
               arr(i+1) = arr(i)
               idx(i+1) = idx(i)
 11         continue
            i = l-1
 2          arr(i+1) = a
            idx(i+1) = b
 12      continue

         if (jstack .eq. 0) return

         ir = istack(jstack)
         l = istack(jstack - 1)
         jstack = jstack - 2
      else
         k = (l + ir) / 2
         temp = arr(k)
         arr(k) = arr(l+1)
         arr(l+1) = temp
         temp_int = idx(k)
         idx(k) = idx(l+1)
         idx(l+1) = temp_int

         if (arr(l) .gt. arr(ir)) then
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp

            temp_int = idx(l)
            idx(l) = idx(ir)
            idx(ir) = temp_int
         endif

         if (arr(l+1) .gt. arr(ir)) then
            temp = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir) = temp

            temp_int = idx(l+1)
            idx(l+1) = idx(ir)
            idx(ir) = temp_int
         endif

         if (arr(l) .gt. arr(l+1)) then
            temp = arr(l)
            arr(l) = arr(l+1)
            arr(l+1) = temp
            temp_int = idx(l)
            idx(l) = idx(l+1)
            idx(l+1) = temp_int
         endif

         i = l + 1
         j = ir
         a = arr(l+1)
         b = idx(l+1)

 3       continue

         i = i+1

         if (arr(i) .lt. a) goto 3

 4       continue

         j = j -1

         if (arr(j) .gt. a) goto 4
         if (j .lt. i) goto 5

         temp = arr(i)
         arr(i) = arr(j)
         arr(j) = temp

         temp_int = idx(i)
         idx(i) = idx(j)
         idx(j) = temp_int
         goto 3

 5       arr(l+1) = arr(j)
         arr(j) = a
         idx(l+1) = idx(j)
         idx(j) = b
         jstack = jstack + 2
         if (jstack .gt. nstack) pause "..." !'NSTACK too small in quicksort2'

         if (ir-i+1 .ge. j-1) then
            istack(jstack) = ir
            istack(jstack-1) = i
            ir = j-1
         else
            istack(jstack) = j-1
            istack(jstack-1) = l
            l = i
         endif
      endif
      goto 1
      end subroutine quicksort2
c-----------------------------------------------------------------------
      subroutine add_noise_scal(qin,fc1,fc2,fc3)
      implicit none
      include 'SIZE'
      include 'TSTEP'           ! TIME, DT
      include 'PARALLEL'        ! LGLEL
      include 'INPUT'           ! if3D
      include 'SOLN'            ! VX, VY, VZ, VMULT
      include 'GEOM'            ! XM1, YM1, ZM1
      real, intent(inout), dimension(lx1*ly1*lz1*lelt) :: qin
      real, intent(in) :: fc1,fc2,fc3
      real, dimension(lx1,ly1,lz1,lelt) :: q
      integer iel,ieg,il,jl,kl,nt
      real xl(ldim),mth_rand,fc(3),nmin,nmax,glmax,glmin
      fc(1)=fc1; fc(2)=fc2; fc(3)=fc3
      nt = nx1*ny1*nz1*nelt
      call copy(q(:,:,:,:),qin(:), nt)
      do iel=1,nelv
         do kl=1,nz1
            do jl=1,ny1
               do il=1,nx1
                  ieg = lglel(iel)
                  xl(1) = xm1(il,jl,kl,iel)
                  xl(2) = ym1(il,jl,kl,iel)
                  if (if3D) xl(ndim) = zm1(il,jl,kl,iel)
                  q(il,jl,kl,iel)=q(il,jl,kl,iel)+mth_rand(il,jl,kl,ieg,xl,fc)
               enddo
            enddo
         enddo
      enddo
      call dssum(q,lx1,ly1,lz1)
      call col2(q, vmult, nt)
      call dsavg(q)
      call bcdirSC(q)
      call copy(qin(:),q(:,:,:,:), nt) ! RESHAPE ARRAY TO 1D
      nmin=glmin(qin(:),nt); nmax=glmax(qin(:),nt)
      if(nid.eq.0)write(6,*)'noise scal min,max',nmin,nmax
      return
      end subroutine add_noise_scal
c-----------------------------------------------------------------------
      subroutine op_add_noise(qx, qy, qz)
   !     input random number to fields
         implicit none
         include 'SIZE'            ! NX1, NY1, NZ1, NELV, NID
         include 'TSTEP'           ! TIME, DT
         include 'PARALLEL'        ! LGLEL
         include 'INPUT'           ! if3D
         include 'SOLN'            ! VX, VY, VZ, VMULT
         include 'GEOM'            ! XM1, YM1, ZM1
   
         real, dimension(lx1,ly1,lz1,lelv) :: qx, qy, qz
         integer iel,ieg,il,jl,kl,nv
         real xl(LDIM),mth_rand,fc(3),nmin,nmax,glmax,glmin
         nv = nx1*ny1*nz1*nelv
   
         do iel=1,NELV
            do kl=1,NZ1
               do jl=1,NY1
                  do il=1,NX1
                     ieg = LGLEL(iel)
                     xl(1) = XM1(il,jl,kl,iel)
                     xl(2) = YM1(il,jl,kl,iel)
                     if (if3D) xl(NDIM) = ZM1(il,jl,kl,iel)
   
                     fc(1)=3.0e4; fc(2)=-1.5e3; fc(3)=0.5e5
                     qx(il,jl,kl,iel)=qx(il,jl,kl,iel)+mth_rand(il,jl,kl,ieg,xl,fc)

                     fc(1)=2.3e4; fc(2)=2.3e3; fc(3)=-2.0e5
                     qy(il,jl,kl,iel)=qy(il,jl,kl,iel)+mth_rand(il,jl,kl,ieg,xl,fc)

                     if (if3D) then
                        fc(1)=2.e4; fc(2)=1.e3; fc(3)=1.e5
                        qz(il,jl,kl,iel)=qz(il,jl,kl,iel)+mth_rand(il,jl,kl,ieg,xl,fc)
                     endif
                     
                  enddo
               enddo
            enddo
         enddo
   
   !     face averaging
         call opdssum(qx(:,:,:,:), qy(:,:,:,:), qz(:,:,:,:))
         call opcolv (qx(:,:,:,:), qy(:,:,:,:), qz(:,:,:,:), VMULT)
   
         call dsavg(qx(:,:,:,:))
         call dsavg(qy(:,:,:,:))
         if (if3D) call dsavg(qz(:,:,:,:))
         
         !Note: v*mask removes points at wall/inflow
         call bcdirVC(qx(:,:,:,:), qy(:,:,:,:), qz(:,:,:,:),v1mask,v2mask,v3mask)

         nmin=glmin(qx,nv);nmax=glmax(qx,nv)
         if(nid.eq.0)write(6,*)'noise vx min,max',nmin,nmax

         nmin=glmin(qy,nv);nmax=glmax(qy,nv)
         if(nid.eq.0)write(6,*)'noise vy min,max',nmin,nmax

         if (if3D) then 
         nmin=glmin(qz,nv);nmax=glmax(qz,nv)
         if(nid.eq.0)write(6,*)'noise vz min,max',nmin,nmax
         endif
         return
         end subroutine op_add_noise
c-----------------------------------------------------------------------
      subroutine add_symmetric_seed(qx, qy, qz, qp) !generate symmetric seed to fields
      implicit none
      include "SIZE"
      include "TOTAL"
      real, dimension(lx1,ly1,lz1,lelv) :: qx, qy, qz, qp
      integer iel,ieg,il,jl,kl,ntot
      real xlx,yly,zlz,alpha,x,y,z
      real glsc3,amp

      ntot = NX1*NY1*NZ1*NELV
      xlx = xmx - xmn
      yly = ymx - ymn
      zlz = zmx - zmn

      alpha = 2*pi/zlz

!     --> Create the initial velocity perturbation.

      do iel=1,NELV
         do kl=1,NZ1
            do jl=1,NY1
               do il=1,NX1

                  ieg = LGLEL(iel)
                  x = XM1(il,jl,kl,iel)
                  y = YM1(il,jl,kl,iel)
                  if (if3D) z = ZM1(il,jl,kl,iel)

!     -> Construct the perturbation. ! Note: Spanwise invariant.
                  qx(il,jl,kl,iel) = cos(alpha*z)*sin(2.*pi*y)
                  qz(il,jl,kl,iel) = -(2.*pi)/(alpha)*cos(alpha*z)*cos(2.*pi*y)
                  qp(il,jl,kl,iel) = cos(alpha*z)*cos(2.*pi*y)

               enddo
            enddo
         enddo
      enddo

      amp = glsc3(qx, bm1, qx, ntot)+ glsc3(qy, bm1, qy, ntot)
      if(if3D) amp = amp + glsc3(qz, bm1, qz, ntot)
      amp = 1e-6/(0.50d0*amp)
      call opcmult(qx, qy, qz, amp)
      call cmult(qp, amp, ntot)

      return
      end subroutine add_symmetric_seed
c-----------------------------------------------------------------------
      real function mth_rand(ix,iy,iz,ieg,xl,fc) !generate random number
      implicit none
      include 'SIZE'
      include 'INPUT' ! if3D
      integer ix,iy,iz,ieg
      real xl(LDIM), fc(3)
      mth_rand = fc(1)*(ieg+xl(1)*sin(xl(2))) + fc(2)*ix*iy+fc(3)*ix
      if (if3D) mth_rand = fc(1)*(ieg +xl(NDIM)*sin(mth_rand))+fc(2)*iz*ix+fc(3)*iz
      mth_rand = cos(1.e3*sin(1.e3*sin(mth_rand)))
      return
      end function mth_rand
c-----------------------------------------------------------------------
      subroutine outpost_vort(ux,uy,uz,name)
         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'
         real, intent(in) :: ux(lv),uy(lv),uz(lv)
         character(len=3)  :: name
         real wo1(lv),wo2(lv),wo3(lv),vort(lv,3)
         logical ifto_sav, ifpo_sav

         if(ifvor)then 
            
            call comp_vort3(vort, wo1, wo2, ux,uy,uz)

            ifto_sav = ifto
            ifpo_sav = ifpo
            ifto = .false.
            ifpo = .false.
            call outpost(vort(1,1), vort(1,2), vort(1,3), pr, t, name)
            ifto = ifto_sav
            ifpo = ifpo_sav
         endif

         return
      end subroutine outpost_vort
c-----------------------------------------------------------------------
      subroutine norm_grad(vx_, vy_, vz_, pr_, t_, norma)

         use krylov_subspace
         implicit none
         include 'SIZE'
         include 'TOTAL'

         real, intent(in), dimension(lv) :: vx_, vy_, vz_
         real, intent(in), dimension(lp) :: pr_
         real, intent(in), dimension(lt,ldimt) :: t_
         real, intent(out) :: norma

         real, dimension(lv) :: dudx, dudy, dudz
         real, dimension(lv) :: dvdx, dvdy, dvdz
         real, dimension(lv) :: dwdx, dwdy, dwdz

         real :: glsc3
         nv = nx1 * ny1 * nz1 * nelv

         ! gradient computation
         call gradm1(dudx, dudy, dudz, vx_, nelv)
         call gradm1(dvdx, dvdy, dvdz, vy_, nelv)
         if (if3D) call gradm1(dwdx, dwdy, dwdz, vz_, nelv)

         ! call dsavg(dudx); call dsavg(dudy); call dsavg(dudz)
         ! call dsavg(dvdx); call dsavg(dvdy); call dsavg(dvdz)
         ! call dsavg(dwdx); call dsavg(dwdy); call dsavg(dwdz)

         norma = 0.0d0

         norma = norma + glsc3(dudx, bm1s, dudx, nv) + glsc3(dudy, bm1s, dudy, nv)
         norma = norma + glsc3(dvdx, bm1s, dvdx, nv) + glsc3(dvdy, bm1s, dvdy, nv)

         if (if3D) then 
            norma = norma + glsc3(dudz, bm1s, dudz, nv)
            norma = norma + glsc3(dvdz, bm1s, dvdz, nv)
            norma = norma + glsc3(dwdx, bm1s, dwdx, nv) 
            norma = norma + glsc3(dwdy, bm1s, dwdy, nv)
            norma = norma + glsc3(dwdz, bm1s, dwdz, nv)
         endif 
      end subroutine norm_grad