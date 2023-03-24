!  Grid Printer Proto-type 01
!       * Copy the value from old grid to new grid
!          using value of same point / interpolation
!       * Bathy data of grid or Nodall attributes
!
!                       Copyleft by Seizo Tanaka
!                                and C.H.Lab, University of Notre Dame
!
!  --- Release
!     2010.07.10. Kick off
!
!
MODULE griddata
     implicit none
     integer :: nnode, nelem, imode
     integer, allocatable :: nc(:,:)
     double precision, allocatable :: xy(:,:), value(:,:)
     integer :: npt
     double precision, allocatable :: xpt(:,:),valpt(:,:)
     character(len=80) :: fort14_old, fort14_new, fort13_old, fort13_new
END MODULE
!
MODULE lattice_tab
     implicit none
      integer, parameter :: ndiv = 1000
      integer :: ne_piece_sum
      integer :: ne_piece_add(ndiv,ndiv), ne_piece(0:ndiv+1,0:ndiv+1)
      integer, allocatable :: ne_piece_list(:), list_ehasn(:,:)
      double precision :: xmax(1:2), xmin(1:2)
END MODULE
!
!
PROGRAM MAIN
     use griddata, only:imode
     implicit none

!     write(6,*) ' Select mode:'
!     write(6,*) '   1. Print Bathydata onto New Grid'
!     write(6,*) "   2. Print Nodal Attributes onto New Grid's"
!     read(5,*) imode
      imode = 2
!
     call read14
     call search
     select case(imode)
       case(1)
         call print_bath
       case(2)
         call print_13
       case default
         write(6,*)  '   !!!! You must select 1 or 2 '
         stop
     end select

END PROGRAM MAIN
!
!--------------------------------------------------------------------------------------------------
  SUBROUTINE read14
!--------------------------------------------------------------------------------------------------
     use griddata
     implicit none
     integer :: i,j,k,n,m
!
     open(10,file='Print13.ctrl',status='old',action='read')
      read(10,'(a)') fort14_old
      read(10,'(a)') fort14_new
!
      if( imode==2 ) then
        read(10,'(a)') fort13_old
        read(10,'(a)') fort13_new
      endif
     close(10)
!
     open(14,file=fort14_old,status='old',action='read')
     write(6,*) ' * READING: ', trim(adjustl(fort14_old))
      read(14,*)
      read(14,*) nelem, nnode
      call showbar(nnode+nelem,0)
      allocate( xy(1:3,1:nnode), nc(1:3,1:nelem) )
      do n = 1, nnode
        read(14,*) i, (xy(j,n),j=1,2)
        if( i /= n) then
          write(6,*) ' Error 501. System requires SEQUENTIAL NODE NUMBER'
          stop
        endif
      call showbar(nnode+nelem,n)
      enddo
      do m = 1, nelem
        read(14,*) i, k, (nc(j,m),j=1,3)
        call showbar(nnode+nelem,m+nnode)
      enddo
      write(6,*) 
     close(14)
!
     open(14,file=fort14_new,status='old',action='read')
     write(6,*) ' * READING:', trim(adjustl(fort14_new))
      read(14,*) 
      read(14,*) n, npt
      call showbar(npt,0)
      allocate( xpt(1:2,1:npt) )
      do n = 1, npt
        read(14,*) i, (xpt(j,n),j=1,2)
        if( i /= n) then
          write(6,*) ' Error 502. System requires SEQUENTIAL NODE NUMBER'
          stop
        endif
      call showbar(npt,n)
      enddo
      write(6,*) 
     close(14)
!
  END SUBROUTINE read14
!
  SUBROUTINE print_bath  ! Now under constructing
     use griddata
     use lattice_tab, only: list_ehasn
     implicit none
     character(len=1000) :: chara

     open(14,file=fort14_new,status='old',action='read')
     open(140,file=fort14_new//'-up',status='replace',action='write')
       read(14,'(a1000)') chara
       write(140,'(a)') trim(adjustl(chara))   !header
       read(14,'(a1000)') chara
       write(140,'(a)') trim(adjustl(chara))   !nelem, nnode
  END SUBROUTINE print_bath
!
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
  SUBROUTINE print_13
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
     use griddata
     use lattice_tab, only: list_ehasn
     implicit none
!
     integer :: i, j, k, n, nondef, nattr, ii(1)
     double precision :: rjunk, r(3), d(3)
     double precision, allocatable :: defval(:,:)
     integer, allocatable :: nvpn(:)
     character(len=200) :: cjunk, aname, attrname
     character(len=20) :: c(13)
!
     integer, allocatable :: nprop(:)
     allocate( nprop(1:npt) )
!
     open(13, file=fort13_old,status='old',    action='read')
     open(130,file=fort13_new,status='replace',action='write')
!
     read(13,'(a)') cjunk   !Agrid
     write(130,'(a)') trim(adjustl(cjunk))
! NNODE
      read(13,*)   i
      write(cjunk,*) npt
      write(130,'(a)') trim(adjustl(cjunk))
! NATTR
      read(13,*) nattr !NATTR
      write(cjunk,*) nattr
      write(130,'(a)') trim(adjustl(cjunk))
!
      allocate(nvpn(1:nattr), defval(1:12,1:nattr))
!
      do i = 1, nattr

        read (13,'(a)') attrname   !AttrName
        write(130,'(x,a)') trim(adjustl(attrname))
        read (13,'(a)') cjunk         !Units
        write(130,'(x,a)') trim(adjustl(cjunk))

        read (13,*) nvpn(i)
        write(cjunk,*) nvpn(i)
        write(130,'(x,a)') trim(adjustl(cjunk))

        read (13,*) ( defval(k,i), k = 1,nvpn(i) )
        do k = 1, nvpn(i)
          write(c(k),'(f20.6)') defval(k,i)
        enddo
        write(130,'(20(x,a))') ( trim(adjustl(c(k))), k=1,nvpn(i) )

      enddo
!
      allocate (value(1:maxval(nvpn(1:nattr)),1:nnode))
      allocate (valpt(1:maxval(nvpn(1:nattr)),1:npt))
!
      do i = 1, nattr
        read(13,'(a)') attrname ! AttrName
        write(130,'(a)') trim(adjustl(attrname))
        write(6,'(a)') trim(adjustl(attrname))
        read(13,*) nondef
!
          do j =1, nvpn(i)
            value(j,1:nnode) = defval(j,i)
          enddo
          do n = 1, nondef
            read(13,*)j, (value(k,j), k=1,nvpn(i))
          enddo
! Print to new grid
          do n = 1, npt
            select case(list_ehasn(1,n))
              case(1)
                do j = 1, nvpn(i)
                  valpt(j,n) = value(j,list_ehasn(2,n))
                enddo
              case(2)
                call inp_ele(nvpn(i), n, list_ehasn(2,n) )
              case(0)
                do j = 1, nvpn(i)
                  valpt(j,n) = -99999.00
                enddo
            end select
          enddo
!    Classification Filter
          aname='average_horizontal_eddy_viscosity_in_sea_water_wrt_depth'
          if( trim(adjustl(aname)) == trim(adjustl(attrname)) ) then
            write(6,*) 'Do you want to classify to two class (20.0 and 2.0)? [Y/N]'
            read(5,'(a)') cjunk
            if( (trim(adjustl(cjunk)) == 'Y').or.(trim(adjustl(cjunk)) == 'y') ) then
              do n = 1, npt
                rjunk = valpt(1,n)
                if ( rjunk < -5000.0d0 ) cycle
                if( dabs(rjunk-20.0d0) < dabs(rjunk-2.0d0) ) then
                   valpt(1,n) = 20.0d0
                else
                   valpt(1,n) =  2.0d0
                endif
              enddo
              write(6,*) 'Classified!'
            endif
          endif
          aname='primitive_weighting_in_continuity_equation'
          if( trim(adjustl(aname)) == trim(adjustl(attrname)) ) then
            write(6,*) 'Do you want to classify to three class (0.03, 0.02 and 0.005)? [Y/N]'
            read(5,'(a)') cjunk
            if( (trim(adjustl(cjunk)) == 'Y').or.(trim(adjustl(cjunk)) == 'y') ) then
              d(1) = 0.03d0
              d(2) = 0.02d0
              d(3) = 0.005d0
              do n = 1, npt
                rjunk = valpt(1,n)
                if ( rjunk < -5000.0d0 ) cycle
                do j = 1, 3
                  r(j) = dabs( rjunk - d(j) )
                enddo
                ii = minloc(r(1:3))
                valpt(1,n) = d(ii(1))
              enddo
              write(6,*) 'Classified!'
            endif
          endif
       
!Count Non Def Value 
          nprop(1:npt) = 0
          NPT_loop:do n = 1, npt
            do j = 1, nvpn(i)
              if( dabs(valpt(j,n)-defval(j,i)) > 1.0d-5 ) then
                nprop(n) = 1
                exit
              endif
            enddo
          enddo NPT_loop
          write(cjunk,*) sum( nprop(1:npt) )
          write(130,'(x,a)') trim(adjustl(cjunk))
          do n = 1, npt
            if( nprop(n) == 0 ) cycle
            write(c(1),*) n
            do j = 1, nvpn(i)
              write(c(j+1),'(f20.6)') valpt(j,n)
            enddo
            write(130,'(20(x,a))') ( trim(adjustl(c(k))), k=1,nvpn(i)+1 )
          enddo
      enddo
!
  END SUBROUTINE print_13
!
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
  SUBROUTINE inp_ele(j, n, m)
     use griddata, only: nc, value, valpt, xy, xpt
     implicit none
     integer, intent(in) :: j, n, m
     integer :: i, n1, n2, n3
     double precision :: xc, yc, x1, x2, x3, y1, y2, y3, a0, a1, a2, a3
!
         xc = xpt(1,n)
         yc = xpt(2,n)
         n1 = nc(1,m)
         n2 = nc(2,m)
         n3 = nc(3,m)
         x1 = xy(1,n1)
         x2 = xy(1,n2)
         x3 = xy(1,n3)
         y1 = xy(2,n1)
         y2 = xy(2,n2)
         y3 = xy(2,n3)
         a0 =  (x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)
         a1 = ((xc - x2) * (yc - y3) - (xc - x3) * (yc - y2) ) / a0
         a2 = ((x1 - xc) * (y1 - y3) - (x1 - x3) * (y1 - yc) ) / a0
         a3 = ((x1 - x2) * (y1 - yc) - (x1 - xc) * (y1 - y2) ) / a0
!
         do i = 1, j
           valpt(i,n) = value(i,n1) * a1 + value(i,n2) * a2 + value(i,n3) * a3
         enddo

  END SUBROUTINE inp_ele
!
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
  SUBROUTINE search
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
     use griddata, only: nnode, nelem, xy, nc, npt, xpt
     use lattice_tab
     implicit none
     call count_list(nnode, nelem, nc, xy, xmax, xmin, ndiv, ne_piece_add, ne_piece, ne_piece_sum)
     allocate( ne_piece_list(ne_piece_sum) )
     call mkelocation( nnode, nelem, nc, xy, &
                       ndiv, xmax, xmin, ne_piece_sum, ne_piece_add, ne_piece, ne_piece_list)
     call search_inelm(nnode,nelem, xy, nc, npt, xpt)
!
  END SUBROUTINE search
!
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
  SUBROUTINE count_list(nnode, nelem, nc, xy, xmax, xmin, ndiv, ne_piece_add, ne_piece, ne_piece_sum)
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
      implicit none
      integer,intent(in) :: nnode, nelem, nc(1:3,1:nelem), ndiv
      double precision, intent(in) :: xy(1:3,1:nnode)
!
      double precision, intent(out) :: xmax(1:2), xmin(1:2)
      integer, intent(out) :: ne_piece_sum, ne_piece(0:ndiv+1,0:ndiv+1), ne_piece_add(1:ndiv,1:ndiv)
!
      integer :: i, j, im , m, n, ix(2,3), istart
      double precision :: dx(2)
!
      do i = 1, 2
         xmin(i) = xy(i,nc(1,1))
         xmax(i) = xy(i,nc(1,1))
      enddo
      do m = 1, nelem
        do i = 1, 2
          do j = 1, 3
            xmin(i) = dmin1(xmin(i), xy(i,nc(j,m)))
            xmax(i) = dmax1(xmax(i), xy(i,nc(j,m)))
          enddo
        enddo
      enddo
!
      write(6,*) '   Domain Info'
      write(6,*) xmin(1), xmin(2)
      write(6,*) xmax(1), xmax(2)
!
      do i = 1, 2
        xmax(i) = xmax(i) + 1.d0
      enddo
      dx(1:2) = ( xmax(1:2)-xmin(1:2) ) / ndiv
      ne_piece(:,:) = 0
      do m = 1, nelem
        do j = 1, 3
          n = nc(j,m)
          do i = 1, 2
             ix(i,j) = int( (xy(i,n)-xmin(i)) / dx(i) ) + 1
          enddo
        enddo
        do i = MINVAL(ix(1,1:3)), MAXVAL(ix(1,1:3))
          do j = MINVAL(ix(2,1:3)), MAXVAL(ix(2,1:3))
            ne_piece(i,j) = ne_piece(i,j) + 1
          enddo
        enddo
      enddo
      istart = 0
      do i = 1, ndiv
         do j = 1, ndiv
            ne_piece_add(i,j) = istart
            istart = istart + ne_piece(i,j)
         enddo
      enddo
      ne_piece_sum = sum( ne_piece(1:ndiv,1:ndiv) )
      write(6,*) nelem, ne_piece_sum
   END SUBROUTINE count_list
!
!
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
   SUBROUTINE mkelocation( nnode, nelem, nc, xy, &
                           ndiv, xmax, xmin, ne_piece_sum, ne_piece_add, ne_piece, ne_piece_list )
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
      implicit none
      integer, intent(in) :: nnode, nelem, nc(1:3,1:nelem),ndiv, ne_piece_sum
      double precision, intent(in) :: xy(1:3,1:nnode), xmax(1:2), xmin(1:2)
      integer, intent(in)  :: ne_piece_add(1:ndiv,1:ndiv)
      integer, intent(out) :: ne_piece(0:ndiv+1,0:ndiv+1), ne_piece_list(1:ne_piece_sum)
!
      integer :: im, n, m, i, j, ix(1:2,1:3)
      double precision :: dx(1:2)
!
      dx(1:2) = ( xmax(1:2)-xmin(1:2) ) / ndiv
!
      ne_piece(1:ndiv,:) = 0
      do m = 1, nelem
        do j = 1, 3
          n = nc(j,m)
          do i = 1, 2
            ix(i,j) = int( (xy(i,n)-xmin(i)) / dx(i) ) + 1
          enddo
        enddo
        do i = MINVAL(ix(1,1:3)), MAXVAL(ix(1,1:3))
          do j = MINVAL(ix(2,1:3)), MAXVAL(ix(2,1:3))
            ne_piece(i,j) = ne_piece(i,j) + 1
            ne_piece_list( ne_piece(i,j)+ne_piece_add(i,j) ) = m
          enddo
        enddo
      enddo
   END SUBROUTINE mkelocation
!
!
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
   SUBROUTINE search_inelm(nnode,nelem, xy, nc, npt, xpt)
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
!$    use OMP_LIB
      use lattice_tab
      implicit none
      integer,intent(in) :: nnode, nelem, nc(1:3,1:nelem), npt
      double precision, intent(in) :: xy(1:3,1:nnode), xpt(1:2,1:npt)
      double precision :: xc, yc
!
      integer :: i, j, k, n, mi, ix, iy, iadd, ipt, iptm
      double precision :: x(1:3), y(1:3), a1, a2, a3, a0, dx, dy, dist
      double precision :: eps = 1.0d-03
      integer :: mythreads, numthreads
      integer :: idiv
      integer,allocatable :: nomp_map(:)
      double precision, allocatable :: dist_near(:)
      character(len=200) :: cjunk
!
      idiv = ndiv
!
      allocate( dist_near(1:npt) )
      dx = ( xmax(1)-xmin(1) )
      dy = ( xmax(2)-xmin(2) )
!
      dist_near(1:npt) = dsqrt(dx*dx + dy*dy) * 2.0d0
!
      allocate( list_ehasn(1:2,1:npt) )
      list_ehasn(1:2,1:npt) = 0
!
      dx = ( xmax(1)-xmin(1) ) / ndiv
      dy = ( xmax(2)-xmin(2) ) / ndiv
!
      write(6,*) '    * Search location of points'
      mythreads = 0
      numthreads = 1
!$    numthreads = omp_get_max_threads()
      j = 0
      allocate (nomp_map(1:npt))
      do i = 1, numthreads
        do n = 1, int(npt/numthreads)
          j = j + 1
          nomp_map(j) = (n-1)*numthreads+i
        enddo
      enddo
      do i = int(npt/numthreads)*numthreads+1, npt
        j = j + 1
        nomp_map(i) = i
      enddo
!
      call showbar(int(npt/numthreads),0)
!
!$ write(*,*)    '        OpenMP mode',omp_get_max_threads(),'Threads'
!$OMP PARALLEL default(none) &
!$OMP          private( ipt, iptm, xc, yc, ix, iy, i, mi, k, dist, x, y, a0, a1, a2, a3, mythreads ) &
!$OMP          shared( npt, xpt, xmin, dx, dy, idiv, ne_piece, ne_piece_list, ne_piece_add,    &
!$OMP                  xy, nc, eps, list_ehasn, dist_near, numthreads, nomp_map)
!$    mythreads = omp_get_thread_num()
!$OMP DO
      POINT_LOOP:do iptm = 1, npt
        ipt = nomp_map(iptm)
!        ipt = iptm
        
!$      if (mythreads==0) then
          call showbar(int(npt/numthreads),iptm)
!$      endif
        xc = xpt(1,ipt)
        yc = xpt(2,ipt)
        ix = int( (xc -xmin(1)) / dx ) + 1
        iy = int( (yc -xmin(2)) / dy ) + 1
        ix = max(0,ix); ix = MIN(idiv+1,ix)
        iy = max(0,iy); iy = MIN(idiv+1,iy)
        do i = 1, ne_piece(ix,iy)
           mi = ne_piece_list(i+ne_piece_add(ix,iy))
           do k = 1, 3
             x(k) = xy(1,nc(k,mi))
             y(k) = xy(2,nc(k,mi))
           enddo
           do k = 1, 3
             dist = dsqrt(( x(k) - xc ) * ( x(k) - xc ) + ( y(k) - yc ) * ( y(k) - yc ))
             if( dist <= dist_near(ipt) ) then
                 dist_near(ipt) = dist
                 list_ehasn(2,ipt) = nc(k,mi)
             endif
             if( dist <= eps*eps ) then
                list_ehasn(1,ipt) = 1            !same node
                list_ehasn(2,ipt) = nc(k,mi)
                cycle POINT_LOOP
             endif
           enddo
          
           a0 =  (x(1) - x(2)) * (y(1) - y(3)) - (x(1) - x(3)) * (y(1) - y(2))
           a1 = (( xc  - x(2)) * ( yc  - y(3)) - ( xc  - x(3)) * ( yc  - y(2)) ) / a0
           a2 = ((x(1) -  xc ) * (y(1) - y(3)) - (x(1) - x(3)) * (y(1) -  yc ) ) / a0
           a3 = ((x(1) - x(2)) * (y(1) -  yc ) - (x(1) -  xc ) * (y(1) - y(2)) ) / a0
           if( (a1 >= -eps ) .and. ( a2 >= -eps ) .and. ( a3 >= -eps ) ) then
              list_ehasn(1,ipt) = 2              !interpolate
              list_ehasn(2,ipt) = mi
              cycle POINT_LOOP
           endif
        enddo
      enddo POINT_LOOP
!$OMP END DO
!$OMP END PARALLEL
      write(6,*)
      i=0; j=0; k=0
      do ipt = 1, npt
        select case(list_ehasn(1,ipt))
          case(0)
            i = i + 1
          case(1)
            j = j + 1
          case(2)
            k = k + 1
        end select
      enddo
      write(6,*) '      ---New Grid Points Data-------------'
      write(6,*) '            Same location:', j
      write(6,*) '            Interpolate  :', k
      write(6,*) '            NOT FOUND    :', i
      write(6,*) '      ------------------------------------'
      
      write(6,*) ' !!! Do you want to use the value of Nearest point for NOT FOUNDED point? [Y/N]'
      read(5,'(a)') cjunk
      if( (trim(adjustl(cjunk)) == 'Y').or.(trim(adjustl(cjunk)) == 'y') ) then
        do ipt = 1, npt
           if ( list_ehasn(1,ipt) == 0 ) then
!              write(6,*) ipt, dist_near(ipt), list_ehasn(2,ipt)
              if( list_ehasn(2,ipt) /=0 ) then
                   list_ehasn(1,ipt) = 1
              else
                write(156,*)  ipt, xpt(1,ipt), xpt(2,ipt)
              endif
           endif
        enddo
        write(6,*)
        i=0; j=0; k=0
        do ipt = 1, npt
          select case(list_ehasn(1,ipt))
            case(0)
              i = i + 1
            case(1)
              j = j + 1
            case(2)
              k = k + 1
          end select
        enddo
        write(6,*) '      ---Re-New Grid Points Data-------------'
        write(6,*) '            Same location:', j
        write(6,*) '            Interpolate  :', k
        write(6,*) '            NOT FOUND    :', i
      endif
      deallocate( dist_near )
  END SUBROUTINE search_inelm
!
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
   subroutine showbar(ntotal, now)
!-----+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+
   implicit none
   integer, intent(in) :: ntotal, now
   integer, save :: ncount, nout
   integer, parameter :: ndiv=50, ishow=10
   integer :: n
!
   if(now == 0) then
     ncount = 0
     nout = int(dble(ntotal/ndiv))
     n = mod(ntotal,ndiv)
     if ( n /= 0 ) nout = nout + 1
     return
   endif
!
!
   if( now >= ncount * nout ) then
      ncount = ncount + 1
      n = mod(ncount,ishow)
      if( n == 0 ) then
         write(6,'(a,$)') '+'
      else
         write(6,'(a,$)') '-'
      endif
      if( ncount==ndiv ) then
         write(6,*)
      end if
   endif
   end subroutine showbar
!
   
