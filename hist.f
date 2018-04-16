c *********************************************************************
c		MAKE A HISTOGRAM OF THE X-POINTS 
c**********************************************************************
	implicit none
	integer nmax 
	parameter (nmax=1000)
	integer ndiv, n, nbin, nx, num(nmax), i
	double precision delx, xmin(nmax)
	double precision p,x
	
	open (unit=2, file='output',form='formatted',status='old')


	write (6,*) 'number of division of x[0,1]'
	read (5,*) ndiv
	write (6,*) 'number of x-points'
	read (5,*) nx
	do n = 1, ndiv
	   num(i) = 0
	enddo
	delx = 1.d0/dfloat(ndiv)
	xmin(1) = 0.d0
	do n = 2, ndiv
	   xmin(n) = xmin(n-1) + delx
	enddo
	do i = 1, nx
	   read (2,*) x, p
	   do n = 1, ndiv-1
	      if (x.ge.xmin(n) .and. x.lt.xmin(n+1)) then
	         nbin = n
	         go to 10
	      endif
	   enddo
10	   continue
	   num(nbin) = num(nbin) + 1
	enddo
	do n = 1, ndiv
c	   write (3,*) n, num(n)
	   write (3,*) xmin(n), num(n)/dfloat(nx)
	enddo
	end
	    
	
