c ******************************************************************
c               plotf.f
c       PROGRAM WRITTEN BY NANDINI TRIVEDI
c       THEORETICAL PHYSICS GROUP, TIFR
c       DEC. 24, 1997
c******************************************************************
c       PROGRAM EVALUATES THE TWO FUNCTIONS G1(x) and G2(x)
c	AND THE INTEGRAND F(x)=G1(x)*G2(x)
c       ON A FINE GRID FOR x=[0,1].
c	G1(x) is a Gaussian with mean a1 and standard deviation sig1 
c	and G2(x) is a Gaussian with mean a2 and standard deviation sig2
c	the file plot_g1 has the data with x in first column and
c					   G1(x) in second column
c	the file plot_g2 has the data with x in first column and
c					   G2(x) in second column
c	the file plot_int has the data with x in first column and
c            the integrand G1(x)*G2(x) in second column


	implicit none
	integer i
	double precision x, f, func
	do i = 0, 1000
c		choose some fine grid in the interval [0,1]
	   x = i/1000.d0
	   f = func(x)
	enddo
	end
c==================================================================
c		DEFINITION OF FUNCTION FUNC
c===================================================================
	double precision function func(x)
	implicit none
	include 'param.dat'
	double precision x
	double precision gauss1, gauss2, gauss

	open (unit=10,file='plot_g1',form='formatted',status='unknown')
	open (unit=11,file='plot_g2',form='formatted',status='unknown')
	open (unit=12,file='plot_int',form='formatted',status='unknown')

	gauss1 = gauss(x,a1,sig1)
	gauss2 = gauss(x,a2,sig2)
	func = gauss1 * gauss2
	write (10,100) x, gauss1
	write (11,100) x, gauss2
	write (12,100) x, func
100	format(2x,2e20.8)
	end
c==================================================================
c		DEFINITION OF GAUSSIAN FUNCTION
c===================================================================
	double precision function gauss(x,a,sig)
	implicit none
	double precision x, a, sig, arg,pi, pre
	pi = 4.d0*datan(1.d0)
	pre = dsqrt(2.d0*pi)*sig
	
	arg = (x-a)*(x-a)/(2.d0*sig*sig)
	gauss = dexp(-arg)/pre
	end
