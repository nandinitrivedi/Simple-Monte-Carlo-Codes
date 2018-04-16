c**************************************************************************
c		quad.f
c       PROGRAM WRITTEN BY NANDINI TRIVEDI
c       THEORETICAL PHYSICS GROUP, TIFR
c       DEC. 24, 1997
c******************************************************************
c       PROGRAM EVALUATES THE INTEGRAL OF F(x) = G1(x)*G2(x)
c       FOR x=[0,1], WHERE G1(x) is a Gaussian with mean a1 and
c       standard deviation sig1 and G2(x) is a Gaussian with mean
c       a2 and standard deviation sig2
c	USING QUADRATURE METHODS TO COMPUTE THE INTEGRAL
c***************************************************************************

	implicit none
	integer m, i
	double precision x,fx,sum,sum2,favg,sig
	double precision func

	write (6,*) 'input number of divisions m'
	read (5,*) m
	sum = 0.d0
	sum2 = 0.d0
	do i=0,m
	   x = dfloat(i)/dfloat(m)
	   fx = func(x)
	   sum = sum + fx
	   sum2 = sum2 + fx*fx
	enddo
	favg = sum/dfloat(m)
	sig = dsqrt((sum2/dfloat(m) - favg*favg)/dfloat(m))
	write (6,*) 'INTEGRAL EVALUATED USING QUADRATURE METHODS'
	write (6,*) 'number of divisions of xrange [0,1]',m
	write (6,*)
	write (6,100) favg, sig
100	format('value of integral=',e20.8,5x, 'error=', e20.8)
	write (6,*)
	end
	   
c==================================================================
c		DEFINITION OF FUNCTION FUNC
c===================================================================
	double precision function func(x)
	implicit none
	include 'param.dat'
	double precision x
	double precision  gauss1, gauss2, gauss
	
	gauss1 = gauss(x,a1,sig1)
	gauss2 = gauss(x,a2,sig2)
	func = gauss1 * gauss2
	end
	   
c==================================================================
c		DEFINITION OF FUNCTION GAUSS
c===================================================================
	double precision function gauss(x,a,sig)
	implicit none
	double precision x, a, sig, arg,pi, pre
	pi = 4.d0*datan(1.d0)
	pre = dsqrt(2.d0*pi)*sig
	
	arg = (x-a)*(x-a)/(2.d0*sig*sig)
	gauss = dexp(-arg)/pre
	end
