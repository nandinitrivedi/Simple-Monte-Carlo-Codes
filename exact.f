c ******************************************************************
c		exact.f
c	PROGRAM WRITTEN BY NANDINI TRIVEDI
c	THEORETICAL PHYSICS GROUP, TIFR
c	DEC. 24, 1997
c******************************************************************
c	PROGRAM EVALUATES THE EXACT INTEGRAL OF F(x) = G1(x)*G2(x)
c	FOR x=[0,1], WHERE G1(x) is a Gaussian with mean a1 and 
c	standard deviation sig1 and G2(x) is a Gaussian with mean
c	a2 and standard deviation sig2
c	The result is evaluated exactly -- see notes -- in terms
c	of Error Functions.

	implicit none
c		read in the values of a1,a2,sig1,sig2
c		through a parameter statement
	include 'param.dat'
	double precision b,c,sigsq,sigma, pre1, pre2,term,arg1,arg2
	double precision pi,erf1,erf2,int
	pi = 4.d0*datan(1.d0)

	sigsq = sig1*sig1 + sig2*sig2
	sigma = sigsq/(sig1*sig1*sig2*sig2)
	b = (a1*sig2*sig2 + a2*sig1*sig1)/sigsq
	c = (a1*a1*sig2*sig2 + a2*a2*sig1*sig1)/sigsq
	pre2 = dsqrt(0.5d0*pi/sigma)
	term = dexp(0.5d0*sigma*(b*b-c))
	arg1 = (1.d0-b)*dsqrt(0.5d0*sigma)
	arg2 = b*dsqrt(0.5d0*sigma)
	write (6,*) arg1, arg2
c		NOTE: erf1 = Erf(arg1) and erf2 = Erf(arg2)
c		However, since in the cases studied here 
c		the arguments are large, we essentially get
c		Erf(arg1)=1 and Erf(arg2)=1
c		In case that is not the case, use mathematica or
c		Numerical Recipes to obtain the values of Erf(arg1) and
c		Erf(arg2) and input those below.
	erf1=1.d0
	erf2 = 1.d0
	pre1 = 2.d0*pi*sig1*sig2
	int=pre2*term*(erf1+erf2)/pre1
	write (6,*) 'PRODUCT OF TWO GAUSSIANS --'
	write (6,100) a1,sig1
	write (6,110) a2,sig2
	write (6,120) int
100	format('G1(x): mean a1=',f8.4, 5x, 'std. dev. sig1=',f8.4)
110	format('G2(x): mean a2=',f8.4, 5x, 'std. dev. sig2=',f8.4)
120	format('exact value of integral = ',f20.8)
	end
