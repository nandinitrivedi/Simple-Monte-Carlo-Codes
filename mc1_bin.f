c ******************************************************************
c               mc1_bin.f
c       PROGRAM WRITTEN BY NANDINI TRIVEDI
c       THEORETICAL PHYSICS GROUP, TIFR
c       DEC. 24, 1997
c******************************************************************

c**************************************************************************
c       PROGRAM EVALUATES THE INTEGRAL OF F(x) = G1(x)*G2(x)
c       FOR x=[0,1], WHERE G1(x) is a Gaussian with mean a1 and
c       standard deviation sig1 and G2(x) is a Gaussian with mean
c       a2 and standard deviation sig2
c       USING MONTE CARLO METHODS.
c
c       This is the simplest form OF MONTE CARLO WITHOUT ANY IMPORTANCE
c       SAMPLING that means that a list of x is generated at random and
c       F(x) is evaluated at those randomly selected points

c==================================================================
c	THE ERROR ON THE INTEGRAL VALUE IS COMPUTED BY BINNING THE DATA
c	WHICH REMOVES THE CORRELATIONS BETWEEN SUCCESSIVELY CHOSEN XPOINTS
c=========================================================================

	implicit none
	integer m, idum, i
	double precision x,sum
	double precision func, ran3
	integer nbin, nb, lbin
	parameter (nbin=10)
	double precision fbin(nbin), average, error 

	write (6,*) 'number of points m'
	read (5,*) m
	write (6,*) 'seed for random number generator'
	read (5,*) idum

	do nb = 1, nbin
	   fbin(nb) = 0.d0
	enddo
c		number of xpoints within a bin=lbin
	lbin = idint((dfloat(m)/dfloat(nbin)) + 0.1d0)
	do i=1,m
c		ran3 is a random number generator in the range [0,1]
c		nb is the bin number
	   nb = idint(dfloat(i-1)/dfloat(lbin)) + 1
	   x = ran3(idum)
	   fbin(nb) = fbin(nb) + func(x)
        enddo
	do nb = 1, nbin
           fbin(nb)  = fbin(nb)/dfloat(lbin)
	enddo
	write (6,*) (fbin(nb),nb=1, nbin)
c		calculate average
	average = 0.d0
	do nb = 1, nbin
	   average = average + fbin(nb)
	enddo
	average = average/dfloat(nbin)
c		compute error of average
	sum = 0.d0
	do nb = 1, nbin
	   sum = sum + (fbin(nb)-average)**2
	enddo
	error = dsqrt(sum/dfloat(nb-1))
        write (6,*) 'INTEGRAL EVALUATED USING MONTE CARLO'
        write (6,*) 'WITHOUT IMPORTANCE SAMPLING'
	write (6,*) 'PROPER ERROR ESTIMATION'
        write (6,*) 'number of xpoints=',m
        write (6,*)
        write (6,*) '=================================================='
        write (6,100) average, error
100     format('value of integral=',e20.8,5x, 'error=', e20.8)
        write (6,*) '=================================================='

        end
c==================================================================
c		DEFINITION OF FUNCTION FUNC
c===================================================================
	double precision function func(x)
	implicit none
	include 'param.dat'
	double precision x
	double precision gauss1, gauss2, gauss
	
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
	pre =dsqrt(2.d0*pi)*sig
	
	arg = (x-a)*(x-a)/(2.d0*sig*sig)
	gauss = dexp(-arg)/pre
	end
C====================================================================== 
C       PROGRAM: ran3.f  
C       TYPE   : main 
C       PURPOSE: 
C       I/O    :
C	
C       VERSION: 8-1-90 
C       COMMENT: TAKEN FROM  NUMERICAL RECIPES
C=============================================
        double precision function ran3(idum)
	implicit none
	integer mbig,mseed,mz,ma,idum,iff,mk,mj
	integer i,k,ii
	integer inext,inextp
	double precision fac
        Parameter (mbig=1000000000,Mseed=161803398,Mz=0,fac=1./Mbig)
        Dimension MA(55)
        common/saved/iff,inext,inextp,ma
        save
        if (idum.lt.0.or.iff.eq.0) then
          iff=1
          mj=mseed-iabs(idum)
          mj=mod(mj,mbig)
          ma(55)=mj
          mk=1
          do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if (mk.lt.mz) mk=mk+mbig
          mj=ma(ii)
11        continue
          do 13 k=1,4
            do 12 i=1,55
              ma(i)=ma(i)-ma(1+mod(i+30,55))
              if (ma(i).lt.mz)ma(i)=ma(i)+mbig
12          continue
13        continue
          inext=0        
          inextp=31
          idum=1
      end if
      inext=inext+1
      if (inext.eq.56) inext=1
      inextp=inextp+1
      if (inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if (mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      ran3=mj*fac
      return
      end	

