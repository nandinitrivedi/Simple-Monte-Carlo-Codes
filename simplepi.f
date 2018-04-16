C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: simplepi.f 
C TYPE   : main
C PURPOSE: example for Monte Carlo (written in Fortran 77)
C VERSION: 05-FEB-96
C COMMENT: This simplest possible Monte Carlo Program generates
C          random x and y between -1 and 1 and checks whether the
C          points fall within the circle.
C========+=========+=========+=========+=========+=========+=========+=$
      program pi
      parameter (N=1000000)
      idum=-91914
      pi=acos(-1.)
      nhit=0
      do i=1,N
         x=2.*(ran3(idum))-1.
         y=2.*(ran3(idum))-1.
         if (x**2+y**2.lt.1) then
            nhit=nhit+1
         end if
         if (mod(i,1000).eq.0) print*,i,4*nhit/real(i)
      end do
      end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: ran3.f  
C TYPE   : function
C PURPOSE: generate random numbers
C I/O    :
C VERSION: 27 MAI 94
C COMMENT: Initialize idum with negative integer
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      real function ran3(idum)
      Parameter (mbig=1000000000,Mseed=161803398,Mz=0,fac=1./Mbig)
      Dimension MA(55)
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
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if (ma(i).lt.mz)ma(i)=ma(i)+mbig
12        continue
13      continue
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
      end	
