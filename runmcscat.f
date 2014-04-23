      program test
      use omp_lib
      use mcscat
      implicit none
      integer lenalt, numphot
      parameter (lenalt=64, numphot=1000)
      integer z,x,y
      double precision, dimension(lenalt) :: parr, aerdens, tarr, dens,
     -     ratm, dtau, alt, ndens, aerext
      double precision, dimension(lenalt) :: singscat
      double precision, dimension(lenalt, numphot) ::  fluxmat
      double precision, dimension(lenalt,numphot,lenalt) :: pathmat
      double precision rstar, dist, phasecumarr(lenalt, lenth)
      double precision, dimension(mxmom, lenalt) :: phmom   
      double precision temp, scale, scatalb
      double precision heightavgout(lenalt)
      do 01 x=1,lenalt
         do 02 y=1,numphot
            do 03 z=1,lenalt
               pathmat(x, y, z) = 0.0
 03         enddo
 02      enddo
 01   enddo
      kb = 1.38D-23
      mp = 1.67D-27
c      mu = 29.0
      mu = 2.2
c      g = 9.8
      g = 9.8*0.91
c      rplan = 6371.
      rplan = 6371*2.678
c     temp = 250.
      temp = 600.0
      pi = acos(-1.0)
      do 10 i=1,lenalt
         parr(i) = 10**(-6+i*0.1)
         tarr(i) = temp
         dens(i) = parr(i)*100000/tarr(i)/kb
 10   enddo
      dens0 = maxval(dens)
      scale = kb*temp/(mp*mu*g)/1000.
      alt = log(dens0/dens)*scale
c      open(unit=10, file='icrccm_62.atm')
c      do 10 i=1, 11
c         read(10, *)
c 10   continue
c      do 20 i=1, lenalt
c         read(10,*) parr(i), tarr(i)
c 20   continue
c      close(10)
c      dens = parr/kb/tarr
c      dens0 = maxval(dens)
c      open(unit=11, file='height.dat')
c      do 30 i=1,lenalt
c         read(11,*) alt(lenalt-i+1)
c 30   continue
c      close(11)
      ratm = alt+rplan
      sig = 10**(-21.)
c      sig = 10**(-40.) ! for testing cases with essentially no scattering
      do 40 i=1, lenalt
         if ( ( parr(i) .gt.  0.01) .and. ( parr(i) .lt. 0.02 ) ) then
            aerdens(i) = 1.5*10**(14.0)
            dtau(i) = aerdens(i)*sig*10**3*(alt(i-1)-alt(i))
            aerext(i) = aerdens(i)*sig*10**3
c     set next layer to a very aerosol optical depth is produce more realistic results
            aerdens(i+1) = 1.5*10**(25.0)
            dtau(i+1) = aerdens(i+1)*sig*10**3*(alt(i)-alt(i+1))
            aerext(i+1) = aerdens(i+1)*sig*10**3
         endif
 40   enddo
c      alpha = 0.000275
c      alpha = 0.000132
      nstp = 1 + 0.000132
      bimp = .175
c      bimp = 0.0
      rstar = 695500.*0.2064
      dist = 1.5*10**8*0.0143
      gpar = 0.70
      scatalb = 1.0
      do 60 i=1, lenth
         tharr(i) = (i-1.0)*pi/float(lenth)
 60   continue
      phasefunc = 0.5*(1-gpar**2)/(1+gpar**2-2*gpar*cos(tharr))**1.5
      call cumtrapz(phasefunc, tharr, 1, lenth, lenth, phasecum)
      phasecum = phasecum/maxval(phasecum)
      thcum(1) = 0.0
      do 70 i=1,lenth
         if (i .eq. 1) then
            thcum(i) = 0.0
         else
            thcum(i) = 0.5*(tharr(i)+tharr(i-1))
         endif
         do 71 j=1,lenalt
            phasecumarr(j,i) = phasecum(i)
 71      enddo
 70   continue

      do 72 j=1,lenalt
         singscat(j) = scatalb
 72   enddo

      wl = 10.0
c     test values for wlarr, u1, u2, u3,u4
      do 80 i=1, lenwlarr
         wlarr(i) = i*.01
         u1(i) = 0; u2(i) = 0; u3(i) = 0; u4(i) = 0
 80   enddo
c      wl = 10.0
c      starang = 0.534*pi/180.0
c      starang = 2.0*pi/180.0
      call mcscatcalc(rplan, lenalt, numphot, alt, dens, nstp,  bimp,
     -     dtau,aerext,singscat, thcum, phasecumarr,rstar,dist, pathmat, 
     -     fluxmat, wl, lenwlarr, wlarr, u1, u2, u3, u4, heightavgout)
      end
