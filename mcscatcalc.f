      MODULE mcscat
      implicit none
      integer i, j, k, lenth, dint, darrlen, maxscat, 
     -     lenrtest, lenwlarr
      integer lenthcum, altind, mxmom
      double precision dl
      parameter ( lenth=10000, darrlen=10000, mxmom=200, maxscat=10,
     -     lenrtest=400, lenwlarr = 100, dl=10)
      double precision kb, mp, mu, g, rplan, nstp
c      double precision, dimension(lenp):: parr, tarr, dens, alt, dtau,
c     -     ratm, slfluxnew, sltaunew, aerdens, ndens, slflux, sltau, 
c     -     refth, refminr, paths, heightavg
      double precision, dimension(lenwlarr) :: wlarr, u1, u2, u3, u4
      double precision, dimension(lenrtest) :: heightavgtest, rangtest
c      double precision, dimension(lenp, numphot, lenp) :: pathtot
c      double precision, dimension(lenp, numphot) :: pathtmp
c      double precision, dimension(lenp, numphot):: thfin, numfin, drfin
c      double precision, dimension(numphot) :: thtmp, numtmp, drtmp
c      double precision, dimension(lenp+4):: rtmp, atmp, ntmp
      double precision, dimension(darrlen):: darr, aerextnew, rparr      
      double precision, dimension(darrlen):: xnew, rnew, taucum, dcum
      double precision tharr(lenth), phasecum(lenth), thcum(lenth)
      double precision phasefunc(lenth), thlim, zp, bimp, rfinang, wl
      double precision dens0, sig, tanh, gpar, pi, ch, densstp, minaer
      double precision ds, dmax, tautot, th, xp, yp, altflux, dflux
      double precision psi, altflux1, angth, starang, slope,
     -     dr
      CONTAINS

      real function refpathcalc(xp, yp, th, rplan, alt, lenalt)
      integer, intent(in) :: lenalt
      double precision, intent(in) :: xp, yp, th, rplan
      double precision, intent(in) :: alt(lenalt)
      double precision dtot(4), a, b, c, rnew
      integer k
      if (sqrt(xp**2+yp**2+zp**2)-rplan .lt. 0.0) then
         refpathcalc = -9999
         return
      endif
      rnew = rplan + maxval(alt)
      a = 1.0
      b = 2*( yp*cos(th)+xp*sin(th) )
      c = xp**2+yp**2-rnew**2
      dtot(1) = (-b+sqrt(b**2-4*a*c))/(2*a)
      dtot(2) = (-b-sqrt(b**2-4*a*c))/(2*a)
      c = xp**2+yp**2-rplan**2
      if (b**2-4*a*c .lt. 0) then
         dtot(3) = -9999.
         dtot(4) = -9999.
      else
         dtot(3) = (-b+sqrt(b**2-4*a*c))/(2*a)
         dtot(4) = (-b-sqrt(b**2-4*a*c))/(2*a)
      endif
      if (maxval(dtot) .le. 0.0) then
         refpathcalc = -9999.
         return
      endif
      refpathcalc = 9999.
      do 01 k=1, 4
         if ( (dtot(k) .lt. refpathcalc) .and. (dtot(k).gt.10E-6) ) then
            refpathcalc = dtot(k)
         endif
 01   continue
      return
      end function
c
      real function pathcalc(xp, yp, zp, th, phi, rplan, alt, lenalt)
      integer, intent(in) :: lenalt
      double precision, intent(in) :: alt(lenalt), xp, yp, zp, th, phi,
     -     rplan
      double precision dtot(4), a, b, c, rnew
      integer k
      if (sqrt(xp**2+yp**2+zp**2)-rplan .lt. 0.0) then
         pathcalc = -9999
         return
      endif
      rnew = rplan+maxval(alt)
      a = 1+tan(th)**2+tan(phi)**2
      b = 2*(xp+yp*tan(th)+zp*tan(phi))
      c = xp**2+yp**2+zp**2-rnew**2
      dtot(1) = (-b+sqrt(b**2-4*a*c))/(2*a)
      dtot(2) = (-b-sqrt(b**2-4*a*c))/(2*a)
      c = xp**2+yp**2+zp**2-rplan**2
      if (b**2-4*a*c .lt. 0) then
         dtot(3) = -9999.
         dtot(4) = -9999.
      else
         dtot(3) = (-b+sqrt(b**2-4*a*c))/(2*a)
         dtot(4) = (-b-sqrt(b**2-4*a*c))/(2*a)
      endif
      if (maxval(dtot) .le. 0.0) then
         pathcalc = -9999
         return
      endif
      pathcalc = 9999.
      do 1 k=1, 4
         if ( (dtot(k) .lt. pathcalc) .and. (dtot(k) .gt. 10E-6) ) then
            pathcalc = dtot(k)
         endif
 1    end do
      pathcalc = pathcalc*sqrt(a)
      return
      end function
c
      subroutine refcalc(xp0, yp0, th0, rtmp, ntmp, atmp, rplan,
     -     thfin, distfin, alt, lenalt, paths)
      integer, intent(in) :: lenalt
      double precision, intent(in) :: xp0, yp0, th0, alt(lenalt), rplan
      double precision, dimension(lenalt+4), intent(in):: rtmp, ntmp, 
     -     atmp
      integer dint, darrlen, q
      parameter (darrlen=10000)
      double precision distfin, th, dmax, dl, thfin
      double precision rp, rnew, n, dndh, eps, xp, yp
      double precision, dimension(darrlen):: darr, aerextnew, rparr
      double precision, dimension(darrlen) :: xparr, yparr, zparr
      double precision phi, beta, val(4), xpnew, ypnew, rpnew
      double precision nnew, maxr
      double precision, dimension(lenalt), intent(out) :: paths
      dmax =  refpathcalc(xp0, yp0, th0, rplan, alt, lenalt)
      dmax = dmax*2. !increase dmax to ensure the entire atmosphere is traversed
      xp = xp0; yp = yp0; th = th0;
      dl = 10.
      if (dmax .lt. dl) then
         thfin = th
         distfin = sqrt(xp**2+yp**2)
         return
      endif
      dint = dmax/dl+1
      do 02 q=1, dint
         darr(q) = (q-1)*dl
 02   continue
      rp = sqrt(xp**2+yp**2)
      n = frefract(rp, rtmp, ntmp, lenalt)
      dndh = 0.0
      eps = 0.0
      phi = atan(xp/yp)
      beta = phi-th+pi/2.0
      val(1) = rp-rplan
      val(2) = phi
      val(3) = beta
      val(4) = eps
      maxr = maxval(alt) + rplan
      do 03 q=1, dint
         if (sqrt(xp**2+yp**2) .gt. maxr+1E-6) exit
         xpnew = xp + sin(th)*dl
         ypnew = yp + cos(th)*dl
         xparr(q) = xpnew
         yparr(q) = ypnew
         zparr(q) = 0
         rpnew = sqrt(xpnew**2 + ypnew**2)
         if (rpnew .lt. rplan) then !photon absorbed
            thfin = -9999. 
            distfin = -1.
            return
         endif
         aerextnew(q) = fr(rpnew, rtmp, atmp, lenalt)
         nnew = frefract(rpnew, rtmp, ntmp, lenalt)
         dndh = (nnew-n)/(rpnew-rp)
         rp = rpnew
         n = nnew      
         xp = xpnew
         yp = ypnew
         rparr(q) = rp
         call rk4(dl, val, dl, n, dndh, rplan)
         th = th0-val(4)
 03   continue
      thfin = th
      distfin = minval(rparr(1:q-1))
      call altpath(xparr, yparr,zparr,q-1,dl, alt, lenalt, rplan, paths)
      return
      end subroutine

      double precision function frefract(rp, rtmp, ntmp, lenalt)
      integer, intent(in) :: lenalt
      double precision, dimension(lenalt+4), intent(in) :: rtmp, ntmp
      double precision rp, dn, dr
      integer k
      do 04 k=1, lenalt+4
         if (rp .gt. rtmp(k)) then
            dn = ntmp(k)-ntmp(k-1)
            dr = rtmp(k)-rtmp(k-1)            
            frefract = ntmp(k) + dn/dr*(rp-rtmp(k))
            exit
         endif
 04   continue
      return
      end function
c
      double precision function fr(rp, rtmp, atmp, lenalt)
      integer, intent(in) :: lenalt
      double precision, dimension(lenalt+4):: rtmp, atmp
      double precision rp, da, dr
      integer k
      do 03 k=1, lenalt+4
         if (rp .gt. rtmp(k)) then
            da = atmp(k)-atmp(k-1)
            dr = rtmp(k)-rtmp(k-1)            
            fr = atmp(k) + da/dr*(rp-rtmp(k))
            exit
         endif
 03   continue
      return
      end function
c
      subroutine cumtrapz(yarr,xarr,start, finish, lenx,  cint)
      integer start, finish, k, lenx
      double precision, dimension(lenx) :: yarr, xarr, cint
      double precision h
      cint(1) = 0.0
      do 04 k=start+1,finish
         h = xarr(k)-xarr(k-1)
c         if(k.eq.1) then 
c            cint(k+1) = 0.5*h*(yarr(k+1)+yarr(k))
c         else 
         cint(k) = cint(k-1) + 0.5*h*(yarr(k)+yarr(k-1))
c         endif
 04   enddo
      return
      end subroutine
c
      subroutine scatcalc(xp, yp,zp, th,rtmp,ntmp,atmp, phasecumarr, 
     -     thcum,rplan, xpfin, ypfin,zpfin, thfin, scatend, alt, lenalt,
     -     paths, psi, numphot, singscat)
      integer, intent(in) :: numphot
      integer dint, scatend, darrlen, tauind, p, lenalt, q, j
      parameter (darrlen=10000)
      double precision, intent(in) :: phasecumarr(lenalt, lenth)
      double precision xp, yp, th, dl, xpfin, ypfin, thfin, xpnew, ypnew
      double precision th0, xp0, yp0, phi, eps, dndh, n, rp, beta, dmax
      double precision, dimension(darrlen) :: darr, aerextnew, taucum
      double precision, dimension(darrlen) :: xparr, yparr, zparr
      double precision nnew, val(4), rpnew, zp, zpfin, zpnew, angth
      double precision, dimension(lenalt+4):: rtmp, ntmp, atmp
      double precision dtau, taurand, tautot, taudiffrand, taudiff
      double precision dnew, dth, rplan, rand, psi, zp0, rpfin
      double precision, dimension(lenth) :: phasecum, thcum
      double precision, dimension(lenalt) :: alt
      double precision, intent(out) :: paths(lenalt)
      double precision scatalb, maxr
      double precision, intent(in) :: singscat(lenalt)
      !initialize variables as zeros
      do p=1, darrlen
         darr(p) = 0.0
         aerextnew(p) = 0.0
         taucum(p) = 0.0
         xparr(p) = 0.0
         yparr(p) = 0.0
         zparr(p) = 0.0
      enddo
      th0 = th; xp0 = xp; yp0 = yp; zp0=zp
      dl = 10.
      dmax = pathcalc(xp, yp, zp, pi/2-th0, psi, rplan, alt,lenalt)
      dmax = dmax*2. !increase dmax to ensure the entire atmosphere is traversed
      if (dmax .lt. dl) then
         xpfin = xp; ypfin=yp; thfin = th0; scatend=1; zpfin=zp
         return
      endif
      dint = dmax/dl
      do 05 p=1, dint
         darr(p) = (p-1)*dl
 05   continue
      rp = sqrt(xp**2+yp**2+zp**2)
      n = frefract(rp, rtmp, ntmp, lenalt)
      dndh = 0.0
      eps = 0.0
      phi = atan(xp/yp)
      beta = phi-th0+pi/2.0
      val(1) = rp-rplan; val(2)=phi; val(3)=beta; val(4)=eps
      maxr = maxval(alt) + rplan
      do 06 p=1, dint
         if (sqrt(xp**2+yp**2) .gt. maxr+1E-6) exit
         dcum(p) = p
         angth = sqrt(1+tan(pi/2-th)**2+tan(psi)**2)
         xpnew = xp+dl/angth
         ypnew = yp+tan(pi/2-th)*dl/angth
         zpnew = zp+tan(psi)*dl/angth
         xparr(p) = xpnew
         yparr(p) = ypnew
         zparr(p) = zpnew
         rpnew = sqrt(xpnew**2+ypnew**2+zpnew**2)
         if (rpnew .lt. rplan) then !photon absorbed
            xpfin = xp; ypfin = yp; thfin = -9999.; scatend=1; zpfin=zp
            return
         endif
         aerextnew(p) = fr(rpnew, rtmp, atmp, lenalt)
         nnew = frefract(rpnew, rtmp, ntmp, lenalt)
         dndh = (nnew-n)/(rpnew-rp)
         rp = rpnew; n = nnew; xp = xpnew; yp = ypnew; zp = zpnew
         call rk4(dl, val, dl, n, dndh, rplan)
         th = th0-val(4)
 06   continue
c      call RANDOM_SEED()
      call RANDOM_NUMBER(dtau)
      taurand = -log(dtau)
c      call cumtrapz(aerextnew*sig*1000.*dl, dcum, 1, darrlen, dint,
c     -     taucum)
      call cumtrapz(aerextnew*dl, dcum, 1, darrlen, p-1, taucum) 
      tautot = maxval(taucum)
      if (taurand .gt. tautot) then
         scatend = 1
         xpfin = xpnew; ypfin = ypnew; thfin = th; zpfin = zpnew
c         return
      else
         scatend = 0
         tauind = 0
         do 07 q=1, p-1
            if (taurand .lt. taucum(q)) then
               tauind = q
               exit
            endif
 07      end do
         taudiffrand = taurand - taucum(tauind-1)
         taudiff = taucum(tauind)-taucum(tauind-1)
c         dnew = darr(tauind-1) + dl*taudiffrand/taudiff         
         dnew = dl*taudiffrand/taudiff  
c         xpfin = xp0 + sin(th)*dnew
c         ypfin = yp0 + cos(th)*dnew
         angth = sqrt(1+tan(pi/2-th)**2+tan(psi)**2)
         xpfin = xparr(tauind-1) + dnew/angth
         ypfin = yparr(tauind-1) + dnew*tan(pi/2-th)/angth
         zpfin = zparr(tauind-1) + dnew*tan(psi)/angth
         rpfin = sqrt(xpfin**2+ypfin**2+zpfin**2)
         do 071 q=1,lenalt
            if (rpfin .gt. (alt(q)+rplan)) then
               do 072 j=1,lenth
                  phasecum(j) = phasecumarr(q,j)
                  scatalb = singscat(q)
 072           enddo
               exit
            endif
 071     enddo                    
         dth = thrand(phasecum, thcum, lenth)
         call RANDOM_NUMBER(rand)
         dth = dth*atan(cos(rand*2*pi))
         psi = psi + dth*atan(sin(rand*2*pi))
c         psi = 0.0 !!!remember to switch this back!
         th = mod(th+dth, 2*pi)
         thfin = th
      endif
      call altpath(xparr,yparr,zparr,p-1, dl, alt, lenalt, rplan, paths)
      return
      end subroutine

      real function thrand(phasecum, thcum, lenarr)
      integer lenarr
      double precision, dimension(lenarr) :: phasecum, thcum
      double precision rand, ang
c      call RANDOM_SEED()
      call RANDOM_NUMBER(rand)  
      if (rand .le. minval(phasecum)) then
         thrand = 0
      else if (rand .ge. maxval(phasecum)) then
         thrand = pi
      else
c         call RANDOM_NUMBER(rand)
         ang = interp1d(phasecum, thcum, lenarr, rand)
         thrand = ang!*cos(rand*2*pi)
c         thrand = sign(ang, rand-0.5)
      endif
      return
      end function
c     
      subroutine pathtrace(k, maxscat, lenalt, lenj, alt, minaer, refth,
     -     rtmp, ntmp, atmp, phasecumarr, thcum, thtmp, numtmp,
     -     pathtmp, drtmp, rplan, singscat, fluxscat)
      integer lenj, numscat, endscat, i, maxscat, lenalt,  scatend
      integer k, q, numphot
      double precision, dimension(lenalt) :: alt, refth, ratm, paths
      double precision, dimension(lenalt+4):: rtmp, atmp, ntmp
      double precision, dimension(lenth) :: thcum
      double precision, dimension(lenalt, lenth) :: phasecumarr
      double precision, dimension(lenj) :: thtmp, numtmp, drtmp,fluxscat
      double precision, dimension(lenalt, lenj) :: pathtmp
      double precision, dimension(lenalt) :: singscat
      double precision thinit, thnew, xp, yp, xpfin, ypfin, thfin
      double precision minaer, zp, zpfin, rplan, thref, minrref
      ratm = alt + rplan
      if (alt(k) .gt. minaer) then        
         yp = ratm(k)
         xp = dble(-sqrt(ratm(1)**2-yp**2))
         thinit = pi/2.0
         call refcalc(xp, yp, th, rtmp, ntmp, atmp, rplan, thref,
     -        minrref, alt, lenalt, paths )
         do 08 j=1, lenj
            thtmp(j) = refth(k)
            drtmp(j) = 0.0
            fluxscat(j) = 1.0
            do 081 q=1,lenalt
               pathtmp(q,j) = paths(q)
 081        enddo
 08      enddo
         return
      endif
!$ call omp_set_num_threads(16)
!$omp parallel do private(q, yp, xp, zp, thinit, numscat, endscat, psi,
!$OMP& fluxscat,  xpfin, ypfin, zpfin, thfin, paths,  pathtmp)
      do 09 j=1, lenj
         do 091 q=1, lenalt
            pathtmp(q,j) = 0.0
 091     end do
         yp = ratm(k)
         xp = dble(-sqrt(ratm(1)**2-yp**2))
         zp = 0
         thinit = pi/2.0
         numscat = 0
         endscat = 0
         psi = 0
         fluxscat(j) = 1.0
         do while ( (endscat .eq. 0) .and. (numscat<maxscat) )            
            call scatcalc(xp, yp, zp,thinit, rtmp,ntmp,atmp,phasecumarr, 
     -           thcum, rplan, xpfin, ypfin, zpfin, thfin, endscat,
     -           alt, lenalt, paths,  psi, numphot,singscat)
            fluxscat(j) = fluxscat(j)*singscat(k)
            if (endscat .eq. 0) then
               numscat = numscat+1
            endif
            xp = xpfin; yp = ypfin; thinit = thfin;
            zp = zpfin
            if (numscat .gt. maxscat) then
               scatend = 1
            endif
            do 092 q=1,lenalt
               pathtmp(q,j) = pathtmp(q,j) + paths(q)
 092        enddo
         enddo
         thtmp(j) = thinit
         drtmp(j) = psi
         numtmp(j) = numscat
 09   end do
!$omp end parallel do
      return
      end subroutine
c
      subroutine deriv(s, val, n, dndh, r0, dval)
      double precision s, n, dndh, r0, val(4), dval(4)
      dval(1) = sin(val(3))
      dval(2) = cos(val(3))/(r0+val(1)) + dndh/n*cos(val(3))
      dval(3) = cos(val(3))/(r0+val(1))
      dval(4) = dndh/n*cos(val(3))
      return 
      end subroutine
      
      subroutine rk4(xn, yn, h, n, dndh, rplan)
      double precision xn, yn(4), h, n, dndh, xnfin, ynfin(4)
      double precision k1(4), k2(4), k3(4), k4(4), rplan
      call deriv(xn, yn, n, dndh, rplan, k1)
      call deriv(xn+0.5*h, yn+0.5*k1, n, dndh, rplan, k2)
      call deriv(xn+0.5*h, yn+0.5*k2, n, dndh, rplan, k3)
      call deriv(xn+h, yn+k3, n, dndh, rplan, k4)
      ynfin = yn + h*(k1+2.0*k2+2.0*k3+k4)/6.0
      yn = ynfin
      return
      end subroutine
c
      real function interp1d(xin, yin, lenx, xtest)
      integer lenx
      double precision, dimension(lenx) :: xin, yin
      double precision xtest, ytest
      do 010 i=1,lenx-1
         if ( (xtest .lt. xin(i+1)) .and. (xtest .gt. xin(i) ) ) then
            interp1d = (yin(i+1)-yin(i))/(xin(i+1)-xin(i))*(xtest-
     -           xin(i))+yin(i)
            exit
         endif
 010  end do
      return
      end function
c
      subroutine altpath(xparr,yparr,zparr,lenxp,dl, alt, lenalt, rplan, 
     -     paths)
      integer lenxp, lenalt, i, j, alti, maxj
      double precision, dimension(lenxp) :: xparr, yparr, zparr
      double precision, dimension(lenalt) :: alt, paths
      double precision rplan, dx, dy, mag, dxnew, dynew, dl, dz, dznew
      double precision, dimension((lenxp-1)*100) :: xpnew, ypnew, zpnew,
     -     rpnew
      double precision, dimension((lenxp-1)*100) :: rpsort
      if (lenxp .le. 1) then
         do i=1, lenalt
            paths(i) = 0
         enddo
         return
      endif
      do 011 i=1, lenxp-1
         dx = xparr(i+1)-xparr(i)
         dy = yparr(i+1)-yparr(i)
         dz = zparr(i+1)-zparr(i)
         mag = sqrt(dx**2+dy**2+dz**2)
         dxnew = dx/mag
         dynew = dy/mag
         dznew = dz/mag
         do 012 j=1, 100
            xpnew((i-1)*100+j) = xparr(i) + dxnew*float(j)*dl*0.01
            ypnew((i-1)*100+j) = yparr(i) + dynew*float(j)*dl*0.01
            zpnew((i-1)*100+j) = zparr(i) + dznew*float(j)*dl*0.01
            rpnew((i-1)*100+j) = sqrt( (xpnew((i-1)*100+j))**2 
     -           + (ypnew((i-1)*100+j))**2 + (zpnew((i-1)*100+j))**2 )
 012     end do
 011  end do
      call sort(rpnew, (lenxp-1)*100, rpsort)
      j = 1
      maxj = (lenxp-1)*100
      do 013 i=1, lenalt
         alti = lenalt-i+1
         paths(alti) = 0.0
         do 014 while ( (rpsort(j) .lt. alt(alti)+rplan) .and. 
     -        (j .lt. maxj) )
            paths(alti) = paths(alti) + dl*0.01
            j=j+1
 014     enddo
 013  end do
      return
      end subroutine
c
c      END MODULE
c
c      program test
c      use mcscat
c      implicit none
      subroutine mcscatcalc(rplan, lenalt, numphot, alt, dens, nstp,
     -     bimp,dtau,aerext,singscat,thcum,phasecumarr,rstar,dist,
     -     pathtot, fluxmat, wl, lenwlarr, wlarr, u1, u2, u3, u4,
     -     heightavgout)
      integer, intent(in) :: lenalt, numphot, lenwlarr
      double precision, dimension(lenalt, numphot), intent(out) :: 
     -     fluxmat
      double precision, dimension(lenalt), intent(in) :: alt, dens, 
     -     dtau, aerext
      double precision, dimension(lenalt), intent(out) :: heightavgout
      double precision, dimension(lenalt) :: refth, refminr, heightavg,
     -     ratm, paths
      double precision, dimension(lenalt) :: singscat
      double precision, intent(in) :: nstp, bimp, rstar, dist, rplan,wl
      double precision, dimension(mxmom, lenalt) :: phmom   
      double precision, dimension(lenwlarr) :: wlarr, u1, u2, u3, u4
      double precision, dimension(lenalt) :: aerdens, ndens
      double precision,  intent(out) :: pathtot(lenalt, numphot, lenalt)
      double precision, dimension(lenalt, numphot) :: pathtmp, thfin, 
     -     drfin, numfin, fluxfin
      double precision, dimension(lenth) :: thcum
      double precision, dimension(lenalt, lenth) :: phasecumarr
      double precision, dimension(lenalt+4) :: rtmp, atmp, ntmp
      double precision, dimension(numphot) :: thtmp, numtmp, drtmp, 
     -     fluxscat
      pi = acos(-1.0)
      minaer = lenalt
      do 41 i=1, lenalt
         if ( dtau(i) .gt. 10D-6) then
            minaer = alt(i)
            exit
         endif
 41   enddo
      ratm = alt+rplan
c     
c      thlim = 0.01
c set up index of refraction profile
c      densstp = 10.0**5/kb/273.15
      densstp = 2.65D25
      ndens = 1.0+dens/densstp*(nstp-1.0)
c
      rtmp(1) = maxval(ratm)+100.0
      rtmp(2) = maxval(ratm)+0.001
      rtmp(lenalt+3) = minval(ratm)-0.001
      rtmp(lenalt+4) = minval(ratm)-100.
      atmp(1) = 0.0
      atmp(2) = 0.0
      atmp(lenalt+3) = 0.0
      atmp(lenalt+4) = 0.0  
      ntmp(1) = 1.0
      ntmp(2) = 1.0
      ntmp(lenalt+3) = 1.0
      ntmp(lenalt+4) = 1.0
      do 80 i=1, lenalt
         rtmp(i+2) = ratm(i)
         ntmp(i+2) = ndens(i)
c         atmp(i+2) = aerdens(i)
         atmp(i+2) = aerext(i) !aerosol opacity/length
 80   end do
c
      do 100 i=1, lenalt
c      do 100 i=23, 23
         yp = ratm(i)
         xp = -sqrt(ratm(1)**2-yp**2)
         th = pi/2.0
         call refcalc(xp, yp, th, rtmp, ntmp, atmp, rplan, refth(i),
     -        refminr(i), alt, lenalt, paths )
 100  enddo
c
      do 200 altind=1, lenalt
         call pathtrace(altind, maxscat, lenalt, numphot, alt, minaer,
     -     refth, rtmp, ntmp, atmp, phasecumarr, thcum, thtmp,numtmp,
     -        pathtmp, drtmp, rplan, singscat, fluxscat)
         do 210 j=1, numphot
            thfin(altind,j) = thtmp(j)
            numfin(altind,j) = numtmp(j)
            drfin(altind, j) = drtmp(j)
            fluxfin(altind, j) = fluxscat(j)
            do 211 i=1,lenalt-1         
               pathtot(altind, j, i) = pathtmp(i,j)
 211        end do
 210     end do
 200  end do
c
      starang = rstar/dist
      call transitref(wl, wlarr, u1, u2, u3, u4, lenwlarr, bimp, 
     -     rangtest, heightavgtest)
      dflux = 1.0/dble(numphot)
      do 220 i=1, lenalt
         heightavg(i) = 0.0
c         altflux = 0.0
c         altflux1 = 0.0
         do 230 j=1, numphot
c            if ( (abs(thfin(i,j)-pi/2.0) .lt. thlim) .and. 
c     -           (abs(drfin(i,j)) .lt. thlim) )then               
c               altflux = altflux + dflux
c            endif
c            if (abs(thfin(i,j)-pi/2.0) .lt. thlim) then
c               altflux1 = altflux1 + dflux
c            endif
c     calculate flux for each photon at each altitude
            rfinang = abs((thfin(i,j)-pi/2)/starang)
            if ((rfinang .lt. 4.0) .and. (rfinang .gt. -1.)) then
               do 221 k=2, lenrtest
                  if (rfinang .lt. rangtest(k) ) then
                     dr = rfinang-rangtest(k-1)
                     slope = (heightavgtest(k)-heightavgtest(k-1) )/
     -                    (rangtest(k)-rangtest(k-1) )
                     heightavg(i) = heightavg(i)+ (heightavgtest(k-1)+
     -                    slope*dr)*dflux*fluxfin(i,j)
                     fluxmat(i,j) = (heightavgtest(k-1) + slope*dr)*
     -                    fluxfin(i,j)
                     exit
                  endif
 221           end do
            endif
 230     end do
         heightavgout(i) = heightavg(i)
c         write(*,"(F10.4, F10.4, F10.4, F10.4, F12.4)") alt(i), altflux, 
c     -        altflux1, slfluxnew(i), heightavg(i) !, thfin(i,1)*180/pi
         write(*,"(F10.4, F10.4, F10.4, F10.4, F12.4)") alt(i),
     -        heightavg(i)
 220  end do
      
      return 
      end subroutine
      end module
