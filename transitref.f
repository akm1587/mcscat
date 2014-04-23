      subroutine transitref(wl, wlarr, u1, u2, u3, u4, lenwlarr, bimp, 
     -     rem, heightavg)
      implicit none
      integer, intent(in) :: lenwlarr
      integer lenth, i, lenx0, j, lenr, q
      parameter (lenth=50, lenx0 = 301, lenr=40)
      double precision, dimension(lenwlarr), intent(in) :: wlarr, u1, u2
     -     ,u3, u4
      double precision, intent(in) :: wl, bimp
      double precision sums, dx, xmax, pi, x0, th, rsq, rtmp
      double precision, dimension(lenth) :: tharr, ftmp
      double precision, dimension(lenx0) :: x0arr
      double precision, dimension(lenr), intent(out) :: rem
      double precision, dimension(lenr,lenx0) :: fluxnew
      double precision, dimension(lenr), intent(out) :: heightavg
      double precision minwl, maxwl, tflux, pflux
      pi = acos(-1.0)
      sums = sum(u1)+sum(u2)+sum(u3)+sum(u4)
      dx = 0.01
      minwl = minval(wlarr)
      maxwl = maxval(wlarr)
      if (bimp .le. 1.0) then
         xmax = -sqrt(1-bimp**2)
      else
         xmax = 0.0
      endif
      do 10 i=1,lenth
         tharr(i) = (i-1)*pi/25.
 10   enddo
      do 11 i=1, lenr
         rem(i) = (i-1)*0.1
 11   enddo
c
      do 20 i=1,lenx0
         x0arr(i) = -3+(i-1)*dx
         x0 = sqrt(x0arr(i)**2+bimp**2)
         do 21 q=1, lenr
c         do 21 q=1, 1
            do 22 j=1, lenth
               th = tharr(j)
               rsq = x0**2+bimp**2+rem(q)**2+2*rem(q)*(x0*cos(th)+
     -           bimp*sin(th))
               if (rsq .lt. 0) then
                  rsq = 0.0
               endif
               rtmp = sqrt(rsq)
               ftmp(j) = limb(rtmp)
 22         enddo
            fluxnew(q, i) = sum(ftmp)/max(1,size(ftmp))
 21      enddo
 20   enddo
      do 23 i=1, lenr       
         tflux = sum(fluxnew(i,202:301))/max(1,size(fluxnew(i,202:301)))
         pflux = sum(fluxnew(i,102:201))/max(1,size(fluxnew(i,102:201)))
c         heightavg(i) = -pflux+tflux !account for different difference in transmitted flux because of refracted flux prior to ingress
         heightavg(i) = tflux
 23   enddo
      return
c    
      contains
      double precision function limb(rang)
      double precision, intent(in) ::  rang
      integer p, wlind
      double precision  mu, dwl, d1, d2, d3, d4
      double precision fu1, fu2, fu3, fu4
      if (abs(rang) .gt. 1.0) then
         limb = 0.0
         return
      else
         if (wl .le. minwl) then
            fu1 = 1.0; fu2=0.0; fu3=0.0; fu4=0.0
         else if (wl .ge. maxwl) then
            fu1 = 0.0; fu2=0.0; fu3=0.0; fu4=0.0
         else
            do 30 p=1, lenwlarr-1
               if (wl .gt. wlarr(p) ) then
                  wlind = p
                  exit
               endif
 30         enddo
            dwl = wlarr(wlind+1)-wlarr(wlind)
            d1 = u1(wlind+1)-u1(wlind)
            d2 = u2(wlind+1)-u2(wlind)
            d3 = u3(wlind+1)-u3(wlind)
            d4 = u4(wlind+1)-u4(wlind)
            fu1 = d1/dwl*(wl-wlarr(wlind)) + u1(wlind)
            fu2 = d2/dwl*(wl-wlarr(wlind)) + u2(wlind)
            fu3 = d3/dwl*(wl-wlarr(wlind)) + u3(wlind)
            fu4 = d4/dwl*(wl-wlarr(wlind)) + u4(wlind)
         endif
         mu = sqrt(1-rang**2)
         limb = 1-fu1*(1-mu**.5)-fu2*(1-mu)-fu3*(1-mu**1.5)-
     -        fu4*(1-mu**2)
       endif
      return
      end function
c
      end
