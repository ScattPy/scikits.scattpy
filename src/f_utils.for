
      SUBROUTINE MAT_A0(NG,NN,Rad,Ang,coef,w,O)
      INTEGER    NG,NN
      COMPLEX*16 Rad(NG,NN)
      COMPLEX*16 Ang(NG,NN),coef(NG)
      COMPLEX*16 O(NN,NN)
      REAL*8     w(NG)
cf2py intent(in) Rad,Ang,coef,w
cf2py intent(out) O
cf2py intent(hide) NG,NN
      INTEGER    l,n,k

      DO l=1,NN
      DO n=1,NN
        O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
      DO l=1,NN
      DO n=1,NN
        O(l,n) = O(l,n) + Rad(k,l)*Ang(k,l)*Ang(k,n)*coef(k)*w(k)
      ENDDO
      ENDDO
      ENDDO
      END

      SUBROUTINE MAT_B(NG,NN,ki,Rad,Radd,Ang,Angd,R,Rdr,Sint,w,O)
      INTEGER    NG,NN
      COMPLEX*16 ki
      COMPLEX*16 Rad(NG,NN),Radd(NG,NN)
      COMPLEX*16 Ang(NG,NN),Angd(NG,NN)
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 R(NG),Rdr(NG),Sint(NG)
      REAL*8     w(NG)
cf2py intent(in) ki,Rad,Radd,Ang,Angd,R,Rdr,Sint,w
cf2py intent(out) O
cf2py intent(hide) NG,NN
      INTEGER    l,n,k

      DO l=1,NN
      DO n=1,NN
        O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
      DO l=1,NN
      DO n=1,NN
        O(l,n) = O(l,n) + ( ki*R(k)*Radd(k,l)*Ang(k,l)
     &                         +Rdr(k)*Sint(k)*Rad(k,l)*Angd(k,l) )
     &                      * Ang(k,n)*Sint(k)*w(k)
      ENDDO
      ENDDO
      ENDDO
      END

      SUBROUTINE MAT_D0(NG,NN,ki,Rad,Radd,Ang,Angd,R,Sint,Cost,coef,w,O)
      INTEGER    NG,NN
      COMPLEX*16 ki
      COMPLEX*16 Rad(NG,NN),Radd(NG,NN)
      COMPLEX*16 Ang(NG,NN),Angd(NG,NN),coef(NG)
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 R(NG),Sint(NG),Cost(NG)
      REAL*8     w(NG)
cf2py intent(in) ki,Rad,Radd,Ang,Angd,R,Sint,Cost,coef,w
cf2py intent(out) O
cf2py intent(hide) NG,NN
      INTEGER    l,n,k

      DO l=1,NN
      DO n=1,NN
        O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
      DO l=1,NN
      DO n=1,NN
        O(l,n) = O(l,n) + ( ki*R(k)*Cost(k)*Radd(k,l)*Ang(k,l)
     &                     +Sint(k)**2*Rad(k,l)*Angd(k,l) )
     &                  * Ang(k,n)*coef(k)*w(k)
      ENDDO
      ENDDO
      ENDDO
      END

      SUBROUTINE MAT_E0(NG,NN,ki,Rad,Radd,Ang,R,coef,w,O)
      INTEGER    NG,NN
      COMPLEX*16 ki
      COMPLEX*16 Rad(NG,NN),Radd(NG,NN)
      COMPLEX*16 Ang(NG,NN),coef(NG)
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 R(NG)
      REAL*8     w(NG)
cf2py intent(in) ki,Rad,Radd,Ang,R,coef,w
cf2py intent(out) O
cf2py intent(hide) NG,NN
      INTEGER    l,n,k

      DO l=1,NN
      DO n=1,NN
        O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
      DO l=1,NN
      DO n=1,NN
        O(l,n) = O(l,n) + ( ki*R(k)*Radd(k,l) +Rad(k,l) )
     &                  * Ang(k,l) * Ang(k,n)*coef(k)*w(k)
      ENDDO
      ENDDO
      ENDDO
      END

      SUBROUTINE MAT_G0(NG,NN,ki,Rad,Radd,Ang,Angd,R,Rd,Sint,coef,w,O)
      INTEGER    NG,NN
      COMPLEX*16 ki
      COMPLEX*16 Rad(NG,NN),Radd(NG,NN)
      COMPLEX*16 Ang(NG,NN),Angd(NG,NN),coef(NG)
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 R(NG),Rd(NG),Sint(NG)
      REAL*8     w(NG)
cf2py intent(in) ki,Rad,Radd,Ang,Angd,R,Rd,Sint,coef,w
cf2py intent(out) O
cf2py intent(hide) NG,NN
      INTEGER    l,n,k

      DO l=1,NN
      DO n=1,NN
        O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
      DO l=1,NN
      DO n=1,NN
        O(l,n) = O(l,n) + ( ki*Rd(k)*Radd(k,l)*Ang(k,l)
     &                         -Sint(k)*Rad(k,l)*Angd(k,l) )
     &                      * Ang(k,n)*coef(k)*w(k)
      ENDDO
      ENDDO
      ENDDO
      END


      SUBROUTINE f1(NG,R,Rd,Rdd,Sint,Cost,f)
      implicit real*8 (r)
      INTEGER NG
      REAL*8 R(NG),Rd(NG),Sint(NG),Cost(NG),f(NG)
cf2py intent(in) R,Rd,Sint,Cost
cf2py intent(out) f
cf2py intent(hide) NG
      INTEGER K
      
      DO k=1,NG
        rs = R(k)
        rsd = Rd(k)
        rsdd = Rdd(k)
        sints = Sint(k)
        costs = Cost(k)

        rdrs = rsd/rs
        r2r2 = rs**2+rsd**2
        r2r22 = r2r2**2
        rsin = rs *sints
        rcos = rs *costs
        rdsin= rsd*sints
        rdcos= rsd*costs
        rdcosrsin = rdcos-rsin
        rdsinrcos = rdsin+rcos

        f(k) = (rs**2 - rs*rsdd + 2*(rsd**2)) 
     &        *( rs *rdcosrsin + rsd*rdsinrcos) 
     &        /r2r22
        ENDDO
      end

      SUBROUTINE f2(NG,R,Rd,Rdd,f)
      implicit real*8 (r)
      INTEGER NG
      REAL*8 R(NG),Rd(NG),Sint(NG),Cost(NG),f(NG)
cf2py intent(in) R,Rd,Sint,Cost
cf2py intent(out) f
cf2py intent(hide) NG
      INTEGER K
      
      DO k=1,NG
        rs = R(k)
        rsd = Rd(k)
        rsdd = Rdd(k)

        rdrs = rsd/rs
        r2r2 = rs**2+rsd**2
        r2r22 = r2r2**2

        f(k) = (rs**4 - 2*(rs**3)*rsdd 
     &          + 2*(rs**2)*(rsd**2) - rsd**4)
     &         /r2r22
        ENDDO
      end

      SUBROUTINE f3(NG,R,Rd,Rdd,Sint,Cost,f)
      implicit real*8 (r)
      INTEGER NG
      REAL*8 R(NG),Rd(NG),Sint(NG),Cost(NG),f(NG)
cf2py intent(in) R,Rd,Sint,Cost
cf2py intent(out) f
cf2py intent(hide) NG
      INTEGER K
      
      DO k=1,NG
        rs = R(k)
        rsd = Rd(k)
        rsdd = Rdd(k)
        sints = Sint(k)
        costs = Cost(k)

        rdrs = rsd/rs
        r2r2 = rs**2+rsd**2
        r2r22 = r2r2**2
        rsin = rs *sints
        rcos = rs *costs
        rdsin= rsd*sints
        rdcos= rsd*costs
        rdcosrsin = rdcos-rsin
        rdsinrcos = rdsin+rcos

        f(k) = 2*(rs**2 - rs*rsdd + 2*(rsd**2)) 
     &          *rdcosrsin*rdsinrcos 
     &          /r2r22
        ENDDO
      end

      SUBROUTINE f4(NG,R,Rd,Rdd,Sint,Cost,f)
      implicit real*8 (r)
      INTEGER NG
      REAL*8 R(NG),Rd(NG),Sint(NG),Cost(NG),f(NG)
cf2py intent(in) R,Rd,Sint,Cost
cf2py intent(out) f
cf2py intent(hide) NG
      INTEGER K
      
      DO k=1,NG
        rs = R(k)
        rsd = Rd(k)
        rsdd = Rdd(k)
        sints = Sint(k)
        costs = Cost(k)

        rdrs = rsd/rs
        r2r2 = rs**2+rsd**2
        r2r22 = r2r2**2
        rsin = rs *sints
        rcos = rs *costs
        rdsin= rsd*sints
        rdcos= rsd*costs
        rdcosrsin = rdcos-rsin
        rdsinrcos = rdsin+rcos

        f(k) = ( ((rs**3)*rsdd+(rs**2)*(rsd**2) 
     &            -rs*(rsd**2)*rsdd+3*(rsd**4)) *sints 
     &          +rdrs*costs
     &           *(rs**4 - 2*(rs**3)*rsdd 
     &             +2*(rs**2)*(rsd**2) - rsd**4)) 
     &         /r2r22
        ENDDO
      end

c---------- EBCM ------------------------------------

      SUBROUTINE axitm (NG,NN,Rad1,Radd1,Rad2,Radd2,Ang,Angd,
     &                        Rs,Rsd,SinT,CosT,k1,k2,ep12,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 Rad1(NG,NN),Radd1(NG,NN),Rad2(NG,NN),Radd2(NG,NN)
      COMPLEX*16 Ang(NG,NN), Angd(NG,NN)
      COMPLEX*16 k1,k2,ep12
      REAL*8     SinT(NG),CosT(NG), Rs(NG),Rsd(NG)
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      integer l,n,k
      complex*16 coef, j1,jd1,j2,jd2,pl,pdl,pn,pdn
      real*8 sn,cs,r,rd
      complex*16 fun1
      
      DO l=1,NN
      DO n=1,NN
        O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
      r  = Rs(k)
      rd = Rsd(k)
      cs = CosT(k)
      sn = SinT(k)

      DO l=1,NN
      DO n=1,NN

      j1  = Rad1(k,l)
      jd1 = Radd1(k,l)
      j2  = Rad2(k,n)
      jd2 = Radd2(k,n)
      pl  = Ang(k,l)
      pdl = Angd(k,l)
      pn  = Ang(k,n)
      pdn = Angd(k,n)

      fun1 = k1**2 * r**2 * ( jd1 * j2 -
     & ep12 * k2/k1 * j1 * jd2) *
     & pl * pn * sn
     & + k1 * rd * sn**2 * (pdl * pn - ep12*
     & pl * pdn) * j1 * j2
     & - (ep12 - 1) * (k1 * r * sn - k1 * rd * cs) *
     & j1 * j2 * pl * pn

      O(l,n) = O(l,n) + fun1*w(k)
      
      ENDDO
      ENDDO
      ENDDO

      DO l=1,NN
       coef = (0d0,1d0) * (2d0 * l + 1) / (2 * l * (l+1))
      DO n=1,NN
       O(l,n) = O(l,n) * coef
      ENDDO
      ENDDO

      end

      SUBROUTINE axite (NG,NN,Rad1,Radd1,Rad2,Radd2,Ang,Angd,
     &                        Rs,Rsd,SinT,CosT,k1,k2,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 Rad1(NG,NN),Radd1(NG,NN),Rad2(NG,NN),Radd2(NG,NN)
      COMPLEX*16 Ang(NG,NN), Angd(NG,NN)
      COMPLEX*16 k1,k2
      REAL*8     SinT(NG),CosT(NG), Rs(NG),Rsd(NG)
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      integer l,n,k
      complex*16 coef, j1,jd1,j2,jd2,pl,pdl,pn,pdn
      real*8 sn,cs,r,rd
      complex*16 fun1
      
      DO l=1,NN
      DO n=1,NN
        O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
      r  = Rs(k)
      rd = Rsd(k)
      cs = CosT(k)
      sn = SinT(k)

      DO l=1,NN
      DO n=1,NN

      j1  = Rad1(k,l)
      jd1 = Radd1(k,l)
      j2  = Rad2(k,n)
      jd2 = Radd2(k,n)
      pl  = Ang(k,l)
      pdl = Angd(k,l)
      pn  = Ang(k,n)
      pdn = Angd(k,n)

      fun1 = k1**2 * r**2 * ( jd1 * j2 -
     & k2/k1 * j1 * jd2) *
     & pl * pn * sn
     & + k1 * rd * sn**2 * (pdl * pn - 
     & pl * pdn) * j1 * j2

      O(l,n) = O(l,n) + fun1*w(k)
      
      ENDDO
      ENDDO
      ENDDO

      DO l=1,NN
       coef = (0d0,1d0) * (2d0 * l + 1) / (2 * l * (l+1))
      DO n=1,NN
       O(l,n) = O(l,n) * coef
      ENDDO
      ENDDO

      end

      SUBROUTINE naxitm (NG,NN,m,Rad1,Radd1,Rad2,Radd2,Ang,Angd,
     &                        Rs,Rsd,SinT,CosT,k1,k2,ep12,ep21,w,O)
      implicit none
      INTEGER NG,NN,m
      COMPLEX*16 Rad1(NG,NN),Radd1(NG,NN),Rad2(NG,NN),Radd2(NG,NN)
      COMPLEX*16 Ang(NG,NN), Angd(NG,NN)
      COMPLEX*16 k1,k2,ep12,ep21
      REAL*8     SinT(NG),CosT(NG), Rs(NG),Rsd(NG)
      COMPLEX*16 O(2*NN,2*NN)
      COMPLEX*16 w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      integer l,n,k,i,l0,n0
      complex*16 coef, j1,jd1,j2,jd2,pl,pdl,pn,pdn
      real*8 sn,cs,r,rd
      integer alm
      complex*16 za1,zb1,za2,zb2
      
      DO l=1,2*NN
      DO n=1,2*NN
        O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
      r  = Rs(k)
      rd = Rsd(k)
      cs = CosT(k)
      sn = SinT(k)

      DO l0=1,NN
      DO n0=1,NN
      l=l0+m-1
      n=n0+m-1

      j1  = Rad1(k,l0)
      jd1 = Radd1(k,l0)
      j2  = Rad2(k,n0)
      jd2 = Radd2(k,n0)
      pl  = Ang(k,l0)
      pdl = Angd(k,l0)
      pn  = Ang(k,n0)
      pdn = Angd(k,n0)

      za1 = k1**2 * r**2 * ( j2 * jd1 -
     & sqrt(ep21) * jd2 * j1) *
     & pn * pl * sn
     & + k1 * rd * j2* j1 * (pn * pdl -
     & pl * pdn) * sn**2
     & - (ep12 - 1) * (sqrt(ep21)* k1**2
     & * r * rd * cs * jd2 * pn +
     & k1 * rd * sn**2 * j2 * pdn) * j1 * pl

      zb1 = (ep12 - 1) * (sqrt(ep21) *
     & k1**2 * r**2 * rd * jd2 * pn +
     & k1 * r * rd * j2 * pn) *
     & j1 * pl
c      need *k(1)
c      zb1 = (ep12 - 1) * (sqrt(ep21) *
c     & k1**3 * r**2 * rd * jd2 * pn +
c     & k1**2 * r * rd * j2 * pn) *
c     & j1 * pl

      za2 = (ep12 - 1) * k1**2 *(rd * cs - r * sn) *
     & (sqrt(ep21) * cs * jd2 * pn + sn**2 / (k1 * r) *
     & j2 * pdn) * j1 * pl
c      need /k(1)
c      za2 = (ep12 - 1) * (k1 * rd * cs - k1 * r * sn) *
c     & (sqrt(ep21) * cs * jd2 * pn + sn**2 / (k1 * r) *
c     & j2 * pdn) * j1 * pl

      zb2 = k1**2 * r**2 * ( j2 * jd1 -
     & sqrt(ep21) * jd2 * j1) *
     & pn * pl * sn
     & + k1 * rd * j2 * j1 * (pn * pdl -
     & pl * pdn) * sn**2
     & + (ep12 - 1) * (k1 * rd * cs - k1 * r * sn) *
     & (sqrt(ep21)* k1 * r * jd2 * pn +
     & j2 * pn) * j1 * pl

      O(l0   ,n0   ) = O(l0   ,n0   ) + za1*w(k)
      O(l0   ,n0+NN) = O(l0   ,n0+NN) + zb1*w(k)
      O(l0+NN,n0   ) = O(l0+NN,n0   ) + za2*w(k)
      O(l0+NN,n0+NN) = O(l0+NN,n0+NN) + zb2*w(k)
      
      ENDDO
      ENDDO
      ENDDO

      DO l0=1,NN
       l=l0+m-1
       alm=1d0
       DO i =1,2*m
         alm = alm*(l-m+i)
       ENDDO
       coef = (0d0,1d0) * (2d0 * l + 1) / (2 * alm)
       DO n0=1,NN
         O(l0   ,n0   ) = O(l0   ,n0   ) * coef
         O(l0   ,n0+NN) =-O(l0   ,n0+NN) * coef
         O(l0+NN,n0   ) = O(l0+NN,n0   ) * coef
         O(l0+NN,n0+NN) = O(l0+NN,n0+NN) * coef
       ENDDO
      ENDDO

      end

      SUBROUTINE naxite (NG,NN,m,Rad1,Radd1,Rad2,Radd2,Ang,Angd,
     &                        Rs,Rsd,Rsdd,SinT,CosT,k1,k2,ep12,ep21,w,O)
      implicit none
      INTEGER NG,NN,m
      COMPLEX*16 Rad1(NG,NN),Radd1(NG,NN),Rad2(NG,NN),Radd2(NG,NN)
      COMPLEX*16 Ang(NG,NN), Angd(NG,NN)
      COMPLEX*16 k1,k2,ep12,ep21
      REAL*8     SinT(NG),CosT(NG), Rs(NG),Rsd(NG),Rsdd(NG)
      COMPLEX*16 O(2*NN,2*NN)
      COMPLEX*16 w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      integer l,n,k,i,l0,n0
      complex*16 coef, j1,jd1,j2,jd2,pl,pdl,pn,pdn
      real*8 sn,cs,r,rd,rdd
      real*8  f1,f2,f3,f4
      complex*16 za1,zb1,za2,zb2
      integer alm
      
      DO l=1,2*NN
      DO n=1,2*NN
        O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
      r  = Rs(k)
      rd = Rsd(k)
      rdd= Rsdd(k)
      cs = CosT(k)
      sn = SinT(k)

      f1 = ((r**2 + rd**2) + (rd**2 - r * rdd)) / (r**2 + rd**2)**2
     & * ( r * (rd*cs - r*sn) + rd * (rd*sn + r*cs))
      f2 = ( (r**2 + rd**2) * (r**2 - rd**2) + 2 * r**2 * (rd**2 -
     & r*rdd) ) / (r**2 + rd**2)**2
      f3 = ( 2 * (rd*cs - r*sn) * (rd*sn + r*cs) *
     & ( (r**2 + rd**2) + (rd**2 - r * rdd) ) ) / (r**2 + rd**2)**2
      f4 = ( (r**3 * rdd + r**2 * rd**2 - r * rd**2 * rdd + 3 * rd**4)
     & * sn + rd/r * (r**4 - 2 * r**3 * rdd + 2 * r**2 * rd**2 - rd**4)
     & * cs ) / (r**2 + rd**2)**2

      DO l=1,NN
      DO n=1,NN

      j1  = Rad1(k,l)
      jd1 = Radd1(k,l)
      j2  = Rad2(k,n)
      jd2 = Radd2(k,n)
      pl  = Ang(k,l)
      pdl = Angd(k,l)
      pn  = Ang(k,n)
      pdn = Angd(k,n)

      za1 =
     & k1**2 * r**2 * ( j2 * jd1 - sqrt(ep21) * jd2 * j1) *
     & pn * pl * sn
     & + k1 * rd * j2 * j1 * (pn * pdl -
     & pl * pdn) * sn**2  -
     & (ep21 - 1) * (
     & (k1 * r**2 * (rd*cs - r*sn) ) / (r**2 + rd**2) *
     & j2 * pn * ( k1 * r * jd1 * pl +
     & rd/r * j1 * pdl * sn ) -
     & ( k1 * r * rd * (rd*cs - r*sn) ) / (r**2 + rd**2) *
     & ( sqrt(ep21) * k1 * rd * jd2 * pn -
     & j2 * pdn * sn ) * j1 * pl +
     & k1 * r * f1 * j2 * pn * j1 * pl )

      zb1 =
     & - (ep21 - 1) * (
     & (k1**2 * r**3 * rd ) / (r**2 + rd**2) *
     & j2 * pn * ( k1 * r * jd1 * pl +
     & rd/r * j1 * pdl * sn ) -
     & ( k1**2 * r**2 * rd**2 ) / (r**2 + rd**2) *
     & ( sqrt(ep21) * k1 * rd * jd2 * pn -
     & j2 * pdn * sn) * j1 * pl +
     & k1**2 * r * rd * f2 * j2 * pn * j1 * pl )
     & /k1

      za2 =
     & (ep21 - 1) * (
     & ( (rd*cs - r*sn) * (rd*sn + r*cs) ) / (r**2 + rd**2) *
     & j2 * pn * ( k1 * r * jd1 * pl +
     & rd/r * j1 * pdl * sn ) -
     & (rd*cs - r*sn)**2 / (r**2 + rd**2) *
     & ( sqrt(ep21) * k1 * rd * jd2 * pn -
     & j2 * pdn * sn) * j1 * pl +
     & f3 * j2 * pn * j1 * pl )
     & *k1

      zb2 =
     & k1**2 * r**2 * ( j2 * jd1 - sqrt(ep21) * jd2 * j1) *
     & pn * pl * sn
     & + k1 * rd * j2 * j1 * (pn * pdl -
     & pl * pdn) * sn**2 +
     & (ep21 - 1) * (
     & (k1 * r * rd * (rd*sn + r*cs) ) / (r**2 + rd**2) *
     & j2 * pn * ( k1 * r * jd1 * pl +
     & rd/r * j1 * pdl * sn ) -
     & ( k1 * r * rd * (rd*cs - r*sn) ) / (r**2 + rd**2) *
     & ( sqrt(ep21) * k1 * rd * jd2 * pn -
     & j2 * pdn * sn ) * j1 * pl +
     & k1 * r * f4 * j2 * pn * j1 * pl )

      O(l   ,n   ) = O(l   ,n   ) + za1*w(k)
      O(l   ,n+NN) = O(l   ,n+NN) + zb1*w(k)
      O(l+NN,n   ) = O(l+NN,n   ) + za2*w(k)
      O(l+NN,n+NN) = O(l+NN,n+NN) + zb2*w(k)
      
      ENDDO
      ENDDO
      ENDDO

      DO l0=1,NN
       l=l0+m-1
       alm=1d0
       DO i =1,2*m
         alm = alm*(l-m+i)
       ENDDO
       coef = (0d0,1d0) * (2d0 * l + 1) / (2 * alm)
       DO n0=1,NN
         O(l0   ,n0   ) = O(l0   ,n0   ) * coef
         O(l0   ,n0+NN) = O(l0   ,n0+NN) * coef
         O(l0+NN,n0   ) = O(l0+NN,n0   ) * coef
         O(l0+NN,n0+NN) = O(l0+NN,n0+NN) * coef
       ENDDO
      ENDDO
      end

c------------- PMM --------------------------------


      SUBROUTINE pmmaxitm (NG,NN,Radj1,Raddj1,Radj2,Raddj2,Radh1,Raddh1,
     &                        Ang,Angd,RadjR,RadhR,
     &                        R,Rd,SinT,CtgT,k1,k2,ep12,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 Radj1(NG,NN),Raddj1(NG,NN),Radj2(NG,NN),Raddj2(NG,NN)
      COMPLEX*16 Radh1(NG,NN),Raddh1(NG,NN),RadjR(NN),RadhR(NN)
      COMPLEX*16 Ang(NG,NN), Angd(NG,NN)
      COMPLEX*16 k1,k2,ep12
      REAL*8     SinT(NG),CtgT(NG), R(NG),Rd(NG)
      COMPLEX*16 O(2*NN,3*NN)
      COMPLEX*16 w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      integer l,n,k
      complex*16 jl1,hl1,jl2,jl1d,hl1d,jl2d,jn1,hn1,jn2,jn1d,hn1d,jn2d
      complex*16 p1l,p1ld,p1n,p1nd, hl1R,jl2R,hn1R,jn2R
      real*8 sints,ctgts,rs,rsd,rsrsd
      complex*16 A1l,A2l,A3l,A1n,A2n,A3n,B1l,B2l,B3l,B1n,B2n,B3n,E2l,E2n
      
      DO l=1,2*NN
      DO n=1,3*NN
        O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
      rs = r (k)
      rsd = rd (k)
c      ts =  athetas (k)
      sints =  SinT (k)
      ctgts =  CtgT (k)

      DO l=1,NN
      DO n=1,NN


          jl1 = Radj1 (k, l)
          hl1 = Radh1 (k, l)
          jl2 = Radj2 (k, l)
          jl1d = Raddj1 (k, l)
          hl1d = Raddh1 (k, l)
          jl2d = Raddj2 (k, l)

          jn1 = Radj1 (k, n)
          hn1 = Radh1 (k, n)
          jn2 = Radj2 (k, n)
          jn1d = Raddj1 (k, n)
          hn1d = Raddh1 (k, n)
          jn2d = Raddj2 (k, n)

          p1l  = Ang  (k, l)
          p1ld = Angd (k, l)
          p1n  = Ang  (k, n)
          p1nd = Angd (k, n)

          hl1R = RadhR (l)
          jl2R = RadjR (l)
          hn1R = RadhR (n)
          jn2R = RadjR (n)
c          hl1R = 1d0
c          jl2R = 1d0
c          hn1R = 1d0
c          jn2R = 1d0

          A1l = jl1 * p1l
          A3l = hl1 * p1l / hl1R
          A2l = jl2 * p1l / jl2R

          A1n = jn1 * p1n
          A3n = hn1 * p1n / hn1R
          A2n = jn2 * p1n / jn2R

          B1l = k1 * rs * jl1d * p1l +
     &          rsd / rs * sints * jl1 * p1ld
          B3l = (k1 * rs * hl1d * p1l +
     &          rsd / rs * sints * hl1 * p1ld) / hl1R
          B2l = (k2 * rs * jl2d * p1l +
     &          rsd / rs * sints * jl2 * p1ld) / jl2R
          E2l = ep12 * B2l - (1 - ep12) *
     &          (1 - rsd / rs * ctgts) * A2l

          B1n = k1 * rs * jn1d * p1n +
     &          rsd / rs * sints * jn1 * p1nd
          B3n = (k1 * rs * hn1d * p1n +
     &          rsd / rs * sints * hn1 * p1nd)  / hn1R
          B2n = (k2 * rs * jn2d * p1n +
     &          rsd / rs * sints * jn2 * p1nd)  / jn2R
          E2n = ep12 * B2n - (1 - ep12) *
     &          (1 - rsd / rs * ctgts) * A2n

          rsrsd = 1 / (rs**2 + rsd**2)
          O(l   ,n   ) = O(l   ,n   ) + (dconjg( A3l ) * A3n +
     &              rsrsd * dconjg( B3l ) * B3n) * w(k)
          O(l   ,n+NN) = O(l   ,n+NN) - (dconjg( A3l ) * A2n +
     &               rsrsd * dconjg( B3l ) * E2n) * w(k)
          O(l+NN,n   ) = O(l+NN,n   ) - (dconjg( A2l ) * A3n +
     &               rsrsd * dconjg( E2l ) * B3n) * w(k)
          O(l+NN,n+NN) = O(l+NN,n+NN) + (dconjg( A2l ) * A2n +
     &              rsrsd * dconjg( E2l ) * E2n) * w(k)

          O(l   ,n+2*NN) = O(l   ,n+2*NN) + (dconjg( A3l ) * A1n +
     &              rsrsd * dconjg( B3l ) * B1n) * w(k)
          O(l+NN,n+2*NN) = O(l+NN,n+2*NN) - (dconjg( A2l ) * A1n +
     &               rsrsd * dconjg( E2l ) * B1n) * w(k)

      ENDDO
      ENDDO
      ENDDO
      end

      SUBROUTINE pmmnaxitm(NG,NN,Radj1,Raddj1,Radj2,Raddj2,Radh1,Raddh1,
     &                        Ang,Angd,RadjR,RadhR,
     &                        R,Rd,SinT,CosT,k1,k2,ep12,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 Radj1(NG,NN),Raddj1(NG,NN),Radj2(NG,NN),Raddj2(NG,NN)
      COMPLEX*16 Radh1(NG,NN),Raddh1(NG,NN),RadjR(NN),RadhR(NN)
      COMPLEX*16 Ang(NG,NN), Angd(NG,NN)
      COMPLEX*16 k1,k2,ep12
      REAL*8     SinT(NG),CosT(NG), R(NG),Rd(NG)
      COMPLEX*16 O(4*NN,6*NN)
      COMPLEX*16 w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      integer l,n,k
      complex*16 jl1,hl1,jl2,jl1d,hl1d,jl2d,jn1,hn1,jn2,jn1d,hn1d,jn2d
      complex*16 pml,pmld,pmn,pmnd, hl1R,jl2R,hn1R,jn2R
      real*8 sints,costs,rs,rsd,rsrsd
      complex*16 A1l,A2l,A3l,A1n,A2n,A3n,B1l,B2l,B3l,B1n,B2n,B3n
      complex*16 C2l,D2l,F2l,C2n,D2n,F2n
      complex*16 rsrsd1,rsrsd2
      
      DO l=1,4*NN
      DO n=1,6*NN
        O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
      rs = r (k)
      rsd = rd (k)
c      ts =  athetas (k)
      sints =  SinT (k)
      costs =  CosT (k)
      rsrsd = 1 / (rs**2 + rsd**2)
      rsrsd1 = (ep12 - 1) * rsd / (rs * sints)
      rsrsd2 = (ep12 - 1) * (rsd * costs - rs * sints) /
     &         (rs**2 * sints)

      DO l=1,NN
      DO n=1,NN


          jl1 = Radj1 (k, l)
          hl1 = Radh1 (k, l)
          jl2 = Radj2 (k, l)
          jl1d = Raddj1 (k, l)
          hl1d = Raddh1 (k, l)
          jl2d = Raddj2 (k, l)

          jn1 = Radj1 (k, n)
          hn1 = Radh1 (k, n)
          jn2 = Radj2 (k, n)
          jn1d = Raddj1 (k, n)
          hn1d = Raddh1 (k, n)
          jn2d = Raddj2 (k, n)

          pml  = Ang  (k, l)
          pmld = Angd (k, l)
          pmn  = Ang  (k, n)
          pmnd = Angd (k, n)

          hl1R = RadhR (l)
          jl2R = RadjR (l)
          hn1R = RadhR (n)
          jn2R = RadjR (n)
c          hl1R = 1d0
c          jl2R = 1d0
c          hn1R = 1d0
c          jn2R = 1d0


           A1l = jl1 * pml
           A3l = hl1 * pml / hl1R
           A2l = jl2 * pml / jl2R
           B1l = k1 * rs * jl1d * pml + rsd / rs * sints * jl1 * pmld
           B3l = (k1 * rs * hl1d * pml + rsd / rs * sints * hl1 * pmld)
     &           / hl1R
           B2l = (k2 * rs * jl2d * pml + rsd / rs * sints * jl2 * pmld)
     &           / jl2R
           C2l = (k2 * rs * costs * jl2d * pml + sints**2 * jl2 * pmld)
     &           / jl2R
           D2l = k1 * rs * ( k2*rs * jl2d + jl2 ) * pml / jl2R
           F2l = B2l + rsrsd1 * C2l

           A1n = jn1 * pmn
           A3n = hn1 * pmn / hn1R
           A2n = jn2 * pmn / jn2R
           B1n = k1 * rs * jn1d * pmn + rsd / rs * sints * jn1 * pmnd
           B3n = (k1 * rs * hn1d * pmn + rsd / rs * sints * hn1 * pmnd)
     &           / hn1R
           B2n = (k2 * rs * jn2d * pmn + rsd / rs * sints * jn2 * pmnd)
     &           / jn2R
           C2n = (k2 * rs * costs * jn2d * pmn + sints**2 * jn2 * pmnd)
     &           / jn2R
           D2n = k1 * rs * ( k2*rs * jn2d + jn2 ) * pmn / jn2R
           F2n = B2n + rsrsd1 * C2n

          O(l,n) = O(l,n)
     &             + (conjg( A3l ) * A3n +
     &             conjg( B3l ) * B3n * rsrsd) * w(k)

          O(l,n+2*NN) = O(l,n+2*NN) -
     &             (conjg( A3l ) * A2n +
     &             conjg( B3l ) * F2n * rsrsd) * w(k)

          O(l,n+3*NN) = O(l,n+3*NN) -
     &             (conjg (B3l) * rsrsd1 * D2n *
     &             rsrsd) * w(k)

          O(l+2*NN,n) = O(l+2*NN,n) -
     &             (conjg( A2l ) * A3n +
     &             conjg( F2l ) * B3n * rsrsd) * w(k)

          O(l+2*NN,n+2*NN) = O(l+2*NN,n+2*NN) +
     &             (conjg( A2l ) * A2n +
     &             conjg( F2l ) * F2n * rsrsd +
     &             rsrsd * ((rsd * costs - rs * sints) / (rs**2 *
     &             sints))**2 * (cdabs(ep12 - 1))**2 *
     &             conjg( C2l ) * C2n) * w(k)

          O(l+2*NN,n+1*NN) = O(l+2*NN,n+1*NN) +
     &             (conjg( rsrsd2 * C2l) * B3n *
     &             rsrsd) * w(k)

          O(l+2*NN,n+3*NN) = O(l+2*NN,n+3*NN) +
     &             (conjg( F2l ) * rsrsd1 * D2n * rsrsd -
     &             conjg( rsrsd2 * C2l ) * (B2n - rsrsd2 * D2n ) *
     &             rsrsd) * w(k)

          O(l+NN,n+2*NN) = O(l+NN,n+2*NN) +
     &             (conjg( B3l ) * rsrsd2 * C2n *
     &             rsrsd) * w(k)

          O(l+NN,n+1*NN) = O(l+NN,n+1*NN) +
     &             (conjg( A3l ) * A3n +
     &             conjg( B3l ) * B3n * rsrsd) * w(k)

          O(l+NN,n+3*NN) = O(l+NN,n+3*NN) -
     &             (conjg( A3l ) * A2n + conjg( B3l ) *
     &             (B2n - rsrsd2 * D2n) * rsrsd ) * w(k)

          O(l+3*NN,n+0*NN) = O(l+3*NN,n+0*NN) -
     &             ( conjg( rsrsd1 * D2l ) * B3n *
     &             rsrsd ) * w(k)

          O(l+3*NN,n+2*NN) = O(l+3*NN,n+2*NN) +
     &             (conjg( rsrsd1 * D2l ) * F2n * rsrsd -
     &             conjg( B2l - rsrsd2 * D2l ) * rsrsd2 * C2n *
     &             rsrsd) * w(k)

          O(l+3*NN,n+1*NN) = O(l+3*NN,n+1*NN) -
     &             ( conjg( A2l ) * A3n + conjg( B2l -
     &             rsrsd2 * D2l ) * B3n * rsrsd ) * w(k)

          O(l+3*NN,n+3*NN) = O(l+3*NN,n+3*NN) +
     &             (conjg( A2l ) * A2n +
     &             (cdabs(ep12 - 1))**2 * (rsd / (rs*sints))**2*
     &             conjg( D2l ) * D2n * rsrsd +
     &             conjg( B2l - rsrsd2 * D2l ) * (B2n - rsrsd2 * D2n)*
     &             rsrsd) * w(k)


          O(l,n+4*NN) = O(l,n+4*NN)
     &           + (conjg(A3l) * A1n + conjg(B3l) * B1n * rsrsd) *w(k)

          O(l+2*NN,n+4*NN) = O(l+2*NN,n+4*NN)
     &           - (conjg(A2l) * A1n + conjg(F2l) * B1n * rsrsd) *w(k)

          O(l+3*NN,n+4*NN) = O(l+3*NN,n+4*NN)
     &           - (conjg(rsrsd1 * D2l) * B1n * rsrsd) * w(k)

      ENDDO
      ENDDO
      ENDDO
      end
