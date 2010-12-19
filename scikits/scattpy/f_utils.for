
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

