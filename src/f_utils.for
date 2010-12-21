
      FUNCTION aa(x)
      implicit none
      REAL*8 aa,x
      aa = x**2
      RETURN
      END

      FUNCTION bb(f,x)
      implicit none
      REAL*8 bb,x
      EXTERNAL f
      REAL*8 f
      bb = 1+ f(x)
      RETURN
      END

      FUNCTION func_alpha(rad,ang)
      implicit none
      COMPLEX*16 func_alpha
      COMPLEX*16 rad,ang
      func_alpha = rad*ang
      RETURN
      END

      FUNCTION func_beta(ki,rad,radd,ang,angd,r,rd,sint)
      implicit none
      COMPLEX*16 func_beta
      complex*16 rad,radd,ang,angd
      real*8     ki,r,rd,sint
      func_beta = ki*r**2*radd*ang+rd*sint*rad*angd
      RETURN
      END

      FUNCTION func_a(rad,ang)
      implicit none
      COMPLEX*16 func_a, func_alpha
      COMPLEX*16 rad,ang
      func_a = func_alpha(rad,ang)
      RETURN
      END

      FUNCTION func_b(ki,rad,radd,ang,angd,r,rd,sint)
      implicit none
      COMPLEX*16 func_b,  func_beta
      complex*16 rad,radd,ang,angd
      real*8     ki,r,rd,sint
      func_b = func_beta(ki,rad,radd,ang,angd,r,rd,sint)
      RETURN
      END

      FUNCTION func_c_tm(ki,rad,radd,ang,angd,r,rd,sint,ctgt,e12)
      implicit none
      COMPLEX*16 func_c_tm
      complex*16 rad,radd,ang,angd,e12
      real*8     ki,r,rd,sint,ctgt
      complex*16 f_beta,f_alpha,  func_beta,func_alpha
      f_beta = func_beta(ki,rad,radd,ang,angd,sint)
      f_alpha= func_alpha(rad,ang)
      func_c_tm = e12*f_beta + (e12-1d0)*(r-rd*ctgt)*f_alpha
      RETURN
      END

      FUNCTION func_c_te(ki,rad,radd,ang,angd,r,rd,sint)
      implicit none
      COMPLEX*16 func_c_te,  func_beta
      complex*16 rad,radd,ang,angd
      real*8     ki,r,rd,sint
      func_c_te = func_beta(ki,rad,radd,ang,angd,r,rd,sint)
      RETURN
      END


      FUNCTION func_aa(rad,ang,r,rd,sint,cost,e21)
      implicit none
      COMPLEX*16 func_aa,  func_alpha
      complex*16 rad,ang,e21
      real*8     r,rd,sint,cost
      func_aa= (1d0-(e21-1d0)*(r*(rd*cost-r*sint))/((r**2+rd**2)*sint))
     &        *func_alpha(rad,ang)
      RETURN
      END

      FUNCTION func_ab(rad,ang,r,rd,sint,cost,e21)
      implicit none
      COMPLEX*16 func_ab,  func_alpha
      complex*16 rad,ang,e21
      real*8     r,rd,sint,cost
      func_ab= -(e21-1d0)*(r**2*rd)/((r**2+rd**2)*sint)
     &        *func_alpha(rad,ang)
      RETURN
      END

      FUNCTION func_ac(rad,ang,r,rd,sint,cost,e21)
      implicit none
      COMPLEX*16 func_ac,  func_alpha
      complex*16 rad,ang,e21
      real*8     r,rd,sint,cost
      func_ac= (e21-1d0)*((rd*sint+r*cost)*(rd*cost-r*sint))
     &                 /(r*(r**2+rd**2)*sint)
     &        *func_alpha(rad,ang)
      RETURN
      END

      FUNCTION func_ad(rad,ang,r,rd,sint,cost,e21)
      implicit none
      COMPLEX*16 func_ad,  func_alpha
      complex*16 rad,ang,e21
      real*8     r,rd,sint,cost
      func_ad= (1d0+(e21-1d0)*(rd*(rd*sint+r*cost))/((r**2+rd**2)*sint))
     &           *func_alpha(rad,ang)
      RETURN
      END

      
      FUNCTION func_gamma(ki,rad,radd,ang,angd,r,rd,sint)
      implicit none
      COMPLEX*16 func_gamma
      complex*16 rad,radd,ang,angd,ki
      real*8     r,rd,sint
      func_gamma = ki*radd*ang - sint*rad*angd
      RETURN
      END


      FUNCTION func_d_te(ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,e21)
      implicit none
      COMPLEX*16 func_d_te
      complex*16 rad,radd,ang,angd,ki,e21
      real*8     r,rd,rdd,sint,cost
      complex*16 f_alpha,f_beta,f_gamma, func_alpha,func_beta,func_gamma
      real*8     f1
      f_alpha = func_alpha(rad,ang)
      f_beta  = func_beta(ki,rad,radd,ang,angd,r,rd,sint)
      f_gamma = func_beta(ki,rad,radd,ang,angd,r,rd,sint)
      func_d_te = f_beta - (e21-1d0)*(
     &                   rd*(rd*cost-r*sint)/((r**2+rd**2)*sint)*f_gamma
     &                  -f1(r,rd,rdd,sint,cost)/sint*f_alpha)
      RETURN
      END

      FUNCTION func_e_te(ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,e21)
      implicit none
      COMPLEX*16 func_e_te
      complex*16 rad,radd,ang,angd,ki,e21
      real*8     r,rd,rdd,sint,cost
      complex*16 f_alpha,f_gamma, func_alpha,func_gamma
      real*8     f2
      f_alpha = func_alpha(rad,ang)
      f_gamma = func_gamma(ki,rad,radd,ang,angd,r,rd,sint)
      func_e_te = -(e21-1d0)*(
     &                   r*rd**2/(r**2+rd**2)/sint*f_gamma
     &                  -rd*f2(r,rd,rdd)/sint*f_alpha)
      RETURN
      END

      FUNCTION func_f_te(ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,e21)
      implicit none
      COMPLEX*16 func_f_te
      complex*16 rad,radd,ang,angd,ki,e21
      real*8     r,rd,rdd,sint,cost
      complex*16 f_alpha,f_gamma, func_alpha,func_gamma
      real*8     f3
      f_alpha = func_alpha(rad,ang)
      f_gamma = func_gamma(ki,rad,radd,ang,angd,r,rd,sint)
      func_f_te = (e21-1d0)*(
     &                   (rd*cost-r*sint)**2/r/(r**2+rd**2)/sint*f_gamma
     &                  -f3(r,rd,rdd,sint,cost)/sint/r*f_alpha)
      RETURN
      END

      FUNCTION func_g_te(ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,e21)
      implicit none
      COMPLEX*16 func_g_te
      complex*16 rad,radd,ang,angd,ki,e21
      real*8     r,rd,rdd,sint,cost
      complex*16 f_alpha,f_beta,f_gamma, func_alpha,func_beta,func_gamma
      real*8     f4
      f_alpha = func_alpha(rad,ang)
      f_beta  = func_beta(ki,rad,radd,ang,angd,r,rd,sint)
      f_gamma = func_gamma(ki,rad,radd,ang,angd,r,rd,sint)
      func_g_te = f_beta + (e21-1d0)*(
     &                   rd*(rd*cost-r*sint)/(r**2+rd**2)/sint*f_gamma
     &                  -f4(r,rd,rdd,sint,cost)/sint*f_alpha)
      RETURN
      END



      FUNCTION func_delta(ki,rad,radd,ang,angd,r,sint,cost)
      implicit none
      COMPLEX*16 func_delta
      complex*16 rad,radd,ang,angd,ki
      real*8     r,sint,cost
      func_delta = ki*r*cost*radd*ang + sint**2*rad*angd
      RETURN
      END

      FUNCTION func_eps(ki,rad,radd,ang,r)
      implicit none
      COMPLEX*16 func_eps
      COMPLEX*16 rad,radd,ang,ki
      REAL*8     r
      func_eps = (ki*r**2*radd+r*rad)*ang
      RETURN
      END


      FUNCTION func_d_tm(ki,rad,radd,ang,angd,r,rd,sint,cost,e12)
      implicit none
      COMPLEX*16 func_d_tm
      complex*16 rad,radd,ang,angd,ki,e12
      real*8     r,rd,sint,cost
      complex*16 f_beta,f_delta, func_beta,func_delta
      f_beta = func_beta(ki,rad,radd,ang,angd,r,rd,sint)
      f_delta = func_delta(ki,rad,radd,ang,angd,r,sint,cost)
      func_d_tm = f_beta+(e12-1d0)*rd/sint*f_delta
      RETURN
      END

      FUNCTION func_e_tm(ki,rad,radd,ang,angd,r,rd,sint,e12)
      implicit none
      COMPLEX*16 func_e_tm
      complex*16 ki,rad,radd,ang,angd,e12
      real*8     r,rd,sint
      complex*16 f_eps,  func_eps
      f_eps = func_eps(ki,rad,radd,ang,r)
      func_e_tm = (e12-1d0)*rd/sint*f_eps
      RETURN
      END

      FUNCTION func_f_tm(ki,rad,radd,ang,angd,r,rd,sint,cost,e12)
      implicit none
      COMPLEX*16 func_f_tm
      complex*16 ki,rad,radd,ang,angd,e12
      real*8     r,rd,sint,cost
      complex*16 f_delta,  func_delta
      f_delta = func_delta(ki,rad,radd,ang,angd,r,sint,cost)
      func_f_tm = -(e12-1d0)*(rd*cost-r*sint)/r/sint*f_delta
      RETURN
      END

      FUNCTION func_g_tm(ki,rad,radd,ang,angd,r,rd,sint,cost,e12)
      implicit none
      COMPLEX*16 func_g_tm
      complex*16 rad,radd,ang,angd,ki,e12
      real*8     r,rd,sint,cost
      complex*16 f_beta,f_eps,  func_beta,func_eps
      f_beta = func_beta(ki,rad,radd,ang,angd,r,rd,sint)
      f_eps = func_eps(ki,rad,radd,ang,angd,r)
      func_g_tm = f_beta-(e12-1d0)*(rd*cost-r*sint)/r/sint*f_eps
      RETURN
      END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc            S V M
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MAT__A(NG,NN,Rad,Ang,Sint,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 Rad(NG,NN),Ang(NG,NN),O(NN,NN)
      REAL*8     Sint(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      INTEGER   k,l,n
      COMPLEX*16 func_a

      DO l=1,NN
      DO n=1,NN
         O(l,n)=0d0
         DO k=1,NG
            O(l,n) = O(l,n) + func_a(Rad(k,l),Ang(k,l))
     &                         *Ang(k,n)*Sint(k)*w(k)
         ENDDO
      ENDDO
      ENDDO
      END

      SUBROUTINE MAT__B(NG,NN,ki,Rad,Radd,Ang,Angd,R,Rd,Sint,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 Rad(NG,NN),Radd(NG,NN),Ang(NG,NN),Angd(NG,NN),O(NN,NN)
      COMPLEX*16 ki
      REAL*8     R(NG),Rd(NG),Sint(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      INTEGER   k,l,n
      COMPLEX*16 func_b

      DO l=1,NN
      DO n=1,NN
         O(l,n)=0d0
         DO k=1,NG
            O(l,n) = O(l,n) 
     &               +func_b(ki,Rad(k,l),Radd(k,l),Ang(k,l),Angd(k,l),
     &                          R(k),Rd(k),Sint(k))
     &                *Ang(k,n)*Sint(k)*w(k)
         ENDDO
      ENDDO
      ENDDO
      END

      SUBROUTINE MAT__C(NG,NN,ki,Rad,Radd,Ang,Angd,R,Rd,
     &                                      Sint,Ctgt,e12,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 Rad(NG,NN),Radd(NG,NN),Ang(NG,NN),Angd(NG,NN)
      COMPLEX*16 ki,e12
      REAL*8     R(NG),Rd(NG),Sint(NG),Ctgt(NG),w(NG)
cf2py intent(hide) NG,NN
cf2py intent(out)  O
      INTEGER   k,l,n
      COMPLEX*16 func_c_tm

      DO l=1,NN
      DO n=1,NN
         O(l,n)=0d0
         DO k=1,NG
            O(l,n) = O(l,n) 
     &              +func_c_tm(ki,Rad(k,l),Radd(k,l),Ang(k,l),Angd(k,l),
     &                          R(k),Rd(k),Sint(k),Ctgt(k),e12)
     &               *Ang(k,n)*Sint(k)*w(k)
         ENDDO
      ENDDO
      ENDDO
      END

      



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


      SUBROUTINE f_f1(NG,R,Rd,Rdd,Sint,Cost,f)
      implicit real*8 (r)
      INTEGER NG
      REAL*8 R(NG),Rd(NG),Sint(NG),Cost(NG),f(NG)
cf2py intent(in) R,Rd,Sint,Cost
cf2py intent(out) f
cf2py intent(hide) NG
      INTEGER K
      real*8  f1
      
      DO k=1,NG
        f(k) = f1(R(k),Rd(k),Rdd(k),Sint(k),Cost(k))
      ENDDO
      end

      SUBROUTINE f_f2(NG,R,Rd,Rdd,f)
      implicit real*8 (r)
      INTEGER NG
      REAL*8 R(NG),Rd(NG),Sint(NG),Cost(NG),f(NG)
cf2py intent(in) R,Rd,Sint,Cost
cf2py intent(out) f
cf2py intent(hide) NG
      INTEGER K
      real*8  f2
      
      DO k=1,NG
        f(k) = f2(R(k),Rd(k),Rdd(k))
      ENDDO
      end

      SUBROUTINE f_f3(NG,R,Rd,Rdd,Sint,Cost,f)
      implicit real*8 (r)
      INTEGER NG
      REAL*8 R(NG),Rd(NG),Sint(NG),Cost(NG),f(NG)
cf2py intent(in) R,Rd,Sint,Cost
cf2py intent(out) f
cf2py intent(hide) NG
      INTEGER K
      real*8  f3
      
      DO k=1,NG
        f(k) = f3(R(k),Rd(k),Rdd(k),Sint(k),Cost(k))
      ENDDO
      end

      SUBROUTINE f_f4(NG,R,Rd,Rdd,Sint,Cost,f)
      implicit real*8 (r)
      INTEGER NG
      REAL*8 R(NG),Rd(NG),Sint(NG),Cost(NG),f(NG)
cf2py intent(in) R,Rd,Sint,Cost
cf2py intent(out) f
cf2py intent(hide) NG
      INTEGER K
      real*8  f4
      
      DO k=1,NG
        f(k) = f4(R(k),Rd(k),Rdd(k),Sint(k),Cost(k))
      ENDDO
      end




      FUNCTION f1(r,rd,rdd,sint,cost)
      implicit real*8 (r)
      real*8 f1,sint,cost
      
      rdr = rd/r
      r2r2 = r**2+rd**2
      r2r22 = r2r2**2
      rsin = r *sint
      rcos = r *cost
      rdsin= rd*sint
      rdcos= rd*cost
      rdcosrsin = rdcos-rsin
      rdsinrcos = rdsin+rcos

      f1 = (r**2 - r*rdd + 2*rd**2) 
     &        *( r *rdcosrsin + rd*rdsinrcos) 
     &        /r2r22
      RETURN
      end

      FUNCTION f2(r,rd,rdd)
      implicit real*8 (r)
      real*8 f2
      
      rdr = rd/r
      r2r2 = r**2+rd**2
      r2r22 = r2r2**2

      f2 = (r**4 - 2*(r**3)*rdd 
     &          + 2*(r**2)*(rd**2) - rd**4)
     &         /r2r22
      RETURN
      end

      FUNCTION f3(r,rd,rdd,sint,cost)
      implicit real*8 (r)
      real*8 f3,sint,cost
      
      rdr = rd/r
      r2r2 = r**2+rd**2
      r2r22 = r2r2**2
      rsin = r *sint
      rcos = r *cost
      rdsin= rd*sint
      rdcos= rd*cost
      rdcosrsin = rdcos-rsin
      rdsinrcos = rdsin+rcos

      f3 = 2*(r**2 - r*rdd + 2*(rd**2)) 
     &          *rdcosrsin*rdsinrcos 
     &          /r2r22
      RETURN
      end

      FUNCTION f4(r,rd,rdd,sint,cost)
      implicit real*8 (r)
      real*8 f4,sint,cost
      
      rdr = rd/r
      r2r2 = r**2+rd**2
      r2r22 = r2r2**2
      rsin = r *sint
      rcos = r *cost
      rdsin= rd*sint
      rdcos= rd*cost
      rdcosrsin = rdcos-rsin
      rdsinrcos = rdsin+rcos

      f4 = ( ((r**3)*rdd + (r**2)*(rd**2) 
     &            -r*(rd**2)*rdd + 3*(rd**4)) *sint 
     &          +rdr*cost
     &           *(r**4 - 2*(r**3)*rdd 
     &             +2*(r**2)*(rd**2) - rd**4)) 
     &         /r2r22
      RETURN
      end
