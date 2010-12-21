
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

      FUNCTION func_alpha()
      implicit none
      COMPLEX*16 func_alpha
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_alpha = rad*ang
      RETURN
      END

      FUNCTION func_beta()
      implicit none
      COMPLEX*16 func_beta
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_beta = ki*r**2*radd*ang + rd*sint*rad*angd
      RETURN
      END

      FUNCTION func_a()
      implicit none
      COMPLEX*16 func_a, func_alpha
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_a = func_alpha(rad,ang)
      RETURN
      END

      FUNCTION func_b()
      implicit none
      COMPLEX*16 func_b,  func_beta
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_b = func_beta(ki,rad,radd,ang,angd,r,rd,sint)
      RETURN
      END

      FUNCTION func_c_tm()
      implicit none
      COMPLEX*16 func_c_tm
      complex*16 func_beta,func_alpha
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_c_tm = e12*func_beta() + (e12-1d0)*(r-rd*ctgt)*func_alpha()
      RETURN
      END

      FUNCTION func_c_te()
      implicit none
      COMPLEX*16 func_c_te,  func_beta
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_c_te = func_beta()
      RETURN
      END


      FUNCTION func_aa()
      implicit none
      COMPLEX*16 func_aa,  func_alpha
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_aa= (1d0-(e21-1d0)*(r*(rd*cost-r*sint))/((r**2+rd**2)*sint))
     &        *func_alpha()
      RETURN
      END

      FUNCTION func_ab()
      implicit none
      COMPLEX*16 func_ab,  func_alpha
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_ab= -(e21-1d0)*(r**2*rd)/((r**2+rd**2)*sint)
     &        *func_alpha()
      RETURN
      END

      FUNCTION func_ac()
      implicit none
      COMPLEX*16 func_ac,  func_alpha
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_ac= (e21-1d0)*((rd*sint+r*cost)*(rd*cost-r*sint))
     &                 /(r*(r**2+rd**2)*sint)
     &        *func_alpha()
      RETURN
      END

      FUNCTION func_ad()
      implicit none
      COMPLEX*16 func_ad,  func_alpha
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_ad= (1d0+(e21-1d0)*(rd*(rd*sint+r*cost))/((r**2+rd**2)*sint))
     &           *func_alpha()
      RETURN
      END

      
      FUNCTION func_gamma()
      implicit none
      COMPLEX*16 func_gamma
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_gamma = ki*radd*ang - sint*rad*angd
      RETURN
      END


      FUNCTION func_d_te()
      implicit none
      COMPLEX*16 func_d_te
      complex*16 func_alpha,func_beta,func_gamma
      real*8     f1
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_d_te = func_beta() - (e21-1d0)*(
     &              rd*(rd*cost-r*sint)/((r**2+rd**2)*sint)*func_gamma()
     &             -f1()/sint*func_alpha())
      RETURN
      END

      FUNCTION func_e_te()
      implicit none
      COMPLEX*16 func_e_te
      complex*16 func_alpha,func_gamma
      real*8     f2
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_e_te = -(e21-1d0)*(
     &                   r*rd**2/(r**2+rd**2)/sint*func_gamma()
     &                  -rd*f2()/sint*func_alpha())
      RETURN
      END

      FUNCTION func_f_te()
      implicit none
      COMPLEX*16 func_f_te
      complex*16 func_alpha,func_gamma
      real*8     f3
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_f_te = (e21-1d0)*(
     &              (rd*cost-r*sint)**2/r/(r**2+rd**2)/sint*func_gamma()
     &             -f3()/sint/r*func_alpha())
      RETURN
      END

      FUNCTION func_g_te()
      implicit none
      COMPLEX*16 func_g_te
      complex*16 func_alpha,func_beta,func_gamma
      real*8     f4
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_g_te = func_beta() + (e21-1d0)*(
     &                rd*(rd*cost-r*sint)/(r**2+rd**2)/sint*func_gamma()
     &               -f4()/sint*func_alpha())
      RETURN
      END



      FUNCTION func_delta()
      implicit none
      COMPLEX*16 func_delta
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_delta = ki*r*cost*radd*ang + sint**2*rad*angd
      RETURN
      END

      FUNCTION func_eps()
      implicit none
      COMPLEX*16 func_eps
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_eps = (ki*r**2*radd+r*rad)*ang
      RETURN
      END


      FUNCTION func_d_tm()
      implicit none
      COMPLEX*16 func_d_tm
      complex*16 func_beta,func_delta
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_d_tm = func_beta()+(e12-1d0)*rd/sint*func_delta()
      RETURN
      END

      FUNCTION func_e_tm()
      implicit none
      COMPLEX*16 func_e_tm
      complex*16 func_eps
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_e_tm = (e12-1d0)*rd/sint*func_eps()
      RETURN
      END

      FUNCTION func_f_tm()
      implicit none
      COMPLEX*16 func_f_tm
      complex*16 func_delta
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_f_tm = -(e12-1d0)*(rd*cost-r*sint)/r/sint*func_delta()
      RETURN
      END

      FUNCTION func_g_tm()
      implicit none
      COMPLEX*16 func_g_tm
      complex*16 func_beta,func_eps
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_g_tm = func_beta()
     &           -(e12-1d0)*(rd*cost-r*sint)/r/sint*func_eps()
      RETURN
      END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc            S V M
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MAT_X(NG,NN,FUNC,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      INTEGER   k,l,n
      EXTERNAL   FUNC
      COMPLEX*16 FUNC
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt

      ki=kki
      e12=e1e2
      e21=e2e1

      DO l=1,NN
      DO n=1,NN
         O(l,n)=0d0
         DO k=1,NG
            rad = Rads(k,l)
            radd= Radds(k,l)
            ang = Angs(k,l)
            angd= Angds(k,l)
            r   = Rs(k)
            rd  = Rds(k)
            rdd = Rdds(k)
            sint= Sints(k)
            cost= Costs(k)
            ctgt= Ctgts(k)
            O(l,n) = O(l,n) + FUNC()*Angs(k,n)*Sints(k)*w(k)
         ENDDO
      ENDDO
      ENDDO
      END


      SUBROUTINE MAT__A(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_a
      COMPLEX*16 func_a
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_a,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__B(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_b
      COMPLEX*16 func_b
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_b,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__C_tm(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_c_tm
      COMPLEX*16 func_c_tm
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_c_tm,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__D_tm(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_d_tm
      COMPLEX*16 func_d_tm
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_d_tm,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__E_tm(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_e_tm
      COMPLEX*16 func_e_tm
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_e_tm,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__F_tm(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_f_tm
      COMPLEX*16 func_f_tm
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_f_tm,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__G_tm(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_g_tm
      COMPLEX*16 func_g_tm
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_g_tm,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__D_te(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_d_te
      COMPLEX*16 func_d_te
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_d_te,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__E_te(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_e_te
      COMPLEX*16 func_e_te
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_e_te,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__F_te(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_f_te
      COMPLEX*16 func_f_te
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_f_te,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__G_te(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_g_te
      COMPLEX*16 func_g_te
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_g_te,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__Aa(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_aa
      COMPLEX*16 func_aa
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_aa,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__Ab(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_ab
      COMPLEX*16 func_ab
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_ab,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__Ac(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_ac
      COMPLEX*16 func_ac
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_ac,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT__Ad(NG,NN,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 kki,e2e1,e1e2
      COMPLEX*16 Rads(NG,NN),Radds(NG,NN),Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_ad
      COMPLEX*16 func_ad
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_X(NG,NN,func_ad,kki,e1e2,e2e1,Rads,Radds,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
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




      FUNCTION f1()
      implicit real*8 (r)
      real*8 f1
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      
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

      FUNCTION f2()
      implicit real*8 (r)
      real*8 f2
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      
      rdr = rd/r
      r2r2 = r**2+rd**2
      r2r22 = r2r2**2

      f2 = (r**4 - 2*(r**3)*rdd 
     &          + 2*(r**2)*(rd**2) - rd**4)
     &         /r2r22
      RETURN
      end

      FUNCTION f3()
      implicit real*8 (r)
      real*8 f3
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      
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

      FUNCTION f4()
      implicit real*8 (r)
      real*8 f4
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      
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
