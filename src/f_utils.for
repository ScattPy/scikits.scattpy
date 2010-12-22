
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
      func_aa= (1d0-(e21-1d0)*r*(rd*cost-r*sint)/(r**2+rd**2)/sint)
     &        *func_alpha()
      RETURN
      END

      FUNCTION func_ab()
      implicit none
      COMPLEX*16 func_ab,  func_alpha
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_ab= -(e21-1d0)*(r**2*rd)/(r**2+rd**2)/sint
     &        *func_alpha()
      RETURN
      END

      FUNCTION func_ac()
      implicit none
      COMPLEX*16 func_ac,  func_alpha
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_ac= (e21-1d0)
     &        *(rd*sint+r*cost)*(rd*cost-r*sint)/(r*(r**2+rd**2))/sint
     &        *func_alpha()
      RETURN
      END

      FUNCTION func_ad()
      implicit none
      COMPLEX*16 func_ad,  func_alpha
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_ad= (1d0+(e21-1d0)*rd*(rd*sint+r*cost)/(r**2+rd**2)/sint)
     &           *func_alpha()
      RETURN
      END

      
      FUNCTION func_gamma()
      implicit none
      COMPLEX*16 func_gamma
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      func_gamma = ki*rd*radd*ang - sint*rad*angd
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
      func_d_te = func_beta() - (e21-1d0)*r*(
     &              rd*(rd*cost-r*sint)/(r**2+rd**2)/sint*func_gamma()
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
      func_e_te = -(e21-1d0)*r*(
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
      func_f_te = (e21-1d0)*r*(
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
      func_g_te = func_beta() + (e21-1d0)*r*(
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
      COMPLEX*16 FUNC,f_sint_w
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt

      ki=kki
      e12=e1e2
      e21=e2e1

      DO n=1,NN
      DO l=1,NN
         O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
            r   = Rs(k)
            rd  = Rds(k)
            rdd = Rdds(k)
            sint= Sints(k)
            cost= Costs(k)
            ctgt= Ctgts(k)
            DO l=1,NN
             rad = Rads(k,l)
             radd= Radds(k,l)
             ang = Angs(k,l)
             angd= Angds(k,l)
             f_sint_w = FUNC()*sint*w(k)
             DO n=1,NN
               O(l,n) = O(l,n) + f_sint_w*Angs(k,n)
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





      FUNCTION f1()
      implicit real*8 (r)
      real*8 f1
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      
      f1 = (r**2 - r*rdd + 2*rd**2) 
     &    *( r*(rd*cost-r*sint) + rd*(rd*sint+r*cost) ) 
     &    /(r**2+rd**2)**2
      RETURN
      end

      FUNCTION f2()
      implicit real*8 (r)
      real*8 f2
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      
      f2 = (r**4 - 2*r**3*rdd + 2*r**2*rd**2 - rd**4) 
     &    /(r**2+rd**2)**2
      RETURN
      end

      FUNCTION f3()
      implicit real*8 (r)
      real*8 f3
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      
      f3 = 2*(r**2 - r*rdd + 2*rd**2) 
     &    *(rd*cost-r*sint)*(rd*sint+r*cost) 
     &    /(r**2+rd**2)**2
      RETURN
      end

      FUNCTION f4()
      implicit real*8 (r)
      real*8 f4
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      
      f4 =((r**3*rdd + r**2*rd**2 - r*rd**2*rdd +3*rd**4)*sint
     &    +(rd/r)*cost*(r**4 - 2*r**3*rdd + 2*r**2*rd**2 - rd**4)) 
     &    /(r**2+rd**2)**2
      RETURN
      end





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc            E B C M
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MAT_XQ(NG,NN,FUNC1,FUNC2,FUNC3,FUNC4,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 k1,k2,e2e1,e1e2
      COMPLEX*16 Rads1(NG,NN),Radds1(NG,NN)
      COMPLEX*16 Rads2(NG,NN),Radds2(NG,NN)
      COMPLEX*16 Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      INTEGER   k,l,n
      EXTERNAL   FUNC1,FUNC2,FUNC3,FUNC4
      COMPLEX*16 FUNC1,FUNC2,FUNC3,FUNC4
      complex*16 fp1,fp3,c1
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt

      c1 = complex(0,1)
      e12=e1e2
      e21=e2e1

      DO n=1,NN
      DO l=1,NN
         O(l,n)=0d0
      ENDDO
      ENDDO

      DO k=1,NG
            r   = Rs(k)
            rd  = Rds(k)
            rdd = Rdds(k)
            sint= Sints(k)
            cost= Costs(k)
            ctgt= Ctgts(k)
            DO n=1,NN
             rad = Rads1(k,n)
             radd= Radds1(k,n)
             ang = Angs(k,n)
             angd= Angds(k,n)
             ki=k1
             fp1 = FUNC1()
             fp3 = FUNC3()
             DO l=1,NN
               rad = Rads2(k,l)
               radd= Radds2(k,l)
               ang = Angs(k,l)
               angd= Angds(k,l)
               ki=k2
               fp1=fp1*FUNC2()
               fp3=fp3*FUNC4()
               O(l,n) = O(l,n) + c1*(fp1-fp3)*sint*w(k)
             ENDDO
            ENDDO
      ENDDO
      END

      FUNCTION czero()
      complex*16 czero
      czero=0d0
      RETURN
      END

      SUBROUTINE MAT_Q_tm(NG,NN,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 k1,k2,e2e1,e1e2
      COMPLEX*16 Rads1(NG,NN),Radds1(NG,NN)
      COMPLEX*16 Rads2(NG,NN),Radds2(NG,NN)
      COMPLEX*16 Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_a,func_b,func_c_tm
      COMPLEX*16 func_a,func_b,func_c_tm
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_XQ(NG,NN,func_b,func_a,func_a,func_c_tm,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT_Q_te(NG,NN,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 k1,k2,e2e1,e1e2
      COMPLEX*16 Rads1(NG,NN),Radds1(NG,NN)
      COMPLEX*16 Rads2(NG,NN),Radds2(NG,NN)
      COMPLEX*16 Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_a,func_b,func_c_te
      COMPLEX*16 func_a,func_b,func_c_te
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_XQ(NG,NN,func_b,func_a,func_a,func_c_te,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT_Q11_te(NG,NN,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 k1,k2,e2e1,e1e2
      COMPLEX*16 Rads1(NG,NN),Radds1(NG,NN)
      COMPLEX*16 Rads2(NG,NN),Radds2(NG,NN)
      COMPLEX*16 Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_a,func_b,func_aa,func_d_te
      COMPLEX*16 func_a,func_b,func_aa,func_d_te
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_XQ(NG,NN,func_b,func_aa,func_a,func_d_te,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT_Q12_te(NG,NN,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 k1,k2,e2e1,e1e2
      COMPLEX*16 Rads1(NG,NN),Radds1(NG,NN)
      COMPLEX*16 Rads2(NG,NN),Radds2(NG,NN)
      COMPLEX*16 Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_a,func_b,func_ab,func_e_te
      COMPLEX*16 func_a,func_b,func_ab,func_e_te
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_XQ(NG,NN,func_b,func_ab,func_a,func_e_te,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT_Q21_te(NG,NN,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 k1,k2,e2e1,e1e2
      COMPLEX*16 Rads1(NG,NN),Radds1(NG,NN)
      COMPLEX*16 Rads2(NG,NN),Radds2(NG,NN)
      COMPLEX*16 Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_a,func_b,func_ac,func_f_te
      COMPLEX*16 func_a,func_b,func_ac,func_f_te
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_XQ(NG,NN,func_b,func_ac,func_a,func_f_te,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT_Q22_te(NG,NN,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 k1,k2,e2e1,e1e2
      COMPLEX*16 Rads1(NG,NN),Radds1(NG,NN)
      COMPLEX*16 Rads2(NG,NN),Radds2(NG,NN)
      COMPLEX*16 Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_a,func_b,func_ad,func_g_te
      COMPLEX*16 func_a,func_b,func_ad,func_g_te
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_XQ(NG,NN,func_b,func_ad,func_a,func_g_te,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END


      SUBROUTINE MAT_Q11_tm(NG,NN,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 k1,k2,e2e1,e1e2
      COMPLEX*16 Rads1(NG,NN),Radds1(NG,NN)
      COMPLEX*16 Rads2(NG,NN),Radds2(NG,NN)
      COMPLEX*16 Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_a,func_d_tm,czero
      COMPLEX*16 func_a,func_d_tm,czero
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_XQ(NG,NN,czero,czero,func_a,func_d_tm,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT_Q12_tm(NG,NN,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 k1,k2,e2e1,e1e2
      COMPLEX*16 Rads1(NG,NN),Radds1(NG,NN)
      COMPLEX*16 Rads2(NG,NN),Radds2(NG,NN)
      COMPLEX*16 Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_a,func_e_tm,czero
      COMPLEX*16 func_a,func_e_tm,czero
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_XQ(NG,NN,czero,czero,func_a,func_e_tm,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT_Q21_tm(NG,NN,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 k1,k2,e2e1,e1e2
      COMPLEX*16 Rads1(NG,NN),Radds1(NG,NN)
      COMPLEX*16 Rads2(NG,NN),Radds2(NG,NN)
      COMPLEX*16 Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_a,func_f_tm,czero
      COMPLEX*16 func_a,func_f_tm,czero
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_XQ(NG,NN,czero,czero,func_a,func_f_tm,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END

      SUBROUTINE MAT_Q22_tm(NG,NN,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      implicit none
      INTEGER NG,NN
      COMPLEX*16 O(NN,NN)
      COMPLEX*16 k1,k2,e2e1,e1e2
      COMPLEX*16 Rads1(NG,NN),Radds1(NG,NN)
      COMPLEX*16 Rads2(NG,NN),Radds2(NG,NN)
      COMPLEX*16 Angs(NG,NN),Angds(NG,NN)
      REAL*8     Rs(NG),Rds(NG),Rdds(NG)
      REAL*8     Sints(NG),Costs(NG),Ctgts(NG),w(NG)
cf2py intent(out) O
cf2py intent(hide) NG,NN
      EXTERNAL   func_a,func_g_tm,czero
      COMPLEX*16 func_a,func_g_tm,czero
      complex*16 ki,e12,e21,rad,radd,ang,angd
      real*8     r,rd,rdd,sint,cost,ctgt
      common /dk/ e12,e21,ki,rad,radd,ang,angd,r,rd,rdd,sint,cost,ctgt
      call MAT_XQ(NG,NN,czero,czero,func_a,func_g_tm,
     &                  k1,k2,e1e2,e2e1,
     &                  Rads1,Radds1,Rads2,Radds2,Angs,Angds,
     &                  Rs,Rds,Rdds,Sints,Costs,Ctgts,w,O)
      END
