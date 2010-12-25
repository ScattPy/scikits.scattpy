      SUBROUTINE get_I(nx,ny,nz,func_R,X,Y,Z,O)
      EXTERNAL func_R
      REAL*8   func_R
      INTEGER  nx,ny,nz, O(nx,ny,nz)
      REAL*8   X(nx),Y(ny),Z(nz)
Cf2py intent(hide) nx,ny,nz
Cf2py intent(out)  O
      INTEGER i,j,k
      REAL*8  r,t,rr
      DO i=1,nx
      DO j=1,ny
      DO k=1,nz
         r = DSQRT(X(i)**2+Y(j)**2+Z(k)**2)
         IF (r.eq.0d0) THEN
            t=0d0
         ELSE
            t=ACOS(Z(k)/r)
         ENDIF
         rr=func_R(t)
         IF(r<rr) THEN
            O(i,j,k)=0
         ELSE 
            O(i,j,k) = 1
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      END


      SUBROUTINE point_C2S ( CP,SP )
      IMPLICIT NONE
      REAL*8 CP(3),SP(3)
      REAL*8 x,y,z,r,t,p
cf2py intent(out) SP
      x = CP(1)
      y = CP(2)
      z = CP(3)
      r = sqrt(x*x + y*y + z*z)
      t = atan2(sqrt(x*x+y*y),z)
      p = atan2(y,x)
      SP(1) = r
      SP(2) = t
      SP(3) = p
      END

      SUBROUTINE point_S2C ( SP,CP )
      IMPLICIT NONE
      REAL*8 CP(3),SP(3)
      REAL*8 x,y,z,r,t,p
cf2py intent(out) CP
      r = SP(1)
      t = SP(2)
      p = SP(3)
      x = r*cos(p)*sin(t)
      y = r*sin(p)*sin(t)
      z = r*cos(t)
      CP(1) = x
      CP(2) = y
      CP(3) = z
      END

      SUBROUTINE vector_S2C ( SV,SP,CV )
      IMPLICIT NONE
      REAL*8 SP(3),TM(3,3)
      COMPLEX*16 SV(3),CV(3)
      REAL*8 r,t,p,sint,cost,sinp,cosp
      INTEGER i,j
cf2py intent(out) CV

      r = SP(1)
      t = SP(2)
      p = SP(3)
      sint = sin(t)
      cost = cos(t)
      sinp = sin(p)
      cosp = cos(p)

      TM(1,1) = sint*cosp
      TM(1,2) = cost*cosp
      TM(1,3) = -sinp
      TM(2,1) = sint*sinp
      TM(2,2) = cost*sinp
      TM(2,3) = cosp
      TM(3,1) = cost
      TM(3,2) = -sint
      TM(3,3) = 0d0

      do i=1,3
         CV(i)=0d0
         do j=1,3
            CV(i) = CV(i) + TM(i,j)*SV(j)
         enddo
      enddo
      END

      SUBROUTINE vector_C2S ( CV,SP,SV )
      IMPLICIT NONE
      REAL*8 SP(3),TM(3,3)
      COMPLEX*16 SV(3),CV(3)
      REAL*8 r,t,p,sint,cost,sinp,cosp
      INTEGER i,j  
cf2py intent(out) SV

      r = SP(1)
      t = SP(2)
      p = SP(3)
      sint = sin(t)
      cost = cos(t)
      sinp = sin(p)
      cosp = cos(p)

      TM(1,1) = sint*cosp
      TM(2,1) = cost*cosp
      TM(3,1) = -sinp
      TM(1,2) = sint*sinp
      TM(2,2) = cost*sinp
      TM(3,2) = cosp
      TM(1,3) = cost
      TM(2,3) = -sint
      TM(3,3) = 0d0

      do i=1,3
         SV(i)=0d0
         do j=1,3
            SV(i) = SV(i) + TM(i,j)*CV(j)
         enddo
      enddo
      END

