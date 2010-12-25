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
