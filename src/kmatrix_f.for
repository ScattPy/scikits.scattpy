
      SUBROUTINE vec2kmatrix(NG,N,VEC,O)
      IMPLICIT NONE
      INTEGER NG,N
      COMPLEX*16 VEC(NG,N),O(NG,1,N)
cf2py intent(in) NG,N
cf2py intent(out) O
cf2py intent(hide) NG,N
      INTEGER i,j
      do i=1,NG
      do j=1,N
        O(i,1,j)=VEC(i,j)
      enddo
      enddo
      END

      SUBROUTINE coef_vec(NG,N,COEF,VEC,O)
      IMPLICIT NONE
      INTEGER NG,N
      COMPLEX*16 COEF(NG),VEC(NG,N),O(NG,N)
cf2py intent(in) NG,N
cf2py intent(out) O
cf2py intent(hide) NG,N
      INTEGER i,j
      do i=1,NG
      do j=1,N
        O(i,j)=COEF(i)*VEC(i,j)
      enddo
      enddo
      END

      SUBROUTINE mat_mat_dot(NG,M,N,MAT1,MAT2,O)
      IMPLICIT NONE
      INTEGER NG,M,N
      COMPLEX*16 MAT1(NG,M,N),MAT2(NG,M,N),O(NG,M,N)
cf2py intent(in) NG,N,M
cf2py intent(out) O
cf2py intent(hide) NG,N,M
      INTEGER i,j,k
      do i=1,NG
      do j=1,M
      do k=1,N
        O(i,j,k)=MAT1(i,j,k)*MAT2(i,j,k)
      enddo
      enddo
      enddo
      END

      SUBROUTINE coef_mat(NG,M,N,COEF,MAT,O)
      IMPLICIT NONE
      INTEGER NG,M,N
      COMPLEX*16 COEF(NG),MAT(NG,M,N),O(NG,M,N)
cf2py intent(in) NG,N,M
cf2py intent(out) O
cf2py intent(hide) NG,M,N
      INTEGER i,j,k
      do k=1,NG
      do i=1,M
      do j=1,N
        O(k,i,j)=COEF(k)*MAT(k,i,j)
      enddo
      enddo
      enddo
      END

      SUBROUTINE mat_mat(NG,L,M,N,MAT1,MAT2,O)
      IMPLICIT NONE
      INTEGER NG,L,M,N
      COMPLEX*16 MAT1(NG,L,M),MAT2(NG,M,N),O(NG,L,N)
cf2py intent(in) NG,L,M,N
cf2py intent(out) O
cf2py intent(hide) NG,L,M,N
      INTEGER i,j,k,t
      do i=1,NG
      do j=1,L
      do k=1,N
        O(i,j,k)=0d0
        do t=1,M
          O(i,j,k)=O(i,j,k) + MAT1(i,j,t)*MAT2(i,t,k)
        enddo
      enddo
      enddo
      enddo
      END

      SUBROUTINE fixmat_mat(NG,L,M,N,FIXMAT,MAT2,O)
      IMPLICIT NONE
      INTEGER NG,L,M,N
      COMPLEX*16 FIXMAT(L,M),MAT2(NG,M,N),O(NG,L,N)
cf2py intent(in) NG,L,M,N
cf2py intent(out) O
cf2py intent(hide) NG,L,M,N
      INTEGER i,j,k,t
      do i=1,NG
      do j=1,L
      do k=1,N
        O(i,j,k)=0d0
        do t=1,M
          O(i,j,k)=O(i,j,k) + FIXMAT(j,t)*MAT2(i,t,k)
        enddo
      enddo
      enddo
      enddo
      END

      SUBROUTINE mat_mfixat(NG,L,M,N,MAT1,FIXMAT,O)
      IMPLICIT NONE
      INTEGER NG,L,M,N
      COMPLEX*16 MAT1(NG,L,M),FIXMAT(M,N),O(NG,L,N)
cf2py intent(in) NG,L,M,N
cf2py intent(out) O
cf2py intent(hide) NG,L,M,N
      INTEGER i,j,k,t
      do i=1,NG
      do j=1,L
      do k=1,N
        O(i,j,k)=0d0
        do t=1,M
          O(i,j,k)=O(i,j,k) + MAT1(i,j,t)*FIXMAT(t,k)
        enddo
      enddo
      enddo
      enddo
      END

      SUBROUTINE sum_fixmat_mat(NG,L,M,N,FIXMAT,MAT2,O)
      IMPLICIT NONE
      INTEGER NG,L,M,N
      COMPLEX*16 FIXMAT(L,M),MAT2(NG,M,N),O(NG,L,N)
cf2py intent(in) NG,L,M,N
cf2py intent(out) O
cf2py intent(hide) NG,L,M,N
      INTEGER i,j,k
      do i=1,NG
      do j=1,L
      do k=1,N
          O(i,j,k)=FIXMAT(j,k)+MAT2(i,j,k)
      enddo
      enddo
      enddo
      END

      SUBROUTINE mat_transpose(NG,M,N,MAT,O)
      IMPLICIT NONE
      INTEGER NG,M,N
      COMPLEX*16 MAT(NG,M,N),O(NG,N,M)
cf2py intent(in) NG,M,N
cf2py intent(out) O
cf2py intent(hide) NG,M,N
      INTEGER i,j,k
      do i=1,NG
      do j=1,N
      do k=1,M
        O(i,j,k)=MAT(i,k,j)
      enddo
      enddo
      enddo
      END

      SUBROUTINE mat_integrate(NG,M,N,MAT,W,O)
      IMPLICIT NONE
      INTEGER M,N,NG
      COMPLEX*16 MAT(NG,M,N), O(M,N)
      REAL*8 W(NG)
cf2py intent(in) NG,M,N
cf2py intent(out) O
cf2py intent(hide) NG,M,N
      INTEGER i,j,k
      do i=1,M
      do j=1,N
        O(i,j)=0d0
      enddo
      enddo

      do k=1,NG
      do i=1,M
      do j=1,N
        O(i,j)=O(i,j)+W(k)*MAT(k,i,j)
      enddo
      enddo
      enddo
      END

