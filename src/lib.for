c************************************************
C Lambda
C
      SUBROUTINE Lambda(K,M,ne,C2,EPS,rlc2,ie)
      parameter (nterms=330, nmmax=40)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      REAL*8 stepr, stepi, imc2, stepi0
      COMPLEX*16 C2, cc2a, cc20
      DIMENSION RLC2(nterms), RLcC2(nterms), rlc20(nmmax,nterms)
      COMMON /K1/ S, S1, AKSI, AK, Kkk, NK, nal
c      COMMON /icc/ ic0, m_max

c      icc = ic0
c-----------------
      icc = 1
      m_max = 1
c-----------------
      cC2 = dreal(c2)
      imC2 = dimag(c2)

      mm = m
      if(m.gt.nmmax.or.m.gt.m_max) icc = 1

c        print *, m, icc
c        print *, c2, cc20
c++++
      if(k.eq.0)   stepi0 = 0.032d0
      if(k.eq.1)   stepi0 = 0.091d0
      if(k.eq.1.and.cc2.gt.22d0)   stepi0 = 0.018d0
      if(k.eq.1.and.cc2.gt.32d0)   stepi0 = 0.0082d0

               stepi0 = 0.0021d0

         if(stepi0.gt.imc2) then
         stepi0 = imc2 / 3d0
         icc = 1
         go to 484
         end if
c++++

c*******************
             if(icc.ne.1.and.imc2.gt.1d-2) then
c             if(icc.ne.1) then

               if(m.le.m_max) rc20 = cc20
               if(m.le.m_max) ne0  = ne00
                 do ij = 1, ne0
                    rlcc2(ij) = rlc20(mm,ij)
                 end do

             if(dabs(cc2-dreal(rc20)).lt.1d-3) go to 786

                if(ne.ne.ne0) then
                do ij = ne0+1, ne
                rlcc2(ij) = s
                end do
                end if

         cc2a = rc20
         stepi = stepi0
         stepr = stepi/(imc2-dimag(cc2a))*(c2-dreal(cc2a))
         go to 787

  786        continue
             icc = 1
  787        continue
             end if
c*******************

c++++
          if(k.eq.0.and.imc2.lt.0.5d0) then
          CALL CDLAMp(RLC2,M,NE,C2,EPS,IE)
          if(ie.ne.0) print *,c2, ie
          if(ie.ne.0) go to 484
          go to 443
          end if
c&-------
                if(k.eq.1) then
                ie = 10
      if(imc2.lt.0.11d0.and.cc2.le.22d0)
     *                  CALL CDLAMo_s(rlc2,M,NE,C2,EPS,IE)
      if(imc2.lt.0.3d0.and.cc2.gt.22d0)
     *                  CALL CDLAMo_l(rlc2,M,NE,C2,EPS,IE)
         if(ie.ne.0) then
         if(ie.ne.10) print *,c2, ie
c                      print *,c2, ie
c         icc = 1
         if(ie.ne.10) icc = 1
         go to 484
         end if
         go to 443
                end if
c&-------

 484        continue

c-------------
       if(icc.eq.1) then
             stepr = 0d0
             stepi = stepi0
                  ac2 = 0.091d0
                  if(imc2.lt.0.1d0) ac2 = 0.0d0
       if(k.eq.1) ac2 = 0d0
       cc2a = DCMPLX(cc2,ac2)

c###
             if(k.eq.0) then
1443    continue
             CALL CDLAMp(RLcC2,M,NE,cc2a,EPS,IE)

         if(ie.ne.0) then
         cc2a = DCMPLX(cc2,0d0)
         go to 1443
         end if
              end if

c###
             if(k.eq.1) then
         if(dreal(cc2a).le.22d0) CALL CDLAMo_s(rlcc2,M,NE,cc2a,EPS,IE)
         if(dreal(cc2a).gt.22d0) CALL CDLAMo_l(rlcc2,M,NE,cc2a,EPS,IE)
             if(cdabs(c2)-cdabs(cc2a).lt.1d-5) then
             do jj = 1, ne
             rlc2(jj) = rlcc2(jj)
             end do
             go to 1445
             end if
             end if
        end if
c-------------

1444    continue
              cc2a = cc2a + dcmplx(stepr,stepi)
              aimc2 = imc2 - dimag(cc2a)
              if(dabs(dreal(cc2a)-cc2).lt.stepr.or.aimc2.lt.0d0)
     *        cc2a = c2

              epss = eps
              if(k.eq.0) epss = 1d-10
              CALL CDLAMn(rlcc2,rlc2,M,NE,cc2a,EPSs,IE)
              if(ie.ne.0) print *,cc2a, ie
c                          print *,cc2a, ie

c            write(7,*)
c            write(7,*) cc2a
c            do jj = 1,5
c            write(7,*) jj, rlc2(jj)
c            end do
*))))))))))))))
         if(ie.ne.0) then
         cc2a = cc2a - 2d0 * dcmplx(stepr,stepi)

         if(icc.eq.1) then
         if(dimag(cc2a).lt.0d0) cc2a = dcmplx(dreal(cc2a),0d0)
         stepi = stepi/2.1d0
ccc         if(stepi.lt.1d-6) stop 000
         if(stepi.lt.1d-6) return
         end if

c..............................
         if(icc.ne.1) then

             if(dreal(cc2a).le.dreal(rc20)) then
             icc = 1
             go to 484
             end if

         stepi = stepi / 2.1d0
         if(stepi.lt.1d-4) then
         icc = 1
         go to 484
         end if
            stepr = stepi/(imc2-dimag(cc2a))*(c2-dreal(cc2a))
         end if
c..............................

         go to 1444
         end if
*))))))))))))))

           do jj=1,ne
             rlcc2(jj)=rlc2(jj)
           end do
              if(cdabs(c2)-cdabs(cc2a).lt.1d-5) go to 1445
              go to 1444
 1445    continue
c++++++++
c++++++++
  443 CONTINUE

               if(m.gt.m_max) ne00 = ne
               if(m.gt.m_max) cc20 = c2
                 do ij = 1, ne
                    rlc20(mm,ij) = rlc2(ij)
                 end do

c            do jj = 1,ne
c            write(*,*) jj, rlc2(jj)
c            write(7,*) jj, rlc2(jj)
c            end do

       return
       end

c************************************************
C BESSJJ0
c
      SUBROUTINE BESSJJ0(A,NUM,BESJ)
      parameter (nterm=330)
      IMPLICIT COMPLEX*16 (A-H,O-Q,T-Z)
      DIMENSION BESJ(NUM+1)
   11 continue
      BESJ(NUM+1)=0D0
      BESJ(NUM)=(1D-300,1d-300)
      N=2*NUM+1
      NUM1=NUM-1
      DO I = 1, NUM1
       N=N-2
       I1=NUM-I
       BESj(I1) = n * a * BESj(I1+1) - BESj(I1+2)
         if(cdabs(besJ(i1)).gt.1d300) then
           num = i - 1
           go to 11
         end if
      end do
      N=2*(NUM/2)
      B=1.2533141373155002D0*CDSQRT(A)
      C=B*BESJ(1)
      DO 12 I=3,N,2
      B=B*(I-0.5D0)*(I-2.0D0)/(I-2.5D0)/(I-1.0D0)
   12 C=C+B*BESJ(I)
      C=1.0D0/C
      DO 13 I=1,NUM
   13 BESJ(I)=C*BESJ(I)
      RETURN
      END

c************************************************
C cdcof4a
C
      SUBROUTINE cdcof4a(VL,AD,IA,BD,IB,ID,VN,IN,M,N,C,IER)
      parameter (nterms=330)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 VL,AD,VN,C,bd
      DIMENSION AD(IA),BD(IB),RL(4*nterms)
      COMMON /EPS3/ EPS3
      COMMON /EPS33/ EPS33
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      common /FACT/ FACT(300)
      AC=C
      IER=0
c      eps33 = eps3
      eps4 = 1d-20
c                 eps4 = 1d-40
C* èPOBEPKA BBOÑàMõX áHAóEHàâ èAPAMETPOB
      IF(M.LT.0.OR.N.LT.M) IER=1
      IF(AC.LT.0D0) IER=3
      IF(ID.EQ.1.AND.IA.LE.((N-M)/2+2)) IER = 4
      IF(ID.EQ.0.AND.IN.NE.0) IER = 6
      IF(IER.NE.0) GO TO 900
      aM = M
      aMM = 2d0 * M
      S2=(1D0,0D0)
      bM = 4D0 * am**2 - 1D0
      RC2=C*C
      RC4=RC2*RC2
      N1=N-M
      IF(ID.EQ.0) GO TO 900
C**** BõóàCãEHàE K-TOB PAáãOÜEHàü
      IF(MOD(N1+1,2)) 21,20,21
   20 I11=1
      GO TO 22
   21 I11=0
   22 CONTINUE
      L1=(N1-I11)/2
      L11=L1+1
      RL(L11)=S2
      IF(L1.EQ.0) GO TO 24
C***  BõóàCãEHàE K-TOB èPà R<N-M
      DO 23 J=1,L1
      JJ1=N1-2*J
      RU=S2
      RB=(am+JJ1)*(am+JJ1+1D0)+AK*RC2/2D0*(1D0-bm/(amm+2D0*JJ1-1D0)/
     *   (amm+2D0*JJ1+3D0))        -VL
      RV=S2/RB
      RW=RV
      JJ2=JJ1
      IF(JJ2.LE.1) GO TO 103
  102 CONTINUE
      II=2*M+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(am+JJ2-2D0)*(am+JJ2-1D0)+AK*RC2/2D0*(1D0-bm/(II1-5D0)/
     *   (II1-1D0))        -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV
      IF(CDABS(RV/rw).LE.eps4.and.j.ge.l1/2) GO TO 103
c      IF(CDABS(RV/rw).LE.eps4) GO TO 103
      JJ2=JJ2-2
      IF(JJ2.LE.1) GO TO 103
      GO TO 102
  103 RL(L11-J)=-(amm+JJ1+1D0)*(amm+JJ1+2D0)*AK*RC2/(amm+2D0*JJ1+3D0)/
     *           (amm+2D0*JJ1+5D0)*RW
   23 RL(L11-J)=RL(L11-J)*RL(L11+1-J)
   24 CONTINUE
C***  BõóàCãEHàE K-TOB èPà R>N-M
      IA1=IA-1
      DO 25 J=L11,IA1
      JJ1=N1+2*(J-L11)+2
      II=2*M+JJ1
      II1=II+JJ1
      RU=S2
      RB=(am+JJ1)*(am+JJ1+1D0)+AK*RC2/2D0*(1D0-bm/(II1-1D0)/(II1+3D0))
     *   -VL
      RV=-JJ1*(JJ1-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RW=RV
      JJ2=JJ1
  104 CONTINUE
      JJ2=JJ2+2
      II=2*m+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(am+JJ2)*(am+JJ2+1D0)+AK*RC2/2D0*(1D0-bm/(II1-1D0)/(II1+3D0))
     *          -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV

c        if(n.eq.1.and.dimag(c).gt.0.1d0.and.j.le.3)  write(7,*) j, rw

c       IF(CDABS(RV/rw).LE.eps4.and.j.ge.ia1/2) GO TO 105
       IF(CDABS(RV/rw).LE.eps4) GO TO 105
c       IF(CDABS(RV/rw).LE.eps3) GO TO 105
      GO TO 104
  105 RL(J+1)=RW*(amm+2D0*JJ1-1D0)*(amm+2D0*JJ1+1D0)/((amm+JJ1-1D0)*
     *      (amm+JJ1)*AK*RC2)
      RL(J+1)=RL(J)*RL(J+1)
       IF(J.GT.L11) RL(J+1)=RL(J+1)/RL(L11)
       IF(CDABS(RL(J+1)).LT.eps3) GO TO  125
   25 CONTINUE
      GO TO 126
  125 IJ=J+2
         DO 127 I=IJ,IA
  127    RL(I)=S
  126 CONTINUE

C *                  BõóàCãEHàE K-TA D(N-M)
      L4=(N+M+I11)/2
c      if(l1+1.gt.170) stop 1701
      F2=1D0/Fact(l1+1)
      DO  I = 1, L4
       F2 = F2 * (L4 + I)
       if(dabs(f2).gt.1d300) go to 27
      end do
      iff = 0
      go to 270
   27 continue
      f2 = f2 * 1d-300
      iff = 1
      DO Ij = i, l4
       F2 = F2 * (J3 + Ij)
      end do
 270  continue
      il1 = 1
      if(mod(l1,2).ne.0) il1 = -1
      D1 = iL1 * F2 / 2D0**N1
       if(iff.eq.1) d1 = d1 * 1d300
      SNM=S

      DO 28 J=1,IA
      J1=2*(J-1)+I11
      J2=(J1-I11)/2
      J3=(J1+2*M+I11)/2
      rF2=RL(J)/2D0**J1
      DO I=1,J3
       rF2 = rF2 * (J3 + I)
       if(cdabs(rf2).gt.1d300) go to 30
      end do
      iff = 0
      go to 330
   30 continue
      rf2 = rf2 * 1d-300
      iff = 1
      DO Ij = i, J3
       rF2 = rF2 * (J3 + Ij)
      end do
 330  continue
      rF1=rF2
      DO 29 I=1,J2
   29 rF1=rF1/I
       SG=(-1D0)**J2*rF1
       if(iff.eq.1) SG = sg * 1d300
       IF(CDABS(SG).LT.eps3.AND.J.GT.ia/2) GO TO 128
   28 SNM=SNM+SG
  128 AD(L11)=D1/SNM

      DO 31 I=1,IA
   31 AD(I)=RL(I)*AD(L11)
      IF(ID.NE.2) GO TO 901
C*       BõóàCãEHàE KOùîîàñàEHTOB PAáãOÜEHàü èPà  R<0
      BD(1)=AD(1)
      DO 223 J=2,IB
      JJ1=I11-2*(J-1)
      RU=S2
      RB=(am+JJ1)*(am+JJ1+1D0)+AK*RC2/2D0*(1D0-bm/(amm+2D0*JJ1-1D0)/
     *   (amm+2D0*JJ1+3D0))        -VL
      RV=1D0/RB
      RW=RV
      JJ2=JJ1
  224 CONTINUE
      II=2*M+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(am+JJ2-2D0)*(am+JJ2-1D0)+AK*RC2/2D0*(1D0-bm/(II1-5D0)/
     *   (II1-1D0))        -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV
c          print *, j, RV/rw
c          pause
c      IF(CDABS(RV/rw).LE.eps4.and.j.ge.ib/2) GO TO 225
      IF(CDABS(RV/rw).LE.eps4) GO TO 225
      JJ2=JJ2-2
      GO TO 224
  225 CONTINUE
      IF(JJ1.EQ.-2*M-2.OR.JJ1.EQ.-2*M-1) GO TO 226
      BD(J)=-(amm+JJ1+1D0)*(amm+JJ1+2D0)*AK*RC2/(amm+2D0*JJ1+3D0)/
     *           (amm+2D0*JJ1+5D0)*RW
      GO TO 227
  226 BD(J)= AK*RC2/(amm+2D0*JJ1+3D0)/(amm+2D0*JJ1+5D0)*RW
      IF(I11.EQ.1) BD(J)=-BD(J)
  227 BD(J)=BD(J)*BD(J-1)
          IF(cDABS(BD(J)).LT.eps33) GO TO 1223
  223 CONTINUE
        GO TO 901
 1223 IJ=J+1
        DO 1224 I=IJ,IB
 1224   BD(I)=s
  901 IF(IN.EQ.0) GO TO 900
C*****  BõóàCãEHàE HOPMàPYûôEÉO MHOÜàTEãü *
      VN=S
         iff = 0
         if(cdabs(ad(1)).gt.1d150) iff = 1

      DO 32 J=1,IA
      JJ1=2*(J-1)+I11
         rf1 = 1d0
         if(iff.eq.1) rf1 = 1d-300
      rF1 = rf1 * AD(J) / (amm+2D0*JJ1+1D0) * ad(j)
      IF(M.EQ.0) GO TO 333
      DO 33 I=1,2*M
   33 rF1=rF1*(JJ1+I)
  333 RVN=rF1
      VN=VN+RVN
      IF(CDABS(RVN/(vn+1d-50)).LT.eps4.and.j.gt.ia/2) GO TO 34
c      IF(CDABS(RVN/(vn+1d-50)).LT.eps4) GO TO 34
   32 continue
   34 VN=2D0*VN
         if(iff.eq.1) vn = vn * 1d300
  900 RETURN
      END

c************************************************
c cdrb12cc
c
      SUBROUTINE CDRB12cc(R,DR,AR,ADR,ID,M,N,C,X,AD,IA,IER)
      parameter (nterm=330)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 AD,C,S1,BESJ,BESY,AC,SIG,
     *           RK,R,DR,F2,RS,RD,S,BR,ABR,ar,adr
      DIMENSION AD (IA), BESJ (4*nterm), BESY (4*nterm)
      DIMENSION iiy (4*nterm)
      COMMON /EPS3/ EPS3
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /pi/ pi
      factor=1d250
c      factor=1d+50
      IER=0
      IF(X.LT.1D-10) IER=10
      IF(IER.NE.0) GO TO 900
      M2=2*M
      IF(MOD((N-M),2)) 1,2,1
    1 JK=1
      GO TO 5
    2 JK=0
    5 CONTINUE
      AM=M
      AM=AM/2D0
      PI2 = pi / 2d0
      X2=X*x
      AC=C*X
      SIG=S
      DO 7 I=1,IA
      J=2*(I-1)+JK
      F2=AD(I)
      DO 8 I1=1,M2
    8 F2=F2*(J+I1)
    7 SIG=SIG+F2
      RK=CDSQRT(PI2/AC)/SIG*((X2-AK)/X2)**AM
       NUM=3*CDABS(AC)+10
c       if(x.gt.2d0) NUM = 3 * CDABS(AC) + 10 * (m + x)
       if(x.gt.2d0) NUM = 3 * CDABS(AC) + 6 * (m + x)
      if(NUM.lt.60) num=60
c      if(NUM-1.ge.4*nterm) write(*,*) num, 4*nterm
      if(NUM-1.ge.4*nterm) num = 4*nterm

      IF (ID.EQ.1) CALL CESSEL0(1D0/AC,NUM,BESJ,BESY,iiy)
      IF (ID.EQ.0) CALL BESSJJ0(1D0/AC,NUM,BESJ)

      R=S
      DR=S
      AR=s
      ADR=s
      F1=AK/X*M/(X2-AK)
      NM2=MIN0(IA,NUM/2-2)

      DO 13 I=1,NM2
      J=2*(I-1)+JK
      F2=AD(I)*S1**(J+M-N)
      DO 18 I1=1,M2
   18 F2=F2*(J+I1)
      IF(I.EQ.1) GO TO 21
      IF(CDABS(RD).LT.eps3.AND.CDABS(RS).LT.eps3.AND.I.GT.(M+5)
     *   .AND.ID.EQ.0)   GO TO 19
   21 RS=F2*BESJ(J+M+1)
      RD=F2*((F1+(M+J)/X)*BESJ(J+M+1)-C*BESJ(J+M+2))
      R=R+RS
      DR=DR+RD
   20 CONTINUE
      IF(ID.EQ.0) GO TO 13
      IF(I.EQ.1) GO TO 22
      IF(CDABS(BR).LT.eps3.AND.CDABS(ABR).LT.eps3.AND.I.GT.(M+5))
     *   GO TO 14
   22 BR = BESY(J+M+1) * F2 * factor**(iiy(j+m+1))
                     AR  = AR  + BR
      ABR=(F1 +(M+J)/X)*BESY(J+M+1) * f2* factor**(iiy(j+m+1))
     *    -C*BESY(J+M+2) * f2 * factor**(iiy(j+m+1))
                     ADR = ADR + ABR
   13 CONTINUE
   14 AR= AR* RK
      ADR= ADR* RK
   19 R=R*RK
      DR=DR*RK
  900 RETURN
      END

c************************************************
c CDRF12cc
c
      SUBROUTINE CDRF12cc(R,DR,AR,ADR,I12,M,N,C,AD,IA,BD,IB,IER)
      parameter (nterms=330)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 AD,C,DR,ar,adr,bd
      DIMENSION AD(IA ),BD(IB+1)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /LEG/ P(2*nterms),PD(2*nterms),Q(2*nterms),QD(2*nterms)
      COMMON /FACT/ FACT(300)
      COMMON /eps3/ eps3
       eps4 = 1d-20
      ier=0
      M2=2*M
      aM2 = 2d0 * M
      nm1=n-m
      nm2=n+m
      IF(MOD((N-M),2)) 1,2,1
    1  CONTINUE
      JK=1
      GO TO 36
    2  CONTINUE
      JK=0
   36 CONTINUE
      IF(MOD(N,2)) 3,4,3
    3  CONTINUE
      Jn=1
      GO TO 5
    4  CONTINUE
      Jn=0
    5 CONTINUE
      IF(MOD(m,2)) 33,34,33
   33  CONTINUE
      Jm=1
      GO TO 35
   34  CONTINUE
      Jm=0
   35 CONTINUE
      SIG=S
      DO 7 I=1,IA
      J=2*(I-1)+JK
       rf2=AD(I)*fact2(m2+j,j)
      SIG=SIG+rf2
      if(cdabs(rf2/sig).lt.eps4) go to 8
    7  CONTINUE
    8  CONTINUE
      if((nm1-jk)/2+1.gt.170) stop 1702
      ff3 =  fact3(nm2+jk,(nm2+jk)/2,1d0/fact((nm1-jk)/2+1))
      RK1=(am2+1d0+jk*2d0)/2d0**nm2 * ff3 /
     *    AD(1)/fact(m+1)*sig/C**(M+JK)
      if(k.eq.1) rk1=rk1/(-s1)**(M+JK)*(-s1)**jn
      R=S
      DR=S
      IF(I12.EQ.0) GO TO 70
      RK2=2d0**nm1/(am2-1d0)*fact(m2+1)/
     *    fact(m+1) / ff3
     *    *bd(m+1)*sig/C**(M-1-JK)
      if(jk.eq.1.and.k.eq.0) rk2=-rk2/(am2-3d0)
      if(k.eq.0) go to 339
      if(jk.eq.1.and.jm.eq.0) rk2=rk2/(am2-3d0)
      if(jk.eq.1.and.jm.eq.1) rk2=-rk2/(am2-3d0)*(-s1)
      if(jk.eq.0.and.jm.eq.0) rk2=rk2*s1
      if(k.eq.1) rk2=rk2/(-s1)**(M-1-JK)
  339 continue
      AR=s
      ADR=s
      DO 12 I=1,M
        J2=2*(M-I)+1+JK
      ADR=ADR+BD(I+1)*QD(J2)
   12 AR=AR+BD(I+1)*Q(J2)           
   70 CONTINUE
      DO 13 I=1,IA
        J1=2*(I-1)+1+JK
        J2=2*(M+I)-1+JK
      RS=AD(I)*P(J1)
      RD=AD(I)*PD(J1)
      R=R+RS
      DR=DR+RD
      IF(I12.EQ.0) GO TO 17
      rBR=AD(I)*Q(J2)
      rABR=AD(I)*QD(J2)
      AR=AR+rBR
      ADR=ADR+rABR
      IF(cDABS(rBR/ar).LT.eps4.AND.cDABS(rABR/adr).LT.eps3) GO TO 14
   17 CONTINUE
      IF(CDABS(RS/r).LT.eps4.AND.CDABS(RD/dr).LT.eps4) GO TO 14
   13 CONTINUE
   14 CONTINUE
      IF(I12.EQ.0) GO TO 71
      MR=M+1
      DO 15 I=MR,IB
        J2=2*(I-M)-JK
      if (k.eq.1.and.mod(m,2).eq.0) go to 116
      rBR=BD(I+1)*P(J2)
      rABR=BD(I+1)*PD(J2)
      go to 115
  116 CONTINUE
      rBR=-BD(I+1)*P(J2)
      rABR=-BD(I+1)*PD(J2)
  115 CONTINUE
      AR=AR+rBR
      ADR=ADR+rABR
      IF(cDABS(rBR/ar).LT.eps4.AND.cDABS(rABR/adr).LT.eps4) GO TO 16
   15 CONTINUE
   16 AR=AR/RK2
      ADR=ADR/RK2
   71 R=R/RK1
      DR=DR/RK1
  900 RETURN
      END

c************************************************
C CESSEL0
c
      SUBROUTINE CESSEL0(A,NUM,BESJ,BESY,iiy)
      IMPLICIT COMPLEX*16 (A-H,O-Q,T-Z)
      DIMENSION BESJ(NUM+1),BESY(NUM+1),iiy(num+1)
      factor=1d250
   11 continue
          do i =1,num+1
            iiy(i) = 0
          end do
      BESJ(NUM+1)=(0D0,0D0)
      BESJ(NUM)=(1D-300,1D-300)
      N=2*NUM+1
      NUM1=NUM-1
      DO I=1,NUM1
       N=N-2
       I1=NUM-I
       BESj(I1) = n * a * BESj(I1+1) - BESj(I1+2)
        if(cdabs(besJ(i1)).gt.1d300) then
           num = i - 1
           go to 11
        end if
      end do
      B1=A*BESJ(1)-BESJ(2)
      B2=-A*B1-BESJ(1)
      N=2*(NUM/2)
      B=1.2533141373155002D0*CDSQRT(A)
      C=B*BESJ(1)
      DO 12 I=3,N,2
      B=B*(I-0.5D0)*(I-2.0D0)/(I-2.5D0)/(I-1.0D0)
   12 C=C+B*BESJ(I)
      C=1.0D0/C
      DO 13 I=1,NUM
      BESJ(I)=C*BESJ(I)
   13 continue
c Y
      BESY(1)=-C*B1
      BESY(2)=C*B2
      DO 14 I=3,NUM
      I2=I-2
         if(iiy(i-1).eq.iiy(i2)) then
          BESY(I) = (2.0D0*I2+1.0D0) * a * BESY(I-1) - BESY(I2)
         else
          BESY(I) = (2.0D0*I2+1.0D0) * a * BESY(I-1) - BESY(I2)/factor
         end if
         iiy(i) = iiy(i-1)
         if(cdabs(besY(i)).gt.1d300) then
           besY(i) = besY(i) / factor
           iiy(i) = iiy(i-1) + 1
         end if
   14 continue
      RETURN
      END
c************************************************
c CDRG1cn
c
      SUBROUTINE CDRG1cn(R,DR,AR,ADR,I12,M,n,C,x,AD,IA,BD,IB,rrl,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 C, dr, drr, drr2, dr2, AD(IA), drr3
      COMPLEX*16 ar, adr, bd
      DIMENSION SG(1000), sa(1000), sb(1000), BD(IB+1),
     *          sam(100)
      COMMON /K1/ S, S1, AKSIi, AK, K, NK, nal
      COMMON /FACT/ FACT(300)
      COMMON /eps3/ eps3
      ier=0
c      if(i12.eq.0) go to 900
c! lambda according Slavyanov
      rl = rrl - c**2
c! aksi!
      if(k.eq.0) aksi = x**2 - 1d0
      if(k.eq.1) aksi = x**2 + 1d0
      eps4 = 1d-20
      aM = 0.5d0 * M
      if(k.eq.0) xx0 = ((x - 1d0) / (x + 1d0)) ** am
      if(k.eq.0) xx =  (x - 1d0) / (x + 1d0)

c      print *,k
c      print *,x, xx
c      pause

      rp = - S1 * c
      re = cdexp(-rp * (x - 1d0))
      rp2 = 2d0 * rp
      rrp2l = (rp2 + 1d0) * (m + 1d0) + rl

c coefficients
      M2=2*M
      nm1=n-m
      nm2=n+m
      IF(MOD((N-M),2)) 1,2,1
    1  CONTINUE
      JK=1
      GO TO 36
    2  CONTINUE
      JK=0
   36 CONTINUE

c       if((nm1-jk)/2+1.gt.170) stop 1707
      ff3 = fact3(nm2+jk,(nm2+jk)/2,1d0/fact((nm1-jk)/2+1))
      rg0=(2d0*m+1d0+jk*2d0)/2d0**(nm2+1d0) * ff3 / AD(1)/C**(M+JK)
      rg0 = 1d0 / rg0

c 2===
      IF(I12.ne.0) then
      IF(MOD(m,2)) 33,34,33
   33  CONTINUE
      Jm=1
      GO TO 35
   34  CONTINUE
      Jm=0
   35 CONTINUE
      SIG=S
      DO 7 I=1,IA
      J=2*(I-1)+JK
      rf2 = AD(I) * fact2(m2+j,j)
      SIG=SIG+rf2
      if(cdabs(rf2/sig).lt.eps4) go to 8
    7  CONTINUE
    8  CONTINUE
      RK2=2d0**nm1/(m2-1d0)*fact(m2+1)/fact(m+1) / ff3
     *    * bd(m+1)*sig/C**(M-1-JK)
      if(jk.eq.1.and.k.eq.0) rk2=-rk2/(m2-3d0)
      if(k.eq.0) go to 339
      if(jk.eq.1.and.jm.eq.0) rk2=rk2/(m2-3d0)
      if(jk.eq.1.and.jm.eq.1) rk2=-rk2/(m2-3d0)*(-s1)
      if(jk.eq.0.and.jm.eq.0) rk2=rk2*s1
      if(k.eq.1) rk2=rk2/(-s1)**(M-1-JK)
  339 continue
      end if
c 2===
c------------
      Sssa0 = S
c. 1s
      DO I = 1, M
        j = - 2 * (M - I + 1) + JK
        if(j.ge.-m) go to 351
      sssa0 = sssa0 + BD(m-I+2) * (-1)**(-j-1) * fact(-j)
     *        * fact(m2+j+1)
      end do

351   continue

c. 2s
      rf1 = s
      DO I = 1, IA
      ir = 2* (i - 1) + jk
      rf3 = s
c=
         do j = 1, m
           ff2 = fact(m+1) / fact(m-j+1)
           rf2 = s
           do jj = 1, j+1
            jjj = jj - 1
            rf2 = rf2 + fact2(m,m-jjj)**2 / fact(jjj+1) /
     *            fact(j-jjj+1) * fact2(ir+m2-jjj,ir+jjj)
           end do
           rf2 = rf2 * (-1)**j / (2d0 * j)  * ff2
           rf3 = rf3 + (-1)**(j-1) / (2d0 * j)  * ff2**2 *
     *           fact2(ir+m2-j,ir+j) + rf2
         end do
c=
         irm = (ir + m - 1 + jk) / 2
         rf4 = s
         if(ir.eq.0) go to 295
         do j = 1, irm
           jj = j - 1
           if(2*jj.gt.ir-1) go to 295
           rf4 = rf4 + (2d0 * ir + m2 - 4d0 * jj - 1d0) / (ir + m - jj)
     *           / (2d0 * jj + 1d0) * fact2(ir+m2-2*jj-1,ir-2*jj-1)
         end do
 295   CONTINUE
      rf11 =  AD(I) * (rf3 - rf4)
      rf1 = rf1 + rf11
      if(cdabs(rf11/rf1).lt.eps4) go to 297
      end do
 297   CONTINUE
      Sssa0 = Sssa0 + rf1
      rs2 = rf1

c. 3s
      rf1 = S
      MR = M + 1
      DO I = MR, IB
        J = 2 * I  - JK
        rf11 = bD(I+1) * fact2(j-1,j-m2-1)
        rf1 = rf1 + rf11
      if(cdabs(rf11/rf1).lt.eps4) go to 298
      end do

  298 continue

      Sssa0 = Sssa0 + rf1

c 1--
      SG0 = (1d0,0d0)
      SG(1) = rrp2l / (m + 1d0) * sg0
      SG(2) = ((2d0 * (m + 2d0 + rp2) + rrp2l) * sg(1) -
     *        (m + 1d0) * sg0) / (2d0 * (m + 2d0))
      r = sg0 + sg(1) * xx + sg(2) * xx**2
c      r = (sg0 + sg(1) * xx + sg(2) * xx**2) * re
      drr = sg(1) * 1d0 +  sg(2) * 2d0 * xx
c      drr = (sg(1) * 1d0 +  sg(2) * 2d0 * xx) * re

c 2--
      if(i12.ne.0) then

      do i = 1, 100
         Sam(i) = s
      end do

      Sam(1) = - 2d0 / m / c / rg0
      if(-m+1.eq.0) go to 352
      Sam(2) = (2d0 * m * (rp2 + 1d0) - rrp2l) / (m - 1d0) * sam(1)
      if(-m+2.eq.0) go to 352
       do i = 1, m-2
         j = - m + i
         Sam(i+2) = ((2d0 * j * (j + rp2 + m + 1d0) + rrp2l)
     *              * Sam(i+1) -
     *           (j * (j + m) * sam(i))) / ((j + 1d0) * (j + m + 1d0))
       end do

  352 continue

c a0--
      s2ic = 2d0 * s1 * c
      sssa1 = s
       do i = 1, m
       j = i - 1
       rf = s
         do ii = 1, m - j
           rf = rf + s2ic**ii / fact(ii+1) * fact(m-j) / fact (ii)
     *          / fact(m-ii-j+1)
         end do
       sssa1 = sssa1 + sam(i) * rf
       end do

      Sa00 = 2d0 / fact(m+1) * sssa0 / rk2 - sssa1

c a0--
c!!!!!!!!          Sa0 = sa00

      Sb0 = (-2d0 * (rp2 + m) + rrp2l) * Sam(m) / m
      if(m.gt.1) sb0 = sb0 + (m - 1d0) * sam(m-1) / m

      sb(1) = rrp2l / (m + 1d0) * sb0
      sb(2) = ((2d0 * (m + 2d0 + rp2) + rrp2l) * sb(1) -
     *        (m + 1d0) * sb0) / (2d0 * (m + 2d0))
      Sa(1) = (rrp2l * sa00 - (m + 2d0) * sb(1)
     *        + 2d0 * (rp2 + m + 1d0) * sb0) / (m + 1d0)
      Sa(2) = ((2d0 * (m + 2d0 + rp2) + rrp2l) * sa(1)
     *        - (m + 1d0) * sa00 - ((m + 4d0) * sb(2)
     *        - 2d0 * (rp2 + m + 3d0) * sb(1)
     *        + (m + 2d0) * sb0)) / (2d0 * (m + 2d0))

      r2 = s
      drr2 = s
       do i = 1, m
         r2 = r2 + sam(m-i+1) / xx**i
         drr2 = drr2 - i * sam(m-i+1) / xx**(i+1)
       end do

      r2 = r2 + sa00 + sa(1) * xx + sa(2) * xx**2
      rr2 = sb0 + sb(1) * xx + sb(2) * xx**2
      drr2 = drr2 + sa(1) * 1d0 +
     *       2d0 * sa(2) * xx
      drr3 = sb(1) * 1d0 + 2d0 * sb(2) * xx
c      r2 = r2 + sa00 + sa(1) * xx + sa(2) * xx**2
c      r2 = r2 * re
c      rr2 = (sb0 + sb(1) * xx + sb(2) * xx**2) * re
c      drr2 = drr2 + sa(1) * 1d0 +
c     *       2d0 * sa(2) * xx
c      drr2 = drr2 * re
c      drr3 = (sb(1) * 1d0 + 2d0 * sb(2) * xx) * re
      end if

c      cosn = cdcos(c*(x-1d0))
c      csin = cdsin(c*(x-1d0))

      DO  I = 2, 998
c 1--
      a = (i + 1d0) * (i + m + 1d0)
      rb = 2d0 * i * (i + m + 1d0 + rp2) + rrp2l
      g = i * (i + m)

      sg(i+1) = (rb * sg(i) - g * sg(i-1)) / a
      rs = sg(i+1) * xx**(i+1d0)
c      drs = dreal(rs)
c      airs = dimag(rs)
c      rsn = drs * cosn - airs * csin + s1 * (drs * csin + airs * cosn)

c      if(dreal(rs).ge.0d0) Rr01 = Rr01 + RS
c      if(dreal(rs).lt.0d0) Rr02 = Rr02 + RS
      R = R + RS
c       R = R + RSn

      rd = sg(i+1) * (i + 1d0) * xx**i
c      drn = dreal(rd)
c      airn = dimag(rd)
c      rdn = drn * cosn - airn * csin + s1 * (drn * csin + airn * cosn)

      dRr = dRr + Rd
c      dRr = dRr + Rdn
c 2--
      rs2 = s
      rrs2 = s
      rd21 = s
      rd23 = s

      if(i12.ne.0) then

        if(cdabs(sb(i)).gt.1d300) return

      sb(i+1) = (rb * sb(i) - g * sb(i-1)) / a
      sa(i+1) = (rb * sa(i) - g * sa(i-1) - ((2d0 * i + m + 2d0) *
     *          sb(i+1) - 2d0 * (2d0 * i + rp2 + m + 1d0) * sb(i) +
     *          (2d0 * i + m) * sb(i-1))) / a
      rs2 = sa(i+1) * xx**(i+1d0)
      rrs2 = sb(i+1) * xx**(i+1d0)

c      drs2 = dreal(rs2)
c      airs2 = dimag(rs2)
c      rs2n = drs2 * cosn - airs2 * csin + s1 * (drs2 * csin
c     *       + airs2 * cosn)

c      drs22 = dreal(rrs2)
c      airs22 = dimag(rrs2)
c      rrs2n = drs22 * cosn - airs22 * csin + s1 * (drs22 * csin
c     *       + airs22 * cosn)

      R2 = R2 + RS2
c      R2 = R2 + RS2n
      Rr2 = Rr2 + RrS2
c      Rr2 = Rr2 + RrS2n

      rd21 = sa(i+1) * (i + 1d0) * xx**i
      rd23 = sb(i+1) * (i + 1d0) * xx**i
c      drd21 = dreal(rd21)
c      aird21 = dimag(rd21)
c      rd21n = drd21 * cosn - aird21 * csin + s1 * (drd21 * csin
c     *       + aird21 * cosn)
c      drd23 = dreal(rd23)
c      aird23 = dimag(rd23)
c      rd23n = drd23 * cosn - aird23 * csin + s1 * (drd23 * csin
c     *       + aird23 * cosn)

      dRr2 = dRr2 + rd21
      dRr3 = dRr3 + rd23
c      dRr2 = dRr2 + rd21n
c      dRr3 = dRr3 + rd23n

c      if(CDABS(RS2n/r2).LT.eps4.AND.CDABS(Rrs2n/rr2).LT.eps4) GO TO 15
      if(CDABS(RS2/r2).LT.eps4.AND.CDABS(Rrs2/rr2).LT.eps4) GO TO 15
      end if
   15 CONTINUE
c      IF(CDABS(RS/r).LT.eps4.AND.CDABS(RDn/drr).LT.eps4) GO TO 14
      IF(CDABS(RS/r).LT.eps4.AND.CDABS(RD/drr).LT.eps4) GO TO 14
      end do
   14 CONTINUE
c 1--

c       R = R + rr01 + rr02

      R = R * xx0 / (x + 1d0) * re * rg0

c      R = R * xx0 / (x + 1d0) * rg0

c      aaR = R
c      R = aaR
      dr = r * (-rp  - 1d0/(x + 1d0) + m / aksi) +
     *     DRr * (2d0 / (x + 1d0)**3) * xx0 * re * rg0
c      dr = r * (-rp  - 1d0/(x + 1d0) + m / aksi) +
c     *     DRr * (2d0 / (x + 1d0)**3) * xx0 * rg0
c      aaR = dR
c      dR = aaR
c 2--
      if(i12.ne.0) then
      R2 = (rr2 * dlog(xx) + R2) * xx0 / (x + 1d0) * re
c      R2 = (rr2 * dlog(xx) + R2) * xx0 / (x + 1d0)
      dr2 = r2 * (-rp -1d0/(x + 1d0) + m / aksi) +
     *     (DRr2 + rr2 / xx + drr3 * dlog(xx)) *
     *     (2d0 / (x + 1d0)**3) * xx0 * re
c      dr2 = r2 * (-rp -1d0/(x + 1d0) + m / aksi) +
c     *     (DRr2 + rr2 / xx + drr3 * dlog(xx)) *
c     *     (2d0 / (x + 1d0)**3) * xx0
      AR = r2
      ADR = dr2
      end if

      if(dimag(c).eq.0d0) then
      r = dreal(r)
      dr = dreal(dr)
      ar = dreal(ar)
      adr = dreal(adr)
      end if

  900 RETURN
      END
c************************************************
C cdrkf45a
c
      SUBROUTINE cDRKF45a(F,NEQN,Y,T,TOUT,
c     1   RELERR,ABSERR,IFLAG,WORK,IWORK)
     1   RELERR,ABSERR,IFLAG)
c      INTEGER IWORK(5),NEQN
      INTEGER NEQN
c      REAL*8 Y(NEQN),T,TOUT,RELERR,ABSERR,
c     1   WORK(27)
      REAL*8 T,TOUT,RELERR,ABSERR
c      complex*16 Y(NEQN), WORK(27)
      complex*16 Y(NEQN)
      EXTERNAL F
c      K1M = NEQN + 1
c      K1 = K1M +1
c      K2 = K1 + NEQN
c      K3 = K2 + NEQN
c      K4 = K3 + NEQN
c      K5 = K4 + NEQN
c      K6=K5+NEQN
c      CALL cDRKFS (F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK(1),
c     1   WORK(K1M),WORK(K1),WORK(K2),WORK(K3),WORK(K4),WORK(K5),
c     2   WORK(K6),WORK(K6+1),IWORK(1),IWORK(2),IWORK(3),
c     3   IWORK(4),IWORK(5))
      CALL cDRKFSa (F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG)
      RETURN
      END

c************************************************
C cdrkfsa
c
      SUBROUTINE cDRKFSa(F,NEQN,Y,T,TOUT,RELERR,ABSERR,
     1   IFLAG)
      LOGICAL  HFAILD, OUTPUT
c      complex*16 Y(NEQN),YP(NEQN),F1(NEQN),
c     *           F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN)
      complex*16 Y(NEQN),YP(4),F1(4),
     *           F2(4),F3(4),F4(4),F5(4)
      REAL*8 T,TOUT,RELERR,ABSERR,h
c      REAL*8 Y(NEQN),T,TOUT,RELERR,ABSERR,H,YP(NEQN),F1(NEQN)
     1   ,SAVRE,SAVAE
      EXTERNAL F
      REAL*8  A,AE,DT,EE,EEOET,ESTTOL,ET,HMIN,REMIN,
     1   RER,S,SCALE,TOL,TOLN,U26,EPSP1,EPS,YPK
      DATA REMIN /1.0D-16/
      DATA MAXNFE /500000/
      IF(NEQN.LT.1) GO TO 10
      IF((RELERR.LT.0D0).OR.(ABSERR.LT.0D0)) GO TO 10
      MFLAG = IABS(IFLAG)
      IF((MFLAG.EQ.0).OR.(MFLAG.GT.8))  GO TO 10
      IF(MFLAG.NE.1) GO TO 20
      EPS = 1.0D0
    5 EPS = EPS/2.0D0
      EPSP1 = EPS + 1.0D0
      IF(EPSP1.GT.1.0D0) GO TO 5
      U26 = 26.0D0 * EPS
      GO TO 50
   10 IFLAG =8
      RETURN
   20 IF((T.EQ.TOUT).AND.(KFLAG.NE.3)) GO TO 10
      IF(MFLAG.NE.2) GO TO 25
      IF((KFLAG.EQ.3).OR.(INIT.EQ.0)) GO TO 45
      IF(KFLAG.EQ.4) GO TO 40
      IF((KFLAG.EQ.5).AND.(ABSERR.EQ.0.0D0)) GO TO 30
      IF((KFLAG.EQ.6).AND.(RELERR.LE.SAVRE).AND.
     1   (ABSERR.LE.SAVAE)) GO TO 30
      GO TO 50
   25 IF(IFLAG.EQ.3) GO TO 45
      IF(IFLAG.EQ.4) GO TO 40
      IF((IFLAG.EQ.5).AND.(ABSERR.GT.0.0D0))  GO TO 45
   30 STOP
   40 NFE = 0
      IF(MFLAG.EQ.2) GO TO 50
   45 IFLAG = JFLAG
      IF(KFLAG.EQ.3)  MFLAG = IABS(IFLAG)
   50 JFLAG = IFLAG
      KFLAG = 0
      SAVRE = RELERR
      SAVAE = ABSERR
      RER = 2.0D0*EPS + REMIN
      IF(RELERR.GE.RER)  GO TO 55
      RELERR = RER
      IFLAG = 3
      KFLAG = 3
      RETURN
   55 DT = TOUT -T
      IF(MFLAG.EQ.1)  GO TO 60
      IF(INIT.EQ.0)  GO TO 65
      GO TO 80
   60 INIT = 0
      KOP = 0
      A = T
      CALL F(A,Y,YP)
      NFE = 1
      IF(T.NE.TOUT) GO TO 65
      IFLAG = 2
      RETURN
   65 INIT = 1
      H = DABS(DT)
      TOLN = 0.0D0
      DO 70 K = 1,NEQN
c         TOL=RELERR*DABS(Y(K)) + ABSERR
         TOL=RELERR*cDABS(Y(K)) + ABSERR
          IF(TOL.LE.0.0D0)  GO TO 70
         TOLN = TOL
c         YPK = DABS(YP(K))
         YPK = cDABS(YP(K))
         IF(YPK *H**5.GT.TOL)  H=(TOL/YPK)**0.2D0
   70 CONTINUE
      IF(TOLN.LE.0.0D0)  H = 0.0D0
      H=1D-4
      H= DMAX1(H, U26*DMAX1(DABS(T),DABS(DT)))
      JFLAG = ISIGN(2,IFLAG)
   80 H =DSIGN(H,DT)
      IF(DABS(H).GE.2.0D0*DABS(DT))  KOP = KOP+1
      IF(KOP.NE.100) GO TO 85
      KOP = 0
      IFLAG = 7
      RETURN
   85 IF(DABS(DT).GT.U26*DABS(T))  GO TO 95
      DO 90  K = 1,NEQN
   90 Y(K) = Y(K) + DT*YP(K)
      A = TOUT
      CALL F(A,Y,YP)
      NFE = NFE+ 1
      GO TO  300
   95 OUTPUT=.FALSE.
      SCALE= 2D0/RELERR
      AE = SCALE * ABSERR
  100 HFAILD = .FALSE.
      HMIN = U26 * DABS(T)
      DT = TOUT -T
      IF(DABS(DT).GE.2.0D0*DABS(H))  GO TO 200
      IF(DABS(DT).GT.DABS(H))  GO TO 150
      OUTPUT= .TRUE.
      H = DT
      GO TO  200
  150 H = 0.5D0 * DT
  200 IF(NFE.LE.MAXNFE)  GO TO 220
      IFLAG = 4
      KFLAG = 4
      RETURN
  220 CALL cFEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,F1)
      NFE = NFE + 5
      EEOET=0D0
      DO 250  K = 1,NEQN
c         ET = DABS(Y(K))+ DABS(F1(K)) + AE
         ET = cDABS(Y(K))+ cDABS(F1(K)) + AE
         IF(ET.GT.0.0D0)  GO TO 240
               IFLAG =5
               RETURN
c  240    EE = DABS((-2.09D3 * YP(K) + (2.197D4 * F3(K) -
c     1       1.5048D4 * F4(K))) + (2.2528D4 *F2(K) -2.736D4 *F5(K)))
  240    EE = cDABS((-2.09D3 * YP(K) + (2.197D4 * F3(K) -
     1       1.5048D4 * F4(K))) + (2.2528D4 *F2(K) -2.736D4 *F5(K)))
  250    EEOET = DMAX1(EEOET,EE/ET)
      ESTTOL = DABS(H)*EEOET*SCALE/7.524D5
      IF (ESTTOL.LE.1.0D0)  GO TO 260
      HFAILD = .TRUE.
      OUTPUT = .FALSE.
      S = 0.1D0
      IF(ESTTOL.LT.5.9049D4)  S = 0.9D0/ESTTOL**0.2D0
      H = S * H
      IF(DABS(H).GT.HMIN)  GO TO  200
      IFLAG = 6
      KFLAG = 6
      RETURN
  260 T=T+H
      DO 270  K=1,NEQN
  270    Y(K) =F1(K)
      A = T
      CALL F(A,Y,YP)
      NFE = NFE + 1
      S = 5.0D0
      IF(ESTTOL.GT.1.889568D-4) S = 0.9D0/ESTTOL**0.2D0
      IF (HFAILD)  S = DMIN1(S,1.0D0)
      H= DSIGN( DMAX1(S*DABS(H),HMIN),H)
      IF(OUTPUT)  GO TO 300
      IF(IFLAG.GT.0) GO TO 100
      IFLAG = -2
      RETURN
  300 T = TOUT
      IFLAG = 2
      RETURN
      END

c************************************************
C cfehl
c
      SUBROUTINE cFEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,S)
c      REAL*8  Y(NEQN),T,H,YP(NEQN),F1(NEQN),F2(NEQN),
c     1   F3(NEQN),F4(NEQN),F5(NEQN),S(NEQN),CH
      REAL*8  T,H,CH
      complex*16  Y(NEQN),YP(NEQN),F1(NEQN),S(NEQN),
     *            F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN)
      CH = H/4.0D0
      DO 221  K = 1,NEQN
  221    F5(K) = Y(K) + CH * YP(K)
      CALL F(T+CH,F5,F1)
      CH = 3.0D0 * H/3.2D1
      DO 222  K = 1,NEQN
  222    F5(K) = Y(K) + CH*(YP(K)+3.0D0*F1(K))
      CALL F(T+3.0D0*H/8.0D0,F5,F2)
      CH = H/2.197D3
      DO 223  K =1,NEQN
  223    F5(K) = Y(K)+CH*(1.932D3*YP(K)+(7.296D3*F2(K)-7.2D3*F1(K)))
      CALL F(T+1.2D1*H/1.3D1,F5,F3)
      CH = H/4.104D3
      DO 224  K= 1,NEQN
  224    F5(K)=Y(K)+CH*((8.341D3*YP(K)-8.45D2*F3(K))
     1         +(2.944D4*F2(K) - 3.2832D4*F1(K)))
      CALL F(T+H,F5,F4)
      CH = H/2.052D4
      DO 225 K= 1,NEQN
  225    F1(K)=Y(K)+CH*((-6.08D3*YP(K)+(9.295D3 *F3(K)-
     1         5.643D3*F4(K)))+(4.104D4*F1(K)-2.8352D4 *F2(K)))
      CALL F(T+H/2.0D0,F1,F5)
      CH = H/7.61805D6
      DO 230  K= 1,NEQN
  230    S(K)=Y(K)+CH*((9.0288D5*YP(K)+(3.855735D6 *F3(K) -
     1        1.371249D6*F4(K)))+(3.953664D6*F2(K)+2.7702D5*F5(K)))
      RETURN
      end
c************************************************
C cDRSF20
C
      SUBROUTINE cDRSF20(R,DR,X,R1,R2,m,c,rl,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C,S1,S,RRA,RRB,R1,R2,R3,RL,cc,rll
c      COMPLEX*16 r, dr, WORK(27),Y(4)
      COMPLEX*16 r, dr, Y(4)
c      DIMENSION IWORK(5)
      COMMON /K1/ S, S1, AKSI, AK, kkk, nk, nal
      COMMON /F2/ RLl, Cc, Mm, i12a
      EXTERNAL  cDRFP2
c   11 FORMAT(1X,I10,2D10.3,4D20.8)
c   31 FORMAT(17H TOLERANCES RESET,2D12.3)
c   41 FORMAT(11H MANY STEPS)
c   71 FORMAT(12H MUCH OUTPUT)
c   81 FORMAT(14H IMPROPER CALL)
          mm = m
          cc = c
          rll = rl
      NEQN = 4
      PL = dreal(RL)
      AC = dreal(C)


c      print *, rL, C
c      print *, 1D0+PL/AC**2, 1D0-PL/AC**2
c      pause

      TKSI=DSQRT((1D0+PL/AC**2)/2D0+DSQRT((1D0-PL/AC**2)**2/4D0+
     *     (M**2-1D0)/AC**2))
      IF(M.EQ.1) TKSI=1D0
       T=TKSI+X/AC*10D0
*       ELERR=1D-5
c      ELERR=1D-12
      ELERR=1D-08
      ABSERR=0D0
      TFINAL=X
      IFLAG=1
      TOUT=T
      CALL cDRFP0(M,C,RL,T,Y)
c  10 CALL cDRKF45(cDRFP2,NEQN,Y,T,TOUT,ELERR,ABSERR,IFLAG,WORK,IWORK)
   10 CALL cDRKF45a(cDRFP2,NEQN,Y,T,TOUT,ELERR,ABSERR,IFLAG)
c      RRA=DCMPLX(Y(1),Y(2))
c      RRB=DCMPLX(Y(3),Y(4))
      RRA = Y(1) + Y(2)
      RRB = Y(3) + Y(4)
      R3=-S1/(C*AKSI*((RRB/RRA+X/AKSI)*R1+R2))
       R=-S1*R3
       DR=-(X/AKSI+RRB/RRA)*R+R1*S1*RRB/RRA
      GO TO (80,20,30,40,50,60,70,80),IFLAG
   20 TOUT=X
      IF(T .GT.TFINAL) GO TO 10
      IER=0
      GO TO  111
c   30 WRITE(7,31)  ELERR,ABSERR
   30 IER=3
      GO TO 10
c   40 WRITE(7,41)
   40 IER = 4
      GO TO 10
   50 ABSERR=1D-9
c      WRITE(7,31)  ELERR,ABSERR
      IER=5
      GO TO 10
   60  ELERR=10D0* ELERR
c      WRITE(7,31)  ELERR,ABSERR
      IFLAG=2
      IER=6
      GO TO 10
c   70 WRITE(7,71)
   70 IFLAG=2
      IER=7
      GO TO 10
c   80 WRITE(7,81)
   80 IER=8
      GO TO 1000
  111 CONTINUE
 1000 RETURN
      END
c************************************************
C cDRFP0
C
      SUBROUTINE cDRFP0(M,C,RL,KSI,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 RL,C, y, c2, gksi, al, an
      REAL*8 KSI
      DIMENSION Y(4)
      aM2 = M * M
      C2=C**2
      AKSI=KSI**2-1D0
      GKSI=C2-(RL-C2)/AKSI+(1D0-aM2)/AKSI**2
      AL=KSI*(2D0*(1D0-aM2)-(RL-C2)*AKSI)/(2D0*GKSI*AKSI**3)
      AN=cDSQRT(1D0+cdabs(GKSI)+AL**2)
      Y(1)=1D0/AN
      Y(2)=0D0
      Y(3)=-AL/AN
      Y(4)=-SQRT(cdabs(GKSI))/AN
      RETURN
      END
c************************************************
C cDRFP2
c
      SUBROUTINE cDRFP2(X,Y,YP)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C,RL, Y(4),YP(4), c2, gksi, akk
      COMMON /F2/ RL,C,M, i1
      AKSI=X**2-1D0
      C2=C**2
      GKSI=C2-(RL-C2)/AKSI+(1D0-M**2)/AKSI**2
      AKK=(Y(1)*Y(3)+Y(2)*Y(4)*(1D0-GKSI)/(Y(1)**2+Y(2)**2+Y(3)**2+
     *     Y(4)**2))
       YP(1)=AKK*Y(1)-Y(3)
       YP(2)=AKK*Y(2)-Y(4)
       YP(3)=AKK*Y(3)+GKSI*Y(1)
       YP(4)=AKK*Y(4)+GKSI*Y(2)
      RETURN
      END
c************************************************
C cDRSF212
C
      SUBROUTINE cDRSF212(R,DR,ar,adr,i12,m,N,c,X,rl,AD,IA,BD,IB,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C,S1,S,AD(IA),RL,cc,rll
c      COMPLEX*16 WORK(27),Y(4),BD(IB), ar, adr, r, dr
      COMPLEX*16 Y(4),BD(IB), ar, adr, r, dr
c      DIMENSION IWORK(5)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /F2/ RLl,Cc,Mm,i12a
      EXTERNAL  cDRF222
c   11 FORMAT(1X,I10,2D10.3,4D20.8)
c   31 FORMAT(17H TOLERANCES RESET,2D12.3)
c   41 FORMAT(11H MANY STEPS)
c   71 FORMAT(12H MUCH OUTPUT)
c   81 FORMAT(14H IMPROPER CALL)

          mm = m
          cc = c
          rll = rl
          i12a = i12

      NEQN=2
      if(i12.eq.2) NEQN=4
       T=0D0
cccc       ELERR=1D-5
       ELERR=1D-15
      ABSERR=0D0
      TFINAL=X
      IFLAG=1
      TOUT=T
      CALL cDRF2K02(M,N,C,AD,IA,BD,IB,i12,Y)
c   10 CALL cDRKF45(cDRF222,NEQN,Y,T,TOUT,ELERR,ABSERR,IFLAG,WORK,IWORK)
   10 CALL cDRKF45a(cDRF222,NEQN,Y,T,TOUT,ELERR,ABSERR,IFLAG)

      if(i12.eq.0) then
      R=Y(1)
      DR=Y(2)
      end if

      if(i12.eq.1) then
      aR=Y(1)
      aDR=Y(2)
      end if

      if(i12.eq.2) then
      R=Y(1)
      DR=Y(2)
      aR=Y(3)
      aDR=Y(4)
      end if

      GO TO (80,20,30,40,50,60,70,80),IFLAG
   20 TOUT=X
      IF(T.LT.TFINAL) GO TO 10
      IER=0
      GO TO  111
c   30 WRITE(7,31)  ELERR,ABSERR
   30 IER=3
      GO TO 10
c   40 WRITE(7,41)
   40 IER = 4
      GO TO 10
   50 ABSERR=1D-9
c      WRITE(7,31)  ELERR,ABSERR
      IER=5
      GO TO 10
   60  ELERR=10D0 * ELERR
c      WRITE(7,31)  ELERR,ABSERR
      IFLAG=2
      IER=6
      GO TO 10
c   70 WRITE(7,71)
   70 IFLAG=2
      IER=7
      GO TO 10
c   80 WRITE(7,81)
   80 IER=8
      GO TO 1000
  111 CONTINUE
 1000 RETURN
      END

c************************************************
C cDRF222
c
      SUBROUTINE cDRF222(X,Y,YP)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C, RL, Y(4),YP(4), ar
      COMMON /F2/ RL,C,M,i12
       AKSI=X**2+1D0
       am = m * m
       ar = RL - (C*X)**2 - aM / AKSI

      if(i12.eq.0.or.i12.eq.1) then
       YP(1)=Y(2)
       YP(2)=(-2D0*X*Y(2)+ aR * Y(1))/AKSI
      end if

      if(i12.eq.2) then
       YP(1)=Y(2)
       YP(2)=(-2D0*X*Y(2) + aR * Y(1))/AKSI
       YP(3)=Y(4)
       YP(4)=(-2D0*X*Y(4) + aR * Y(3))/AKSI
      end if

      RETURN
      END

c************************************************
C cDRF2K02
C
      SUBROUTINE cDRF2K02(M,N,C,AD,IA,BD,IB,i12,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 AD,C,R1,S1,r,s,s2,t,r2,ss,ss2, bd, y
      DIMENSION AD(IA),BD(IB),FACT(300),Y(4)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /FACT/ FACT
      COMMON /eps3/ eps3
      COMMON /pi/ pi
       eps4 = 1d-20
      M2=2*M
      IF(MOD((N-M),2)) 1,2,1
    1 JK=1
      GO TO 3
    2 JK=0
    3 CONTINUE
      S2 = s
      DO 4 I=1,IA
      J=2*(I-1)+JK
      F=1D0
      DO 5 JJ=1,M2
    5 F=F*(J+JJ)
      T=AD(I)*F
      S2=S2+T
      IF(cDABS(T/s2).LT.eps4) GO TO 6
    4 CONTINUE
    6 CONTINUE
      MF=N+M+1+JK
      F=FACT(MF)
c      IF(MF.LE.170) F=FACT(MF)
c      IF(MF.gt.170) stop 1706

      if(i12.eq.1.or.i12.eq.2) then
c----
       jmn = mod(N-m-jk,4)
       if(jmn.eq.0) ss = 1d0
       if(jmn.eq.1) ss = s1
       if(jmn.eq.2) ss = -1d0
       if(jmn.eq.3) ss = -s1
c       jmn2 = dabs(mod(m-n-jk,4))
       jmn2 = iabs(mod(m-n-jk,4))
       if(jmn2.eq.0) ss2 = 1d0
       if(jmn2.eq.1) ss2 = s1
       if(jmn2.eq.2) ss2 = -1d0
       if(jmn2.eq.3) ss2 = -s1
c-----
      R=Ss*(M2-1D0)*FACT(M+1)*PI/2D0**(2*N-M+1)*C**(M-1-JK)/
     *  FACT(M2+1)/BD(M+1)/S2*
     *  (F/FACT((N-M-JK)/2+1)/FACT((N+M+JK)/2+1))**2
      R1=Ss2*S2/2D0**M/C**(M+1+JK)/AD(1)/FACT(M+1)
      end if

      if(i12.eq.0.or.i12.eq.2) then
          R2=Ss*2D0**M*C**(M+JK)*AD(1)*FACT(M+1)/s2
      end if

      IF(JK.EQ.1) GO TO 10

      if(i12.eq.0) then
      Y(1)=R2/(m2+1d0)
      Y(2)=0D0
      end if

      if(i12.eq.1) then
      Y(1)=R
      Y(2)=R1*(M2+1D0)
      end if

      if(i12.eq.2) then
      Y(1)=R2/(m2+1d0)
      Y(2)=0D0
      Y(3)=R
      Y(4)=R1*(M2+1D0)
      end if

      GO TO 900

   10 CONTINUE

      if(i12.eq.0) then
      Y(1)=0D0
      Y(2)=R2/(M2+3D0)
      end if

      if(i12.eq.1) then
      Y(1)=R1*(M2+3D0)
      Y(2)=R*(M2-3D0)
      end if

      if(i12.eq.2) then
      Y(1)=0D0
      Y(2)=R2/(M2+3D0)
      Y(3)=R1*(M2+3D0)
      Y(4)=R*(M2-3D0)
      end if

  900 CONTINUE

c      write(*,*) n
c      write(*,22) y(1)
c      write(*,22) y(2)
c      write(7,22) y(1)
c      write(7,22) y(2)
c      write(*,22) y(3)
c      pause
c      write(*,22) y(4)
c      write(7,22) y(3)
c      write(7,22) y(4)
c  22  format(1x,2d24.15)

      RETURN
      END
C**********************************************************************
C CDLAMo_l - Oblate eigenvalues (large values)
C**********************************************************************
C
      SUBROUTINE CDLAMo_l(VL,M,NN,C,EPS,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
c      COMPLEX*16 VL,C,dx,alam0
      COMPLEX*16 VL,C,dx
      DIMENSION VL(NN),RL(1550),SL(1550),RDA(1550)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /FACT/ FACT(300)

*****************************************
c      COMMON /r0/ alam0(500)
*****************************************
      inull = 0
 991  continue
      eps1=eps
      AC=C
      CC=CDABS(C)
      AF=0D0
      IF(DABS(AC-CC).GT.1D-10) AF=0.1D0
      SC=DCMPLX(0.1D0,AF)
      IER=0
c -------------------------------------
C**  Parameters check-up
c -------------------------------------
      IF(AC.LT.0D0) IER=3
      IF(IER.NE.0) GO TO 900
      DO 112 IN=1,NN
 112      vl(in)=s
      DO 113 IN = 1,550
          rl(in) = s
          sl(in) = s
 113      rda(in) = s

      sdelta=s

      M2=M*M
      aM =  M
      aMM = 2d0 * M
      AM2 = 4D0 * M2 - 1D0
      RC2=C*C
      RC4=RC2*RC2
cc      NC=idint(aC)
             jcn=0
             DO 111 IN=1,NN

c       print *, in
c       pause


              inout = 4
c              if(ac.gt.75d0)  inout = 8
              nt = 0
              iv = 0

*mmm             if(ac.gt.27d0.and.m.gt.10) sc=-sc
c             dx=dcmplx(0.1d0,0.1d0*af)
             dx=dcmplx(0.25d0,0.25d0*af)
c             if(ac.gt.70d0)  dx = dcmplx(0.3d0,0.3d0*af)
c             if(ac.gt.25d0) dx=dcmplx(0.25d0,0.25d0*af)
             if(ac*m.gt.400d0)
     *          dx=dcmplx(0.35d0+(ac*m-400d0)/1000d0,0.35d0*af)

             if(in.gt.6.and.dreal(vl(in-1)).lt.-4*ac) then
             dx=dcmplx(0.1d0,0.1d0*af)
             if(ac.gt.60d0) dx=dcmplx(0.25d0,0.25d0*af)
             if(ac.gt.70d0) dx=dcmplx(0.3d0,0.3d0*af)
             end if

             if(in.gt.6.and.ac.gt.40d0.and.dreal(vl(in-2)).lt.0d0.and.
     *          dreal(vl(in-1)).gt.-4*ac)
     *          dx=dcmplx(0.35d0+(ac-40d0)/200d0,0.35d0*af)

             if(in.gt.6.and.ac.gt.70d0.and.
     *          dabs(dreal(vl(in-2)) + dreal(vl(in-1))).lt.cn)
     *          dx=dcmplx(0.35d0,0.35d0*af)


c       print *, '1'
c       pause

             n1001=0
             if(dreal(sc).lt.0d0.and.n1001.eq.0) sc = -sc
      N=M+IN-1
      N1=N-M
      NM2=N-M+2
      NMAX=NM2+2*CC+2
 1113 continue
             JN=0
c     IF(NMAX.GT.500) WRITE(7,2000) NMAX
c 2000 FORMAT(14X,'******* CDLAMo_l ****** NMAX=',I4)
      IF(MOD((NMAX-NM2),2))250,251,250
  250 NMAX=NMAX+1
  251 CONTINUE
      IF(aC.LT.5.50D0) GO TO 401

      IF(n.lT.7.or.in.le.3) GO TO 401
      IF(in.lT.7) GO TO 401
      IF(dreal(vl(in-1)).lt.0d0.and.m.ge.10) GO TO 401

      IF(dreal(vl(in-1)).lt.3d0*ac.and.m.ge.15) then
      if(dreal(sc).lt.0d0.and.n1001.eq.0) sc = - sc
      GO TO 401
      end if

        cn=ac
        an=dfloat(n)
        if(an.lt.cn) cn=an

                   IF(dabs(dreal(vl(in-1))).lT.cn.
     *                and.dreal(vl(in-2)).lT.cn.
     *                or.dabs(dreal(vl(in-2)+vl(in-1))).lT.ac+m) then
           IF(dreal(vl(in-1)).lT.0d0)  RL(1) = S
           IF(dreal(vl(in-1)).gT.0d0)  RL(1) = 2d0*vl(in-1)
           IF(dreal(vl(in-1)).gT.cn)   RL(1) = vl(in-1) + ac/10d0
             r0=rl(1)
c             if(mod(n-m+1,2).ne.0) nt = 1
             if(n1001.gt.1.and.m.lt.5) dx = dx / 1.25d0
             if(n1001.ne.0.and.m.ge.5)  dx = 2d0 * dx
ccccc             if(n1001.eq.0.and.dreal(dx).lt.ac/100d0)  dx = ac / 100d0
ccccc             if(dreal(sc).gt.0d0.and.ac.lt.60d0.and.m.lt.5) sc = - sc

             if(n1001.eq.0) then
             if(dreal(dx).lt.ac/100d0)  dx = ac / 100d0
             if(dreal(sc).gt.0d0.and.ac.lt.60d0.and.m.lt.5) sc = - sc
             end if

c             write(*,*) in,r0

*             print *, n, '   0'
c             write(7,*) n, '   0'
             iv = 100

             GO TO 2
                    end if

         IF(dreal(vl(in-1)).lt.cn.
     *      and.dabs(dreal(vl(in-2)+vl(in-1))).lT.cn.
     *      or.dreal(vl(in-1)).lt.cn.
     *      and.dabs(dreal(vl(in-2)+vl(in-3))).lT.cn) then
            nt = 2
            GO TO 401
         end if


c~~~~~~~~~~~~~~~~~
         IF(dreal(vl(in-1)).gt.cn.and.dreal(vl(in-1)+vl(in-2)).gt.ac+m.
     *      or.dreal(vl(in-2)).gt.0d0) then

c         IF(dreal(dx).gt.0.1d0.and.dabs(VL(IN-1)-VL(IN-2)).gt.cn+2*m+5)

                 IF(dreal(dx).gt.0.1d0.and.
     *           cdabs(VL(IN-1)-VL(IN-2)).gt.ac+2*m+5)
     *           dx = dcmplx(0.1d0,0.1d0*af)


c===============
       if(cdabs(VL(IN-2)-VL(IN-3)).lt.cn.or.
     *    cdabs(VL(IN-1)-VL(IN-2)).lt.cn) then

       if(mod(n-m+1,2).eq.0) r0 = VL(IN-1) + cdabs(VL(IN-2)-VL(IN-3))
       if(mod(n-m+1,2).ne.0) r0 = VL(IN-2) + 2d0*ac
*             print *, n, '   ++++'
c            write (7,*) n, '   ++++'
             iv = 1111
c----------------
             if(n1001.eq.0.or.n1001.eq.1) then
             dx = dcmplx(0.25d0+ac*m/10000d0,0.25d0*af)
             if(ac.gt.60d0) dx = dcmplx(0.3d0+ac*m/10000d0,0.3d0*af)
             if(ac.gt.75d0) dx = dcmplx(0.4d0+ac*m/10000d0,0.4d0*af)
             end if
c----------------

       go to 223
       end if
c===============

       r0=VL(IN-1)+cdabs(VL(IN-1)-VL(IN-2))

c        write(*,*) in,r0

c^^^^^^^^^^^^^^^
      if(cdabs(VL(IN-1)-VL(IN-2)).lt.10d0-m/15d0) then

c----------------
         if(ac.ge.40d0) then
         if(ac*m.ge.600d0.or.dreal(vl(in-1)).gt.cn+2*m)
     *      r0 = r0 + VL(IN-2)/2d0
         if(ac*m.lt.600d0.or.dreal(vl(in-1)).lt.cn+2*m)
     *      r0 = r0 + VL(IN-2)
         end if
c----------------

c----------------
          if(ac.lt.40d0.and.m.le.15) then
          if(dreal(vl(in-1)).gt.ac+m) r0 = r0 + VL(IN-2)/2d0
          if(dreal(vl(in-1)).le.ac+m) r0 = r0 + VL(IN-2)
          end if
c----------------

      end if
c^^^^^^^^^^^^^^^

 223    continue
      if(ac.ge.30d0.and.m.gt.3) r0=r0+m/5d0
      if(ac.ge.16d0.and.ac.lt.29d0) r0=r0+2.086d0
      if(ac.ge.30d0.and.ac.lt.40d0) r0=r0+(ac-30d0)/2d0
      if(ac.ge.30d0.and.ac.lt.40d0.and.dreal(vl(in-2)).lt.cn) r0=r0-m
      if(ac.ge.40d0.and.ac.lt.50d0)  r0=r0+n/10d0-(ac-40d0)/3d0-m/2.5d0
c      if(ac.ge.40d0.and.ac.lt.50d0)  r0=r0+n/10d0-(ac-40d0)/3d0
      if(ac.ge.50d0.and.ac.le.60d0) r0=r0+(ac-50d0)/3d0
      if(ac.gt.60d0.and.ac.le.70d0) r0=r0+(ac-60d0)/3d0
      if(ac.gt.70d0.and.ac.le.80d0) r0=r0+(ac-70d0)/4d0
      if(ac.gt.80d0.and.ac.le.90d0) r0=r0-(ac-80d0)/8d0+3.5d0
      if(ac.gt.90d0.and.ac.le.100d0) r0=r0-(ac-90d0)/8d0+3.5d0
      if(ac.gt.100d0) r0=r0-(ac-100d0)/8d0+3.5d0
      if(ac.gt.110d0) r0=r0+(ac-110d0)/8d0+0.5d0

      if(ac.gt.50d0.and.m.gt.15)  r0=r0-n/40d0
             jcn=jcn+1

c        write(*,*) in,r0

c             print *, n, '   +'
c            write(7,*) n, '   +'
             iv = 11
      GO TO 2
          endif
c~~~~~~~~~~~~~~~~~

                    nt = 0
cccc                    nt = 3
                    GO TO 401

*             write(*,*) in,r0

c                   GO TO 2
  400 CONTINUE

cccccc         sdelta=s

             JN=JN+1
c 3000        FORMAT(1X,'CDLAMo_l N,RL(1)=',I5,5D20.10)
             RL(1)=R0-JN*dx*SC
             IF(JN.EQ.1001) GO TO 1112
c            GO TO 2
             GO TO 22
  401 CONTINUE

c             print *, n, '   401'
c             write(7,*) n, '   401'

             iv = 401

      RL(1)=S
      r0=rl(1)
      IF(CC.LE.0.5D0) GO TO 2
      N3=2*N-1
      N4=N3+4
      IF(CC.GT.(N+3D0)) GO TO 1
c -------------------------------------
C** Initial approximation (small 'c')
c -------------------------------------
      N2=N+M
      RL(1)=N*(N+1D0)-AK*RC2/2D0*(am2/N3/N4-1D0)+RC4/2D0*
     *((N1-1D0)*N1*(N2-1D0)*N2/((N3+2D0)*N3**3*N4)-(N1+1D0)*
     *(N1+2D0)*(N2+1D0)*(N2+2D0)/((N4+2D0)*N4**3*(N3+2D0)))
      r0=rl(1)
c                                 WRITE(7,3000) N,RL(1)
c                                 WRITE(*,3000) N,RL(1)
      GO TO 2
    1 CONTINUE
c -------------------------------------
C** Initial approximation (large 'c')
c -------------------------------------
      IF(MOD(N1,2)) 3,4,3
    3 q=N
      as = 1d0
      GO TO 16
    4 q=N+1d0
      as = - 1d0
   16 continue
      q2=q*q
      q4=q2*q2
      q6=q4*q2
      d=q2+1d0-M2
      RL(1)=-RC2+2D0*C*q-d/2D0-q*d/8D0/C
     *      -(4D0*q2*(d+1D0)+d**2)/64D0/RC2
     *      -q*(33d0*q4+114d0*q2+37d0-2d0*m2*(23d0*q2+25d0)
     *      +13*m2**2)/512d0/rc2/c
     *      -(63d0*q6+340d0*q4+239d0*q2+14d0
     *      -10d0*m2*(10d0*q4+23d0*q2+3d0)+m2*m2*(39d0*q2-18d0)
     *      -2d0*m2**3)/1024d0/rc4
     *      -q*(527d0*q6+4139d0*q4+5221d0*q2+1009d0
     *      -m2*(939d0*q4+3750d0*q2+1591d0)+m2*m2*(465d0*q2+635d0)
     *      -53d0*m2**3)/8192d0/rc4/c
c     *      +as*(4d0*c)**(q+1d0)/cdexp(2d0*c)/fact(n+1)/fact(n+m+1)
c     *      *(1d0-(3d0*q2+4d0*q+1d0-m2)/8d0/c
c     *      +((9d0*q4+4d0*q2*q-18d0*q2-12d0*q-7d0)
c     *      -m2*(6d0*q2-4d0*q-6d0)+m2*m2)/128d0/rc2
c     *      -((23d0*q6-44d0*q4*q-171d0*q4+400d0*q2*q+257d0*q2+444d0*q
c     *      +51d0)-m2*(19d0*q4-56d0*q2*q-50d0*q2+240d0*q+as*55d0)
c     *      +m2*m2*(5d0*q2-12d0*q+5d0)-m2*m2*m2)/3d0/1024d0/rc2/c)

c       print *, '2'
c       pause


       IF(in.GT.2) then

               IF(dreal(vl(in-1)).lT.0d0.and.dreal(vl(in-2)).lT.0d0.
     *             and.dreal(rl(1)).gT.0d0) then

               if(inull.eq.0.or.ac.lt.83d0.or.dreal(rl(1)).ge.ac/2d0)
     *            then
                 rl(1) = s
                 go to 121
                 end if

               if(inull.ne.0) then
               if(ac.ge.83d0.and.dreal(rl(1)).lT.ac/2d0)
     *            rl(1) = vl(in-1) / 2d0
               sdelta = s
               end if

  121          CONTINUE
               end if
       end if

c       print *, '3'
c       pause

      if(ac.gt.24.92d0.and.ac.lt.24.93d0.and.m.eq.1.and.n.eq.15) then
      rl(1)=-2.78d0
      inout = 32
      end if

      if(ac.gt.55.32d0.and.ac.lt.55.33d0.and.m.eq.1.and.n.eq.34) then
      rl(1)=-52.247d0
      inout = 32
      sdelta = s
      end if

      if(ac.gt.63.69d0.and.ac.lt.63.70d0.and.m.eq.1.and.n.eq.37) then
      rl(1)=-144.2d0
      inout = 32
      sdelta = s
      end if

      r0=rl(1)
    2 CONTINUE


*****************************************

cc      IF(JN.EQ.0.and.n1001.eq.0) then
      IF(JN.EQ.0) then


      delta=sdelta
c      if(dreal(r0).le.0d0.or.nt.ne.0.or.dreal(vl(in-1)).gt.5d0*ac)
c       if(dreal(r0).le.0d0.or.nt.ne.0)
       if(dreal(r0).lt.0d0.or.nt.ne.0)
     *                   r0 = r0 + delta
         if(m.gt.6.and.dabs(dreal(r0)).lt.cn.
     *      and.dreal(dx).lt.0.7d0) dx = 0.7d0
c         if(m.gt.6.and.n1001.eq.0)  then
         if(m.gt.6.and.n1001.eq.0.and.in.gt.1)  then
         if(dabs(dreal(r0)).lt.cn.or.dabs(dreal(vl(in-1))).lt.cn) then
c     *                   dx = 0.5d0
                        dx = 0.7d0
             if(m.gt.11.and.ac*m.gt.400d0.and.mod(n-m+1,2).ne.0) then
c             if(m.gt.11.and.ac*m.gt.600d0) then
             if(dreal(sc).gt.0d0) sc = - sc
             dx = 2d0*(ac+m-11)/10d0
             end if
             end if
         end if

c       print *, '4'
c       pause

      if(in.gt.1.and.dreal(r0).lt.dreal(vl(in-1)))  r0 = vl(in-1)
      rl(1) = r0

c         write(*,*) r0
      end if

c      alam0(in)=r0
c         write(7,*) in, r0
*****************************************
   22 CONTINUE
c      RDEL=dx
c      if(dreal(dx).gt.0.2d0) RDEL = dcmplx(0.2d0,0.2d0*af)
        RDEL0 = dcmplx(0.1d0,0.1d0*af)
        RDEL = rdel0
  222 CONTINUE
      NUM=1

  100 CONTINUE
c -------------------------------------
C*       Iteration: `from below to top'
c -------------------------------------
      IF(MOD(N1,2)) 201,200,201
  200 RDA(2)=(amm+3D0)*(amm+5D0)/((amm+2D0)*(amm+1D0)*AK*RC2)*
     *       (RL(NUM)-aM*(aM+1D0)-AK*RC2/(amm+3D0))

      IF(NM2.EQ.2) GO TO 101
      I1=4
      GO TO  1102
  201 RDA(3)=(amm+5D0)*(amm+7D0)/((amm+3D0)*(amm+2D0)*AK*RC2)*
     *       (RL(NUM)-(aM+1D0)*(aM+2D0)-(6D0*aM+3D0)*AK*RC2/((amm+1D0)*
     *       (amm+5D0)))
      IF(NM2.EQ.3) GO TO 101
      I1=5
 1102 CONTINUE
      DO 11 I=I1,NM2,2
      aJJ = I - 2D0
      aII = 2d0 * aJJ + aMM

      RDA(I)=(aII+3D0)*(aII+5D0)/((amm+aJJ+2D0)*(amm+aJJ+1D0))*
     *       ((RL(NUM)-(aM+aJJ)*(aM+aJJ+1D0)-(2D0*(aM+aJJ)*(aM+aJJ+1D0)
     *       -2D0*M2-1D0)*AK*RC2/((aII-1D0)*(aII+3D0)))/
c     *       RC2*AK-aJJ*(aJJ-1D0)/((aII-3D0)*(aII-1D0))/RDA(I-2))
     *       RC2*AK-aJJ*(aJJ-1D0)/((aII-3D0)*(aII-1D0))/
     *       (RDA(I-2)+1d-50))
   11 CONTINUE
  101 CONTINUE
      RAA=RDA(NM2)
c -------------------------------------
C*       Iteration: `from top to below'
c -------------------------------------
      RDA(NMAX)=S
      NMX1=NMAX-2
      DO 12 I=2,NMX1,2
      J = NMAX - I
      aJ = NMAX - I
      aJJ = aMM + 2d0 * aJ
      RDA(J)=aj*(aj-1D0)*(ajJ+3D0)*(ajJ+5D0)/((ajJ-3D0)*(ajJ-1D0)*
     *       (amm+aj+2D0)*(amm+aj+1D0))/((aJJ+3D0)*(aJJ+5D0)/
     *       ((amm+aJ+2D0)*(amm+aJ+1D0)*AK*RC2)*(RL(NUM)-(aM+aJ)*
     *       (aM+aJ+1D0)-(2D0*(aM+aJ)*(M+aJ+1D0)-2D0*M2-1D0)*AK*RC2/
     *       ((aJJ-1D0)*(aJJ+3D0)))-RDA(J+2))
   12 CONTINUE
      RBB=RDA(NM2)
      SL(NUM)=RAA-RBB

c      if(in.eq.1)   write(7,*) 'sl  ', num, sl(num)

      NUM=NUM+1

      if(num.gt.540) eps1 = 1d-10
*      if(num.gt.500) eps1 = 1d-08
*      if(num.gt.540)  print *, num

c -------------------------------------
C***       Iteration scheme      ***
c -------------------------------------
      IF(MOD(NUM,2)) 13,14,13
   13 CONTINUE

*         print *,num,num-2

      RL(NUM)=RL(NUM-2)
      IF(CDABS(SL(NUM-1)-SL(NUM-2)).LE.EPS1) GO TO 15
      RL(NUM)=RL(NUM-2)+SL(NUM-2)*RDEL/(SL(NUM-2)-SL(NUM-1))
      IF(CDABS(RL(NUM)-RL(NUM-2)).LE.EPS1) GO TO 15
      GO TO 100
   14 CONTINUE
      IF(NUM.EQ.2) GO TO 18
      RDEL=RDEL*SL(NUM-1)/(SL(NUM-3)-SL(NUM-2))
   18 RL(NUM)=RL(NUM-1)+RDEL
      GO TO 100
c -------------------------------------
C***    Eigenvalue is found !!! (YPA!)
c -------------------------------------
   15 VL(IN)=RL(NUM)
      IF(aC.LT.2.5D0) GO TO 110
      IF(in.le.2) GO TO 110

        Cc2=dreal(rl(1))

***        if(m.ge.10.and.dreal(vl(in-1)).gt.0d0)  Cc2=dreal(r0)

        Cc1=dreal(VL(IN))
        CFF=dabs(cc1-cc2)
c                if(cff.gt.10d0*ac.and.cc1.lt.10d0*ac) then
                if(cff.gt.20d0*ac.and.cc1.lt.-ac) then
                rdel = rdel/2d0
                rl(1) = r0
*                write(*,*) in, rdel
*                write(*,*) cff, cc1
                go to 222
                end if
        cn=ac
        an=dfloat(n)
        if(an.lt.cn) cn=an
                   IF(n1001.eq.2.
     *                and.dreal(vl(in)).gT.0d0.
     *                and.dreal(vl(in-1)).lT.0d0.
     *                and.dreal(vl(in-2)).lT.0d0) then
                cn=ac+an
                end if
         cbb=dreal(VL(IN))-dreal(VL(IN-1))

        IF(CFF.gt.cn) go to 400
        IF(cbb.le.0d0.and.dreal(vl(in)).gT.-5*ac) go to 400
        IF(cbb.lt.0d0.and.dreal(vl(in)).le.-5*ac) go to 400

 1112  CONTINUE

c 3001  FORMAT(I2,3D15.8,1x,2d10.3)
c 3002  FORMAT(1x,'L',2f7.3,2I5,2x,3(f10.3,f10.3))

        ch=r0
        cr=rl(1)
        crn=vl(in-2)
        cc0=vl(in-1)
        ck=vl(in)
c       if(jn.gt.1) WRITE(7,3002) c,N,jn,ch,cr,crn,cc0,ck
c       if(jn.gt.1) WRITE(*,3002) c,N,jn,ch,cr,crn,cc0,ck
            if(jn.lt.1001) go to 110
            sc = -sc
            n1001=n1001+1

c            if(n1001.ge.4) then
            if(n1001.ge.inout) then
c       write(7,*) ac
c       write(*,*) ac
c       write(7,*) '**   Problems with computations of eigenvalues!'
c       write(*,*) '**   Problems with computations of eigenvalues!'
            ier = 8
            go to 900
            end if

  109  CONTINUE
ccc            if(mod(n1001,2).eq.0) dx=dx/1.5d0
            if(mod(n1001,2).eq.0.and.nt.eq.0) dx=dx/3d0
            go to 1113
  110  CONTINUE
c        sdelta=vl(in)-r0
       if(dreal(vl(in)).gt.0d0) sdelta=vl(in)-r0
       if(dreal(vl(in)).lt.0d0.and.mod(in,2).ne.0) sdelta=vl(in)-r0

                if(cdabs(sdelta).gt.5d0*ac.and.in.lt.3) then
c                rdel = rdel/200d0
                rdel = rdel/2d0
                eps1 =  eps1 * 10d0
c                write(*,*) in, rdel
c                write(*,*) r0, vl(in)
                go to 222
                end if

c       print *, '5'
c       pause

      IF(in.le.2) go to 815

         IF(dreal(vl(in)-vl(in-1)).gt.in.or.
c         IF(dreal(vl(in)-vl(in-1)).gt.in.and.
c     *      dreal(vl(in-1)-vl(in-2)).gt.in) go to 111
     *      dreal(vl(in-1)-vl(in-2)).gt.in-1-ac/10d0) go to 111

  815 continue
c       print *, '6'
c       pause


          IF(CC.GT.N+3D0.and.iv.ne.401) then

      IF(MOD(N1,2)) 83,84,83
   83 q=N
      as = 1d0
      GO TO 816
   84 q=N+1d0
      as = - 1d0
  816 continue
      q2=q*q
      q4=q2*q2
      q6=q4*q2
      d=q2+1d0-M2
      RL(1)=-RC2+2D0*C*q-d/2D0-q*d/8D0/C
     *      -(4D0*q2*(d+1D0)+d**2)/64D0/RC2
     *      -q*(33d0*q4+114d0*q2+37d0-2d0*m2*(23d0*q2+25d0)
     *      +13*m2**2)/512d0/rc2/c
     *      -(63d0*q6+340d0*q4+239d0*q2+14d0
     *      -10d0*m2*(10d0*q4+23d0*q2+3d0)+m2*m2*(39d0*q2-18d0)
     *      -2d0*m2**3)/1024d0/rc4
     *      -q*(527d0*q6+4139d0*q4+5221d0*q2+1009d0
     *      -m2*(939d0*q4+3750d0*q2+1591d0)+m2*m2*(465d0*q2+635d0)
     *      -53d0*m2**3)/8192d0/rc4/c
c     *      +as*(4d0*c)**(q+1d0)/cdexp(2d0*c)/fact(n+1)/fact(n+m+1)
     *      +as*(4d0*c)**(q+1d0)/cdexp(2d0*c)/fact2(n+m,n)
     *      *(1d0-(3d0*q2+4d0*q+1d0-m2)/8d0/c
     *      +((9d0*q4+4d0*q2*q-18d0*q2-12d0*q-7d0)
     *      -m2*(6d0*q2-4d0*q-6d0)+m2*m2)/128d0/rc2
     *      -((23d0*q6-44d0*q4*q-171d0*q4+400d0*q2*q+257d0*q2+444d0*q
     *      +51d0)-m2*(19d0*q4-56d0*q2*q-50d0*q2+240d0*q+as*55d0)
     *      +m2*m2*(5d0*q2-12d0*q+5d0)-m2*m2*m2)/3d0/1024d0/rc2/c)

               IF(dreal(vl(in-1)).lT.0d0.and.dreal(vl(in-2)).lT.0d0.
     *             and.dreal(rl(1)).gT.0d0) then
                           rl(1) = s
               end if

          r0=rl(1)
          IF(dabs(dreal(vl(in)-r0)).gT.cn) then
          iv = 401
          if(dreal(sc).lt.0d0) sc = - sc
          go to 2
          end if
          end if
c             write(*,*) in,r0,vl(in)

  111  CONTINUE
      N = nn
      N = M + nN - 1
      N1 = N - M
      N3=2*N-1
      N4=N3+4
      N2=N+M
      R_cont = N*(N+1D0)-AK*RC2/2D0*(am2/N3/N4-1D0)+RC4/2D0*
     *((N1-1D0)*N1*(N2-1D0)*N2/((N3+2D0)*N3**3*N4)-(N1+1D0)*
     *(N1+2D0)*(N2+1D0)*(N2+2D0)/((N4+2D0)*N4**3*(N3+2D0)))

      if(cdabs(r_cont - vl(nn)).gt.3d0*ac) then
      inull = inull + 1
c      write(7,*) '**   ', r_cont
c      write(7,*) '**   ', vl(nn)
c      write(*,*) '**   ', r_cont
c      write(*,*) '**   ', vl(nn)
      if(inull.eq.2) then
      ier = 8
      go to 900
      end if
      go to 991
      end if

  900 CONTINUE
      return
      END

C**********************************************************************
C CDLAMo_s  - Oblate eigenvalues (small values)
C**********************************************************************
C
      SUBROUTINE CDLAMo_s(VL,M,NN,C,EPS,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
c      COMPLEX*16 VL,C,dx,alam0
      COMPLEX*16 VL,C,dx
      DIMENSION VL(NN),RL(550),SL(550),RDA(550)
c      DIMENSION VL(NN),RL(1550),SL(1550),RDA(1550)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal

*****************************************
c      COMMON /r0/ alam0(500)
*****************************************

      eps1=eps
      AC=C
      CC=CDABS(C)
      AF=0D0
      IF(DABS(AC-CC).GT.1D-10) AF=0.1D0
      SC=DCMPLX(0.1D0,AF)
      IER=0
c       print *, c, eps1, nn
c       pause
c -------------------------------------
C**  Parameters check-up
c -------------------------------------
      IF(AC.LT.0D0) IER=3
      IF(IER.NE.0) GO TO 900
      DO 112 IN=1,NN
 112      vl(in)=s

      DO 113 IN = 1,550
          rl(in) = s
          sl(in) = s
 113      rda(in) = s
c      sdelta=s

      M2=M*M
      MM=2*M
      AM=4D0*M2-1D0
      RC2=C*C
      RC4=RC2*RC2
      NC=idint(aC)
             jcn=0
             DO 111 IN=1,NN

c       print *, in
c       pause
              inout = 8

             dx=dcmplx(0.1d0,0.1d0*af)
             IF(aC.gt.20D0.and.m.gt.10)
     *             dx = dcmplx(0.1d0+(m-10d0)/100d0,0.1d0*af)
             n1001=0
 1113 continue
             JN=0
      N=M+IN-1
      N1=N-M
      NM2=N-M+2
      NMAX=NM2+2*CC+2
c      IF(NMAX.GT.500) WRITE(7,2000) NMAX
c 2000 FORMAT(14X,'******* CDLAMo_s ****** NMAX=',I4)
      IF(MOD((NMAX-NM2),2))250,251,250
  250 NMAX=NMAX+1
  251 CONTINUE

c       print *, in-1
c       pause
      IF(aC.LT.5.50D0) GO TO 401
c      IF(in.lT.3) GO TO 401
      IF(in.lT.7) GO TO 401
ccccc      IF(in.lT.7.and.m.lt.12) GO TO 401
      IF(n.lT.Nc-nc/3) GO TO 401
      IF(dreal(vl(in-1)).lT.0d0) GO TO 401
c      IF(dreal(vl(in-1)-vl(in-2)).lT.n.and.mod(n-m,2).ne.0) GO TO 401
      r0=VL(IN-1)+2d0*N+1d0
c        IF(aC.gt.22D0) r0 = r0 - 2.5d0
        IF(aC.gt.18D0) r0 = r0 - 2.532d0
             jcn=jcn+1

         IF(m.eq.1.and.aC.gt.15.998D0.and.aC.lt.16.003D0.and.in.eq.13)
     *      then
            r0 = 66.1923d0 - (ac - 16d0) * 10d0
            dx = dx / 100d0
         end if

         IF(m.eq.7.and.aC.gt.18.690D0.and.aC.lt.18.695D0.and.in.eq.10)
     *      then
            r0 = 134.161d0 - (ac - 18.691d0) * 13.69744d0
            dx = dx / 5.5d0
         end if
         IF(m.eq.8.and.aC.gt.18.535D0.and.aC.lt.18.540D0.and.in.eq.8)
     *      then
            r0 = 116.313d0 - (ac - 18.535d0) * 13.69744d0
            dx = dx / 5.5d0
         end if

      GO TO 2
  400 CONTINUE
         JN=JN+1
           RL(1)=R0-JN*dx*SC
         IF(JN.EQ.1001) GO TO 1112
         GO TO 2
  401 CONTINUE
      RL(1)=S
      IF(CC.LE.0.5D0) GO TO 2
      N3=2*N-1
      N4=N3+4
      IF(CC.GT.(N+3D0)) GO TO 1
c      IF(CC.GT.(N+3D0).and.cc.le.21d0) GO TO 1
c      IF(dabs(CC-(N+3D0)).gt.0.5d0.and.cc.gt.21d0) GO TO 1
c -------------------------------------
C** Initial approximation (small 'c')
c -------------------------------------
      N2=N+M
      RL(1)=N*(N+1D0)-AK*RC2/2D0*(AM/N3/N4-1D0)+RC4/2D0*
     *((N1-1D0)*N1*(N2-1D0)*N2/((N3+2D0)*N3**3*N4)-(N1+1D0)*
     *(N1+2D0)*(N2+1D0)*(N2+2D0)/((N4+2D0)*N4**3*(N3+2D0)))
      r0=rl(1)

c       print *, '2'
c       pause

      GO TO 2
    1 CONTINUE
c -------------------------------------
C** Initial approximation (large 'c')
c -------------------------------------
      IF(MOD(N1,2)) 3,4,3
    3 q=N
      GO TO 16
    4 q=N+1d0
   16 continue
      q2=q*q
      q4=q2*q2
      q6=q4*q2
      d=q2+1d0-M2
      RL(1)=-RC2+2D0*C*q-d/2D0-q*d/8D0/C
     *      -(4D0*q2*(d+1D0)+d**2)/64D0/RC2
     *      -q*(33d0*q4+114d0*q2+37d0-2d0*m2*(23d0*q2+25d0)
     *      +13*m2**2)/512d0/rc2/c
     *      -(63d0*q6+340d0*q4+239d0*q2+14d0
     *      -10d0*m2*(10d0*q4+23d0*q2+3d0)+m2*m2*(39d0*q2-18d0)
     *      -2d0*m2**3)/1024d0/rc4
     *      -q*(527d0*q6+4139d0*q4+5221d0*q2+1009d0
     *      -m2*(939d0*q4+3750d0*q2+1591d0)+m2*m2*(465d0*q2+635d0)
     *      -53d0*m2**3)/8192d0/rc4/c
      r0=rl(1)
    2 CONTINUE
c       print *, '3'
c       print *, in
c       pause
      IF(in.GT.1) then
        IF(dreal(vl(in-1)).lT.0d0.and.jn.eq.1) then
           r0=r0-2d0
           IF(dreal(r0).gt.0d0) r0=r0-2d0
        end if
      end if
               RDEL=dx
      NUM=1

ccc          write(*,*) n,r0

c       print *, '4'
c       pause
c -------------------------------------
c      alam0(in)=r0
c -------------------------------------

  100 CONTINUE
c -------------------------------------
C*       Iteration: 'from below to top'
c -------------------------------------
      IF(MOD(N1,2)) 201,200,201
  200 RDA(2)=(MM+3D0)*(MM+5D0)/((MM+2D0)*(MM+1D0)*AK*RC2)*
     *       (RL(NUM)-M*(M+1D0)-AK*RC2/(MM+3D0))
      IF(NM2.EQ.2) GO TO 101
      I1=4
      GO TO  1102
  201 RDA(3)=(MM+5D0)*(MM+7D0)/((MM+3D0)*(MM+2D0)*AK*RC2)*
     *       (RL(NUM)-(M+1D0)*(M+2D0)-(6D0*M+3D0)*AK*RC2/((MM+1D0)*
     *       (MM+5D0)))
      IF(NM2.EQ.3) GO TO 101
      I1=5
 1102 CONTINUE
      DO 11 I=I1,NM2,2
      JJ=I-2D0
      II=2D0*JJ+MM
      RDA(I)=(II+3D0)*(II+5D0)/((MM+JJ+2D0)*(MM+JJ+1D0))*((RL(NUM)-
     *       (M+JJ)*(M+JJ+1D0)-(2D0*(M+JJ)*(M+JJ+1D0)-2D0*M2-1D0)*AK*
     *   RC2/((II-1D0)*(II+3D0)))/RC2*AK-JJ*(JJ-1D0)/((II-3D0)*(II-1D0))
     *      /RDA(I-2))
   11 CONTINUE
  101 CONTINUE
      RAA=RDA(NM2)
c -------------------------------------
C*       Iteration: 'from top to below'
c -------------------------------------
      RDA(NMAX)=S
      NMX1=NMAX-2
      DO 12 I=2,NMX1,2
      J=NMAX-I
      JJ=MM+2*J
      RDA(J)=J*(J-1D0)*(JJ+3D0)*(JJ+5D0)/((JJ-3D0)*(JJ-1D0)*
     *       (MM+J+2D0)*(MM+J+1D0))/((JJ+3D0)*(JJ+5D0)/((MM+J+2D0)*
     *       (MM+J+1D0)*AK*RC2)*(RL(NUM)-(M+J)*
     *       (M+J+1D0)-(2D0*(M+J)*(M+J+1D0)-2D0*M2-1D0)*AK*RC2/
     *       ((JJ-1D0)*(JJ+3D0)))-RDA(J+2))
   12 CONTINUE
      RBB=RDA(NM2)
      SL(NUM)=RAA-RBB
      NUM=NUM+1

      if(num.gt.540) eps1 = 1d-10
c      if(num.gt.540) eps1 = eps1 * 10d0
c      if(num.gt.540)  print *, num, eps1
c       print *, num, eps1
c       pause

c -------------------------------------
C***       Iteration scheme      ***
c -------------------------------------
      IF(MOD(NUM,2)) 13,14,13
   13 CONTINUE
      RL(NUM)=RL(NUM-2)
      IF(CDABS(SL(NUM-1)-SL(NUM-2)).LE.EPS1) GO TO 15
      RL(NUM)=RL(NUM-2)+SL(NUM-2)*RDEL/(SL(NUM-2)-SL(NUM-1))
      IF(CDABS(RL(NUM)-RL(NUM-2)).LE.EPS1) GO TO 15
      GO TO 100
   14 CONTINUE
      IF(NUM.EQ.2) GO TO 18
      RDEL=RDEL*SL(NUM-1)/(SL(NUM-3)-SL(NUM-2))
   18 RL(NUM)=RL(NUM-1)+RDEL
      GO TO 100
c -------------------------------------
C***    Eigenvalue is found !!! (YPA!)
c -------------------------------------
   15 VL(IN)=RL(NUM)
      IF(aC.LT.2.5D0) GO TO 111
      IF(in.le.2) GO TO 111
        Cc2=dreal(rl(1))
        Cc1=dreal(VL(IN))
        CFF=dabs(cc1-cc2)
        cn=ac
        an=dfloat(n)
        if(an.lt.cn) cn=an
        if(an.gt.cn.and.ac.lt.10d0) cn=an
         cbb=dreal(VL(IN))-dreal(VL(IN-1))
c        IF(CFF.gt.cn.or.cbb.lt.0d0) go to 400
        IF(CFF.gt.cn.or.cbb.le.0d0) go to 400
 1112  CONTINUE
c 3001  FORMAT(I2,3D15.8,1x,2d10.3)
c 3002  FORMAT(1x,'s',2f7.3,2I5,2x,3(f10.3,f10.3))
        ch=r0
        cr=rl(1)
        crn=vl(in-2)
        cc0=vl(in-1)
        ck=vl(in)
c       if(jn.gt.1) WRITE(7,3002) c,N,jn,ch,cr,crn,cc0,ck
c       if(jn.gt.1) WRITE(*,3002) c,N,jn,ch,cr,crn,cc0,ck

            if(jn.lt.1001) then
            go to 111
            end if

            sc=-sc
            n1001=n1001+1

            if(n1001.ge.inout) then
c       write(7,*) '**   Problems with computations of eigenvalues!'
c       write(*,*) '**   Problems with computations of eigenvalues!'
            ier = 8
            go to 900
            end if

  878 continue

  119  CONTINUE
            if(mod(n1001,2).eq.0) dx=dx/3d0
            go to 1113

  111  CONTINUE
  900 CONTINUE
      return
      END

C**********************************************************************
C CDLAMP  - Prolate eigenvalues
C
      SUBROUTINE CDLAMp(VL,M,NN,C,EPS,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 VL,C,dx
      DIMENSION VL(NN),RL(1550),SL(1550),RDA(1550)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal

*****************************************
*      COMMON /r0/ alam0(500)
*****************************************

      eps1=eps
      AC=C
      CC=CDABS(C)
      AF=0D0
      IF(DABS(AC-CC).GT.1D-10) AF=0.1D0
      SC=DCMPLX(0.1D0,AF)
*****      if(ac.lt.32.5+m-3.or.ac.gt.47.4d0) SC=-sc
      if(ac.lt.29d0) SC=-sc
      IER=0
c -------------------------------------
C*  Parameters check-up
c -------------------------------------
      IF(AC.LT.0D0) IER=3
      IF(IER.NE.0) GO TO 900
      DO 112 IN=1,NN
 112      vl(in)=s
      nt=0
      sdelta=s
      M2=M*M
      m4=m2*m2
      MM=2*M
      AM = 4D0*M2-1D0
      RC2=C*C
      RC4=RC2*RC2
      NC=idint(aC)
             jcn=0
             DO 111 IN=1,NN

              inout = 16
c              inout = 4
****             write(*,*) in,nt

             dx=dcmplx(0.1d0,0.1d0*af)
             if(ac.gt.60d0) dx=dcmplx(0.2d0,0.2d0*af)
             n1001=0
 1113 continue
             JN=0
      N=M+IN-1
      N1=N-M
      NM2=N-M+2
      NMAX=NM2+2*CC+2
c      IF(NMAX.GT.3090) WRITE(7,2000) NMAX
c      IF(NMAX.GT.3090) WRITE(*,2000) NMAX
c 2000 FORMAT(14X,'******* CDLAMp ****** NMAX=',I4)
      IF(MOD((NMAX-NM2),2))250,251,250
  250 NMAX=NMAX+1
  251 CONTINUE
      IF(aC.LT.11.90D0) GO TO 401
      IF(in.lT.7) GO TO 401
**      IF(n.lT.Nc-nc/3-1.and.nt.eq.0) GO TO 401
      IF(n.lT.Nc-nc/3-m.and.nt.eq.0) GO TO 401
       r0=VL(IN-1)+dabs(dreal(VL(IN-1))-dreal(VL(IN-2)))
*       if(ac.ge.30d0.and.m.gt.3) r0=r0+m/2d0
       if(ac.ge.30d0.and.m.gt.3) r0=r0+m/5d0
      if(ac.ge.16d0.and.ac.lt.29d0) r0=r0+1.5d0
      if(ac.ge.30d0.and.ac.lt.40d0) r0=r0+(ac-30d0)/2d0
      if(ac.ge.40d0.and.ac.lt.50d0) r0=r0+n/10d0-(ac-40d0)/3d0
      if(ac.ge.50d0.and.ac.le.60d0) r0=r0+(ac-50d0)/3d0
      if(ac.gt.60d0.and.ac.le.70d0) r0=r0+(ac-60d0)/3d0
      if(ac.gt.70d0.and.ac.le.80d0) r0=r0+(ac-70d0)/4d0
      if(ac.gt.80d0.and.ac.le.90d0) r0=r0-(ac-80d0)/8d0+3.5d0
      if(ac.gt.90d0.and.ac.le.100d0) r0=r0-(ac-90d0)/8d0+4d0
      if(ac.gt.100d0) r0=r0-(ac-100d0)/8d0+4.5d0
      if(ac.gt.110d0) r0=r0+(ac-110d0)/8d0+1.5d0

**********      sdelta=s

             jcn=jcn+1
      GO TO 2
  400 CONTINUE

         sdelta=s

         JN=JN+1
         RL(1)=R0-JN*dx*SC*AK
         IF(JN.EQ.1001) GO TO 1112
         GO TO 2
  401 CONTINUE
      RL(1)=S
      R0=S
      IF(CC.LE.4.3D0) GO TO 2
      IF(CC.GT.(N+3D0)) GO TO 1
c -------------------------------------
C** Initial approximation (small 'c')
c -------------------------------------
      an3=2d0*N-1d0
      an4=an3+4d0
      N2=N+M
      RL(1)=N*(N+1D0)-AK*RC2/2D0*(AM/an3/an4-1D0)+RC4/2D0
     *      *((N1-1D0)*N1*(N2-1D0)*N2/((an3+2D0)*an3**3*an4)-(N1+1D0)
     *      *(N1+2D0)*(N2+1D0)*(N2+2D0)/((an4+2D0)*an4**3*(an3+2D0)))
     *      +rc4*rc2*(4d0*m2-1d0)*((N1+1D0)*(N1+2D0)*(N2+1D0)*(N2+2D0)
     *      /(an3*(2d0*N+1D0)*an4**5*(an4+2D0)*(an4+4D0))
     *      -(N1-1D0)*N1*(N2-1D0)*N2
     *      /((an3-4D0)*(an3-2D0)*an3**5*(an3+2D0)*an4))
      r0=rl(1)
      GO TO 2
    1 CONTINUE
c -------------------------------------
C** Initial approximation (large 'c')
c -------------------------------------
       q=2d0*N1+1d0
       q2=q*q
       q3=q2*q
       q4=q3*q
      RL(1)=C*q+M2-(q2+5D0)/8D0-q*(q2+11D0-32D0*M2)/C/64D0
     *      -(5D0*(q2**2+26D0*q2+21D0)-384D0*M2*(q2+1D0))/1024D0/RC2
     *      -1d0/rc2/c*(1d0/128d0**2*(33d0*q4*q+1594d0*q3
     *      +5621d0*q)-m2/128d0*(37d0*q3+167d0*q)+m4/8d0*q)
     *      -1d0/rc2/rc2*(1d0/256d0**2*(63d0*q4*q2+4940d0*q4
     *      +43327d0*q2+22470d0)-m2/512d0*(115d0*q4
     *      +1310d0*q2+735d0)+3d0*m4/8d0*(q2+1d0))
     *      -1d0/rc2/rc2/c*(1d0/1024d0**2*(527d0*q4*q3+61529d0*q4*q
     *      +1043961d0*q3+2241599d0*q)-m2/32d0/1024d0*(5739d0*q4*q
     *      +127550d0*q3+298951d0*q)+m4/512d0*(355d0*q3+1505d0*q)
     *      -m4*m2*q/16d0)
      r0=rl(1)
    2 CONTINUE
               RDEL=0.1D0*SC
      NUM=1

*****************************************
      delta=sdelta
      IF(delta.gT.1d0) delta=1d0
      IF(delta.lT.-1d0) delta=-1d0
**      if(n1001.eq.0)  r0=r0+sdelta
      if(n1001.eq.0)  r0=r0+delta
ccc         r0=r0+delta
*      alam0(in)=r0
c      write(*,*) in, r0,sdelta

        IF(cdabs(sdelta).gT.0.8d0.and.nt.eq.0) then
        sdelta=s
        nt=1
        end if
*****************************************

  100 CONTINUE
c -------------------------------------
C*            Iteration: 'from bottom to top'
c -------------------------------------
      IF(MOD(N1,2)) 201,200,201
  200 RDA(2)=(MM+3D0)*(MM+5D0)/((MM+2D0)*(MM+1D0)*AK*RC2)*
     *       (RL(NUM)-M*(M+1D0)-AK*RC2/(MM+3D0))
      IF(NM2.EQ.2) GO TO 101
      I1 = 4
      GO TO  1102
  201 RDA(3)=(MM+5D0)*(MM+7D0)/((MM+3D0)*(MM+2D0)*AK*RC2)*
     *       (RL(NUM)-(M+1D0)*(M+2D0)-(6D0*M+3D0)*AK*RC2/((MM+1D0)*
     *       (MM+5D0)))
      IF(NM2.EQ.3) GO TO 101
      I1=5
 1102 CONTINUE
      DO 11 I=I1,NM2,2
      JJ=I-2D0
      II=2D0*JJ+MM
      RDA(I)=(II+3D0)*(II+5D0)/((MM+JJ+2D0)*(MM+JJ+1D0))*((RL(NUM)-
     *       (M+JJ)*(M+JJ+1D0)-(2D0*(M+JJ)*(M+JJ+1D0)-2D0*M2-1D0)*AK*
     *   RC2/((II-1D0)*(II+3D0)))/RC2*AK-JJ*(JJ-1D0)/((II-3D0)*(II-1D0))
     *      /RDA(I-2))
   11 CONTINUE
  101 CONTINUE
      RAA=RDA(NM2)
c -------------------------------------
C*            Iteration: 'from top to bottom'
c -------------------------------------
      RDA(NMAX)=S
      NMX1=NMAX-  2
      DO 12 I=  2,NMX1,2
      J=NMAX-I
      JJ=MM+2*J
      RDA(J)=J*(J-1D0)*(JJ+3D0)*(JJ+5D0)/((JJ-3D0)*(JJ-1D0)*
     *       (MM+J+2D0)*(MM+J+1D0))/((JJ+3D0)*(JJ+5D0)/((MM+J+2D0)*
     *       (MM+J+1D0)*AK*RC2)*(RL(NUM)-(M+J)*
     *       (M+J+1D0)-(2D0*(M+J)*(M+J+1D0)-2D0*M2-1D0)*AK*RC2/
     *       ((JJ-1D0)*(JJ+3D0)))-RDA(J+2))
   12 CONTINUE
      RBB=RDA(NM2)
      SL(NUM)=RAA-RBB
      NUM=NUM+1

      if(num.gt.540) eps1 = 1d-10
*      if(num.gt.500) eps1 = 1d-08
*      if(num.gt.540)  print *, num

c -------------------------------------
C***       Iteration scheme      ***
c -------------------------------------
      IF(MOD(NUM,2)) 13,14,13
   13 CONTINUE
      IF(NUM.gE.1009) GO TO 15
      RL(NUM)=RL(NUM-2)
      IF(CDABS(SL(NUM-1)-SL(NUM-2)).LE.EPS1) GO TO 15
      RL(NUM)=RL(NUM-2)+SL(NUM-2)*RDEL/(SL(NUM-2)-SL(NUM-1))
      IF(CDABS(RL(NUM)-RL(NUM-2)).LE.EPS1) GO TO 15
      GO TO 100
   14 CONTINUE
      IF(NUM.EQ.2) GO TO 18
      RDEL=RDEL*SL(NUM-1)/(SL(NUM-3)-SL(NUM-2))
   18 RL(NUM)=RL(NUM-1)+RDEL
      GO TO 100
c -------------------------------------
C***    Eigenvalue is found !!! (YPA!)
c -------------------------------------
   15 VL(IN)=RL(NUM)
       sdelta=vl(in)-r0
      IF(aC.LT.11.90D0) GO TO 111
      IF(in.lT.3) GO TO 111
        Cc2=r0
        Cc1=VL(IN)
        CFF=dabs(cc1-cc2)
        cn=ac
        an=dfloat(n)
        if(an.lt.cn) cn=an
        IF(CFF.ge.cn) go to 400
 1112  CONTINUE
c 3002  FORMAT(1x,'p',2f6.2,2I5,2x,3(f10.3,f10.3))
        if(jn.gt.1) then
        cr=rl(1)
        crn=vl(in-2)
        cc0=vl(in-1)
c        WRITE(7,3002) c,N,jn,cc2,cr,crn,cc0,cc1
c        WRITE(*,3002) c,N,jn,cc2,cr,crn,cc0,cc1
        end if
            if(jn.lt.1001) go to 111
            sc=-sc
            n1001=n1001+1
            if(mod(n1001,2).eq.0) dx=dx/1.5d0

            if(n1001.ge.inout) then
c       write(7,*) ac
c       write(*,*) ac
c       write(7,*) '**   Problems with computations of eigenvalues!'
c       write(*,*) '**   Problems with computations of eigenvalues!'
            ier = 8
            go to 900
            end if

c            if(mod(n1001,48).eq.0) go to 900
cyy            if(mod(n1001,2).eq.0) go to 900
            go to 1113
  111  CONTINUE
  900 CONTINUE
      return
      END
C**********************************************************************
C CDLAMn  - ... eigenvalues
C
      SUBROUTINE CDLAMn(vl0,VL,M,NN,C,EPS,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 VL,C,dx,vl0(nn)
c      DIMENSION VL(NN),RL(550),SL(550),RDA(550)
      DIMENSION VL(NN),RL(1550),SL(1550),RDA(1550)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      eps1=eps
      AC=C
      CC=CDABS(C)
      AF=0D0
      IF(DABS(AC-CC).GT.1D-10) AF=0.1D0
      SC=DCMPLX(0.1D0,AF)
      if(ac.lt.32+m-3.or.ac.gt.47.4d0) SC=-sc
      IER=0
      in0 = 0
C*  Parameters check-up
      IF(AC.LT.0D0) IER=3
      IF(IER.NE.0) GO TO 900
      DO 112 IN=1,NN
 112      vl(in)=s
      M2=M*M
      MM=2*M
      RC2=C*C
             DO 111 IN=1,NN
             dx=dcmplx(0.1d0,0.1d0*af)
             n1001=0
c             if(dreal(sc).lt.0d0) sc = -sc
 1113 continue
             JN=0
      N=M+IN-1
      N1=N-M
      NM2=N-M+2
      NMAX=NM2+2*CC+2
c      IF(NMAX.GT.1550) WRITE(7,2000) NMAX
c      IF(NMAX.GT.1550) WRITE(*,2000) NMAX
c 2000 FORMAT(14X,'******* CDLAMn ****** NMAX=',I4)
      IF(MOD((NMAX-NM2),2))250,251,250
  250 NMAX=NMAX+1
  251 CONTINUE

          if(cdabs(vl0(in)).lt.1d-10) then
                 in0 = in
                 vl0(in)=VL0(IN-1)+dabs(dreal(VL0(IN-1)-VL0(IN-2)))

      IF(CC.le.(N+3D0)) then
      N2=N+M
      N3=2*N-1
      N4=N3+4
      AM=4D0*M2-1D0
c -------------------------------------
C** Initial approximation (small 'c')
c -------------------------------------
      vl0(in) = N*(N+1D0)-AK*RC2/2D0*(AM/N3/N4-1D0)+RC2**2/2D0*
     *          ((N1-1D0)*N1*(N2-1D0)*N2/((N3+2D0)*N3**3*N4)-(N1+1D0)*
     *          (N1+2D0)*(N2+1D0)*(N2+2D0)/((N4+2D0)*N4**3*(N3+2D0)))
      end if

c                    print *,in,vl0(in)

                  end if

                 r0=vl0(in)
                 rl(1)=r0
                 go to 2

  400 CONTINUE
         JN=JN+1
         RL(1)=R0-JN*dx*SC*AK

             IF(JN.EQ.1001.and.in.gt.2) GO TO 1112
             IF(JN.EQ.1001.and.in.le.2) GO TO 1131

    2 CONTINUE
               RDEL=dx
      NUM=1
  100 CONTINUE
c -------------------------------------
C*            Iteration: 'from below to top'
c -------------------------------------
      IF(MOD(N1,2)) 201,200,201
  200 RDA(2)=(MM+3D0)*(MM+5D0)/((MM+2D0)*(MM+1D0)*AK*RC2)*
     *       (RL(NUM)-M*(M+1D0)-AK*RC2/(MM+3D0))
      IF(NM2.EQ.2) GO TO 101
      I1 = 4
      GO TO  1102
  201 RDA(3)=(MM+5D0)*(MM+7D0)/((MM+3D0)*(MM+2D0)*AK*RC2)*
     *       (RL(NUM)-(M+1D0)*(M+2D0)-(6D0*M+3D0)*AK*RC2/((MM+1D0)*
     *       (MM+5D0)))
      IF(NM2.EQ.3) GO TO 101
      I1=5
 1102 CONTINUE
      DO 11 I=I1,NM2,2
      JJ=I-2D0
      II=2D0*JJ+MM
      RDA(I)=(II+3D0)*(II+5D0)/((MM+JJ+2D0)*(MM+JJ+1D0))*((RL(NUM)-
     *       (M+JJ)*(M+JJ+1D0)-(2D0*(M+JJ)*(M+JJ+1D0)-2D0*M2-1D0)*AK*
     *   RC2/((II-1D0)*(II+3D0)))/RC2*AK-JJ*(JJ-1D0)/((II-3D0)*(II-1D0))
     *      /RDA(I-2))
   11 CONTINUE
  101 CONTINUE
      RAA=RDA(NM2)
c -------------------------------------
C*            Iteration: 'from top to below'
c -------------------------------------
      RDA(NMAX)=S
      NMX1=NMAX-  2
      DO 12 I=  2,NMX1,2
      J=NMAX-I
      JJ=MM+2*J
      RDA(J)=J*(J-1D0)*(JJ+3D0)*(JJ+5D0)/((JJ-3D0)*(JJ-1D0)*
     *       (MM+J+2D0)*(MM+J+1D0))/((JJ+3D0)*(JJ+5D0)/((MM+J+2D0)*
     *       (MM+J+1D0)*AK*RC2)*(RL(NUM)-(M+J)*
     *       (M+J+1D0)-(2D0*(M+J)*(M+J+1D0)-2D0*M2-1D0)*AK*RC2/
     *       ((JJ-1D0)*(JJ+3D0)))-RDA(J+2))
   12 CONTINUE
      RBB=RDA(NM2)
      SL(NUM)=RAA-RBB
      NUM=NUM+1

      if(num.gt.1540) eps1 = eps1 * 10d0
*      if(num.gt.540) eps1 = eps1 * 10d0 ** (num - 540)
*      if(num.gt.545) eps1 = 1d-08
*      if(num.gt.540)  print *, num, eps1

            if(num.eq.1550) then
            ier = 8
            go to 900
            end if

c -------------------------------------
C***       Iteration scheme      ***
c -------------------------------------
      IF(MOD(NUM,2)) 13,14,13
   13 CONTINUE
      IF(NUM.gE.1009) GO TO 15
      RL(NUM)=RL(NUM-2)
      IF(CDABS(SL(NUM-1)-SL(NUM-2)).LE.EPS1) GO TO 15
      RL(NUM)=RL(NUM-2)+SL(NUM-2)*RDEL/(SL(NUM-2)-SL(NUM-1))
      IF(CDABS(RL(NUM)-RL(NUM-2)).LE.EPS1) GO TO 15
      GO TO 100
   14 CONTINUE
      IF(NUM.EQ.2) GO TO 18
      RDEL=RDEL*SL(NUM-1)/(SL(NUM-3)-SL(NUM-2))
   18 RL(NUM)=RL(NUM-1)+RDEL
      GO TO 100
c -------------------------------------
C***    Eigenvalue is found !!! (YPA!)
c -------------------------------------
   15 VL(IN)=RL(NUM)

cc                    print *,in,vl0(in)-vl(in)

      IF(aC.LT.2.5D0) GO TO 110
        Cc2=dreal(r0)
        Cc1=dreal(VL(IN))
        CFF=dabs(cc1-cc2)

        cn=ac
        an=dfloat(n)
        if(ac.gt.13d0) an = an + 1
        if(an.lt.cn) cn = an
        if(an.gt.cn.and.ac.lt.10d0) cn = an

cc        cn=cc
cc        an=dfloat(n)
cc        if(cc.gt.13d0) an = an + 1
cc        if(an.lt.cn) cn = an
cc        if(an.gt.cn.and.cc.lt.10d0) cn = an

      IF(in.le.2) GO TO 1121
      IF(in0.ne.0) GO TO 110

         cbb = dreal(VL(IN))-dreal(VL(IN-1))
c         cbb2 = cdabs(VL(IN) - VL(IN-2))

                  IF(k.eq.0.or.k.eq.1.and.
c     *                MOD(in,2).eq.0.and.dimag(c).gt.0.5d0) then
     *                 dimag(c).gt.0.5d0) then
        IF(cbb.lt.0d0)  go to 110
                   end if

c        IF(CFF.gt.cn.or.cbb.le.0d0) then
c        write(7,*) jn
c        write(7,*) cn, ac, an
c        write(7,*) jn, cbb, cbb2
c        write(7,*) jn, vl(in),cff,cbb
c        go to 400
c        end if

        IF(CFF.gt.cn.or.cbb.le.0d0) go to 400
c        IF(Cc-ac.gt.1d0.and.cbb2.le.1d-6) go to 400

 1112  CONTINUE
c 3001  FORMAT(I2,3D15.8,1x,2d10.3)
c 3002  FORMAT(1x,'*n*',2f6.2,2I5,2x,3(f10.3,f10.3),'*')
        ch=r0
        cr=rl(1)
        crn=vl(in-2)
        cc0=vl(in-1)
        ck=vl(in)

c        if(jn.gt.1) WRITE(7,3002) c,N,jn,cbb,cbb2,vl(in),vl(in-2)

c       if(jn.gt.1) WRITE(7,3002) c,N,jn,ch,cr,crn,cc0,ck
c       if(jn.gt.1) WRITE(*,3002) c,N,jn,ch,cr,crn,cc0,ck
            if(jn.lt.1001) go to 110
            sc=-sc
            n1001=n1001+1
            eps1 = eps1 * 10d0

cccc            if(n1001.ge.8) then
            if(n1001.ge.4) then
            ier = 8
            go to 900
            end if

c            if(n1001.ge.28) go to 110
            if(n1001.le.3) go to 119
             jn=0
             dx=dx/2d0
c            if(mod(n1001,2).eq.0) dx=dx/5d0
            go to 400
  119  CONTINUE
            if(mod(n1001,2).eq.0) dx=dx/1.5d0
            go to 1113
 1121  CONTINUE
        IF(CFF.gt.ac) go to 400
 1131  CONTINUE
        ch=r0
        cr=rl(1)
        ck=vl(in)
c       if(jn.gt.1) WRITE(7,3002) c,N,jn,ch,cr,ck
c       if(jn.gt.1) WRITE(*,3002) c,N,jn,ch,cr,ck
            if(jn.lt.1001) go to 110
            sc=-sc
            n1001=n1001+1
            eps1 = eps1 * 10d0

ccc*            if(n1001.ge.8) then
            if(n1001.ge.4) then
            ier = 8
            go to 900
            end if

c            if(n1001.ge.28) go to 110
            if(n1001.le.3) go to 119
             jn=0
             dx=dx/2d0
c            if(mod(n1001,2).eq.0) dx=dx/5d0
            go to 400
  110  CONTINUE

  111  CONTINUE
  900 CONTINUE
      return
      END
c************************************************
C fact2
C
      real*8 FUNCTION fact2(n,k)
      REAL*8 f
      f = 1d0
       do i = k + 1, n
        f = f * i
       end do
      fact2 = f
      RETURN
      END

c************************************************
C fact3
C
      real*8 FUNCTION fact3(n,k,f1)
      REAL*8 f, f1
      f = f1
       do i = k + 1, n
        f = f * i
       end do
      fact3 = f
      RETURN
      END

c******************************************
      SUBROUTINE DUGOLM1(X,M,LMAX,PI,TAU)
      REAL*8 X,PI(LMAX),TAU(LMAX),A
      M2=2*M
      A=1.D0
      DO I=3,M2-1
       IF(MOD(I,2).EQ.1) A=A*I
      ENDDO
      PI(1)=A
      PI(2)=A*X*(M2+1.D0)
      TAU(1)=A*M*X
      TAU(2)=A*(M2+1.D0)*((M+1.D0)*X*X-1.D0)
      DO L1=3,LMAX
        PI(L1)=((2.D0*L1+M2-3.D0)*X*PI(L1-1)-(L1+M2-2.D0)*PI(L1-2))
     *          /(L1-1.D0)
      TAU(L1)=(L1+M-1.D0)*X*PI(L1)-(L1+M2-1.D0)*PI(L1-1)
      ENDDO
      RETURN
      END

c************************************************
C funlegnn
C
      SUBROUTINE funlegnn(M,X,NF)
      parameter (nterms=330)
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 s, s1
      DIMENSION P(2*nterms),PD(2*nterms),Q(2*nterms),QD(2*nterms)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /LEG/ P, PD, Q, QD
c      IF(NF.GT.2*nterms) WRITE(7,1) NF
c    1 FORMAT(1X,'funlegnn  NF > 2*nterms, NF=',2I10)
      IF(NF.GT.2*nterms) GO TO 900
         DO I  = 1, 2*nterms
         P (I) = 0D0
         PD(I) = 0D0
         Q (I) = 0D0
         QD(I) = 0D0
         end do
      if(k.eq.0) then
        CALL DUGOLM1(X,M,Nf,P,PD)
      IF(X.LT.1.D0) THEN
        A1=DSQRT(1.D0-X*X)
      ELSE
        A1=DSQRT(X*X-1.D0)
      ENDIF
      DO I = 1, Nf
        P(I)=P(I)*A1**M
        PD(I)=-PD(I)*A1**(M-2)
        IF(X.GT.1.D0) PD(I)=-PD(I)
      ENDDO
      end if
       ier1 = 0
      if(k.eq.1) CALL DLEGF1(P(1),D,X,D,d,M,M,K,0,IER1)
      CALL DLEGF2(Q(1),D,X,D,-M,M,K,0,IER2)
c       IF((IER1+IER2).NE.0) WRITE(7,2) M,IER1,IER2
c       IF(IER2.NE.0) WRITE(7,2) M,IER2
c    2 FORMAT(1X,'funlegnn  M,I,IER1,2,3,4=',6I10)
       ier3 = 0
       DO 3 I=2,NF
      if(k.eq.1) then
       if(i.gt.2) d=P(I-2)
      CALL DLEGF1(P(I),PD(I-1),X,P(I-1),d,M+I-1,M,K,1,IER3)
      end if
      CALL DLEGF2(Q(I),QD(I-1),X,Q(I-1),-M+I-1,M,K,1,IER4)
c       IF((IER3+IER4).NE.0) WRITE(7,2) M,I,IER1,IER2,IER3,IER4
c       IF(IER4.NE.0) WRITE(7,2) M,I,IER4
    3 CONTINUE
  900 RETURN
      END
c************************************************
C DLEGf00
C
      subroutine DLEGF00(res,red,n,M,id,iER)
      IMPLICIT REAL*8(A-H,O-Z)                                        
      COMMON  /FACT/ FACT(300)
      ier=0
      res=0d0
      red=0d0
      if(n.eq.0.and.m.eq.0) res=1d0
      if(n.eq.0.and.m.eq.0)  go to 100
      if(n.le.0.or.m.lt.0) ier=1
      if(ier.ne.0) go to 100
c**************************************   << x1 = 0 >>
      IF(MOD(iabs(N-M),2)) 50,1,50
    1 continue
      k1=(m+n)/2
      e=1d0
      do 2 kk=1,k1
    2 e=e*(k1+kk)
      k2=(n-m)/2
      if(k2+1.gt.170) stop 1703
      res=(-1)**k2*e/fact(k2+1)/2d0**n
      if(id) 50,100,50
   50 continue
c=========================             <<     Ø‡Æ®ß¢Æ§≠†Ô    >>
      if(id.eq.0) go to 100
c************************************  << x1 = 0 >>
      IF(MOD(iabs(N-M),2)) 71,100,71
   71 continue
      k1=(m+n+1)/2
      e=1d0
      do 72 kk=1,k1
   72 e=e*(k1+kk)
      k2=(n-m-1)/2
      red=(-1)**k2*e/fact(k2+1)/2d0**n   
  100 CONTINUE                                                      
      RETURN                                                        
      END      

c************************************************
C DLEGf1
C   
      subroutine DLEGF1(rs,rd,X1,x2,x3,n,M,ic,id,iER)
      IMPLICIT REAL*8(A-H,O-Z) 
      complex*16 s, s1
      COMMON /FACT/ FACT(300)
      COMMON /K1/ S, S1, AKSI, AK, Kkk, NK, nal
      ier=0
      res=0d0
      red=0d0
      eps4 = 1d-20
      xx2=dabs(x2)
      nm=n-m
      j=mod(n,4)
      jm=mod(m,4)
      if(n.eq.0.and.m.eq.0) res=1d0
      if(n.eq.0.and.m.eq.0)  go to 100
      if(n.le.0.or.m.lt.0) ier=1
      if(id.ne.0.and.dabs(x1-1d0).lt.1d-6.and.ic.eq.0) ier=2
      if(ier.ne.0) go to 100
      if(dabs(x1).gt.1d-6)  go to 20
c**************************************   << x1 = 0 >>
      IF(MOD(iabs(nm),2)) 100,1,100
    1 continue
      k1=(m+n)/2
      e=1d0
      do 2 k=1,k1
    2 e=e*(k1+k)
      if(nm/2+1.gt.170) stop 1703
      res=e/fact(nm/2+1)/2d0**n
      if(id) 50,102,50
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   << x1 /= 0 >>
   20 continue
      if(m.gt.n) go to 100
      if(nm.ge.2.and.x1.gt.1d0.and.ic.eq.0) go to 48
      z=1d0/(2d0**n)
      k1=(nm+2)/2
      d=0d0
      if(ic.ne.0) go to 21
c.....................................   << ic = 0 >>
c      if(nm.gt.169) go to 291
      if(nm.gt.169) stop 170

c       print *, m, n, k1
c       pause

      do 22 k=1,k1
      k3=nm-2*(k-1)
      if(k3+1.gt.170) stop 1704
      h=1d0/fact(k3+1)
      k4=n-k+1

c       print *,  k, k4
c       pause

      do 27 i=1,k4
      if(dabs(h).gt.1d300) print *,  n, h

   27 h=h*(k4+i)
      id = 1
       if(mod(k-1,2).ne.0) id = -1
       dd = id * h / fact(k) * x1**k3
c       print *,  dd, d
      d = d + dd
      if(dabs(dd/d).lt.eps4) go to 23
   22 continue
   23 continue
       res=d*z*dsqrt(dabs(1d0-x1**2)**m)
      if(id) 50,102,50     
c<<<<<<<<<++++++++++++++++++++++++++++   << ic /= 0 >>
   21 continue
c      if(nm.gt.169) go to 29
      if(nm.gt.169) stop 1705
      do 28 k=1,k1
      k3=nm-2*(k-1)
      h=1d0/fact(k3+1)
      k4=n-k+1
c      if(k4.gt.169) stop 17055
      do 33 i=1,k4
   33 h=h*(k4+i)
   28 d=d+h/fact(k)*x1**k3
      res=d*z*dsqrt((1+x1**2)**m)
      if(id) 50,102,50     
   48 continue
      res=((2d0*n-1d0)*x1*x2-(n+m-1d0)*x3)/nm
      if(id) 50,102,50     
   50 continue
c=========================             <<  derivature  >>
c*****************************************   << x1 /= + - 1 >>
      if(ic.ne.0) go to 101
      red=((nm)*res-x1*n*x2)/(x1**2-1d0)
       rs=res
       rd=red
       return
  101 continue
      red=((nm)*res-x1*n*xx2)/(x1**2-ak)
  100 CONTINUE        
      if(jm.eq.2.or.jm.eq.3) go to 103
      rs=res
      if(j.eq.2.or.j.eq.3) rs=-dabs(res)
      go to 104
  103 continue
      rs=-dabs(res)
      if(j.eq.0.or.j.eq.1) rs=-rs
  104 continue
      j=mod(n-1,4)                          
      rd=red
      if(j.eq.2.or.j.eq.3) rd=-dabs(red)
      RETURN      
  102 continue
       rs=res
      if(ic.eq.0.and.jm.eq.3.and.x1.gt.1d0) rs=dabs(res)
      if(ic.ne.0.and.jm.eq.2) rs=-dabs(res)
      if(ic.ne.0.and.jm.eq.3) rs=-dabs(res)
      RETURN                                               
      END

c************************************************
C DLEGf2
C   
      subroutine DLEGF2(res,red,X1,x2,n,M,ic,id,iER) 
      IMPLICIT REAL*8(A-H,O-Z)  
      complex*16 s1, s
      dimension r(1500),q(1500)
      COMMON /K1/ S, S1, AKSI, AK, Kkk, NK, nal
      COMMON /FACT/ FACT(300)
      COMMON /PI/ PI
      ier=0
      res=0d0
      red=0d0
      if(n+m.lt.0) ier=1
      if(ic.eq.0.and.dabs(x1).le.1d0) ier=2
      if(ic.ne.0.and.dabs(x1).lt.0d0) ier=2
      if(n+m-1.lt.0.and.id.ne.0) ier=3
      if(ier.ne.0) go to 100
      if(dabs(x1).gt.1d-6.or.ic.eq.0)  go to 20
c**************************************   << x1 = 0 & ic /= 0 >>
      if(m.gt.n) go to 3
c                             m  < = n
      k1=n-m
      k2=n+m-1
      d=1d0
      e=1d0
      IF(MOD(k1,2)) 2,1,2
    1 continue
      if(k1.eq.0) go to 6
      do 4 i=1,k2,2
    4 d=d*i
      do 5 i=2,k1,2
    5 e=e*i
    6 res=-pi/2d0*d/e
      if(id) 50,100,50
    2 continue
      do 7 i=2,k2,2
    7 d=d*i
      do 8 i=1,k1,2
    8 e=e*i
      res=d/e
      if(id) 50,100,50
c                                m  >  n 
    3 continue
      IF(MOD(iabs(N-M),2)) 10,9,10
    9 go to 100
   10 continue
      k1=(n-m-1)/2
      k2=(m+n-1)/2
      d=1d0
      do 11 i=1,k1
   11 d=d*(k1+i)
      if(k1+1.gt.170) stop 1706
      res=2d0**n*d*fact(k1+1)
      if(id) 50,100,50
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   << x1 /= 0 >>
   20 continue
      if(m.eq.0)  go to 61
      if(m.eq.1)  go to 71
c *********************************     m  > 1
      if(ic.eq.0) res=dlog((x1+1d0)/(x1-1d0))/2d0
      if(ic.ne.0.and.x1.lt.0.5d0) res=(datan(x1)-pi/2d0)
      if(ic.ne.0.and.x1.ge.0.5d0) res=-datan(1d0/x1)
      r(1)=-ak/(dsqrt(x1**2-ak))
      r(2)=-2d0*x1*ak/(dsqrt(x1**2-ak))*r(1)
      if(m.eq.2) go to 221
      do 21 i=3,m
   21 r(i)=-2d0*(i-1d0)*x1*ak/(dsqrt(x1**2-ak))*r(i-1)+
     *     (2d0-i)*(i-1d0)*r(i-2)
  221 continue
      q0=r(m)
      if(n.eq.0) res=q0
      if(n.eq.0) go to 26
      j=n
      if(n.lt.0) j=-n
      nmax=30*x1+j+50
      if(nmax.gt.500) nmax=500
      r(nmax)=x1-dsqrt(x1**2-ak)
      xer=dabs(x1-1d0/dsqrt(3d0))
      do 22 i=nmax-1,1,-1
   22 r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
      q(1)=r(1)*q0*ak
      if(xer.lt.1d-10.and.m.eq.3) q(1)=3d0
      if(xer.lt.1d-10.and.m.eq.6) q(1)=-144d0
      if(xer.lt.1d-10.and.m.eq.9) q(1) = 4.536d4
      if(xer.lt.1d-10.and.m.eq.12) q(1) = -4.35456d+07
      if(xer.lt.1d-10.and.m.eq.15) q(1) = 9.340529643126767d+10
      if(xer.lt.1d-10.and.m.eq.18) q(1) = -3.76610155207262d+14
      if(xer.lt.1d-10.and.m.eq.21) q(1) = 2.554546682742534d+18
      if(xer.lt.1d-10.and.mod(m,3).eq.0) q0=0d0
      if(n.eq.1) go to 24
      if(n.lt.0) go to 225
      do 23 i=1,j-1
   23 q(i+1)=r(i+1)*q(i)*ak**(i+1)
       go to 24
c******************************************      n < 0
  225 continue 
      if(n.eq.-1) res=(ak*x1*q0+(m-1d0)*q(1))/m
      if(n.eq.-1) go to 26
      r(1)=q0
      r(2)=(ak*x1*q0+(m-1d0)*q(1))/m
      do 25 i=-2,n,-1
      r(1-i)=(ak**i*(2d0*i+3d0)*x1*r(-i)-
     *       (i-m+2d0)*r(-i-1))/(i+m+1d0)
   25 continue 
      res=r(-n+1)
      go to 26
   24 continue
      res=q(n)
   26 continue
      if(id) 50,100,50
c***********************************    m = 0
   61 continue
      if(ic.eq.0) res=dlog((x1+1d0)/(x1-1d0))/2d0
      if(ic.ne.0.and.x1.lt.0.5d0) res=(datan(x1)-pi/2d0)
      if(ic.ne.0.and.x1.ge.0.5d0) res=-datan(1d0/x1)
      if(n.eq.0) go to 65
      q0=res
      nmax=1.5*x1+n+10
      r(nmax)=x1-dsqrt(x1**2-ak)
      do 62 i=nmax-1,1,-1
      r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
   62 continue
      q(1)=r(1)*q0*ak
      if(n.eq.1) GO TO 64
      do 63 i=1,n-1
      q(i+1)=r(i+1)*q(i)*ak**(i+1)
   63 continue
   64 continue
      res=q(n)
   65 if(id) 50,100,50
c***********************************    m = 1
   71 continue
      res=-ak/(dsqrt(x1**2-ak))
      if(n.eq.0) go to 75
      q0=res
      if(n.eq.-1) res=ak*x1*q0
      if(n.eq.-1) go to 75
      nmax=30*x1+n+50
      r(nmax)=x1-dsqrt(x1**2-ak)
      do 72 i=nmax-1,1,-1
      r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
   72 continue
      q(1)=r(1)*q0*ak
      if(n.eq.1) GO TO 74
      do 73 i=1,n-1
   73 q(i+1)=r(i+1)*q(i)*ak**(i+1)
   74 continue
      res=q(n)
   75 if(id) 50,100,50
   50 continue
c=========================             <<     Ø‡Æ®ß¢Æ§≠†Ô    >>
      if(id.eq.0) go to 100
      if(ic.eq.0.or.mod(n,2).eq.0)  go to 174
      red=-((n-m)*res+x1*n*x2)/(x1**2-ak)
      go to 100
  174 continue
      red=((n-m)*res-x1*n*x2)/(x1**2-ak)
  100 CONTINUE                                                      
      RETURN                                                        
      END

c************************************************
C dreal ---- Attention ---> for Lahey (only)
C
      FUNCTION dreal(r)
      REAL*8 dreal
      complex*16 r
      dreal = real(r)
      RETURN
      END

