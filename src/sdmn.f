C       Shanjie Zhang and Jianming Jin
C
C       Copyrighted but permission granted to use code in programs.
C       Buy their book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
        SUBROUTINE SDMN(M,N,C,CV,KD,DF)
C
C       =====================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions, dk
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                cv --- Characteristic value
C                KD --- Function code
C                       KD=1 for prolate; KD=-1 for oblate
C       Output:  DF(k) --- Expansion coefficients dk;
C                          DF(1), DF(2), ... correspond to
C                          d0, d2, ... for even n-m and d1,
C                          d3, ... for odd n-m
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(200),D(200),G(200),DF(200)
cf2py intent(in) M, N, C, CV, KD
cf2py intent(out) DF
        NM=25+INT(0.5*(N-M)+C)
        IF (C.LT.1.0D-10) THEN
           DO 5 I=1,NM
5             DF(I)=0D0
           DF((N-M)/2+1)=1.0D0
           RETURN
        ENDIF
        CS=C*C*KD
        IP=1
        K=0
        IF (N-M.EQ.2*INT((N-M)/2)) IP=0
        DO 10 I=1,NM+2
           IF (IP.EQ.0) K=2*(I-1)
           IF (IP.EQ.1) K=2*I-1
           DK0=M+K
           DK1=M+K+1
           DK2=2*(M+K)
           D2K=2*M+K
           A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
           D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &          *(DK2+3.0))*CS
           G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
10      CONTINUE
        FS=1.0D0
        F1=0.0D0
        F0=1.0D-100
        KB=0
        DF(NM+1)=0.0D0
        FL=0.0D0
        DO 30 K=NM,1,-1
           F=-((D(K+1)-CV)*F0+A(K+1)*F1)/G(K+1)
           IF (DABS(F).GT.DABS(DF(K+1))) THEN
              DF(K)=F
              F1=F0
              F0=F
              IF (DABS(F).GT.1.0D+100) THEN
                 DO 12 K1=K,NM
12                  DF(K1)=DF(K1)*1.0D-100
                 F1=F1*1.0D-100
                 F0=F0*1.0D-100
              ENDIF
           ELSE
              KB=K
              FL=DF(K+1)
              F1=1.0D-100
              F2=-(D(1)-CV)/A(1)*F1
              DF(1)=F1
              IF (KB.EQ.1) THEN
                 FS=F2
              ELSE IF (KB.EQ.2) THEN
                 DF(2)=F2
                 FS=-((D(2)-CV)*F2+G(2)*F1)/A(2)
              ELSE
                 DF(2)=F2
                 DO 20 J=3,KB+1
                    F=-((D(J-1)-CV)*F2+G(J-1)*F1)/A(J-1)
                    IF (J.LE.KB) DF(J)=F
                    IF (DABS(F).GT.1.0D+100) THEN
                       DO 15 K1=1,J
15                        DF(K1)=DF(K1)*1.0D-100
                       F=F*1.0D-100
                       F2=F2*1.0D-100
                    ENDIF
                    F1=F2
20                  F2=F
                 FS=F
              ENDIF
              GO TO 35
           ENDIF
30      CONTINUE
35      SU1=0.0D0
        R1=1.0D0
        DO 40 J=M+IP+1,2*(M+IP)
40         R1=R1*J
        SU1=DF(1)*R1
        DO 45 K=2,KB
           R1=-R1*(K+M+IP-1.5D0)/(K-1.0D0)
45           SU1=SU1+R1*DF(K)
        SU2=0.0D0
        SW=0.0D0
        DO 50 K=KB+1,NM
           IF (K.NE.1) R1=-R1*(K+M+IP-1.5D0)/(K-1.0D0)
           SU2=SU2+R1*DF(K)
           IF (DABS(SW-SU2).LT.DABS(SU2)*1.0D-14) GOTO 55
50         SW=SU2
55      R3=1.0D0
        DO 60 J=1,(M+N+IP)/2
60         R3=R3*(J+0.5D0*(N+M+IP))
        R4=1.0D0
        DO 65 J=1,(N-M-IP)/2
65         R4=-4.0D0*R4*J
        S0=R3/(FL*(SU1/FS)+SU2)/R4
        DO 70 K=1,KB
70         DF(K)=FL/FS*S0*DF(K)
        DO 75 K=KB+1,NM
75         DF(K)=S0*DF(K)
        RETURN
        END