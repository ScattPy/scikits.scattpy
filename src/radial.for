                                                                     
                                                                     
                                                                     
                                             
c****************************************************************************
c                    Radial spheroidal functions
c
c parameters:
c    kob  - switch between prolate (0) and oblate (1) functions
c    m    - index m
c    ne   - maximum index n
c    c2   - argument c (complex*16)
c    ksi0 - argument xi (real*8)
c    eps  - accuracy of getting eigenvalues
c
c results (all complex*16):
c    R1f(ne) - functions R^(1)_mn (c,xi),  n = 1, ne
c    R1d(ne) - derivatives dR^(1)_mn (c,xi)/dxi, n = 1, ne
c    R2f(ne) - functions R^(2)_mn (c,xi),  n = 1, ne
c    R2d(ne) - derivatives dR^(2)_mn (c,xi)/dxi, n = 1, ne
c
c 2005, Nov
c****************************************************************************

      SUBROUTINE rad_fun (kob, m, ne, C2, KSI0, EPS, R1f, R1d, R2f,R2d)
      parameter (nterms=330)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z), COMPLEX*16 (R-S)
      REAL*8 C2, ksi0
      complex*16 R1f , R1d , R2f , R2d 
      COMPLEX*16 bdc2
      DIMENSION RLC2(nterms), rDC2(4*nterms), bDC2(4*nterms)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /EPS1/ EPS1
      COMMON /EPS33/ EPS33
      COMMON /EPS3/ EPS3
      COMMON /FACT/ FACT(300)
      COMMON /PI/ PI
cf2py intent(in) kob, m, ne, C2, KSI0, EPS
cf2py intent(out) R1f, R1d, R2f, R2d
c      open(unit=07,file='sph_fun.out',status='unknown',access='append')

      NMAX = ne
      NMXE = nmax
      NK = NMAX+40
      IF (NK.LT.60) NK = 60
      IF (Ksi0.gt.1.5d0) NK = nk+40
c      NK=4*nterms
      IF (NK.GT.4*nterms) NK = 4*nterms

      k = kob
      IF(K) 5005,6,5005
 5005 CONTINUE
       AK=-1D0
       GO TO 7
    6 CONTINUE
       AK= 1D0
    7 CONTINUE
      x = ksi0
      IF(K.eq.0) AKSI=X**2-1D0
      IF(K.eq.1) AKSI=X**2+1D0

c      if (k.eq.0 .and. x.lt.1.00001d0) then
c      write (*,2114) x
c      write (7,2114) x
c      end if

c      if (k.eq.1 .and. x.lt.0.00001d0) then
c      write (*,2115) x
c      write (7,2115) x
c      end if
c 2114 FORMAT(1X,3('!'), '    Prolate x > 1  x = ',f5.3)
c 2115 FORMAT(1X,3('!'), '    Oblate x > 0   x = ',f5.3)

      S=(0D0,0D0)
      S1=(0D0,1D0)
      nal = 1
      EPS=1D-15
      EPS1=1D-12
*       eps3 = 1d-80
*       eps3 = 1d-100
       eps3 = 1d-200
       eps33 = eps3
      pi = 4d0 * datan(1d0)

c factorial
       FACT(1)=1D0
       FACT(2)=1D0
       DO 130 I=3,170
  130  FACT(I)=FACT(I-1)*(I-1D0)
       FACT(170)=FACT(170)*(1.D-300)
       DO J = 170, 297
        FACT(J+1)=J*FACT(J)
       end do

c-----------------

c  112 FORMAT(1X,61('.'))
c  202 FORMAT(1X,'homFUNq NN>nterms',5X,'NN=',I5)
c  210 FORMAT(1X,'L=',I4,5X,'IER1,2,3,4,5=',5I5)
c  212 FORMAT(1X,'I,L=',2I5,5X,'IER=',I5)

      W1=1D0/(C2*(ksi0**2-1D0+2*K))
      ncc = real(c2)
      IF (NE-nterms)  40,40,41
c   41   WRITE(7,202) NE
c        WRITE(*,202) NE
   41 RETURN
   40 CONTINUE
      RC1=c2

c calc of lambda

      call lambda(K,M,ne,C2,EPS,rlc2,ie)
      if(ie.ne.0) return

      inum1 = nk
      if(nk.lt.inum1) inum1=nk+10

      if(k.eq.0) then
        ifun1 = 44
        IF(ksi0.GE.1.5D0) IFUN1 = 22
      end if

      if(k.eq.1) then
        ifun1 = 33
        IF(ksi0.GE.1.5D0) IFUN1 = 22
      end if

c calc of Leg. functions

      CALL funlegnn (m, ksi0, inum1)

c loop over l

      L = NE
c      DO 1 L=1,NE
        LM=L+M-1

        RL=RLC2(L)

        CALL cdcof4a(RL,RDC2,NK,BDC2,NK,2,rvn,1,M,LM,RC1,IER1)
       if(abs(bdc2(m+1)).le.eps33*1d20) eps33 = eps33 / 1d20

      if(k.eq.1.and.l.gt.ncc.and.ksi0.gt.0.25d0.
     *    and.ksi0.lt.0.8d0) ifun1 = 11
      if(ifun1.eq.11) go to 11
      if(ifun1.eq.22) go to 22
      if(ifun1.eq.33) go to 33
      if(ifun1.eq.44) go to 44

c LEG/LEG    ifun1 = 11
   11 CONTINUE
      CALL CDRF12cc(R1,R2,r3,r4,1,M,LM,RC1,RDC2,NK,
     *            BDC2,NK,IER2)
      W = r1 * r4 - r2 * r3
      W2a = dabs(W / W1 - 1d0)
c      if(W2a.gt.eps1)  write(*,99999) ifun1,l,W,W1,W2a
      if(W2a.gt.W0.and.k.eq.0) go to 13
      go to 5
c99999 format(1x,'ifun=',i3,2x,'L=',i3,2x,'W=',1pd13.6,2x,
c     *       'W1=',d13.6,2x,'W2=',d13.6)

c LEG/DIFF  (JAF/DIFF)  ifun1 = 13
   13  CONTINUE
      IF(K.EQ.0) CALL cDRSF20(rr3,rr4,ksi0,R1,R2,m,rc1,rl,IER4)
      W = r1 * rr4 - r2 * rr3
      W2b = dabs(W / W1 - 1d0)
      ifu = 13
      if(ifun1.eq.22) ifu = 23
      if(ifun1.eq.44) ifu = 43
c      if(W2b.gt.eps1)  write(*,99999) ifu,l,W,W1,W2b
      IF(W2b.lt.w0) then
        r3 = rr3
        r4 = rr4
      end if
      go to 5

c BES/BES    ifun1 = 22
   22  CONTINUE
      CALL CDRB12cc (R1,R2,r3,r4,1,M,LM,RC1,ksi0,RDC2,NK,IER3)
      W = r1 * r4 - r2 * r3
      W2a = dabs(W / W1 - 1d0)
      if (W2a.le.eps1) go to 5
      if (k.eq.0) go to 13
      go to 5

c DIFF/DIFF   ifun1 = 33
   33 CONTINUE
      CALL cDRSF212(R1,R2,r3,r4,2,m,LM,rc1,ksi0,rl,
     *             RDC2,NK,BDC2,m+1,IER4)
      W = r1 * r4 - r2 * r3
      W2a = dabs(W / W1 - 1d0)
      if (W2a.gt.eps1) then
c        write(*,99999) ifun1,l,W,W1,W2a
        go to 11
      end if
      go to 5

c JAF/JAF    ifun1 = 44
 44   CONTINUE
      CALL CDRG1cn(r1,r2,r3,r4,1,M,lm,rC1,ksi0,
     *              RDC2,NK,bdc2,nk,rl,IER0)
      W = r1 * r4 - r2 * r3
      W0 = dabs(W / W1 - 1d0)
      ifun1 = 44
c      if (W0.gt.eps1)  write(*,99999) ifun1,l,W,W1,W0
      if(W0.gt.eps1)  go to 13

c final results
c    5 continue
    5 R1f = r1
      R1d = r2
      R2f = r3
      R2d = r4

c       write(7,*) 'n ', l
c       write(7,*) 'R^(1) ',r1
c       write(7,*) 'R^(1)" ',r2
c       write(7,*) 'R^(2) ',r3
c       write(7,*) 'R^(2)" ', r4
c      W = r1 * r4 - r2 * r3
c       W2a = dabs(W / W1 - 1d0)
c       write(7,*) 'W, log W ',w2a, dlog10(W2a+1d-100)
c       write(*,*) 'W, log W ',w2a, dlog10(W2a+1d-100)
c       write(7,112)

    1 CONTINUE
      RETURN
      END


