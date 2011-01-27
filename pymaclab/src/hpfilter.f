C ----------------------------------------------------------------------
C  SR: hpfilt
C  Kalman smoothing routine for HP filter written by E Prescott.
C   y=data series, d=deviations from trend, t=trend, n=no. obs,
C   s=smoothing parameter (eg, 1600 for std HP).
C   Array v is scratch area and must have dimension at least 3n.
C   If IOPT=1 and n and s are the same as for the previous call,
C   the numbers in v are not recomputed.  This reduces execution
C   time by about 30 percent.  Note that if this option is exercised,
C   v cannot be used for other purposes between calls.
C   This version does NOT release the trend in order to save memory.
C ----------------------------------------------------------------------
      SUBROUTINE HPFILT(Y,D,V,N,S,IOPT)
      INTEGER*2 IOPT,NN,I,I1,IB,N
      REAL*8 Y(N),T(N),V(N,3),D(N),SS
      REAL*8 M1,M2,V11,V12,V22,X,Z,B11,B12,B22,DET,E1,E2,S
      DATA SS,NN/0.D0,0/
Cf2py intent(in) y,n,s,v,iopt
Cf2py intent(out) d
Cf2py depend(n) y,d,v
C
C     compute sequences of covariance matrix for f[x(t),x(t-1) | y(<t)]
C
      IF(IOPT.NE.1.OR.NN.NE.N.OR.S.NE.SS)  THEN
        SS=S
        NN=N
        V11=1.D0
        V22=1.D0
        V12=0.D0
       DO 5 I=3,N
        X=V11
        Z=V12
        V11=1.D0/S + 4.D0*(X-Z) + V22
        V12=2.D0*X - Z
        V22=X
        DET=V11*V22-V12*V12
        V(I,1)=V22/DET
        V(I,3)=V11/DET
        V(I,2)=-V12/DET
        X=V11+1.D0
        Z=V11
        V11=V11-V11*V11/X
        V22=V22-V12*V12/X
        V12=V12-Z*V12/X
  5    CONTINUE
                                       ENDIF
C
C     this is the forward pass
C
      M1=Y(2)
      M2=Y(1)
      DO 10 I=3,N
        X=M1
        M1=2.0*M1-M2
        M2=X
        T(I-1)= V(I,1)*M1+V(I,2)*M2
        D(I-1)= V(I,2)*M1+V(I,3)*M2
        DET=V(I,1)*V(I,3)-V(I,2)*V(I,2)
        V11=V(I,3)/DET
        V12=-V(I,2)/DET
        Z=(Y(I)-M1)/(V11+1.D0)
        M1=M1+V11*Z
        M2=M2+V12*Z
 10     CONTINUE
      T(N)=M1
      T(N-1)=M2
C
C       this is the backward pass
C
         M1=Y(N-1)
         M2=Y(N)
        DO 15 I=N-2,1,-1
           I1=I+1
           IB=N-I+1
         X=M1
         M1=2.D0*M1 - M2
         M2=X
C
C           combine info for y(.lt.i) with info for y(.ge.i)
C
         IF(I.GT.2)                 THEN
           E1=V(IB,3)*M2 + V(IB,2)*M1 + T(I)
           E2=V(IB,2)*M2 + V(IB,1)*M1 + D(I)
           B11=V(IB,3)+V(I1,1)
           B12=V(IB,2)+V(I1,2)
           B22=V(IB,1)+V(I1,3)
           DET=B11*B22-B12*B12
           T(I)=(-B12*E1+B11*E2)/DET
                                    ENDIF
C
C           end of combining
C
        DET=V(IB,1)*V(IB,3)-V(IB,2)*V(IB,2)
        V11=V(IB,3)/DET
        V12=-V(IB,2)/DET
        Z=(Y(I)-M1)/(V11+1.D0)
         M1=M1+V11*Z
         M2=M2+V12*Z
 15     CONTINUE
       T(1)=M1
       T(2)=M2
        DO 20 I=1,N
 20      D(I)=Y(I)-T(I)
        RETURN
        END
C***********************************************************************
