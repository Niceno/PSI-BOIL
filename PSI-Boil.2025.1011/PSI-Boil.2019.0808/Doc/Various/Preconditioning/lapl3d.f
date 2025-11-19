C#######################################################################
      SUBROUTINE ICCG(FI,NS)
C#######################################################################
      PARAMETER (NX=41,NY=41,NZ=41,NXYZ=NX*NY*NZ)
      COMMON /INDI/ NI,NJ,NK,NIM,NJM,NKM,NIJ,NIJK,LI(NX),LK(NZ),
     *       LTEST
      COMMON /INDR/ RESMAX,ALFA
      COMMON /COEF/ AE(NXYZ),AW(NXYZ),AN(NXYZ),AS(NXYZ),
     *       AT(NXYZ),AB(NXYZ),AP(NXYZ),Q(NXYZ)
      DIMENSION FI(NXYZ),ZK(NXYZ),RES(NXYZ),D(NXYZ),PK(NXYZ)
C
C.....INITALIZE WORKING ARRAYS
C
      DO IJK=1,NIJK
        PK(IJK)=0.
        ZK(IJK)=0.
        D(IJK)=0.
        RES(IJK)=0.
      END DO
C
C.....CALCULATE INITIAL RESIDUAL VECTOR AND THE NORM
C
      RES0=0.
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            RES(IJK)=Q(IJK)-AE(IJK)*FI(IJK+NJ)-AW(IJK)*FI(IJK-NJ)-
     *        AN(IJK)*FI(IJK+1)-AS(IJK)*FI(IJK-1)-AT(IJK)*FI(IJK+NIJ)-
     *        AB(IJK)*FI(IJK-NIJ)-AP(IJK)*FI(IJK)
            RES0=RES0+ABS(RES(IJK))
          END DO
        END DO
      END DO
C
C.....IF LTEST=1, PRINT THE NORM 
C
      IF(LTEST.EQ.1) WRITE(6,*) 0,' SWEEP, RES0 = ',RES0
C
C.....CALCULATE ELEMENTS OF DIAGONAL PRECONDITIONING MATRIX
C
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            D(IJK)=1./(AP(IJK)-AW(IJK)**2*D(IJK-NJ)-AS(IJK)**2*D(IJK-1)
     *             -AB(IJK)**2*D(IJK-NIJ))
          END DO
        END DO
      END DO
C
      S0=1.E20
C
C....START INNER ITERATIONS
C
      DO L=1,NS
C
C.....SOLVE FOR ZK(IJK) -- FORWARD SUBSTITUTION
C
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=(RES(IJK)-AW(IJK)*ZK(IJK-NJ)-AS(IJK)*ZK(IJK-1)-
     *              AB(IJK)*ZK(IJK-NIJ))*D(IJK)
          END DO
        END DO
      END DO
C
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=ZK(IJK)/(D(IJK)+1.E-30)
          END DO
        END DO
      END DO
C
C..... BACKWARD SUBSTITUTION; CALCULATE SCALAR PRODUCT SK
C
      SK=0.
      DO K=NKM,2,-1
        DO I=NIM,2,-1
          DO J=NJM,2,-1
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=(ZK(IJK)-AE(IJK)*ZK(IJK+NJ)-AN(IJK)*ZK(IJK+1)-
     *               AT(IJK)*ZK(IJK+NIJ))*D(IJK)
            SK=SK+RES(IJK)*ZK(IJK)
          END DO
        END DO
      END DO
C
C.....CALCULATE BETA
C
      BET=SK/S0
C
C.....CALCULATE NEW SEARCH VECTOR PK
C
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            PK(IJK)=ZK(IJK)+BET*PK(IJK)
          END DO
        END DO
      END DO
C
C.... CALCULATE SCALAR PRODUCT (PK.A PK) AND ALPHA (OVERWRITE ZK)
C
      PKAPK=0.
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=AP(IJK)*PK(IJK)+AE(IJK)*PK(IJK+NJ)+
     *        AW(IJK)*PK(IJK-NJ)+AN(IJK)*PK(IJK+1)+AS(IJK)*PK(IJK-1)+
     *        AT(IJK)*PK(IJK+NIJ)+AB(IJK)*PK(IJK-NIJ)
            PKAPK=PKAPK+PK(IJK)*ZK(IJK)
          END DO
        END DO
      END DO
C
      ALF=SK/PKAPK
C
C.....CALCULATE VARIABLE CORRECTION, NEW RESIDUAL VECTOR, AND NORM
C
      RESL=0.
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            FI(IJK)=FI(IJK)+ALF*PK(IJK)
            RES(IJK)=RES(IJK)-ALF*ZK(IJK)
            RESL=RESL+ABS(RES(IJK))
          END DO
        END DO
      END DO
C
      S0=SK
C
C.....CHECK CONVERGENCE
C
      RSM=RESL/(RES0+1.E-30)
      IF(LTEST.EQ.1) WRITE(6,*) L,' SWEEP, RESL = ',RESL,' RSM = ',RSM
      IF(RSM.LT.RESMAX) RETURN
C
C.....END OF ITERATION LOOP
C
      END DO
C
      RETURN
      END
