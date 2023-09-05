      SUBROUTINE NCONC(P0,P5,P1,A,KW,MSO)
!     APEX1905
!     THIS SUBPROGRAM COMPUTES PARAMETERS OF AN EQUATION DESCRIBING THE
!     N AND P RELATIONS TO CROP MATURITY.
      DIMENSION KW(MSO)
      A=5.0
      DO I=1,10
          A5=A*.5
          EA=EXP(A)
          EA1=EA-1.
          EG=EXP(-A5)
          P0G=P0*EG
          EG1=EXP(A5)
          PEG=P1*(EA-EG1)
          P01=P0*(1.-EG)
          X1=PEG-P01
          PG5=.5*P0G
          FU=X1/EA1+P0G-P5
          IF (ABS(FU)<1.E-7) GO TO 3
          DFDA=(EA1*(P1*(EA-.5*EG1)-PG5)-EA*X1)/(EA1*EA1)-PG5
          A=A-FU/DFDA
      END DO
      WRITE (KW(1),4) A,FU
    3 P5=(P1*EA-P0)/EA1
      P0=P0-P5
      RETURN
    4 FORMAT (//T10,'NCONC DID NOT CONVERGE',2E16.6)
      END