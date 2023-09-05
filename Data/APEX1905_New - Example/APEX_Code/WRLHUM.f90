      SUBROUTINE WRLHUM(I)
!     APEX1905
!     THIS SUBPROGRAM SIMULATES DAILY RELATIVE HUMIDITY FROM TRIANGULAR
!     DISTRIBUTION.
      USE PARM 
      Q1=RHM-1.
      UPLM=RHM-Q1*EXP(Q1)
      BLM=RHM*(1.-EXP(-RHM))
      RHD(I)=ATRI(BLM,RHM,UPLM,IDG(7))
      RETURN
      END