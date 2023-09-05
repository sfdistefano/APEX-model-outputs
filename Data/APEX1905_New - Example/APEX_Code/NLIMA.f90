      SUBROUTINE NLIMA(SB,DSB,C1,PH,ALS,OC,BSA)
!     APEX1905
!     THIS SUBPROGRAM ESTIMATES ALUMINUM SATURATION USING BASE SATURA-
!     TION, ORGANIC C, AND PH.
      SB=SB-DSB
      IF(SB<.02)SB=.02
      BSA=C1*SB
      IF(PH>5.6)THEN
          ALS=0.
      ELSE
          ALS=154.2-1.017*BSA-3.173*OC-14.23*PH
          IF(ALS<0.)THEN
              ALS=0.
          ELSE
              IF(ALS>95.)ALS=95.
          END IF
      END IF    
      RETURN
      END