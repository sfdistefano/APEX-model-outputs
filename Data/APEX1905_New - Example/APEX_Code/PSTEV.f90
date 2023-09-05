      SUBROUTINE PSTEV
!     APEX1905
!     THIS SUBPROGRAM ESTIMATES UPWARD PESTICIDE MOVEMENT CAUSED BY SOIL 
!     EVAPORATION.
      USE PARM 
      IF(NEV==1)RETURN
      LD1=LID(1,ISA)
      DO K=1,NDP
          XX=0.
          VPST=0.
          KK=NEV
          DO KK=NEV,2,-1
              ISL=LID(KK,ISA)
              Y1=PSTZ(K,ISL,ISA)
              IF(Y1<.001)CYCLE
              V=SEV(ISL,ISA)
              IF(V<=0.)CYCLE
              VPST=Y1*(1.-EXP(-V/(PO(ISL,ISA)-S15(ISL,ISA)+.001*PKOC(K)*&
              WT(ISL,ISA)*WOC(ISL,ISA))))
              XX=XX+VPST
              PSTZ(K,ISL,ISA)=Y1-VPST
          END DO
          PSTZ(K,LD1,ISA)=PSTZ(K,LD1,ISA)+XX
      END DO
      RETURN
      END