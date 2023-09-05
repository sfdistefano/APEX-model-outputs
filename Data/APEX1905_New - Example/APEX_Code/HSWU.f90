      SUBROUTINE HSWU(CPWU,RGS)
!     APEX1905
!     THIS SUBPROGRAM DISTRIBUTES PLANT EVAPORATION THROUGH THE ROOT
!     ZONE AND CALCULATES ACTUAL PLANT WATER USE BASED ON SOIL WATER
!     AVAILABILITY.
      USE PARM
      BLM=S15(ISL,ISA)
      IF(Z(ISL,ISA)<=.5)BLM=PRMT(5)*S15(ISL,ISA)
      IF(ISL/=LID(1,ISA))THEN
          CALL CRGBD(RGS)
          CPWU=CPWU*RGS
      END IF
      SUM=EP(JJK,ISA)*(1.-EXP(-UB1(ISA)*GX/RD(JJK,ISA)))/UOB(ISA)
      TOS=36.*ECND(ISL,ISA)
      XX=LOG10(S15(ISL,ISA))
      X1=3.1761-1.6576*(LOG10(SWST(ISL,ISA))-XX)/(LOG10(FC(ISL,ISA))-XX)
      IF(X1<4.)THEN
          WTN=MAX(5.,10.**X1)
          XX=TOS+WTN
          IF(XX<5000.)THEN
              IF(SCLM(22)>0.)XX=MIN(XX,SCLM(22))
              F=1.-XX/(XX+EXP(SCRP(22,1)-SCRP(22,2)*XX))
              UW(ISL,ISA)=MIN(SUM-CPWU*AEP(JJK,ISA)-(1.-CPWU)*UX,SWST(ISL,ISA)-BLM)*F*RGS
              UW(ISL,ISA)=MAX(0.,UW(ISL,ISA))*SALF
          END IF    
      END IF
      UX=SUM
      RETURN
      END