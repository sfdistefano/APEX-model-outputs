      SUBROUTINE HUSE
!     APEX1905
!     THIS SUBPROGRAM IS THE MASTER WATER AND NUTRIENT USE SUBROUTINE.
!     CALLS HSWU AND NUPPO FOR EACH SOIL LAYER.
      USE PARM
      LRD(ISA)=0
      L1=0
      UX=0.
      SEP=0.
      SUM=0.
      TOT=0.
      CPWU=1.
      RGS=1.
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          SUM=SUM+SWST(ISL,ISA)-FC(ISL,ISA)
          TOT=TOT+PO(ISL,ISA)-FC(ISL,ISA)
          SEP=Z(ISL,ISA)
          IF(L1>0)CYCLE
          IF(RD(JJK,ISA)>Z(ISL,ISA))THEN
              GX=Z(ISL,ISA)
          ELSE
              GX=RD(JJK,ISA)
              LRD(ISA)=MAX(LRD(ISA),J)
              L1=J
          END IF
          CALL HSWU(CPWU,RGS)
          AEP(JJK,ISA)=AEP(JJK,ISA)+UW(ISL,ISA)
      END DO 
      IF(LRD(ISA)==0)LRD(ISA)=NBSL(ISA)
      RTO=MIN(1.,SUM/TOT)
      F=100.*(RTO-CAF(JJK))/(1.0001-CAF(JJK))
      IF(F>0.)THEN
          IF(SCLM(7)>0.)F=MIN(F,SCLM(7))
          SAT=1.-F/(F+EXP(SCRP(7,1)-SCRP(7,2)*F))
      ELSE
          SAT=1.
      END IF
      RETURN
      END