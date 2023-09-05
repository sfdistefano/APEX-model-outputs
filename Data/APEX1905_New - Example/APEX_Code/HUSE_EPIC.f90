      SUBROUTINE HUSE_EPIC
!     EPIC0810
!     THIS SUBPROGRAM IS THE MASTER WATER AND NUTRIENT USE SUBPROGRAM.
!     CALLS HSWU AND NUPPO FOR EACH SOIL LAYER.
      USE PARM
	  LRD(ISA)=1
      IAR=0
      UX=0.
      XX=0.
      TOT=0.
      SAT=0.
      CPWU=1.
      RGS=1.
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          IF(Z(ISL,ISA)<1.)THEN
              SAT=SAT+SWST(ISL,ISA)
              TOT=TOT+PO(ISL,ISA)
              XX=Z(ISL,ISA)
          ELSE
              IF(IAR==0)THEN
                  IAR=1
                  X3=1.-XX
                  X4=Z(ISL,ISA)-XX
                  RTO=X3/X4
                  IF(WTBL(ISA)<=Z(ISL,ISA))THEN
                      X1=PO(ISL,ISA)*(Z(ISL,ISA)-WTBL(ISA))/X4
                      X2=SWST(ISL,ISA)-X1
                      IF(WTBL(ISA)>1.)THEN
                          SAT=SAT+X2*X3/(WTBL(ISA)-XX)
                      ELSE
                          SAT=SAT+X2+PO(ISL,ISA)*(1.-WTBL(ISA))/X4
                      END IF
                  ELSE    
                      SAT=SAT+RTO*SWST(ISL,ISA)
                  END IF    
                  TOT=TOT+RTO*PO(ISL,ISA)
              END IF  
          END IF
          IF(LRD(ISA)>1)CYCLE
          IF(RD(JJK,ISA)>Z(ISL,ISA))THEN
              GX=Z(ISL,ISA)
          ELSE
              GX=RD(JJK,ISA)
              LRD(ISA)=MAX(LRD(ISA),J)
          END IF
          CALL HSWU(CPWU,RGS)
          AEP(JJK,ISA)=AEP(JJK,ISA)+UW(ISL,ISA)
      END DO
      IF(LRD(ISA)==0)LRD(ISA)=NBSL(ISA)
      IF(RZSW(ISA)>PAW(ISA))THEN
          RTO=MIN(1.,SAT/TOT)
          F=100.*(RTO-CAF(JJK))/(1.0001-CAF(JJK))
          IF(F>0.)THEN
              IF(SCLM(7)>0.)F=MIN(F,SCLM(7))
              SAT=1.-F/(F+EXP(SCRP(7,1)-SCRP(7,2)*F))
          ELSE
              SAT=1.
          END IF
      ELSE
          SAT=1.
      END IF
      RETURN
      END