      SUBROUTINE HWTBL
!     APEX1905
!     THIS SUBPROGRAM SIMULATES WATER TABLE DYNAMICS AS A FUNCTION
!     OF RAIN AND EVAP.
      USE PARM
      RTO=(SMRF(ISA)-SMEO(ISA))/SMEO(ISA)
      IF(RTO>0.)THEN
          XX=WTMN(ISA)
          X1=1.
      ELSE
          XX=WTMX(ISA)
          X1=PRMT(87)*(REAL(IDA)/REAL(ND))**PRMT(89)
      END IF
      X2=MIN(P88MX,ABS(RTO)*X1)
      WTBP(ISA)=(WTBP(ISA)-(WTBP(ISA)-XX)*X2)
      SUM=0.
      TOT=0.
      ISL=LID(NBSL(ISA),ISA)
      IF(WTBP(ISA)<=Z(ISL,ISA))THEN
          SWST(ISL,ISA)=SWST(ISL,ISA)+QSSI
          DO K=NBSL(ISA),2,-1
              ISL=LID(K,ISA)
              L1=LID(K-1,ISA)
              X1=SWST(ISL,ISA)-PO(ISL,ISA)
              IF(X1>0.)THEN
                  SWST(L1,ISA)=SWST(L1,ISA)+X1
                  SWST(ISL,ISA)=PO(ISL,ISA)
              ELSE
                  WTBL(ISA)=MAX(WTBP(ISA),Z(ISL,ISA)-SWST(ISL,ISA)*(Z(ISL,ISA)-Z(L1,ISA))/PO(ISL,ISA))
                  EXIT
              END IF
          END DO    
          SMM(19,MO,ISA)=SMM(19,MO,ISA)+QSSI
          VAR(19,ISA)=QSSI
          QIN(MO,ISA)=QIN(MO,ISA)+QSSI
      ELSE
          WTBL(ISA)=WTBP(ISA)
      END IF
      RETURN
      END