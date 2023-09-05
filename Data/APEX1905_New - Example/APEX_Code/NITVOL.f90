      SUBROUTINE NITVOL(Z5)
!     APEX1905
!     THIS SUBPROGRAM SIMULATES THE TRANSFORMATION FROM NH3 TO NO3, AND
!     THE VOLATILIZATION OF NH3 USING MODIFIED METHODS OF REDDY AND OF
!     THE CERES MODEL
      USE PARM
      X1=.41*((STMP(ISL,ISA)-5.)/10.)
      WSAX=WSA(ISA)
      IF(X1>0.)THEN
          IF(ISL==LID(1,ISA))THEN
              FAF=.335+.16*LOG(U10(IRF(ISA))+.2)
              AKAV=X1*FAF
          ELSE
              FCEC=MAX(PRMT(27),1.-.038*CEC(ISL,ISA))
              X2=Z5
              IF(SCLM(12)>0.)X2=MIN(X2,SCLM(12))
              FZ=1.-.5*(X2/(X2+EXP(SCRP(12,1)-SCRP(12,2)*X2)))
              AKAV=X1*FCEC*FZ
          END IF
          IF(PH(ISL,ISA)>7.)THEN
              IF(PH(ISL,ISA)>7.4)THEN
                  FPH=5.367-.599*PH(ISL,ISA)
              ELSE
                  FPH=1.
              END IF
          ELSE
              FPH=.307*PH(ISL,ISA)-1.269
          END IF    
          AKAN=X1*SUT(ISL,ISA)*FPH
          AKAV=AKAV*SUT(ISL,ISA)
          XX=AKAV+AKAN
          IF(XX>0.)THEN
              F=MIN(PRMT(80),1.-EXP(-XX))
              X1=F*WNH3(ISL,ISA)
              AVOL=X1*PRMT(72)
              SVOL=SVOL+AVOL
              IF(NTV==0)THEN
                  RNIT=X1-AVOL
                  WNH3(ISL,ISA)=WNH3(ISL,ISA)-AVOL-RNIT
                  WNO3(ISL,ISA)=WNO3(ISL,ISA)+RNIT
                  SNIT=SNIT+RNIT
              ELSE
                  F=MIN(1.,PRMT(80)*AKAN)
                  GNO2=F*WNH3(ISL,ISA)
                  WNH3(ISL,ISA)=WNH3(ISL,ISA)-GNO2-AVOL
                  WNO2(ISL,ISA)=WNO2(ISL,ISA)+GNO2
                  IF(PH(ISL,ISA)>5.5)THEN
                      IF(PH(ISL,ISA)>7.2)THEN
                          FPH=4.367-.5324*PH(ISL,ISA)
                      ELSE
                          FPH=1.
                      END IF
                  ELSE
                      FPH=.307*PH(ISL,ISA)-1.269
                  END IF
                  AKAN=X1*SUT(ISL,ISA)*FPH
                  F=MIN(1.,PRMT(80)*AKAN)
                  GNO3=F*WNO2(ISL,ISA)
                  WNO2(ISL,ISA)=WNO2(ISL,ISA)-GNO3
                  WNO3(ISL,ISA)=WNO3(ISL,ISA)+GNO3
                  SNIT=SNIT+GNO2
              END IF    
          END IF    
      END IF    
      RETURN
      END