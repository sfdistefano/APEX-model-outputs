      SUBROUTINE ERHEM
!     APEX1905
!     THIS SUBPROGRAM PREDICTS DAILY SOIL LOSS FOR RANGELANDS USING 
!     EROSION EQS TAKEN FROM RHEM.  THIS MODIFICATION OF RHEM IS CALLED
!     REMX
      USE PARM
      DATA H2OC/9807./
      LD1=LID(1,ISA)
      DTS=3600.*DTHY
      IF(INFL==0.AND.IRFT==0)CALL HRFDT
      X2=1.-EXP(-.03*ROK(LD1,ISA))
      XX=1.-FRSD      
      UBK=7.747E-6
      X1=MAX(FGC(ISA),FGSL(ISA))
      IF(IREM==0)THEN
          IF(FBAR(ISA)>.5)THEN
              SSKG=1.3*10**(3.2727-.347*X1-.5226*XX)*RSF(ISA)
              SSKB=1.3*10**(4.1438-.347*X1-.5226*XX)*RSF(ISA)
              SSK=FBAR(ISA)*SSKB+(1.-FBAR(ISA))*SSKG
          ELSE
              SSK=1.3*10**(3.2727-.347*X1-.5226*XX)*RSF(ISA)
          END IF
      ELSE
          SSK=.1*YSSK/(DTS*QVOL(IDO))
      END IF
      PKQ=0.
      XX=0.
      BRFX=0.
      TOT=0.
      AD1=0.
      T1=0.
      SPB=H2OC*STP(ISA) 
      W0=2.46/STP(ISA)**.4
      DTSL=DTS/SPLG(ISA)
      DO I=2,NRF
          T1=T1+DTHY
          IF(INFL==1.OR.INFL==4)THEN
              QMM=QGA(I)
          ELSE
              X1=RFDT(I)-SCN2
              IF(X1>0.)THEN
                  QV=X1*X1/(RFDT(I)+.8*SCN)
              ELSE
                  QV=0.
              END IF
              QMM=QV-XX
          END IF
          IF(QMM>0.)THEN
              PKQ=MAX(QMM,PKQ)
              !QM3=10.*QMM
              QM3=QMM*SPLG(ISA)/1000.
              QCMS=QM3/DTS
              W=MIN(1.,W0*QCMS**.39)
              QCMSR=QCMS/W
              SP=SPB*QCMSR
              X2=EXP(.845+.412*LOG10(1000.*SP))
              TFC=.1*W*10**(-34.47+38.61*X2/(1.+X2))
              IF(IREM==0)THEN
                  XR=(RFDT(I)-RFDT(I-1))/(1000.*DTS)
                  X2=.001*QMM/DTS
                  DSS=SSK*XR**1.052*X2**.592
              ELSE
                  DSS=SSK*QMM
              END IF
              DCF=.5*SP*UBK
              DTOT=DSS+DCF
              X3=DTOT*DTS*W
              X4=TFC*DTSL
              DY=MIN(X3,X4)
              AD1=AD1+DY
              XX=QV
              BRFX=MAX(QCMS,BRFX)
              TOT=TOT+QMM
          END IF
      END DO
      PKQ=PKQ/DTHY
      YSD(6,IDO)=MAX(0.,10.*AD1)*WSA(ISA)
      !WRITE(KW(1),'(T10,A,F8.3,A,2X,A,F8.3,A)')'PKQ=',PKQ,'mm/h','YSD=',&
      !YSD(6,IDO),'t/ha'                               
      RETURN
      END