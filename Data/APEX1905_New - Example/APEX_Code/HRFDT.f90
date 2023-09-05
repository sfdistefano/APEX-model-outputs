      SUBROUTINE HRFDT
!     APEX1905
!     THIS SUBPROGRAM DISTRIBUTES DAILY RAINFALL EXPONENTIALLY &
!     FURNISHES THE RAINFALL DISTRIBUTION AT TIME INTERVALS DTHY.
      USE PARM
      UPLM=.95
      QMN=.25
      BLM=.05
      R1=ATRI(BLM,QMN,UPLM,10)
      XK1=R1/9.21
      XK2=XK1*(1.-R1)/R1
      DURG=RFV(IRF(ISA))/(REP*(XK1+XK2))
      IF(DURG>24.)THEN
          DURG=24.
          REP=RFV(IRF(ISA))/(24.*(XK1+XK2))
      END IF
      XKP1=XK1*DURG
      XKP2=XK2*DURG
      TP=R1*DURG
      XX1=REP*XKP1
      XX2=REP*XKP2
      T1=0.
	  RFDT(1)=0.
	  I=1
      DO WHILE(T1<DURG)
          I=I+1
          T1=T1+DTHY
          IF(T1<TP)THEN
              RFDT(I)=XX1*EXP((T1-TP)/XKP1)
          ELSE
              RFDT(I)=XX2*(1.-EXP((TP-T1)/XKP2))+XX1
              X1=MIN(1.,T1/DURG)
              RFDT(I)=X1*RFV(IRF(ISA))+(1.-X1)*RFDT(I)
          END IF
      END DO
      NRF=I
      IF(KFL(26)>0)THEN
          IF(RFDT(I)>50.)WRITE(KW(26),8)(RFDT(J),J=1,NRF)
	      WRITE(KW(26),10)RFV(IRF(ISA)),RFDT(I)
      END IF
      RETURN
    8 FORMAT(20X,7F7.2)
   10 FORMAT(5X,'RFV= ',F8.3,' mm',2X,'SUM RF DST= ',F8.3,' mm')
      END