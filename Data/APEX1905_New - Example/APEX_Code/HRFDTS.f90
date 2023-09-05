      SUBROUTINE HRFDTS
      ! THIS SUBPROGRAM GENERATES RAINFALL DISTRIBUTIONS AT DTHY TIME 
      ! INTERVALS USING AN S CURVE
      USE PARM
      T1=0.
      !DTX=DTHY/.12
      X1=RFV(IRF(ISA))
      RFDT(1)=0.
      I=1
      DO WHILE(T1<=DUR)
          I=I+1
          T1=T1+DTHY
          T2=1000.*T1/DUR
          IF(SCLM(26)>0.)T2=MIN(T2,SCLM(26))
          RFF=T2/(T2+EXP(SCRP(26,1)-SCRP(26,2)*T2))
          !I=I+1
          RFDT(I)=X1*RFF
          !IF(T1>100.)EXIT
      END DO
      NRF=I     
      RTO=X1/RFDT(NRF)
      DO J=1,NRF
          RFDT(J)=RFDT(J)*RTO
      END DO
      RETURN
      END