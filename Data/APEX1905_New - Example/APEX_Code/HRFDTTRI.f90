      SUBROUTINE HRFDTTRI
      ! THIS SUBPROGRAM GENERATES RAINFALL DISTRIBUTIONS AT DT=DTHY TIME 
      ! INTERVALS USING A TRIANGULAR DISTRIBUTION
      USE PARM
      T1=0.
      X1=RFV(IRF(ISA))
      PRT=2.*X1/DUR
      TP=.5*DUR
      RTO1=PRT/TP
      RTO2=PRT/(DUR-TP)
      RFDT(1)=0.
      Y1=0.
      I=1
      DO WHILE(T1<=DUR)
          I=I+1
          I1=I-1
          T1=T1+DTHY
          IF(T1<=TP)THEN
              Y2=T1*RTO1
              RFDT(I)=RFDT(I1)+.05*(Y1+Y2)
          ELSE
              Y2=(DUR-T1)*RTO2
              RFDT(I)=RFDT(I1)+.05*(Y1+Y2)
          END IF
          Y1=Y2
      END DO
      NRF=I     
      RTO=X1/RFDT(NRF)
      DO J=1,NRF
          RFDT(J)=RFDT(J)*RTO
      END DO
      RETURN
      END