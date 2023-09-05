      SUBROUTINE SOIL_05(KWX)
      ! APEX1307
	  ! THIS SUBPROGRAM OUTPUTS SOIL PROPERTIES IN THE TOP 0.05 M OF SOIL
	  USE PARM 
      DIMENSION STOT(27)
      STOT=0.
      XX=0.
      JJ=0
      DO I=1,NBSL(ISA)
          ISL=LID(I,ISA)
          IF(Z(ISL,ISA)<.05)THEN
              DXZ=Z(ISL,ISA)-XX
              RTO=1.
          ELSE
              DXZ=.05-XX
              RTO=DXZ/(Z(ISL,ISA)-XX)
              JJ=1
          END IF 
          STOT(1)=STOT(1)+PO(ISL,ISA)*RTO
          STOT(2)=STOT(2)+FC(ISL,ISA)*RTO
          STOT(3)=STOT(3)+S15(ISL,ISA)*RTO
          STOT(4)=STOT(4)+SWST(ISL,ISA)*RTO
          STOT(5)=STOT(5)+SATC(ISL,ISA)*DXZ
          STOT(6)=STOT(6)+HCL(ISL,ISA)*DXZ
          STOT(7)=STOT(7)+BDP(ISL,ISA)*DXZ
          STOT(9)=STOT(9)+SAN(ISL,ISA)*DXZ
          STOT(10)=STOT(10)+SIL(ISL,ISA)*DXZ
          STOT(11)=STOT(11)+CLA(ISL,ISA)*DXZ
          STOT(12)=STOT(12)+ROK(ISL,ISA)*DXZ
          STOT(13)=STOT(13)+RSD(ISL,ISA)*RTO
          STOT(14)=STOT(14)+PH(ISL,ISA)*DXZ
          STOT(15)=STOT(15)+SMB(ISL,ISA)*DXZ
          STOT(16)=STOT(16)+CEC(ISL,ISA)*DXZ
          STOT(17)=STOT(17)+ALS(ISL,ISA)*DXZ
          STOT(18)=STOT(18)+CAC(ISL,ISA)*DXZ
          STOT(19)=STOT(19)+PSP(ISL,ISA)*DXZ
          STOT(20)=STOT(20)+WPML(ISL,ISA)*RTO
          STOT(21)=STOT(21)+WPMA(ISL,ISA)*RTO
          STOT(22)=STOT(22)+WPMS(ISL,ISA)*RTO
          STOT(23)=STOT(23)+WPO(ISL,ISA)*RTO
          STOT(24)=STOT(24)+WNO3(ISL,ISA)*RTO
          STOT(25)=STOT(25)+WON(ISL,ISA)*RTO
          STOT(26)=STOT(26)+WOC(ISL,ISA)*RTO
          STOT(27)=STOT(27)+WT(ISL,ISA)*RTO
          IF(JJ>0)EXIT
          XX=Z(ISL,ISA)
      END DO
      STOT(8)=STOT(27)/500.
      DO K=1,4
          STOT(K)=STOT(K)/50.
      END DO   
      DO K=5,19
          IF(K==8.OR.K==13)CYCLE
          STOT(K)=STOT(K)/.05
      END DO
      DO K=20,25
          STOT(K)=1000.*STOT(K)/STOT(27)
      END DO 
      STOT(26)=.1*STOT(26)/STOT(27)
      Z5=.05
      WRITE(KWX,901)ISA,NBSA(ISA),IYR,IY,MO1,(LORG(ISL,ISA),ISL=1,2),&
      Z5,(STOT(K),K=1,26)
	  RETURN
  901 FORMAT(2I8,1X,3I4,I2,'/',I2,F8.2,4F8.3,4F8.2,4F8.1,F8.2,5F8.1,F8.2,6F8.0,F8.2)            
      END     