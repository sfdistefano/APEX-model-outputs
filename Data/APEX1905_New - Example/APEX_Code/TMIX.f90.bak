      SUBROUTINE TMIX(EE,DMX,NMIX,IBMX)
!     APEX1905
!     THIS SUBPROGRAM MIXES N,P, AND CROP RESIDUE WITHIN THE PLOW DEPTH
!     ACCORDING TO THE MIXING EFFICIENCY OF THE IMPLEMENT, CALCULATES
!     THE CHANGE IN BULK DENSITY, CONVERTS STANDING RESIDUE TO FLAT
!     RESIDUE, AND ESTIMATES THE IMPLEMENT'S EFFECT ON RIDGE HEIGHT AND
!     INTERVAL.
      USE PARM 
      DIMENSION TST(100),DUM(MSL,MSA),XTZ(9),YTP(8)
      ISM=NDP+34
      LD1=LID(1,ISA)
      !IF(IYR==1972.AND.MO==10.AND.KDA==13)THEN
          !XTZ=0.
          !DO ISL=1,NBSL(ISA)
              !XTZ(1)=XTZ(1)+WNO3(ISL,ISA)
              !XTZ(2)=XTZ(2)+WNH3(ISL,ISA)
              !XTZ(3)=XTZ(3)+WLSN(ISL,ISA)
              !XTZ(4)=XTZ(4)+WLMN(ISL,ISA)
              !XTZ(5)=XTZ(5)+WBMN(ISL,ISA)
              !XTZ(6)=XTZ(6)+WHSN(ISL,ISA)
              !XTZ(7)=XTZ(7)+WHPN(ISL,ISA)
              !XTZ(8)=XTZ(8)+WNMU(ISL,ISA)
              !XTZ(9)=XTZ(9)+WNOU(ISL,ISA)
          !END DO
          !AD1=XTZ(1)+XTZ(2)+XTZ(3)+XTZ(4)+XTZ(5)+XTZ(6)+XTZ(7)+XTZ(8)+XTZ(9)
      !END IF 
      !TOT=0.
      IF(IBMX==0)THEN
          IF(NMIX==0)THEN
              II=IHC(JT1)
              IF(DMX<0..AND.II/=NHC(2).AND.II/=NHC(19).AND.II/=NHC(22))THEN
                  DO J=1,LC
                      IF(IDC(J)==NDC(7).OR.IDC(J)==NDC(8).OR.IDC(J)==NDC(10))CYCLE
                      IF(CPHT(J,ISA)+DMX<0..OR.STL(J,ISA)<.001)CYCLE
                      XX=MIN(1.,(CPHT(J,ISA)+DMX)/CPHT(J,ISA))
                      ZZ=XX*STD(J,ISA)
                      ZL=XX*STL(J,ISA)
                      STD(J,ISA)=MAX(0.,STD(J,ISA)-ZZ)
                      STL(J,ISA)=MAX(0.,STL(J,ISA)-ZL)
                      X1=1.-XX
                      SLAI(J,ISA)=SLAI(J,ISA)*X1
                      HU(J,ISA)=HU(J,ISA)*X1
                      STDL(J,ISA)=STDL(J,ISA)*X1
                      DX=DM(J,ISA)+1.E-10
                      X1=ZZ+ZL+ZO
                      ZZ=XX*STDN(J,ISA)
                      ZN=ZL*UN1(J,ISA)/DX
                      STDN(J,ISA)=MAX(0.,STDN(J,ISA)-ZZ)
                      X2=ZZ+ZN+ZX
                      CALL NCNSTD(X1,X2,LD1)
                      ZZ=XX*STDP(J,ISA)
                      STDP(J,ISA)=MAX(1.E-5,STDP(J,ISA)-ZZ)
                      ZP=ZL*UP1(J,ISA)/DX
                      FOP(LD1,ISA)=FOP(LD1,ISA)+ZZ+ZP+X1
                      DM(J,ISA)=DM(J,ISA)-ZL
                      UN1(J,ISA)=MAX(1.E-5,UN1(J,ISA)-ZN)
                      UP1(J,ISA)=UP1(J,ISA)-ZP
                      CPHT(J,ISA)=-DMX
                  END DO    
                  RETURN
              END IF    
          END IF
          IF(Z(LD1,ISA)>=DMX)RETURN
          RCF(ISA)=1.
          IF(RHT(JT1)<RHT(JT2).AND.TLD(JT2)>0.)THEN
              RHTT(ISA)=RHT(JT1)+(RHT(JT2)-RHT(JT1))*EXP(-DMX/TLD(JT2))
          ELSE
              RHTT(ISA)=RHT(JT1)
              RINT(ISA)=MAX(RIN(JT1),RINT(ISA))
          END IF
          F=1.-EXP(-56.9*DMX*EE)
          SUM=0.
          TOT=0.
          DO JJ=1,NCP(IRO(ISA),ISA)
              X1=STD(JJ,ISA)*F
              STD(JJ,ISA)=MAX(1.E-5,STD(JJ,ISA)-X1)
              SUM=SUM+X1
              XX=F*STDN(JJ,ISA)
              STDN(JJ,ISA)=MAX(0.,STDN(JJ,ISA)-XX)
              TOT=TOT+XX
              XX=F*STDL(JJ,ISA)
              STDL(JJ,ISA)=MAX(.1*STD(JJ,ISA),STDL(JJ,ISA)-XX)
              XX=F*STDP(JJ,ISA)
              STDP(JJ,ISA)=MAX(1.E-10,STDP(JJ,ISA)-XX)
              FOP(LD1,ISA)=FOP(LD1,ISA)+XX
          END DO
          CALL NCNSTD(SUM,TOT,LD1)
          FOP(LD1,ISA)=FOP(LD1,ISA)+XX
          RRUF(ISA)=RR(JT1)
          IF(EE>1.E-5)TLMF(ISA)=0 !Jaehak to improve wind erosion calculation, 2020

      END IF
      DO I=1,ISM
          TST(I)=0.
      END DO
      XX=0.
      YTP(1)=WLS(LD1,ISA)
      YTP(2)=WLM(LD1,ISA)
      YTP(3)=WLSL(LD1,ISA)
      YTP(4)=WLSC(LD1,ISA)
      YTP(5)=WLMC(LD1,ISA)
      YTP(6)=WLSLC(LD1,ISA)
      YTP(7)=WLSN(LD1,ISA)
      YTP(8)=WLMN(LD1,ISA)
      DO J=1,NBSL(ISA)
          K=LID(J,ISA)
          IF(IBMX>0)EE=WBMX
          UN(K,ISA)=ROK(K,ISA)
          ZZ=Z(K,ISA)-XX
          IF(Z(K,ISA)>=DMX)GO TO 8
          IF(NMIX==0)THEN
              BDP(K,ISA)=BDP(K,ISA)-(BDP(K,ISA)-.6667*BD(K,ISA))*EE
              CLA(K,ISA)=CLA(K,ISA)*ZZ
              SIL(K,ISA)=SIL(K,ISA)*ZZ
              ROK(K,ISA)=ROK(K,ISA)*ZZ
          END IF
          PMA=WPMA(K,ISA)+WPML(K,ISA)
          DUM(K,ISA)=PSP(K,ISA)*PMA
          UP(K,ISA)=PMA-DUM(K,ISA)
          TST(1)=EAJL(WNO3(K,ISA),EE)+TST(1)
          TST(2)=EAJL(WHPN(K,ISA),EE)+TST(2)
          TST(3)=EAJL(WHSN(K,ISA),EE)+TST(3)
          TST(4)=EAJL(WBMN(K,ISA),EE)+TST(4)
          TST(5)=EAJL(WLSN(K,ISA),EE)+TST(5)
          TST(6)=EAJL(WLMN(K,ISA),EE)+TST(6)
          TST(7)=EAJL(WHPC(K,ISA),EE)+TST(7)
          TST(8)=EAJL(WHSC(K,ISA),EE)+TST(8)
          TST(9)=EAJL(WBMC(K,ISA),EE)+TST(9)
          TST(14)=EAJL(WLS(K,ISA),EE)+TST(14)
          TST(15)=EAJL(WLM(K,ISA),EE)+TST(15)
          TST(16)=EAJL(WLSL(K,ISA),EE)+TST(16)
          TST(10)=EAJL(WLSC(K,ISA),EE)+TST(10)
          TST(11)=EAJL(WLMC(K,ISA),EE)+TST(11)
          TST(12)=EAJL(WLSLC(K,ISA),EE)+TST(12)
          IF(J==1)THEN
              YTP(1)=WLS(LD1,ISA)
              YTP(2)=WLM(LD1,ISA)
              YTP(3)=WLSL(LD1,ISA)
              YTP(4)=WLSC(LD1,ISA)
              YTP(5)=WLMC(LD1,ISA)
              YTP(6)=WLSLC(LD1,ISA)
              YTP(7)=WLSN(LD1,ISA)
              YTP(8)=WLMN(LD1,ISA)
          END IF
          TST(17)=EAJL(WPO(K,ISA),EE)+TST(17)
          TST(19)=EAJL(WPMA(K,ISA),EE)+TST(19)
          TST(20)=EAJL(WPOU(K,ISA),EE)+TST(20)
          TST(21)=EAJL(FOP(K,ISA),EE)+TST(21)
          TST(22)=EAJL(WPMS(K,ISA),EE)+TST(22)
          IF(NMIX==0)THEN
              TST(23)=EAJL(CLA(K,ISA),EE)+TST(23)
              TST(24)=EAJL(SIL(K,ISA),EE)+TST(24)
              TST(27)=EAJL(ROK(K,ISA),EE)+TST(27)
          END IF
          TST(25)=EAJL(DUM(K,ISA),EE)+TST(25)
          TST(26)=EAJL(UP(K,ISA),EE)+TST(26)
          TST(28)=EAJL(WNH3(K,ISA),EE)+TST(28)
          TST(29)=EAJL(WPML(K,ISA),EE)+TST(29)
          TST(30)=EAJL(WNOU(K,ISA),EE)+TST(30)
          TST(31)=EAJL(RSDM(K,ISA),EE)+TST(31)
          TST(32)=EAJL(WCOU(K,ISA),EE)+TST(32)
	      TST(33)=EAJL(WPMU(K,ISA),EE)+TST(33)
	      TST(34)=EAJL(WNMU(K,ISA),EE)+TST(34)
          I1=35
          DO I=1,NDP
              TST(I1)=EAJL(PSTZ(I,K,ISA),EE)+TST(I1)
              I1=I1+1
          END DO
          XX=Z(K,ISA)
      END DO          
      J=NBSL(ISA)
      DMX=Z(LID(NBSL(ISA),ISA),ISA)
      GO TO 10
    8 RTO=(DMX-XX)/ZZ
      RE=RTO*EE
      IF(NMIX==0)THEN
          BDP(K,ISA)=BDP(K,ISA)-(BDP(K,ISA)-.6667*BD(K,ISA))*RE
          CLA(K,ISA)=CLA(K,ISA)*ZZ
          SIL(K,ISA)=SIL(K,ISA)*ZZ
          ROK(K,ISA)=ROK(K,ISA)*ZZ
      END IF
      PMA=WPMA(K,ISA)+WPML(K,ISA)
      DUM(K,ISA)=PSP(K,ISA)*PMA
      UP(K,ISA)=PMA-DUM(K,ISA)
      TST(1)=EAJL(WNO3(K,ISA),RE)+TST(1)
      TST(2)=EAJL(WHPN(K,ISA),RE)+TST(2)
      TST(3)=EAJL(WHSN(K,ISA),RE)+TST(3)
      TST(4)=EAJL(WBMN(K,ISA),RE)+TST(4)
      TST(5)=EAJL(WLSN(K,ISA),RE)+TST(5)
      TST(6)=EAJL(WLMN(K,ISA),RE)+TST(6)
      TST(7)=EAJL(WHPC(K,ISA),RE)+TST(7)
      TST(8)=EAJL(WHSC(K,ISA),RE)+TST(8)
      TST(9)=EAJL(WBMC(K,ISA),RE)+TST(9)
      TST(10)=EAJL(WLSC(K,ISA),RE)+TST(10)
      TST(11)=EAJL(WLMC(K,ISA),RE)+TST(11)
      TST(12)=EAJL(WLSLC(K,ISA),RE)+TST(12)
      TST(14)=EAJL(WLS(K,ISA),RE)+TST(14)
      TST(15)=EAJL(WLM(K,ISA),RE)+TST(15)
      TST(16)=EAJL(WLSL(K,ISA),RE)+TST(16)
      TST(17)=EAJL(WPO(K,ISA),RE)+TST(17)
      TST(19)=EAJL(WPMA(K,ISA),RE)+TST(19)
      TST(20)=EAJL(WPOU(K,ISA),RE)+TST(20)
      TST(21)=EAJL(FOP(K,ISA),RE)+TST(21)
      TST(22)=EAJL(WPMS(K,ISA),RE)+TST(22)
      IF(NMIX==0)THEN
          TST(23)=EAJL(CLA(K,ISA),RE)+TST(23)
          TST(24)=EAJL(SIL(K,ISA),RE)+TST(24)
          TST(27)=EAJL(ROK(K,ISA),RE)+TST(27)
      END IF
      TST(25)=EAJL(DUM(K,ISA),RE)+TST(25)
      TST(26)=EAJL(UP(K,ISA),RE)+TST(26)
      TST(28)=EAJL(WNH3(K,ISA),RE)+TST(28)
      TST(29)=EAJL(WPML(K,ISA),RE)+TST(29)
      TST(30)=EAJL(WNOU(K,ISA),RE)+TST(30)
      TST(31)=EAJL(RSDM(K,ISA),RE)+TST(31)
      TST(32)=EAJL(WCOU(K,ISA),RE)+TST(32)
	  TST(33)=EAJL(WPMU(K,ISA),RE)+TST(33)
	  TST(34)=EAJL(WNMU(K,ISA),RE)+TST(34)
      I1=35
      DO I=1,NDP
          TST(I1)=EAJL(PSTZ(I,K,ISA),RE)+TST(I1)
          I1=I1+1
      END DO
   10 K1=J-1
      DO I=1,ISM
          TST(I)=TST(I)/DMX
      END DO
      XX=0.
      DO J=1,K1
          ISL=LID(J,ISA)
          ZZ=Z(ISL,ISA)-XX
          !RT1=MIN(1.,WCMU(ISL,ISA)/WBMC(ISL,ISA))    ! Commented by Luca Doro - 20191018
          WNO3(ISL,ISA)=TST(1)*ZZ+WNO3(ISL,ISA)
          WHPN(ISL,ISA)=TST(2)*ZZ+WHPN(ISL,ISA)
          WHSN(ISL,ISA)=TST(3)*ZZ+WHSN(ISL,ISA)
          WBMN(ISL,ISA)=TST(4)*ZZ+WBMN(ISL,ISA)
          WLSN(ISL,ISA)=TST(5)*ZZ+WLSN(ISL,ISA)
          WLMN(ISL,ISA)=TST(6)*ZZ+WLMN(ISL,ISA)
          WHPC(ISL,ISA)=TST(7)*ZZ+WHPC(ISL,ISA)
          WHSC(ISL,ISA)=TST(8)*ZZ+WHSC(ISL,ISA)
          WBMC(ISL,ISA)=TST(9)*ZZ+WBMC(ISL,ISA)
          WLSC(ISL,ISA)=TST(10)*ZZ+WLSC(ISL,ISA)
          WLMC(ISL,ISA)=TST(11)*ZZ+WLMC(ISL,ISA)
          WLSLC(ISL,ISA)=TST(12)*ZZ+WLSLC(ISL,ISA)
          WLS(ISL,ISA)=TST(14)*ZZ+WLS(ISL,ISA)
          WLM(ISL,ISA)=TST(15)*ZZ+WLM(ISL,ISA)
          WLSL(ISL,ISA)=TST(16)*ZZ+WLSL(ISL,ISA)
          IF(J==1)THEN
              IF(DMX>.01)THEN
                  IF(WLS(ISL,ISA)>YTP(1))CALL TMXL1(DMX,TST(14),WLS(ISL,ISA),&
                  YTP(1),YTP(1))
                  IF(WLM(ISL,ISA)>YTP(2))CALL TMXL1(DMX,TST(15),WLM(ISL,ISA),&
                  YTP(2),YTP(2))
                  IF(WLSL(ISL,ISA)>YTP(3))CALL TMXL1(DMX,TST(16),WLSL(ISL,ISA),&
                  YTP(3),YTP(3))
                  IF(WLSC(ISL,ISA)>YTP(4))CALL TMXL1(DMX,TST(10),WLSC(ISL,ISA),&
                  YTP(4),YTP(4))
                  IF(WLMC(ISL,ISA)>YTP(5))CALL TMXL1(DMX,TST(11),WLMC(ISL,ISA),&
                  YTP(5),YTP(5))
                  IF(WLSLC(ISL,ISA)>YTP(6))CALL TMXL1(DMX,TST(12),WLSLC(ISL,ISA),&
                  YTP(6),YTP(6))
                  IF(WLSN(ISL,ISA)>YTP(7))CALL TMXL1(DMX,TST(5),WLSN(ISL,ISA),&
                  YTP(7),YTP(7))
                  IF(WLMN(ISL,ISA)>YTP(8))CALL TMXL1(DMX,TST(6),WLMN(ISL,ISA),&
                  YTP(8),YTP(8))
              END IF
          END IF
          WLSLNC(ISL,ISA)=WLSC(ISL,ISA)-WLSLC(ISL,ISA)
          RSD(ISL,ISA)=.001*(WLS(ISL,ISA)+WLM(ISL,ISA))
          WPO(ISL,ISA)=TST(17)*ZZ+WPO(ISL,ISA)
          WPMA(ISL,ISA)=TST(19)*ZZ+WPMA(ISL,ISA)
          WPOU(ISL,ISA)=TST(20)*ZZ+WPOU(ISL,ISA)
          FOP(ISL,ISA)=TST(21)*ZZ+FOP(ISL,ISA)
          WPMS(ISL,ISA)=TST(22)*ZZ+WPMS(ISL,ISA)
          DUM(ISL,ISA)=TST(25)*ZZ+DUM(ISL,ISA)
          UP(ISL,ISA)=TST(26)*ZZ+UP(ISL,ISA)
          IF(NMIX==0)THEN
              ROK(ISL,ISA)=TST(27)+ROK(ISL,ISA)/ZZ
              CLA(ISL,ISA)=TST(23)+CLA(ISL,ISA)/ZZ
              SIL(ISL,ISA)=TST(24)+SIL(ISL,ISA)/ZZ
          END IF
          WNH3(ISL,ISA)=TST(28)*ZZ+WNH3(ISL,ISA)
          WPML(ISL,ISA)=TST(29)*ZZ+WPML(ISL,ISA)
	      WNOU(ISL,ISA)=TST(30)*ZZ+WNOU(ISL,ISA)
          RSDM(ISL,ISA)=TST(31)*ZZ+RSDM(ISL,ISA)
          WCOU(ISL,ISA)=TST(32)*ZZ+WCOU(ISL,ISA)
	      WPMU(ISL,ISA)=TST(33)*ZZ+WPMU(ISL,ISA)
	      WNMU(ISL,ISA)=TST(34)*ZZ+WNMU(ISL,ISA)
          I1=35
          DO I=1,NDP
              PSTZ(I,ISL,ISA)=TST(I1)*ZZ+PSTZ(I,ISL,ISA)
              I1=I1+1
          END DO
          PSP(ISL,ISA)=DUM(ISL,ISA)/(UP(ISL,ISA)+DUM(ISL,ISA))
          RX=MIN(1.,(100.-ROK(ISL,ISA))/(100.-UN(ISL,ISA)))
          FC(ISL,ISA)=FC(ISL,ISA)*RX
          S15(ISL,ISA)=S15(ISL,ISA)*RX
          PO(ISL,ISA)=PO(ISL,ISA)*RX
          CALL SPOFC(ISL)
          SAN(ISL,ISA)=100.-CLA(ISL,ISA)-SIL(ISL,ISA)
          WT(ISL,ISA)=BD(ISL,ISA)*ZZ*1.E4
          RT1=MIN(1.,WCMU(ISL,ISA)/WBMC(ISL,ISA))    ! Moved by Luca Doro - 20191018
          WCMU(ISL,ISA)=MAX(1.E-5,WBMC(ISL,ISA)*RT1)
          XX=Z(ISL,ISA)
      END DO          
      XX=DMX-Z(LID(K1,ISA),ISA)
      RT1=MIN(1.,WCMU(K,ISA)/WBMC(K,ISA))
      WNO3(K,ISA)=WNO3(K,ISA)+TST(1)*XX
      WHPN(K,ISA)=WHPN(K,ISA)+TST(2)*XX
      WHSN(K,ISA)=WHSN(K,ISA)+TST(3)*XX
      WBMN(K,ISA)=WBMN(K,ISA)+TST(4)*XX
      WLSN(K,ISA)=WLSN(K,ISA)+TST(5)*XX
      WLMN(K,ISA)=WLMN(K,ISA)+TST(6)*XX
      WHPC(K,ISA)=WHPC(K,ISA)+TST(7)*XX
      WHSC(K,ISA)=WHSC(K,ISA)+TST(8)*XX
      WBMC(K,ISA)=WBMC(K,ISA)+TST(9)*XX
      WLSC(K,ISA)=WLSC(K,ISA)+TST(10)*XX
      WLMC(K,ISA)=WLMC(K,ISA)+TST(11)*XX
      WLSLC(K,ISA)=WLSLC(K,ISA)+TST(12)*XX
      WLS(K,ISA)=WLS(K,ISA)+TST(14)*XX
      WLM(K,ISA)=WLM(K,ISA)+TST(15)*XX
      WLSL(K,ISA)=WLSL(K,ISA)+TST(16)*XX
      WLSLNC(K,ISA)=WLSC(K,ISA)-WLSLC(K,ISA)
      RSD(K,ISA)=.001*(WLS(K,ISA)+WLM(K,ISA))
      WPO(K,ISA)=WPO(K,ISA)+TST(17)*XX
      WPMA(K,ISA)=WPMA(K,ISA)+TST(19)*XX
      WPOU(K,ISA)=WPOU(K,ISA)+TST(20)*XX
      FOP(K,ISA)=FOP(K,ISA)+TST(21)*XX
      WPMS(K,ISA)=WPMS(K,ISA)+TST(22)*XX
      DUM(K,ISA)=DUM(K,ISA)+TST(25)*XX
      UP(K,ISA)=UP(K,ISA)+TST(26)*XX
      IF(NMIX==0)THEN
          ROK(K,ISA)=ROK(K,ISA)+TST(27)*XX
          CLA(K,ISA)=CLA(K,ISA)+TST(23)*XX
          SIL(K,ISA)=SIL(K,ISA)+TST(24)*XX
      END IF
      WNH3(K,ISA)=WNH3(K,ISA)+TST(28)*XX
      WCMU(K,ISA)=MAX(1.E-5,WBMC(K,ISA)*RT1)
      WPML(K,ISA)=WPML(K,ISA)+TST(29)*XX
	  WNOU(K,ISA)=WNOU(K,ISA)+TST(30)*XX
      RSDM(K,ISA)=RSDM(K,ISA)+TST(31)*XX
      WCOU(K,ISA)=WCOU(K,ISA)+TST(32)*XX
	  WPMU(K,ISA)=WPMU(K,ISA)+TST(33)*XX
	  WNMU(K,ISA)=WNMU(K,ISA)+TST(34)*XX
      I1=35
      DO I=1,NDP
          PSTZ(I,K,ISA)=PSTZ(I,K,ISA)+TST(I1)*XX
          I1=I1+1
      END DO
      PSP(K,ISA)=DUM(K,ISA)/(UP(K,ISA)+DUM(K,ISA))
      ZZ=Z(K,ISA)-Z(LID(K1,ISA),ISA)
      IF(NMIX==0)THEN
          ROK(K,ISA)=ROK(K,ISA)/ZZ
          CLA(K,ISA)=CLA(K,ISA)/ZZ
          SIL(K,ISA)=SIL(K,ISA)/ZZ
      END IF
      IF(UN(K,ISA)>0.)THEN
          RX=MIN(1.,(100.-ROK(K,ISA))/(100.-UN(K,ISA)))
          FC(K,ISA)=FC(K,ISA)*RX
          S15(K,ISA)=S15(K,ISA)*RX
          PO(K,ISA)=PO(K,ISA)*RX
          CALL SPOFC(K)
      END IF
      SAN(K,ISA)=100.-CLA(K,ISA)-SIL(K,ISA)
      WT(K,ISA)=BD(K,ISA)*ZZ*1.E4
      !XTX=0.
      !DO I=1,NBSL(ISA)
          !ISL=LID(I,ISA)
          !XTX(1)=XTX(1)+WPML(ISL,ISA)
          !XTX(2)=XTX(2)+WPO(ISL,ISA)
          !XTX(3)=XTX(3)+FOP(ISL,ISA)
          !XTX(4)=XTX(4)+WPOU(ISL,ISA)
          !XTX(5)=XTX(5)+WPMU(ISL,ISA)
          !XTX(6)=XTX(6)+WPMA(ISL,ISA)
          !XTX(7)=XTX(7)+WPMS(ISL,ISA)
          !XTX(8)=XTX(8)+WNOU(ISL,ISA)
      !END DO
      !XTX(3)=XTX(3)+STDP(JJK,ISA)+UP1(JJK,ISA)
      !DO I=1,8
          !DF=XTZ(I)-XTX(I)
          !IF(ABS(DF)>.001)WRITE(KW(1),'(A,4I4,20E16.6)')'#####',IY,MO,KDA,I,XTZ(I),XTX(I),DF
      !END DO
!     IF(EE<.9)RETURN
!     LD2=LID(2,ISA)
!     RSDM(LD2,ISA)=RSDM(LD2,ISA)+RSDM(LD1,ISA)
!     RSDM(LD1,ISA)=1.E-5
!     X1=STL(JJK,ISA)+STD(JJK,ISA)
!     DM(JJK,ISA)=DM(JJK,ISA)-STL(JJK,ISA)
!     XX=STL(JJK,ISA)/(DM(JJK,ISA)+1.E-10)
!     X2=XX*UN1(JJK,ISA)
!     X3=XX*UP1(JJK,ISA)
!     X4=STDN(JJK,ISA)+X2
!     CALL NCNSTD(X1,X4,LD1)
!     WLMN(LD2,ISA)=WLMN(LD2,ISA)+WLMN(LD1,ISA)
!     WLMN(LD1,ISA)=1.E-5
!     WLSN(LD2,ISA)=WLSN(LD2,ISA)+WLSN(LD1,ISA)
!     WLSN(LD1,ISA)=1.E-5
!     WLS(LD2,ISA)=WLS(LD2,ISA)+WLS(LD1,ISA)
!     WLS(LD1,ISA)=1.E-5
!     WLM(LD2,ISA)=WLM(LD2,ISA)+WLM(LD1,ISA)
!     WLM(LD1,ISA)=1.E-5
!     WLSL(LD2,ISA)=WLSL(LD2,ISA)+WLSL(LD1,ISA)
!     WLSL(LD1,ISA)=1.E-5
!     WLSC(LD2,ISA)=WLSC(LD2,ISA)+WLSC(LD1,ISA)
!     WLSC(LD1,ISA)=1.E-5
!     WLMC(LD2,ISA)=WLMC(LD2,ISA)+WLMC(LD1,ISA)
!     WLMC(LD1,ISA)=1.E-5
!     WLSLC(LD2,ISA)=WLSLC(LD2,ISA)+WLSLC(LD1,ISA)
!     WLSLC(LD1,ISA)=1.E-5
!     WLSLNC(LD2,ISA)=WLSLNC(LD2,ISA)+WLSLNC(LD1,ISA)
!     WLSLNC(LD1,ISA)=1.E-5
!     WNO3(LD2,ISA)=WNO3(LD2,ISA)+WNO3(LD1,ISA)
!     WNO3(LD1,ISA)=1.E-5
!     WNH3(LD2,ISA)=WNH3(LD2,ISA)+WNH3(LD1,ISA)
!     WNH3(LD1,ISA)=1.E-5
!     WCMU(LD2,ISA)=WCMU(LD2,ISA)+WCMU(LD1,ISA)
!     WCMU(LD1,ISA)=1.E-5
!     WNMU(LD2,ISA)=WNMU(LD2,ISA)+WNMU(LD1,ISA)
!     WNMU(LD1,ISA)=1.E-5
!     WPMA(LD2,ISA)=WPMA(LD2,ISA)+WPMA(LD1,ISA)
!     WPMA(LD1,ISA)=1.E-5
!     STL(JJK,ISA)=0.
!     STD(JJK,ISA)=0.
!     STDN(JJK,ISA)=0.
!     STDL(JJK,ISA)=0.
!     STDP(JJK,ISA)=0.
!     UN1(JJK,ISA)=MAX(1.E-5,UN1(JJK,ISA)-X2)
!     UP1(JJK,ISA)=UP1(JJK,ISA)-X3
!     DM(JJK,ISA)=DM(JJK,ISA)-STL(JJK,ISA)
!     STL(JJK,ISA)=0.
!     STD(JJK,ISA)=0.
!     X1=STL(JJK,ISA)/(DM(JJK,ISA)+1.E-10)
!     X2=X1*UN1(JJK,ISA)
!     X3=X1*UP1(JJK,ISA)
!     WNOU(LD2,ISA)=WNOU(LD2,ISA)+WNOU(LD1,ISA)
!     WNOU(LD1,ISA)=1.E-5
!     WCOU(LD2,ISA)=WCOU(LD2,ISA)+WCOU(LD1,ISA)
!     WCOU(LD1,ISA)=1.E-5
!     STDN(JJK,ISA)=0.
!     FOP(LD2,ISA)=FOP(LD2,ISA)+FOP(LD1,ISA)+STDP(JJK,ISA)+X3
!     FOP(LD1,ISA)=1.E-5
!     WPOU(LD2,ISA)=WPOU(LD2,ISA)+WPOU(LD1,ISA)
!     WPOU(LD1,ISA)=1.E-5
!     WPML(LD2,ISA)=WPML(LD2,ISA)+WPML(LD1,ISA)
!     WPML(LD1,ISA)=1.E-5
!     WPMU(LD2,ISA)=WPMU(LD2,ISA)+WPMU(LD1,ISA)
!     WPMU(LD1,ISA)=1.E-5
!     STDP(JJK,ISA)=0.
!     UN1(JJK,ISA)=MAX(1.E-5,UN1(JJK,ISA)-X2)
!     UP1(JJK,ISA)=UP1(JJK,ISA)-X3
      !IF(IYR==1972.AND.MO==10.AND.KDA==13)THEN
          !WRITE(KW(1),'(5X,10E16.6)')XTZ,AD1
          !XTZ=0.
          !DO ISL=1,NBSL(ISA)
              !XTZ(1)=XTZ(1)+WNO3(ISL,ISA)
              !XTZ(2)=XTZ(2)+WNH3(ISL,ISA)
              !XTZ(3)=XTZ(3)+WLSN(ISL,ISA)
              !XTZ(4)=XTZ(4)+WLMN(ISL,ISA)
              !XTZ(5)=XTZ(5)+WBMN(ISL,ISA)
              !XTZ(6)=XTZ(6)+WHSN(ISL,ISA)
              !XTZ(7)=XTZ(7)+WHPN(ISL,ISA)
              !XTZ(8)=XTZ(8)+WNMU(ISL,ISA)
              !XTZ(9)=XTZ(9)+WNOU(ISL,ISA)
          !END DO
          !AD2=XTZ(1)+XTZ(2)+XTZ(3)+XTZ(4)+XTZ(5)+XTZ(6)+XTZ(7)+XTZ(8)+XTZ(9)
          !WRITE(KW(1),'(5X,10E16.6)')XTZ,AD2
          !DF=AD1+TOT-AD2
          !IF(ABS(DF)>.1)WRITE(KW(1),'(A,3I4,20E16.6)')'#####',IY,MO,KDA,&
          !AD1,TOT,AD2,DF
      !END IF        
      RETURN
      END