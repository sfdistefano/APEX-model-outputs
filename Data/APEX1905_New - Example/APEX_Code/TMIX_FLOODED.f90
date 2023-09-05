      SUBROUTINE TMIX_FLOODED(EE,DMX,NMIX,JNT)
!     APEX0806
!     THIS SUBPROGRAM MIXES N,AND P WITHIN THE PLOW DEPTH ACCORDING TO THE MIXING EFFICIENCY OF THE IMPLEMENT, 
!     CALCULATES THE CHANGE IN BULK DENSITY, AND INCREASE SEDIMENT CONCENTRAION IN THE PONDING WATER

      USE PARM 
      DIMENSION TST(100),DUM(MSL,MSA),YTP(8)

      !PADDY_STO(1,ISA) Water storage volume, mm
      !PADDY_STO(2,ISA) Sediment storage, t/ha
      !PADDY_STO(3,ISA) Min N storage, kg/ha
      !PADDY_STO(4,ISA) Min P storage, kg/ha
      !PADDY_STO(5,ISA) Org N storage, kg/ha
      !PADDY_STO(6,ISA) Org P storage, kg/ha
      !PADDY_STO(7,ISA) Ammonia N storage, kg/ha

      !Assume 10,000 mg/l sediment concentration after plowing (sort of mud...)
      !which is equivalent to 100 kg/ha-mm sediment.
      
      ISM=NDP+34
      LD1=LID(1,ISA)
      II=IHC(JT1)
      IF(Z(LD1,ISA)>=DMX)RETURN
      RCF(ISA)=1.
      TLMF(ISA)=0.
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
      
      !Sediment concentration increase due to web ploughing
      PADDY_STO(2,ISA) = PADDY_STO(1,ISA) * 100. ! 10,000 mg/l 
      
      TST(1) = EAJL(PADDY_STO(3,ISA),EE) !soluble N
      TST(29) = EAJL(PADDY_STO(4,ISA),EE) !mineral P
      TST(17) = EAJL(PADDY_STO(6,ISA),EE) !organic P
      TST(18) = EAJL(PADDY_STO(7,ISA),EE) !organic N
      
      DO J=1,NBSL(ISA)
         K=LID(J,ISA)
         UN(K,ISA)=ROK(K,ISA)
         ZZ=Z(K,ISA)-XX
          
         IF(Z(K,ISA)>=DMX)THEN 
            RTO=(DMX-XX)/ZZ
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
            !     TST(2)=EAJL(WHPN(K,ISA),RE)+TST(2)
            !     TST(3)=EAJL(WHSN(K,ISA),RE)+TST(3)
            TST(4)=EAJL(WBMN(K,ISA),RE)+TST(4)
            TST(5)=EAJL(WLSN(K,ISA),RE)+TST(5)
            TST(6)=EAJL(WLMN(K,ISA),RE)+TST(6)
            !     TST(7)=EAJL(WHPC(K,ISA),RE)+TST(7)
            !     TST(8)=EAJL(WHSC(K,ISA),RE)+TST(8)
            TST(9)=EAJL(WBMC(K,ISA),RE)+TST(9)
            TST(10)=EAJL(WLSC(K,ISA),RE)+TST(10)
            TST(11)=EAJL(WLMC(K,ISA),RE)+TST(11)
            TST(12)=EAJL(WLSLC(K,ISA),RE)+TST(12)
            TST(14)=EAJL(WLS(K,ISA),RE)+TST(14)
            TST(15)=EAJL(WLM(K,ISA),RE)+TST(15)
            TST(16)=EAJL(WLSL(K,ISA),RE)+TST(16)
            TST(17)=EAJL(WPO(K,ISA),RE)+TST(17)
            TST(18)=EAJL(WON(K,ISA),EE)+TST(18) ! Jaehak 2014
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
            GOTO 10
         ELSE
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
         !     TST(2)=EAJL(WHPN(K,ISA),EE)+TST(2)
         !     TST(3)=EAJL(WHSN(K,ISA),EE)+TST(3)
            TST(4)=EAJL(WBMN(K,ISA),EE)+TST(4)
            TST(5)=EAJL(WLSN(K,ISA),EE)+TST(5)
            TST(6)=EAJL(WLMN(K,ISA),EE)+TST(6)
         !     TST(7)=EAJL(WHPC(K,ISA),EE)+TST(7)
         !     TST(8)=EAJL(WHSC(K,ISA),EE)+TST(8)
            TST(9)=EAJL(WBMC(K,ISA),EE)+TST(9)
            TST(14)=EAJL(WLS(K,ISA),EE)+TST(14)
            TST(15)=EAJL(WLM(K,ISA),EE)+TST(15)
            TST(16)=EAJL(WLSL(K,ISA),EE)+TST(16)
            TST(10)=EAJL(WLSC(K,ISA),EE)+TST(10)
            TST(11)=EAJL(WLMC(K,ISA),EE)+TST(11)
            TST(12)=EAJL(WLSLC(K,ISA),EE)+TST(12)
            TST(17)=EAJL(WPO(K,ISA),EE)+TST(17)
            TST(18)=EAJL(WON(K,ISA),EE)+TST(18) ! Jaehak 2014
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
         ENDIF
      END DO          
      J=NBSL(ISA)
      DMX=Z(LID(NBSL(ISA),ISA),ISA)

10    K1=J-1
      ZZ = PADDY_STO(1,ISA) / 1000.
      DO I=1,ISM
          TST(I)=TST(I)/(DMX+ZZ)
      END DO
      
      !Update nutrient concentration in the ponding water
      PADDY_STO(3,ISA) = TST(1)*ZZ + PADDY_STO(3,ISA) !mineral N
      PADDY_STO(4,ISA) = TST(29)*ZZ + PADDY_STO(4,ISA) !mineral P
      PADDY_STO(6,ISA) = TST(17)*ZZ + PADDY_STO(6,ISA) !organic P
      PADDY_STO(7,ISA) = TST(18)*ZZ + PADDY_STO(7,ISA) !organic N
      
      XX=0.
      DO J=1,K1
         ISL=LID(J,ISA)
         ZZ=Z(ISL,ISA)-XX
         RT1=MIN(1.,WCMU(ISL,ISA)/WBMC(ISL,ISA))
         WNO3(ISL,ISA)=TST(1)*ZZ+WNO3(ISL,ISA)
         !     WHPN(ISL,ISA)=TST(2)*ZZ+WHPN(ISL,ISA)
         !     WHSN(ISL,ISA)=TST(3)*ZZ+WHSN(ISL,ISA)
         WBMN(ISL,ISA)=TST(4)*ZZ+WBMN(ISL,ISA)
         WLSN(ISL,ISA)=TST(5)*ZZ+WLSN(ISL,ISA)
         WLMN(ISL,ISA)=TST(6)*ZZ+WLMN(ISL,ISA)
         !     WHPC(ISL,ISA)=TST(7)*ZZ+WHPC(ISL,ISA)
         !     WHSC(ISL,ISA)=TST(8)*ZZ+WHSC(ISL,ISA)
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
               &YTP(1),YTP(1))
               IF(WLM(ISL,ISA)>YTP(2))CALL TMXL1(DMX,TST(15),WLM(ISL,ISA),&
               &YTP(2),YTP(2))
               IF(WLSL(ISL,ISA)>YTP(3))CALL TMXL1(DMX,TST(16),WLSL(ISL,ISA),&
               &YTP(3),YTP(3))
               IF(WLSC(ISL,ISA)>YTP(4))CALL TMXL1(DMX,TST(10),WLSC(ISL,ISA),&
               &YTP(4),YTP(4))
               IF(WLMC(ISL,ISA)>YTP(5))CALL TMXL1(DMX,TST(11),WLMC(ISL,ISA),&
               &YTP(5),YTP(5))
               IF(WLSLC(ISL,ISA)>YTP(6))CALL TMXL1(DMX,TST(12),WLSLC(ISL,ISA),&
               &YTP(6),YTP(6))
               IF(WLSN(ISL,ISA)>YTP(7))CALL TMXL1(DMX,TST(5),WLSN(ISL,ISA),&
               &YTP(7),YTP(7))
               IF(WLMN(ISL,ISA)>YTP(8))CALL TMXL1(DMX,TST(6),WLMN(ISL,ISA),&
               &YTP(8),YTP(8))
            END IF
         END IF
         WLSLNC(ISL,ISA)=WLSC(ISL,ISA)-WLSLC(ISL,ISA)
         RSD(ISL,ISA)=.001*(WLS(ISL,ISA)+WLM(ISL,ISA))
         WPO(ISL,ISA)=TST(17)*ZZ+WPO(ISL,ISA)
         Won(ISL,ISA)=TST(18)*ZZ+WON(ISL,ISA)
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
         WCMU(ISL,ISA)=MAX(1.E-5,WBMC(ISL,ISA)*RT1)
         XX=Z(ISL,ISA)
      END DO          
      XX=DMX-Z(LID(K1,ISA),ISA)
      RT1=MIN(1.,WCMU(K,ISA)/WBMC(K,ISA))
      WNO3(K,ISA)=WNO3(K,ISA)+TST(1)*XX
!     WHPN(K,ISA)=WHPN(K,ISA)+TST(2)*XX
!     WHSN(K,ISA)=WHSN(K,ISA)+TST(3)*XX
      WBMN(K,ISA)=WBMN(K,ISA)+TST(4)*XX
      WLSN(K,ISA)=WLSN(K,ISA)+TST(5)*XX
      WLMN(K,ISA)=WLMN(K,ISA)+TST(6)*XX
!     WHPC(K,ISA)=WHPC(K,ISA)+TST(7)*XX
!     WHSC(K,ISA)=WHSC(K,ISA)+TST(8)*XX
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
      WON(K,ISA)=WON(K,ISA)+TST(18)*XX
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
      RETURN
      END