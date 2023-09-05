      SUBROUTINE CAGRO
!     APEX1905
!     THIS SUBPROGRAM CALCULATES THE DAILY INCREASE IN PLANT BIOMASS,
!     ROOT WEIGHT, AND YIELD BY ADJUSTING THE POTENTIAL VALUES WITH THE
!     ACTIVE STRESS CONSTRAINT.
      USE PARM
      REAL:: DFL(MNC,MSA)
 
      DFL(JJK,ISA)=0.				
 
      LD1=LID(1,ISA)
	  IPH(JJK,ISA)=1				
      II=IDC(JJK)
      X1=REG(JJK,ISA)*AJWA(ISA)*SHRL(ISA)
      RWL=RW(JJK,ISA)
      STLL=STL(JJK,ISA)
      RGD(JJK,ISA)=DDM(JJK,ISA)*X1
      X2=100.*HUI(JJK,ISA)
      AJH0=AJH1(JJK,ISA)
      IF(SCLM(3)>0.)X2=MIN(X2,SCLM(3))
      AJH1(JJK,ISA)=HI(JJK)*X2/(X2+EXP(SCRP(3,1)-SCRP(3,2)*X2))
      DHI=AJH1(JJK,ISA)-AJH0
      X3=TOPC(JJK)-TX
      RTO=RZSW(ISA)/PAW(ISA)
      IF(X3<0..AND.HUI(JJK,ISA)>.7.AND.RTO<.8)THEN    
          XX=EXP(PRMT(30)*X3/TOPC(JJK))
          DHI=DHI*XX
      END IF    
      AJHI(JJK,ISA)=AJHI(JJK,ISA)+DHI
      X3=DM(JJK,ISA)-DDM(JJK,ISA)
      X4=MAX(1.E-5,X3+RGD(JJK,ISA))
      DM(JJK,ISA)=X4
      DM1(JJK,ISA)=DM1(JJK,ISA)+RGD(JJK,ISA)
      VARC(18,JJK,ISA)=DM1(JJK,ISA)   
      SMMC(18,JJK,MO,ISA)=DM1(JJK,ISA)   
      X1=HUI(JJK,ISA)
      RF=MAX(.2,RWPC(1,JJK)*(1.-X1)+RWPC(2,JJK)*X1)
      IF(HUI(JJK,ISA)>.3)RF=RF*(1.1-.1*WS(ISA))
      RF=MIN(RF,.99)
      SMMC(20,JJK,MO,ISA)=SMMC(20,JJK,MO,ISA)+RF
      VARC(20,JJK,ISA)=RF
      IF((II==NDC(3).OR.II==NDC(6)).AND.TMN(IRF(ISA))<TBSC(JJK)+4..AND.HR1<0.)THEN
          FALLRS=MAX(0.,(TMN(IRF(ISA))-TBSC(JJK))/4.)
      ELSE
          FALLRS=0.
      END IF
      IF(SLAI(JJK,ISA)/XLAI(JJK,ISA)<.03*XLAI(JJK,ISA).AND.TGX(JJK,ISA)>0..AND.MATX(JJK,ISA)<1)THEN
          X1=TGX(JJK,ISA)/PHU(JJK,IHU(JJK,ISA),ISA)
          RW(JJK,ISA)=(PRMT(97)-X1)*RW(JJK,ISA)+(1.-PRMT(97)+X1)*RF*DM(JJK,ISA)
          IR2S=1
      ELSE
          RFDF=RW(JJK,ISA)/(DM(JJK,ISA)+1.E-5)-RF
          IF(ABS(RFDF)<.02)THEN !WE ARE VERY CLOSE TO THE DESIRED ROOT FRACTION, SO ALLOCATE NEW BIOMASS ACCORDING TO DESIRED ROOT FRACTION
			  RF2=RF
          ELSE
              IF(RFDF>0.)THEN !ROOT:BIOMASS RATIO IS WAY TOO HIGH. ALLOCATE ALL NEW BIOMASS TO SHOOTS.
			      RF2=0.
              ELSE !ROOT:BIOMASS RATIO IS WAY TOO LOW.  ALLOCATE ALL NEW BIOMASS TO ROOTS.
                  RF2=1.
              END IF
          END IF    
          RW(JJK,ISA)=RW(JJK,ISA)+RF2*RGD(JJK,ISA)
          IR2S=0
      END IF    
      RWLOS=MIN(MAX(RDEAT(JJK),0.),0.5)*RW(JJK,ISA)
      DM(JJK,ISA)=MAX(1.E-5,DM(JJK,ISA)-RWLOS)
      IF(RW(JJK,ISA)-RWLOS<RWPC(2,JJK)*DM(JJK,ISA))THEN
          RW(JJK,ISA)=RWPC(2,JJK)*DM(JJK,ISA)
          RWLOS=MAX(0.,RWL-RW(JJK,ISA))
      ELSE
          RW(JJK,ISA)=MAX(0.,RW(JJK,ISA)-RWLOS)
      END IF
      RW(JJK,ISA)=MIN(RW(JJK,ISA),DM(JJK,ISA))
      DRW=RW(JJK,ISA)-RWL
      X1=STL(JJK,ISA)
      STL(JJK,ISA)=DM(JJK,ISA)-RW(JJK,ISA)
      IF(IR2S>0)THEN
          VARC(19,JJK,ISA)=RWLOS-DRW
          SMMC(19,JJK,MO,ISA)=SMMC(19,JJK,MO,ISA)+RWLOS-DRW
      ELSE
          VARC(19,JJK,ISA)=0.
      END IF    
      DGRO(JJK,ISA)=STL(JJK,ISA)-X1
      ! CALCULATE TDN OF FORAGE
      IF(HUI(JJK,ISA)<.3.AND.MATX(JJK,ISA)<1)THEN
          TDNF(JJK,ISA)=TDNX(JJK)
      ELSE
          IF(STL(JJK,ISA)>0.)THEN
              X1=(TDNX(JJK)-TDNN(JJK))/PHU(JJK,IHU(JJK,ISA),ISA)
              TD1=STL(JJK,ISA)-DGRO(JJK,ISA)
              TD2=TDNF(JJK,ISA)-TGX(JJK,ISA)*X1
              TD3=DGRO(JJK,ISA)*TDNX(JJK)
              TDNF(JJK,ISA)=(TD1*TD2+TD3)/STL(JJK,ISA)
              TDNF(JJK,ISA)=MIN(MAX(TDNF(JJK,ISA),TDNN(JJK)),TDNX(JJK))
          ELSE
              TDNF(JJK,ISA)=TDNN(JJK)
          END IF    
      END IF    
      IF(II==NDC(7).OR.II==NDC(8).OR.II==NDC(10))THEN
          X1=100.*IDA/ND
          IF(MO>8.AND.(X1>75..OR.TMN(IRF(ISA))<0.))THEN
              IF(TREF0(JJK,ISA)<1.E-10)THEN
                  DFL(JJK,ISA)=FTO(JJK)*(STL(JJK,ISA)-STL0(JJK,ISA))
                  TREF0(JJK,ISA)=X1
                  TREFB(JJK,ISA)=100.-TREF0(JJK,ISA)
              END IF
              XX=100.*(X1-TREF0(JJK,ISA))/TREFB(JJK,ISA)
              F=XX/(XX+EXP(SCRP(30,1)-SCRP(30,2)*XX))
              FF=F-FLF0(JJK,ISA)
              W1=MAX(1.,U10(IRF(ISA))/4.)     !MAX(1.,U10(ISA)/4.)  Per Luca Doro. 6/30/2023
              FALF=MIN(DFL(JJK,ISA)*FF*W1,.1*STL(JJK,ISA))
              SMFL(JJK,ISA)=SMFL(JJK,ISA)+FALF
              IF(SMFL(JJK,ISA)>DFL(JJK,ISA))THEN
                  FALF=0.
                  SMFL(JJK,ISA)=DFL(JJK,ISA)
              END IF    
              FLF0(JJK,ISA)=F
              SMM(109,MO,ISA)=SMM(109,MO,ISA)+FALF
              DM(JJK,ISA)=DM(JJK,ISA)-FALF
              STL(JJK,ISA)=STL(JJK,ISA)-FALF
              X11=FALF*1000.
              X5=CNY(JJK)*X11
              X6=CPY(JJK)*X11
              X10=CKY(JJK)*X11
              IF(FALF>1.E-5)CALL NCNSTD(FALF,X5,LD1)
              FOP(LD1,ISA)=FOP(LD1,ISA)+X6
              UN1(JJK,ISA)=UN1(JJK,ISA)-X5
              UP1(JJK,ISA)=UP1(JJK,ISA)-X6
              UK1(JJK,ISA)=UK1(JJK,ISA)-X10
              SOLK(LD1,ISA)=SOLK(LD1,ISA)+X10
          END IF
      ELSE
          STLNEW=STL(JJK,ISA)-STLOLD(JJK,ISA)
          X7NEW=MIN(.99,.01*(HUI(JJK,ISA)+.01)**10)
          X7NEW=X7NEW*STLNEW
          X7OLD=MIN(.99,.01*(HUI(JJK,ISA)+1.)**10)
          X7OLD=X7OLD*STLOLD(JJK,ISA)
          X7=X7NEW+X7OLD
          IF(II==NDC(3).OR.II==NDC(6).OR.II==NDC(11))THEN
              IF(HUI(JJK,ISA)>1.)THEN
                  HU(JJK,ISA)=0.
                  MATX(JJK,ISA)=1
                  STLOLD(JJK,ISA)=STL(JJK,ISA)
              ELSE
                  IF(HUI(JJK,ISA)>.6.AND.HUI(JJK,ISA)<.65.AND.SLAI(JJK,ISA)/XLAI(JJK,ISA)<(BLAI(JJK,ISA)))THEN
                      HU(JJK,ISA)=0.
                      STLOLD(JJK,ISA)=STL(JJK,ISA)
                  END IF    
              END IF    
          END IF
          XX=MIN(X7,STL(JJK,ISA))
          STL(JJK,ISA)=STL(JJK,ISA)-XX
          STLOLD(JJK,ISA)=MAX(0.,STLOLD(JJK,ISA)-MAX(0.,X7OLD))
          VARC(21,JJK,ISA)=STL(JJK,ISA)-STLL
          DM(JJK,ISA)=DM(JJK,ISA)-XX
          STD(JJK,ISA)=STD(JJK,ISA)+XX
          STDL(JJK,ISA)=STDL(JJK,ISA)+CLG(JJK,ISA)*XX
          SMMC(6,JJK,MO,ISA)=SMMC(6,JJK,MO,ISA)+XX
          VARC(6,JJK,ISA)=XX
          X8=XX*BN(3,JJK)
          XUN=UN1(JJK,ISA)
          IF(XUN-X8<.01)X8=XUN-.01
          X9=XX*BP(3,JJK)
          XUP=UP1(JJK,ISA)
          IF(XUP-X9<.01)X9=XUP-.01
          STDN(JJK,ISA)=STDN(JJK,ISA)+X8
          STDP(JJK,ISA)=STDP(JJK,ISA)+X9
          UN1(JJK,ISA)=XUN-X8
          UP1(JJK,ISA)=XUP-X9
      END IF          
      SUM=0.
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          RTO=WNMU(ISL,ISA)/(WNMU(ISL,ISA)+WNO3(ISL,ISA)+1.E-5)
          UU=UN(ISL,ISA)*RTO
          UN(ISL,ISA)=UN(ISL,ISA)-UU
          IF(WNO3(ISL,ISA)<UN(ISL,ISA)) THEN
               UN(ISL,ISA)=WNO3(ISL,ISA)
          ENDIF
          IF(WNMU(ISL,ISA)<UU)UU=WNMU(ISL,ISA)
          WNO3(ISL,ISA)=WNO3(ISL,ISA)-UN(ISL,ISA)
          WNMU(ISL,ISA)=WNMU(ISL,ISA)-UU
          UTO=RWT(ISL,JJK,ISA)/RWL
          X1=RWLOS*UTO
          X2=MIN(UN1(JJK,ISA),1000.*BN(3,JJK)*X1)
          CALL NCNSTD(X1,X2,ISL)
          UN1(JJK,ISA)=UN1(JJK,ISA)-X2
          IF(DRW>0.)UTO=UW(ISL,ISA)/(AEP(JJK,ISA)+1.E-20)
          SWST(ISL,ISA)=MAX(1.E-10,SWST(ISL,ISA)-UW(ISL,ISA))
          X1=WPML(ISL,ISA)+WPMU(ISL,ISA)
          IF(X1>UP(ISL,ISA))THEN
              X2=WPML(ISL,ISA)/X1
              WPML(ISL,ISA)=WPML(ISL,ISA)-UP(ISL,ISA)*X2
              WPMU(ISL,ISA)=WPMU(ISL,ISA)-UP(ISL,ISA)*(1.-X2)
          ELSE
              UP(ISL,ISA)=X1
              WPML(ISL,ISA)=0.
              WPMU(ISL,ISA)=0.
          END IF
          X1=SOLK(ISL,ISA)+WKMU(ISL,ISA)
          IF(X1>UK(ISL,ISA))THEN
              X2=SOLK(ISL,ISA)/X1
              SOLK(ISL,ISA)=SOLK(ISL,ISA)-UK(ISL,ISA)*X2
              WKMU(ISL,ISA)=WKMU(ISL,ISA)-UK(ISL,ISA)*(1.-X2)
          ELSE
              UK(ISL,ISA)=X1
              SOLK(ISL,ISA)=0.
              WKMU(ISL,ISA)=0.
          END IF
          RWT(ISL,JJK,ISA)=RWT(ISL,JJK,ISA)+DRW*UTO
          SUM=SUM+RWT(ISL,JJK,ISA)
      END DO
      RW(JJK,ISA)=SUM
      !AD2=0.
      !AD3=0.
      !DO K=1,NBSL(ISA)
          !ISL=LID(K,ISA)
          !AD2=AD2+ST(ISL,ISA)
          !AD3=AD3+UW(ISL)
      !END DO
      !DF=AD1-AD2-AD3
      !IF(ABS(DF)>.001)WRITE(KW(1),1)IY,MO,KDA,AD1,AD2,AD3,DF
    !1 FORMAT(5X,'CAGRO',3I4,10E13.5)  
      RETURN
      END