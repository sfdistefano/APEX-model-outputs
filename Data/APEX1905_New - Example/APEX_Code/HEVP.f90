      SUBROUTINE HEVP
!     APEX1905
!     THIS SUBPROGRAM ESTIMATES DAILY EVAPOTRANSPIRATION.  THERE ARE
!     FOUR OPTIONS FOR COMPUTING POTENTIAL EVAP(PENMAN-MONTEITH, PENMAN,
!     PRIESTLEY-TAYLOR, & HARGREAVES)
      USE PARM
      SWL1=SWLT(ISA)
      AD1=SNO(ISA)+SWLT(ISA)+XRFI(ISA)
      DO K=1,NBSL(ISA)
          ISL=LID(K,ISA)
          AD1=AD1+SWST(ISL,ISA)
      END DO
      SMLA(ISA)=.01
      EPP=0.
      ADD=0.
      SEV=0.
      NN=NCP(IRO(ISA),ISA)
      CHMX(ISA)=0.
      DO K=1,NN
          K1=JE(K,ISA)
          SMLA(ISA)=SMLA(ISA)+SLAI(K1,ISA)
          IF(CPHT(K1,ISA)>CHMX(ISA))CHMX(ISA)=CPHT(K1,ISA)
          EP(K1,ISA)=0.
      END DO
      X1=MAX(.4*SMLA(ISA),PRMT(17)*(CV(ISA)+.1))
      IF(X1>10.)THEN
          EAJ=.00005
      ELSE
          EAJ=EXP(-X1)
      END IF
      IF(SNO(ISA)>5.)THEN
          ALB=.6
          EAJ=.5
      ELSE
          ALB=SALB(ISA)
          IF(SMLA(ISA)>0.)ALB=.23*(1.-EAJ)+SALB(ISA)*EAJ
      END IF
      TK=TX+273.
      RSO=RAMX
      XL=2.501-2.2E-3*TX
      EA=ASVP(TK)
      ED=EA*RHD(IRF(ISA))
      VPD=EA-ED
      SMM(9,MO,ISA)=SMM(9,MO,ISA)+VPD
      VAR(9,ISA)=VPD
      RALB1=SRAD(IRF(ISA))*(1.-ALB)
      DLT=EA*(6790.5/TK-5.029)/TK
      XX=DLT+GMA(ISA)
      TK4=TK**4
      RBO=(.34-.14*SQRT(ED))*4.9E-9*TK4
      RTO=MIN(.99,SRAD(IRF(ISA))/(RSO+.1))
      RN=RALB1-RBO*(.9*RTO+.1)
      X2=RN*DLT
      SELECT CASE(IET)
          CASE(5)
              !BAIER-ROBERTSON PET METHOD
              EO=MAX(0.,.288*TMX(IRF(ISA))-.144*TMN(IRF(ISA))+.139*RSO-4.931)
          CASE(4)
              !HARGREAVES PET METHOD
              RAMM=RSO/XL
              EO=MAX(0.,PRMT(23)*RAMM*(TX+17.8)*(TMX(IRF(ISA))-TMN(IRF(ISA)))&
              **HGX)
          CASE(3)
              !PRIESTLEY-TAYLOR PET METHOD
              RAMM=RALB1/XL
              EO=1.28*RAMM*DLT/XX
          CASE(2)
              !PENMAN PET METHOD
              FWV=2.7+1.63*U10(IRF(ISA))
              X3=GMA(ISA)*FWV*VPD
              X1=X2/XL+X3
              EO=X1/XX
          CASE(1)
              !PENMAN-MONTEITH PET METHOD
              RHO=.01276*PB/(1.+.00367*TX)
              IF(IGO(ISA)>0)THEN
                  IF(CPHT(JJK,ISA)<8.)THEN
                      UZZ=U10(IRF(ISA))
                      ZZ=10.
                  ELSE
                      ZZ=CHMX(ISA)+2.
                      UZZ=U10(IRF(ISA))*LOG(ZZ/.0005)/9.9035
                  END IF
                  X1=LOG10(CPHT(JJK,ISA)+.01)
                  Z0=10.**(.997*X1-.883)
                  ZD=10.**(.979*X1-.154)
                  RV=6.25*(LOG((ZZ-ZD)/Z0))**2/UZZ
                  X3=VPD-VPTH(JJK)
                  IF(X3>0.)THEN
                      FVPD=MAX(.1,1.-VPD2(JJK)*X3)
                  ELSE
                      FVPD=1.
                  END IF
                  G1=GSI(JJK)*FVPD
                  RC=PRMT(1)/((SMLA(ISA)+.01)*G1*EXP(.00155*(330.-CO2)))
                  EPP=(X2+86.66*RHO*VPD/RV)/(XL*(DLT+GMA(ISA)*(1.+RC/RV)))
              END IF
              RV=350./U10(IRF(ISA))
              EO=(X2+86.66*RHO*VPD/RV)/(XL*XX)
              IF(EPP>EO)EO=EPP
          CASE DEFAULT
              !HARGREAVES PET METHOD
              RAMM=RSO/XL
              EO=MAX(0.,PRMT(23)*RAMM*(TX+17.8)*(TMX(IRF(ISA))-TMN(IRF(ISA)))&
              **HGX)                
      END SELECT
      VARW(10)=EO    
      IF(LUN(ISA)==35)RETURN
      SALF=SALA(ISA)/WSA(ISA)
      EOR=EO
      IF(IET>1.OR.IET==0)EPP=MIN(SMLA(ISA)*EOR/3.,EOR)
      DO
          IF(XRFI(ISA)>EOR)THEN
              ES=EOR
              XRFI(ISA)=XRFI(ISA)-ES
              SWLT(ISA)=SWLT(ISA)+XRFI(ISA)
              EPP=0.
              EXIT
          ELSE
              EOR=EOR-XRFI(ISA)
              ADD=XRFI(ISA)
              EPP=MIN(EOR,EPP)
          END IF
          IF(IGO(ISA)>0)THEN
              XX=EPP/SMLA(ISA)
              DO K=1,NN
                  K1=JE(K,ISA)
                  EP(K1,ISA)=SLAI(K1,ISA)*XX
              END DO
          END IF
          !Actual ET in flooded/dry field
          IF(PADDY_STO(1,ISA)>EO)THEN   !Paddy ponding condition Jaehak 2016
              IF(SMLA(ISA)<=4.)THEN
                  ES=.6*(1.-SMLA(ISA)/4.)*EO !Sakaguchi et al. 2014
              ELSE
                  ES=0.
              END IF
              EO=ES+EPP
              EOR=EO
              VAR(63,ISA)=EPP
              SMM(63,MO,ISA)=SMM(63,MO,ISA)+EPP
              RETURN
          ELSE !unsaturated soil
              ES=EOR*EAJ
              ST0(ISA)=RALB1
              ES=MIN(ES,ES*EOR/(ES+EPP+1.E-10))
              SNO1=SNO(ISA)
              IF(SNO(ISA)>=ES)THEN
                  NEV=1
                  SNO(ISA)=SNO(ISA)-ES
              ELSE
                  XX=ES-SNO(ISA)
                  ES=SNO(ISA)
                  IF(SWLT(ISA)<XX)THEN
                      XX=XX-SWLT(ISA)
                      ES=ES+SWLT(ISA)
                      SWLT(ISA)=0.
                      TOT=0.
                      DO J=1,NBSL(ISA)
                          ISL=LID(J,ISA)
                          RTO=1000.*Z(ISL,ISA)
                          IF(SCLM(2)>0.)RTO=MIN(RTO,SCLM(2))
                          SUM=XX*RTO/(RTO+EXP(SCRP(2,1)-SCRP(2,2)*RTO))
                          XZ=FC(ISL,ISA)-S15(ISL,ISA)
                          IF(SWST(ISL,ISA)<FC(ISL,ISA))THEN
                              F=EXP(PRMT(12)*(SWST(ISL,ISA)-FC(ISL,ISA))/XZ)
                          ELSE
                              F=1.
                          END IF
                          ZZ=SUM-TOT
                          XY=PRMT(5)*S15(ISL,ISA)
                          SEV(ISL,ISA)=MAX(0.,MIN(ZZ*F,SWST(ISL,ISA)-XY))
                          IF(Z(ISL,ISA)>.2)EXIT       !0.2 M FOR SOIL EVAPORATION,HARD CODED, LIWANG MA
                          ES=ES+SEV(ISL,ISA)
                          TOT=SUM
                      END DO
                      IF(J<NBSL(ISA))THEN
                          Z1=Z(LID(J-1,ISA),ISA)
                          RTO=(.2-Z1)/(Z(ISL,ISA)-Z1)
                          X1=RTO*SWST(ISL,ISA)
                          X2=RTO*XY
                          SEV(ISL,ISA)=MAX(0.,MIN(SEV(ISL,ISA),X1-X2))
                          ES=ES+SEV(ISL,ISA)
                      ELSE
                          J=NBSL(ISA)
                      END IF
                      NEV=J
                  ELSE
                      SWLT(ISA)=SWLT(ISA)-XX
                      ES=ES+XX
                      NEV=1
                  END IF
              END IF 
          END IF    
          SUM=SNO1-SNO(ISA)+SWL1-SWLT(ISA)
          DO J=1,NBSL(ISA)
              SEV(J,ISA)=SEV(J,ISA)*SALF
              SUM=SUM+SEV(J,ISA)
              SWST(J,ISA)=SWST(J,ISA)-SEV(J,ISA)
          END DO
          ES=SUM
          EXIT
      END DO    
      XX=EOR-ES
      IF(EPP>XX)THEN
          X1=XX/EPP
          EPP=XX
          DO K=1,NN
              K1=JE(K,ISA)
              EP(K1,ISA)=EP(K1,ISA)*X1
          END DO 
      END IF
      VAR(63,ISA)=EPP
      SMM(63,MO,ISA)=SMM(63,MO,ISA)+EPP
      ES=ES+ADD
      AD2=SNO(ISA)+SWLT(ISA)
      DO K=1,NBSL(ISA)
          ISL=LID(K,ISA)
          AD2=AD2+SWST(ISL,ISA)
      END DO
      DF=AD1-ES-AD2
      IF(ABS(DF)>.001)WRITE(KW(1),1)IY,MO,KDA,AD1,ES,AD2,DF
    1 FORMAT(5X,'HEVP',3I4,10E13.5)  
      RETURN
      END