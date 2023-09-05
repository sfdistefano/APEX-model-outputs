      SUBROUTINE TLOP(CSTX,COX,JRT)
!     APEX1905
!     THIS SUBPROGRAM CONTROLS ALL TILLAGE OPERATIONS INCLUDING PLANTING
!     & HARVESTING.
      USE PARM
!      DIMENSION YLSD(MNC,MSA)
      FNPP(X)=DMLA(JJK)*X/(X+EXP(PPCF(1,JJK)-PPCF(2,JJK)*X))
      JRT=0
      II=IHC(JT1)
      LD1=LID(1,ISA)
      NN=NCP(IRO(ISA),ISA)
      N1=MAX(1,NN)
      WSAX=WSA(ISA)
      X1=CND(IRO(ISA),KT(ISA),ISA)
      IF(ABS(X1-CN0(ISA))>0.)THEN
          X2=SMX(ISA)
          CALL HCNSLP(X1,X3)
          CN0(ISA)=X1
          CN2(ISA)=X1
          SCI(ISA)=SMX(ISA)*SCI(ISA)/X2
      END IF
      SELECT CASE(II)
          CASE(1,2,3,19,22)
              GO TO 10
          CASE(5)
              IDRL(ISA)=0
          CASE(6)
              IDRL(ISA)=1
          CASE(7,8)
              GO TO 57
          CASE(10)
              CSTX=-CSTX*YLD1(JJK,ISA)/(1.-WCY(JJK))
              COX=CSTX
              GO TO 57
          CASE(11)
              CSTX=-CSTX*YLD(JJK,ISA)/(1.-WCY(JJK))
              COX=CSTX
              GO TO 57
          CASE(12,13)
              IF(ICUS(JT1)==0)GO TO 57
              CSTX=-CSTX*YLD(JJK,ISA)/(1.-WCY(JJK))
              COX=CSTX
              GO TO 57
          CASE(14)
              CALL TBURN
              GO TO 7
          CASE(15) !Puddle paddy fields Jaehak Jeong 2016 paddy
	          SATC(LID(2,ISA),ISA)=PRMT(39)
              GO TO 57
          CASE(16) !Destroy Puddle Jaehak Jeong 2016 paddy
              SATC(LID(2,ISA),ISA)=SATK(ISA)
              GO TO 57
          CASE(21)
              KOMP(KT(ISA),ISA)=0
              ISCP(ISA)=ISCP(ISA)+1
              IF(ISCP(ISA)<MSCP)THEN
                  JRT=1
                  RETURN
              END IF
              ISCP(ISA)=0
              XX=1.-ORHI(JT1)
              X4=RSDM(LD1,ISA)*ORHI(JT1)
              RSDM(LD1,ISA)=RSDM(LD1,ISA)-X4
              X1=WCMU(LD1,ISA)*ORHI(JT1)
              WCMU(LD1,ISA)=WCMU(LD1,ISA)-X1
              X2=WCOU(LD1,ISA)*ORHI(JT1)
              WCOU(LD1,ISA)=WCOU(LD1,ISA)-X2
              X6=WNOU(LD1,ISA)*ORHI(JT1)
              WNOU(LD1,ISA)=WNOU(LD1,ISA)-X6
              X7=WNMU(LD1,ISA)*ORHI(JT1)
              WNMU(LD1,ISA)=WNMU(LD1,ISA)-X7
              WLM(LD1,ISA)=WLM(LD1,ISA)*XX
              WLS(LD1,ISA)=WLS(LD1,ISA)*XX
	          WLSL(LD1,ISA)=WLSL(LD1,ISA)*XX
              X2=WLMC(LD1,ISA)*ORHI(JT1)
              WLMC(LD1,ISA)=WLMC(LD1,ISA)-X2
              X1=WLSC(LD1,ISA)*ORHI(JT1)
              WLSC(LD1,ISA)=WLSC(LD1,ISA)-X1
              X101=WSAX*(X1+X2)
              SMM(101,MO,ISA)=SMM(101,MO,ISA)+X101
              VAR(101,ISA)=X101
              WLSLC(LD1,ISA)=WLSLC(LD1,ISA)*XX
              WLSLNC(LD1,ISA)=WLSLNC(LD1,ISA)*XX
              X2=ORHI(JT1)*WNO3(LD1,ISA)
              X3=ORHI(JT1)*WNH3(LD1,ISA)
              WNO3(LD1,ISA)=MAX(1.E-5,WNO3(LD1,ISA)-X2)
              WNH3(LD1,ISA)=WNH3(LD1,ISA)-X3
              X1=ORHI(JT1)*WLMN(LD1,ISA)
              X5=ORHI(JT1)*WLSN(LD1,ISA)
              WLMN(LD1,ISA)=WLMN(LD1,ISA)-X1
              WLSN(LD1,ISA)=WLSN(LD1,ISA)-X5
              X9=WSAX*(X1+X2+X3+X5+X6+X7)
              SMM(89,MO,ISA)=SMM(89,MO,ISA)+X9
              VAR(89,ISA)=X9
              X3=WPOU(LD1,ISA)*ORHI(JT1)
              WPOU(LD1,ISA)=WPOU(LD1,ISA)-X3
              X1=ORHI(JT1)*FOP(LD1,ISA)
              FOP(LD1,ISA)=FOP(LD1,ISA)-X1
              X2=ORHI(JT1)*WPML(LD1,ISA)
              WPML(LD1,ISA)=WPML(LD1,ISA)-X2
              X5=WPMU(LD1,ISA)*ORHI(JT1)
              WPMU(LD1,ISA)=WPMU(LD1,ISA)-X5
              X90=WSAX*(X1+X2+X3+X5)
              SMM(90,MO,ISA)=SMM(90,MO,ISA)+X90
              VAR(90,ISA)=X90
              X1=ORHI(JT1)*SOLK(LD1,ISA)
              SOLK(LD1,ISA)=SOLK(LD1,ISA)-X1
              X2=ORHI(JT1)*EXCK(LD1,I)
              EXCK(LD1,I)=EXCK(LD1,I)-X2
              X3=ORHI(JT1)*FIXK(LD1,I)
              FIXK(LD1,I)=FIXK(LD1,I)-X3
              SMM(153,MO,ISA)=SMM(153,MO,ISA)+X1+X2+X3
              VAR(153,ISA)=X1+X2+X3
              SMNU(IDON(ISA))=SMNU(IDON(ISA))+X4*WSAX
              IF((NOP>0.OR.NBSA(ISA)==ISAP).AND.KFL(1)>0)WRITE(KW(1),44)ISA,&
              NBSA(ISA),IYR,MO,KDA,TIL(JT1),X4,ORHI(JT1),RSD(LD1,ISA),&
              RSDM(LD1,ISA),SMNU(IDON(ISA)),XHSM(ISA)
              RETURN
          CASE(23)
              ICV=1
              GO TO 57
          CASE(24)
              ICV=0
              GO TO 57
          CASE(27)
              CALL NMULCH(KT(ISA))
              RETURN
          CASE DEFAULT
              GO TO 6
      END SELECT
      ISL=LID(2,ISA)
      DO K=1,LC
          !I2=LY(IRO(ISA),K,ISA)
          IF(JH(IRO(ISA),KT(ISA),ISA)==KDC(K))EXIT
      END DO
      IF(K>LC.OR.KGO(K,ISA)>0)GO TO 26
      ZX=0.
      DO I=1,NBSL(ISA)
          ISL=LID(I,ISA)
          Z1=Z(ISL,ISA)
          ZZ=.5*(ZX+Z1)
          IF(ZZ>=.075)EXIT
   	      ZX=Z1 
      END DO
      IF(I>NBSL(ISA))THEN
          ISL=LID(NBSL(ISA),ISA)
          LRD(ISA)=NBSL(ISA)
      ELSE
          LRD(ISA)=I
      END IF
      IF(STMP(ISL,ISA)<TBSC(K)+2.)THEN
          KOMP(KT(ISA),ISA)=0
          JRT=1
          RETURN 
      END IF
      AWC(JJK,ISA)=RZSW(ISA)
      IGO(ISA)=IGO(ISA)+1
!	  TCPA(I2)=TCPA(I2)+WSA(ISA)
	  KC(ISA)=1
      DO J=1,NN
          IF(JE(J,ISA)>=MNC)EXIT
      END DO
      KC(ISA)=MIN(J,NN)
      JE(KC(ISA),ISA)=K
      JJK=K
      KGO(JJK,ISA)=K
      JP(JJK,ISA)=0
	  IYH(JJK,ISA)=1
      SWH(JJK,ISA)=0.
      SWP(JJK,ISA)=0.
      ACET(JJK,ISA)=0.
      XDLAI(JJK)=DLAI(JJK)
      IF(PADDY_STO(1,ISA)>=.1.AND.PADDY_HWEIR(ISA)>0)THEN
          SLAI(JJK,ISA)=LAI_INIT
          BIR(ISA)=0.
      END IF    
	  IF(CPNM(JJK)=='FALW')NCR(JJK,ISA)=NCR(JJK,ISA)+1
      WCYD=.3
      RD(JJK,ISA)=TLD(JT1)
      HU(JJK,ISA)=0.
      DM(JJK,ISA)=SDW(JJK)*5.E-4
      DM1(JJK,ISA)=DM(JJK,ISA)
      RW(JJK,ISA)=.4*DM(JJK,ISA)
      RWT(ISL,JJK,ISA)=RW(JJK,ISA)
      ROSP(ISA)=RIN(JT1)
      PPL0(JJK,ISA)=POP(JJK,IHU(JJK,ISA),ISA)
      XLAI(JJK,ISA)=FNPP(PPL0(JJK,ISA))
      DMLX(JJK,ISA)=XLAI(JJK,ISA)
	  X1=SDW(JJK)*CSTS(JJK)
      COST(ISA)=COST(ISA)+X1
      LRD(ISA)=MAX(2,LRD(ISA))
      IPLD(JJK,ISA)=IYR*10000+MO*100+KDA
	  IHVD(JJK,ISA)=0
      IMTU(JJK,ISA)=0
      JPL(JJK,ISA)=1
      IF((NOP>0.OR.NBSA(ISA)==ISAP).AND.KFL(1)>0)WRITE(KW(1),32)ISA,&
      NBSA(ISA),IYR,MO,KDA,CPNM(JJK),CV(ISA)
      IF(KFL(31)>0)WRITE(KW(31),49)ISA,NBSA(ISA),IYR,MO,KDA,TIL(JT1),&
      KDC(JJK),II,NBE(JT1),NBT(JT1),X1,X1,SDW(JJK)
    6 EE=EMX(JT1)
      IF(II==NHC(19))THEN
          DO KHD=1,NHRD(IOW)
              IHD=KOW(KHD,IOW)
              IF(IGZX(IHD,IOW)==ISA)EXIT
          END DO
          IF(KHD<=NHRD(IOW))EE=EE*GZNB(IHD,ISA)/WSAX
      END IF
      STRC(ISA)=STRC(ISA)+STIR(JT1)   !Liwang Ma, Stir factor (e.g. 10), make STRC infinity 10^10, huge number
      PPL0(JJK,ISA)=(1.-FPOP(JT1))*PPL0(JJK,ISA)
      XLAI(JJK,ISA)=FNPP(PPL0(JJK,ISA))
      DMLX(JJK,ISA)=XLAI(JJK,ISA)
      DMX=TLD(JT1)
      !Mix soil properties in flooded/dry paddy fields Jaehak Jeong 2016
      IF(PADDY_STO(1,ISA)>=.1)THEN
          CALL TMIX_FLOODED(EE,DMX,0,0) 
      ELSE
          CALL TMIX(EE,DMX,0,0)
      ENDIF     
	  IF(DMX>BIG(ISA))TLD(JT1)=BIG(ISA)
      IF(II==NHC(22).OR.II==NHC(19))GO TO 26
	  IF(II==NHC(15))THEN
	      SATC(LID(2,ISA),ISA)=PRMT(39)
	  ELSE
	      IF(II==NHC(16))SATC(LID(2,ISA),ISA)=SATK(ISA)
      END IF
   57 IF(IDR(ISA)>0)THEN
          IF(II==NHC(25))THEN
              HCL(IDR(ISA),ISA)=HCLN(ISA)
          ELSE
              IF(II==NHC(26))HCL(IDR(ISA),ISA)=HCLD(ISA)
          END IF
      END IF
      IF(KFL(31)>0)WRITE(KW(31),50)ISA,NBSA(ISA),IYR,MO,KDA,TIL(JT1),&
      KDC(JJK),II,NBE(JT1),NBT(JT1),CSTX,COX,FULU(JT1)
      SMFU(ISA)=SMFU(ISA)+FULU(JT1)
      SMST(ISA)=SMST(ISA)+STIR(JT1)
    7 XX=TLD(JT1)*1000.
      IF((NOP>0.OR.NBSA(ISA)==ISAP).AND.KFL(1)>0)WRITE(KW(1),28)ISA,&
      NBSA(ISA),IYR,MO,KDA,TIL(JT1),XX,XHSM(ISA)
      IF(II/=NHC(17).AND.II/=NHC(18))GO TO 26
      IF(II/=NHC(18))THEN
          DHT(ISA)=DKH(JT1)
          DKHL(ISA)=DHT(ISA)
          DKIN(ISA)=DKI(JT1)
          IF((NOP>0.OR.NBSA(ISA)==ISAP).AND.KFL(1)>0)WRITE(KW(1),30)ISA,&
          NBSA(ISA),IYR,MO,KDA,DHT(ISA),DKIN(ISA),XHSM(ISA)
          GO TO 26
      END IF
      DHT(ISA)=0.
      DKHL(ISA)=0.
      IF((NOP>0.OR.NBSA(ISA)==ISAP).AND.KFL(1)>0)WRITE(KW(1),30)ISA,&
      NBSA(ISA),IYR,MO,KDA,DHT(ISA),DKIN(ISA),XHSM(ISA)
      GO TO 26
   10 CALL PESTF
      X6=PSTF(ISA)
      NPSF(ISA)=NPSF(ISA)+1
      IF(II==NHC(19))THEN
          KOMP(KT(ISA),ISA)=0
          CALL TGRAZ
          GO TO 6
      END IF
      IF(ORHI(JT1)>0..OR.II==NHC(22))THEN
          CALL THVOR(X6,JRT)
             IF(IHD>0)THEN
                 IF(HAY_supply(IHD,IY,IOW)>0)THEN   !it works before giving birth or after weaning
                  Hay_amount(IHD,IY,IOW)=Hay_amount(IHD,IY,IOW)+Cut_hay(ISA)*1000*WSA(ISA)*0.75   !!!calculate cut hay for total hay application
                  Cut_hay(ISA)=0.
                 ENDIF
              ENDIF
	      IF(JRT==0)GO TO 6
	      RETURN
	  END IF
      DO K=1,LC
          !J1=JE(K,ISA)
          IF(KGO(K,ISA)==0)CYCLE
          IF(JH(IRO(ISA),KT(ISA),ISA)==KDC(K))EXIT
      END DO
      IF(K>LC)THEN
           CSTX=0.
           COX=0.
           JRT=1
           RETURN 
      ELSE
          JJK=KDC1(KDC(K))
          DMF(JJK,ISA)=DM1(JJK,ISA)
          IF(JP(JJK,ISA)==0)THEN
              JP(JJK,ISA)=1
              NCR(JJK,ISA)=NCR(JJK,ISA)+1
          END IF
          IF(II/=NHC(1))THEN
              HUF(JJK,ISA)=MAX(HUF(JJK,ISA),HU(JJK,ISA))
              TRA(JJK,ISA)=SRA(JJK,ISA)+TRA(JJK,ISA)
              IF(RD(JJK,ISA)>RDF(JJK,ISA))RDF(JJK,ISA)=RD(JJK,ISA)
              XX=DM(JJK,ISA)+.001
              X2=UN1(JJK,ISA)/XX
              X7=X2*.001
              X3=UP1(JJK,ISA)/XX
              XX=STD(JJK,ISA)+1.E-10
              RNR=STDN(JJK,ISA)/XX
              RPR=STDP(JJK,ISA)/XX
              RKR=STDK(JJK,ISA)/XX
              STDL(JJK,ISA)=CLG(JJK,ISA)*XX
              RLR=STDL(JJK,ISA)/XX
              XX=100.*SWH(JJK,ISA)/(SWP(JJK,ISA)+1.E-5)
              IF(SCLM(10)>0.)XX=MIN(XX,SCLM(10))
              FWS=XX/(XX+EXP(SCRP(10,1)-SCRP(10,2)*XX))
              !FNS=1.-SFCP(2,JJK,ISA)/REAL(NGD(JJK,ISA))
              XX=MAX(AJHI(JJK,ISA)-WSYF(JJK),0.)
              FT=MAX(.1,1.+PRMT(81)*(IYR-2000))
              X1=MIN(FWS*XX+WSYF(JJK),.9*DM(JJK,ISA)/(STL(JJK,ISA)+1.E-10))*FT
              X1=MAX(X1,WSYF(JJK))
              X2=1000.*CNY(JJK)*(X7/BN(3,JJK))**.1
              X3=1000.*CPY(JJK)*(.001*X3/BP(3,JJK))**.1
              X8=1000.*CKY(JJK)
              XZ=X1*STL(JJK,ISA)
              YZ=X1*STD(JJK,ISA)
              ZZ=MAX(.01,1.-X1)
              CPHT(JJK,ISA)=CPHT(JJK,ISA)*ZZ
              HU(JJK,ISA)=MAX(.1,HU(JJK,ISA)*ZZ)
              SLAI(JJK,ISA)=SLAI(JJK,ISA)*ZZ
              STL(JJK,ISA)=STL(JJK,ISA)*ZZ
              STD(JJK,ISA)=STD(JJK,ISA)*ZZ
              STDL(JJK,ISA)=STDL(JJK,ISA)*ZZ
              YLD(JJK,ISA)=XZ*HE(JT1)*X6
              YLSD(JJK,ISA)=YZ*HE(JT1)
              Y4=YZ*RNR
              Y5=YZ*RPR
              Y6=YZ*RKR
              STDN(JJK,ISA)=MAX(0.,STDN(JJK,ISA)-Y4)
              STDP(JJK,ISA)=MAX(0.,STDP(JJK,ISA)-Y5)
              STDK(JJK,ISA)=MAX(0.,STDK(JJK,ISA)-Y6)
              STDL(JJK,ISA)=MAX(STDL(JJK,ISA)-YZ*RLR,.1*STD(JJK,ISA))
              X4=MIN(XZ*X2,UN1(JJK,ISA))
              X5=MIN(XZ*X3,UP1(JJK,ISA))
              X11=XZ-YLD(JJK,ISA)+YZ-YLSD(JJK,ISA)
              Z2=YLSD(JJK,ISA)*RNR
              Z3=YLSD(JJK,ISA)*RPR
              Z4=YLSD(JJK,ISA)*RKR
              YLN(JJK,ISA)=MIN(.9*(UN1(JJK,ISA)+STDN(JJK,ISA)),YLD(JJK,ISA)*X2+Z2)
              YLP(JJK,ISA)=MIN(.9*(UP1(JJK,ISA)+STDP(JJK,ISA)),YLD(JJK,ISA)*X3+Z3)
              YLK(JJK,ISA)=MIN(.9*(UK1(JJK,ISA)+STDK(JJK,ISA)),YLD(JJK,ISA)*X8+Z4)
              X10=X4-YLN(JJK,ISA)+Y4
              CALL NCNSTD(X11,X10,LD1)
              FOP(LD1,ISA)=MAX(.01,FOP(LD1,ISA)+X5-YLP(JJK,ISA)+Y5)
              YY=YLD(JJK,ISA)+YLSD(JJK,ISA)
              YLD(JJK,ISA)=YY
	          IF(IDC(JJK)/=NDC(9))THEN
                  YLD1(JJK,ISA)=YLD1(JJK,ISA)+YY
              ELSE
                  YLD1(JJK,ISA)=YLD1(JJK,ISA)+FTO(JJK)*YY
                  YLD2(JJK,ISA)=YLD2(JJK,ISA)+YLD1(JJK,ISA)*(1./FLT(JJK)-1.)
              END IF
              X2=DM(JJK,ISA)
              X3=RW(JJK,ISA)
              JD(ISA)=JJK
              SRA(JJK,ISA)=0.
              UN1(JJK,ISA)=UN1(JJK,ISA)-X4
              UP1(JJK,ISA)=UP1(JJK,ISA)-X5
              DM(JJK,ISA)=DM(JJK,ISA)-XZ
              YLNF(JJK,ISA)=YLNF(JJK,ISA)+YLN(JJK,ISA)
              YLPF(JJK,ISA)=YLPF(JJK,ISA)+YLP(JJK,ISA)
              YLKF(JJK,ISA)=YLKF(JJK,ISA)+YLK(JJK,ISA)
              TYN(ISA)=TYN(ISA)+YLN(JJK,ISA)
              TYP(ISA)=TYP(ISA)+YLP(JJK,ISA)
              TYK(ISA)=TYK(ISA)+YLK(JJK,ISA)
              HIF(JJK,ISA)=X1
              IHVD(JJK,ISA)=IYR*10000+MO*100+KDA
	          IF(ICUS(JT1)/=0.AND.CSTX<1.E-10)THEN
                  CSTX=-CSTX*YLD(JJK,ISA)
                  COX=CSTX
              END IF
              IF((NOP>0.OR.NBSA(ISA)==ISAP).AND.KFL(1)>0)WRITE(KW(1),29)ISA,&
              NBSA(ISA),IYR,MO,KDA,TIL(JT1),CPNM(JD(ISA)),YY,X2,X3,X1,X6,X7,XHSM(ISA),AJHI(JJK,ISA)
          ELSE
              IF(IPD==5)THEN
                  CALL SPRNT
                  IF(KFL(1)>0)THEN
                      WRITE(KW(1),31)
                      CALL SOLIOP
                      CALL SOLIOC
                  END IF    
              END IF
              IF(IDC(JJK)==NDC(7).OR.IDC(JJK)==NDC(8).OR.IDC(JJK)==NDC(10))THEN
                  I1=LYR(IRO(ISA),KT(ISA),ISA)
                  IF(IYH(JJK,ISA)/=I1.AND.I1/=1)THEN
                      KOMP(KT(ISA),ISA)=0
                      JRT=1
                      RETURN
                  END IF
              END IF    
              CALL TRDST
              TPSF(JJK,ISA)=TPSF(JJK,ISA)+X6
          END IF
          GO TO 6
      END IF    
   26 RETURN
   28 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,2X,'DPTH = ',F7.0,' mm',2X,'HUSC = ',F6.2)
   29 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,2X,A4,2X,'YLD=',F6.2,'t/ha',2X,&
      'BIOM=',F6.2,'t/ha',2X,'RW=',F5.2,'t/ha',2X,'HI=',F6.2,2X,'PSTF=',&
      F5.2,2X,'NCN=',F6.3,'G/G',2X,'HUSC=',F5.2,2X,'AJHI=',F6.3)
   30 FORMAT(1X,2I8,1X,I4,2I2,2X,'DKH = ',F6.0,' mm',3X,'DKI = ',F7.2,&
      ' m',2X,'HUSC= ',F5.2)
   31 FORMAT(T5,'SOIL DATA')
   32 FORMAT(1X,2I8,1X,I4,2I2,2X,A4,2X,'RSD = ',F5.1,'t')
   44 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,2X,'MNU SCRP=',F8.2,'t',2X,'SCRP EF=&
      ',F6.3,2X,'RSD REMAIN=',E13.5,'t/ha',2X,'MNU REMAIN=',E13.5,'t/ha'&
      ,2X,'MNU STK PL=',E13.5,'t',2X,'HUSC=',F5.2)
   49 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,8X,I6,6X,3I4,F10.2,10X,3F10.2)
   50 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,8X,I6,6X,3I4,2F10.2,20X,F10.2)
      END