      SUBROUTINE RESRT
!     APEX1905
!     THIS SUBPROGRAM ROUTES WATER AND SEDIMENT THROUGH RESEVOIRS.
!     COMPUTES EVAPORATION AND SEEPAGE FROM THE RESERVOIR.
      USE PARM
      DIMENSION VBK(5)
      IDO=IDOT(ICMD)
      IDN1=IDN1T(ICMD)
      IDN2=IDN2T(ICMD)
      IDRO(IDN2)=IDO
	  IDX=IDOA(IDN2)
	  IDNB(IDO)=NBSA(IDN2)
	  LD1=LID(1,IDN2)
	  RWSA(IDO)=RWSA(IDN1)
      !IF(WYLD(IDN1)<1.E-10)RETURN
      YB=RSYS(IDN2)
      RSSF(IDO)=0.
	  SST(IDO)=0.
      QDR(IDO)=0.
      QRF(IDO)=0.
      CPVH(IDO)=0.
	  DEP=0.
	  Y2=0.
      Y2O=0.
      RSPK=0.
      VBK=0.
      Q1M3=WYLD(IDN1)
      IF(IEXT(IDN2)>0)THEN
	      A1=WSA(IDN2)
          II=IDN2
	  ELSE
	      A1=RWSA(IDN1)
          II=IDN2-1
      END IF
      II0=II
      A10=10.*A1
      V0=RSV(IDN2)
      TC(IDO)=TC(IDN1)
      SST(IDO)=SST(IDN1)
	  !XX=MAX(.1,V0+RVP0(IDN2)-RSVP(IDN2)+RVE0(IDN2)-RSVE(IDN2))
      !X2=BR1(IDN2)*XX**BR2(IDN2)
   	  !RSSA(IDN2)=MIN(X2,RSAE(IDN2))
      STV(15,MO,IDN2)=RSSA(IDN2)
      VARS(15)=RSSA(IDN2)
      X1=RSSA(IDN2)
      EVRS(IDN2)=10.*EO*X1
      SP=RSHC(IDN2)*X1*240.
      PRCP=RFV(IRF(IDN2))-VAR(6,IDN2)
      RFRA=PRCP*X1*10.
	  IF(IDX==IDN1)THEN
          QIX=0.
      ELSE
          QIX=Q1M3
      END IF    
	  QSA=WYLD(IDX)
      VARH(37,IDX)=QSA
	  SMM(125,MO,IDN2)=SMM(125,MO,IDN2)+QIX
	  SMM(126,MO,IDN2)=SMM(126,MO,IDN2)+QSA
	  SMM(127,MO,IDN2)=SMM(127,MO,IDN2)+RFRA
	  Q1=QIX+QSA+RFRA
	  Y1=YSD(NDRV,IDN1)
	  YU1=YMNU(IDN1)
	  RSM1=RSDM(LD1,IDN2)
	  XX=SP+EVRS(IDN2)
      X1=Q1/A10
      SMM(64,MO,IDN2)=SMM(64,MO,IDN2)+Q1
      VAR(64,IDN2)=Q1
      V0=V0+Q1
      IF(V0<=XX)THEN
          X1=V0/(XX+1.E-10)
          EVRS(IDN2)=EVRS(IDN2)*X1
          SP=SP*X1
          RSPK(IDN2)=.1*SP/WSA(IDN2)
          RSV(IDN2)=0.
          QVOL(IDO)=0.
          YSD(NDRV,IDO)=0.
          YN(IDO)=0.
          YP(IDO)=0.
          QN(IDO)=0.
          QP(IDO)=0.
          TSFN(IDO)=0.
          RSFN(IDO)=0.
	      RSYS(IDN2)=0.
	      RSOP(IDN2)=.0001
          RSON(IDN2)=.0001
	      RSSP(IDN2)=.0001
	      RSO3(IDN2)=.0001
	      DEP=YB
	  ELSE
          V0=V0-XX
          VRR=V0-RSVE(IDN2)
          OFLO=0.
	      OFP=0.
	      I1=1
	      A3=A1
          IF(VRR>0.)OFLO=VRR
          VVR=V0-RSVP(IDN2)
          IF(VVR>0.)THEN
	          X1=MIN(RSRR(IDN2),RSVE(IDN2)-RSVP(IDN2),VVR)
	          IF(ISAO(IDN2)==0)THEN
	              OFLO=OFLO+X1
	          ELSE
	              OFP=X1
	              I2=NISA(ISAO(IDN2))
	              I1=IDOA(I2)
	              A3=WSA(I2)
                  QRP(I1)=QRP(I1)+OFP
 	              SMM(111,MO,IDN2)=SMM(111,MO,IDN2)+OFP
	              VAR(111,IDN2)=VAR(111,IDN2)+OFP
	              V0=V0-OFP
	          END IF    
              V0=V0-OFLO
          END IF
          RSV(IDN2)=V0
          STV(13,MO,IDN2)=V0
          VARS(13)=V0
          !AET=(10.*SALA(ISA)*AET+EVRS(IDN2))/(10.*WSA(IDN2))
          !SMM(11,MO,ISA)=SMM(11,MO,ISA)+AET
          !VAR(11,ISA)=AET
          XX=MAX(.1,V0+RVP0(IDN2)-RSVP(IDN2)+RVE0(IDN2)-RSVE(IDN2))
          X2=BR1(IDN2)*XX**BR2(IDN2)
          RSSA(IDN2)=MIN(X2,RSAE(IDN2))
          SALA(IDN2)=MAX(0.,WSA(IDN2)-RSSA(IDN2))
          VV=MAX(1.E-5,V0+OFLO+OFP)
          SMM(65,MO,IDN2)=SMM(65,MO,IDN2)+OFLO
          VAR(65,IDN2)=OFLO
          WYLD(IDO)=OFLO
          VARH(34,IDO)=OFLO
          SMMH(34,MO,IDO)=SMMH(34,MO,IDO)+OFLO
          VARH(35,IDO)=OFLO/86400.
          SMMH(35,MO,IDO)=SMMH(35,MO,IDO)+VARH(35,IDO)
          QVOL(IDO)=OFLO
	      IF(IEXT(IDN2)==0.AND.V0>RSVP(IDN2))THEN
              !COMPUTE BACKWATER VOLUME
              SUM=WSA(IDN2)
              SUM1=SUM
              X1=24.*RSHC(IDN2)
              IF(RSSA(IDN2)<SUM)THEN
                  RTO=RSSA(IDN2)/WSA(IDN2)
                  RSPK(IDN2)=X1*RTO
                  IQT=1
              ELSE
                  RSPK(IDN2)=X1
                  IQT=0
              END IF    
              !DISTRIBUTE BACKWATER TO UPSTREAM SA'S
              DO
                  IF(IQT>0)EXIT
                  SUM=SUM+WSA(II)
                  IF(RSSA(IDN2)<SUM)THEN
                      RTO=(RSSA(IDN2)-SUM1)/WSA(II)
                      RSPK(II)=X1*RTO
                      IQT=1
                  ELSE
                      RSPK(II)=X1
                  END IF
                  VBK(II)=10.*RSPK(II)*WSA(II)
                  SMM(86,MO,II)=SMM(86,MO,II)+VBK(II)
                  VAR(86,II)=VBK(II)
                  II=II-1
                  SUM1=SUM
              END DO 
          ELSE
              RSPK(IDN2)=240.*RSHC(IDN2)*RSSA(IDN2)
              SMM(86,MO,IDN2)=SMM(86,MO,IDN2)+RSPK(IDN2)
          END IF
          YY=RSYS(IDN2)+Y1
	      IF(YY>0.)THEN
              SMM(68,MO,IDN2)=SMM(68,MO,IDN2)+Y1
              VAR(68,IDN2)=Y1
              CY=YY/VV
              CD=MAX(RSYN(IDN2),(CY-RSYN(IDN2))*RSDP(IDN2)+RSYN(IDN2))
	          X1=MIN(CD,CY)
              Y2=OFLO*X1
	          Y2O=MIN(YY,OFP*X1)
	          YSD(NDRV,I1)=YSD(NDRV,I1)+Y2O    !/A3    area division removed to have YSD reported as mass - 2020-03-02
	          SMM(112,MO,IDN2)=SMM(112,MO,IDN2)+Y2O
	          VAR(112,IDN2)=Y2O
	          YSD(NDRV,IDO)=Y2    !/A1    area division removed to have YSD reported as mass - 2020-03-02
	          RSM1=WSA(IDN2)*RSM1
	          YYU=RSM1+YU1
	          CUN=YYU/VV
	          CDU=CD*CUN/CY
	          X1=MIN(CDU,CUN)
	          YU2=MIN(OFLO*X1,.5*YU1+.1*RSM1)
	          YMNU(IDO)=MAX(0.,YU2)
	          RSDM(LD1,IDN2)=(YYU-YU2)/WSA(IDN2)
	          SMM(69,MO,IDN2)=SMM(69,MO,IDN2)+Y2
              VAR(69,IDN2)=Y2
              DEP=MAX(0.,V0*(CY-CD))
              SRCH(13,IDO)=SRCH(13,IDO)+DEP
              SMM(70,MO,IDN2)=SMM(70,MO,IDN2)+DEP
              VAR(70,IDN2)=DEP
              X1=YY-Y2-Y2O-DEP
              RSYS(IDN2)=MAX(1.E-10,X1)
              XX=X1/A1
              STV(14,MO,IDN2)=X1
              VARS(14)=XX
              TOT=0.
              DRTO=MIN(1.,(.001+Y2)/(Y1+.001))
              B1=LOG(DRTO)/4.47
              DO I=1,NSZ
                  PCTH(I,IDO)=PCT(I,IDN1)
                  X1=MAX(-10.,B1*PSZY(I))
                  PCT(I,IDO)=PCT(I,IDN1)*EXP(X1)
                  TOT=TOT+PCT(I,IDO)
              END DO
              DO I=1,NSZ
                  PCT(I,IDO)=PCT(I,IDO)/(TOT+1.E-10)
              END DO
              YX=YY+1.E-5
	          RSON(IDN2)=RSON(IDN2)+YN(IDN1)*A1
              CON=RSON(IDN2)/YX
	          RSOP(IDN2)=RSOP(IDN2)+YP(IDN1)*A1
              COP=RSOP(IDN2)/YX
	          X2=Y2+DEP+Y2O
              X1=CON*X2
              RTO=X1/RSON(IDN2)
              IF(RTO>.1)THEN
                  X1=.1*RSON(IDN2)
                  CON=X1/X2
              END IF
              YN(IDO)=CON*Y2/A1
              DRTP=MIN(1.,(.001+Y2O)/(Y1+.001))
	          X3=MIN(RSON(IDN2),CON*Y2O)
	          SMM(113,MO,IDN2)=SMM(113,MO,IDN2)+X3
              VAR(113,IDN2)=X3
              YN(I1)=YN(I1)+X3/A3
              RSON(IDN2)=MAX(.01,RSON(IDN2)-X1)
              X1=COP*X2
              RTO=X1/RSOP(IDN2)
              IF(RTO>.1)THEN
                  X1=.1*RSOP(IDN2)
                  COP=X1/X2
              END IF
              YP(IDO)=COP*Y2/A1
	          X3=MIN(RSOP(IDN2),COP*Y2O)
	          SMM(114,MO,IDN2)=SMM(114,MO,IDN2)+X3
              VAR(114,IDN2)=X3
              YP(I1)=YP(I1)+X3/A3
              RSOP(IDN2)=MAX(.01,RSOP(IDN2)-X1)
              RSO3(IDN2)=RSO3(IDN2)+A1*QN(IDN1)
              CO3=RSO3(IDN2)/VV
              RSSP(IDN2)=RSSP(IDN2)+A1*QP(IDN1)
              CSP=RSSP(IDN2)/VV
              QN(IDO)=OFLO*CO3
	          X3=OFP*CO3
	          SMM(115,MO,IDN2)=SMM(115,MO,IDN2)+X3
	          QN(I1)=QN(I1)+X3/A3
              RSO3(IDN2)=MAX(.01,RSO3(IDN2)-QN(IDO)-X3)
              QP(IDO)=OFLO*CSP
	          X3=OFP*CSP
	          SMM(116,MO,IDN2)=SMM(116,MO,IDN2)+X3
	          QP(I1)=QP(I1)+X3/A3
              RSSP(IDN2)=MAX(.01,RSSP(IDN2)-QP(IDO)-X3)
              QN(IDO)=QN(IDO)/A1
              QP(IDO)=QP(IDO)/A1
          END IF    
      END IF    
      SMM(66,MO,IDN2)=SMM(66,MO,IDN2)+EVRS(IDN2)
      VAR(66,IDN2)=EVRS(IDN2)
      TREV=TREV+EVRS(IDN2)
      SAET=10.*SALA(IDN2)*AET
      TSAE=TSAE+SAET
      SMM(67,MO,IDN2)=SMM(67,MO,IDN2)+SP
      VAR(67,IDN2)=SP
      X1=DEP/RSBD(IDN2)
      RSVP(IDN2)=MAX(0.,RSVP(IDN2)-X1)
	  RSVE(IDN2)=MAX(0.,RSVE(IDN2)-X1)
	  DF=YB+Y1-Y2-Y2O-DEP-RSYS(IDN2)
	  IF(KFL(13)==0)RETURN
	  !IF(ABS(DF)>1.E-5)WRITE(KW(13),13)DF
	  WRITE(KW(13),12)IDN2,NBSA(IDN2),IY,MO,KDA,RFRA,QSA,QIX,EVRS(IDN2),SP,&
      OFLO,RSV(IDN2),RSVP(IDN2),RSVE(IDN2),Y1,Y2O,DEP,RSSA(IDN2),VBK(2),VBK(1)
      RETURN     
   12 FORMAT(1X,2I4,1X,3I4,2X,4F10.0,F10.3,4F10.0,4F10.2,2F10.3)
   13 FORMAT(1X,'!!!!!',E16.6)	 
	  END 
      