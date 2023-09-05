      SUBROUTINE HVOLQ(IVR)
!     APEX1905
!     THIS SUBPROGRAM PREDICTS DAILY RUNOFF VOLUME AND PEAK RUNOFF RATE
!     GIVEN DAILY PRECIPITATION AND SNOW MELT.
      USE PARM
      REAL::WEIRL,PONDINGH,PADDYQ
      WSAX1=10.*WSA(ISA)
      IF(IVR==0)THEN
          LD1=LID(1,ISA)
          IF(LUN(ISA)==35)THEN
              SCN=25400./CN0(ISA)-254.
              CN=CN0(ISA)
          ELSE
              SUM=0.
              ADD=0.
              IJ=NVCN(ISA)
              SELECT CASE(IJ)
                  CASE(0)
                      XX=0.
                      DO JJ=1,NBSL(ISA)
                          ISL=LID(JJ,ISA)
                          IF(Z(ISL,ISA)>1.)EXIT
                          ZZ=(Z(ISL,ISA)-XX)/Z(ISL,ISA)
                          SUM=SUM+(SWST(ISL,ISA)-S15(ISL,ISA))*ZZ/(FC(ISL,ISA)-S15(ISL,ISA))
                          ADD=ADD+ZZ
                          XX=Z(ISL,ISA)
                      END DO
                      IF(JJ<=NBSL(ISA))THEN 
                          ZZ=1.-XX
                          SUM=SUM+(SWST(ISL,ISA)-S15(ISL,ISA))*ZZ/(FC(ISL,ISA)-S15(ISL,ISA))
                          ADD=ADD+ZZ
                      END IF
                      SUM=SUM/ADD
                      IF(SUM>0.)THEN
                          SUM=100.*SUM
                          IF(SUM<1000.)THEN
                              IF(SCLM(30)>0.)SUM=MIN(SUM,SCLM(30))
                              SCN=SMX(ISA)*(1.-SUM/(SUM+EXP(CNSC(1,ISA)-CNSC(2,ISA)*SUM)))
                          ELSE
                              SCN=3.
                          END IF
                      ELSE   
                          SCN=SMX(ISA)*(1.-SUM)**2
                      END IF
                      SCN=SCNX(ISA)*SCN   
                  CASE(1)                           
                      DO JJ=1,NBSL(ISA)
                          ISL=LID(JJ,ISA)
                          IF(Z(ISL,ISA)>1.)EXIT
                          SUM=SUM+SWST(ISL,ISA)-S15(ISL,ISA)
                          ADD=ADD+FC(ISL,ISA)-S15(ISL,ISA)
                          L1=ISL
                      END DO
                      IF(JJ<=NBSL(ISA))THEN
                          RTO=(1.-Z(L1,ISA))/(Z(ISL,ISA)-Z(L1,ISA))
                          SUM=SUM+(SWST(ISL,ISA)-S15(ISL,ISA))*RTO
                          ADD=ADD+(FC(ISL,ISA)-S15(ISL,ISA))*RTO 
                      END IF
                      SUM=SUM/ADD
                      IF(SUM>0.)THEN
                          SUM=100.*SUM
                          IF(SUM<1000.)THEN
                              SCN=SMX(ISA)*(1.-SUM/(SUM+EXP(CNSC(1,ISA)-CNSC(2,ISA)*SUM)))
                          ELSE
                              SCN=3.
                          END IF
                      ELSE   
                          SCN=SMX(ISA)*(1.-SUM)**2
                      END IF
                      SCN=SCNX(ISA)*SCN   
                  CASE(2)
                      DO JJ=1,NBSL(ISA)
                          ISL=LID(JJ,ISA)
                          SUM=SUM+SWST(ISL,ISA)-S15(ISL,ISA)
                          ADD=ADD+FC(ISL,ISA)-S15(ISL,ISA)
                          IF(Z(ISL,ISA)>1.)EXIT
                      END DO
                      SUM=SUM/ADD
                      IF(SUM<0.)THEN
                          SCN=SMX(ISA)*(1.-SUM)**2
                      ELSE    
                          RTO=MIN(.98,SUM)
                          SCN=SMX(ISA)*(1.-RTO)
                      END IF    
                  CASE(3)
                      SCN=25400./CN0(ISA)-254.
                      CN=CN0(ISA)     
                  CASE(4)
                      SCN=SCI(ISA)
              END SELECT
              IF(IJ/=3)THEN        
                  X1=PRMT(15)*(PRMT(35)-RSD(LD1,ISA))
                  IF(X1>-10.)THEN
                      SCN=MAX(3.,(SCN-SMX(ISA))*EXP(X1)+SMX(ISA))  
                  ELSE
                      SCN=SMX(ISA)
                  END IF
                  IF(STMP(LID(2,ISA),ISA)<-1.)THEN
                      SCN=SCN*PRMT(22)
                  ELSE
                      IF(RFV(IRF(ISA))>0.)SCN=SCN*EXP(PRMT(25)*(.2-AL5))
	              END IF
                  CN=25400./(SCN+254.)
                  IF(ISCN==0)THEN
                      UPLM=MIN(99.5,CN+5.)
                      BLM=MAX(1.,CN-5.)
                      CN=ATRI(BLM,CN,UPLM,IDG(8))
                  END IF
                  SCN=25400./CN-254.
                  IF(SCN<3.)THEN
                      SCN=3.
                      CN=25400./(SCN+254.)
                  END IF    
              END IF
          END IF
          SCN2=PRMT(20)*SCN    
          TOT=100.
          DO I=1,9
              TOT=TOT-5.
              IF(CN>TOT)EXIT
          END DO
          CNDS(I,ISA)=CNDS(I,ISA)+1.
          IF(RFV(IRF(ISA))>0.)THEN
              SELECT CASE(INFL)
	              CASE(0)
	                  X1=RFV(IRF(ISA))-SCN2
                      IF(X1>0.)THEN
                          QVOL(IDO)=X1*X1/(RFV(IRF(ISA))+(1.-PRMT(20))*SCN)
                      ELSE
                          QVOL(IDO)=0.
                      END IF        
	              CASE(1,2)
	                  CALL HREXP
	              CASE(3)
	                  CALL HRUNF
                  CASE(4)
                      CALL HRFDTQ
              END SELECT
          END IF    
      ELSE
          QVOL(IDO)=EFI(ISA)*RFV(IRF(ISA))
      END IF
      IF(IFD(ISA)>0)THEN
          IF(DHT(ISA)>.1)THEN
              CALL HFURD
              X1=DVOL-MAX(0.,SWST(LD1,ISA)-PO(LD1,ISA))
              IF(QVOL(IDO)<X1)THEN
                  QVOL(IDO)=0.
                  IF(DHT(ISA)/DKHL(ISA)<.7)DHT(ISA)=DKHL(ISA)
              ELSE
                  DHT(ISA)=0.
                  IF(IDRL(ISA)==0.AND.CPHT(JJK,ISA)<1.)DHT(ISA)=DKHL(ISA)   
              END IF
              IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),17)ISA,NBSA(ISA),IYR,&
              MO,KDA,DHT(ISA),RHTT(ISA),QVOL(IDO),DVOL,X1,XHSM(ISA)
          END IF
      END IF
      !Paddy field water storage Jaehak jeong 2014
      IF(PADDY_HWEIR(ISA)>0)THEN 
          QVOL(IDO)=0.
          IF(PADDY_STO(1,ISA)>PADDY_HWEIR(ISA))THEN
              DO I=1,96
                  PONDINGH=(PADDY_STO(1,ISA)-PADDY_HWEIR(ISA))/1000.
                  IF(PONDINGH<=0.)EXIT
                  ! RECTANGULAR WEIR FOR DISCHARGE CONTROL
                  PADDYQ=1.84*PADDY_LWEIR(ISA)*PONDINGH**2.6
                  PADDYQ=900.*PADDYQ/WSAX1
                  QVOL(IDO)=QVOL(IDO)+PADDYQ
                  PADDY_STO(1,ISA)=PADDY_STO(1,ISA)-PADDYQ
              END DO    
          ELSE
              QVOL(IDO)=0.
          END IF
          IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),17)ISA,NBSA(ISA),IYR,&
          MO,KDA,DHT(ISA),RHTT(ISA),QVOL(IDO),DVOL,X1,XHSM(ISA)
      ELSE
          IF(PADDY_STO(1,ISA)>0)THEN
              QVOL(IDO)=QVOL(IDO)+PADDY_STO(1,ISA)+RFV(IRF(ISA))
              PADDY_STO(1,ISA)=0.
          END IF
      END IF
      !Urban impervious cover
      IF(URBF(ISA)>0.)THEN
          X1=RFV(IRF(ISA))-2.67
          IF(X1>0.)QURB(IDO)=X1*X1/(RFV(IRF(ISA))+10.69)
          QVOL(IDO)=QVOL(IDO)*(1.-URBF(ISA))+QURB(IDO)*URBF(ISA)
      END IF
      IF(RFV(IRF(ISA))>0.)THEN
          X2=QVOL(IDO)/DUR
          IF(X2>1.)THEN
              X2=X2**.25
          ELSE
              X2=1.
          END IF
          X4=MIN(SPLG(ISA)/360.,TCS(ISA)/X2)
          TC(IDO)=X4+TCC(ISA)/X2
          SELECT CASE(ITYP)
              CASE(1,2,3,4)
                  TC(IDO)=TCC(ISA)+TCS(ISA)/SQRT(RFV(IRF(ISA)))
                  QRB=MIN(.95,SCN2/RFV(IRF(ISA)))
                  CALL HTR55
                  ALTC=MAX(ALMN,1.-EXP(-TC(IDO)*PRFF))
                  X1=ALTC*QVOL(IDO)
              CASE(0,-1)
                  ALTC=MAX(ALMN,1.-EXP(-TC(IDO)*PRFF))
                  X1=ALTC*QVOL(IDO)
                  RQRB(IDO)=X1/TC(IDO)
              CASE(5)
                  IF(QVOL(IDO)>0.)THEN
                      CALL HRFDTTRI
                      !CALL HRFDTS
                      !CALL HRFDT
                      QV=0.
                      S4=0.
                      XX=0.
                      QMAX=0.
                      DO I=2,NRF
                          X1=RFDT(I)-SCN2
                          IF(X1>0.)THEN
                              QV=X1*X1/(RFDT(I)+.8*SCN)
                          ELSE
                              QV=0.
                          END IF
                          QMM=QV-XX
                          IF(QMM>QMAX)QMAX=QMM
                          S4=S4+QMM*QMM 
                          XX=QV
                      END DO
                      WQRT=S4/(QVOL(IDO)*DTHY)                                                                 
                      DURE=QVOL(IDO)/WQRT
                      !RTO=QVOL(IDO)/RFV(IRF(ISA))
                      !DURE=DUR*RTO
                      TPBT=.5*DURE+.6*TC(IDO)
                      RQRB(IDO)=2.*BETA(ISA)*QVOL(IDO)/TPBT
                  END IF        
                  X1=2.*BETA(ISA)*QVOL(IDO)
          END SELECT        
          QPR(IDO)=X1/X4
          IF(DALG(ISA)>0.)QVOL(IDO)=QVOL(IDO)*(X1-DALG(ISA))/X1
          IF(LUN(ISA)==35)ES=RFV(IRF(ISA))-QVOL(IDO)
      ELSE
          ES=0.
      END IF
      QVOL(IDO)=WSAX1*QVOL(IDO)
      QURB(IDO)=WSAX1*QURB(IDO)
      IF(NVCN(ISA)/=4)SCI(ISA)=SCN
      RETURN
   17 FORMAT(1X,2I8,1X,3I4,12X,'DHT= ',F6.0,'mm',2X,'RHTT= ',F6.0,&
      'mm',2X,'Q= ',F6.1,'mm',2X,'DV= ',F7.1,'mm',2X,'ADV= ',F7.1,'mm',&
      2X,'HUSC = ',F5.2)
      END