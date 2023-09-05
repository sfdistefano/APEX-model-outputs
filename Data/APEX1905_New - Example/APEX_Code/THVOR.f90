      SUBROUTINE THVOR(X6,JRT)
!     APEX1905
!     THIS SUBPROGRAM HARVESTS CROPS WITH OVER RIDE HARVEST INDEX--
!     FORAGE(HAY, SILAGE, ETC) AND GRAZING.
      USE PARM 
!      DIMENSION YLSD(MNC,MSA)
      JRT=0
      II=IHC(JT1)
      LD1=LID(1,ISA)
      NN=NCP(IRO(ISA),ISA)
      YY=0.
      YLK=0.
      YLN=0.
      YLP=0.
      YLR=0.
      YLRN=0.
      YLRP=0.
 	  IOW=IDON(ISA)
      JJ=JJK
      DO
          IF(KGO(JJ,ISA)==0)EXIT
          IF(IDC(JJ)==NDC(7).OR.IDC(JJ)==NDC(8).OR.IDC(JJ)==NDC(10).OR.II==NHC&
          (3))THEN
              IF(JH(IRO(ISA),KT(ISA),ISA)/=KDC(JJ))EXIT
              I1=LYR(IRO(ISA),KT(ISA),ISA)
              IF(IYH(JJ,ISA)<I1.AND.I1/=1)THEN
                  KOMP(KT(ISA),ISA)=0
                  JRT=1
                  RETURN
              END IF
          END IF
          IF(II==NHC(3))THEN
              IF(IHT(JT1,ISA)>0)RETURN
              IHT(JT1,ISA)=1
          END IF
          IF(JP(JJ,ISA)<=0)THEN
              JP(JJ,ISA)=1
              NCR(JJ,ISA)=NCR(JJ,ISA)+1
          END IF
          HUF(JJ,ISA)=MAX(HUF(JJ,ISA),HU(JJ,ISA))
          DMF(JJ,ISA)=DM1(JJ,ISA)
          TRA(JJ,ISA)=SRA(JJ,ISA)+TRA(JJ,ISA)
          IF(RD(JJ,ISA)>RDF(JJ,ISA))RDF(JJ,ISA)=RD(JJ,ISA)
          XX=DM(JJ,ISA)+.001
          X3=RW(JJ,ISA)
          X12=STL(JJ,ISA)
          X2=MAX(.001,UN1(JJ,ISA)/XX)
          X7=X2*.001
          X10=MAX(.0001,UP1(JJ,ISA)/XX)
          X8=STD(JJ,ISA)+1.E-10
          RNR=STDN(JJ,ISA)/X8
          RPR=STDP(JJ,ISA)/X8
          X1=ORHI(JT1)
          HEX=HE(JT1)
          IF(HEX<0.)THEN
              HEX=-HEX
              CRSN=(WLSN(LD1,ISA)+WLMN(LD1,ISA))/RSD(LD1,ISA)
              CRSP=FOP(LD1,ISA)/RSD(LD1,ISA)
              YLR=X1*RSD(LD1,ISA)*HEX 
              YLRN=YLR*CRSN
              YLRP=YLR*CRSP
              FOP(LD1,ISA)=FOP(LD1,ISA)-YLRP
              CALL NCNSTD(-YLR,-YLRN,LD1)
              AJHI(JJ,ISA)=0.
          ELSE
              IF(TLD(JT1)>0.)THEN
!                  CALL THVRT(YLSD,X6,X1,X3,JJ,LD1)
                  CALL THVRT(X6,X1,X3,JJ,LD1)
                  GO TO 8
              END IF
              IF(II==NHC(22))THEN
                  KOMP(KT(ISA),ISA)=0	
                  IF(AGPM(ISA)<PRMT(79).OR.NMW(ISA)<IMW(ISA))EXIT
                  NMW(ISA)=0
              END IF
          END IF
          ZZ=MAX(.01,1.-X1)
          YZ=X1*X8      
          XZ=X1*X12
          IF(II/=NHC(3))THEN
              IF(IDC(JJ)/=NDC(7).AND.IDC(JJ)/=NDC(8).AND.IDC(JJ)/=NDC(10))THEN
                  CPHT(JJ,ISA)=MAX(.001,CPHT(JJ,ISA)*ZZ)
                  HU(JJ,ISA)=MAX(.1,HU(JJ,ISA)*ZZ)
                  SLAI(JJ,ISA)=MAX(.01,SLAI(JJ,ISA)*ZZ)
              END IF
          END IF
          STD(JJ,ISA)=MAX(.01,STD(JJ,ISA)*ZZ)
          YLD(JJ,ISA)=XZ*HEX
          YLSD(JJ,ISA)=YZ*HEX
          Y4=YZ*RNR
          Y5=YZ*RPR
          STDN(JJ,ISA)=MAX(0.,STDN(JJ,ISA)-Y4)
          STDP(JJ,ISA)=MAX(0.,STDP(JJ,ISA)-Y5)
          X4=MIN(XZ*X2,UN1(JJ,ISA))
          X5=MIN(XZ*X3,UP1(JJ,ISA))
          Z2=YLSD(JJ,ISA)*RNR
          Z3=YLSD(JJ,ISA)*RPR
          YLN(JJ,ISA)=MIN(.9*(UN1(JJ,ISA)+STDN(JJ,ISA)),YLD(JJ,ISA)*X2)
          YLP(JJ,ISA)=MIN(.9*(UP1(JJ,ISA)+STDP(JJ,ISA)),YLD(JJ,ISA)*X3)
          X11=XZ-YLD(JJ,ISA)+YZ-YLSD(JJ,ISA)
          X10=X4-YLN(JJ,ISA)+Y4-Z2
          JJK=JJ
          CALL NCNSTD(X11,X10,LD1)
          FOP(LD1,ISA)=MAX(.01,FOP(LD1,ISA)+X5-YLP(JJ,ISA)+Y5-Z3)
          YY=YLD(JJ,ISA)+YLSD(JJ,ISA)+YLR
          YLD2(JJ,ISA)=YLD2(JJ,ISA)+YY
          JD(ISA)=JJ
          SRA(JJ,ISA)=0.
          UN1(JJ,ISA)=UN1(JJ,ISA)-X4
          UP1(JJ,ISA)=UP1(JJ,ISA)-X5
          DM(JJ,ISA)=DM(JJ,ISA)-XZ
          IF(DM(JJK,ISA)<RW(JJK,ISA))RW(JJK,ISA)=DM(JJK,ISA)
          STL(JJ,ISA)=DM(JJ,ISA)-RW(JJ,ISA)
          YLN(JJ,ISA)=YLN(JJ,ISA)+Z2+YLRN
          YLP(JJ,ISA)=YLP(JJ,ISA)+Z3+YLRP
          YLNF(JJ,ISA)=YLNF(JJ,ISA)+YLN(JJ,ISA)
          YLPF(JJ,ISA)=YLPF(JJ,ISA)+YLP(JJ,ISA)
          IYH(JJ,ISA)=0
        8 TYN(ISA)=TYN(ISA)+YLN(JJ,ISA)
          TYP(ISA)=TYP(ISA)+YLP(JJ,ISA)
          TYK(ISA)=TYK(ISA)+YLK(JJ,ISA)
          IF (TIL(JT1)=='BALEPTOL'.OR.TIL(JT1)=='BALERMTD'.OR.TIL(JT1)=='BALEPTOS') then
		  Cut_hay(ISA)=Cut_hay(ISA)+YLD(JJ,ISA)!+YLSD may not needed  !add a variable to contain how much hay can be avalialbe for each site  T/ha
          ENDIF  
          IF(NOP>0.OR.NBSA(ISA)==ISAP)WRITE(KW(1),7)ISA,NBSA(ISA),IYR,MO,KDA,TIL&
          (JT1),CPNM(JD(ISA)),YLD(JJ,ISA),YLSD(JJ,ISA),AGPM(ISA),X12,X3,X1,X7,XHSM(ISA),YLN,YLP
          EXIT
      END DO    
      RETURN
    3 FORMAT(1X,2I8,1X,I4,2I2,3I4,1X,A8,1X,A4,6F8.4)
    7 FORMAT(1X,2I8,1X,3I4,2X,A8,2X,A4,2X,'YLD=',F7.2,'t/ha',2X,'YSD=',&
      F7.2,'t/ha',2X,'AGPM=',F7.2,'t/ha',2X,'STL=',F7.2,'t/ha',2X'RWT=',&
      F7.2,'t/ha',2X,'HI=',F7.3,2X,'NCN=',F7.3,'G/G',2X,'HUSC=',F5.2,2X,&
      'YLN=',F5.0,'kg/ha',2X,'YLP=',F5.0,'kg/ha')
      END