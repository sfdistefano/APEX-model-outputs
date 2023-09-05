      SUBROUTINE CGROW(JRT)
!     APEX1905
!     THIS SUBPROGRAM CALCUALTES LEAF AREA INDEX, HEAT UNITS, ROOT DEPTH
!     AND TEMPERATURE STRESS FOR THE CROP.
      USE PARM
!	  DIMENSION SLA0(MNC,MSA),XDLA0(MNC,MSA) 
      JRT=0
      II=IDC(JJK)
	  X1=DM(JJK,ISA)+1.E-10 
      CPR=UP1(JJK,ISA)/X1
      CNR=UN1(JJK,ISA)/X1
      AJWA(ISA)=1.
      XPHU=PHU(JJK,IHU(JJK,ISA),ISA)
      TGX(JJK,ISA)=MAX(0.,TX-TBSC(JJK))
      IF(MATX(JJK,ISA)>0)THEN
          XX=MAX(0.,1.-HUI(JJK,ISA))
      ELSE
          XX=1.
      END IF
      SLLG=MAX(SLAI(JJK,ISA),BLAI(JJK,ISA))
      HU(JJK,ISA)=HU(JJK,ISA)+XX*TGX(JJK,ISA)
      IF(DETIN(JJK)>0..AND.HUI(JJK,ISA)>DETIN(JJK))HU(JJK,ISA)=DETIN(JJK)*XPHU
      IF(JDA==JDHU.AND.II/=NDC(7).AND.II/=NDC(8).AND.II/=NDC(10))THEN
          HU(JJK,ISA)=XPHU*PRMT(26)
          PSTS(ISA)=MIN(0.,PSTS(ISA))
          IPST(ISA)=0
          MATX(JJK,ISA)=0
      END IF    
      HUI(JJK,ISA)=HU(JJK,ISA)/XPHU
      IF(HU(JJK,ISA)>XPHU)THEN
          WCYD=MAX(WCY(JJK),WCYD-EO*.002)
          IF(II==NDC(3).OR.II==NDC(6))THEN
              HU(JJK,ISA)=0.
              JRT=2
          ELSE
              IF(II==NDC(7).OR.II==NDC(8).OR.II==NDC(10))THEN
                  JRT=0
                  HU(JJK,ISA)=.9*HU(JJK,ISA)
              ELSE
                  JRT=1
              END IF
          END IF
          RETURN
      END IF
      NGD(JJK,ISA)=NGD(JJK,ISA)+1
      X1=HUI(JJK,ISA)
      IF(SCLM(36)>0.)X1=MIN(X1,SCLM(36))
      F2=X1/(X1+EXP(DLAP(1,JJK)-DLAP(2,JJK)*X1))
      IF(II==NDC(8).OR.II==NDC(10))THEN
          X1=HSM(ISA)/AHSM
          IF(SCLM(36)>0.)X1=MIN(X1,SCLM(36))
          F1=X1/(X1+EXP(DLAP(1,JJK)-DLAP(2,JJK)*X1))
	      F=F1
          XLAI(JJK,ISA)=MAX(.1,DMLX(JJK,ISA)*F2)
      ELSE
          F=F2    
      END IF
      IF(II==NDC(8).OR.II==NDC(7).OR.II==NDC(10))THEN
          F3=HUI(JJK,ISA)
	  ELSE
          F3=SQRT(F2+1.E-10) 
	  END IF
      FF=F-WLV(JJK,ISA)
      XX=FF*XLAI(JJK,ISA)
      X2=1.
      SLAX=0.
      X3=(SLAI(JJK,ISA)+.001)*CPHT(JJK,ISA)
      IF(IGO(ISA)==1)THEN
          SUM=X3
          TOT=STL(JJK,ISA)
          ADD=SLAI(JJK,ISA)
      ELSE
          SUM=0.
	      TOT=0.
	      ADD=0.
          DO I=1,IGO(ISA)
              K1=JE(I,ISA)
              IF(SLAI(K1,ISA)>SLAX)SLAX=SLAI(K1,ISA)
	          TOT=TOT+STL(K1,ISA)
              SUM=SUM+SLAI(K1,ISA)*CPHT(K1,ISA)
              ADD=ADD+SLAI(K1,ISA)
          END DO
          IF(SLAX>2.)X2=X3/SUM
	  END IF
      IF(XX>0.)THEN
          X1=XX*X2*(1.+HR1)**3
          IF(II/=NDC(7).AND.II/=NDC(8).AND.II/=NDC(10))&
          X1=X1*SQRT(REG(JJK,ISA))*SHRL(ISA)
          SLAI(JJK,ISA)=MIN(XLAI(JJK,ISA),SLAI(JJK,ISA)+X1)
      END IF
      WLV(JJK,ISA)=F
      IF(SLAI(JJK,ISA)<BLAI(JJK,ISA))THEN
          SLAI(JJK,ISA)=BLAI(JJK,ISA)
      ELSE    
          !CPHT(JJK,ISA)=MAX(CPHT(JJK,ISA),HMX(JJK)*F3*REG(JJK,ISA))
          FF=MAX(0.,F3-WCHT(JJK,ISA))
          CPHT(JJK,ISA)=MIN(HMX(JJK),CPHT(JJK,ISA)+FF*HMX(JJK)*REG(JJK,ISA))
          IF((HUI(JJK,ISA)>XDLAI(JJK).OR.MATX(JJK,ISA)==1).AND.HRLT>HLMN(ISA).AND.XDLA0(JJK,ISA)>0.)THEN
              IF(MATX(JJK,ISA)==1)THEN
                  XXCJZ=.95
              ELSE
                  XXCJZ=HUI(JJK,ISA)
              END IF
              XX=(1.-XXCJZ)/XDLA0(JJK,ISA)
              IF(XX>1.E-5)THEN
                  XX=LOG10(XX)
              ELSE
                  XX=-5.
              END IF 
              IF(II/=NDC(7).AND.II/=NDC(8).AND.II/=NDC(10))THEN
                  RTO=RLAD(JJK)*XX
                  IF(RTO<-10.)RTO=-10.
                  SLAINEW=SLA0(JJK,ISA)*10.**RTO
				  SLAILOSS=SLLG-SLAINEW  
			      IF(MATX(JJK,ISA)==1)THEN				  
				      SLAILOSS=MAX(0.,MIN(SLAILOSS,SLAIOLD(JJK,ISA)))
                      SLAI(JJK,ISA)=SLAI(JJK,ISA)-SLAILOSS
                  ELSE           
                      SLAI(JJK,ISA)=MIN(SLAI(JJK,ISA),SLAINEW)    
                      SLAI(JJK,ISA)=MAX(BLAI(JJK,ISA),SLAI(JJK,ISA))			  
                      SLLN=SLLG-SLAI(JJK,ISA)
                      IF(SLLN<1.E-10)THEN
                          SLAI(JJK,ISA)=SLLG*(1.-SLAL(JJK,ISA)/(SLLG+SLAL(JJK,ISA)))
                          SLAI(JJK,ISA)=MAX(BLAI(JJK,ISA),SLAI(JJK,ISA))
                      END IF 
                  END IF    
				  SLAL(JJK,ISA)=SLLG-SLAI(JJK,ISA)
                  SLAIOLD(JJK,ISA)=MAX(SLAIOLD(JJK,ISA)-SLAILOSS,0.)				  
              END IF 
              RTO=RBMD(JJK)*XX
              IF(RTO<-10.)RTO=-10.
              AJWA(ISA)=10.**RTO  
          END IF
          IF(HUI(JJK,ISA)<XDLAI(JJK))THEN    
              XDLA0(JJK,ISA)=1.-XDLAI(JJK)
              SLA0(JJK,ISA)=SLAI(JJK,ISA)
              IF(MATX(JJK,ISA)== 0)SLAIOLD(JJK,ISA)=SLA0(JJK,ISA)
          END IF
      END IF
      WCHT(JJK,ISA)=F3
      XX=MAX(CPHT(JJK,ISA),RD(JJK,ISA),2.5*RDMX(JJK)*HUI(JJK,ISA))
      RD(JJK,ISA)=MIN(RDMX(JJK),Z(LID(NBSL(ISA),ISA),ISA),XX)
      IF(SCLM(23)>0.)ADD=MIN(ADD,SCLM(23))
	  FGC(ISA)=ADD/(ADD+EXP(SCRP(23,1)-SCRP(23,2)*ADD))
      IF(SCLM(24)>0.)TOT=MIN(TOT,SCLM(24))
	  FGSL(ISA)=TOT/(TOT+EXP(SCRP(24,1)-SCRP(24,2)*TOT))
      X1=HUI(JJK,ISA)
      IF(SCLM(32)>0.)X1=MIN(X1,SCLM(32))
      CLG(JJK,ISA)=BLG(3,JJK)*X1/(X1+EXP(BLG(1,JJK)-BLG(2,JJK)*X1))
      X1=TOPC(JJK)-TBSC(JJK)
      RTO=TGX(JJK,ISA)/X1
      IF(RTO<2..AND.TGX(JJK,ISA)>0.)THEN
          REG(JJK,ISA)=SIN(1.5707*RTO)
      ELSE
          REG(JJK,ISA)=0.
      END IF
      SHRL(ISA)=1.
      FHR=0.
      ! WINTER DORMANCY & FROST
      IF(HRLT+1.E-5<WDRM(ISA))THEN
          SHRL(ISA)=0.
          FHR=1.-HRLT/WDRM(ISA)
      ELSE
          SRA(JJK,ISA)=SRA(JJK,ISA)+SRAD(IRF(ISA))
      END IF
      ! NO EFFECT ON EVERGREEN TREES
      IF(II/=NDC(7))THEN
          IF(TMN(IRF(ISA))<-1.)THEN
              XX=ABS(TMN(IRF(ISA)))
              IF(SCLM(37)>0.)XX=MIN(XX,SCLM(37))
              F=XX/(XX+EXP(FRST(1,JJK)-FRST(2,JJK)*XX))
              F=MAX(F,FHR)
          ELSE    
              F=FHR
          END IF   
          IF(F>0.)THEN
              ! ONLY LAI REDUCED FOR DECIDUOUS TREES
              IF(II/=NDC(8).AND.II/=NDC(10))THEN
                  IF(STL(JJK,ISA)>1.E-5)THEN
                      XX=F*STL(JJK,ISA)
                      STL(JJK,ISA)=STL(JJK,ISA)-XX
                      DM(JJK,ISA)=DM(JJK,ISA)*SQRT(1.-F)
                      !DM(JJK,ISA)=DM(JJK,ISA)-XX
                      RW(JJK,ISA)=DM(JJK,ISA)-STL(JJK,ISA)
                      STD(JJK,ISA)=STD(JJK,ISA)+XX
                      STDL(JJK,ISA)=STDL(JJK,ISA)+CLG(JJK,ISA)*F
                      XY=XX*CNR
                      XZ=XX*CPR
                      XUN=UN1(JJK,ISA)
                      IF(XUN-XY<.01)XY=XUN-.01
                      XUP=UP1(JJK,ISA)
                      IF(XUP-XZ<.01)XZ=XUP-.01
                      STDN(JJK,ISA)=STDN(JJK,ISA)+XY
                      STDP(JJK,ISA)=STDP(JJK,ISA)+XZ
                      UN1(JJK,ISA)=XUN-XY
                      UP1(JJK,ISA)=XUP-XZ
                  END IF
                  IF(F*(1.-SNOF)>.99)THEN
                      IF(II/=NDC(3).AND.II/=NDC(6))CALL TRDST
                  END IF
              END IF
              SLAI(JJK,ISA)=SLAI(JJK,ISA)*(1.-F)    
          END IF    
      END IF    
      IF(REG(JJK,ISA)>0.)RETURN
      SFCP(5,JJK,ISA)=SFCP(5,JJK,ISA)+1.
      SFMO(5,JJK,ISA)=SFMO(5,JJK,ISA)+1.
      SMMC(14,JJK,MO,ISA)=SMMC(14,JJK,MO,ISA)+1.
      SMMC(17,JJK,MO,ISA)=SMMC(17,JJK,MO,ISA)+1.
      VARC(17,JJK,ISA)=0.
      JRT=1
      RETURN
      END