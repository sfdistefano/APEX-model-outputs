      SUBROUTINE NFERT(IRC,JFT,KKG,KTX,JRT)
!     APEX1905
!     THIS SUBPROGRAM APPLIES N AND P FERTILIZER AT SPECIFIED DATES,
!     RATES, AND DEPTH OR AUTOMATICALLY.
      USE PARM
      JRT=0
      IF(LUN(ISA)==35)RETURN
!     IF(ANA(IRO(ISA),ISA)>=FNMX(IRO(ISA),ISA).AND.IRC/=1)RETURN
      X1=0.
      X8=0.
      I=LID(1,ISA)
      ZFT=TLD(JFT)
      WSAX=WSA(ISA)
      SELECT CASE(IRC)
	      CASE(1,2)
	          I1=IABS(IAPL(ISA))
              KF=IDFT(IRC,I1)
	      CASE(3)
	          KF=IDMU(KKG)
	      CASE(4,5)
	          KF=IDFT(IRC,ISA)
	      CASE(6)
	          KF=LFT(IRO(ISA),KTX,ISA)
              X1=WFA(IRO(ISA),KTX,ISA)
              GO TO 5
          CASE(7)
              KF=IDFT(3,ISA)
              X1=APMU(ISA)
	          X3=0.
	          GO TO 12              
	      CASE(8)
              KF=IDFT(6,ISA)
              X1=APMU(ISA)
	          X3=0.
	          GO TO 12
          CASE(9)
              KF=MFT
              X1=APMU(ISA)
              X3=X1*FN(KF)
              FTNM(KF)='FECE'
              GO TO 12
      END SELECT
      IF(IRC/=3.AND.IRC/=4.AND.IRC/=7.AND.MNUL/=0)THEN
          XX=FP(KF)+FPO(KF)
          XX1=XMAP(ISA)/XX
          XX2=SAMA(ISA)/XX
          IF(XX2>=XX1)THEN
              JRT=1
              RETURN 
          END IF
          APMU(ISA)=MIN(APMU(ISA),XX1-XX2)
      END IF
      X1=APMU(ISA)
    5 IF(FOC(KF)<1.E-10)THEN
          DO J=1,NBSL(ISA)
              I=LID(J,ISA)
              IF(ZFT<Z(I,ISA))EXIT
          END DO

          IF(X1<1.E-10)THEN
              X3=MAX(UNA(JJK,ISA)-TNOR(ISA),0.)
              X1=X3/FN(KF)
              IF(X1>0.)GO TO 1
              RETURN
          END IF
      END IF
      X3=X1*FN(KF)
      GO TO 12
    1 X8=FNMX(IRO(ISA),ISA)-ANA(IRO(ISA),ISA)
      IF(X8<1.E-10)THEN
          JRT=1
          RETURN 
      END IF
      IF(X3<=X8)GO TO 12
      X3=X8
      X1=X3/FN(KF)
   12 X2=X1*FP(KF)
      X9=X1*FK(KF)
      SOLK(I,ISA)=SOLK(I,ISA)+X9
      SMM(148,MO,ISA)=SMM(148,MO,ISA)+X9
      VAR(148,ISA)=X9
      X4=X3*FNMA(KF)
      X5=X1*FNO(KF)
      X6=X1*FPO(KF)
      X7=X3-X4
      X8=X1*FOC(KF)
      X11=X1*FSLT(KF)
      XX8=X8*WSAX
      SMM(99,MO,ISA)=SMM(99,MO,ISA)+XX8
      VAR(99,ISA)=XX8
      IF(X8<.1)THEN
          WPML(I,ISA)=WPML(I,ISA)+X2
          WNO3(I,ISA)=WNO3(I,ISA)+X7
          WNH3(I,ISA)=WNH3(I,ISA)+X4
      ELSE
          WHSC(I,ISA)=WHSC(I,ISA)+X8
          WPMU(I,ISA)=WPMU(I,ISA)+X2
          WNMU(I,ISA)=WNMU(I,ISA)+X3
          Z1=.001*X1
          RSDM(I,ISA)=RSDM(I,ISA)+Z1
          WCOU(I,ISA)=WCOU(I,ISA)+X8
          WNOU(I,ISA)=WNOU(I,ISA)+X5
          WPOU(I,ISA)=WPOU(I,ISA)+X6
          SMM(78,MO,ISA)=SMM(78,MO,ISA)+Z1
          X53=WSAX*X5
          SMM(53,MO,ISA)=SMM(53,MO,ISA)+X53
          VAR(53,ISA)=VAR(53,ISA)+X53
          X56=WSAX*X6
          SMM(56,MO,ISA)=SMM(56,MO,ISA)+X56
          VAR(56,ISA)=VAR(56,ISA)+X56
      END IF
      ANA(IRO(ISA),ISA)=ANA(IRO(ISA),ISA)+X3
      X54=X7*WSAX
      SMM(54,MO,ISA)=SMM(54,MO,ISA)+X54
      VAR(54,ISA)=VAR(54,ISA)+X54
      X55=X4*WSAX
      SMM(55,MO,ISA)=SMM(55,MO,ISA)+X55
      VAR(55,ISA)=VAR(55,ISA)+X55
      X57=WSAX*X2
      SMM(57,MO,ISA)=SMM(57,MO,ISA)+X57
      VAR(57,ISA)=VAR(57,ISA)+X57
      IF(ZFT<.01)THEN
          FSFN(ISA)=FSFN(ISA)+X5+X7+X4
          FSFP(ISA)=FSFP(ISA)+X6+X2
      END IF
      WSLT(I,ISA)=WSLT(I,ISA)+X11
      SMM(132,MO,ISA)=SMM(132,MO,ISA)+X11
      IF(IRC/=3)THEN
          FRTN(JJK,ISA)=FRTN(JJK,ISA)+X3+X5
          FRTP(JJK,ISA)=FRTP(JJK,ISA)+X2+X6
          FRTK(JJK,ISA)=FRTK(JJK,ISA)+X9
      END IF
      VAR(132,ISA)=X11
      IF(IRC/=9)THEN
	      XX=X1*FCST(KF)
          COST(ISA)=COST(ISA)+XX
          IF(IRC==4)THEN
              Y1=COTL(JFT)
              Y2=Y1-COOP(JFT)
              COST(ISA)=COST(ISA)+Y1
              CSFX=CSFX+Y2
          END IF
      ELSE
          XX=0.
      END IF
      IF(KFL(31)>0)THEN
          WRITE(KW(31),34)ISA,NBSA(ISA),IYR,MO,KDA,FTNM(KF),KDC(JJK),&
          KDF(KF),IHC(JFT),NBE(JFT),NBT(JFT),XX,XX,X1
          IF(IRC==4)WRITE(KW(31),50)ISA,NBSA(ISA),IYR,MO,KDA,TIL(JFT),&
          KDC(JJK),IHC(JFT),NBE(JFT),NBT(JFT),Y1,Y2,FULU(JFT)
      END IF
      IF((NOP>0.OR.NBSA(ISA)==ISAP).AND.X3>0..AND.IRC==1)WRITE(KW(1),25)ISA,&
      NBSA(ISA),IYR,MO,KDA,IDON(ISA),NBSA(I1),FTNM(KF),X1,ZFT,X3,X4,&
      X5,X2,X6,X9,XHSM(ISA)
      IF(KFL(19)>0.AND.KKG>0)WRITE(KW(19),22)ISA,NBSA(ISA),IYR,MO,KDA,IY,&
      IDON(ISA),GNAM(KKG),FTNM(KF),X1,X3,X4,X5,X2,X6
      APMU(ISA)=X1*WSA(ISA)
      IF((IRC==6.AND.X8<.1).OR.IRC==4.OR.IRC==7)THEN
          FCMN(ISA)=FCMN(ISA)+X3+X5
          FCMP(ISA)=FCMP(ISA)+X2+X6
      ELSE
          SMM(81,MO,ISA)=SMM(81,MO,ISA)+X1
          SAMA(ISA)=SAMA(ISA)+X2+X6
          TMAP=TMAP+APMU(ISA)
	      OMAP(IDON(ISA))=OMAP(IDON(ISA))+APMU(ISA)
          IF(IRC==3)RETURN
	  END IF
      IF(IRC/=9)NDFA(ISA)=0
      RETURN
   11 FORMAT(1X,2I8,1X,I4,2I2,2X,'IDON=',I4,2X,'HRD#=',I3,2X,A8,2X,&
      'RATE=',F8.0,'kg/ha',1X,'DPTH=',F8.2,'mm ELEM WT(kg/ha)--',1X,'MN='&
      ,F7.0,1X,'NH3=',F7.0,1X,'ON=',F7.0,1X,'MP=',F7.0,1X,'OP=',F7.0,1X,&
      'MK=',F7.0,1X,'HUSC=',F5.2)
   22 FORMAT(1X,2I8,1X,2I8,3I4,1X,A16,3X,A8,F8.0,6F8.2)      
   25 FORMAT(1X,2I8,1X,I4,2I2,2X,'IDON=',I4,2X,'FA#=',I4,2X,A8,2X,&
      'RATE=',F8.0,'kg/ha',1X,'DPTH=',F8.2,'mm ELEM WT(kg/ha)--',1X,'MN='&
      ,F7.0,1X,'NH3=',F7.0,1X,'ON=',F7.0,1X,'MP=',F7.0,1X,'OP=',F7.0,1X,&
      'MK=',F7.0,1X,'HUSC=',F5.2)
   34 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,8X,I6,2X,4I4,F10.2,10X,2F10.2)
   50 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,8X,I6,6X,3I4,2F10.2,20X,F10.2)
      END