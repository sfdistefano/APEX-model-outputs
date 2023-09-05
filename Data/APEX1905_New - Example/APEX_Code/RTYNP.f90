      SUBROUTINE RTYNP(QAV,SUM,TOT,IZ)
!     APEX1905
!     THIS SUBPROGRAM ROUTES SEDIMENT THRU A REACH USING THE TIME STEP 
!     OF THE RTVSC. POTENTIAL SEDIMENT CONC IS ESTIMATED AS A FUNCTION 
!     OF VELOCITY. NET DEPOSITION OR DEGRADATION MAY RESULT DEPENDING ON
!     CONC OF REACH INFLOW, VEGETATIVE COVER, AND CHANNEL ERODIBILITY.
      USE PARM
      LD1=LID(1,IDN2)
      ZCH=RCHD(IDN2)
      CBW=RCBW(IDN2)
	  Y1=.5*(YHY(IZ,IDN1)+YHY(IZ-1,IDN1))
      Q1=QAV/CMS
	  CIN=.1*Y1/Q1
	  DEP=0.
      DEG=0.
      IF(QCAP(IDN2)<QAV)THEN
          QQ=QCAP(IDN2)
          QFP=QAV-QCAP(IDN2)
          DFP=(QFP/RFPX(IDN2))**.6
          VFP=QFP/(DFP*RFPW(IDN2))
          X1=QFP*VFP
          CYFP=PRMT(19)*VFP**PRMT(18)
          Q3=QFP/CMS
          SMM(79,MO,IDN2)=SMM(79,MO,IDN2)+Q3
          VAR(79,IDN2)=Q3
          VFPA(IDN2)=VFPA(IDN2)+VFP
          IF(VFP>VFPB(IDN2))VFPB(IDN2)=VFP
          NBFX=NBFX+1
          TRT=RFPL(IDN2)/(3.6*VFP)
          XX=FPSC(IDN2)
          F=MIN(Q3,XX*TRT)
          X2=F*RWSA(IDN1)/WSA(IDN2)
	      FPF(IDN2)=FPF(IDN2)+X2
          SMM(98,MO,IDN2)=SMM(98,MO,IDN2)+X2
          Q2=Q1-F
          QO=Q3-F
          IF(QO<1.E-5)THEN
              DEP=10.*CIN*Q3
              GO TO 34
          END IF
          X7=PRMT(45)*TRT*PSZM(IDN1)
          IF(X7<10.)THEN
              F=1.-EXP(-X7)
          ELSE
              F=1.
          END IF   
          X5=10.*(CYFP*QO-CIN*Q3)*F
          IF(X5>0.)THEN
              DEG=EK(IDN2)*CVF(IDN2)*X5
          ELSE
              DEP=-X5
          END IF
          GO TO 34
      END IF
      X1=0.
      QQ=QAV
   34 IF(QQ<1.E-5)THEN
          VCH=0.
          DEP=Y1
          GO TO 10
      END IF
      DCH=MAX(.1,ZCH*LOG10(QQ+1.)/LOG10(QCAP(IDN2)+1.))
      X3=SQRT(RCSS(IDN2)*RCSS(IDN2)+1.)
      DO IT=1,10
          A=DCH*(CBW+DCH*RCSS(IDN2))
          P=CBW+2.*X3*DCH
          Q=A**1.66667*RCHX(IDN2)/P**.66667
          FU=Q-QQ
          X6=MAX(1.,QQ)
          IF(ABS(FU/X6)<.01)GO TO 14
          R=A/P
          DFDD=RCHX(IDN2)*(1.6667*R**.6667*(CBW+2.*DCH*RCSS(IDN2))-1.333*R**&
         &1.66667*X3)
          DCH=DCH-FU/DFDD
      END DO
      WRITE(KW(1),16)IDN2,NBSA(IDN2),IYR,MO,KDA,QQ,Q,DFDD,DCH,ZCH
   14 VCH=QQ/A
      TRT=RCHL(IDN2)/(3.6*VCH)
      CYCH=PRMT(19)*VCH**PRMT(18)
      VCHA(IDN2)=VCHA(IDN2)+VCH
      IF(VCH>VCHB(IDN2))VCHB(IDN2)=VCH
      NBCX=NBCX+1
      Q3=QQ/CMS
      X7=PRMT(45)*TRT*PSZM(IDN1)
      IF(X7<10.)THEN
          F=1.-EXP(-X7)
      ELSE
          F=1.
      END IF
      X5=10.*(CYCH-CIN)*Q3
      IF(X5>0.)THEN
          X5=RCHK(IDN2)*X5
          DEG=DEG+X5
      ELSE
          DEP=DEP-X5*F
          DPMT(IDN2)=DPMT(IDN2)-X5
          ZCH=ZCH+.1*X5/BD(LD1,IDN2)
      END IF
      X7=2.*PRMT(45)*TRT
      IF(X7<10.)THEN
          F=1.-EXP(-X7)
      ELSE
          F=1.
      END IF   
      DO I=1,NBSL(IDN2)
          ISL=LID(I,IDN2)
          IF(ZCH<Z(ISL,IDN2))EXIT
      END DO
      ZCH=ZCH+.1*X5/BD(ISL,IDN2)
   10 SRCH(14,IDO)=SRCH(14,IDO)+DEG
      TDEG=TDEG+DEG
      DEP=MIN(Y1,DEP)
      SRCH(13,IDO)=SRCH(13,IDO)+DEP
      TDEP=TDEP+DEP
	  Y2=Y1-DEP+DEG
	  YHY(IZ,IDO)=Y2
	  YSD(NDRV,IDO)=YSD(NDRV,IDO)+Y2
	  SUM=SUM+Y2
	  TOT=TOT+Y1
	  RETURN
   12 FORMAT(1X,2I4,1X,I4,2I2,'FLOOD QI=',F7.2,' QO=',F7.2,' QAV=',F7.2,&
      ' QQ=',F7.2,' QFP=',F7.2,' DCH=',F7.2,' DFP=',F7.3,' VCH=',F7.3,&
      ' VFP=',F7.3,' CIN=',F8.5,' CCH=',F8.5,' CFP=',F8.5,' DEP=',F7.3,&
      ' DEG=',F7.3,' DR=',F7.3,' ER=',F7.2,' YI=',F7.3,' YO=',F7.3)
   16 FORMAT(1X,'_______________PROBLEM DID NOT CONVERGE______________',&
      2I8,1X,I4,2I2,5E13.5)
      END