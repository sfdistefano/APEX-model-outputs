      SUBROUTINE SOLT
      ! APEX1905
      ! THIS SUBPROGRAM ESTIMATES DAILY AVEAGE TEMPERATURE AT THE CENTER
      ! OF EACH SOIL LAYER.
      USE PARM 
      DATA XLAG/.8/,P98/2./,P99/.5/,P100/1./
      XLG1=1.-XLAG
      F=ABD(ISA)/(ABD(ISA)+686.*EXP(-5.63*ABD(ISA)))
      DP=1.+2.5*F
      WW=.356-.144*ABD(ISA)
      B=LOG(.5/DP)
      WC=.001*SW(ISA)/(WW*Z(LID(NBSL(ISA),ISA),ISA))
      F=MAX(.2,EXP(B*((1.-WC)/(1.+WC))**2))
      DD=F*DP
      X4=.5*AMPX
      X8=(IDA-200)/PIT
      DST0(ISA)=AVT+X4*COS(X8)
      DSTC=DST0(ISA)
      ZZ=2.*DD
      XX=0.
      IF(ISLT==0)THEN
          XZ=.5*(TMX(IRF(ISA))-TMN(IRF(ISA)))*ST0(ISA)/10.
          X2=TX
          X3=(1.-BCV(ISA))*X2+BCV(ISA)*STMP(LID(2,ISA),ISA)
          TG=.5*(X2+X3)
          X1=TG-DST0(ISA)
          DO J=1,NBSL(ISA)
              ISL=LID(J,ISA)
              ZD=(XX+Z(ISL,ISA))/ZZ
              X5=EXP(-ZD)
              STMP(ISL,ISA)=AVT+X5*(X4*COS(X8-ZD)+X1)
              XX=Z(ISL,ISA)
          END DO
      ELSE        
          IF(ST0(ISA)>P99)THEN
              B2=P100/(25.-P99)
              XZ=1.+B2*(ST0(ISA)-P99)
          ELSE
              XZ=1.
          END IF
          IF(TX<0.)THEN
              X1=XZ*ABS(TX)+TX
              X2=TX+X1
          ELSE    
              X2=TX*XZ
          END IF    
          DST0(ISA)=AVT+X4*COS(X8)
          DSTC=(1.-BCV(ISA))*X2+BCV(ISA)*DST0(ISA)
          !DSTC=(1.-ALB(ISA))*(TX+(TMX(IRF(ISA))-TMN(IRF(ISA)))*SQRT(.03*SRAD(IRF(ISA))))
          X1=AVT-DSTC
          DO J=1,NBSL(ISA)
              ISL=LID(J,ISA)
              ZD=(XX+Z(ISL,ISA))/ZZ
              F=1.-EXP(-P98*ZD)
              STMP(ISL,ISA)=XLAG*STMP(ISL,ISA)+XLG1*(F*X1+DSTC)
              XX=Z(ISL,ISA)
          END DO
      END IF
      IF(KFL(45)>0)WRITE(KW(45),107)ISA,NBSA(ISA),IYR,MO,KDA,DD,TX,&
      SRAD(IRF(ISA)),DST0(ISA),DSTC,(STMP(LID(K,ISA),ISA),K=1,NBSL(ISA)),&
      SNO(ISA),CV(ISA),BCV(ISA)               
      RETURN
  107 FORMAT(1X,2I8,1X,3I4,50F10.3)                                              
      END