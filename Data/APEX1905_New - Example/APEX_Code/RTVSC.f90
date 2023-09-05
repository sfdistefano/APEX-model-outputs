      SUBROUTINE RTVSC
!     APEX1905
!     THIS SUBPROGRAM ROUTES FLOW THROUGH REACHES USING THE VARIABLE
!     STORAGE COEF METHOD
      USE PARM
      DIMENSION QMS(MHP),QMSI(MHP) 
      CMS=RWSA(IDN1)/(DTHY*360.)
      IWH=0
      NIT=0
      IF(HYDV(IDN1)>QTH)IWH=1
      IF(KFL(26)>0.AND.IWH>0)THEN
          WRITE(KW(26),4)IY,MO,KDA,IDN1,IDN2,NBSA(IDN2),IDO,RWSA(IDN1),WSA&
          (IDN2),RFPL(IDN2),HYDV(IDN1),QCAP(IDN2),YSD(NDRV,IDN1),&
          STY(IDN1)
          WRITE(KW(26),7)
      END IF
!     ROUTE THROUGH THE REACH
      QMSI=0.
      QMS=0.
      ADD=0.
      TOT=0.
      TM=0.
      QPK=0.
      TPK=0.
      NBCX=0
	  NBFX=0
	  I=1
	  DO J=1,MHX
	      I=I-1
          DO K=1,NPD
              I=I+1
              QMSI(I)=QHY(K,IDN1,IHX(J))
          END DO
	  END DO
      QI1=QHY(1,IDN1,IHX(1))
      QO1=QHY(1,IDO,IHX(1))
      ADI=0.
      DO L=1,NHY(IDN1)
          IF(QMSI(L)>0.)EXIT
          QMS(L)=0.
          IF(KFL(26)>0.AND.IWH>0)WRITE(KW(26),6)NIT,TM,TOT,TOT,TOT,TOT,&
          TOT,TOT,TOT,QMS(L),TOT,TOT,TOT,TOT,TOT
          TM=TM+DTHY
      END DO
      ZCH=RCHD(IDN2)
      CBW=RCBW(IDN2)
      XL3=RFPL(IDN2)/3.6
      XLT=XL3/DTHY
      XLS=1000.*RFPL(IDN2)*RFPS(IDN2)
      RCHX(IDN2)=SQRT(RCHS(IDN2))/RCHN(IDN2)
      SSS=SQRT(RCSS(IDN2)*RCSS(IDN2)+1.)
      STHY=0.
      SIA=0.
      QO2=0.
      !     ROUTE UNTIL OUTFLOW IS GREATER THAN 0.                                                                                
      DO I=L,NHY(IDN1)
          I1=I-1                                                                    
          QI2=QMSI(I)
          ADI=ADI+QI2
          CALL HQDAV(AI,CBW,QI2,SSS,ZCH,ZI,CHW,FPW,IDN2)
          DD=SQRT((XLS+ZI)/XLS)                                                          
          XFLO=QI2/DD 
          CALL HQDAV(AII,CBW,XFLO,SSS,ZCH,ZII,CHW,FPW,IDN2)
          SMO=2.*ADI-QI2                                                         
          XFLO=SMO-XLT*AII
          !O2 = SUM OF INFLOW - STORAGE                                                                         
          IF(XFLO>0.)EXIT
          QMS(I)=0.
          STHY=STHY+QI2                                                                                                                                        
          ZOO=0.                                                                         
          AOO=0.                                                                         
          V=.5*QMSI(I)/AII                                                           
          TT=XL3/(V+1.E-5)                                                                
          CVSC=MIN(.99,2.*DTHY/(2.*TT+DTHY))
          QI1=QI2
          TM=TM+DTHY
          IF(KFL(26)>0)WRITE(KW(26),6)IT,TM,RFDT(I),QI2,V,TT,CVSC,SIA,STHY,QO2
          !WRITE(KW(26),6)NIT,TM,ZII,ZOO,AII,AOO,V,TT,CVSC,STHY,QI2,QO2
      END DO
      ISM=I
      IIY=0
      DO 
          I1=I-1
          QI2=QMSI(I)
          SIA=.5*(QI1+QI2)+STHY
          G1=(QI1+QI2+QO1)/3.
          GB=MAX(QI1,QI2,QO1)                      
          GL=0.
          CALL HQDAV(AI,CBW,QI2,SSS,ZCH,ZI,CHW,FPW,IDN2)
          IF(IIY==0)THEN
              V=MAX(.1,QI2/AI)
              T1=XL3/V
              IIY=1
          END IF
          DO IT=1,20
              CALL HQDAV(AO,CBW,G1,SSS,ZCH,ZO,CHW,FPW,IDN2)
              XX=MAX(.75,(XLS+ZI-ZO)/XLS)
              DD=SQRT(XX)      
              V=MAX(.05,(QI2+G1)*DD/(AI+AO))
              T2=XL3/V
              TT=.5*(T1+T2)
              C=MIN(.99,2.*DTHY/(2.*TT+DTHY))
              !Q1=C*SIA
              Q1=C*(SIA-QO1/2.)
              GQ=Q1-G1                                                                       
                  IF(ABS(GQ/G1)<.001)EXIT
                  IF(GQ>0.)THEN
                      GL=G1
                  ELSE
                      GB=G1
                  END IF
                  G1=.5*(GB+GL)
                  IF(GB-GL<.001)EXIT                                                         
          END DO
          IF(IT>20)THEN
              IF(KFL(26)>0)WRITE(KW(26),47)QI1,QI2,QO1,DD,G1,Q1,GB,GL
          END IF
          NIT=NIT+IT
          QO2=MAX(.01,G1)                                                  
          QMS(I)=QO2 
          STHY=SIA-.5*(QO2+QO1)
          T1=T2
          QAV=.25*(QI1+QI2+QO2+QO1)
	      IF(QAV>0..AND.IHY==0)CALL RTYNP(QAV,ADD,TOT,I)
	      DF=TOT+TDEG-TDEP-ADD
	      IF(KFL(26)>0.AND.IWH>0)WRITE(KW(26),6)IT,TM,RFDT(I),QI2,&
          V,TT,C,SIA,STHY,QO2,TOT,TDEG,TDEP,ADD,DF
          TM=TM+DTHY
          !IF(STHY<1.E-10)EXIT
          IF(I>NHY(IDN1))THEN
              QMSI(I)=.9*QMSI(I1)                                              
              IF(QMS(I1)<=1.)EXIT
          END IF   
          I=I+1
          QI1=QI2
          QO1=QO2
      END DO
      NHY(IDO)=I-1
      SUM=0.
      TM=0.
      I=0
      DO K=1,NPD
          I=I+1
          QHY(K,IDO,IHX(1))=QMS(I)
          SUM=SUM+QMS(I)
          TM=TM+DTHY
          IF(QMS(I)>QPK)THEN
              QPK=QMS(I)
              TPK=TM
          END IF
      END DO
      SUM=(SUM-.5*(QHY(1,IDO,IHX(1))+QHY(NPD,IDO,IHX(1))))/CMS
      HYDV(IDO)=SUM
	  RQRB(IDO)=QPK/CMS
	  TC(IDO)=TPK
	  TCAV(IDO)=TCAV(IDO)+TPK
	  !PRAV(IDO)=PRAV(IDO)+RQRB(IDO)
      !IF(RQRB(IDO)>PRB(IDO))PRB(IDO)=RQRB(IDO)
      !NQRB(IDO)=NQRB(IDO)+1
      IF(TC(IDO)>TCMX(IDO))THEN
          TCMX(IDO)=TC(IDO)
      ELSE
          IF(TC(IDO)<TCMN(IDO))TCMN(IDO)=TC(IDO)
      END IF
      DO J=2,MHX
          I=I-1
          DO K=1,NPD
              I=I+1
              QHY(K,IDO,IHX(J))=QMS(I)
          END DO
	  END DO
      XO=I-ISM
      XIT=NIT
      XO=XIT/XO
      IF(NBCX>0)NBCF(IDN2)=NBCF(IDN2)+1
	  IF(NBFX>0)NBFF(IDN2)=NBFF(IDN2)+1
	  NBCT(IDN2)=NBCT(IDN2)+NBCX
	  NBFT(IDN2)=NBFT(IDN2)+NBFX
      II=I+1
      IF(KFL(26)>0.AND.IWH>0)WRITE(KW(26),8)IDN1,IDN2,NBSA(IDN2),SUM,&
      QPK,TPK,XO,TOT,TDEG,TDEP,ADD
      RETURN
    4 FORMAT(//T10,'ROUTE'/T10,'DATE(Y-M-D)=',3I4,2X,'IDN1= ',I8,2X,&
      'ISA= ',I8,2X,'IDSA= ',I8,2X,'IDO= ',I8/T10,'WSAI= ',F12.3,' ha',&
      2X,'WSA= ',F12.3,' ha',2X,'RFPL= ',F12.3,' km'/T10,'HYD VOL= ',F8.3,&
      ' mm',2X,'QCAP= ',F8.3,' m3/s',2X,'YI= 'F8.3,' t/ha'/T10,'STY= ',&
      F8.3,' t/ha')
    6 FORMAT(6X,I4,20F10.3)
    7 FORMAT(//8X,'IT',5X,'Th',8X,'RFmm',5X,'QIm3/s',4X,'Vm/s',7X,'TTh',&
      7X,'C',5X,'QI+STm3/s',4X,'STm3/s',4X,'QOm3/s',4X,'YIt/s',5X,&
      'DEGt/s',4X'DEPt/s',4X,'YOt/s',5X,'DF')
    8 FORMAT(T10,'IDN1= ',I8,2X,'ISA= ',I8,2X,'IDSA= ',I8,2X,'HYD VOL= ',&
      F7.3,' mm',2X,'PEAK RATE= ',F7.2,' m3/s',2X,'TP= ',F7.2,' h',2X,&
      'AVE IT= ',F6.3/T10,'YI= ',F8.3,' t/ha',2X,'DEG= ',F8.3,' t/ha',&
      2X,'DEP= ',F8.3,' t/ha',2X,'YO= ',F8.3,'t/ha')
   15 FORMAT(20X,7F7.2)
   47 FORMAT(// T2,'CONVERGENCE PROBLEM--QI1=',F10.4,' QI2=',F10.4,&
      ' QO1=',F10.4,' QO2=',F10.4,'  G=',F10.4,' G1=',F10.4//)                            
      END