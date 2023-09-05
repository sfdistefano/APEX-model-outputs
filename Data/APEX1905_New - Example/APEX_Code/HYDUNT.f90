      SUBROUTINE HYDUNT                                                              
!     THIS SUBPROGRAM DEVELOPS A UNIT HYDROGRAPH,CONVERTS MASS RAINFALL              
!     TO POINT RUNOFF,AND COMPUTES STORM HYDROGRAPHS BY SUMMATION.
      USE PARM
      DIMENSION CFS(MHP),QMS(MHP),QVX(MHP)
      CMM=WSA(ISA)/(DTHY*360.)
      IWH=0
      QI=QVOL(IDO)
      IF(QVOL(IDO)>QTH)IWH=1
      XLDW=100.*CHL(ISA)**2/WSA(ISA)                                                                                 
      !XK=.00282*WSA(ISA)**.325*CHS(ISA)**(-.729)*XLDW**.296                                
      !TPU=.0158*WSA(ISA)**.351*CHS(ISA)**(-.353)*XLDW**.45 
      XK=7.
      TPU=10.
      IF(KFL(26)>0.AND.IWH>0)THEN
          WRITE(KW(26),7)IY,MO,KDA,IDO,ISA,NBSA(ISA),WSA(ISA),CN,&
          QVOL(IDO),TC(IDO),XK,TPU
      END IF
      QVX=0.
      QMS=0.
      XX=0.
      BRFX=0.
      TOT=0.
      T1=0. 
      DO I=2,NRF
          T1=T1+DTHY
          IF(INFL==4)THEN
              QVX(I)=QGA(I)
          ELSE
              X1=RFDT(I)-SCN2
              IF(X1>0.)THEN
                  QV=X1*X1/(RFDT(I)+.8*SCN)
              ELSE
                  QV=0.
              END IF
              QVX(I)=QV-XX
              XX=QV
          END IF
          BRFX=MAX(QVX(I),BRFX)
          TOT=TOT+QVX(I)
          QGA(I)=0.
      END DO                               
      BRFX=BRFX*CMM
      ! COMPUTE N BY ITERATION.
      XKTP=XK/TPU                                                                    
      XN=5.                                                                          
      DO I=1,20                                                                    
          XNM=XN-1.                                                                      
          T=1.+SQRT(1./XNM)                                                              
          TN=T**XNM                                                                      
          EE=EXP(-XNM*(T-1.))                                                            
          FU=EE*(XNM*(T**(XN-2.)-TN)+TN/XKTP)                                            
          IF(ABS(FU)<1.E-4)EXIT                                                    
          IF(I>1)THEN
              DFDN=(FU-FU1)/(XN-XN1)                                                         
              DXN=FU/DFDN                                                                    
              ADX=ABS(DXN)                                                                   
              IF(ADX>.5*XNM)DXN=.5*XNM*ADX/DXN                                            
          ELSE                                                             
              DXN=-.1                                                                        
          END IF
          FU1=FU                                                                         
          XN1=XN                                                                         
          XN=XN-DXN
      END DO
      IF(KFL(26)>0.AND.I>20)WRITE(KW(26),'(T2,A)')'N DID NOT CONVERGE'                                                                                                                        
      !     DETERMINE C1. 
      TINF=T
      DELT=T/50.                                                                 
      TC1=0.                                                                         
      SUM=0.
      DO I=2,51 
          TC1=TC1+DELT                                                                   
          X1=(TC1**XNM)*EXP(XNM*(1.-TC1))
          SUM=SUM+X1
      END DO 
      SUM=SUM-X1/2.                                                                
      C1=SUM*DELT
      CFSII=X1                                                                  
      TTINF=TINF*TPU                                                                 
      TREC1=TTINF+2.*XK                                                              
      EEE=EXP((TTINF-TREC1)/XK)                                                      
      XK1=XK                                                                     
      ! COMPUTE B,QPK,AND CFSI.                                                         
      B=360.*(C1+CFSII*(XKTP*(1.-EEE)+EEE*(XK1/TPU)))                                                                                                     
      QPK=WSA(ISA)/(B*TPU)
      CFSI=QPK*X1                                                                
      CFSR1=CFSI*EEE                                                                                                                               
      SUM=0.                                                                         
!     COMPUTE UNIT HYDROGRAPH.                                                       
      T2=0.                                                                          
      CFS(1)=0.
      DO I=2,MHP 
          T2=T2+DTHY
          IF(T2<=TTINF)THEN
              CFS(I)=QPK*((T2/TPU)**XNM)*EXP(XNM *( 1.-T2/TPU))                               
          ELSE                  
              IF(T2<=TREC1)THEN
                  CFS(I)=CFSI*EXP((TTINF-T2)/XK)                                                 
              ELSE
                  CFS(I)=CFSR1*EXP((TREC1-T2)/XK1)                                               
                  IF(CFS(I)<=.0001)EXIT
              END IF
          END IF                                                                              
          SUM=SUM+CFS(I)
      END DO
      ICND=I                                                                         
      SUM=SUM/CMM
      IF(KFL(26)>0)THEN
          WRITE(KW(26),30)XK,TPU,XN,QPK,SUM,TOT,BRFX                              
          WRITE(KW(26),9)
      END IF    
!     COMPUTE STORM HYDROGRAPH.                                                      
      DO J=2,NRF                                                                  
          N=MIN(MHP,J+ICND-2)                                                                     
          I=2                                                                            
          DO K=J,N                                                                    
              QMS(K)=QMS(K)+QVX(J)*CFS(I)                                           
              I=I+1
          END DO
      END DO                                                                          
      NHY(IDO)=N
      QPK=0.
      TPK=0.
      T1=0.
      SUM=0.
      I=0
      DO K=1,NPD
          I=I+1
          CFS(I)=QHY(K,IDO,IHX(1))
          QHY(K,IDO,IHX(1))=QHY(K,IDO,IHX(1))+QMS(I)
          SUM=SUM+QHY(K,IDO,IHX(1))
          IF(QHY(K,IDO,IHX(1))>QPK)THEN
              QPK=QHY(K,IDO,IHX(1))
              TPK=T1
          END IF
          IF(KFL(26)>0.AND.IWH>0)WRITE(KW(26),6)T1,RFDT(I),QVX(I),QMS(I),&
          CFS(I),QHY(K,IDO,IHX(1))
          T1=T1+DTHY
      END DO
      SUM=(SUM-.5*(QHY(1,IDO,IHX(1))+QHY(NPD,IDO,IHX(1))))/CMM
	  HYDV(IDO)=SUM
	  IF(KFL(26)>0.AND.IWH>0)WRITE(KW(26),8)IDO,ISA,NBSA(ISA),QI,SUM,&
      QPK,TPK
      RQRB(IDO)=QPK*360./WSA(ISA)
      !NQRB(IDO)=NQRB(IDO)+1
      !PRAV(IDO)=PRAV(IDO)+RQRB(IDO)
      !PRSD(ISA)=PRSD(ISA)+RQRB(IDO)*RQRB(IDO)
      !IF(RQRB(IDO)>PRB(IDO))PRB(IDO)=RQRB(IDO)
      !X1=RQRB(IDO)/QI
      !IF(X1>QRQB(ISA))QRQB(ISA)=X1
      !QRBQ(ISA)=QRBQ(ISA)+X1
      TCAV(IDO)=TCAV(IDO)+TC(IDO)
      IF(TC(IDO)<TCMX(IDO))THEN
          IF(TC(IDO)<TCMN(IDO))TCMN(IDO)=TC(IDO)
      ELSE
          TCMX(IDO)=TC(IDO)
      END IF
      DO J=2,MHX
          I=I-1
          T1=T1-DTHY
          DO K=1,NPD
              I=I+1
              QHY(K,IDO,IHX(J))=QHY(K,IDO,IHX(J))+QMS(I)
              IF(KFL(26)>0.AND.IWH>0)WRITE(KW(26),6)T1,RFDT(I),QVX(I),QHY(K,IDO,IHX(J))
              T1=T1+DTHY
          END DO
          IF(I>NHY(IDO))EXIT
	  END DO
      IF(SUM<.01)NHY(IDO)=0
            IF(INFL/=4)NRF=0
      RETURN
    6 FORMAT(6X,3F10.2,3E13.5)  
    7 FORMAT(//T10,'COMPUTE HYD'/T10,'DATE(Y-M-D)=',3I4,2X,'IDO= ',I8&
      ,2X,'ISA= ',I8,2X,'IDSA= ',I8/T10,'WSA= ',F12.1,' ha',2X,'CN= ',&
      F5.1,2X,'QVOL= ',F7.3,' mm',2X,'TC= ',F5.2,' h',2X,'XK= ',F7.3,&
      ' h',2X,'TPU= ',F7.3,' h')
    8 FORMAT(5X,'IDO= ',I8,2X,'ISA= ',I8,2X,'IDSA= ',I8,2X,'QVOL=',F8.3,&
      ' mm',2X,'HYD VOL= ',F8.3,' mm'/5X,'PEAK RATE= ',E13.5,' m3/s',2X,&
      'TP= ',F7.2,' h')
    9 FORMAT(//11X,'Th',9X,'RFmm',6X,'QVmm',5X,'DFLOm3/s',6X,'PFLOm3/s',&
      6X,'TFLOm3/s')  
   30 FORMAT(T10,'K=',F7.3,' h'/T10,'TP=',F7.3,' h'/T10,'SHAPE CONSTANT,N=',&
      F6.3/T10,'UNIT PEAK=',F10.3,' m3/s'/T10,'UNIT VOLUME=',F6.3,' mm'/T10,&
      'RF EXCESS'/T15,'VOL=',F7.2,' mm'/T15,'QPK=',F10.1,' m3/s'/)                        
      END          