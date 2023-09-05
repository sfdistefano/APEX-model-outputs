      SUBROUTINE NYPA
!     APEX1905
!     THIS SUBPROGRAM PREDICTS DAILY P & K LOSS, GIVEN SOIL LOSS AND
!     ENRICHMENT RATIO.
      USE PARM
      LD1=LID(1,ISA)
      QPU(IDO)=0.
      WSAX=WSA(ISA)
      QMM=.1*QVOL(IDO)/WSAX
      X2=WPML(LD1,ISA)
      !AD1=X2+WPMA(LD1,ISA)+WPO(LD1,ISA)
      YP(IDO)=YEW*WPO(LD1,ISA)
      YPWN(IDO)=YEWN*WPO(LD1,ISA)
      X4=MIN(WPOU(LD1,ISA),YMNU(IDO)*WPOU(LD1,ISA)/(RSDM(LD1,ISA)+1.E-10))
      YPOU(IDO)=X4
      WPOU(LD1,ISA)=WPOU(LD1,ISA)-X4
      WPO(LD1,ISA)=WPO(LD1,ISA)-YP(IDO)-YPWN(IDO)
      YP(IDO)=YP(IDO)+X4
      QAPY=YEW*X2
      APWN=YEWN*X2
      X2=X2-QAPY-APWN
      X1=EXCK(LD1,ISA)*YEW
      EXCK(LD1,ISA)=EXCK(LD1,ISA)-X1
      X3=FIXK(LD1,ISA)*YEW
      FIXK(LD1,ISA)=FIXK(LD1,ISA)-X3
      VAR(149,ISA)=X1+X3
      SMM(149,MO,ISA)=SMM(149,MO,ISA)+X1+X3
      VH=QSF(LD1,ISA)+CPFH(LD1,ISA)
      V=QMM+PKRZ(LD1)+SSF(LD1,ISA)+VH
      IF(V>0.)THEN
          V5=5.*V
          X1=MAX(V5,WT(LD1,ISA)*PRMT(8))
          XZ=MAX(V5,WT(LD1,ISA))
          VPU(ISA)=MIN(.5*WPMU(LD1,ISA),WPMU(LD1,ISA)*PKRZ(LD1)/XZ)
	      WPMU(LD1,ISA)=MAX(0.,WPMU(LD1,ISA)-VPU(ISA))
	      IF(QMM>0.)THEN
	          ! SOLUBLE P RUNOFF
	          XX=QMM/X1
	          IF(LBP==0)THEN
	              ! GLEAMS LINEAR EQ
                  QP(IDO)=X2*XX
              ELSE
                  ! LANGMUIR EQ SOLUTION
                  QQ=MIN(.1*QMM,5.)
                  AD1=0.
                  AD2=0.
                  QT=0.
                  IND=0
                  X5=X2
                  !WRITE(KW(1),'(T10,A,3I4,A,F8.3)')'Y-M-D',IY,MO,KDA,' QVOL=',QVOL(IDO)
                  DO 
                      AD1=AD1+QQ
                      IF(AD1>QMM)THEN
                          AD1=AD1-QQ
                          QQ=QMM-AD1
                          AD1=QMM
                          IND=1
                      END IF    
                      CS=1000.*X5/WT(LD1,ISA)
                      X1=MAX(1.,CPMX(ISA)-CS)
	                  CL=10.*CS/(PRMT(8)*X1)    
	                  QT=.01*CL*QQ
	                  AD2=AD2+QT
	                  !WRITE(KW(1),2)X5,CS,CL,QT,AD1,AD2
	                  IF(IND>0)EXIT
                      X5=X5-QT
	              END DO    
                  QP(IDO)=AD2
	          END IF
	          ! SOLUBLE P RUNOFF FROM MANURE
              QPU(IDO)=WPMU(LD1,ISA)*XX*(1.-URBF(ISA))
              QP(IDO)=QP(IDO)*(1.-URBF(ISA))+URBF(ISA)*QURB(IDO)*.00015
              X2=X2-QP(IDO)
          END IF
          VSS=PKRZ(LD1)+SSF(LD1,ISA)+VH
          IF(VSS>0.)THEN
              X3=MIN(.5*X2,X2*VSS/WT(LD1,ISA))
              WPML(LD1,ISA)=X2-X3
              VAP(ISA)=X3*PKRZ(LD1)/VSS
              QSFP=X3*QSF(LD1,ISA)/VSS
          END IF
          WPMU(LD1,ISA)=MAX(0.,WPMU(LD1,ISA)-QPU(IDO))
          QP(IDO)=QP(IDO)+QPU(IDO)
      END IF 
      QP(IDO)=WSAX*QP(IDO)
      YMP=WPMA(LD1,ISA)*YEW
	  PMWN=YEWN*WPMA(LD1,ISA)
      WPMA(LD1,ISA)=WPMA(LD1,ISA)-YMP-PMWN
      YP(IDO)=WSAX*(YP(IDO)+YMP+QAPY)
      YPWN(IDO)=WSAX*(YPWN(IDO)+PMWN+APWN)
	  TPRK=TPRK+PKRZ(LD1)
      TPSP=TPSP+VAP(ISA)
    !2 FORMAT(T10,' AP=',F8.2,' CS=',F8.3,' CL=',F8.3,' QT=',&
      !F8.3,' ADQ=',F8.3,' ADP=',F8.4)
      !AD2=WPML(LD1,ISA)+WPMA(LD1,ISA)+WPO(LD1,ISA)
      !DF=AD1-QP(IDO)-VAP(ISA)-YP(IDO)-YPWN(IDO)-AD2
      !IF(ABS(DF)>1.E-5)WRITE(KW(1),1)IY,MO,KDA,AD1,QP(IDO),VAP(ISA),YP&
      !(IDO),YPWN(IDO),AD2,DF
    !1 FORMAT(1X,'!!!!!',3I4,10E13.5)  
      RETURN
      END