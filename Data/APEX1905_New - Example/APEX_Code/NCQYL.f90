      SUBROUTINE NCQYL
!     APEX1905
!     THIS SUBPROGRAM PREDICTS DAILY C LOSS, GIVEN SOIL LOSS AND
!     ENRICHMENT RATIO.
      USE PARM
      LD1=LID(1,ISA)
      WSAX=WSA(ISA)
      Y1=WBMC(LD1,ISA)
      YY=YSD(NDRV,ISA)
      Y4=PKRZ(LD1)
      V=0.
      VBC=0.
      QC(IDO)=0.
      YBC=0.
      QMM=.1*QVOL(IDO)/WSAX
      TOT=WSAX*(WHPC(LD1,ISA)+WHSC(LD1,ISA)+WLMC(LD1,ISA)+WLSC(LD1,ISA))
      X1=MAX(1.E-5,1.-YEW-YEWN)
      YC(IDO)=YEW*TOT
      YCWN(IDO)=YEWN*TOT
      WHSC(LD1,ISA)=WHSC(LD1,ISA)*X1
      WHPC(LD1,ISA)=WHPC(LD1,ISA)*X1
      WLS(LD1,ISA)=WLS(LD1,ISA)*X1
      WLM(LD1,ISA)=WLM(LD1,ISA)*X1
      WLSL(LD1,ISA)=WLSL(LD1,ISA)*X1
      WLSC(LD1,ISA)=WLSC(LD1,ISA)*X1
      WLMC(LD1,ISA)=WLMC(LD1,ISA)*X1
      WLSLC(LD1,ISA)=WLSLC(LD1,ISA)*X1
      WLSLNC(LD1,ISA)=WLSC(LD1,ISA)-WLSLC(LD1,ISA)
      IF(Y1>.01)THEN
          DK=.0001*PRMT(21)*WOC(LD1,ISA)
          X1=PO(LD1,ISA)-S15(LD1,ISA)
          XX=X1+DK
          V=QMM+Y4
          X3=0.
          IF(V>0.)THEN
              X3=Y1*MIN(.5,(1.-EXP(-V/XX)))
              CO=X3/(Y4+PRMT(16)*QMM)
              CS=PRMT(16)*CO
              VBC=CO*Y4
              Y1=Y1-X3
              QC(IDO)=CS*QMM
          END IF
          ! COMPUTE WBMC LOSS WITH SEDIMENT
          IF(YEW>0.)THEN
              CS=DK*Y1/XX
              YBC=MIN(YEW*CS,.5*Y1)
          END IF
      END IF
      WBMC(LD1,ISA)=Y1-YBC      
      DO L=2,NBSL(ISA)
          ISL=LID(L,ISA)
          Y1=WBMC(ISL,ISA)+VBC
          WCMU(ISL,ISA)=WCMU(ISL,ISA)+VBC
          VBC=0.
          IF(Y1>.01)THEN
              V=PKRZ(ISL)
              IF(V>0.)THEN
                  XX=MIN(.75,PKRZ(ISL)/WT(ISL,ISA))
                  VBC=XX*Y1
              END IF
          END IF
          WBMC(ISL,ISA)=Y1-VBC
          WCMU(ISL,ISA)=WCMU(ISL,ISA)-VBC
      END DO
      VBC=VBC*WSAX
      SMM(75,MO,ISA)=SMM(75,MO,ISA)+VBC
      VAR(75,ISA)=VBC
      QC(IDO)=QC(IDO)*WSAX
      SMM(76,MO,ISA)=SMM(76,MO,ISA)+QC(IDO)
      VAR(76,ISA)=QC(IDO)
      SQC(IDO)=SQC(IDO)+QC(IDO)
      VQC(IDO)=QC(IDO)
      X4=MIN(.5*WHSC(LD1,ISA),.5*WCOU(LD1,ISA),YMNU(IDO)*WCOU(LD1,ISA)/(WSAX*RSDM(LD1,ISA)+1.E-10))
      WCOU(LD1,ISA)=WCOU(LD1,ISA)-X4
      WHSC(LD1,ISA)=WHSC(LD1,ISA)-X4
      X4=WSAX*X4
      YCOU(IDO)=X4
      YC(IDO)=YC(IDO)+X4
      SMM(77,MO,ISA)=SMM(77,MO,ISA)+YC(IDO)
      VAR(77,ISA)=YC(IDO)
      VYC(IDO)=YC(IDO)
      SMM(136,MO,ISA)=SMM(136,MO,ISA)+YCWN(IDO)
      VAR(136,ISA)=YCWN(IDO)
      RETURN
      END