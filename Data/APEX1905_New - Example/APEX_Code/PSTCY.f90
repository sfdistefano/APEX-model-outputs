      SUBROUTINE PSTCY(PQPX,PYPX,IXP)
      ! APEX1905
      ! THIS SUBPROGRAM SIMULATES PESTICIDE TRANSPORT & DEGRADATION.
      USE PARM 
      LD1=LID(1,ISA)
      Y3=QVOL(IDO)
      Y5=SSF(LD1,ISA)
      Y6=QSF(LD1,ISA)
      Y7=CPFH(LD1,ISA)
      QQ=Y3+Y5+Y6+Y7
      YY=YSD(NDRV,IDO)
      AD1=0.
      DO I=1,NBSL(ISA)
          ISL=LID(I,ISA)
          AD1=AD1+PSTZ(1,ISL,ISA)
      END DO
      SP2=0.
      SP4=0.
      DO K=1,NDP
          ADD=0.
          SUM=0.
          TOT=0.
          X3=0.
          DGF=0.
          WO=0.
          QPST(K,IDO)=0.
          YPST(K,IDO)=0.
          TSPS(K,IDO)=0.
          PSSF(K,LD1,ISA)=0.
	      Y1=PSTZ(K,LD1,ISA)
          Y2=PFOL(K,ISA)
          Y4=PKRZ(LD1)
          IF(IGO(ISA)>0)THEN
              IF(Y2<1.E-5)THEN
	              PFOL(K,ISA)=0.
              ELSE	          
	              IF(RFV(IRF(ISA))>2.54)THEN
                      ! COMPUTE PESTICIDE WASH OFF FROM FOLIAGE
                      WO=PWOF(K)*Y2
                      Y2=Y2-WO
                      Y1=Y1+WO
                  END IF
                  ! COMPUTE PESTICIDE DEGRADATION FROM FOLIAGE
                  DGF=Y2*PHLF(K)
                  PFOL(K,ISA)=Y2-DGF
                  SMMP(6,K,MO,IDO)=SMMP(6,K,MO,IDO)+DGF
                  VARP(6,K,IDO)=DGF
                  VARH(31,ISA)=DGF
              END IF
          END IF
    !     COMPUTE PESTICIDE LOSS FROM TOP SOIL LAYER IN RUNOFF,
    !     LATERAL SUBSURFACE FLOW, & PERCOLATION
          IF(Y1>1.E-5)THEN
              DK=.0001*PKOC(K)*WOC(LD1,ISA)
              X1=PO(LD1,ISA)-S15(LD1,ISA)
              XX=X1+DK
              V=Y3+Y5+Y4+Y6+Y7
              IF(V>0.)THEN
                  VPST=Y1*(1.-EXP(-V/XX))
                  CO=MIN(PSOL(K),VPST/(Y4+PRMT(24)*(Y3+Y5+Y6+Y7)))
                  CS=PRMT(24)*CO
                  X3=CO*Y4
                  QPST(K,IDO)=CS*Y3
                  SMMP(2,K,MO,IDO)=SMMP(2,K,MO,IDO)+QPST(K,IDO)
                  VARP(2,K,IDO)=QPST(K,IDO)
                  SUM=CS*(Y5+Y6+Y7)
                  PSSF(K,1,ISA)=SUM
	              Y1=Y1-X3-QPST(K,IDO)-SUM
                  ! COMPUTE PESTICIDE LOSS WITH SEDIMENT
                  IF(YEW>0.)THEN
                      CS=DK*Y1/XX
                      YPST(K,IDO)=YEW*CS
                      SMMP(5,K,MO,IDO)=SMMP(5,K,MO,IDO)+YPST(K,IDO)
                      VARP(5,K,IDO)=YPST(K,IDO)
                      Y1=Y1-YPST(K,IDO)
                  END IF
              END IF
              ! COMPUTE PESTICIDE DEGRADATION IN TOP SOIL LAYER
              DGS=Y1*PHLS(K)
              Y1=Y1-DGS
              TOT=DGS
              ADD=Y1
          ELSE
              Y1=0.
          END IF
          PSTZ(K,LD1,ISA)=Y1
          IF(IPTS(ISA)>0)THEN
              IF(KDP(K)==KPSN(IPSO(ISA)))THEN
		          QPST(K,IDO)=QPST(K,IDO)+PQPX
		          YPST(K,IDO)=YPST(K,IDO)+PYPX
		      END IF
	      END IF
    !     COMPUTE PESTICIDE MOVEMENT THRU SOIL LAYERS BY LATERAL
    !     SUBSURFACE FLOW & PERCOLATION
          X2=0.
          DO L1=2,NBSL(ISA)
              ISL=LID(L1,ISA)
              Y1=PSTZ(K,ISL,ISA)
              Y1=Y1+X3
              X3=0.
              PSSF(K,ISL,ISA)=0.
              IF(Y1>.01)THEN
                  VH=SSF(ISL,ISA)+QSF(ISL,ISA)+CPFH(ISL,ISA)
                  V=PKRZ(ISL)+VH
                  IF(V>0.)THEN
                      VPST=Y1*(1.-EXP(-V/(PO(ISL,ISA)-S15(ISL,ISA)+.0001*PKOC(K)&
                      *WOC(ISL,ISA))))
                      CO=MIN(PSOL(K),VPST/(PKRZ(ISL)+PRMT(24)*VH))
                      CS=PRMT(24)*CO
                      X5=CS*VH
                      PSSF(K,ISL,ISA)=X5
                      IF(ISL==IDR(ISA))THEN
                          SMMP(10,K,MO,IDO)=SMMP(10,K,MO,IDO)+X5
                          VARP(10,K,IDO)=X5
                      END IF
                      SUM=SUM+X5
                      X3=CO*PKRZ(ISL)
                      !	WRITE(KW(1),32)CS,CO,V,X5,X3
                      IF(ISL==LID(NBSL(ISA),ISA))X2=X3
                      Y1=Y1-X5-X3
                  END IF
                  ! COMPUTE PESTICIDE DEGRADATION IN SOIL LAYERS
                  DGS=Y1*PHLS(K)
                  Y1=Y1-DGS
                  TOT=TOT+DGS
              END IF
              ADD=ADD+Y1
              PSTZ(K,ISL,ISA)=Y1
          END DO 
          SMMP(3,K,MO,IDO)=SMMP(3,K,MO,IDO)+X2
          VARP(3,K,IDO)=X2
          GWPS(K,ISA)=GWPS(K,ISA)+X2
          SMMP(4,K,MO,IDO)=SMMP(4,K,MO,IDO)+SUM
          VARP(4,K,IDO)=SUM
          SMMP(7,K,MO,IDO)=SMMP(7,K,MO,IDO)+TOT
          VARP(7,K,IDO)=TOT
          SMMP(8,K,MO,IDO)=PFOL(K,ISA)
          VARP(8,K,IDO)=PFOL(K,ISA) 
          SMMP(9,K,MO,IDO)=ADD
          VARP(9,K,IDO)=ADD
          PLCH(K)=X2
          SSPS(K)=SUM
          SP2=SP2+QPST(K,IDO)
	      SP4=SP4+YPST(K,IDO)
	      WSAX=WSA(ISA)
          TSPS(K,IDO)=SUM
          VARH(31,IDO)=WSAX*(TOT+DGF)
          SMMH(31,MO,IDO)=SMMH(31,MO,IDO)+VARH(31,IDO)
          PQST=QPST(K,IDO)+SUM
      END DO
      ! CALL PSTFRQ(QQ,YY,SSPS,ISA,IXP)
      VARH(27,IDO)=WSAX*SP2
	  SMMH(27,MO,IDO)=SMMH(27,MO,IDO)+VARH(27,IDO)
	  VARH(29,IDO)=WSAX*SP4	
	  SMMH(29,MO,IDO)=SMMH(29,MO,IDO)+VARH(29,IDO)
	  AD2=0.
      DO I=1,NBSL(ISA)
          ISL=LID(I,ISA)
          AD2=AD2+PSTZ(1,ISL,ISA)
      END DO
      !DF=AD1+WO-TOT-SSPS(1)-QPST(1,IDO)-PLCH(1)-YPST(1,IDO)-AD2
      !IF(ABS(DF/(AD1+1.))>.001)WRITE(KW(1),20)ISA,IY,MO,KDA,AD1,WO,QPST(1,IDO),&
      !SSPS(1),PLCH(1),YPST(1,IDO),AD2,DF
      RETURN
      !20 FORMAT(1X,'#####',4I4,20E13.5)           
      END