      SUBROUTINE HPURK
!     APEX1905
!     THIS SUBPROGRAM IS THE MASTER PERCOLATION COMPONENT.  IT MANAGES
!     THE ROUTING PROCESS
      USE PARM
      ADD=0.
      SUM=0.
      TOT=0.
      XX=0.
      CPVV=0.
      QDR(IDO)=0.
      QMM=.1*QVOL(IDO)/WSA(ISA)
      VAR(98,ISA)=FPF(ISA)
      AD1=SWLT(ISA)
      DO K=1,NBSL(ISA)
          ISL=LID(K,ISA)
          AD1=AD1+SWST(ISL,ISA)
      END DO
      SEP=MAX(0.,(10.*SALA(ISA)*(RFV(IRF(ISA))-QMM)+RSPK(ISA))/(10.*WSA(ISA)))
      SP1=SEP
      LD1=LID(1,ISA)
	  STLT(ISA)=PRMT(51)*RSD(LD1,ISA)
      SWLT(ISA)=SWLT(ISA)+SEP
      IF(SWLT(ISA)>STLT(ISA))THEN
          SEP=SWLT(ISA)-STLT(ISA)
          SWLT(ISA)=STLT(ISA)
      ELSE
          SEP=0.
      END IF
      !Percolation while flooded/dry Jaehak paddy 2016
	  IF(PADDY_STO(1,ISA)>0)THEN
          SEP=.2*PRMT(39)*PADDY_STO(1,ISA) !Total head is referenced at 120mm ponding depth
          IF(SEP>PADDY_STO(1,ISA))SEP=PADDY_STO(1,ISA) 
      END IF
      IF(IPRK==2)THEN
          CALL HPERC2(3)
      ELSE    
          DO KK=1,NBSL(ISA)
              ISL=LID(KK,ISA)
              DZ=Z(ISL,ISA)-XX
              XX=Z(ISL,ISA)
              SWST(ISL,ISA)=SWST(ISL,ISA)+SEP
              IF(WTBL(ISA)<=Z(ISL,ISA))THEN
                  SSF(ISL,ISA)=0.
                  PKRZ(ISL)=0.
                  CPFH(ISL,ISA)=0.
                  SEP=0.
              ELSE    
                  CPVV=SEP*CPRV(ISL,ISA)
                  X1=SEP-CPVV
                  CPVH(IDO)=X1*CPRH(ISL,ISA)
                  SWST(ISL,ISA)=MAX(1.E-5,SWST(ISL,ISA)-CPVV-CPVH(IDO))
	              !IF(RSAE(ISA)>0.)THEN
	                  !SATX=MAX(1.E-10,RSHC(ISA))
	              !ELSE
	                  SATX=SATC(ISL,ISA)
                  !END IF
                  IF(IPRK==0)THEN
                      CALL HPERC(DZ,SATX)
                  ELSE    
                      CALL HPERC1(DZ,SATX)
                  END IF    
                  SWST(ISL,ISA)=MAX(1.E-5,SWST(ISL,ISA)-SEP-SST(IDO)-QRF(IDO))
                  IF(ISL/=IDR(ISA))THEN
                      SUM=SUM+QRF(IDO)
                  ELSE
                      QDR(IDO)=QRF(IDO)
                  END IF
                  ADD=ADD+SST(IDO)
                  TOT=TOT+CPVH(IDO)
                  SSF(ISL,ISA)=SST(IDO)
                  QSF(ISL,ISA)=QRF(IDO)
                  CPFH(ISL,ISA)=CPVH(IDO)
                  SEP=SEP+CPVV
                  PKRZ(ISL)=SEP
              END IF    
          END DO
          SST(IDO)=ADD
          QRF(IDO)=SUM
          CPVH(IDO)=TOT
          L1=LD1
          DO K=NBSL(ISA),2,-1
              ISL=LID(K,ISA)
              L1=LID(K-1,ISA)
              XX=SWST(ISL,ISA)-PO(ISL,ISA)
              IF(XX>0.)THEN
                  SWST(L1,ISA)=SWST(L1,ISA)+XX
                  PKRZ(L1)=MAX(0.,PKRZ(L1)-XX)
                  SWST(ISL,ISA)=PO(ISL,ISA)
              END IF
              XX=SWST(ISL,ISA)-FC(ISL,ISA)
              IF(XX>0.)THEN
                  X1=VGN(L1,ISA)/(VGN(L1,ISA)-1.)
                  X2=SWST(L1,ISA)-S15(L1,ISA)
                  IF(X2>0.)THEN
                      RTO=(PO(L1,ISA)-S15(L1,ISA))/X2
                      X3=RTO**X1-1.
                      IF(X3>0.)THEN
                          T1=(X3/VGA(L1,ISA)**VGN(L1,ISA))**(1./VGN(L1,ISA))
                          T1=T1/10.19
                      ELSE
                          T1=1.
                      END IF
                  ELSE
                      T1=1500.
                  END IF    
                  ZH=10.*(Z(ISL,ISA)-Z(L1,ISA))
	              X1=VGN(ISL,ISA)/(VGN(ISL,ISA)-1.)
                  X2=SWST(ISL,ISA)-S15(ISL,ISA)
                  IF(X2>0.)THEN
                      RTO=(PO(ISL,ISA)-S15(ISL,ISA))/X2
                      X3=RTO**X1-1.
                      IF(X3>0.)THEN
                          T2=(X3/VGA(ISL,ISA)**VGN(ISL,ISA))**(1./VGN(ISL,ISA))
                          T2=T2/10.19
                      ELSE
                          T2=1.
                      END IF
                  ELSE
                      T2=1500.
                  END IF    
                  T2=T2+ZH
                  IF(T1<T2)CYCLE
                  X1=XX*MIN(PRMT(61),(T1-T2)/(T1+T2),PKRZ(L1))
                  SWST(L1,ISA)=SWST(L1,ISA)+X1
                  PKRZ(L1)=PKRZ(L1)-X1
                  SWST(ISL,ISA)=SWST(ISL,ISA)-X1
              END IF    
          END DO
      END IF
      FPF(ISA)=0.
      AD2=SWLT(ISA)
      DO K=1,NBSL(ISA)
          ISL=LID(K,ISA)
          AD2=AD2+SWST(ISL,ISA)
      END DO
      !DF=AD1+SP1-SST(IDO)-QRF(IDO)-PKRZ(LID(NBSL(ISA),ISA))-AD2
      !IF(ABS(DF)>.001)WRITE(KW(1),1)IY,MO,KDA,AD1,AD2,SP1,SST(IDO),QRF(IDO),&
      !PKRZ(LID(NBSL(ISA),ISA)),DF
      XX=10.*WSA(ISA)
      SST(IDO)=SST(IDO)*XX
      QRF(IDO)=QRF(IDO)*XX
      CPVH(IDO)=CPVH(IDO)*XX
      QDR(IDO)=QDR(IDO)*XX
      SMM(17,MO,ISA)=SMM(17,MO,ISA)+QDR(IDO)
      VAR(17,ISA)=QDR(IDO)
    !1 FORMAT(5X,'PURK',3I4,10E13.5)  
      RETURN
      END

