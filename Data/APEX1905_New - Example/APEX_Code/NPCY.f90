      SUBROUTINE NPCY
!     APEX1905
!     THIS SUBPROGRAM IS THE MASTER NUTRIENT CYCLING SUBROUTINE.
!     CALLS NPMIN, NYNIT, NLCH, NCNMI, AND NDNIT FOR EACH SOIL
!     LAYER.
      USE PARM 
      WSAX=WSA(ISA)
      STDX=0.
      DO K=1,LC
          STDX=STDX+STD(K,ISA)
      END DO          
      SGMN=0.
      SN2=0.
      SMP=0.
      SNIT=0.
      TRSP=0.
      TSFS=0.
      WBMX=0.
      IRTO=0
      QDRN=0.
      STRC(ISA)=MAX(1.E-10,STRC(ISA)*PRMT(52))   !should be **, changed by Liwang
      LD1=LID(1,ISA)
      XX=0.
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          RSPC(ISL,ISA)=0.
          RNMN(ISL,ISA)=0.
          WDN=0.
          ZX=Z(ISL,ISA)-XX
          X1=SWST(ISL,ISA)-S15(ISL,ISA)
          IF(X1<0.)THEN
              SUT(ISL,ISA)=.1*(SWST(ISL,ISA)/S15(ISL,ISA))**2
          ELSE
              SUT(ISL,ISA)=.1+.9*SQRT(X1/(FC(ISL,ISA)-S15(ISL,ISA)))
          END IF    
          CALL NPMIN
          CALL NKMIN
          IF(ISL/=LID(1,ISA))THEN
              L1=LID(J-1,ISA) 
              CALL NLCH(L1)
          ELSE
              ZZ=MIN(.5,PRMT(76)*(1.+.1*RFV(IRF(ISA))))
              IF(STDX>.001)CALL NSTDFAL(ZZ)
              CALL NYNIT
          END IF    
          IF(ISL/=IDR(ISA))THEN
              QRFN(IDO)=QRFN(IDO)+QSFN*WSAX
              SMM(84,MO,ISA)=SMM(84,MO,ISA)+QRFN(IDO)
              QRFP(IDO)=QRFP(IDO)+QSFP*WSAX
              SMM(142,MO,ISA)=SMM(142,MO,ISA)+QRFP(IDO)
          ELSE
              QDRN(IDO)=QSFN*WSAX  !kg
              VAR(47,ISA)=QDRN(IDO)
              SMM(47,MO,ISA)=SMM(47,MO,ISA)+QDRN(IDO)  !kg
              QDRP(IDO)=QSFP*WSAX
              VAR(143,ISA)=QDRP(IDO) !kg
              SMM(143,MO,ISA)=SMM(143,MO,ISA)+QDRP(IDO)
          END IF
          TSFN(IDO)=TSFN(IDO)+SSFN
          TSFK(IDO)=TSFK(IDO)+SSFK
          TSFS=TSFS+SSST
          Z5=500.*(Z(ISL,ISA)+XX)
          IF(WNH3(ISL,ISA)>.01)CALL NITVOL(Z5)
          IF(STMP(ISL,ISA)>0..AND.LUNS(ISA)/=35)THEN
              X1=STMP(ISL,ISA)
              IF(SCLM(14)>0.)X1=MIN(X1,SCLM(14))
              CDG(ISL,ISA)=X1/(X1+EXP(SCRP(14,1)-SCRP(14,2)*X1))
              IF(Z(ISL,ISA)>BIG(ISA).OR.STRC(ISA)<10.)THEN
                  X1=1.
              ELSE    
                  X1=STRC(ISA)**PRMT(67)
              END IF    
              X2=Z5
              IF(SCLM(20)>0.)X2=MIN(X2,SCLM(20))
              OX=1.-PRMT(53)*X2/(X2+EXP(SCRP(20,1)-SCRP(20,2)*X2))
              IF(IDNT>2.AND.CGO2(LID(1,ISA),ISA)>0..AND.IOX==2)OX=CGO2(ISL,ISA)&
              /CGO2(LID(1,ISA),ISA)
              X2=SQRT(CDG(ISL,ISA)*SUT(ISL,ISA))*OX
              CS=MIN(10.,X2*PRMT(70)*X1)
	          CALL NRSPC(CS)
              CALL NPMN(CS)
              X3=1.
              IF(SWST(ISL,ISA)>1.1*FC(ISL,ISA))X3=0.
              IF(ISL==LD1)THEN
                  X4=RSD(LD1,ISA)
                  IF(SCLM(27)>0.)X4=MIN(X4,SCLM(27))
                  F=1.+10.*X4/(X4+EXP(SCRP(27,1)-SCRP(27,2)*X4))
              ELSE
                  F=1.
              END IF
              IF(Z(ISL,ISA)<ZMIX)THEN
                  WBMX=WBMX+ZX*X2*X3*F
              ELSE
                  IF(IRTO==0)THEN
                      X1=ZMIX-XX
                      RTO=X1/ZX
                      WBMX=WBMX+X1*X2*RTO*X3
                      IRTO=1
                  END IF    
              END IF    
              SMP=SMP+WMP
              IF(IDNT<3)THEN
	              DZA=1000.*(Z(ISL,ISA)-XX)
	              IF(IDNT==2)THEN
	                  CALL NDNITAK(DZA)
	              ELSE
	                  WDN=0.
	                  DN2=0.
	                  DN2O=0.
	                  CALL NDNIT
	                  SN2=SN2+DN2
	              END IF  
	              SDN=SDN+WDN
	              SN2O=SN2O+DN2O
	          END IF
          END IF
          XX=Z(ISL,ISA)
      END DO
      WBMX=WBMX*PRMT(29)/ZMIX
      GWSN(ISA)=GWSN(ISA)+VNO3(LNS,ISA)*WSAX
      TSFN(IDO)=TSFN(IDO)*WSAX
      VAR(80,ISA)=RSFN(IDO)
      SMM(80,MO,ISA)=SMM(80,MO,ISA)+RSFN(IDO)
      X51=WSAX*(VAP(ISA)+VPU(ISA))
      SMM(51,MO,ISA)=SMM(51,MO,ISA)+X51
      VAR(51,ISA)=X51
      SMM(133,MO,ISA)=SMM(133,MO,ISA)+VSLT(ISA)
      SMM(152,MO,ISA)=SMM(152,MO,ISA)+VSK(ISA) 
      RETURN
      END