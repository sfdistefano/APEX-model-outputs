      SUBROUTINE ROUTE
!     APEX1905
!     THIS SUBPROGRAM CONTROLS ROUTING OPERATIONS.
      USE PARM 
      IDO=IDOT(ICMD)
      IDN1=IDN1T(ICMD)
      IDN2=IDN2T(ICMD)
      IDNB(IDO)=NBSA(IDN2)
      QVOL(IDO)=QVOL(IDN1)
      RSSF(IDO)=RSSF(IDN1)
      SST(IDO)=SST(IDN1)
      VARH(15,IDN1)=SST(IDN1)
      VARH(15,IDO)=SST(IDN1)
      QDR(IDO)=QDR(IDN1)
      QRF(IDO)=QRF(IDN1)
      WYLD(IDO)=WYLD(IDN1)
      IDNF(IDN2)=IDN1
      SMIO(IDN1)=SMIO(IDN1)+WYLD(IDN1)
      VARH(36,IDN1)=WYLD(IDN1)
      VSSN(IDN1)=TSFN(IDN1)
      SSN(IDN1)=SSN(IDN1)+TSFN(IDN1)
      RWSA(IDO)=RWSA(IDN1)
      TC(IDO)=TC(IDN1)
      RSFN(IDO)=RSFN(IDN1)
      QDRN(IDO)=QDRN(IDN1)
      QRFN(IDO)=QRFN(IDN1)
      RQRB(IDO)=RQRB(IDN1)
      DO I=1,NSZ
          PCT(I,IDO)=PCT(I,IDN1)
      END DO
      DO K=1,NDP
          QPST(K,IDO)=QPST(K,IDN1)
          TSPS(K,IDO)=TSPS(K,IDN1)
	      YPST(K,IDO)=YPST(K,IDN1)
      END DO
      X2=RWSA(IDN1)/WSA(IDN2)
      !X6=SST(IDN1)*X2
      X5=TSFN(IDN1)
      II=IDOA(IDN2)
      SST(II)=SST(IDN1)
      TSFN(IDO)=X5
      TSFN(II)=X5
      SSIN(IDN2)=SSIN(IDN2)+X5
      X5=X5/WSA(IDN2)
      X2=X5/Z(LID(NBSL(IDN2),IDN2),IDN2)
      !SSFI(IDN2)=SSFI(IDN2)+X6
      !SMM(95,MO,IDN2)=SMM(95,MO,IDN2)+X6
      !VAR(95,IDN2)=X6
      X6=.1*SST(IDN1)/WSA(IDN2)
      XX=X6/Z(LID(NBSL(IDN2),IDN2),IDN2)
      Z1=0.
      DO I=1,NBSL(IDN2)
          ISL=LID(I,IDN2)
          X1=Z(ISL,IDN2)-Z1
          X3=XX*X1
          X4=X2*X1
          WNO3(ISL,IDN2)=WNO3(ISL,IDN2)+X4
          SWST(ISL,IDN2)=SWST(ISL,IDN2)+X3
          Z1=Z(ISL,IDN2)
      END DO
      YSD(NDRV,IDO)=0.
      YN(IDO)=0.
      YP(IDO)=0.
      YC(IDO)=0.
      QC(IDO)=0.
      QN(IDO)=0.
      QP(IDO)=0.
	  QPU(IDO)=0.
      YMNU(IDO)=0.
      YCOU(IDO)=0.
      YNOU(IDO)=0.
      YPOU(IDO)=0.
      CALL EYCC(IDN2)
      IF(IHY==0.AND.WYLD(IDO)>0.)THEN
          CALL RTSED
      ELSE    
          IF(NHY(IDN1)>0.)THEN
              SELECT CASE(IHY)
                  CASE(1)
                      CALL RTVSC
                  CASE(2)
                      CALL RTSVS
                  CASE(3)
                      CALL RTM_CVC
                  CASE(4)
                      CALL RTM_CVC4
                  CASE DEFAULT
                      CALL RTVSC
              END SELECT                                      
          END IF
      END IF    
      RETURN
      END