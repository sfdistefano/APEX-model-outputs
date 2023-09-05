      SUBROUTINE HQDAV(A,CBW,QQ,SSS,ZCH,ZX,CHW,FPW,IDY)
!     APEX1905
!     THIS SUBPROGRAM COMPUTES FLOW AREA AND DEPTH GIVEN RATE
      USE PARM
      ZX=.5*ZCH 
      DO IT=1,10
          IF(QQ>QCAP(IDY))THEN
              ZX=MAX(ZX,ZCH)
              ZZ=ZX-ZCH
              !COMPUTE CH FLOW ABOVE QCAP
              ACH=CHXA(IDY)+ZZ*RCTW(IDY)
              R=ACH/CHXP(IDY)
              QCH=ACH*R**.66667*RCHX(IDY)
              CHW=RCTW(IDY)
              !COMPUTE FP FLOW
              AFPX=ZZ*(RFPW(IDY)-RCTW(IDY))
              QFP=AFPX*ZZ**.66667*RFPX(IDY)/RFPW(IDY)
              Q=QCH+QFP
              A=ACH+AFPX
              FPW=RFPW(IDY)
              NBCF=1
          ELSE
              X1=ZX*RCSS(IDY)
              A=ZX*(CBW+X1)
              P=CBW+2.*SSS*ZX
              Q=A**1.66667*RCHX(IDY)/P**.66667
              CHW=CBW+2.*X1
              FPW=0.
              NBCX=1     
          END IF
          FU=Q-QQ
          !WRITE(KW(26),2)IT,QQ,Q,ZX,FU
          !2 FORMAT(1X,I4,4E16.6)
          X6=MAX(1.,QQ)
          IF(ABS(FU/X6)<.001)EXIT
          IF(IT==1)THEN
              DFQ=-.1*ZX
          ELSE
              DFDZ=(FU-FU1)/(ZX-ZX1)
              DFQ=FU/DFDZ
          END IF
          FU1=FU
          ZX1=ZX
          ZX=ZX-DFQ
      END DO
      IF(KFL(26)>0.AND.IT>20)WRITE(KW(26),1)IT,ZX,FU,QQ,Q
      RETURN
    1 FORMAT(1X,'QDAV XCONV',I4,4E16.6)
      END