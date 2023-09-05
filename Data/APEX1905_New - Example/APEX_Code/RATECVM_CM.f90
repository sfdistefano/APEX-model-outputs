      SUBROUTINE RATECVM_CM(CBW,DZRT,SSS,ZCH)
!     APEX1905
!     THIS SUBPROGRAM COMPUTES A RATING CURVE (FLOW RATE AND FLOW AREA 
!     GIVEN DEPTH) AND THE M_C EXPONENT M
      USE PARM
      DIMENSION AX(30)
      ZX=0.
      DO I=1,NPRC
          ZX=ZX+DZRT
          ZZ=ZX-ZCH
          IF(ZZ>0.)THEN
              !COMPUTE CH FLOW ABOVE QCAP
              ACH=CHXA(ISA)+ZZ*RCTW(ISA)
              R=ACH/CHXP(ISA)
              QCH=ACH*R**.66667*RCHX(ISA)
              CHW=RCTW(ISA)
              !COMPUTE FP FLOW
              AFPX=ZZ*(RFPW(ISA)-RCTW(ISA))
              QFP=AFPX*ZZ**.66667*RFPX(ISA)/RFPW(ISA)
              SCFS(I,ISA)=QCH+QFP
              AX(I)=ACH+AFPX
          ELSE
              X1=ZX*RCSS(ISA)
              AX(I)=ZX*(CBW+X1)
              P=CBW+2.*SSS*ZX
              SCFS(I,ISA)=AX(I)**1.66667*RCHX(ISA)/P**.66667
          END IF
      END DO
      XMS(1,ISA)=1.66667
      XMS(2,ISA)=LOG10(SCFS(2,ISA)/SCFS(1,ISA))/LOG10(AX(2)/AX(1))
      TOT=XMS(2,ISA)*SCFS(2,ISA)
      DO J=3,NPRC
          J1=J-1          
          X1=LOG10(SCFS(J,ISA)/SCFS(J1,ISA))/LOG10(AX(J)/AX(J1))
          TOT=TOT+X1*(SCFS(J,ISA)-SCFS(J1,ISA))
          XMS(J,ISA)=TOT/SCFS(J,ISA)
      END DO
      IF(KFL(26)>0)THEN
          ZX=0.
          DO J=1,NPRC
              ZX=ZX+DZRT
              WRITE(KW(26),1)ZX,AX(J),SCFS(J,ISA),XMS(J,ISA)
          END DO                                                            
      END IF    
      RETURN
    1 FORMAT(5X,4F10.3)      
      END