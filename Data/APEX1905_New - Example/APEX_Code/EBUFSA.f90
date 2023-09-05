      SUBROUTINE EBUFSA
!     THIS SUBPROGRAM ROUTES SEDIMENT THROUGH UPLAND BUFFERS WITHIN A 
!     SUBAREA.  USED IN SIMULATING LARGE SUBAREAS WHERE ONLY THE 
!     FRACTION OF THE SA CONTROLLED BY BUFFERS IS KNOWN--NOT THE 
!     LOCATION OF INDIVIDUAL BUFFERS.
      USE PARM
      IDX=IDOA(ISA)
      WSA1=WSA(ISA)*BCOF(ISA)
      Q1=QVOL(IDO)
      Y1=YSD(NDRV,IDO)
      CIN=.1*Y1/Q1
      DEP=0.
      QFP=Q1*WSA1/(360.*(DUR+TC(ISA)))
      BFW=50.*WSA1
      DFP=(QFP/(BFSN(ISA)*BFW))**.6
      VFP=QFP/(DFP*BFW)
      CYFP=PRMT(19)*VFP**PRMT(18)
      TRT=.000278*BFFL(ISA)/VFP
      XX=FPSC(ISA)
      F=MIN(Q1,XX*TRT)
      QO=Q1-F
      IF(QO<1.E-5)THEN
          DEP=Y1
      ELSE
          X7=PRMT(45)*TRT*PSZM(ISA)
          IF(X7<10.)THEN
              F=1.-EXP(-X7)
          ELSE
              F=1.
          END IF
          X5=10.*(CYFP*QO-CIN*Q1)*F
          IF(X5>0.)THEN
              DEG=EK(ISA)*CVF(ISA)*X5
          ELSE
              DEP=-X5
          END IF
      END IF      
      Y2=MAX(1.E-10,Y1-DEP)
      QVOL(IDO)=QO*BCOF(ISA)+Q1*(1.-BCOF(ISA))
      YSD(NDRV,IDO)=Y2*BCOF(ISA)+YSD(NDRV,IDO)*(1.-BCOF(ISA))
      CY=.1*YSD(NDRV,IDO)/QVOL(IDO)
      IF(IERT>0)THEN
          ERTO=MIN(1.25,.78*(CY+1.E-4)**(-.2468))
      ELSE 
          B2=ALOG10(DRTO)/2.301
          B1=1./.1**B2
          ERTO=MIN(1.25,B1*(CY+1.E-4)**B2)
      END IF
      RETURN
      END 