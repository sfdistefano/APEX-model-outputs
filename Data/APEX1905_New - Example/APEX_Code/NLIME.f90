      SUBROUTINE NLIME
!     APEX1905
!     THIS SUBPROGRAM APPLIES LIME WHEN THE SUM OF THE SOIL LIME REQUIRE
!     MENT AND ACCUMULATED LIME REQUIREMENT CAUSED BY N FERTILIZER EX-
!     CEED 4 t/ha.
      USE PARM
      OC=0.
      TOT=0.
      XZ=0.
      ZZ=0.
      XX=0.
      DO J=1,NBSL(ISA)
          L=LID(J,ISA)
          IF(Z(L,ISA)>BIG(ISA))EXIT
          XY=WT(L,ISA)
          OC=OC+WOC(L,ISA)*XY
          XZ=XZ+CEC(L,ISA)*XY
          TOT=TOT+PH(L,ISA)*XY
          ZZ=ZZ+SMB(L,ISA)*XY
          XX=XX+XY
      END DO
      IF(J>NBSL(ISA))THEN
          J=NBSL(ISA)
          L=LID(J,ISA)
      ELSE    
          L1=LID(J-1,ISA)
          W3=Z(L,ISA)-Z(L1,ISA)
          W2=BIG(ISA)-Z(L1,ISA)
          RTO=W2*WT(L,ISA)/W3
          OC=OC+RTO*WOC(L,ISA)
          TOT=TOT+RTO*PH(L,ISA)
          ZZ=ZZ+RTO*SMB(L,ISA)
          XZ=XZ+RTO*CEC(L,ISA)
          XX=XX+RTO
      END IF    
      XZ=XZ/XX
      OC=OC/XX
      TOT=TOT/XX
      ZZ=ZZ/XX
      XY=.001*XX
      X1=SMY(54,ISA)+SMY(55,ISA)
      DSB=.036*(SMY(43,ISA)+X1-SMFN(ISA))/XX
      SMFN(ISA)=X1
      BS=100./XZ
      TOT=TOT-.05*DSB*BS
      CALL NLIMA(ZZ,DSB,BS,TOT,ALSX,OC,BSA)
      IF(TOT>6.5)GO TO 6
      IF(IDS(ISA)/=4)GO TO 5
      EAL=.01*ALSX*XZ
      TLA=EAL*XY
      IF(TLA<1.)GO TO 6
      TOT=5.4
      CALL NLIMA(ZZ,-EAL,BS,TOT,ALSX,OC,BSA)
      GO TO 7
    5 DBS=MIN((6.5-TOT)/.023,90.-BSA)
      ALN=0.
      TLA=DBS*XY/BS
      IF(TLA<2..AND.TOT>5.)GO TO 6
      PHN=6.5
      BSA=(BSA+DBS)/BS
      GO TO 8
    6 TLA=0.
    7 ALN=ALSX
      PHN=TOT
      BSA=ZZ
    8 DO K=1,J
          L=LID(K,ISA)
          TOT=SMB(L,ISA)
          XZ=PH(L,ISA)
          ALSX=ALS(L,ISA)
          SMB(L,ISA)=BSA
          PH(L,ISA)=PHN
          ALS(L,ISA)=ALN
      END DO          
      IF(J==NBSL(ISA))RETURN
      L=LID(J,ISA)
      W1=Z(L,ISA)-BIG(ISA)
      SMB(L,ISA)=(W1*TOT+W2*SMB(L,ISA))/W3
      PH(L,ISA)=(W1*XZ+W2*PH(L,ISA))/W3
      ALS(L,ISA)=(W1*ALSX+W2*ALS(L,ISA))/W3
      RETURN
      END