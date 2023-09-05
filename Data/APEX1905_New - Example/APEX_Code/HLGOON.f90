      SUBROUTINE HLGOON(JRT)
!     APEX1905
!     THIS SUBPROGRAM ESTIMATES INFLOW, STORAGE, & IRRIGATION FROM
!     LAGOONS. RUNOFF IS ESTIMATED WITH 90 CN.
      USE PARM
      II=IDON(ISA)
      ALGI(II)=0.
      X2=10.*DALG(ISA)
      DP=.1677*VLG(ISA)**.3333
      TW=18.*DP
      SURA=.0001*TW*TW
      EV=6.*EO*SURA
      VLG(ISA)=MAX(.01,VLG(ISA)-EV+COWW(ISA))
      SMM(20,MO,ISA)=SMM(20,MO,ISA)+EV
      VAR(20,ISA)=EV
      XX=COWW(ISA)
      SMM(21,MO,ISA)=SMM(21,MO,ISA)+XX
      VAR(21,ISA)=XX
      X1=RFV(IRF(ISA))-5.64
      IF(X1>0.)THEN
          QLG=X1*X1/(RFV(IRF(ISA))+22.6)
          QLG=10.*(QLG*(DALG(ISA)-SURA)+RFV(IRF(ISA))*SURA)
          VLG(ISA)=VLG(ISA)+QLG
          SMM(22,MO,ISA)=SMM(22,MO,ISA)+QLG
          VAR(22,ISA)=QLG
      END IF
      IF(VLG(ISA)<=VLGN(ISA))THEN
          JRT=1
          RETURN
      END IF
      IF(VLG(ISA)>VLGM(ISA))THEN
          X1=(VLG(ISA)-VLGM(ISA))
          !QVOL(IDO)=QVOL(IDO)+X1/(10.*WSA(ISA))
          QVOL(IDO)=QVOL(IDO)+X1 !m^3, Jaehak Oct 2020. Unit error fixed
          CFNP(ISA)=WTMU(ISA)/VLG(ISA)
          X3=CFNP(ISA)*X1
          SOFU=SOFU+X3
          YMNU(IDO)=YMNU(IDO)+X3
          SMM(23,MO,ISA)=SMM(23,MO,ISA)+X1
          VAR(23,ISA)=X1
          X1=X1/X2
          WRITE(KW(1),4)ISA,NBSA(ISA),IYR,MO,KDA,X1,X3
          VLG(ISA)=VLGM(ISA)
          WTMU(ISA)=WTMU(ISA)-X3
      END IF
      IF(JDA==LPD)GO TO 1
      IF((VLG(ISA)-VLGN(ISA))/(VLGM(ISA)-VLGN(ISA))>.75)GO TO 1
      IF(IPMP(ISA)==0)THEN
          JRT=1
          RETURN
      END IF
!     IF(RZSW(ISA)>=PAW(ISA))GO TO 7
    1 XX=VLGI(ISA)
      IPMP(ISA)=1
      X3=.95*(VLG(ISA)-VLGN(ISA))
      IF(XX>=X3)THEN
          IPMP(ISA)=0
          XX=X3
      END IF
      SMM(52,MO,ISA)=SMM(52,MO,ISA)+XX
      VAR(52,ISA)=XX
      ALGI(II)=XX
      VLG(ISA)=MAX(.01,VLG(ISA)-XX)
      JRT=0
      RETURN
    4 FORMAT(1X,2I8,1X,I4,2I2,2X,'______LAGOON OVERFLOWED VOL = ',E13.5,&
      ' mm',5X,'MANURE = ',E13.5,' T')
      END