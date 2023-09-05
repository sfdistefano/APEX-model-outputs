      SUBROUTINE EWNINT
!     APEX1905
!     THIS SUBPROGRAM ESTIMATES DAILY POTENTIAL WIND EROSION FOR A BARE
!     SOIL BY INTEGRATING EROSION RATES WITH TIME GIVEN THE WIND SPEED
!     DISTRIBUTION
      USE PARM
      WW=ATRI(.1,.35,.6,IDG(9))
      SUM=.6424*WW**(-.1508)*EXP(.4336*WW)
      DU10=USTRT/.0408
      Y1=DU10/U10(IRF(ISA))
      X1=EXP(-Y1**(1./WW))
      X11=X1
      DX=.1
      YW(IDO)=0.
      Z1=0.
      DO WHILE(DX>1.E-4)
          X2=X1-DX
          IF(X2>0.)THEN
              Y2=(-LOG(X2))**WW
              Z2=EROWN(Y2)
              XY=(Y1+Y2)*DX
              XZ=(Z2+Z1)*DX
              YW(IDO)=YW(IDO)+XZ
              X1=X2
              Y1=Y2
              Z1=Z2
              IF(XY<.1)CYCLE
          END IF
          DX=DX*.5
      END DO
      YW(IDO)=.5*YW(IDO)*X11/SUM
      RETURN
      END