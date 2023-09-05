      SUBROUTINE AEXINT(WW,SUM)
!     APEX1905
!     THIS SUBPROGRAM INTEGRATES THE MODIFIED EXPONENTIAL EQ.
      X1=1.
      DX=.1
      SUM=0.
      Y1=0.
      DO WHILE(DX>1.E-4)
          XY=0.
          DO WHILE(XY<.1)
              X2=X1-DX
              IF(X2<=0.)EXIT
              Y2=(-LOG(X2))**WW
              XY=(Y2+Y1)*DX
              SUM=SUM+XY
              X1=X2
              Y1=Y2
          END DO
          DX=DX*.5
      END DO
      SUM=SUM*.5
      RETURN
      END