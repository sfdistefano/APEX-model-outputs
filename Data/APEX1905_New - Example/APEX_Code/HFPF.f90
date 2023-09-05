      SUBROUTINE HFPF(QQ,TRT,SC,Q2)
!     THIS SUBPROGRAM CALCULATES FLOODPLAIN INFILTRATION USING AN ITERATIVE
!     TECHNIQUE WITH THE GREEN & AMPT INFILTRATION EQ 
      USE PARM
      QR=QQ/TRT
      ZI=SC*(SCI(IDN2)/QQ+1.)
      ZI1=ZI
      IF(QR<=ZI)THEN
          Q2=0.
          RETURN
      END IF
      QB=TRT*(QR-ZI)
      GB=0.
      GL=QB
      ZI=SC*(SCI(IDN2)+1.)
      QL=MAX(0.,TRT*(QR-ZI))
      DO IR=1,10
          X1=GL-GB
          IF(ABS(X1)<.01)THEN
              Q2=GL
              EXIT
          END IF
          B2=(QL-QB)/X1
          B1=QL-B2*GL
          X2=1.-B2
          IF(ABS(X2)<1.E-5.OR.B1<1.E-5)THEN
              G1=.5*(GL+GB)
          ELSE    
              G1=B1/X2
          END IF    
          FF=QQ-G1
          IF(FF>0.)THEN
              ZI=SC*(SCI(IDN2)/FF+1.)
          ELSE
              Q2=QQ
              EXIT
          END IF
          Q2=TRT*(QR-ZI)
          GQ=Q2-G1
          IF(ABS(GQ/G1)<.001)EXIT
          IF(GQ>0.)THEN
              GL=G1
              QL=Q2
          ELSE
              GB=G1
              QB=Q2
          END IF
      END DO
      IF(Q2>QQ)Q2=TRT*(QR-ZI1)
      RETURN
      END       