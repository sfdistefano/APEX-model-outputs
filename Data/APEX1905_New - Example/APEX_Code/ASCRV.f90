      SUBROUTINE ASCRV(X1,X2,X3,X4,SCLM,I,KWX)
!     APEX1905
!     THIS SUBPROGRAM COMPUTES S CURVE PARMS GIVEN 2(X,Y)POINTS.
      DIMENSION SCLM(38)
      XX=LOG(X3/X1-X3)
      X2=(XX-LOG(X4/X2-X4))/(X4-X3)
      X1=XX+X3*X2
      X=X4
      IF(X2<0.)THEN
          DO IT=1,20
              Y1=X1*EXP(-X2*X)
              FU=X+Y1-X*(1.-X2*Y1)
              IF(ABS(FU)<.001)EXIT
              DFDX=-X2*X2*X*Y1
              DF=FU/DFDX
              XX=ABS(DF)
              X5=.5*X
              IF(XX>X5)DF=X5*XX/DF
              X=X-DF
          END DO
          IF(IT>20)WRITE(KWX,'(5X,A)')'S-CURVE LIMIT NOT CONVERGE'
          SCLM(I)=X
          WRITE(KWX,'(5X,A,I4,4E12.5)')'!!!!!',IT,X1,X2,X,FU
      ELSE
          SCLM(I)=0.
      END IF    
      RETURN
      END