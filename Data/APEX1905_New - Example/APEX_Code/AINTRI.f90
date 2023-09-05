      SUBROUTINE AINTRI(X,Y,N1,N2)
      ! APEX1905
      ! THIS SUBPROGRAM INTERPOLATES SOIL PROPERTIES FROM LAYERS WITH UN
      ! EQUAL THICKNESS TO LAYERS OF EQUAL THICKNESS USED IN DIFFERENTIAL
      ! EQUATIONS OF GAS DIFFUSION.
      USE PARM
      DIMENSION X(MSL,MSA),Y(MSC,MSA)
      ZZ=0.
      Z1=0.
      TOT=0.
      J=1
      DO K=1,N1
          L=LID(K,ISA)
          DO 
              IF(ZC(J,ISA)>Z(L,ISA))EXIT
              Y(J,ISA)=TOT+X(L,ISA)*(ZC(J,ISA)-ZZ)/(Z(L,ISA)-Z1)
              ZZ=ZC(J,ISA)
              J=J+1
              IF(J>N2)EXIT
              TOT=0.
          END DO
          TOT=TOT+X(L,ISA)*(Z(L,ISA)-ZZ)/(Z(L,ISA)-Z1)
          IF(J>N2)EXIT
          Z1=Z(L,ISA) 
          ZZ=Z(L,ISA)
      END DO
      I=MIN(J,N2)
      Y(I,ISA)=TOT
      RETURN
      END