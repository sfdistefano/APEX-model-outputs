      SUBROUTINE AINTRIC(X,Y,N1,N2)
      ! APEX1905
      ! THIS SUBPROGRAM INTERPOLATES CONCENTRATIONS FROM LAYERS WITH UN
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
              Y(J,ISA)=TOT+X(L,ISA)*(ZC(J,ISA)-ZZ)
              ZZ=ZC(J,ISA)
              J=J+1
              IF(J>N2)RETURN
              TOT=0.
          END DO
          IF(J>N2)RETURN
          TOT=TOT+X(L,ISA)*(Z(L,ISA)-ZZ)
          Z1=Z(L,ISA) 
          ZZ=Z(L,ISA)
      END DO
      DO I=1,N2-1
          Y(I,ISA)=Y(I,ISA)/DZDN
      END DO 
      Y(N2,ISA)=MAX(X(LID(N1,ISA),ISA),TOT/DZDN)
      RETURN
      END