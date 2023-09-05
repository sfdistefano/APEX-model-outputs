      SUBROUTINE AINTRX(X,Y,N1,N2)
      ! APEX1905
      ! THIS SUBPROGRAM INTERPOLATES SOIL PROPERTIES FROM LAYERS WITH   
      ! EQUAL THICKNESS (OUTPUT FROM DIF EQ SOLN OF GAS DIFFUSION EQS) TO
      ! LAYERS OF UNEQUAL THICKNESS (INPUT SOIL LAYERS).
      USE PARM
      DIMENSION X(MSC,MSA),Y(MSL,MSA)
      Z1=0.
      TOT=0.
      J=1
      DO K=1,N2
          DO WHILE(J<=N1)
              L=LID(J,ISA)
              IF(Z(L,ISA)>ZC(K,ISA))EXIT
              Y(L,ISA)=TOT+X(K,ISA)*(Z(L,ISA)-Z1)
              Z1=Z(L,ISA)
              J=J+1
              TOT=0.
          END DO
          IF(J<=N1)THEN
              TOT=TOT+X(K,ISA)*(ZC(K,ISA)-Z1)
              Z1=ZC(K,ISA)
          ELSE
              EXIT  
          END IF
      END DO
      Z1=0.
      DO J=1,N1
          L=LID(J,ISA)     
          Y(L,ISA)=Y(L,ISA)/(Z(L,ISA)-Z1)
          Z1=Z(L,ISA)
      END DO
      RETURN
      END