      SUBROUTINE AINTRO(X,Y,N1,N2)
      ! APEX1905
      ! THIS SUBPROGRAM INTERPOLATES SOIL PROPERTIES FROM LAYERS WITH   
      ! EQUAL THICKNESS (OUTPUT FROM DIF EQ SOLN OF GAS DIFFUSION EQS) TO
      ! LAYERS OF UNEQUAL THICKNESS (INPUT SOIL LAYERS).
      USE PARM
      DIMENSION X(MSC,MSA),Y(MSC,MSA)
      Z1=0.
      TOT=0.
      AD1=0.
      AD2=0.
      J=1
      DO K=1,N2
          AD1=AD1+X(K,ISA)
          DO WHILE(J<=N1)
              L=LID(J,ISA)
              IF(Z(L,ISA)>ZC(K,ISA))EXIT
              Y(L,ISA)=TOT+X(K,ISA)*(Z(L,ISA)-Z1)/DZDN
              AD2=AD2+Y(L,ISA)
              Z1=Z(L,ISA)
              J=J+1
              TOT=0.
          END DO
          IF(J<=N1)THEN
              TOT=TOT+X(K,ISA)*(ZC(K,ISA)-Z1)/DZDN
              Z1=ZC(K,ISA)
          ELSE
              EXIT  
          END IF
      END DO
      L=LID(N1,ISA)     
      Y(L,ISA)=Y(L,ISA)+AD1-AD2
      RETURN
      END