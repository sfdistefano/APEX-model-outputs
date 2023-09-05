      SUBROUTINE HGWST
      ! APEX1905
      ! THIS SUBPROGRAM TRANSFERS GROUNDWATER AMONG SUBAREAS
      USE PARM
      DIMENSION GWS1(MSA)
      DO I=MSA,2,-1
          II=JSA(I)
          RN=AUNIF(IDG(13))
          K=I*RN+1
          JSA(I)=JSA(K)
          JSA(K)=II
          GWS1(I)=GWST(I)
      END DO
      GWS1(1)=GWST(1)
      M1=MSA-1
      DO I=1,M1
          I1=JSA(I)
          X1=GWE0(I1)-(1.-GWST(I1)/GWMX(I1))*5.
          DO J=I+1,MSA
              J1=JSA(J)
              X2=GWE0(J1)-(1.-GWST(J1)/GWMX(J1))*5.
              IF(X1>X2)THEN
                  SGN=1.
                  II=I1
              ELSE
                  SGN=-1.
                  II=J1
              END IF
              XX=MIN(.1,1.-EXP(SGN*(X2-X1)*RFTT(II)/(SADST(I1,J1)+.01)))
              X3=SGN*XX*GWST(II)
              X4=SGN*(GWST(I1)*GWMX(J1)-GWST(J1)*GWMX(I1))/(GWMX(I1)+GWMX(J1))
              X5=MIN(ABS(X3),ABS(X4))
              GWTR=SGN*X5    
              GWST(I1)=GWST(I1)-GWTR
              GWST(J1)=GWST(J1)+GWTR
          END DO
      END DO
      DO I=1,MSA
         GWEL(I)=GWE0(I)-(1.-GWST(I)/GWMX(I))*5. 
         VAR(157,I)=GWST(I)-GWS1(I)
         SMM(157,MO,I)=SMM(157,MO,I)+VAR(157,I)
      END DO   
      RETURN
      END
      