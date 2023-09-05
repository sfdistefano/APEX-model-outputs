      SUBROUTINE AUAC(MY)
      ! APEX1905
      ! THIS SUBPROGRAM CONVERTS FROM MASS TO MASS/UNIT AREA
      USE PARM
      DO J=1,NSM                                                                     
          SELECT CASE(MASA(J))
              CASE(0)
                  X1=1.
              CASE(1)
                  X1=WSA(ISA)
              CASE(2)
                  X1=10.*WSA(ISA)
          END SELECT
          SELECT CASE(MY)
              CASE(0) ! DAILY
                  VARUA(J,ISA)=VAR(J,ISA)/X1
              CASE(1) ! MONTHLY
                  SMMUA(J,MO1,ISA)=SMM(J,MO1,ISA)/X1
                  !SMYUA(J,ISA)=SMY(J,ISA)/X1
              CASE(2) ! ANNUAL
                  SMYUA(J,ISA)=SMY(J,ISA)/X1
              CASE(3) ! TOTAL
                  SMUA(J,ISA)=SM(J,ISA)/X1
          END SELECT
      END DO
      RETURN
      END