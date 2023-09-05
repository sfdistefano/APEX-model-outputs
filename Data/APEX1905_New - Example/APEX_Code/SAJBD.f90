      SUBROUTINE SAJBD
!     APEX1905
!     THIS SUBPROGRAM SIMULATES THE CHANGE IN BULK DENSITY WITHIN THE
!     PLOW LAYER CAUSED BY INFILTRATION SETTLING.
      USE PARM 
      XX=RFV(IRF(ISA))-QVOL(IDO)
      DO J=1,NBSL(ISA)
          ISL=LID(J,ISA)
          IF(XX>0.)THEN
              XX=XX*.2*(1.+2.*SAN(ISL,ISA)/(SAN(ISL,ISA)+EXP(8.597-.075*&
              SAN(ISL,ISA))))/Z(ISL,ISA)**.6
              IF(SCLM(6)>0.)XX=MIN(XX,SCLM(6))
              F=XX/(XX+EXP(SCRP(6,1)-SCRP(6,2)*XX))
              BDP(ISL,ISA)=BDP(ISL,ISA)+F*(BD(ISL,ISA)-BDP(ISL,ISA))
              XX=PKRZ(ISL)
          ELSE
              IF(Z(ISL,ISA)>BIG(ISA))RETURN
          END IF
      END DO
      RETURN
      END