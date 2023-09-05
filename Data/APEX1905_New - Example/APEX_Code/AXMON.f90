      SUBROUTINE AXMON(JDX,MOX)
!     APEX1905
!     THIS SUBPROGRAM DETERMINES THE MONTH, GIVEN THE DAY OF THE YEAR.
      USE PARM
      IF(JDX>NC(2))THEN
          M=MOX
          DO MOX=M,12
              M1=MOX+1
              NDA=NC(M1)-NYD
              IF(JDX<=NDA)RETURN
          END DO
      END IF
      MOX=1
      RETURN
      END