      SUBROUTINE NUSE
!     APEX1905
!     THIS SUBPROGRAM CALCULATES THE DAILY POTENTIAL SOIL SUPPLY OF P
!     FOR EACH LAYER.
      USE PARM
      XX=1.5*UPM/RW(JJK,ISA)
      DO J=1,LRD(ISA)
          ISL=LID(J,ISA)
          UN(ISL,ISA)=WNO3(ISL,ISA)*UW(ISL,ISA)/(SWST(ISL,ISA)+.001)
          SUN=SUN+UN(ISL,ISA)
          F=1000.*WPML(ISL,ISA)/WT(ISL,ISA)
          IF(F>30.)THEN
              F=1.
          ELSE
              IF(SCLM(11)>0.)F=MIN(F,SCLM(11))
              F=F/(F+EXP(SCRP(11,1)-SCRP(11,2)*F))
          END IF
          UP(ISL,ISA)=XX*F*RWT(ISL,JJK,ISA)
          X1=WPML(ISL,ISA)+WPMU(ISL,ISA)
          IF(UP(ISL,ISA)>=X1)UP(ISL,ISA)=.9*X1
          SUP=SUP+UP(ISL,ISA)
      END DO
      RETURN
      END