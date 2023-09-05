      SUBROUTINE WHLRMX(II)
!     APEX1905
!     THIS SUBPROGRAM COMPUTES DAY LENGTH & MAX SOLAR RADIATION AT THE
!     EARTHS SURFACE.
      USE PARM 
      XI=II
      SD=.4102*SIN((XI-80.25)/PIT)
      CH=-YTN(ISA)*TAN(SD)
      IF(CH>=1.)THEN
          H=0.
      ELSE
          IF(CH<=-1.)THEN
              H=3.1416
          ELSE
              H=ACOS(CH)
          END IF
      END IF  
      DD=1.+.0335*SIN((XI+88.2)/PIT)
      HRLT=7.72*H
      HR1=HRLT-HR0(ISA)
      HR0(ISA)=HRLT
      RAMX=30.*DD*(H*YLS(ISA)*SIN(SD)+YLC(ISA)*COS(SD)*SIN(H))
      RETURN
      END