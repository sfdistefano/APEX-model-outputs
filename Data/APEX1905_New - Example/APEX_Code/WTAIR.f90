      SUBROUTINE WTAIR(I)
!     APEX1905
!     THIS SUBPROGRAM SIMULATES MAXIMUM & MINIMUM DAILY AIR TEMPERATURE
!     FROM A NORMAL DISTRIBUTION.
      USE PARM  
      TMX(I)=TXXM+SDTMX(IWI,MO)*WX(1)
      TMN(I)=OBMN(IWI,MO)+SDTMN(IWI,MO)*WX(2)
      IF(TMN(I)>TMX(I))TMN(I)=TMX(I)-.2*ABS(TMX(I))
      RETURN
      END