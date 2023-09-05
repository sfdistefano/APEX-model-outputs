      SUBROUTINE WTAIX(I)
!     APEX1905
!     THIS SUBPROGRAM GENERATES TMX AND TMN FOR MISSING RECORDS (999)
      USE PARM  
      IF(TMX(I)<100.)THEN
          TMN(I)=MIN(OBMN(IWI,MO)+SDTMN(IWI,MO)*WX(2),TMX(I)-.2*ABS(TMX(I)))
      ELSE
          TMX(I)=MAX(TXXM+SDTMX(IWI,MO)*WX(1),TMN(I)+.2*ABS(TMN(I)))
      END IF
      RETURN
      END