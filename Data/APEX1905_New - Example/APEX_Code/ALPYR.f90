      SUBROUTINE ALPYR(IYR,NYD,LPYR)
!     APEX1905
      NYD=1
      IF(MOD(IYR,100)/=0)THEN
          IF(MOD(IYR,4)/=0)RETURN
      ELSE
          IF(MOD(IYR,400)/=0)RETURN
      END IF
      NYD=LPYR
      RETURN 
      END