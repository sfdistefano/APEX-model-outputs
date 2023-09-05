      SUBROUTINE SWRTNR(CL,SA,OC,WP,FC)
!     APEX1905
!     THIS SUBPROGRAM USES WALTER RAWLS'S METHOD FOR ESTIMATING SOIL
!     WATER CONTENT AT 33 AND 1500 KPA.
      WP=.026+.005*CL+.0158*OC
      FC=.2576-.002*SA+.0036*CL+.0299*OC
      RETURN
      END