      SUBROUTINE TMXL1(DMX,T,W,X,Y)
!     THIS SUBPROGRAM PREVENTS TILLAGE FROM INCREASING TOP LAYER CONTENT
!     BY MIXING
      T=(T*DMX-X+Y)/(DMX-.01)
      W=X
      RETURN
      END