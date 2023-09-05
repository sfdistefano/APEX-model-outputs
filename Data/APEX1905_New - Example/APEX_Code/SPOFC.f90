      SUBROUTINE SPOFC(I)
!     APEX1905
!     THIS SUBPROGRAM CHECKS SOIL WATER STORAGE CAPACITY AT WILTING
!     POINT, FIELD CAPACITY, & POROSITY & ADJUSTS INCONSISTANCIES AS
!     NEEDED
      USE PARM 
      X1=.95*PO(I,ISA)
      IF(FC(I,ISA)>X1)THEN
          X2=MAX(1.,FC(I,ISA)-S15(I,ISA))
          FC(I,ISA)=X1
          S15(I,ISA)=FC(I,ISA)-X2
          IF(S15(I,ISA)<=0.)S15(I,ISA)=.01*FC(I,ISA)
      END IF
      RETURN
      END