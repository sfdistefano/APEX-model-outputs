      SUBROUTINE EWIK
!     THIS SUBPROGRAM ESTIMATES THE SOIL ERODIBILITY FACTOR FOR THE WIND
!     EROSION EQ.
      USE PARM
      LD1=LID(1,ISA)
      IF(SAN(LD1,ISA)>85.+.5*CLA(LD1,ISA))THEN
          WK=1.
          RETURN
      END IF
      IF(SAN(LD1,ISA)>70.+CLA(LD1,ISA))THEN
          WK=.43
          RETURN
      END IF
      IF(SIL(LD1,ISA)>80..AND.CLA(LD1,ISA)<12.)THEN
          WK=.12
          RETURN
      END IF
      IF(CAC(LD1,ISA)>0.)THEN
          IF(SAN(LD1,ISA)<45..OR.CLA(LD1,ISA)<20..OR.SIL(LD1,ISA)>28.)THEN
              WK=.28
              RETURN
          ELSE
              WK=.18
              RETURN
          END IF
      END IF  
      IF(CLA(LD1,ISA)<7.)THEN
          IF(SIL(LD1,ISA)<50.)THEN
              WK=.28
              RETURN
          ELSE
              WK=.18
              RETURN
          END IF
      END IF        
      IF(CLA(LD1,ISA)<20.)THEN
          IF(SAN(LD1,ISA)>52.)THEN
              WK=.28
              RETURN
          ELSE
              WK=.18
              RETURN
          END IF
      END IF        
      IF(CLA(LD1,ISA)<27.)THEN
          IF(SIL(LD1,ISA)<28.)THEN
              WK=.18
              RETURN
          ELSE
              WK=.16
              RETURN
          END IF
      END IF        
      IF(CLA(LD1,ISA)<35..AND.SAN(LD1,ISA)<20.)THEN
          WK=.12
          RETURN
      END IF
      IF(CLA(LD1,ISA)<35.)THEN
          IF(SAN(LD1,ISA)<45.)THEN        
              WK=.16
              RETURN
          ELSE
              WK=.18
              RETURN
          END IF
      END IF        
      IF(SAN(LD1,ISA)>45.)THEN
          WK=.18
          RETURN
      ELSE
          WK=.16 !Susan to improve wind erosion, 2020
          RETURN
      END IF
      RETURN
      END