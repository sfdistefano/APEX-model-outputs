      SUBROUTINE SOIL_BDCD(BDCD,I)
      ! THIS SUBPROGRAM ESTIMATES THE SOIL BULK DENSITY CODE GIVEN TEXTURE
      ! USED TO CALCULATE PARM(2)
      USE PARM
      IF(SAN(I,ISA)>85.)THEN
          ! SANDY
          BDCD=1.65
          RETURN
      END IF
      IF(SAN(I,ISA)>70..AND.CLA(I,ISA)<15.)THEN
          ! SANDY
          BDCD=1.65
          RETURN
      END IF
      IF(SAN(I,ISA)<15.)THEN
          IF(CLA(I,ISA)<18.)THEN
              ! COURSE SILTY
              BDCD=1.35
              RETURN
          END IF    
          IF(CLA(I,ISA)<35.)THEN
              !FINE SILTY
              BDCD=1.3
              RETURN
          END IF
      ELSE
          IF(CLA(I,ISA)<18.)THEN
              ! COURSE LOAMY 
              BDCD=1.55
              RETURN
          END IF
          IF(CLA(I,ISA)<35.)THEN
              !FINE LOAMY
              BDCD=1.45
              RETURN
          END IF
      END IF    
      IF(CLA(I,ISA)>60.)THEN
          ! VERY FINE
          BDCD=1.15
          RETURN
      END IF    
      BDCD=1.15
      RETURN
      END