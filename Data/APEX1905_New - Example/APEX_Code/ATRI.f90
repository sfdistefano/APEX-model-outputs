      FUNCTION ATRI(BLM,QMN,UPLM,KK)
!     APEX1905
!     THIS SUBPROGRAM GENERATES NUMBERS FROM A TRIANGULAR DISTRIBUTION
!     GIVEN X AXIS POINTS AT START & END AND PEAK Y VALUE.
      USE PARM
      U3=QMN-BLM
      RN=AUNIF(IDG(KK))
	  IF(RN>=1.)THEN
	      ATRI=1.
	  ELSE
	      Y=2./(UPLM-BLM)
          B2=UPLM-QMN
          B1=RN/Y
          X1=Y*U3/2.
          IF(RN>X1)THEN
              ATRI=UPLM-SQRT(B2*B2-2.*B2*(B1-.5*U3))
          ELSE
              ATRI=SQRT(2.*B1*U3)+BLM
          END IF
      END IF
      IF(KK==7.OR.KK==4)THEN
          AMN=(UPLM+QMN+BLM)/3.
          ATRI=ATRI*QMN/AMN
          IF(ATRI>=1.)ATRI=.99
      END IF
      RETURN
      END