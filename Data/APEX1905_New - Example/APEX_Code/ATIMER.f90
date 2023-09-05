      SUBROUTINE ATIMER(ITR)
!     APEX1905
!     THIS SUBPROGRAM SETS DATE AND TIME FOR OUTPUT AND CALCULATES
!     ELAPSED TIME.
      USE PARM 
      INTEGER :: values(8)
      IF(ITR==0)THEN
          CALL DATE_AND_TIME(VALUES=values)
          IYER = values(1)
          IMON = values(2)
          IDAY = values(3)
          IT1 = values(5)
          IT2 = values(6)
          IT3 = values(7)
          I100 = values(8)
          RETURN
      END IF
      CALL DATE_AND_TIME(VALUES=values)
      IY = values(1)
      MO = values(2)
      IDA = values(3)
      IEH = values(5)
      IEM = values(6)
      IES = values(7)
      I100 = values(8)
      CALL ALPYR(IY,NYD,LPYR)
      CALL ADAJ(NC,IEDT,MO,IDA,NYD)
      CALL ALPYR(IYER,NYD,LPYR)
      CALL ADAJ(NC,IBDX,IMON,IDAY,NYD)
      I1=86400*((IY-IYER)*(366-NYD)+IEDT-IBDX)
      IBT=IT1*3600+IT2*60+IT3
      IEX=IEH*3600+IEM*60+IES
      WRITE(KW(1),5000)IT1,IT2,IT3
      WRITE(KW(1),5100)IEH,IEM,IES
      WRITE(KW(1),5200)
      II=I1+IEX-IBT
      ITS=MOD(II,60)
      II=(II-ITS)/60
      ITM=MOD(II,60)
      ITH=(II-ITM)/60
      WRITE(KW(1),5320)ITH,ITM,ITS
      RETURN
 5000 FORMAT(10X,'BEGINNING TIME: ',I2,2(':',I2))
 5100 FORMAT(10X,'ENDING    TIME: ',I2,2(':',I2))
 5200 FORMAT(10X,'----------------------------')
 5320 FORMAT(10X,'TOTAL RUN TIME: ',I2,2(':',I2))
      END