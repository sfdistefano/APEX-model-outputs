      SUBROUTINE WREAD(IU,KK)
!     APEX1905
!     THIS SUBPROGRAM READS THE DAILY WEATHER FILE TO THE DAY BEFORE THE
!     SIMULATION BEGINS.
      USE PARM 
      DIMENSION MOFD(12)
      DATA MOFD/31,29,31,30,31,30,31,31,30,31,30,31/
      IF(KK==3)THEN
          READ(KRST(IU),*)JDA,I3
          MO=1
          CALL AXMON(JDA,MO)
		  CALL AICL
		  I1=KDA
		  I2=MO
	  ELSE
	 	  READ(KRST(IU),2)I3,I2,I1,(XTP(L),L=1,5)
	  END IF
      J3=10000*I3
      J1=100*I2+J3
      II=I1+J1
      K1=I2
	  IF(KK==1)THEN
	      IF(II<IBDT)GO TO 9
          WRITE(KW(36),5)IBDT,II,FWTH(NDWT)
          IBDT=II
	      IYR0=I3
	      GO TO 9
	  END IF
      IF(II==IBDT)GO TO 9
      DO
	      CALL ALPYR(I3,NYD,LPYR)
          DO I2=K1,12
              N1=MOFD(I2)
              IF(I2==2)N1=N1-NYD
              J2=100*I2
              J1=J2+J3
              DO WHILE(I1<N1)
                  I1=I1+1
                  II=J1+I1
                  IF(II==IBDT)RETURN
                  READ(KRST(IU),'()',IOSTAT=NFL)
                  IF(NFL/=0)THEN
                      IF(KK==3)THEN
                          WRITE(KW(1),*)'START DATE EXCEEDS PSO FILE'
                          RETURN
                      ELSE    
                          WRITE(KW(1),*)'START DATE EXCEEDS WTH FILE--&
                          WEATHER GENERATED'
                          NGN=-1
                          KGN=0
                          RETURN
                      END IF    
                  END IF
              END DO
              I1=0
          END DO
          K1=1
          I2=1
          I3=I3+1
          J3=10000*I3
      END DO
    9 REWIND KRST(IU)
      RETURN
    1 FORMAT(14X,7F6.0)
    2 FORMAT(2X,3I4,6F6.0)
    5 FORMAT(T10,'BEGINNING DATE CHANGED FROM',I10,' TO',I10,1X,A80)       
      END