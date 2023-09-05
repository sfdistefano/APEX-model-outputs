      SUBROUTINE COMP2EXP(KW50)
!     APEX1905
!     THIS SUBPROGRAM COMPARE EXPERIMENTAL RESULTS WITH MODELED RESULTS
      USE PARM
      CHARACTER(80):: STRING,COPYCAG1,NAME,CTIME 
      CHARACTER(4):: ECPNM(MNC),SCPNM,VORCP(MNC),SCPNM1
      CHARACTER(8):: values
	  LOGICAL::XMIS1,XMIS2,XMIS3,XMIS4
      INTEGER::EYEAR,EMON,EDAY,SYEAR,SMON,SDAY,EHERD,SHERD,EOWNER,SOWNER,EPAID,SPAID,KW50
      INTEGER::K,L,J, VALCAG1(13), IREAD, LCC, NOBS
      REAL::ESTLD(MNC),SSTL,SSTD,EWTG,SWTG,SGZSL,SGZSD,ESTLD1,EXPOBS(100),SIMDATA(100)
      REAL::RSQ,SLOPE,SINTER,RMSECAGE,RMSEVOR,RMSEWTG,TTEST,T5
!      CALL DATE_AND_TIME(VALUES=values)
!      CALL TIME(VALUES)
      ASTN=ASTN1
      OPEN(KW50,FILE=TRIM(ADJUSTL(ASTN))//'.CMP')
	  INQUIRE(FILE=TRIM(ADJUSTL(ESTN))//'.CAG',EXIST=XMIS1)
      NOBS=0
	  IF(XMIS1==.TRUE.)THEN
	      OPEN(KR(64),FILE=TRIM(ADJUSTL(ESTN))//'.CAG')
          OPEN(KR(62),FILE=TRIM(ADJUSTL(ASTN))//'.SAD')        !READ CAGE BIOMASS DATA, KW(49), KW(11) FOR .SAD FILE
          DO L=1,9
          READ (KR(62),'(A80)') STRING
          ENDDO
          READ (KR(64),'(A80)') STRING
          WRITE(KW50,*) 'CAGE_B'
          WRITE(KW50,*) TRIM(STRING)//'  S_BIOMAS'
          LCC=6
          DO I=1,LCC
          READ (KR(64),157) EPAID,EYEAR,EMON,EDAY,ECPNM(I),ESTLD(I)
          END DO
          DO WHILE (.NOT.EOF(KR(62)))   !OR DO WHILE (.TRUE.)
          READ (KR(62),154) SPAID,SYEAR,SMON,SDAY,SCPNM,SSTL,SSTD,SGZSL,SGZSD
124       IF (SYEAR>EYEAR) THEN
          if (.NOT.EOF(KR(64))) THEN
          DO I=1,LCC
          READ (KR(64),157) EPAID,EYEAR,EMON,EDAY,ECPNM(I),ESTLD(I)
          END DO
          GO TO 124
          ELSE
          GO TO 325
          END IF
          END IF
          IF (SPAID==EPAID .AND. SYEAR==EYEAR .AND. SMON==EMON .AND. SDAY==EDAY) then
              DO J=1,LCC
              IF (SCPNM==ECPNM(J)) THEN
                WRITE(KW50,155) EPAID,EYEAR,EMON,EDAY,ECPNM(J),ESTLD(J),(SSTL+SSTD)*1000+SGZSL+SGZSD
                NOBS=NOBS+1
                EXPOBS(NOBS)=ESTLD(J)
                SIMDATA(NOBS)=(SSTL+SSTD)*1000+SGZSL+SGZSD
              ENDIF
              END DO
              IREAD=1
          ELSE
          if (.NOT.EOF(KR(64)).AND.IREAD==1) THEN
          DO I=1,LCC
              READ (KR(64),157) EPAID,EYEAR,EMON,EDAY,ECPNM(I),ESTLD(I)
          ENDDO
          IREAD=0
          ENDIF
          ENDIF
          END DO
	  ELSE
          WRITE(*,'(/A/)')'File '//TRIM(ADJUSTL(ESTN)//'.CAG')//' IS MISSING. NO CAGE BIOMASS COMPARED'
	  END IF	
325   CONTINUE
      CALL STATS(EXPOBS,SIMDATA,NOBS,RSQ,SLOPE,SINTER,RMSECAGE,TTEST,T5)
      WRITE(KW50,FMT=174) RMSECAGE,RSQ,SLOPE,SINTER,TTEST,T5
      Close (KR(62))
      NOBS=0
	  INQUIRE(FILE=TRIM(ADJUSTL(ESTN))//'.VOR',EXIST=XMIS2)
	  IF(XMIS2==.TRUE.)THEN
	      OPEN(KR(66),FILE=TRIM(ADJUSTL(ESTN))//'.VOR')
          OPEN(KR(67),FILE=TRIM(ADJUSTL(ASTN))//'.SAD')        !READ CAGE BIOMASS DATA, KW(49), KW(11) FOR .SAD FILE
          REWIND(KR(67))
          DO L=1,9
          READ (KR(67),'(A80)') STRING
          END DO
          READ (KR(66),'(A80)') STRING
          K=INUMB(STRING)
          BACKSPACE(KR(66))
          READ (KR(66),*) STRING,(VORCP(L),L=1,K-1)
          READ (KR(66),'(A80)') STRING
          WRITE(KW50,*)
          WRITE(KW50,*) 'VOR_B'
          WRITE(KW50,*) TRIM(STRING)//'       S_BIOMAS'
          READ (KR(66),*) EPAID,EYEAR,EMON,EDAY,ESTLD1
          DO WHILE (.NOT.EOF(KR(67)))   !OR DO WHILE (.TRUE.)
          READ (KR(67),154) SPAID,SYEAR,SMON,SDAY,SCPNM,SSTL,SSTD
224       IF (SYEAR>EYEAR) THEN
          if (.NOT.EOF(KR(66))) THEN
          READ (KR(66),*) EPAID,EYEAR,EMON,EDAY,ESTLD1
          GO TO 224
          ELSE
          GO TO 525
          END IF
          END IF
          TSTLD=0.
          IF (SPAID==EPAID .AND. SYEAR==EYEAR .AND. SMON==EMON .AND. SDAY==EDAY) THEN
324       DO L=1,K-1
          IF (SCPNM==VORCP(L)) TSTLD=TSTLD+(SSTL+SSTD)*1000.
          END DO
          if (.NOT.EOF(KR(67))) READ (KR(67),154) SPAID,SYEAR,SMON,SDAY,SCPNM,SSTL,SSTD
          IF (SPAID==EPAID .AND. SYEAR==EYEAR .AND. SMON==EMON .AND. SDAY==EDAY) THEN
          GOTO 324
          ELSE
          WRITE(KW50,156) EPAID,EYEAR,EMON,EDAY,ESTLD1,TSTLD
          NOBS=NOBS+1
          EXPOBS(NOBS)=ESTLD1
          SIMDATA(NOBS)=TSTLD
          BACKSPACE(KR(67))
          END IF
          if (.NOT.EOF(KR(66))) READ (KR(66),*) EPAID,EYEAR,EMON,EDAY,ESTLD1
          ENDIF
          END DO
	  ELSE
          WRITE(*,'(/A/)')'File '//TRIM(ADJUSTL(ESTN)//'.VOR')//' IS MISSING. NO VOR BIOMASS COMPARED'
	  END IF		  
525   CONTINUE
      CALL STATS(EXPOBS,SIMDATA,NOBS,RSQ,SLOPE,SINTER,RMSEVOR,TTEST,T5)
      WRITE(KW50,174) RMSEVOR,RSQ,SLOPE,SINTER,TTEST,T5
      Close (KR(67))
      NOBS=0
      INQUIRE(FILE=TRIM(ADJUSTL(ESTN))//'.WTG',EXIST=XMIS3)
	  IF(XMIS3==.TRUE.)THEN
	      OPEN(KR(65),FILE=TRIM(ADJUSTL(ESTN))//'.WTG')
          OPEN(KR(63),FILE=TRIM(ADJUSTL(ASTN))//'.AGZ')        !READ WEIGT GAIN DATA, KW(33)
          REWIND(KR(63))
          DO I=1,9
          READ (KR(63),'(A80)') STRING
          END DO
          READ (KR(65),'(A80)') STRING
          WRITE(KW50,*)
          WRITE(KW50,*) 'WTG'
          WRITE(KW50,'(13X,A40)') TRIM(STRING)//'      S_WTG'
          READ (KR(65),*) EYEAR,EOWNER,EHERD,EWTG
          READ (KR(63),254) SYEAR,SCPNM1,SOWNER,SHERD,SWTG
          BACKSPACE KR(63)
          DO WHILE (.NOT.EOF(KR(63)))   !OR DO WHILE (.TRUE.)
!          DO I=1,LC
          READ (KR(63),254) SYEAR,SCPNM,SOWNER,SHERD,SWTG
!          END DO
123       IF (SYEAR>EYEAR) THEN
          if (.NOT.EOF(KR(65))) THEN
          READ (KR(65),*) EYEAR,EOWNER,EHERD,EWTG
          GO TO 123
          ELSE
          GO TO 625
          END IF
          END IF
          IF (SYEAR==EYEAR .AND. SOWNER==EOWNER .AND. SHERD==EHERD .AND. SCPNM==SCPNM1) THEN
          WRITE(KW50,255) EYEAR,EOWNER,EHERD,EWTG,SWTG
          NOBS=NOBS+1
          EXPOBS(NOBS)=EWTG
          SIMDATA(NOBS)=SWTG
          if (.NOT.EOF(KR(65))) READ (KR(65),*) EYEAR,EOWNER,EHERD,EWTG
          END IF
          END DO
	  ELSE
          WRITE(*,'(/A/)')'File '//TRIM(ADJUSTL(ESTN)//'.WTG')//' IS MISSING. NO WEIGHT GAIN COMPARED'
	  END IF
      
625   CONTINUE
      CALL STATS(EXPOBS,SIMDATA,NOBS,RSQ,SLOPE,SINTER,RMSEWTG,TTEST,T5)
      WRITE(KW50,174) RMSEWTG,RSQ,SLOPE,SINTER,TTEST,T5
      WRITE(KW50,175) RMSECAGE,RMSEVOR,RMSEWTG
      CLOSE(ALL)
      
      RETURN
  154 FORMAT(9X,I8,1X,3I4,6X,A4,60X,1X,F11.5,12X,1X,F11.5,12X,2(1X,F11.5))      
  155 FORMAT(1X,I8,1X,3(4X,I4),4X,A4,2(1X,F9.2))      
  156 FORMAT(1X,I8,1X,3(4X,I4),3X,2(1X,F9.2))      
  254 FORMAT(18X,I4,7X,A4,2(1X,I4),81X,2F10.3)                                                                                          
  255 FORMAT(14X,I4,1X,I5,4X,I5,2(2X,F8.3))      
  157 FORMAT(4X,I4,3(4X,I4),4X,A4,1X,F9.2)      
  174 FORMAT (/10X,'RMSE',11X,'R2',11X,'SLOPE',6X,'INTERCEPT',5X, &
              'T-SLOPE',7X,'T-5% CONFIDENCE'/,5X,6(1X,F12.4)/)
  175 FORMAT (/10X,'RMSECAGE',5X,'RMSEVOR',8X,'RMSEWTG'/,5X,6(1X,F12.4)/)
      END
      
      INTEGER FUNCTION INUMB(STRING)
!
      CHARACTER STRING*(*)
      INTEGER J,I
!
      J=itrim(STRING)
       inumb=0
       do i=1, j
           if ((((iachar(string(i:i)).ne.9).and. &
            (iachar(string(i:i)).ne.32)).and. &
            ((iachar(string(i+1:i+1)).eq.9).or. &
            (iachar(string(i+1:i+1)).eq.32).or. &
            (iachar(string(i+1:i+1)).eq.13)))) then
!        if ((string(i:i).ne.' ').and.(string(i+1:i+1).eq.' ')) then
	    inumb=inumb+1
        endif
!	  write (*,*) iachar(string(i:i)),string(i:i)
	enddo
      END
      
      INTEGER FUNCTION ITRIM(STRING)
!
      CHARACTER STRING*(*)
      INTEGER J
!
      J=LEN(STRING)
   10 IF(J.GT.1.AND.STRING(J:J).EQ.' ') THEN
        J=J-1
        GOTO 10
      ENDIF
      IF(J.EQ.1.AND.STRING(1:1).EQ.' ') THEN
        ITRIM=0
      ELSE
        ITRIM=J
      ENDIF
    END
!
      SUBROUTINE STATS(YY,YYP,NOB,RSQ,SLOPE,SINTER,RMSE,TTEST,TVAR) 
      IMPLICIT REAL (A-H,O-Z)
      DIMENSION YY(nob),YYP(nob)
      PARAMETER (NP=2)
      YSUM=0.0
      FSUM=0.0
      YFSUM=0.0
      YSSQ=0.0
      FSSQ=0.0
      YFSSQ=0.0
      NOBREAL = 0
      RESSSQ = 0.0
      RSQ = -99.0
      SLOPE = -99.0
      SINTER = -99.0
      RMSE = -99.0
      TTEST = -99.0
      TVAR = -99.0
      DO 501 I=1,NOB
        IF (YY(I).GT.0.0.AND.YYP(I).GT.0.0) THEN
          NOBREAL = NOBREAL + 1
          YSUM=YSUM+YY(I)
          FSUM=FSUM+YYP(I)
          YSSQ=YSSQ+YY(I)*YY(I)
          FSSQ=FSSQ+YYP(I)*YYP(I)
          YFSSQ=YFSSQ+YY(I)*YYP(I)
          RESSSQ=RESSSQ + (YY(I)-YYP(I))**2
        ENDIF
501   CONTINUE
      IF (NOBREAL.EQ.0) RETURN
      RMSE=SQRT(RESSSQ/NOBREAL)
      SXX=YSSQ-YSUM*YSUM/NOBREAL
      SYY=FSSQ-FSUM*FSUM/NOBREAL
      SXY=YFSSQ-YSUM*FSUM/NOBREAL
      YBAR=YSUM/NOBREAL
      FBAR=FSUM/NOBREAL
      IF (SXX.EQ.0.0) RETURN
      SLOPE = SXY/SXX
      SINTER = FBAR - SLOPE * YBAR
      IF (SYY.EQ.0.0) RETURN
      RSQ = SXY**2/(SXX*SYY)
      IF (NOBREAL.EQ.NP) RETURN
      Z=1.0/(NOBREAL-NP)
      TVAR=1.96+Z*(2.3779+Z*(2.7135+Z*(3.187936+2.4666*Z**2)))
      SYX = SQRT(abs(SYY-SLOPE*SXY)/(NOBREAL-NP))
      IF (SYX.NE.0.0) TTEST = (SLOPE-1.0)/(SYX/SQRT(SXX))
      RETURN
      END
      