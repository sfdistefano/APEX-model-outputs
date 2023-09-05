      SUBROUTINE CANOPY(CANOPYD,CANOPYV,MCANOPY,NREPLACE,NBASAL,XMIS,ISWITCH)
!     APEX1905
!     THIS SUBPROGRAM READS CANOPY CROP DATA, LIWANG MA
      USE PARM
      INTEGER CANOPYD(50,LC,MSA1),MCANOPY,ISWITCH,NREPLACE,NBASAL
      REAL CANOPYV(50,LC,MSA1)   !CANOPYV(50,LC,MSA1)
      CHARACTER(80):: ADUM
	  LOGICAL::XMIS	
      FNPP(X)=DMLA(JJK)*X/(X+EXP(PPCF(1,JJK)-PPCF(2,JJK)*X))
      IF (ISWITCH==0) THEN
      !ADUM=ESTN//'.CAN'
	  INQUIRE(FILE=TRIM(ADJUSTL(ESTN))//'.CAN',EXIST=XMIS)
	  IF(XMIS==.TRUE.)THEN
	      OPEN(KR(61),FILE=TRIM(ADJUSTL(ESTN))//'.CAN')
	  ELSE
          WRITE(*,'(/A/)')'File '//TRIM(ADJUSTL(ESTN//'.CAN'))//' IS MISSING. NO BASAL AREA UPDATED'
          RETURN
	  END IF	
      IBDT1=0
      I=0
      MCANOPY=0
      READ(KR(61),*) NREPLACE
      READ(KR(61),*) NBASAL
      IF (NREPLACE.GT.0) THEN
      DO WHILE (.NOT.EOF(KR(61)))   !OR DO WHILE (.TRUE.)
      READ(KR(61),2470,IOSTAT=NFL)JX1,JX2,JX3,JX4,JX5,OPV5                                            
      IF (IBDT1>(JX3+100*JX2+10000*JX1)) THEN
            I=0
            IBDT1=0
      END IF
      IF (IBDT1/=(JX3+100*JX2+10000*JX1)) THEN
            I=I+1
            IBDT1=JX3+100*JX2+10000*JX1
            MCANOPY=MAX(I,MCANOPY)
      ENDIF      
      IF (KDC1(JX5).GT.0) THEN
      CANOPYD(I,KDC1(JX5),JX4)=JX3+100*JX2+10000*JX1
      CANOPYV(I,KDC1(JX5),JX4)=OPV5
      ENDIF
 2470 FORMAT(I4,2I3,2I5,F10.0)
      END DO
      ENDIF
      CLOSE(KR(61))
      ELSE IF (ISWITCH==1) THEN
      DO I=1,LC ! FORAGE LOOP													
      DO J=1,MSA  !SUBAREA LOOP
      JJK=I
      Y1=PPLP(1,JJK)
      Y2=PPLP(2,JJK)
      IF(EXTC(LC)<1.E-3)EXTC(LC)=0.65
      IF(Y2>Y1)THEN
          X4=Y1
          X5=Y2
      ELSE
          X4=Y2
          X5=Y1
      END IF 
      X1=ASPLT(X4)
      X2=ASPLT(X5)
      CALL ASCRV(X4,X5,X1,X2,SCLM,31,KW(1))
      PPCF(1,JJK)=X4
      PPCF(2,JJK)=X5
      IF (NREPLACE==1 .AND. NBASAL==0) THEN
      DO K=1,MCANOPY
      IF (CANOPYD(K,JJK,J)==KDA+100*MO+10000*IYR) THEN
	  IF(CANOPYV(K,JJK,J)>0..AND. NREPLACE==1)THEN   !if it is zero, should omit the whole upate canopy program, Liwang Ma
          X3=CANOPYV(K,JJK,J)
      ELSE 
          cycle   !added to jump over, Liwang Ma
	      G1=X2
	      DO IT=1,10
              Z1=EXP(X4-X5*G1)
	          Z2=G1+Z1
              FU=G1/Z2-.9
	          IF(ABS(FU)<1.E-5)EXIT
	          DFDU=Z1*(1.+X5*G1)/(Z2*Z2)
		      X1=FU/DFDU
		      X2=ABS(X1)
	          X3=.1*G1
	          IF(X2>X3)X1=X3*X1/X2
		      G1=G1-X1
          END DO
	      IF(IT>10)WRITE(KW(1),5)
          X3=G1
      END IF
      IF(SCLM(31)>0.)X3=MIN(X3,SCLM(31))
      PPLA(JJK,IHU(JJK,J),J)=DMLA(JJK)*X3/(X3+EXP(X4-X5*X3))
      POP(JJK,IHU(JJK,J),J)=X3
	  PLAX(JJK)=PPLA(JJK,IHU(JJK,J),J)
      PPL0(JJK,J)=POP(JJK,IHU(JJK,J),J)
      XLAI(JJK,J)=FNPP(PPL0(JJK,J))
      DMLX(JJK,J)=XLAI(JJK,J)
      ENDIF
      END DO
      END IF
      IF (NBASAL==1 .AND. PBAS(JJK,J)>0 .AND. NREPLACE==0) THEN
      X3=PBAS(JJK,J)
      IF(SCLM(31)>0.)X3=MIN(X3,SCLM(31))
      PPLA(JJK,IHU(JJK,J),J)=DMLA(JJK)*X3/(X3+EXP(X4-X5*X3))
      POP(JJK,IHU(JJK,J),J)=X3
	  PLAX(JJK)=PPLA(JJK,IHU(JJK,J),J)
      PPL0(JJK,J)=POP(JJK,IHU(JJK,J),J)
      XLAI(JJK,J)=FNPP(PPL0(JJK,J))
      DMLX(JJK,J)=XLAI(JJK,J)
      ENDIF
      END DO
      END DO
      END IF
	  RETURN
    5 FORMAT(5X,'!!!!! PLANT POP DID NOT CONVERGE')
      END