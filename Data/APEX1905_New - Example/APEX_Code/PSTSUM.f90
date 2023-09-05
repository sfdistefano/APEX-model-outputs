      SUBROUTINE PSTSUM(XWS,XYR,NX,ISX)
!     APEX1905
!     THIS SUBPROGRAM OUTPUTS THE PESTICIDE SUMMARY TABLES TO THE .OUT FILE.
      USE PARM
      DIMENSION NX(100),NXX(3,5),NYY(3,5),NY(3)
      DIMENSION XTX(100),XTY(100),AD1(NDP),AD2(NDP)
      DIMENSION PSSM(3,50)
      IF(ISX<0)THEN
          NY(1)=1
          NY(2)=.1*XYR+1.5
          NY(3)=.5*XYR+1.5
          DO K=1,NDP
              WRITE(KW(1),462)PSTN(K)
              DO I=1,5
                  DO J=1,NBYR
                      XTY(J)=APY(I,K,J,ISX)
                      IF(XTY(J)<1.E-4)THEN
                          APY(I,K,J,ISX)=0.
                          AYB(I,K,J,ISX)=0.
                      END IF
                      XTX(J)=APQ(I,K,J,ISX)
                      IF(XTX(J)>1.E-4)CYCLE
                      APQ(I,K,J,ISX)=0.
                      AQB(I,K,J,ISX)=0.
                  END DO
                  CALL ASORT1(XTX,NX,NBYR)
                  NXX(1,I)=NX(NY(1))
                  NXX(2,I)=NX(NY(2))
                  NXX(3,I)=NX(NY(3))
                  CALL ASORT1(XTY,NX,NBYR)
                  NYY(1,I)=NX(NY(1))
                  NYY(2,I)=NX(NY(2))
                  NYY(3,I)=NX(NY(3))
              END DO
              DO N2=1,3
                  !     PRINTOUT PESTICIDE FREQ SUMMARY
                  WRITE(KW(1),463)(APQ(I,K,NXX(N2,I),ISX),I=1,5)
                  WRITE(KW(1),464)(AQB(I,K,NXX(N2,I),ISX),I=1,5)
                  WRITE(KW(1),465)(APY(I,K,NYY(N2,I),ISX),I=1,5)
                  WRITE(KW(1),466)(AYB(I,K,NYY(N2,I),ISX),I=1,5)
                  IF(N2==1)WRITE(KW(1),474)
                  IF(N2==2)WRITE(KW(1),473)
              END DO
          END DO
      END IF
      WRITE(KW(1),3940)
      IF(KFL(30)>0)THEN
          IF(ISX==NCMD)THEN
              WRITE(KW(30),'(//T10,A)')'WATERSHED OUTLET'
          ELSE
              WRITE(KW(30),2232)ISX,NBSA(ISX)
          END IF
      END IF
      IF(ISX==NCMD)THEN 
          PSSM=0.
		  DO I=1,MSA
			  RTO=WSA(I)/RWSA(NCMD)
			  DO K=1,NDP
				  PSSM(1,K)=PSSM(1,K)+PFOL(K,I)*RTO
				  PSSM(2,K)=PSSM(2,K)+GWPS(K,I)*RTO
				  SUM=0.
				  DO J=1,NBSL(I) 
					  SUM=SUM+PSTZ(K,J,I)
				  END DO
				  PSSM(3,K)=PSSM(3,K)+SUM*RTO
			  END DO
		  END DO
		  II=ISX
      ELSE
          II=IDOA(ISX)		  
	  END IF
      DO K=1,NDP
          AD1(K)=0.
          AD2(K)=0.
          DO I=1,12
              AD1(K)=AD1(K)+SMRP(3,K,I)
              AD2(K)=AD2(K)+SMRP(4,K,I)
          END DO
          DO L=1,7
              SMAP(L,K,II)=SMAP(L,K,II)/XWS
          END DO
          DO L=10,13
			  SMAP(L,K,II)=SMAP(L,K,II)/XWS
		  END DO
		  IF(ISX==NCMD)THEN
			  WRITE(KW(1),'(/1X,A,A16)')'-----PESTICIDE BALANCE(g) ',&
              PSTN(K)
              AD1(K)=AD1(K)-SMAP(4,K,ISX)-SMAP(11,K,ISX)
              DF=SMAP(1,K,ISX)-AD1(K)-SMAP(4,K,ISX)-AD2(K)-SMAP(6,K,ISX)&
              -SMAP(7,K,ISX)-SMAP(11,K,ISX)-SMAP(12,K,ISX)-PSSM(1,K)-&
              PSSM(2,K)-PSSM(3,K)
              PER=DF/(SMAP(1,K,ISX)+1.E-10)
              WRITE(KW(1),472)PER,DF,SMAP(1,K,ISX),AD1(K),SMAP(4,K,ISX),&
              AD2(K),SMAP(6,K,ISX),SMAP(7,K,ISX),SMAP(11,K,ISX),SMAP(12,K,&
              ISX),PSSM(1,K),PSSM(2,K),PSSM(3,K)
          END IF 
		  DO L=1,7
			  SMAP(L,K,II)=SMAP(L,K,II)/XYR
		  END DO
		  DO L=10,13
		      SMAP(L,K,II)=SMAP(L,K,II)/XYR
		  END DO
          AD1(K)=AD1(K)/XYR
          AD2(K)=AD2(K)/XYR
      END DO 
      IF(KFL(30)>0)WRITE(KW(30),3910)(PSTN(K),K=1,NDP)
      I1=0
      K1=0
      N1=0
      DO WHILE(N1<NDP)
          I1=I1+10
          N1=MIN(I1,NDP)
          K2=K1+1
          N2=MIN(10,NDP-K1)
    !     PRINTOUT PESTICIDE SUMMARY MONTHLY
          WRITE(KW(1),3910)(PSTN(K),K=K2,N1)
          IF(ISX==NCMD)THEN
              WRITE(KW(1),470)HEDP(1),(SMAP(1,K,ISX),K=K2,N1)
              WRITE(KW(1),470)HEDP(2),(AD1(K),K=K2,N1)
              DO L=3,4
                  WRITE(KW(1),470)HEDP(L),(SMAP(L,K,ISX),K=K2,N1)
              END DO
              WRITE(KW(1),470)HEDP(5),(AD2(K),K=K2,N1)
          ELSE
              DO L=1,7
                  WRITE(KW(1),470)HEDP(L),(SMAP(L,K,II),K=K2,N1)
              END DO
          END IF
          DO L=10,13
              WRITE(KW(1),470)HEDP(L),(SMAP(L,K,II),K=K2,N1)
              IF(KFL(30)>0)WRITE(KW(30),470)HEDP(L),(SMAP(L,K,II),K=K2,N1)
          END DO
          K1=I1
      END DO
      RETURN 
  462 FORMAT(5X,A16,T35,'MAXIMUM')
  463 FORMAT(8X,'SOL  ',5E13.5)
  464 FORMAT(8X,'Q+SSF',5E13.5)
  465 FORMAT(8X,'ADSRB',5E13.5)
  466 FORMAT(8X,'SED Y',5E13.5)
  470 FORMAT(5X,A4,10E16.6)
  472 FORMAT(5X,'PER =',E13.6,2X,'DF  =',E13.6,2X,'PAPL=',E13.6,2X,&
      'PSRO=',E13.6,2X,'PSSF=',E13.6,2X,'PSED=',E13.6/5X,'PDGF=',E13.6,&
      2X,'PDGS=',E13.6,2X,'PRSF=',E13.6,2X,'PDPK=',E13.6,2X,'PFOL=',&
      E13.6,2X,'PGW =',E13.6/5X,'PSOL=',E13.6/)
  473 FORMAT(T35,'50 % EXCEED')
  474 FORMAT(T35,'10 % EXCEED')
 2232 FORMAT(//T10,'SA#=',I8,1X,'ID=',I8)  
 3910 FORMAT(/11X,10A16)
 3940 FORMAT(/1X,'-----AVE ANNUAL VALUES(g/ha)')
      END