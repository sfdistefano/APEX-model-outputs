      SUBROUTINE WDLYSTA
!     APEX1905      
!     THIS SUB PROGRAM ASSIGNS EACH SUBAREA A DAILY WEATHER FILE USING
!     IWTH IF GIVEN OR NEAREST STATION TO SUBAREA CENTROID IF IWTH IS
!     NOT GIVEN.  IT ALSO OPENS THE WEATHER FILES.
      USE PARM  
      CHARACTER(20)::ADUM
      IF(IWTH(ISA)==0)THEN
          XX=YCT(ISA)/CLT
          SIN1=SIN(XX)
          COS1=COS(XX)
          D0=1.E20
          DO 
              READ(KR(25),*,IOSTAT=NFL)JJ,ADUM,Y,X
              IF(NFL/=0)EXIT
	          RY=Y/CLT
	          XX=SIN1*SIN(RY)+COS1*COS(RY)*COS((X-XCT(ISA))/CLT)
              D=6378.8*ACOS(XX)
              IF(D>=D0)CYCLE
              D0=D
              IWTH(ISA)=JJ
              FWTH(NDWT)=ADUM
          END DO
      ELSE
          JJ=-1
          DO WHILE(JJ/=IWTH(ISA)) 
              READ(KR(25),*,IOSTAT=NFL)JJ,FWTH(NDWT)
              IF(NFL/=0)THEN
                  WRITE(*,*)'FWTH NO = ',IWTH(ISA),' NOT IN DAILY &
                  WEATHER LIST FILE     SAID = ',NBSA(ISA)
                  STOP
              END IF
	      END DO
	  END IF
      REWIND KR(25) 
	  IF(NDWT/=1)THEN
          DO L=1,MXW
              IF(NBW(L)==IWTH(ISA))GO TO 801
          END DO
      END IF
      MXW=NDWT
      L=MXW!+KND Jaehak 2019
      NBW(MXW)=IWTH(ISA)
	  IRF(ISA)=NDWT
	  CALL OPENV(KRST(L),FWTH(NDWT),IDIR)
	  IF(KFL(1)>0)WRITE(KW(1),'(T10,A,A80)')'DAILY WEATHER FILE = ',FWTH(NDWT)
      IF(NGN0<=0)THEN
          CALL WREAD(L,2)
      ELSE
          !KK=1 TO CONTINUE WEATHER FILE
          !  =2 TO REWIND WEATHER FILE
          KK=2
          CALL WREAD(L,KK)
      END IF
	  IYR=IYR0
      CALL ALPYR(IYR,NYD,LPYR)
	  NDWT=NDWT+1
      NWTH=MXW
      RETURN
  801 IRF(ISA)=L
      RETURN  
	  END