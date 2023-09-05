      SUBROUTINE SOLIOC
!     APEX1905
!     THIS SUBPROGRAM OUTPUTS THE CHEMICAL SOIL TABLE.
      USE PARM 
      WRITE(KW(1),2)
      DO I=1,NBSL(ISA)
          ISL=LID(I,ISA)
          WRITE(KW(1),40)LORG(ISL,ISA),Z(ISL,ISA),PH(ISL,ISA),SMB(ISL,ISA),&
          CEC(ISL,ISA),ALS(ISL,ISA),CAC(ISL,ISA),PSP(ISL,ISA),SOIL(1,ISL,&
          ISA),SOIL(2,ISL,ISA),SOIL(3,ISL,ISA),SOIL(4,ISL,ISA),SOIL(5,ISL,&
          ISA),SOIL(6,ISL,ISA),SOIL(10,ISL,ISA),SOIL(11,ISL,ISA),&
          SOIL(15,ISL,ISA),SOIL(7,ISL,ISA),SOIL(16,ISL,ISA)
          X1=.001*ZOC(ISA)
      END DO
      WRITE(KW(1),41)ZPML(ISA),ZPMA(ISA),ZPMS(ISA),ZPO(ISA),ZNO3(ISA),&
      ZON(ISA),ZSK(ISA),ZEK(ISA),ZFK(ISA),X1,ZSLT(ISA)
      RETURN
    2 FORMAT(T27,'SUM',13X,'AL',12X,'P SORP',3X,'_______MIN P_______',5X&
      ,'ORG',5X,'NO3',5X,'ORG',T138,'ORG'/2X,'LAY',4X,'DEPTH',5X,'PH',5X,&
      'BASE',4X,'CEC',4X,'SAT',5X,'CACO3',4X,'RTO',4X,'LAB',5X,'ACT',5X,&
      'STB',6X,'P',7X,'N',7X,'N',4X,'SOLK',4X,'EXCK',4X,'FIXK',7X,'C',&
      6X,'SALT'/3X,'NO',5X,'(m)',12X,'-(cmol/kg)--',4X,'-----(%)-----',&
      11X,'---------------------------------(g/t)-----------------------&
      --------',3X,'(%)',5X,'(ppm)')
   40 FORMAT(T2,I4,T7,F8.2,5F8.1,F8.2,9F8.0,F8.2,F8.0)
   41 FORMAT(T2,'TOTAL',T63,9F8.0,F8.2,F8.0)
      END