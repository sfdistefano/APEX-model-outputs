      SUBROUTINE SOLIOP
!     APEX1905
!     THIS SUBPROGRAM OUTPUTS THE SOIL TABLE.
      USE PARM 
      WRITE(KW(1),2)
      DO I=1,NBSL(ISA)
          ISL=LID(I,ISA)
          WRITE(KW(1),40)LORG(ISL,ISA),Z(ISL,ISA),SOIL(8,ISL,ISA),SOIL(9,ISL&
          ,ISA),SSF(ISL,ISA),SOIL(12,ISL,ISA),SATC(ISL,ISA),HCL(ISL,ISA),BDP&
          (ISL,ISA),SOIL(13,ISL,ISA),SAN(ISL,ISA),SIL(ISL,ISA),CLA(ISL,ISA),&
          ROK(ISL,ISA),RSD(ISL,ISA)
      END DO
      WRITE(KW(1),41)UW(4,ISA),UW(1,ISA),UW(3,ISA),UW(2,ISA),TRSD(ISA)
      RETURN
    2 FORMAT(T18,'_______WATER CONTENT________   ___SAT COND__    __BULK &
      DEN___    ___________TEXTURE__________   CROP'/2X,'LAY',4X,&
      'DEPTH',4X,'POR',5X,'FC',6X,'WP',4X,'CURNT',4X,'VERT',4X,'HORZ',&
      4X,'33KPA  OV DRY',4X,'SAND',4X,'SILT',4X,'CLAY',4X,'ROCK',4X,&
      'RSD'/3X,'NO',5X,'(m)',4X,'------------(m/m)-----------',3X,&
      '---(mm/h)----',5X,'--(t/m3)---',5X,'------------(%)-------------'&
      ,2X,'(t/ha)')
   40 FORMAT(T2,I4,T7,F8.2,4F8.3,4F8.2,4F8.1,2F8.2)
   41 FORMAT(T2,'TOTAL',T15,4F8.3,T111,F8.2)
      END