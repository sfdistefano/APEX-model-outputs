      SUBROUTINE HIRG(AVIR,EFD,ZX,JRT,IRX,IRY)
!     APEX1307
!     THIS SUBPROGRAM IS USED TO SIMULATE AUTOMATIC OR USER SPECIFIED
!     IRRIGATION APPLICATIONS.  COMPUTES THE AMOUNT OF IRRIGATION WATER
!     NEEDED TO BRING THE ROOT ZONE WATER CONTENT TO FIELD CAPACITY, FOR
!     AUTOMATIC OPTION.  USER SPECIFIED AMOUNT IS APPLIED FOR
!     MANUAL OPTION.  EROSION AND RUNOFF ARE ESTIMATED.
      USE PARM
      JRT=0
      IF(VIRT(ISA)>=VIMX(ISA).OR.NII(ISA)<IRI(ISA))THEN
          AVIR=0.
          JRT=1
          RETURN
      END IF
      WSAX1=10.*WSA(ISA)
      IF(IRY==0)THEN
          ! AUTO
          SELECT CASE(IAC(ISA))
              CASE(0)
                  ! FLEXIBLE
                  XX=0.
                  DO J=1,NBSL(ISA)
                      ISL=LID(J,ISA)
                      XX=XX+FC(ISL,ISA)-SWST(ISL,ISA)
                  END DO
                  XX=FIRG(ISA)*XX/((1.-EFI(ISA))*EFD)
                  X4=MIN(VIMX(ISA)-VIRT(ISA),XX,ARMX(ISA))
              CASE(1)
                  ! RIGID
                  X4=MIN(VIMX(ISA)-VIRT(ISA),ARMX(ISA))
              CASE(2)    
                  !Irrigate when ponding depth goes below the minimum depth in the paddy jaehak jeong 2016, revised 2020
                  IF(PADDY_STO(1,ISA)<PADDY_HMIN(ISA)) THEN
                     X4=AVIR
                  ELSE
                     X4=0
                  ENDIF
          END SELECT
      ELSE    
          ! MANUAL OR LAGOON
          IF(IAC(ISA)==0)THEN
              ! FLEXIBLE
              XX=0.
              DO J=1,NBSL(ISA)
                  ISL=LID(J,ISA)
                  XX=XX+FC(ISL,ISA)-SWST(ISL,ISA)
              END DO
              XX=FIRG(ISA)*XX/((1.-EFI(ISA))*EFD)
              X4=MIN(VIMX(ISA)-VIRT(ISA),XX,AVIR)
          ELSE
              ! RIGID
              X4=AVIR
          END IF
      END IF 
      IF(IRRS(ISA)>0)THEN
          I1=NISA(IRRS(ISA))
          X3=X4*WSAX1
          RSM3=MIN(RSV(I1),X3)
          X4=RSM3/WSAX1
      ELSE           
          IF(IRRW(ISA)>0)THEN
              I1=NISA(IRRW(ISA))
              X3=X4*WSAX1
              GWM3=MIN(GWST(I1),X3)
              X4=GWM3/WSAX1
          END IF           
      END IF    
      AVIR=X4*EFD
	  X3=COIR*X4
      IF(AVIR>ARMN(ISA))THEN      
          NII(ISA)=0
          REPI(ISA)=AVIR/24.
	      XEF=WSAX1*(X4-AVIR)
          SMM(110,MO,ISA)=SMM(110,MO,ISA)+XEF
          VAR(110,ISA)=XEF
          X1=RZSW(ISA)-PAW(ISA)
          IF((NOP>0.OR.NBSA(ISA)==ISAP).AND.KFL(1)>0)WRITE(KW(1),9)ISA,&
          NBSA(ISA),IYR,MO,KDA,TIL(IRX),X4,WS(ISA),WTN,X1,BIR(ISA),&
          EFI(ISA),FIRG(ISA),XHSM(ISA)
          VIR(JJK,ISA)=VIR(JJK,ISA)+X4
          VIRT(ISA)=VIRT(ISA)+X4
	      COST(ISA)=COST(ISA)+X3
          IF(VIRT(ISA)>0.)THEN
              IF(IRY==0)THEN
                  X1=COTL(IRX)
                  X2=X1-COOP(IRX)
                  COST(ISA)=COST(ISA)+X1
                  CSFX=CSFX+X2
              END IF
              IF(KFL(31)>0)THEN
                  WRITE(KW(31),14)ISA,NBSA(ISA),IYR,MO,KDA,TIL(IRX),KDC(JJK),IHC&
                  (IRX),NBE(IRX),X3,X3,X4
                  IF(IRY==0)WRITE(KW(31),50)ISA,NBSA(ISA),IYR,MO,KDA,TIL(IRX),KDC&
                  (JJK),IHC(IRX),NBE(IRX),NBT(IRX),X1,X2,FULU(IRX)
              END IF
          END IF
          VSLT(ISA)=.01*X4*CSLT
          SMM(129,MO,ISA)=SMM(129,MO,ISA)+VSLT(ISA)
          VAR(129,ISA)=VSLT(ISA)
          IF(IRRS(ISA)>0)THEN
              VAR(146,I1)=RSM3
              SMM(146,MO,I1)=SMM(146,MO,I1)+RSM3
              RSV(I1)=RSV(I1)-RSM3
          ELSE
              IF(IRRW(ISA)>0)THEN
                  VAR(156,I1)=GWM3
                  SMM(156,MO,I1)=SMM(156,MO,I1)+GWM3
                  GWST(I1)=GWST(I1)-GWM3
              END IF
          END IF
          X4=X4*WSAX1
          SMM(18,MO,ISA)=SMM(18,MO,ISA)+X4
          VAR(18,ISA)=X4
      ELSE
          AVIR=0.
      END IF
      JRT=0
      IF(IRR(ISA)/=5)RETURN
	  DO J=1,NBSL(ISA)
          I=LID(J,ISA)
          IF(ZX<Z(I,ISA))EXIT
      END DO
      SWST(I,ISA)=SWST(I,ISA)+AVIR
      AVIR=0.
	  REPI(ISA)=0.
      RETURN
    9 FORMAT(1X,2I8,1X,3I4,2X,A8,2X,'VOL=',F5.0,' mm',2X,'WS=',F5.2,2X,&
      'WTN=',F6.0,'KPA',2X,'PWDF=',F6.0,' mm',2X,'TRGR=',F7.2,2X,'Q/VIR='&
      ,F5.2,2X,'FIRG=',F5.2,2X,'HUSC=',F5.2)      
   14 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,8X,I6,6X,2I4,4X,F10.2,10X,2F10.2)
   50 FORMAT(1X,2I8,1X,I4,2I2,2X,A8,8X,I6,6X,3I4,2F10.2,20X,F10.2)
      END