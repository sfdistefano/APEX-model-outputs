      SUBROUTINE DUSTDST
!     APEX1905
!     THIS SUBPROGRAM SETTLES AND DISTRIBUTES PARTICULATE MATTER(DUST) 
!     FROM FEEDYARDS TO SURROUNDING SUBAREAS.
      USE PARM
      DIMENSION NX(100)
      DIMENSION DPY(100),WDRX(16)
      DATA A90/1.5708/,A180/3.1416/,A270/4.7124/,A360/6.2832/,WDRX/.1963&
      ,.589,.9817,1.374,1.767,2.16,2.553,2.945,3.338,3.731,4.123,4.516,&
      4.909,5.301,5.694,6.087/
      IF(IDFH(ISA)==0)RETURN
!     IF(IFED(IDFH(ISA),IDON(ISA))/=ISA)RETURN
      TTH=TAN(THW)
      GRDX=.0005*GRDL(ISA)
      DO I=1,16
          IF(THW<WDRX(I))EXIT
      END DO
      IF(I>16)I=1
      NWDR(I)=NWDR(I)+1
      CALL DUSTEM
      IF(DEMR<1.E-5)RETURN
      SUM=0.
      NDX=0
      DO I=1,MSA
          DPY(I)=0.
          NX(I)=I
          IF(I==ISA)CYCLE
          JDR=10
          ! COMPUTE BASE Y COORDINATE--X COORDINATE = WIND DIRECTION (THW)
          YB=(XCT(ISA)-XCT(I))*TTH+YCT(ISA)
          ! DETERMINE WHICH SUBAREAS ARE DOWN WIND FROM FEEDYARD
          IF(THW<A90.OR.THW>A270)THEN
              ! CASE 1 OR 2 (N WIND)
              IF(YCT(I)>YB)CYCLE
              Y0=YCT(I)-YCT(ISA)
              X0=XCT(I)-XCT(ISA)
              DIST=SQRT(X0*X0+Y0*Y0)
              IF(THW<A90)THEN
                  ! CASE 1 (NE WIND) 
                  IF(Y0>0.)THEN
                      ! CASE 1.1
                      ALF=ATAN(Y0/(X0+.001))
                      BTA=THW+ALF
                      JDR=9
                  ELSE
                      IF(X0>0.)THEN
                          ! CASE 1.4
                          ALF=ATAN(X0/(Y0+.001))
                          BTA=A90-THW+ALF
                          JDR=9
                      ELSE 
                          IF(Y0>X0/TTH)THEN
                              ! CASE 1.2
                              ALF=ATAN(Y0/(X0+.001))
                              BTA=A90-THW-ALF
                          ELSE    
                              ! CASE 1.3
                              ALF=ATAN(X0/(Y0+.001))
                              BTA=THW-ALF
                          END IF   
                      END IF    
                  END IF    
              ELSE    
                  ! CASE 2 (NW WIND)      
                  IF(X0<0.)THEN
                      ! CASE 2.1
                      ALF=ATAN(X0/(Y0+.001))
                      BTA=THW-A270-ALF
                      JDR=9
                  ELSE
                      IF(Y0>0.)THEN
                          ! CASE 2.4
                          ALF=ATAN(Y0/(X0+.001))
                          BTA=A360-THW-ALF
                          JDR=9
                      ELSE
                          IF(Y0<X0*TTH)THEN
                              ! CASE 2.2
                              ALF=ATAN(X0/(Y0+.001))
                              BTA=A360-THW+ALF
                          ELSE    
                              ! CASE 2.3
                              ALF=ATAN(Y0/(X0+.001))
                              BTA=THW-A270+ALF
                          END IF
                      END IF    
                  END IF    
              END IF    
          ELSE    
              IF(YCT(I)<YB)CYCLE
              ! CASE 3 OR 4 (S WIND)
              Y0=YCT(I)-YCT(ISA)
              X0=XCT(I)-XCT(ISA)
              DIST=SQRT(X0*X0+Y0*Y0)
              IF(THW<A180)THEN
                  ! CASE 3 (SE WIND)
                  IF(Y0<0.)THEN
                      ! CASE 3.1
                      ALF=ATAN(Y0/(X0+.001))
                      BTA=A180-THW-ALF
                      JDR=9
                  ELSE
                      IF(X0>0.)THEN
                          ! CASE 3.4
                          ALF=ATAN(X0/(Y0+.001))
                          BTA=THW-A90-ALF
                          JDR=9
                      ELSE    
                          IF(Y0<X0*TTH)THEN
                              ! CASE 3.2
                              ALF=ATAN(Y0/(X0+.001))
                              BTA=THW-A90+ALF
                          ELSE    
                              ! CASE 3.3
                              ALF=ATAN(X0/(Y0+.001))
                              BTA=A180-THW+ALF
                          END IF
                      END IF
                  END IF    
              ELSE    
                  ! CASE 4 (SW WIND)
                  IF(X0<0.)THEN
                      ! CASE 4.1
                      ALF=ATAN(X0/(Y0+.001))
                      BTA=A270-THW+ALF
                      JDR=9
                  ELSE
                      IF(Y0<0.)THEN
                          ! CASE 4.4
                          ALF=ATAN(Y0/(X0+.001))
                          BTA=THW-A180+ALF
                          JDR=9
                      ELSE    
                          IF(Y0>X0/TTH)THEN
                              ! CASE 4.2
                              ALF=ATAN(X0/(Y0+.001))
                              BTA=THW-A180-ALF
                          ELSE    
                              ! CASE 4.3
                              ALF=ATAN(Y0/(X0+.001))
                              BTA=A270-THW-ALF
                          END IF
                      END IF
                  END IF    
              END IF
          END IF
          ! COMPUTE NEW X AND Y COORDINATES ON BASE (XB,YB)
          IF(JDR==9)THEN
              XT=ABS(DIST*COS(BTA))
              YT=ABS(DIST*SIN(BTA))
          ELSE    
              XT=ABS(DIST*SIN(BTA))
              YT=ABS(DIST*COS(BTA))
          END IF    
          IF(XT<GRDX)THEN
              FX=1.
          ELSE
              X1=YT/XT
              IF(X1<.01)CYCLE
              X2=ATAN(X1)
              FX=(SIN(X2))**PRMT(112)
          END IF
          TRT=YT/(3.6*U10(IRF(ISA)))
          X1=PRMT(111)*TRT*10.
          IF(X1>10.)CYCLE
          FY=EXP(-X1)
          DPY(I)=FX*FY*WSA(I)
          SUM=SUM+DPY(I)
          NDX=NDX+1
          ! WRITE(KW(1),1)I,NBSA(I),IY,MO,KDA,THW,ALF,BTA,X0,Y0,XT,YT,DEMR,Q,
          ! 1U10(IRF(ISA)),TRT,CIN,CW,FX,FY,DPY(I)
      END DO
      IF(SUM>0.)THEN
          B1=DEMR/SUM
      ELSE
          B1=1.
      END IF   
      TOT=0.
      DO I=1,MSA
          DPY(I)=DPY(I)*B1
          TOT=TOT+DPY(I)
          PM10(I)=PM10(I)+DPY(I)
          SMM(94,MO,I)=SMM(94,MO,I)+DPY(I)
      END DO
      CALL ASORT1(DPY,NX,MSA)
      IF(KFL(20)>0)WRITE(KW(20),32)IY,IYR,MO,KDA,THW,U10(IRF(ISA)),DEMR,TOT
      SUM=0.
      DO I=1,NDX
          X1=DPY(NX(I))/TOT
          SUM=SUM+X1
          IF(SUM>.999)CYCLE
          ! WRITE(KW(1),30)I,NX(I),NBSA(NX(I)),DPY(NX(I)),X1,SUM
          IF(KFL(20)>0)WRITE(KW(20),31)I,NX(I),NBSA(NX(I)),DPY(NX(I)),X1,&
          SUM
      END DO
      ! WRITE(KW(1),29)DEMR,TOT
      RETURN
!   1 FORMAT(1X,3I4,2I2,20E13.5)
!  29 FORMAT(1X,2E13.5)
!  30 FORMAT(1X,3I4,3E13.5)
   31 FORMAT(5X,I8,1X,2I8,F10.2,2F10.4)
   32 FORMAT(2X,'YR#=',I4,1X,'Y M D=',I4,2I2,1X,'W DIR=',F5.1,' RAD',2X,&
      'W SPD=',F6.2,' m/s',2X,'DEMR=',F6.0,' KG',2X,'DDPR=',F6.0,' KG')
      END