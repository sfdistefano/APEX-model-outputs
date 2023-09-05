      SUBROUTINE WGN
!     APEX1905
!     THIS IS THE MASTER WEATHER GENERATING SUBPROGRAM. IT CALCULATES
!     MEAN DAILY MAXIMUM AND MINUMUM AIR TEMPERATURE, SOLAR RADIATION
!     AND RELATIVE HUMIDITY BASED ON WET OR DRY DAY. ALSO CALCULATES THE
!     STANDARD NORMAL DEVIATES.
      USE PARM 
      DIMENSION A(3,3),B(3,3),XX(3),SNDV(3)
      DATA A/.594,.454,-.004,.076,.261,-.037,-.018,-.129,.222/,B/.767,&
      .304,.274,0.,.692,-.33,0.,0.,.873/
      XXX=.5*(OBMX(IWI,MO)-OBMN(IWI,MO))
      III=1
      Z2=WFT(IWI,MO)
      YY=.9*Z2
      TXXM=OBMX(IWI,MO)+XXX*Z2
      RHM=(RH(IWI,MO)-YY)/(1.-YY)
      IF(RHM<.05)RHM=.5*RH(IWI,MO)
      RM=OBSL(IWI,MO)/(1.-.25*Z2)
      IF(RFV0(1)>2.)THEN
          TXXM=TXXM-XXX
          RM=.5*RM
          RHM=RHM*.1+.9
      END IF
      DO I=1,3
          V2(ISA)=AUNIF(IDG(2))
          SNDV(I)=ADSTN(VV1(ISA),V2(ISA))
          VV1(ISA)=V2(ISA)
      END DO
      DO I=1,3
          WX(I)=0.
          XX(I)=0.
          DO J=1,3
              WX(I)=WX(I)+B(I,J)*SNDV(J)
              XX(I)=XX(I)+A(I,J)*XIM(J)
          END DO
      END DO        
      DO I=1,3
          WX(I)=WX(I)+XX(I)
          XIM(I)=WX(I)
      END DO
      RETURN
      END