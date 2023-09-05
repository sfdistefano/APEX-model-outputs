      SUBROUTINE CPTBL
!     APEX1905
!     THIS SUBPROGRAM READS CROP TABLE TO DETERMINE CROP PARAMETERS
      USE PARM
      DIMENSION YTP(65)
      L=1
      IF(LC>0)THEN
          DO L=1,LC
              IF(KDC(L)==JX(6))EXIT
          END DO
      END IF
      IF(L>LC)THEN
          LC=LC+1
          KDC1(JX(6))=LC
          KDC(LC)=JX(6)
          XMTU(LC)=JX(7)
          ! READ CROP DATA
          ! CPNM = NAMES OF CROPS IN CROP PARAMETER TABLE
          ! 1  WA   = BIOMASS/ENERGY(kg/ha/MJ)(FOR CO2 = 330 ppm)
          ! 2  HI   = HARVEST INDEX(CROP YIELD/ABOVE GROUND BIOMASS)
          ! 3  TB   = OPTIMAL TEMP FOR PLANT GROWTH(c)
          ! 4  TG   = MIN TEMP FOR PLANT GROWTH(c)
          ! 5  DMLA = MAX LEAF AREA INDEX
          ! 6  DLAI = FRACTION OF GROWING SEASON WHEN LEAF AREA STARTS DECLINING
          ! 7  DLAP(= LAI DEVELOPMENT PARMS--NUMBERS BEFORE DECIMAL = %
          ! 8   1,2)  OF GROWING SEASON .  NUMBERS AFTER DDECIMAL = FRACTION OF
          !           DMLA AT GIVEN %.
          ! 9  RLAD = LAI DECLINE RATE FACTOR.
          !10  RBMD = WA  DECLINE RATE FACTOR.
          !11  ALT  = ALUMINUM TOLERANCE INDEX(1-5)1=SENSITIVE, 5=TOLERANT
          !12  GSI  = MAX STOMATAL CONDUCTANCE(DROUGTH TOLERANT PLANTS HAVE
          !           LOW VALUES.)
          !13  CAF  = CRITICAL AERATION FACTOR(SW/POR > CAF CAUSES AIR STRESS)
          !14  SDW  = SEEDING RATE(kg/ha)
          !15  HMX  = MAX CROP HEIGHT(m)
          !16  RDMX = MAX ROOT DEPTH(m)
          !17  WAC2 = NUMBER BEFORE DECIMAL = CO2 CONC IN FUTURE ATMOSPHERE
          !           (ppm).  NUMBER AFTER DECIMAL = RESULTING WA VALUE.
          !18  CNY  = N FRACTION OF YIELD(KG/KG)
          !19  CPY  = P FRACTION OF YIELD(KG/KG)
          !20  CKY  = K FRACTION OF YIELD(KG/KG)
          !21  WSYF = LOWER LIMIT OF HARVEST INDEX
          !22  PST  = PEST(INSECTS,WEEDS,AND DISEASE)FACTOR(0-1)
          !23  CSTS = SEED COST($/KG)
          !24  PRYG = PRICE FOR SEED YIELD ($/t)
          !25  PRYF = PRICE FOR FORAGE YIELD ($/t)
          !26  WCY  = FRACTION WATER IN YIELD.
          !27-29BN   = N FRACTION IN PLANT WHEN GROWTH IS 0.,.5,1.0
          !30-32BP   = P FRACTION IN PLANT WHEN GROWTH IS 0.,.5,1.
          !33-35BK   = K FRACTION IN PLANT WHEN GROWTH IS 0.,.5,1.
          !36-38BW   = WIND EROSION FACTORS FOR STANDING LIVE, STANDING DEAD,
          !            AND FLAT RESIDUE
          !39  IDC  = CROP ID NUMBERS. USED TO PLACE CROPS IN CATEGORIES
          !           1 FOR WARM SEASON ANNUAL LEGUME
          !           2 FOR COLD SEASON ANNUAL LEGUME
          !           3 FOR PERENNIAL LEGUME
          !           4 FOR WARM SEASON ANNUAL
          !           5 FOR COLD SEASON ANNUAL
          !           6 FOR PERENNIAL
          !           7 FOR EVERGREEN TREES
          !           8 FOR DECIDUOUS TREES
          !           9 FOR COTTON
          !          10 FOR LEGUME TREES
          !40/ FRST(= FROST DAMAGE PARMS--NUMBERS BEFORE DECIMAL = MIN TEMP(DEG C
          !41   1,2) NUMBERS AFTER DECIMAL = FRACTION YLD LOST WHEN GIVEN MIN TE
          !          IS EXPERIENCED.
          !42  WAVP = PARM RELATING VAPOR PRESSURE DEFFICIT TO WA
          !43  VPTH = THRESHOLD VPD (KPA)(F=1.
          !44  VPD2 = NUMBER BEFORE DECIMAL = VPD VALUE (KPA).  NUMBER AFTER
          !           DECIMAL = F2 < 1.
          !45  RWPC1= ROOT WEIGHT/BIOMASS PARTITIONING COEF
          !46  RWPC2= ROOT WEIGHT/BIOMASS PARTITIONING COEF
          !47  GMHU = HEAT UNITS REQUIRED FOR GERMINATION
          !48  PPLP1= PLANT POP PARM--NUMBER BEFORE DECIMAL = # PLANTS
          !           NO AFTER DEC = FRACTION OF MAX LAI
          !49  PPLP2= 2ND POINT ON PLANT POP-LAI CURVE. PPLP1<PPLP2--PLANTS/M2
          !                                             PPLP1>PPLP2--PLANTS/HA
          !50  STX1 = YLD DECREASE/SALINITY INCREASE ((t/ha)/(mmho/cm))
          !51  STX2 = SALINITY THRESHOLD (mmho/cm)
          !52  BLG1 = LIGNIN FRACTION IN PLANT AT .5 MATURITY
          !53  BLG2 = LIGNIN FRACTION IN PLANT AT MATURITY
          !54  WUB  = WATER USE CONVERSION TO BIOMASS(t/mm)
          !55  FTO  = FRACTION TURNOUT (COTTON LINT/STRIPPER YLD)
          !56  FLT  = FRACTION LINT (COTTON LINT/PICKER YLD)
          !57  EXTC = EXTINCTION COEFFICIENT  
          !58  
          !59  TDNX = MAX TOTAL DIGESTABLE NUTRIENTS (TDN) %
          !60  TDNN = MIN TDN %
          !61  ANTQX= MAXIMUM ANTI TOLERANCE FOR GRAZERS
          !62  ANTQN= MINIMUM ANTI TOLERANCE FOR GRAZERS
          !63  DETIN= (0 INDICATES NORMAL FUNCTION.  GREATER THAN 0 IS THE MAX ACHIEVED HUI)
		  !64  RDEAT= DAILY RATE OF ROOT DEATH
          READ(KR(4),'()')
          READ(KR(4),'()')
          J2=-10
          DO WHILE(J2/=JX(6))
              READ(KR(4),3420,IOSTAT=NFL)J2,CPNM(LC),(YTP(L),L=1,65)
              IF(NFL/=0)THEN
                  WRITE(*,*)'CROP NO = ',JX(6),' NOT IN CROP LIST FILE     &
                  SAID = ',NBSA(ISA)
                  PAUSE
	              STOP       
	          END IF
          END DO
          WA(LC)=YTP(1)
          WAI(LC)=WA(LC)
          HI(LC)=YTP(2)
          TOPC(LC)=YTP(3)
          TBSC(LC)=YTP(4)
          DMLA(LC)=YTP(5)
          DLAI(LC)=YTP(6)
          DLAP(1,LC)=YTP(7)
          DLAP(2,LC)=YTP(8)
          RLAD(LC)=YTP(9)
          RBMD(LC)=YTP(10)
          ALT(LC)=YTP(11)
          GSI(LC)=YTP(12)
          CAF(LC)=YTP(13)
          SDW(LC)=YTP(14)
          HMX(LC)=YTP(15)
          RDMX(LC)=YTP(16)
          WAC2(2,LC)=YTP(17)      
          CNY(LC)=YTP(18)
          CPY(LC)=YTP(19)
          CKY(LC)=YTP(20)
          WSYF(LC)=YTP(21)
          PST(LC)=YTP(22)
          CSTS(LC)=YTP(23)
          PRYG(LC)=YTP(24)
          PRYF(LC)=YTP(25)
          WCY(LC)=YTP(26)
          BN(1,LC)=YTP(27)
          BN(2,LC)=YTP(28)
          BN(3,LC)=YTP(29)
          BP(1,LC)=YTP(30)
          BP(2,LC)=YTP(31)
          BP(3,LC)=YTP(32)
          BK(1,LC)=YTP(33)
          BK(2,LC)=YTP(34)
          BK(3,LC)=YTP(35)
          BWN(1,LC)=YTP(36)
          BWN(2,LC)=YTP(37)
          BWN(3,LC)=YTP(38)
          IDC(LC)=YTP(39)
          FRST(1,LC)=YTP(40)
          FRST(2,LC)=YTP(41)
          WAVP(LC)=YTP(42)
          VPTH(LC)=YTP(43)
          VPD2(LC)=YTP(44)
          RWPC(1,LC)=YTP(45)
          RWPC(2,LC)=YTP(46)
          GMHU(LC)=YTP(47)
          PPLP(1,LC)=YTP(48)
          PPLP(2,LC)=YTP(49)
          STX(1,LC)=YTP(50)
          STX(2,LC)=YTP(51)
          BLG(1,LC)=YTP(52)
          BLG(2,LC)=YTP(53)
	      FTO(LC)=YTP(55)
	      FLT(LC)=YTP(56)
	      EXTC(LC)=YTP(57)
          TDNX(LC)=.01*YTP(59)
	      TDNN(LC)=.01*YTP(60)
          ANTQX(LC)=YTP(61)
          ANTQN(LC)=YTP(62)
		  DETIN(LC)=YTP(63)
          RDEAT(LC)=YTP(64)
          BASAL(LC)=YTP(65)
          !EXTC(LC)=YTP(63)
          REWIND KR(4)
          X1=.6*BK(1,LC)
          IF(BK(2,LC)>X1)BK(2,LC)=X1
          X1=.8*BK(2,LC)
          IF(BK(3,LC)>X1)BK(3,LC)=X1
      END IF    
      JJK=KDC1(JX(6))
      CNT(JJK,ISA)=BN(1,JJK)
      IHU(JJK,ISA)=IHU(JJK,ISA)+1
      PHU(JJK,IHU(JJK,ISA),ISA)=OPV(1)
	  BLAI(JJK,ISA)=.01*DMLA(JJK)							 
      IF(XMTU(JJK)>0)PHU(JJK,IHU(JJK,ISA),ISA)=CAHU(1,365,TBSC(JJK),0)&
      *XMTU(JJK)
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
	  IF(OPV(5)>0.)THEN
          X3=OPV(5)
      ELSE          
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
      PPLA(JJK,IHU(JJK,ISA),ISA)=DMLA(JJK)*X3/(X3+EXP(X4-X5*X3))
      POP(JJK,IHU(JJK,ISA),ISA)=X3
	  PHUX(JJK)=PHU(JJK,IHU(JJK,ISA),ISA)
	  PLAX(JJK)=PPLA(JJK,IHU(JJK,ISA),ISA)
	  POPX(JJK)=X3
      TDNF(JJK,ISA)=.5*(TDNX(JJK)+TDNN(JJK))
	  RETURN
    5 FORMAT(5X,'!!!!! PLANT POP DID NOT CONVERGE')
 3420 FORMAT(1X,I4,1X,A4,65F8.0)
      END