      SUBROUTINE TGRAZ
!     APEX1605
!     THIS SUBPROGRAM SIMULATES ANIMAL GRAZING USING CONCEPTS FROM PHY-GROW
!     GRAZING PREFERENCE COMPONENT.
      USE PARM
	  DIMENSION IPFD(KOMX(IOW),MNC),IPFL(KOMX(IOW),MNC)
      DIMENSION AD1(KOMX(IOW)),AD2(KOMX(IOW)),AD3(KOMX(IOW)),ADH2O(KOMX(IOW)),&
      ADNX(KOMX(IOW)),ADP(KOMX(IOW)),DESD(MNC),DESL(MNC),DESX(KOMX(IOW)),&
      DMDF(KOMX(IOW)),DEF(KOMX(IOW)),SPFD(KOMX(IOW)),SPGZ(KOMX(IOW))
      DIMENSION DSP(1000,MSA)
      DIMENSION RTD(KOMX(IOW),MNC),RTL(KOMX(IOW),MNC),TSD(MNC),TSL(MNC),YLSD(MNC)
      DIMENSION FSEL(3,KOMX(IOW)),YGZL(KOMX(IOW),MNC),YGZD(KOMX(IOW),MNC),YTP(3,KOMX(IOW))
      DIMENSION GZSD(3,KOMX(IOW),MNC),GZSL(3,KOMX(IOW),MNC),SPSL(3,KOMX(IOW),MNC),SPSD(3,KOMX(IOW),MNC)
      DATA FCA,FCO0/2*0./
      REAL::INKG, NEmcs, NEcs, NEhs, NEmpa, Slope, LoCTemp, InsTot, MEcs, PREV_PLANE, EMNR, HeatProd,&
	   SurfArea, HIDE, HairCoat, HairDepth, ExtIns, TissIns, EMN, EMNE, DOM2CP, Indicator, CPWG,PSWG,DWT
      LD1=LID(1,ISA)
      KT2=KT(ISA)
      KRT=0
	  KGZ=0
      LGZ=1
	  N1=KT(ISA)
      DMIN=0.
      DMDFN=0.
      GZNL=0.
      SPFD=0.
	  YGZL=0.
      YGZD=0.
	  INKG=0.
	  NEmcs=0.
	  NEcs=0.
	  NEhs=0.
	  NEmpa=0.
	  Slope=0.
	  LoCTemp=0.
	  InsTot=0.
	  MEcs=0.
	  Indicator=0.
	  DWT=0.
      IOW=IDON(ISA)
      KDG0=2*MSA+MSO       	 
      NN=NCP(IRO(ISA),ISA)	!number of forages
      TSLSD=0.
      TOT1=0.
      DO K=1,NN
          KKF=JE(K,ISA)
          IF(KKF==MNC)CYCLE
          TSLSD=TSLSD+STL(KKF,ISA)+STD(KKF,ISA)
          TOT1=TOT1+UN1(KKF,ISA)+STDN(KKF,ISA)
      END DO    
      ! DETERMINE STOCKING RATE, INITIAL HERD WEIGHT, POTENTIAL DAILY GRAZING RATE 
      DO KHD=1,NHRD(IOW) ! HERD LOOP
          IHD=KOW(KHD,IOW)
          IF(IHRD==0)THEN
              IF(RSTK(IRO(ISA),KT2,ISA)<1.E-10)CYCLE
              IGZR(IHD,ISA)=1
          ELSE    
              IF(IGZX(IHD,IOW)/=ISA.OR.IHD/=JGRZ(IRO(ISA),KT2,ISA))CYCLE
              IGZR(IHD,ISA)=1
          END IF    
          IF(RSTK(IRO(ISA),KT2,ISA)>0.)THEN
              GZNB(IHD,ISA)=WSA(ISA)/RSTK(IRO(ISA),KT2,ISA)
          ELSE    
              GZNB(IHD,ISA)=NHD(IHD,IOW)
          END IF    
          KKG=NGZ(IHD,ISA)
          IF(KKG==0)CYCLE
          IF(GZWX(IRO(ISA),KT2,ISA)>0.)THEN
              X1=GZWX(IRO(ISA),KT2,ISA)
          ELSE
              X1=GZWI(KKG)
          END IF    
          GZWT(IHD,IOW)=MAX(GZWT(IHD,IOW),X1)
          IF(INWT(IHD,IOW)==0)THEN
              WTBG(IHD,IOW)=GZWT(IHD,IOW)
              IGZB(IHD,IOW)=IDA
          END IF    
          INWT(IHD,IOW)=1
          IF(GZSM(IHD)<1.E-10)GZSM(IHD)=0. !Cody added this line to replace the line below.  The below line caused trouble when hay was not being fed and forage quantity was low.  The animals wouldn't lose weight because the grazed sum was set to full diet.
!          IF(GZSM(IHD)<1.E-10)GZSM(IHD)=.001*GZIN(KKG)*GZWT(IHD,IOW)*GZNB(IHD,ISA)/WSA(ISA)!GZSM = intake of grazed and supplemental forage.  that is to say, everything eaten.  Units are Mg/ha.  Here it is just calculated in case the value was zero.  
          DO K=1,NN	! FORAGE LOOP			
              KKF=JE(K,ISA)
              IF(KKF==MNC)CYCLE
              IF(JP(KKF,ISA)<=0)THEN
                  JP(KKF,ISA)=1														
                  NCR(KKF,ISA)=NCR(KKF,ISA)+1											
              END IF
          END DO ! FORAGE LOOP    
      END DO ! HERD LOOP
      YLN=0.
      YLP=0.
      SPSL=0.
	  SPSD=0.
	  TDMD=0.
      AD1=0.
      AD3=0.
      YMLK=0.
      DO KHD=1,NHRD(IOW) ! HERD DEMAND LOOP
          IHD=KOW(KHD,IOW)
          IF(IGZR(IHD,ISA)==0)CYCLE
          KKG=NGZ(IHD,ISA)
          IF(KKG==0)CYCLE
          RSTX=WSA(ISA)/GZNB(IHD,ISA) !RSTX = hectares/animal  (stocking density)
          INKG=1000.*GZSM(IHD)*RSTX !INKG is intake in kg/hd or kg/animal (hd stands for "head"; synonymous with animal)
		  Intake_kg_ha=1000.*GZSM(IHD) !Intake_kg_ha is intake in kg/ha
          DOM2CP=(ATDN(IHD)/1.05)/(ADN(IHD)*6.25)

		  SELECT CASE(IGZD(KKG))
              CASE(1) !GENERAL. Intake demand (DMI) is a percentage (GZIN) of body weight (GZWT)
                  DMI=GZIN(KKG)*GZWT(IHD,IOW) 														
                  DMDF(IHD)=.001*DMI/RSTX											
                  CH4=79.87+9.95*INKG
              CASE(2) !GENERAL using the forage value index (FVI).  This option is not recommended
                  DMIP=GZIN(KKG)*(.769+.231*FVI(IHD))
                  DMI=DMIP*GZWT(IHD,IOW)
                  DMDF(IHD)=.001*DMI/RSTX
                  CH4=79.87+9.95*INKG
              CASE(3) ! BEEF YEARLINGS.  Based on the NRC (2000).  Yearling intake demand and weight gain are based on the net energy system. Gain can be limited by protein or energy
			  !DE=Digestible energy; EM=Metabolizable energy; EMN=Net energy for maintenance; EGN=Net energy for gain; 
			  !SBW=Shrunk body weight; EMNR=Required net energy for maintenance; SRW=Standard reference weight
			  !FSBW=Final shrunk body weight; GZWM=Grazer weight at maturity; RE=Retained energy; EQSBW=Equivalent shrunk body weight;
			  !PSWG=Potential shrunk weight gain; ADCP=crude protein (CP) intake; CPD=Crude protein digestibility; CPWG=Crude protein-based weight gain prediction;
			  !DWT=Change in weight gain
                  DE=4.4*TDN(IHD,ISA)
                  EM=.82*DE
                  EMN=EM*(1.37-EM*(.138-EM*.0105))-1.12
                  IF(INKG>.001)THEN
                      SBW=.96*GZWT(IHD,IOW)
                      EGN=EM*(1.42-EM*(.174-EM*.0122))-1.65
                      EMNR=.077*SBW**.75
                      SRW=478.
                      FSBW=.96*GZWM(KKG)
                      RE=MAX(.001,(INKG-EMNR/EMN)*EGN)
                      EQSBW=SBW*SRW/FSBW
                      PSWG=13.91*RE**.9116*EQSBW**(-.6837)
                      ADCP=6250.*ADN(IHD)*INKG  !check units on this.  ADN is g N/g dry matter.  
                      EQEBW=.891*EQSBW
                      CPD=.72
                      CPWG=1.
                      XX2=3.8*EQEBW**.75
                      DO IT=1,10
                          FU=(3.8*SBW**.75+544.7*CPWG-XX2*(.956*CPWG)**1.097)/CPD-ADCP
                          IF(ABS(FU)<.001)EXIT
                          DFDC=(544.7-XX2*1.044*CPWG**.097)/CPD
                          !CPWG=MAX(.1,CPWG-FU/DFDC)
                          X1=FU/DFDC
                          X3=ABS(X1)
                          X2=.25*CPWG
                          IF(X3>X2)X1=X2*X1/X3
                          CPWG=CPWG-X1
                      END DO
                      DWT=MIN(CPWG,PSWG)/.96
                      GZWT(IHD,IOW)=MIN(GZWM(KKG),GZWT(IHD,IOW)+DWT)
                      !WRITE(KW(1),'(5X,A,5I4,5F10.3)')'!!!!!',IY,MO,KDA,IHD,IT,FU,CPWG,DWT,GZWT(IHD,IOW),ADN(IHD)
                  END IF
                  IF(IMPL(KKG)>0)THEN
                      ADTV=1.0 !dry matter intake modifier
                  ELSE    
                      ADTV=0.94
				  END IF
                  EMNI=(GZWT(IHD,IOW)*.96)**.75*(.2435*EMN-.0466*EMN*EMN-.0869) !EMNI is Net Energy for Maintenance Intake (potential value based on forage quality.  Potential is reached if adequate forage is available)
                  DMI=EMNI/MAX(.95,EMN)*ADTV
                  DMDF(IHD)=.001*DMI/RSTX
                  CH4=79.87+9.95*INKG
              CASE(4) ! COW-CALF    
                  ICVW=0
                  IF(ICWD(KKG)>ICVB(KKG))THEN
                      IF(IDA>=ICVB(KKG).AND.IDA<=ICWD(KKG))ICVW=1
                  ELSE
                      IF(IDA>=ICVB(KKG).AND.(IDA<=ND.OR.IDA<=ICWD(KKG)))ICVW=1
                  END IF
                  IF(ICVW==0)THEN
                      ! DRY COWS
                      JJ=GZNB(IHD,ISA)+1
                      DO J=1,JJ
                          DSP(J,ISA)=0.
                      END DO    
                      CAWT=0.
                      FCA=0.
                      DE=4.4*TDN(IHD,ISA)
                      EM=.82*DE
                      EMN=EM*(1.37-EM*(.138-EM*.0105))-1.12
				      EMNI=(GZWT(IHD,IOW)*.96)**.75*(.04997*EMN*EMN+.04631)
                      DMID=EMNI/MAX(.95,EMN)
                      DMI=DMID
                      DMDF(IHD)=.001*DMID/RSTX
                      CH4=79.87+9.95*INKG
                  ELSE
                      ! CALVING WINDOW
                      ICVP=0
                      IF(ICVF(KKG)>ICVB(KKG))THEN
                          IF(IDA>=ICVB(KKG).AND.IDA<=ICVF(KKG))ICVP=1
                      ELSE
                          IF(IDA>=ICVB(KKG).AND.(IDA<=ND.OR.IDA<=ICVF(KKG)))ICVP=1
                      END IF    
                      IF(ICVP>0)THEN
                          FCW=100.*(REAL(IDA)-CV0(KKG))/CVW(KKG)
                          FCO=FCW/(FCW+EXP(SCRP(29,1)-SCRP(29,2)*FCW))
                          DFCO=FCO-FCO0
                          DFCO=DFCO*ATRI(.95,1.,1.05,13)
                          FCA=FCA+DFCO
                          FCO0=FCO
                          GZNL=FCA*GZNB(IHD,ISA)
                      ELSE
                          GZNL=GZNB(IHD,ISA)
                      END IF 
                      NCLV=GZNL
                      DO J=1,NCLV
                          DSP(J,ISA)=DSP(J,ISA)+1.
                      END DO    
                      ! DRY COWS
                      DE=4.4*TDN(IHD,ISA)
                      EM=.82*DE
                      EMN=EM*(1.37-EM*(.138-EM*.0105))-1.12
				      EMNI=(GZWT(IHD,IOW)*.96)**.75*(.04997*EMN*EMN+.04631)
                      DMID=EMNI/MAX(.95,EMN)
                      DMI=DMID
                      DMDFD=.001*DMID/RSTX
                      CH4D=79.87+9.95*INKG
                      ! LACTATING COWS
                      DE=4.4*TDN(IHD,ISA)
                      EM=.82*DE
                      EMN=EM*(1.37-EM*(.138-EM*.0105))-1.12
                      SM1=0.
                      DO J=1,NCLV
                          IF(DSP(J,ISA)>1.)THEN
                              CAWT(J)=CAWT(J)+.0016*GZWT(IHD,IOW)
                          ELSE
                              CAWT(J)=.074*GZWT(IHD,IOW)
                          END IF
                          AD3(IHD)=AD3(IHD)+CAWT(J)
                          X3=1./(.3197*PMLK(KKG))
                          YMLK(IHD)=(DSP(J,ISA)/7.+.3571)/(X3*EXP(.0168*DSP(J,ISA)+.042))
                          AD1(IHD)=AD1(IHD)+YMLK(IHD)
                          EMNI=(GZWT(IHD,IOW)*.96)**.75*(.04997*EMN*EMN+.0384)+.2*YMLK(IHD) 
                          DMIL=EMNI/MAX(.95,EMN)
                          DMI=DMIL
                          CH4L=5.39+20.85*INKG
                          SM1=SM1+DMIL
                      END DO
                      IF(NCLV>0)THEN
                          X4=NCLV
                          AD1(IHD)=AD1(IHD)/X4
                          AD3(IHD)=AD3(IHD)/X4
                          WTGL(IHD,IOW)=AD3(IHD)
                          SM1=SM1/X4
                          DMDFL=.001*SM1/RSTX
                          IF(IDA/=ICWD(KKG).AND.GZNL>0.)THEN
                              ! NURSING BEEF CALVES
                              SM2=0.
                              DO J=1,NCLV
                                  DMIN=(.0033+.000109*DSP(J,ISA))*CAWT(J)
                                  SM2=SM2+DMIN
                              END DO
                              DMIN=SM2/X4
                              DMI=DMIN
                              DMDFN=.001*DMI/RSTX
                              CH4N=79.87+9.95*DMIN
                          END IF    
                      ELSE
                          GZNL=0.
                          FCO0=0.
                          FCA=0.
                          CAWT=0.
                      END IF    
                      RTO=GZNL/GZNB(IHD,ISA)
                      RTO1=1.-RTO
                      DMI=RTO*(DMIL+DMIN)+RTO1*DMID
                      DMDF(IHD)=RTO*(DMDFL+DMDFN)+RTO1*DMDFD
                      CH4=RTO*(CH4L+CH4N)+RTO1*CH4D
                  END IF
              CASE(5)  ! SHEEP
                  DMI=GZIN(KKG)*GZWT(IHD,IOW) 														
                  DMDF(IHD)=.001*DMI/RSTX											
                  YM1=.1083-.0667*TDN(IHD,ISA)
                  CH4=.2875*DMI*YM1    
              CASE(6)  ! GOATS
                  DMI=GZIN(KKG)*GZWT(IHD,IOW) 										
                  DMDF(IHD)=.001*DMI/RSTX											
                  CH4=.000342*GZWT(IHD,IOW)   
			  CASE(7) !Dry (non-lactating) cows that gain/lose weight based on diet quality.  Intake is a constant % of BW
			      DMI=GZIN(KKG)*GZWT(IHD,IOW) 	!calculate intake demand/animal (kg/hd)						
				  TEMPERATURE=TX !NRC does not define temperature as a min, max, mean, etc.  I'm using the APEX daily mean (TX).
				  IF(TEMPERATURE > 25.) THEN !Cody guessed at a minimum critical temp of 22 C to start loss of performance and intake due to heat stress
					!West et al. (2003; J. Dairy Science 86:232-242) found linear decrease of intake.  See Fig. 2. 
					DMI=DMI-(TEMPERATURE-25.)*0.85
				  END IF							
                  DMDF(IHD)=.001*DMI/RSTX		!calculate intake on area basis, converting from kg/hd to Mg/ha
                  DE=4.4*TDN(IHD,ISA)
				  IF(DOM2CP>7.) THEN !Limit Digestible Energy (DE) based on DOM:CP ratio.  Cody came up with this equation to allow intake to remain constant (Lori's request) while limiting growth when DOM:CP is very high (Above 7.0)
					DE=4.4 * (7.*1.05*(ADNX(IHD)*6.25))
					Indicator=1.
				  END IF				  
              
				  EM=.82*DE
                  EMN=EM*(1.37-EM*(.138-EM*.0105))-1.12 !Net Energy for Maintenance of forage concentration Mcal/kg
                  IF(INKG>=0.)THEN
                      SBW=.96*GZWT(IHD,IOW)
                      EGN=EM*(1.42-EM*(.174-EM*.0122))-1.65
					  PREV_PLANE=0.8 + ((BCS(KKG)-1)*0.05) !Adjustment for previous plane of nutrition.  See page 114 of NRC (2001), variable COMP. Nutrient Requirements of Beef Cattle
                      EMNR=.077*PREV_PLANE*SBW**.75 !Net energy maintenance requirement

					  WindSpeed=Min(32.,U10(IRF(ISA))*3.6) !wind speed (km/h). bring in from weather.  Must be <= 32 km/h. Variable WS in NRC.  U10 is in m/s.
					  HeatProd=EM*INKG - EMNR !Heat Production (Mcal/d). Variable HE in NRC
					  SurfArea=0.09*GZWT(IHD,IOW)**0.67 !surface area (m^2). Variable SA in NRC (2016; p.356)
					  !Evap_HL=0.15*HeatProd/SurfArea 
					  HIDE=1. !0.8 for thin hide (Brahman), 1.0 for average, 1.2 for thick (Hereford); Could be put in input file
					  HairCoat=1. !condition of hair coat.  1=no mud; 0.2=heavy mud;  could be set to interact with weather and/or put in input file
					  HairDepth=5. !could potentially be input in grazers file and/or interact with weather
					  ExtIns=MAX(0.,(6.1816 - 0.5575*WindSpeed + 0.0152*WindSpeed**2 + 5.298*HairDepth - 0.4297*HairDepth**2 - 0.1029*WindSpeed*HairDepth)*HairCoat*HIDE) !external insulation (C * M^2 * d/Mcal); Variable EI in NRC (2016; p.356)
					  TissIns=5.25+0.75*BCS(KKG) !Tissue insulation (p. 356 of NRC, 2016). Variable TI in NRC
					  InsTot=TissIns+ExtIns !Total insulation; variable IN in NRC
					  !EATemp=0.0 !Effective Ambient temperature (degrees C) adjusted for thermal radiation.  Set to 0 but should be calculated based on temperature.  see pag 177 of NRC (2016), variable EAT
					  LoCTemp=39.-0.85*InsTot*(HeatProd/SurfArea) !Lower Critical Temp; Variable LCT in NRC (2016; p.356)
					  MEcs=0. !Metabolizable energy required due to cold stress
					  IF(LoCTemp > TEMPERATURE) THEN
					  	  MEcs=SurfArea*(LoCTemp-TEMPERATURE)/InsTot
                      END IF
					  NEhs=0. !Net Energy requirement due to heat stress (Mcal/d)
					  IF(TEMPERATURE > 25.) THEN !Cody guessed at a minimum critical temp of 25 C to start loss of performance due to heat stress.  See also "NUTBAL Technical Support" p. 22
					      !West et al. (2003; J. Dairy Science 8232-242) found linear decrease of intake with temps from 22.5 to 34.4 C (see fig. 2).
						  !NRC p. 356 lists loss of 0.07*NEm for rapid shallow panting and 0.18*NEm at open-mouthed panting. 32 C - 25 C = 9 C.  0.18NEm/9C = 0.02NEm/C
						  NEhs=(TEMPERATURE-25.)*0.02*EMN*INKG
					  END IF
					  NEcs=0.576*MEcs !0.576 is km ion page 177 of NRC(2016)
					  NEmcs=EMNR+NEcs+NEhs !energy for maintenance plus energy for cold and heat stress
					  Slope=0.01 !slope. 0.01=1%. This should be tied to subarea info
					  NEmpa= (0.10 + 0.10 * (MIN(Slope, 0.20)/0.20)) * EMNR !energy for activity

                      EMNE=(INKG*EMN-NEmcs-NEmpa) !Excess intake energy (Mcal)
                      IF(EMNE<0) THEN !check for energy deficit
					    EMNE=EMNE/0.8 !mobilized energy is less efficient
					  END IF
					  DWT = EMNE / 5.82 !weight change based on excess (deficit) energy
					  GZWT(IHD,IOW)=MIN(1.44*GZWM(KKG),MAX(.77*GZWM(KKG),GZWT(IHD,IOW)+DWT)) !adjust weight, within upper and lower limits
					  PMCW = GZWT(IHD,IOW)/GZWM(KKG) !current weight relative to mature weight
					  BCS(KKG) = -20.03 + 36.084*PMCW -11.099*PMCW**2				      !body condition score (BCS)
                  END IF
				  EMNI=(GZWT(IHD,IOW)*.96)**.75*(.04997*EMN*EMN+.04631) !Net Energy for Maintenance intake
				  CH4=79.87+9.95*INKG	!methane production			  
              CASE DEFAULT ! GENERAL
                  DMI=GZIN(KKG)*GZWT(IHD,IOW) 										
                  DMDF(IHD)=.001*DMI/RSTX											
                  CH4=79.87+9.95*INKG    
          END SELECT
          FECE(IHD,ISA)=DMI*(1.-TDN(IHD,ISA))				                
          TDMD=TDMD+DMDF(IHD)													
	  END DO ! HERD DEMAND LOOP 
	  ! DETERMINE SUPPLY OF STL & STD FOR EACH FORAGE
      ! FOR EACH HERD
      TSL=0.
      TSD=0.
      IPFL=0
      IPFD=0
      DO KKP=1,3 ! PREFERENCE LOOP
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              IF(DMDF(IHD)<1.E-10)CYCLE
              DESX(IHD)=0.
              DO I=1,NN ! FORAGE LOOP
                  KKF=JE(I,ISA)
                  IF(KKF==MNC)CYCLE
                  IF(HUI(KKF,ISA)<.2)THEN
                      ANTQF=ANTQN(KKF)
                  ELSE
                      ANTQF=ANTQN(KKF)*(1.-HUI(KKF,ISA))+HUI(KKF,ISA)*ANTQX(KKF)
                  END IF    
                  DESL(KKF)=MIN(.032,CNT(KKF,ISA))*MIN(.9,TDNF(KKF,ISA))*&
                  MIN(1.,1.+ANTQG(KKG)-ANTQF)
                  DESD(KKF)=MIN(.032,BN(3,KKF))*MIN(.9,TDNN(KKF))*&
                  MIN(1.,1.+ANTQG(KKG)-ANTQF)
                  DESX(IHD)=MAX(DESX(IHD),DESL(KKF),DESD(KKF))
              END DO ! FORAGE LOOP   
              DO I=1,NN ! FORAGE LOOP
                  KKF=JE(I,ISA)
                  IF(KKF==MNC)CYCLE
                  RTL(IHD,KKF)=(DESX(IHD)-DESL(KKF))/DESX(IHD)		             
                  RTD(IHD,KKF)=(DESX(IHD)-DESD(KKF))/DESX(IHD)													
                  ! STANDING LIVE        
                  IF(IPFL(IHD,KKF)==0.AND.RTL(IHD,KKF)<PDUX(KKP,KKG))THEN							
                      X2=STL(KKF,ISA)
                      IF(X2>.001)THEN												
                          IF(IDC(KKF)==NDC(8))X2=X2*FTO(KKF)*SLAI(KKF,ISA)/XLAI(KKF,ISA)
                          SPSL(KKP,IHD,KKF)=X2*DMDF(IHD)/TDMD									
                          TSL(KKF)=TSL(KKF)+SPSL(KKP,IHD,KKF)
                          IPFL(IHD,KKF)=1					
                      END IF    
                  END IF
                  ! STANDING DEAD    
                  IF(IPFD(IHD,KKF)==0.AND.RTD(IHD,KKF)<PDUX(KKP,KKG))THEN							
                      X2=STD(KKF,ISA)
                      IF(X2>.001)THEN
                          SPSD(KKP,IHD,KKF)=X2*DMDF(IHD)/TDMD
                          TSD(KKF)=TSD(KKF)+SPSD(KKP,IHD,KKF)
                          IPFD(IHD,KKF)=1
                      END IF
                  END IF  
              END DO ! FORAGE LOOP
          END DO ! HERD LOOP
      END DO ! PREFERENCE LOOP
      SPGZ=0.														
      ! SUPPLY MUST NOT EXCEED STANDING LIVE AND STANDING DEAD
      DO KKP=1,3 ! PREFERENCE LOOP													
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              ! FORAGE LOOP
              DO I=1,NN												
                  KKF=JE(I,ISA)
                  IF(KKF==MNC)CYCLE
                  X2=HE(JT1)*STL(KKF,ISA)							
                  IF(TSL(KKF)>X2.AND.X2>0.)THEN						
                      RTO=X2/TSL(KKF)								
                      SPSL(KKP,IHD,KKF)=SPSL(KKP,IHD,KKF)*RTO		
                  END IF
                  X2=HE(JT1)*STD(KKF,ISA)    						
                  IF(TSD(KKF)>X2)THEN
                      RTO=X2/TSD(KKF)
                      SPSD(KKP,IHD,KKF)=SPSD(KKP,IHD,KKF)*RTO
                  END IF
                  SPGZ(IHD)=SPGZ(IHD)+SPSL(KKP,IHD,KKF)+SPSD(KKP,IHD,KKF)    
              END DO
          END DO
      END DO 
      ! CALCULATE POTENTIAL GRAZING RATE BY PREFERENCE
      GZSL=0.
      GZSD=0.
      YLD=0.
      YLSD=0.
      GZSM=0.
      YTP=0.
      ! CALCULATE FORAGE SUPPLY BY PREFERENCE  
      DO KKP=1,3 ! PREFERENCE LOOP
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              DO I=1,NN ! FORAGE LOOP          
                  KKF=JE(I,ISA)
                  IF(KKF==MNC)CYCLE
                  IF(TSLSD<GZLM(IHD,ISA))CYCLE
                  YTP(KKP,IHD)=YTP(KKP,IHD)+SPSL(KKP,IHD,KKF)+SPSD(KKP,IHD,KKF)
              END DO ! FORAGE LOOP          
          END DO ! HERD LOOP  
      END DO ! PREFERENCE LOOP
      ! CALCULATE SELECTION COEFFICIENTS
      FSEL=0.
      AD2=0.
      DO KKP=1,3 ! PREFERENCE LOOP
          DO KHD=1,NHRD(IOW)
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              IF(YTP(KKP,IHD)<1.E-10)CYCLE
              SELECT CASE(KKP)										
                  CASE(1) ! PREFERRED FORAGE
                      RTO=YTP(1,IHD)/SPGZ(IHD)											
                      PFL=1.-EXP(-3.65*RTO)							
                      FSEL(1,IHD)=PFL
                  CASE(2) ! DESIRABLE FORAGE
                      IF(YTP(3,IHD)>0.)THEN
                          RTO=YTP(3,IHD)/SPGZ(IHD)								
                          UFL=.03197*EXP(2.89*RTO)
                      ELSE
                          UFL=0.
                      END IF                                                              						
                      FSEL(2,IHD)=MAX(0.,1.-UFL-AD2(IHD))                     	
                  CASE(3) ! UNDESIRABLE FORAGE                                          
                      RTO=YTP(3,IHD)/SPGZ(IHD)											
                      FSEL(3,IHD)=.03197*EXP(2.89*RTO)
                      IF(FSEL(3,IHD)>FSEL(2,IHD).AND.YTP(3,IHD)<YTP(2,IHD))THEN
                          X1=1.-AD2(IHD)
                          FSEL(2,IHD)=.8*X1
                          FSEL(3,IHD)=.2*X1
                      END IF                                                           					
              END SELECT
              AD2(IHD)=AD2(IHD)+FSEL(KKP,IHD)
          END DO  
      END DO                  
      ! CALCULATE POTENTIAL GRAZING RATES
      DO KKP=1,3 ! PREFERENCE LOOP														
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              DO I=1,NN ! FORAGE LOOP													
                  KKF=JE(I,ISA)
                  IF(KKF==MNC)CYCLE
                  IF(TSLSD<GZLM(IHD,ISA))CYCLE                                                				
                  IF(YTP(KKP,IHD)<1.E-10)CYCLE
                  RTX=SPSL(KKP,IHD,KKF)/YTP(KKP,IHD)           
                  GZSL(KKP,IHD,KKF)=MIN(SPSL(KKP,IHD,KKF),FSEL(KKP,IHD)*RTX*DMDF(IHD)) 
                  RTX=SPSD(KKP,IHD,KKF)/YTP(KKP,IHD)           
                  GZSD(KKP,IHD,KKF)=MIN(SPSD(KKP,IHD,KKF),FSEL(KKP,IHD)*RTX*DMDF(IHD))
              END DO    
          END DO
      END DO
      ! SUM GRAZING
      DO KKP=1,3 ! PREFERENCE LOOP                                 
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              ! FORAGE LOOP
              DO I=1,NN
                  KKF=JE(I,ISA)
                  IF(KKF==MNC)CYCLE
                  GZSM(IHD)=GZSM(IHD)+GZSL(KKP,IHD,KKF)+GZSD(KKP,IHD,KKF) !add actual grazed standing live and dead to total grazed for this herd (preliminary calculation.  Will be recalculated below a couple times.)
              END DO
          END DO
      END DO    
      ! ADJUST GRAZING RATES TO TRY TO MATCH DEMAND
      DO KKP=1,3 ! PREFERENCE LOOP
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              IF(GZSM(IHD)<1.E-10)CYCLE
              RTO=DMDF(IHD)/GZSM(IHD) !ratio of demand to actual grazing intake
              IF(RTO<1.)THEN !if there is a deficit of intake, consume some more forage if available
                  DO I=1,NN ! FORAGE LOOP
                      KKF=JE(I,ISA)
                      IF(KKF==MNC)CYCLE
                      IF(TSLSD<GZLM(IHD,ISA))CYCLE
                      GZSL(KKP,IHD,KKF)=RTO*GZSL(KKP,IHD,KKF)
                      GZSD(KKP,IHD,KKF)=RTO*GZSD(KKP,IHD,KKF)
                  END DO    
              END IF    
          END DO
      END DO    
      ! RECALCULATE TOTAL GRAZING AND REDUCE SUPPLY
      GZSM=0.
      SPGZ=0.
      DO KKP=1,3 ! PREFERENCE LOOP
          DO KHD=1,NHRD(IOW)
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              DO I=1,NN ! FORAGE LOOP
                  KKF=JE(I,ISA)
                  IF(KKF==MNC)CYCLE
                  GZSM(IHD)=GZSM(IHD)+GZSL(KKP,IHD,KKF)+GZSD(KKP,IHD,KKF) !add actual grazed standing live and dead to total grazed for this herd (preliminary calculation.  Will be recalculated below once more.)
                  SPSL(KKP,IHD,KKF)=SPSL(KKP,IHD,KKF)-GZSL(KKP,IHD,KKF) !update supply of standing live forage
                  SPSD(KKP,IHD,KKF)=SPSD(KKP,IHD,KKF)-GZSD(KKP,IHD,KKF)
                  SPGZ(IHD)=SPGZ(IHD)+SPSL(KKP,IHD,KKF)+SPSD(KKP,IHD,KKF)
              END DO
          END DO
      END DO
      ! CALCULATE DEFICIT
      DO KHD=1,NHRD(IOW) ! HERD LOOP
          IHD=KOW(KHD,IOW)
          IF(IGZR(IHD,ISA)==0)CYCLE
          KKG=NGZ(IHD,ISA)
          IF(KKG==0)CYCLE
          DEF(IHD)=DMDF(IHD)-GZSM(IHD) !the difference between demand and grazed forage is the deficit
      END DO    
      ! IF DEFICIT EXISTS TRY TO MEET DEMAND
      AD2=0.
      DO KKP=1,3 ! PREFERENCE LOOP
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              IF(DEF(IHD)>0.)THEN
                  DO I=1,NN ! FORAGE LOOP
                      KKF=JE(I,ISA)
                      IF(KKF==MNC)CYCLE
                      IF(TSLSD<GZLM(IHD,ISA))CYCLE
                      IF(ABS(AD2(IHD)-DEF(IHD))<1.E-10.OR.SPGZ(IHD)<1.E-10)CYCLE
                      RTO=SPSL(KKP,IHD,KKF)/SPGZ(IHD)
                      X1=MIN(DEF(IHD)*RTO,SPSL(KKP,IHD,KKF))
                      AD2(IHD)=AD2(IHD)+X1
                      GZSL(KKP,IHD,KKF)=GZSL(KKP,IHD,KKF)+X1
                      SPSL(KKP,IHD,KKF)=SPSL(KKP,IHD,KKF)-X1
                      RTO=SPSD(KKP,IHD,KKF)/SPGZ(IHD)
                      X1=MIN(DEF(IHD)*RTO,SPSD(KKP,IHD,KKF))
                      AD2(IHD)=AD2(IHD)+X1
                      GZSD(KKP,IHD,KKF)=GZSD(KKP,IHD,KKF)+X1
                      SPSD(KKP,IHD,KKF)=SPSD(KKP,IHD,KKF)-X1
                  END DO
              END IF
          END DO
      END DO 
      ! CALCULATE YIELD FOR STANDING LIVE AND STANDING DEAD
      YTP=0.
      GZSM=0.
      ATDN=0.
      ADNX=0.
      ADP=0.
      ADH2O=0.
      DO KKP=1,3 ! PREFERENCE LOOP
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              DO I=1,NN ! FORAGE LOOP
                  KKF=JE(I,ISA)
                  IF(KKF==MNC)CYCLE
                  YLD(KKF)=YLD(KKF)+GZSL(KKP,IHD,KKF)
                  YLSD(KKF)=YLSD(KKF)+GZSD(KKP,IHD,KKF)
                  YGZL(IHD,KKF)=YGZL(IHD,KKF)+GZSL(KKP,IHD,KKF)
                  YGZD(IHD,KKF)=YGZD(IHD,KKF)+GZSD(KKP,IHD,KKF)
                  YTP(KKP,IHD)=YTP(KKP,IHD)+GZSL(KKP,IHD,KKF)+GZSD(KKP,IHD,KKF)
                  GZSM(IHD)=GZSM(IHD)+GZSL(KKP,IHD,KKF)+GZSD(KKP,IHD,KKF) !add actual grazed standing live and dead to total grazed for this herd (final calculation except for adding supplement if needed.)
                  ATDN(IHD)=ATDN(IHD)+TDNF(KKF,ISA)*GZSL(KKP,IHD,KKF)+TDNN(KKF)*GZSD(KKP,IHD,KKF) 
                  ADNX(IHD)=ADNX(IHD)+CNT(KKF,ISA)*GZSL(KKP,IHD,KKF)+BN(3,KKF)*GZSD(KKP,IHD,KKF) !sum of N consumed here.  kg/ha?
                  ADP(IHD)=ADP(IHD)+CPT(KKF)*GZSL(KKP,IHD,KKF)+BP(3,KKF)*GZSD(KKP,IHD,KKF)  !sum of P consumed here.  kg/ha?
                  XX=100.*HUI(KKF,ISA)
                  WC=.5*(1.-XX/(XX+EXP(SCRP(28,1)-SCRP(28,2)*XX)))+.25
                  ADH2O(IHD)=ADH2O(IHD)+GZSL(KKP,IHD,KKF)*WC/(1.-WC)+.111*GZSD(KKP,IHD,KKF) !sum of water consumed here. kg/ha?
              END DO !FORAGE LOOP
          END DO ! HERD LOOP
      END DO ! PREFERENCE LOOP
      ! UPDATE STANDING LIVE & STANDING DEAD
      DO I=1,NN
          KKF=JE(I,ISA)
          IF(KKF==MNC)CYCLE
          STL(KKF,ISA)=MAX(0.,STL(KKF,ISA)-YLD(KKF))
          STD(KKF,ISA)=STD(KKF,ISA)-YLSD(KKF)
      END DO                   
      ! ESTIMATE SUPPLEMENTAL FEEDING = DEFICIT
      IF(IHAY>0)THEN 
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              DEF(IHD)=DEF(IHD)-AD2(IHD)
              SPFD(IHD)=MAX(0.,DEF(IHD))
              YTP(2,IHD)=YTP(2,IHD)+SPFD(IHD)
              GZSM(IHD)=GZSM(IHD)+SPFD(IHD) !add supplemental feed (e.g. hay) to the grazed quantity
              ATDN(IHD)=ATDN(IHD)+.5*SPFD(IHD)
              ADNX(IHD)=ADNX(IHD)+.015*SPFD(IHD)
              ADP(IHD)=ADP(IHD)+.0021*SPFD(IHD)
              ADH2O(IHD)=ADH2O(IHD)+.1*SPFD(IHD)
              FVI(IHD)=(YTP(1,IHD)+.6*YTP(2,IHD))/GZSM(IHD)
              SMMF(IHD,MO,ISA)=SMMF(IHD,MO,ISA)+SPFD(IHD)
              IF(SPFD(IHD)>.0001)SMMFD(IHD,MO,ISA)=SMMFD(IHD,MO,ISA)+1.
          END DO
      END IF
      ! OUTPUT GRAZING RESULTS TO .OUT & .DGZ    
      DO KHD=1,NHRD(IOW) ! HERD LOOP
          IHD=KOW(KHD,IOW)
          IF(IGZR(IHD,ISA)==0)CYCLE
          KKG=NGZ(IHD,ISA)
          IF(KKG==0)CYCLE
          IF(GZSM(IHD)>0.)THEN
              ATDN(IHD)=ATDN(IHD)/GZSM(IHD)
              ADNX(IHD)=ADNX(IHD)/GZSM(IHD)
              ADP(IHD)=ADP(IHD)/GZSM(IHD)
              ADN(IHD)=ADNX(IHD)
          END IF    
          RSTX=WSA(ISA)/GZNB(IHD,ISA)
          ! CALCULATE DRINKING WATER INTAKE l/(hd*d).
          XX=MAX(TX,2.)
          H2OD=.076*GZWT(IHD,IOW)**.7086*XX**.7146
          SMMH2OD(IHD,MO,ISA)=SMMH2OD(IHD,MO,ISA)+H2OD
          ! WATER CONSUMPTION = DRINKING AND GRAZING FORAGE WATER CONTENT.
          ADH2O(IHD)=1000.*ADH2O(IHD)*RSTX
          X10=H2OD+ADH2O(IHD)
          ! URINE VOLUME
          VURN(IHD,ISA)=.31*X10
          SMMH2O(IHD,MO,ISA)=SMMH2O(IHD,MO,ISA)+X10
          X5=1.E-4*VURN(IHD,ISA)/RSTX
          RFV(IRF(ISA))=RFV(IRF(ISA))+X5
          SMMURN(IHD,MO,ISA)=SMMURN(IHD,MO,ISA)+VURN(IHD,ISA)
          ! N & P IN URINE 
          FNU=11.25*ADNX(IHD)+.288
          XMLKN=.0054*AD1(IHD)
          Intake_kg_ha=1000.*GZSM(IHD)*RSTX
          XXN=Intake_kg_ha*ADNX(IHD)-XMLKN
          URNN=FNU*XXN
          XMLKP=.001*YMLK(IHD)
          XXP=Intake_kg_ha*ADP(IHD)-XMLKP
          URNP=.04*XXP
          ! N & P IN MANURE
          FECN=XXN-URNN
          FECP=XXP-URNP
          ! ADD N & P ANIMAL DEPOSITION TO SOIL
          ! URINE APPLIED TO SOIL
          APMU=VURN(IHD,ISA)/RSTX
          XCZ=.02*URNN
          FNO(MFT)=XCZ/APMU
          FN(MFT)=(URNN-XCZ)/APMU
          FNMA(MFT)=.99
          XCZ=.99*URNP
          FP(MFT)=XCZ/APMU
          FPO(MFT)=(URNP-XCZ)/APMU
          FOC(MFT)=.43*FN(MFT)
          FSLT(MFT)=.001
          SMMMU(IHD,MO,ISA)=SMMMU(IHD,MO,ISA)+FECE(IHD,ISA)
          CALL NFERT(APMU,8,IAMF(ISA),IHD,KT2,JRT)
          ! FECES APPLIED TO SOIL
          APMU=FECE(IHD,ISA)/RSTX
          XCZ=.96*FECN
          FNO(MFT)=XCZ/APMU
          FN(MFT)=(FECN-XCZ)/APMU
          FNMA(MFT)=.99
          XCZ=.48*FECP
          FP(MFT)=XCZ/APMU
          FPO(MFT)=(FECP-XCZ)/APMU
          FOC(MFT)=.35
          CALL NFERT(APMU,8,IAMF(ISA),IHD,KT2,JRT)
          SMMTDN(IHD,MO,ISA)=SMMTDN(IHD,MO,ISA)+ATDN(IHD)
          SMMGD(IHD,MO,ISA)=SMMGD(IHD,MO,ISA)+1.
          !SUM=0.
		  TDN(IHD,ISA)=ATDN(IHD)
          ! CALCULATE HERD TDN AND OUTPUT TO .DGZ 
          DO I=1,NN ! FORAGE LOOP
              KKF=JE(I,ISA)
              IF(KKF==MNC)CYCLE
              !Cody replaced the following lines because the TDN of an individual forage can exceed TDN of the average diet.  Thus, TDN would not be the actual overall diet TDN.
			  !Cody's new line is TDN(IHD,ISA)=ATDN(IHD) which is placed above the DO loop
			  !TDN(IHD,ISA)=MAX(TDNN(KKF),ATDN(IHD))
              !SUM=SUM+TDN(IHD,ISA)
              Y1=1000.*YGZL(IHD,KKF)
              Y2=1000.*YGZD(IHD,KKF)
              SMMSL(IHD,KKF,MO,ISA)=SMMSL(IHD,KKF,MO,ISA)+Y1
              SMMSD(IHD,KKF,MO,ISA)=SMMSD(IHD,KKF,MO,ISA)+Y2
              !IF(Y1<1.E-5.AND.Y2<1.E-5)CYCLE
              X1=1000.*DMDF(IHD)
              X2=1000.*SPSL(1,IHD,KKF)
              X3=1000.*SPSL(2,IHD,KKF)
              X4=1000.*SPSL(3,IHD,KKF)
              X5=1000.*SPSD(1,IHD,KKF)
              X6=1000.*SPSD(2,IHD,KKF)
              X7=1000.*SPSD(3,IHD,KKF)
              X8=1000.*SPFD(IHD)
              K1=GZNB(IHD,ISA)
              K2=GZNL
              IF(KFL(18)>0)THEN
                  I1=GZNB(IHD,ISA)+.5
                  IF(NDGZ==1)THEN
                      KDG=KDG0+ISA
                      WRITE(KW(KDG),3)ISA,NBSA(ISA),IOW,IYR,MO,KDA,&
                      WSA(ISA),GNAM(KKG),GZNB(IHD,ISA),RSTX,GZWT(IHD,IOW),TX,H2OD,&
                      CPNM(KKF),STL(KKF,ISA),HUI(KKF,ISA),CNT(KKF,ISA),STD(KKF,ISA),&
                      BN(3,KKF),X1,FVI(IHD),TDNF(KKF,ISA),X2,X3,X4,X5,X6,X7,X8,&
                      Y1,Y2,ATDN(IHD),ADNX(IHD),DOM2CP,Indicator,FECE(IHD,ISA),FECN,FECP,&
                      VURN(IHD,ISA),URNN,URNP,ADH2O(IHD),X10,AD1(IHD),K1,K2,&
                      AD3(IHD),DMIN,CH4,DE,EM,EMN,EMNE,EMNR,NEmcs,NEmpa,LoCTemp,MECS,HeatProd,&
					  SurfArea,NEhs,TEMPERATURE,WindSpeed,InsTot,ExtIns,TissIns,DWT,CPWG,PSWG,BCS(KKG),GZSM(IHD),INKG,PREV_PLANE
                  ELSE 
                      WRITE(KW(18),3)ISA,NBSA(ISA),IOW,IYR,MO,KDA,&
                      WSA(ISA),GNAM(KKG),GZNB(IHD,ISA),RSTX,GZWT(IHD,IOW),TX,H2OD,&
                      CPNM(KKF),STL(KKF,ISA),HUI(KKF,ISA),CNT(KKF,ISA),STD(KKF,ISA),&
                      BN(3,KKF),X1,FVI(IHD),TDNF(KKF,ISA),X2,X3,X4,X5,X6,X7,X8,&
                      Y1,Y2,ATDN(IHD),ADNX(IHD),DOM2CP,Indicator,FECE(IHD,ISA),FECN,FECP,&
                      VURN(IHD,ISA),URNN,URNP,ADH2O(IHD),X10,AD1(IHD),K1,K2,&
                      AD3(IHD),DMIN,CH4,DE,EM,EMN,EMNE,EMNR,NEmcs,NEmpa,LoCTemp,MEcs,HeatProd,&
					  SurfArea,NEhs,TEMPERATURE,WindSpeed,InsTot,ExtIns,TissIns,DWT,CPWG,PSWG,BCS(KKG),GZSM(IHD),INKG,PREV_PLANE
                  END IF
              END IF    
          END DO    
          HDTM(IHD,ISA)=HDTM(IHD,ISA)+1.
          !IF(SUM>0.)TDN(IHD,ISA)=SUM/REAL(NN)
      END DO
      ! UPDATE PLANT STATUS.
      TOT2=0.
      TOT3=0.
      DO I=1,NN ! FORAGE LOOP
          KKF=JE(I,ISA)
          IF(KKF==MNC)CYCLE
          HUF(KKF,ISA)=MAX(HUF(KKF,ISA),HU(KKF,ISA))
          DMF(KKF,ISA)=DM1(KKF,ISA)
          TRA(KKF,ISA)=SRA(KKF,ISA)+TRA(KKF,ISA)
          IF(RD(KKF,ISA)>RDF(KKF,ISA))RDF(KKF,ISA)=RD(KKF,ISA)
          AJHI(KKF,ISA)=0.
          YHE=YLD(KKF)/HE(JT1)
	      X1=MIN(YHE/(STL(KKF,ISA)+1.E-5),.9)
          ZZ=MAX(.01,1.-X1)
          ZZ2=ZZ**5
          ZZ3=MAX(.01,1.-.5*X1)
          ZZ4=MAX(.01,1.-ZZ*HUI(KKF,ISA)**2)
          IF(HUI(KKF,ISA)<.3)THEN
              ZZ4=0.
              ZZ5=1.
          ELSE
              IF(HUI(KKF,ISA)<.5)THEN
                  ZZ4=HU(KKF,ISA)*MAX(.01,.2*X1)
				  ZZ4=MIN(ZZ4,TGX(KKF,ISA)*.5)
                  ZZ5=1.-1.5*(HUI(KKF,ISA)-.3)
              ELSE
                  ZZ4=0.
                  ZZ5=.7
              END IF
          END IF
          IF(DETIN(KKF)>0.)THEN
              CZZ4=MAX(.01,MIN(.3,X1))
              ZZ4=HU(KKF,ISA)*CZZ4*.5
          END IF    
          ZZ5=MAX(.01,1.-ZZ5*YLD(KKF)/1.5)
		  CPHT(KKF,ISA)=MAX(.001,CPHT(KKF,ISA)*ZZ2)
          HUMIN=MIN(HU(KKF,ISA),.4*PHU(KKF,IHU(KKF,ISA),ISA))
		  HU(KKF,ISA)=MAX(HUMIN,HU(KKF,ISA)-ZZ4)
          SLAI(KKF,ISA)=MAX(.01,SLAI(KKF,ISA)*(1.-X1*0.5))
          Y4=1000.*YLSD(KKF)*BN(3,KKF)
          Y5=1000.*YLSD(KKF)*BP(3,KKF)
          XY1=STDN(KKF,ISA)-Y4
          IF(XY1>0.)THEN
              STDN(KKF,ISA)=XY1
          ELSE
              Y4=STDN(KKF,ISA)
              STDN(KKF,ISA)=0.
          END IF
          XY1=STDP(KKF,ISA)-Y5
          IF(XY1>0.)THEN
              STDP(KKF,ISA)=XY1
          ELSE
              Y5=STDP(KKF,ISA)
              STDP(KKF,ISA)=0.
          END IF
          X1=DM(KKF,ISA)+1.E-5
          CNLV=UN1(KKF,ISA)/X1
          CPLV=UP1(KKF,ISA)/X1
          X5=MIN(YLD(KKF)*CPLV,UP1(KKF,ISA))
          YLN=MIN(.9*(UN1(KKF,ISA)+STDN(KKF,ISA)),YLD(KKF)*CNLV)
          YLP=MIN(.9*(UP1(KKF,ISA)+STDP(KKF,ISA)),YLD(KKF)*CPLV)
          X1=1./HE(JT1)-1.
          XZ1=YLD(KKF)*X1
          XZ2=MIN(YLSD(KKF)*X1,STD(KKF,ISA))
          X11=XZ1+XZ2
          STD(KKF,ISA)=STD(KKF,ISA)-XZ2
          XZ1=YLN*X1
          XZ2=Y4*X1
          X10=XZ1+XZ2
          STDN(KKF,ISA)=STDN(KKF,ISA)-XZ2
          UN1(KKF,ISA)=UN1(KKF,ISA)-XZ1
          TOT3=TOT3+X10
	      JJK=KKF
	      ! ADD TRAMPLED STL AND NUTRIENT CONTENT TO LAYER 1 RESIDUE  
          CALL NCNSTD(X11,X10,LD1)
          FOP(LD1,ISA)=MAX(.01,FOP(LD1,ISA)+YLP*X1+Y5*X1)
          YLD2(KKF,ISA)=YLD2(KKF,ISA)+YLD(KKF)+YLSD(KKF)
          JD(ISA)=KKF
          SRA(KKF,ISA)=0.
          UN1(KKF,ISA)=UN1(KKF,ISA)-YLN
          UP1(KKF,ISA)=UP1(KKF,ISA)-X5
          DM(KKF,ISA)=DM(KKF,ISA)-YHE
	      IF(DM(KKF,ISA)<RW(KKF,ISA))RW(KKF,ISA)=DM(KKF,ISA)
          STL(KKF,ISA)=MAX(0.,DM(KKF,ISA)-RW(KKF,ISA))
          YLN=YLN+Y4
          TOT3=TOT3+YLN
          YLP=YLP+Y5
          YLNF(KKF,ISA)=YLNF(KKF,ISA)+YLN
          YLPF(KKF,ISA)=YLPF(KKF,ISA)+YLP
          TYN(ISA)=TYN(ISA)+YLN
          TYP(ISA)=TYP(ISA)+YLP
          !IF(NOP>0.OR.NBSA(ISA)==ISAP)THEN
              !IF(YLD(KKF)>0..OR.YLSD(KKF)>0.)THEN
                  !Y1=1000.*YLD(KKF)
                  !Y2=1000.*YLSD(KKF)
                  !KGZ=1
                  !WRITE(KW(1),29)ISA,NBSA(ISA),IYR,MO,KDA,IOW,TIL(JT1),&
                  !CPNM(KKF),STL(KKF,ISA),STD(KKF,ISA),Y1,Y2,HUI(KKF,ISA),YLN
              !END IF 
          !END IF    
          TOT2=TOT2++UN1(KKF,ISA)+STDN(KKF,ISA)
      END DO
      DF=TOT1-TOT3-TOT2
      IF(ABS(DF)>.001)WRITE(KW(1),1)IY,MO,KDA,TOT1,TOT2,TOT3,DF
      SMM(141,MO,ISA)=SMM(141,MO,ISA)+1.                         
      RETURN
    1 FORMAT(1X,'!!!!!',3I4,4E12.5)       
    2 FORMAT(1X,2I8,1X,3I4,2X,'GNAM= ',A16,2X,'RSTK= ',F8.3,' ha/hd',2X,&
      'CPNM= ',A4,2X,'STD= ',F7.3,' t/ha',2X,'STL= ',F7.3,' t/ha',2X,&
      'DMD= ',F7.4,' t/ha',2X,'SUPLSL(P-D-U)= ',3F7.3,' t/ha',2X,&
      'SUPLSD(P-D-U)= ',3F7.3,' t/ha',2X,'HAY= ',F7.1,' kg/ha',2X,'YLD= ',&
      F7.1,' kg/ha',2X,'YLSD= ',F7.1,' kg/ha')
    3 FORMAT(1X,2I8,1X,I4,1X,3I4,1X,F10.2,A16,F10.2,F10.3,F10.0,F10.2,&
      F10.1,1X,A4,30F10.3,2I10,F10.0,F10.3,F10.0,F10.2,F10.2,&
	  22F10.3,F10.1,F10.1)      
   29 FORMAT(1X,2I8,1X,3I4,2X,'IDON=',I4,2X,A8,2X,A4,2X,'STL=',F7.3,&
      't/ha',2X,'STD=',F7.3,'t/ha',2X,'GZSL=',F7.1,' kg/ha',2X,'GZSD=',&
      F7.1,'kg/ha',2X,'HUSC = ',F6.2,2X,'YLN=',F7.4,' kg/ha')
      END