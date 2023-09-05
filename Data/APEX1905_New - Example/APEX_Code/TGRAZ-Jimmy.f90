      SUBROUTINE TGRAZ
!     APEX1905
!     THIS SUBPROGRAM SIMULATES ANIMAL GRAZING USING CONCEPTS FROM PHY-GROW
!     GRAZING PREFERENCE COMPONENT.
!=====================================================================================================
!       IHD     Herd ID
!       IHRD    = 0 FOR PHY-GROW GRAZING PREFERENCE MODEL
!               = 1 FOR LEVEL 1(HYBRID) GRAZING MODE(HERD FILE REQUIRED) 
!               = 2 FOR LEVEL 2(AUTOMATIC) GRAZING MODE(HERD FILE REQUIRED)  
!       FECE    feces weight
!       WC      water content of standing live vegetation, 10% water was assumed for standing dead
!       H2OD    water drank 
!       CNT,BN  N FRACTION IN PLANT WHEN GROWTH IS 0.,.5,1.0
!       IDON    Owner identification number (cols. 1-4)
!       EMN     net energy for maintenance,  EMN=EM*(1.37-EM*(.138-EM*.0105))-1.12
!       EM      metabolizable energy, EM=.82*DE
!       DE      digestable energy,    DE=4.4*TDN(IHD,ISA)
!       TDN     TOTAL DIGESTABLE NUTRIENTS
!       DMI     daily dry matter intake
!       HUI     Heat Unit Index
!       NHD     NUMBER OF GRAZERS PER HERD BY OWNER(HD)
!       TDNX    MAX TOTAL DIGESTABLE NUTRIENTS (TDN) %
!       TDNN    MIN TDN %
!       IDA     DAY OF MONTH SIMULATION BEGINS
!       STL     standing plant live biomass (kg/ha)
!       STD     standing dead residue (kg/ha)
!       UN1     N in plant biomass (kg N/ha)
!       STDN    N in standing dead (kg N/ha)
!       RSTK    stocking density (animals/ha)
!       WSA     subarea (ha)
!       ANTQX  MAXIMUM ANTIQUALITY TOLERANCE FOR GRAZERS
!       ANTQN  MINIMUM ANTIQUALITY TOLERANCE FOR GRAZERS
!       GZNB    # of grazing animals in each subarea
!       IRO     rotation no.
!       KT2     TILLAGE SEQUENCE # IN .OPC.  START GRAZE & STOP GRAZE
!       ISA     subarea index
!       GZWX    animal weigth before entering a grazing area (for rotation among pastures) (kg)
!       GZWI    initial animal weight (kg)
!       RSTX    area per animal (ha/animal)
!       GZIN    intake as % of body weight (?)
!       FVI     forage value index
!       FTO     FRACTION TURNOUT (COTTON LINT/STRIPPER YLD)
!       ICWD    WEANING DAY OF YEAR
!       ICVB    BEGINNING DAY OF YEAR OF CALVING
!       ICVF    FINAL DAY OF YEAR OF CALVING
!       IGZD    = 1 FOR GRAZING INTAKE RATE CONSTANT
!               = 2 FOR VARIABLE GRAZING INTAKE RATE BASED ON FORAGE QUALITY
!               = 3 FOR VARIABLE GRAZING INTAKE RATE BASED ON NET ENERGY EQ FROM NRC BEEF YEARLINGS
!               = 4 FOR VARIABLE GRAZING INTAKE RATE BASED ON NET ENERGY EQ FROM NRC DRY COWS
!               = 5 FOR VARIABLE GRAZING INTAKE RATE BASED ON NET ENERGY EQ FROM NRC LACTATING COWS
!               = 6  Goat?
!               = 7  If I want to add new weight algorithm. 
!       DWT     daily weight gain    
!       TDN     total digestable nutrient
!       GZWM    naturity weight (kg/head)
!       CAWT    weight of nursing calf (kg/head)
!       GZLM    grazing limit (Mg/ha)
!       GZSL    STANDING LIVE FORAGE THAT WAS GRAZED (KG/HA). DOES NOT INCLUDE TRAMPLED FORAGE.
!       GZSD    STANDING DEAD FORAGE THAT WAS GRAZED (KG/HA). DOES NOT INCLUDE TRAMPLED FORAGE.      
!       NGZ     1=non feeding area; 2=feed lot  (KKG)
!       PDUX    THERESHOLD RATIOS LIGNIN/N FOR GRAZING PREFERENCES
!               (PREFERED, DESIRABLE, UNDESIRABLE)
!       VURN    Urine volume produced
!       GZSM    daily dry matter requirement per animal T/ha/animal (calculated from GZIN, GZWT and GZNB)
!       IGZR    >0 IF IHD IS GRAZING ON ISA THAT DAY
!       IGZX    ISA FOR IHD,IOW (USED IN AUTO GRAZING ROTATION)
!       INWT    SETS INITIAL WEIGHT ON 1ST DAY OF GRAZING TO GZWX>0 OR GZWI
!       MNC     MAX # CROPS USED IN SIMULATION FROM APEXDIM
!       NCR     # TIMES CROP JJK IS HARVESTED IN ISA
!       JP      0 WHEN CROP IS PLANTED; 1 AT HARVEST
!       ADN     N consumed, N concentration of forage diet for the herd
!       ADCP    Crude protein consumed
!       YTP     (here is two dimensional, but in MAIN_1605.f90, it is one dimensional as grazing limit (t/ha))
!       TDNF    TDN OF FORAGE KKF ON ISA
!       DESL    DESIRABILITY OF STL FORAGE KKF 
!       DESD    DESIRABILITY OF STD FORAGE KKF
!       DESX    MAX DESIRABILITY (DESL,DESD)
!       SLAI    LAI OF FORAGE
!       XLAI    MAX LAI CONSIDERING PLANT POP
!       TSL     TOTAL SUPPLY OF Live Biomass
!       TSD     TOTAL SUPPLY OF Dead Biomass
!       FSEL    GRAZING SELECTION OF PREFERRED, DESIRABLE, UNDESIRABLE (PDU) BY FORAGE BY HERD 
!       SPGZ    TOTAL SUPPLY OF STL + STD FORAGE BY HERD
!       SPSL    SUPPLY OF STL BY PREFERENCE BY HERD BY FORAGE
!       SPSD    SUPPLY OF STD BY PREFERENCE BY HERD BY FORAGE
!       RTO     RATIO USED TEMPORARILY IN MANY PLACES
      USE PARM
      DIMENSION IPFD(KOMX(IOW),MNC),IPFL(KOMX(IOW),MNC)
      DIMENSION AD1(KOMX(IOW)),AD2(KOMX(IOW)),AD3(KOMX(IOW)),ADH2O(KOMX(IOW)),&
      ADNX(KOMX(IOW)),ADP(KOMX(IOW)),DESD(MNC,MSA),DESL(MNC,MSA),DESX(KOMX(IOW)),&
      SPFD(KOMX(IOW)),SPGZ(KOMX(IOW)),DEF(MHD,MOW)
      !DIMENSION DSP(1000,MOW),DMDF(KOMX(IOW)),GZSX(KOMX(IOW)),DEF(KOMX(IOW))
      DIMENSION RTD(KOMX(IOW),MNC),RTL(KOMX(IOW),MNC),TSD(MNC,MSA),TSL(MNC,MSA),YLSD(MNC,MSA)
      DIMENSION FSEL(3,KOMX(IOW)),YGZL(KOMX(IOW),MNC),YGZD(KOMX(IOW),MNC),YTP(3,KOMX(IOW))
      DIMENSION GZSD(3,KOMX(IOW),MNC,MSA),GZSL(3,KOMX(IOW),MNC,MSA),SPSL(3,KOMX(IOW),MNC,MSA), &
                SPSD(3,KOMX(IOW),MNC,MSA),GZSLT(3,KOMX(IOW),MNC),GZSDT(3,KOMX(IOW),MNC), &
                GZSS(MNC,MSA),SPSLT(3,KOMX(IOW),MNC),SPSDT(3,KOMX(IOW),MNC)
      DIMENSION TSLSD(MSA),TTOT1(MSA),TTOT2(MSA),TTOT3(MSA),SUMTDN(MSA)
      DATA FCA,FCO0/2*0./
      REAL ME_calf, SUMFOR(MHD,MSA)
      LD1=LID(1,ISA) 
      !SUM1=WPML(LD1,ISA)+WPMS(LD1,ISA)+WPMA(LD1,ISA)+WPMU(LD1,ISA)+&
      !WPOU(LD1,ISA)+WPO(LD1,ISA)+FOP(LD1,ISA)+UP1(1,ISA)
      !WRITE(KW(1),'(10E16.6)')WPML(LD1,ISA),WPMS(LD1,ISA),WPMA(LD1,ISA),&
      !WPMU(LD1,ISA),WPOU(LD1,ISA),WPO(LD1,ISA),FOP(LD1,ISA),UP1(1,ISA)
      WSAX=WSA(ISA)
      KT2=KT(ISA)
      KRT=0
	  KGZ=0
      N1=KT(ISA)
      DMIN=0.
      DMDFN=0.
      GZNL=0.
      SPFD=0.
	  YGZL=0.
      YGZD=0.
      GZSL=0.
      GZSD=0.
      GZSLT=0.
      GZSDT=0.
      SPSLT=0.
      SPSDT=0.
      SPSL=0.
      SPSD=0.
      SPGZ=0.
      GZSS=1.0
      CH4=0.
      IOW=IDON(ISA)
      KDG0=2*MSA+MSO       	 
      NN=NCP(IRO(ISA),ISA)	!number of forages
      TSLSD(ISA)=0.
      TTOT1(ISA)=0.
      ME_calf=0.
! ADDED BY LIWANG MA 7-10-2-22, FROM CODY PAPER HTTP://dx.doi.org/10.1016/j.ecolmodel.2017.02.004
      IF (RSTK(IRO(ISA),KT2,ISA)>=2.0) THEN
            HE(JT1)=0.4
      ELSE IF (RSTK(IRO(ISA),KT2,ISA)<0.8) THEN
            HE(JT1)=0.7
      ELSE
            HE(JT1)=0.9-0.25*RSTK(IRO(ISA),KT2,ISA)
      ENDIF
      DO K=1,NN
          KKF=JE(K,ISA)
          IF(KKF==MNC)CYCLE
          TSLSD(ISA)=TSLSD(ISA)+STL(KKF,ISA)+STD(KKF,ISA)
          TTOT1(ISA)=TTOT1(ISA)+UN1(KKF,ISA)+STDN(KKF,ISA)
      END DO   
      ! DETERMINE STOCKING RATE, INITIAL HERD WEIGHT, POTENTIAL DAILY GRAZING RATE 
      DO KHD=1,NHRD(IOW) ! HERD LOOP
          IHD=KOW(KHD,IOW)
          IF(IHRD<2)THEN
              IF((RSTK(IRO(ISA),KT2,ISA)<1.E-10.AND.GZNB(IHD,ISA)<1.E-10).OR.IHD/=JGRZ(IRO(ISA),KT2,ISA))CYCLE
              IGZR(IHD,ISA)=1
              IHD0=IHD
          ELSE    
              IF(IGZX(IHD,IOW)/=ISA.OR.IHD/=JGRZ(IRO(ISA),KT2,ISA))CYCLE
              IF(IBSL(IHD,IOW)>0)IGZR(IHD,ISA)=0
              IHD0=IHD
          END IF    
          IF(RSTK(IRO(ISA),KT2,ISA)>0.)THEN
              GZNB(IHD,ISA)=WSAX/RSTK(IRO(ISA),KT2,ISA)
          ELSE    
              GZNB(IHD,ISA)=NHD(IHD,IOW)
              RSTK(IRO(ISA),KT2,ISA)=WSAX/GZNB(IHD,ISA)
          END IF    
          KKG=NGZ(IHD,ISA)
          IF(KKG==0)CYCLE
          IF(GZWX(IRO(ISA),KT2,ISA)>0.)THEN 
              X1=GZWX(IRO(ISA),KT2,ISA)
          ELSE
              X1=GZWI(KKG)
          END IF    
          GZWT(IHD,IOW)=MAX(GZWT(IHD,IOW),X1)
          !IF(KFL(37)>0.AND.LGZ(ISA)==0)WRITE(KW(37),101)ISA,NBSA(ISA),&
          !IYR,MO,KDA,TIL(JT1),IHD,GZNB(IHD,ISA),RSTK(IRO(ISA),KT2,ISA),GZWT(IHD,IOW)
          LGZ(ISA)=1
          IF(INWT(IHD,IOW)==0)THEN
              WTBG(IHD,IOW)=GZWT(IHD,IOW)
              IGZB(IHD,IOW)=IDA
          END IF    
          INWT(IHD,IOW)=1

          
          DO K=1,NN	! FORAGE LOOP			
              KKF=JE(K,ISA)
              IF(KKF==MNC)CYCLE
              IF(JP(KKF,ISA)<=0)THEN
                  JP(KKF,ISA)=1														
                  NCR(KKF,ISA)=NCR(KKF,ISA)+1											
              END IF
          END DO !  FORAGE LOOP    
      END DO ! HERD LOOP
      IHD=IHD0
      DO KHD=1,NHRD(IOW) ! HERD DEMAND LOOP
          IHD=KOW(KHD,IOW)
          IF(IGZR(IHD,ISA)==0)CYCLE   !cycling the 11 herd needed
          KKG=NGZ(IHD,ISA)
          IF(KKG==0)CYCLE
          !IF(IHRD>0)CYCLE   !no rotation grazing
		  IF (IHOUR(IY,KT2,ISA,IHD).LT. 12.0 .AND. IHOUR(IY,KT2,ISA,IHD) .GT.0.) THEN
		    GZDH(IHD,ISA)=IHOUR(IY,KT2,ISA,IHD)/12.0
            TBOS(ISA)= TSLSD(ISA)
          ELSE
			GZDH(IHD,ISA)=1.0
          END IF
       END DO
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
          IF(GZWT(IHD,IOW)<1.E-10)GZWT(IHD,IOW)=GZWI(NGZ(IHD,ISA))
          RSTX=WSAX/GZNB(IHD,ISA)
          !PRINT 1002,IY,MO,KDA,IHD,ISA,KDAN(IHD,IOW),GZDH(IHD,ISA),TDNH(IHD,IOW),GZSMH(IHD,IOW),GZSMD(IHD,IOW),RSTX,TTBS(IHD),'test0'
          1002 format (6I4,6(f10.4),A8)
          IF(KDAN(IHD,IOW).NE. KDA)THEN
              IF(GZSMH(IHD,IOW)>0.)THEN
                XX1=GZSMH(IHD,IOW)
                TDN(IHD,ISA)= TDNH(IHD,IOW)/GZSMH(IHD,IOW)    
                GZSMH(IHD,IOW)=0.
                TDNH(IHD,IOW)=0.
                !KDAN(IHD,IOW)=KDA
              ELSEIF(GZSMD(IHD,IOW)>0.)THEN
                 XX1=GZSMD(IHD,IOW) 
                 TDN(IHD,ISA)=TDND(IHD,IOW)
                 GZSMD(IHD,IOW)=0.
                 !KDAN(IHD,IOW)=KDA
                 TDND(IHD,IOW)=0. 
              ELSE
                 GZSM(IHD,IOW)=0. 
              ENDIF   
              Total_B(IHD)=TTBS(IHD)
              KDAN(IHD,IOW)=KDA
              DAYSUB(IHD,IOW)=1  
              TTBS(IHD)=0.
          ELSE
            GZSM(IHD,IOW)=0. 
			DAYSUB(IHD,IOW)=0
            XX1=0.
          ENDIF
          CALL SWN1530 
          IF(KFL(37)>0. .AND. DAYSUB(IHD,IOW)==1)WRITE(KW(37),101)ISA,NBSA(ISA),&   !OUTPUT DALIY GRAZING MANAGEMENT  AND.LGZ(ISA)==1
          IYR,MO,KDA,TIL(JT1),IOW,IHD,GZNB(IHD,ISA),RSTK(IRO(ISA),KT2,ISA),GZWT(IHD,IOW),&
          IHOUR(IY,KT2,ISA,IHD),GZDH(IHD,ISA),WTGL(IHD,IOW),Hays(IHD,IOW),Hay_total(IHD,IOW),&
          ICVW,ICVP,SW15,SW30,SNN15+SNA15,SNN30+SNA30
          !XX1=1000.*GZSM(IHD,IOW)*RSTX          !change this to caculate the XX1 in each pastures at the end
          
		IF(DAYSUB(IHD,IOW)==1)THEN   !only daily update  																   
          SELECT CASE(IGZD(KKG))
              CASE(1) !GENERAL
                  DMI=GZIN(KKG)*GZWT(IHD,IOW) 														
                  DMDF(IHD,IOW)=.001*DMI											
                  CH4=79.87+9.95*XX1
              CASE(2) ! BEEF YEARLINGS
                  IF(IDA==ICVB(KKG))THEN
                     YEARK=IY
                     Hays(IHD,IOW)=Hay_amount(IHD,IY,IOW)
                     !SHays(IHD)=SHay_amount(IHD,IY)
                     Hay_Q(IHD,IOW)=MAX(Hay_quality(IHD,IY,IOW),0.5)  !quality of hay percent of TDN 
                     Hay_pri(IHD,IOW)=Hay_P(IHD,IY,IOW)      
                   ENDIF				   
                  IF(IDA==ICWD(KKG))THEN
                    Hay_total(IHD,IOW)=0.
                  END IF                      
				  
                  DE=4.4*TDN(IHD,ISA)
                  EM=.82*DE
                  EMN=EM*(1.37-EM*(.138-EM*.0105))-1.12
                  IF(XX1>.001)THEN
                      SBW=.96*GZWT(IHD,IOW)
                      EGN=EM*(1.42-EM*(.174-EM*.0122))-1.65
                      EMNR=.077*SBW**.75
                      SRW=478.
                      FSBW=.96*GZWM(KKG)
                      RE=MAX(.001,(XX1-EMNR/EMN)*EGN)
                      EQSBW=SBW*SRW/FSBW
                      PSWG=13.91*RE**.9116*EQSBW**(-.6837)
                      ADCP=6250.*ADN(IHD)*XX1
                      EQEBW=.891*EQSBW
                      CPD=.72
                      CPWG=PSWG
                      XX2=3.8*EQEBW**.75
                      IF(CPWG>.0001)THEN
                          DO IT=1,20
                              FU=(3.8*SBW**.75+544.7*CPWG-XX2*(.956*CPWG)**1.097)/CPD-ADCP
                              IF(ABS(FU)<.001)EXIT
                              DFDC=(544.7-XX2*1.044*CPWG**.097)/CPD
													
                              X1=FU/DFDC
                              X3=ABS(X1)
                              X2=.25*CPWG
                              IF(X3>X2)X1=X2*X1/X3
                              CPWG=CPWG-X1
                          END DO
                      END IF    
                      DWT=MIN(CPWG,PSWG)/.96
                      WTGL(IHD,IOW)=DWT   !FANG ADDED
                      GZWT(IHD,IOW)=MIN(GZWM(KKG),GZWT(IHD,IOW)+DWT)
                      IF(IT>20)WRITE(KW(1),'(5X,A,5I4,5F10.3)')'!!!!!',IY,MO,KDA,IHD,IT,FU,CPWG,DWT,GZWT(IHD,IOW),ADN(IHD)
                  END IF
                  IF(IMPL(KKG)>0)THEN
                      ADTV=1.0
                  ELSE    
                      ADTV=0.94
                  END IF
                  X1=.2435*EMN-.0466*EMN*EMN-.0869
                  EMNI=(GZWT(IHD,IOW)*.96)**.75*X1*ADTV
                  DMI=MAX(EMNI/MAX(.95,EMN),GZIN(KKG)*GZWT(IHD,IOW)) 														
                  DMDF(IHD,IOW)=.001*DMI   !make DMDF no related to area (demand per animal)
                  CH4=79.87+9.95*XX1
                  CH4E(IHD,IOW)=CH4
              CASE(3) ! COW-CALF  
                  IF(IDA==ICWD(KKG))THEN
        	        DSP=0
                    Hay_total(IHD,IOW)=0.
                  END IF   
				  
                  ICVW=0
                  !IF(ICWD(KKG)>ICVB(KKG))THEN
                  !    IF(IDA>=ICVB(KKG).AND.IDA<=ICWD(KKG))ICVW=1
                  !ELSE
                  !    IF(IDA<=ICWD(KKG).AND.NCLV(IHD)>0)ICVW=1
                  !END IF

                  IF(ICWD(KKG)>ICVB(KKG))THEN
                      IF(IDA>=ICVB(KKG).AND.IDA<=ICWD(KKG))ICVW=1
                  ELSE
                      IF(IDA>=ICVB(KKG))THEN
                          IF(IDA<=ND)ICVW=1
                      ELSE
                          IF(IDA<=ICWD(KKG).AND.NCLV(IHD)>0)ICVW=1
                          
                      END IF    
                  END IF                  
                  
                  IF(ICVW==0)THEN
                     YEARK=IY
                     Hays(IHD,IOW)=Hay_amount(IHD,IY,IOW)
                     !SHays(IHD)=SHay_amount(IHD,IY)
                     Hay_Q(IHD,IOW)=MAX(Hay_quality(IHD,IY,IOW),0.5)  !quality of hay percent of TDN 
                     Hay_pri(IHD,IOW)=Hay_P(IHD,IY,IOW)     
					 
                      ! DRY COWS
                      JJ=GZNB(IHD,ISA)+1
                      DO J=1,JJ
                          DSP(J,IHD,IOW)=0.
                      END DO    
                      CAWT=0.
                      FCA=0.
                      DE=4.4*TDN(IHD,ISA)
                      EM=.82*DE
                      EMN=EM*(1.37-EM*(.138-EM*.0105))-1.12
					  !using the steer gain model to estimate pregnant cow gain
                  IF(XX1>.001)THEN
                      SBW=.96*GZWT(IHD,IOW)                     !shrunk body weight
                      EGN=EM*(1.42-EM*(.174-EM*.0122))-1.65     !net energy for gain content of the jth feedstuff, Mcal/kg
                      EMNR=.077*SBW**.75                        !feed for maintainance. NEm in the book?
                      SRW=478.                                  !standard reference weight
                      FSBW=.96*GZWM(KKG)                        !final shrunk body weight
                      RE=MAX(.001,(XX1-EMNR/EMN)*EGN)           !retained energy, Mcal/day, different than in the reference manual 
                      EQSBW=SBW*SRW/FSBW                        !equivalent shrunk body weight, kg
                      PSWG=13.91*RE**.9116*EQSBW**(-.6837)      !predicted daily gain in shrunk body weight
                      ADCP=6250.*ADN(IHD)*XX1
                      EQEBW=.891*EQSBW                          !equivilant empty body weight
                      CPD=.72                                   !crude protein content in diet?
                      XX2=3.8*EQEBW**.75                        !metabolizable protein requirement for maintenance, g/day, in the manual, it was SBW, not EQEBW
                      CPWG=PSWG
                      IF(CPWG>.0001)THEN
                          DO IT=1,20
                              FU=(3.8*SBW**.75+544.7*CPWG-XX2*(.956*CPWG)**1.097)/CPD-ADCP
                              IF(ABS(FU)<.001)EXIT
                              DFDC=(544.7-XX2*1.044*CPWG**.097)/CPD
                              X1=FU/DFDC
                              X3=ABS(X1)
                              X2=.25*CPWG
                              IF(X3>X2)X1=X2*X1/X3
                              CPWG=CPWG-X1
                              IF(ABS(X1)<.001)EXIT
                          END DO
                      END IF    					  
                      DWT=MIN(CPWG,PSWG)/.96                      !converted back to animal weight. 
                      GZWT(IHD,IOW)=MIN(GZWM(KKG),GZWT(IHD,IOW)+DWT)                      
                   ENDIF
				      EMNI=(GZWT(IHD,IOW)*.96)**.75*(.04997*EMN*EMN+.04631)
                      !DMID=EMNI/MAX(.95,EMN)
					  DMID=MAX(EMNI/MAX(.95,EMN),GZIN(KKG)*GZWT(IHD,IOW))													 
                      DMI=DMID
                      DMDF(IHD,IOW)=.001*DMID           !demand/animal
                      CH4=79.87+9.95*XX1
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
                          FCO=FCW/(FCW+EXP(SCRP(31,1)-SCRP(31,2)*FCW))
                          DFCO=FCO-FCO0
                          DFCO=DFCO*ATRI(.95,1.,1.05,13)
                          FCA=FCA+DFCO
                          FCO0=FCO
                          GZNL=FCA*GZNB(IHD,ISA)
                      ELSE
                          GZNL=GZNB(IHD,ISA)+.5
                      END IF 
                      NCLV(IHD)=GZNL
                      DO J=1,NCLV(IHD)
						  DSP(J,IHD,IOW)=DSP(J,IHD,IOW)+1.						  
                          !DSP(J,IHD,IOW)=DSP(J,IHD,IOW)+1.
                      END DO    
                      ! DRY COWS
                      DE=4.4*TDN(IHD,ISA)
                      EM=.82*DE
                      EMN=EM*(1.37-EM*(.138-EM*.0105))-1.12
				      EMNI=(GZWT(IHD,IOW)*.96)**.75*(.04997*EMN*EMN+.04631)
                      DMID=MAX(EMNI/MAX(.95,EMN),GZIN(KKG)*GZWT(IHD,IOW))
                      DMI=DMID
                      DMDFD=.001*DMID
                      CH4D=79.87+9.95*XX1
                      ! LACTATING COWS
                      DE=4.4*TDN(IHD,ISA)
                      EM=.82*DE
                      EMN=EM*(1.37-EM*(.138-EM*.0105))-1.12

                      SBW=.96*GZWT(IHD,IOW)                     !shrunk body weight
                      EGN=EM*(1.42-EM*(.174-EM*.0122))-1.65     !net energy for gain content of the jth feedstuff, Mcal/kg
                      EMNR=.077*SBW**.75                        !feed for maintainance. NEm in the book?
                      SRW=478.                                  !standard reference weight
                      FSBW=.96*GZWM(KKG)                        !final shrunk body weight
                      RE=MAX(.001,(XX1-EMNR/EMN))               !retained energy, kg/day                      
                      SM1=0.
                      DO J=1,NCLV(IHD)
                          !IF(DSP(J,IHD,IOW)>1.)THEN
                              !CAWT(IHD,J)=CAWT(IHD,J)+.0016*GZWT(IHD,IOW)
                          !ELSE
                              !CAWT(IHD,J)=.074*GZWT(IHD,IOW)
                          !END IF
                      IF(MO==1 .AND. KDA==1)THEN   !adapting for SR change from year to year
                        IF(CAWT(IHD,J).LE. .085*GZWT(IHD,IOW))THEN
                            CAWT(IHD,J)=WTGL(IHD,IOW)*1.0
                            DSP(J,IHD,IOW)=DSP(INT(NCLV(IHD)/2.0),IHD,IOW)
                        ENDIF
                      ENDIF    
                      IF(DSP(J,IHD,IOW)>1.)THEN
                        IF(CAWT(IHD,J)==0.)CAWT(IHD,J)= WTGL(IHD,IOW)   !adapting for the SR change from year to your or from SA to SA
                        IF(XX1>.001)THEN
                          X3=1./(.3197*PMLK(KKG))
                          YMLK(IHD)=(DSP(J,IHD,IOW)/7.+.3571)/(X3*EXP(.0168*DSP(J,IHD,IOW)+.042))
                          ME_calf=0.077*(0.96*CAWT(IHD,J))**0.75                    ! from NRBC  Mkal/kg
                          REEG=MAX(.0,(RE-ME_calf/EMN)*EGN)
                          ADCP=6250.*ADN(IHD)*RE 
                          IF(DSP(J,IHD,IOW)<60.)THEN
                              REEG=(YMLK(IHD)*0.751-ME_calf)                        !assuming same EMN in MILK and maintainance
                          ENDIF
						  !RE=MIN(RE, MAX(YMLK(IHD)*0.749,RE))                       !RE not change? using RE instead of YMLKsuppose all the produced milk feed to calf  1kg milk =0.763 Mkcal/kg
						  !ME_calf=0.291*(CAWT(IHD,J)**0.75)*0.23866                !Based on Jouven et al,m 2008
                          SBW=.96*CAWT(IHD,J)                                       !shrunk body weight
                          SRW=478.                                                  !standard reference weight from cow
                          FSBW=.96*GZWI(KKG)                                        !final shrunk body weight from calf weight
                          EQSBW=SBW*SRW/FSBW                                        !equivalent shrunk body weight, kg
						  EQEBW=.891*EQSBW    
                          IF (REEG>0.) THEN
                             PSWG=13.91*(REEG)**.9116*EQSBW**(-.6837)        !predicted calf daily gain in shrunk body weight
                          ELSE
                             PSWG=0.
                          END IF						  
                          CPD=.72                                   !crude protein content in diet?
                          CPWG=PSWG
                          XX2=3.8*EQEBW**.75  !metabolizable protein requirement for maintenance, g/day, in the manual, it was SBW, not EQEBW
                          IF((CPWG>.0001).AND.(DSP(J,IHD,IOW)>60))THEN
                              ADCP=6250.*ADN(IHD)*(REEG/EGN+ME_calf/EMN)                    !potentiall ADCP for calf: not used for now
                              DO IT=1,20
                                 FU=(3.8*SBW**.75+544.7*CPWG-XX2*(.956*CPWG)**1.097)/CPD-ADCP
                                 IF(ABS(FU)<.001)EXIT
                                 DFDC=(544.7-XX2*1.044*CPWG**.097)/CPD
                                 X1=FU/DFDC
                                 X3=ABS(X1)
                                 X2=.25*CPWG
                                 IF(X3>X2)X1=X2*X1/X3
                                 CPWG=CPWG-X1
                                 IF(ABS(X1)<.001)EXIT
                                 !PRINT*,CPWG,X1,ME_calf
                              END DO
                          END IF    						  
                          CAWT(IHD,J)=MIN(CAWT(IHD,J)+.0021*GZWT(IHD,IOW),CAWT(IHD,J)+PSWG/0.96)  !add this ratio to reduce weight gain
                          END IF
                        ELSE
                          CAWT(IHD,J)=.074*GZWT(IHD,IOW)
                        END IF
                          AD3(IHD)=AD3(IHD)+CAWT(IHD,J)
                          X3=1./(.3197*PMLK(KKG))
                          YMLK(IHD)=(DSP(J,IHD,IOW)/7.+.3571)/(X3*EXP(.0168*DSP(J,IHD,IOW)+.042))
                          AD1(IHD)=AD1(IHD)+YMLK(IHD)
                          EMNI=((GZWT(IHD,IOW)*.96)**.75*(.04997*EMN*EMN+.0384)+.2*YMLK(IHD))/MAX(.95,EMN) 
                          EMNI=((GZWT(IHD,IOW)*.96)**.75*(.04000*EMN*EMN+.0300)+.2*YMLK(IHD))/MAX(.95,EMN) 
                          EMNI_calf=(CAWT(IHD,J)*.96)**.75*(.2435*EMN-0.0466*EMN*EMN-0.1128)
                          EMNR=.077*(GZWT(IHD,IOW)*.96)**.75
                          DMIL=EMNI + EMNI_calf/MAX(.95,EMN)                                           ! need to further test
                          DMIL=EMNR/MAX(.95,EMN) +.2*YMLK(IHD)/MAX(.95,EMN) + EMNI_calf/MAX(.95,EMN)   ! need to further test
                          !DMIL=EMNI/MAX(.95,EMN)
                          DMI=DMIL
                          CH4L=5.39+20.85*XX1
                          SM1=SM1+DMIL
                      END DO
                      IF(NCLV(IHD)>0)THEN
                          X4=NCLV(IHD)
                          AD1(IHD)=AD1(IHD)/X4
                          AD3(IHD)=AD3(IHD)/X4
                          WTGL(IHD,IOW)=AD3(IHD)
                          MILK(IHD,IOW)=AD1(IHD)
                          SM1=SM1/X4
                          DMDFL=.001*SM1
                          IF(IDA/=ICWD(KKG).AND.GZNL>0.)THEN
                              ! NURSING BEEF CALVES
                              SM2=0.
                              DO J=1,NCLV(IHD)
                                  DMIN=(.0033+.000109*DSP(J,IHD,IOW))*CAWT(IHD,J)
                                  !DMIN=(.0033+.000109*DSP(J,IHD,IOW))*CAWT(IHD,J)
                                  SM2=SM2+DMIN
                              END DO
                              DMIN=SM2/X4
                              DMI=DMIN
                              DMDFN=.001*DMI
                              CH4N=79.87+9.95*DMIN
                          END IF    
                      ELSE
                          GZNL=0.
                          FCO0=0.
                          FCA=0.
                          CAWT=0.
                      END IF    
                      RTO=GZNL/GZNB(IHD,ISA)
                      IF(RTO>=1.)RTO=1.0
                      RTO1=1.-RTO
                      DMI=RTO*(DMIL+DMIN)+RTO1*DMID
                      DMDF(IHD,IOW)=(RTO*(DMDFL+DMDFN)+RTO1*DMDFD)      ! DMDF as demand/animal
                      CH4=RTO*(CH4L+CH4N)+RTO1*CH4D
                      CH4E(IHD,IOW)=CH4
                  END IF
              CASE(4)  ! SHEEP
                  DMI=GZIN(KKG)*GZWT(IHD,IOW) 														
                  DMDF(IHD,IOW)=.001*DMI											
                  YM1=.1083-.0667*TDN(IHD,ISA)
                  CH4=.2875*DMI*YM1    
              CASE(5)  ! GOATS
                  DMI=GZIN(KKG)*GZWT(IHD,IOW) 										
                  DMDF(IHD,IOW)=.001*DMI											
                  CH4=.000342*GZWT(IHD,IOW)    
              CASE DEFAULT ! GENERAL
                  DMI=GZIN(KKG)*GZWT(IHD,IOW) 										
                  DMDF(IHD,IOW)=.001*DMI											
                  CH4=79.87+9.95*XX1    
          END SELECT
          !FECE(IHD,ISA)=DMI*(1.-TDN(IHD,ISA))				                
        ENDIF
		  DMDT(IHD,IOW)=DMDF(IHD,IOW)/RSTX
          TDMD=TDMD+DMDT(IHD,IOW)													
          !PRINT 1002,IY,MO,KDA,IHD,ISA,KDAN(IHD,IOW),GZDH(IHD,ISA),TSLSD,GZSX(IHD,IOW),DMDT(IHD,IOW)*RSTX,TDN(IHD,ISA),RSTX,'test1'      
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
              IF(DMDT(IHD,IOW)<1.E-10)CYCLE
              DESX(IHD)=0.
              DO I=1,NN ! FORAGE LOOP
                  KKF=JE(I,ISA)
                  IF(KKF==MNC)CYCLE
                  IF(HUI(KKF,ISA)<.2)THEN
                      ANTQF(ISA)=ANTQN(KKF)
                  ELSE
                      ANTQF(ISA)=ANTQN(KKF)*(1.-HUI(KKF,ISA))+HUI(KKF,ISA)*ANTQX(KKF)
                  END IF    
                  DESL(KKF,ISA)=MIN(.032,CNT(KKF,ISA))*MIN(.9,TDNF(KKF,ISA))*&
                  MIN(1.,1.+ANTQG(KKG)-ANTQF(ISA))
                  DESD(KKF,ISA)=MIN(.032,BN(3,KKF))*MIN(.9,TDNN(KKF))*&
                  MIN(1.,1.+ANTQG(KKG)-ANTQF(ISA))
                  DESX(IHD)=MAX(DESX(IHD),DESL(KKF,ISA),DESD(KKF,ISA),.01)
              END DO ! FORAGE LOOP   
              DO I=1,NN ! FORAGE LOOP
                  KKF=JE(I,ISA)
                  IF(KKF==MNC)CYCLE
                  RTL(IHD,KKF)=(DESX(IHD)-DESL(KKF,ISA))/DESX(IHD)		             
                  RTD(IHD,KKF)=(DESX(IHD)-DESD(KKF,ISA))/DESX(IHD)													
                  ! STANDING LIVE        
                  IF(IPFL(IHD,KKF)==0.AND.RTL(IHD,KKF)<PDUX(KKP,KKG))THEN							
                      X2=STL(KKF,ISA)
                      IF(X2>.001)THEN												
                          IF(IDC(KKF)==NDC(8))X2=X2*FTO(KKF)*SLAI(KKF,ISA)/XLAI(KKF,ISA)
                          SPSL(KKP,IHD,KKF,ISA)=X2*DMDT(IHD,IOW)/TDMD									
                          TSL(KKF,ISA)=TSL(KKF,ISA)+SPSL(KKP,IHD,KKF,ISA)
                          IPFL(IHD,KKF)=1					
                      END IF    
                  END IF
                  ! STANDING DEAD    
                  IF(IPFD(IHD,KKF)==0.AND.RTD(IHD,KKF)<PDUX(KKP,KKG))THEN							
                      X2=STD(KKF,ISA)
                      IF(X2>.001)THEN
                          SPSD(KKP,IHD,KKF,ISA)=X2*DMDT(IHD,IOW)/TDMD
                          TSD(KKF,ISA)=TSD(KKF,ISA)+SPSD(KKP,IHD,KKF,ISA)
                          IPFD(IHD,KKF)=1
                      END IF
                  END IF  
              END DO ! FORAGE LOOP
          END DO ! HERD LOOP
      END DO ! PREFERENCE LOOP
!      SPGZ=0.														
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
                  IF(TSL(KKF,ISA)>X2.AND.X2>0.)THEN						
                      RTO=X2/TSL(KKF,ISA)								
                      SPSL(KKP,IHD,KKF,ISA)=SPSL(KKP,IHD,KKF,ISA)*RTO		
                  END IF
                  X2=HE(JT1)*STD(KKF,ISA)    						
                  IF(TSD(KKF,ISA)>X2)THEN
                      RTO=X2/TSD(KKF,ISA)
                      SPSD(KKP,IHD,KKF,ISA)=SPSD(KKP,IHD,KKF,ISA)*RTO
                  END IF
!                  SPGZ(IHD)=SPGZ(IHD)+SPSL(KKP,IHD,KKF,ISA)+SPSD(KKP,IHD,KKF,ISA) 
              END DO
          END DO
      END DO 
      ! CALCULATE POTENTIAL GRAZING RATE BY PREFERENCE
      GZSL=0.
      GZSD=0.
      YLD=0.
      YLSD=0.
      GZSX=0.
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
                  IF(TSLSD(ISA)<GZLM(IHD,ISA)/HE(JT1))CYCLE
                  YTP(KKP,IHD)=YTP(KKP,IHD)+SPSL(KKP,IHD,KKF,ISA)+SPSD(KKP,IHD,KKF,ISA)
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
                      IF (SPGZ(IHD)>0) THEN
                      RTO=YTP(1,IHD)/SPGZ(IHD)	
                      ELSE
                      RTO=1.
                      ENDIF
                      PFL=1.-EXP(-3.65*RTO)							
                      FSEL(1,IHD)=PFL
                  CASE(2) ! DESIRABLE FORAGE  is this correct? why not YTP(2,IHD)? Liwang Ma
                      IF(YTP(2,IHD)>0.)THEN
                      IF (SPGZ(IHD)>0) THEN
                          RTO=YTP(2,IHD)/SPGZ(IHD)								
                      ELSE
                      RTO=1.
                      ENDIF
                          UFL=.03197*EXP(2.89*RTO)
                      ELSE
                          UFL=0.
                      END IF                                                              						
                      FSEL(2,IHD)=MAX(0.,1.-UFL-AD2(IHD))                     	
                  CASE(3) ! UNDESIRABLE FORAGE                                          
                      IF (SPGZ(IHD)>0) THEN
                      RTO=YTP(3,IHD)/SPGZ(IHD)											
                      ELSE
                      RTO=1.
                      ENDIF
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
                  IF(TSLSD(ISA)<GZLM(IHD,ISA)/HE(JT1))CYCLE                                                				
                  IF(YTP(KKP,IHD)<1.E-10)CYCLE
                  RTX=SPSL(KKP,IHD,KKF,ISA)/YTP(KKP,IHD)           
                  !GZSL(KKP,IHD,KKF)=MIN(SPSL(KKP,IHD,KKF),FSEL(KKP,IHD)*RTX*DMDT(IHD,IOW)) 
				  GZSL(KKP,IHD,KKF,ISA)=MIN(SPSL(KKP,IHD,KKF,ISA),FSEL(KKP,IHD)*RTX*DMDT(IHD,IOW)*GZDH(IHD,ISA))																						 
                  RTX=SPSD(KKP,IHD,KKF,ISA)/YTP(KKP,IHD)           
                  !GZSD(KKP,IHD,KKF)=MIN(SPSD(KKP,IHD,KKF),FSEL(KKP,IHD)*RTX*DMDT(IHD,IOW))
				  GZSD(KKP,IHD,KKF,ISA)=MIN(SPSD(KKP,IHD,KKF,ISA),FSEL(KKP,IHD)*RTX*DMDT(IHD,IOW)*GZDH(IHD,ISA))																						
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
                  GZSX(IHD,IOW)=GZSX(IHD,IOW)+GZSL(KKP,IHD,KKF,ISA)+GZSD(KKP,IHD,KKF,ISA)
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
              IF(GZSX(IHD,IOW)<1.E-10)CYCLE
              !RTO=DMDT(IHD,IOW)/GZSX(IHD,IOW)
			  RTO=DMDT(IHD,IOW)*GZDH(IHD,ISA)/GZSX(IHD,IOW)											   
              IF(RTO<1.)THEN
                  DO I=1,NN ! FORAGE LOOP
                      KKF=JE(I,ISA)
                      IF(KKF==MNC)CYCLE
                      IF(TSLSD(ISA)<GZLM(IHD,ISA)/HE(JT1))CYCLE
                      GZSL(KKP,IHD,KKF,ISA)=RTO*GZSL(KKP,IHD,KKF,ISA)
                      GZSD(KKP,IHD,KKF,ISA)=RTO*GZSD(KKP,IHD,KKF,ISA)
                  END DO    
              END IF    
          END DO
      END DO    
      ! RECALCULATE TOTAL GRAZING AND REDUCE SUPPLY
      GZSX=0.
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
                  GZSX(IHD,IOW)=GZSX(IHD,IOW)+GZSL(KKP,IHD,KKF,ISA)+GZSD(KKP,IHD,KKF,ISA)
                  SPSL(KKP,IHD,KKF,ISA)=SPSL(KKP,IHD,KKF,ISA)-GZSL(KKP,IHD,KKF,ISA)
                  SPSD(KKP,IHD,KKF,ISA)=SPSD(KKP,IHD,KKF,ISA)-GZSD(KKP,IHD,KKF,ISA)
                  SPGZ(IHD)=SPGZ(IHD)+SPSL(KKP,IHD,KKF,ISA)+SPSD(KKP,IHD,KKF,ISA)
              END DO
          END DO
      END DO
      ! CALCULATE DEFICIT
      DO KHD=1,NHRD(IOW) ! HERD LOOP
          IHD=KOW(KHD,IOW)
          IF(IGZR(IHD,ISA)==0)CYCLE
          KKG=NGZ(IHD,ISA)
          IF(KKG==0)CYCLE
          !DEF(IHD,IOW)=DMDT(IHD,IOW)-GZSX(IHD,IOW)
		  DEF(IHD,IOW)=DMDT(IHD,IOW)*GZDH(IHD,ISA)-GZSX(IHD,IOW)
          !PRINT*,MO,KDA,IHD,DEF(IHD),GZSX(IHD,IOW)
      END DO    
      ! IF DEFICIT EXISTS TRY TO MEET DEMAND
      AD2=0.
      DO KKP=1,3 ! PREFERENCE LOOP
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              IF(DEF(IHD,IOW)>0.)THEN
                  DO I=1,NN ! FORAGE LOOP
                      KKF=JE(I,ISA)
                      IF(KKF==MNC)CYCLE
                      IF(TSLSD(ISA)<GZLM(IHD,ISA)/HE(JT1))CYCLE
                      IF(ABS(AD2(IHD)-DEF(IHD,IOW))<1.E-10.OR.SPGZ(IHD)<1.E-10)CYCLE
                      RTO=SPSL(KKP,IHD,KKF,ISA)/SPGZ(IHD)
                      X1=MIN(DEF(IHD,IOW)*RTO,SPSL(KKP,IHD,KKF,ISA))
                      AD2(IHD)=AD2(IHD)+X1
                      GZSL(KKP,IHD,KKF,ISA)=GZSL(KKP,IHD,KKF,ISA)+X1
                      SPSL(KKP,IHD,KKF,ISA)=SPSL(KKP,IHD,KKF,ISA)-X1
                      RTO=SPSD(KKP,IHD,KKF,ISA)/SPGZ(IHD)
                      X1=MIN(DEF(IHD,IOW)*RTO,SPSD(KKP,IHD,KKF,ISA))
                      AD2(IHD)=AD2(IHD)+X1
                      GZSD(KKP,IHD,KKF,ISA)=GZSD(KKP,IHD,KKF,ISA)+X1
                      SPSD(KKP,IHD,KKF,ISA)=SPSD(KKP,IHD,KKF,ISA)-X1
                  END DO
              END IF
          END DO
      END DO 
      ! CALCULATE YIELD FOR STANDING LIVE AND STANDING DEAD
      YTP=0.
      ATDN=0.
      ADNX=0.
      ADP=0.
      ADH2O=0.
      GZSX=0.
      DO KKP=1,3 ! PREFERENCE LOOP
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
              DO I=1,NN ! FORAGE LOOP
                  KKF=JE(I,ISA)
                  IF(KKF==MNC)CYCLE
                  YLD(KKF,ISA)=YLD(KKF,ISA)+GZSL(KKP,IHD,KKF,ISA)
                  YLSD(KKF,ISA)=YLSD(KKF,ISA)+GZSD(KKP,IHD,KKF,ISA)
                  YGZL(IHD,KKF)=YGZL(IHD,KKF)+GZSL(KKP,IHD,KKF,ISA)
                  YGZD(IHD,KKF)=YGZD(IHD,KKF)+GZSD(KKP,IHD,KKF,ISA)
                  YTP(KKP,IHD)=YTP(KKP,IHD)+GZSL(KKP,IHD,KKF,ISA)+GZSD(KKP,IHD,KKF,ISA)
                  GZSX(IHD,IOW)=GZSX(IHD,IOW)+GZSL(KKP,IHD,KKF,ISA)+GZSD(KKP,IHD,KKF,ISA)
                  ATDN(IHD)=ATDN(IHD)+TDNF(KKF,ISA)*GZSL(KKP,IHD,KKF,ISA)+TDNN(KKF)*GZSD(KKP,IHD,KKF,ISA)
                  ADNX(IHD)=ADNX(IHD)+CNT(KKF,ISA)*GZSL(KKP,IHD,KKF,ISA)+BN(3,KKF)*GZSD(KKP,IHD,KKF,ISA)
                  ADP(IHD)=ADP(IHD)+CPT(KKF)*GZSL(KKP,IHD,KKF,ISA)+BP(3,KKF)*GZSD(KKP,IHD,KKF,ISA)
                  XX=100.*HUI(KKF,ISA)
                  WC=.5*(1.-XX/(XX+EXP(SCRP(34,1)-SCRP(34,2)*XX)))+.25
                  ADH2O(IHD)=ADH2O(IHD)+GZSL(KKP,IHD,KKF,ISA)*WC/(1.-WC)+.111*GZSD(KKP,IHD,KKF,ISA)
              END DO !FORAGE LOOP
          END DO ! HERD LOOP
      END DO ! PREFERENCE LOOP
      ! UPDATE STANDING LIVE & STANDING DEAD
      DO I=1,NN
          KKF=JE(I,ISA)
          IF(KKF==MNC)CYCLE
          STL(KKF,ISA)=MAX(0.,STL(KKF,ISA)-YLD(KKF,ISA))
          STD(KKF,ISA)=MAX(0.,STD(KKF,ISA)-YLSD(KKF,ISA))
      END DO                   
      ! ESTIMATE SUPPLEMENTAL FEEDING = DEFICIT
      IF(IHAY>0)THEN 
       DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              IF(IGZR(IHD,ISA)==0)CYCLE
              KKG=NGZ(IHD,ISA)
              IF(KKG==0)CYCLE
        IF(DEF(IHD,IOW)>0. .AND.Hay_total(IHD,IOW)<Hays(IHD,IOW))THEN  
          IF ((MO .GT. Hay_Mo(IHD,YEARK,IOW)) .OR.(MO==Hay_Mo(IHD,YEARK,IOW).AND.KDA.GE.Hay_Day(IHD,YEARK,IOW))&
              .OR.(YEARK.LE.IY .AND. MO.LT.SHay_Mo(IHD,IY,IOW)).OR.&
              (YEARK.LE.IY.AND.MO==SHay_Mo(IHD,IY,IOW).AND.KDA.LE.SHay_Day(IHD,IY,IOW)))THEN
              DEF(IHD,IOW)=DEF(IHD,IOW)-AD2(IHD)
              SPFD(IHD)=MAX(0.,DEF(IHD,IOW))
              Hay_cost(IHD,IOW)=Hay_pri(IHD,IOW)*SPFD(IHD)*WSA(ISA)
              Hay_total(IHD,IOW)=Hay_total(IHD,IOW)+SPFD(IHD)*1000*WSA(ISA)                     !kg/ha----kg/total              
              YTP(2,IHD)=YTP(2,IHD)+SPFD(IHD)
              GZSX(IHD,IOW)=GZSX(IHD,IOW)+SPFD(IHD)
              ATDN(IHD)=ATDN(IHD)+Hay_Q(IHD,IOW)*SPFD(IHD)
              ADNX(IHD)=ADNX(IHD)+.015*SPFD(IHD)
              ADP(IHD)=ADP(IHD)+.0021*SPFD(IHD)
              ADH2O(IHD)=ADH2O(IHD)+.1*SPFD(IHD)
              !FVI(IHD)=(YTP(1,IHD)+.6*YTP(2,IHD))/GZSX(IHD,IOW)
              SMMF(IHD,MO,ISA)=SMMF(IHD,MO,ISA)+SPFD(IHD)
              IF(KFL(31)>0)WRITE(KW(31),49)ISA,NBSA(ISA),IYR,MO,KDA,TIL(JT1),&
                KDC(JJK),II,NBE(JT1),NBT(JT1),Hay_cost(IHD,IOW),SPFD(IHD)*WSA(ISA)
              IF(SPFD(IHD)>.0001)SMMFD(IHD,MO,ISA)=SMMFD(IHD,MO,ISA)+1.
          END IF
        END IF
       END DO
      END IF

      ! OUTPUT GRAZING RESULTS TO .OUT & .DGZ    
      DO KHD=1,NHRD(IOW) ! HERD LOOP
          IHD=KOW(KHD,IOW)
          IF(IGZR(IHD,ISA)==0)CYCLE
          KKG=NGZ(IHD,ISA)
          IF(KKG==0)CYCLE
          RSTX=WSAX/GZNB(IHD,ISA)
          IF(GZSX(IHD,IOW)>0.)THEN
               ATDN(IHD)=ATDN(IHD)/GZSX(IHD,IOW)
               ADNX(IHD)=ADNX(IHD)/GZSX(IHD,IOW)
               ADP(IHD)=ADP(IHD)/GZSX(IHD,IOW)
               ADN(IHD)=ADNX(IHD)
               FECE(IHD,ISA)=MAX(0.01,1000.*GZSX(IHD,IOW)*RSTX*(1.-ATDN(IHD)))   !kg/animal			
             IF (GZDH(IHD,ISA)==1.)THEN
                GZSMD(IHD,IOW)=GZSX(IHD,IOW)*RSTX*1000.                          !kg/animal
                KDAN(IHD,IOW)=KDA
			    TDND(IHD,IOW)=ATDN(IHD)
                !KISA(IHD,IOW)=ISA
                !PRINT*,IY,MO,KDA,IHD,GZSX(IHD,IOW)*RSTX*1000,'***'
             ELSE               
                GZSMH(IHD,IOW)=GZSMH(IHD,IOW)+GZSX(IHD,IOW)*RSTX*1000.            !kg/animal
                TDNH(IHD,IOW)=TDNH(IHD,IOW)+ATDN(IHD)*GZSX(IHD,IOW)*RSTX*1000.    !kg TDN/animal
                TTBS(IHD)=TTBS(IHD) + TBOS(ISA)
                KDAN(IHD,IOW)=KDA 
                !PRINT*,IY,MO,KDA,IHD,GZSX(IHD,IOW)*RSTX*1000,GZSMH(IHD,IOW)
             ENDIF 
          ELSE
              FECE(IHD,ISA)=0.01                                                    !set a low value when GZSX=0
          END IF
          
          ! CALCULATE DRINKING WATER INTAKE l/(hd*d).
          XX=MAX(TX,2.)
          !H2OD=.076*GZWT(IHD,IOW)**.7086*XX**.7146
          H2OD=.076*GZWT(IHD,IOW)**.7086*XX**.7146*GZDH(IHD,ISA)            !for each fasture
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
          !XMLKN=.0054*AD1(IHD)
		  !XMLKN=.0054*AD1(IHD)*GZDH(IHD,ISA)																							  
		  XMLKN=.0054*MILK(IHD,IOW)*GZDH(IHD,ISA)																							  
          XX1=1000.*GZSX(IHD,IOW)*RSTX
          !GZSM(IHD,IOW)=GZSX(IHD,IOW)
          XXN=XX1*ADNX(IHD)-XMLKN
          URNN=FNU*XXN
          !XMLKP=.001*YMLK(IHD)
		  XMLKP=.001*MILK(IHD,IOW)*GZDH(IHD,ISA)																										  
          XXP=XX1*ADP(IHD)-XMLKP
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
          FK(MFT)=FP(MFT)
          SMMMU(IHD,MO,ISA)=SMMMU(IHD,MO,ISA)+FECE(IHD,ISA)
          CALL NFERT(APMU,9,IAMF(ISA),IHD,KT2,JRT)
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
          FK(MFT)=FP(MFT)
          CALL NFERT(APMU,9,IAMF(ISA),IHD,KT2,JRT)
          SMMTDN(IHD,MO,ISA)=SMMTDN(IHD,MO,ISA)+ATDN(IHD)
          SMMGD(IHD,MO,ISA)=SMMGD(IHD,MO,ISA)+1.
          !TDN(IHD,ISA)=ATDN(IHD)
          SUMTDN=0.
          CROPN=0           !add a variable for exact crop Number for each ISA
          ! CALCULATE HERD TDN AND OUTPUT TO .DGZ
          DO I=1,NN ! FORAGE LOOP
              KKF=JE(I,ISA)
              IF(KKF==MNC)CYCLE
              TDN(IHD,ISA)=MAX(TDNN(KKF),ATDN(IHD))
              SUMTDN(ISA)=SUMTDN(ISA)+TDN(IHD,ISA)
              CROPN=CROPN+1
              Y1=1000.*YGZL(IHD,KKF)
              Y2=1000.*YGZD(IHD,KKF)
              SMMSL(IHD,KKF,MO,ISA)=SMMSL(IHD,KKF,MO,ISA)+Y1
              SMMSD(IHD,KKF,MO,ISA)=SMMSD(IHD,KKF,MO,ISA)+Y2
              !IF(Y1<1.E-5.AND.Y2<1.E-5)CYCLE
              X1=1000.*DMDT(IHD,IOW)
              X2=1000.*SPSL(1,IHD,KKF,ISA)
              X3=1000.*SPSL(2,IHD,KKF,ISA)
              X4=1000.*SPSL(3,IHD,KKF,ISA)
              X5=1000.*SPSD(1,IHD,KKF,ISA)
              X6=1000.*SPSD(2,IHD,KKF,ISA)
              X7=1000.*SPSD(3,IHD,KKF,ISA)
              X8=1000.*SPFD(IHD)
              K1=GZNB(IHD,ISA)
              K2=NCLV(IHD)!GZNL
              IF(KFL(18)>0. .and. isa==1)THEN
                  I1=GZNB(IHD,ISA)+.5
                  IF(NDGZ==1)THEN
                      KDG=KDG0+ISA
                      WRITE(KW(KDG),3)ISA,NBSA(ISA),IOW,IYR,MO,KDA,&
                      WSAX,GNAM(KKG),IOW,IHD,GZNB(IHD,ISA),RSTX,GZWT(IHD,IOW),TX,H2OD,&
                      CPNM(KKF),STL(KKF,ISA),HUI(KKF,ISA),CNT(KKF,ISA),STD(KKF,ISA),&
                      BN(3,KKF),X1,TDNF(KKF,ISA),X2,X3,X4,X5,X6,X7,X8,&
                      Y1,Y2,ATDN(IHD),ADNX(IHD),FECE(IHD,ISA),FECN,FECP,&
                      VURN(IHD,ISA),URNN,URNP,ADH2O(IHD),X10,MILK(IHD,IOW),K1,MAX(0,K2),&
                      WTGL(IHD,IOW),DMIN,MAX(0.,CH4E(IHD,IOW)), &
                      SW15,SW30,SNN15+SNA15,SNN30+SNA30,PPL0(KKF,ISA)
                  ELSE 
                      WRITE(KW(18),3)ISA,NBSA(ISA),IOW,IYR,MO,KDA,&
                      WSAX,GNAM(KKG),IOW,IHD,GZNB(IHD,ISA),RSTX,GZWT(IHD,IOW),TX,H2OD,&
                      CPNM(KKF),STL(KKF,ISA),HUI(KKF,ISA),CNT(KKF,ISA),STD(KKF,ISA),&
                      BN(3,KKF),X1,TDNF(KKF,ISA),X2,X3,X4,X5,X6,X7,X8,&
                      Y1,Y2,ATDN(IHD),ADNX(IHD),FECE(IHD,ISA),FECN,FECP,&
                      VURN(IHD,ISA),URNN,URNP,ADH2O(IHD),X10,MILK(IHD,IOW),K1,MAX(0,K2),&
                      WTGL(IHD,IOW),DMIN,MAX(0.,CH4E(IHD,IOW)), &
                      SW15,SW30,SNN15+SNA15,SNN30+SNA30,PPL0(KKF,ISA)
                  END IF
                  PSTL(KKF,ISA)=MAX(PSTL(KKF,ISA),STL(KKF,ISA)*1000.)   !PEAK STL KG/HA
                  PSTD(KKF,ISA)=MAX(PSTD(KKF,ISA),STD(KKF,ISA)*1000.)   !PEAK STD KG/HA
              END IF    
            END DO 
          
          3 FORMAT(1X,2I8,1X,I4,1X,3I4,1X,F10.2,A16,2I4,F10.2,F10.3,F10.0,F10.2,&
      F10.1,3X,A4,27F10.3,2(I8,2X),F10.0,F10.3,F10.0,5(2X,F10.4))      
          HDTM(IHD,ISA)=HDTM(IHD,ISA)+1.
          IF(SUMTDN(ISA)>0.)TDN(IHD,ISA)=SUMTDN(ISA)/REAL(CROPN)
      END DO
      ! UPDATE PLANT STATUS.
      TTOT2(ISA)=0.
      TTOT3(ISA)=0.
      DO I=1,NN ! FORAGE LOOP
          KKF=JE(I,ISA)
          IF(KKF==MNC)CYCLE
          HUF(KKF,ISA)=MAX(HUF(KKF,ISA),HU(KKF,ISA))
          DMF(KKF,ISA)=DM1(KKF,ISA)
          TRA(KKF,ISA)=SRA(KKF,ISA)+TRA(KKF,ISA)
          IF(RD(KKF,ISA)>RDF(KKF,ISA))RDF(KKF,ISA)=RD(KKF,ISA)
          AJHI(KKF,ISA)=0.
          YHE=YLD(KKF,ISA)/HE(JT1)
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
          ZZ5=MAX(.01,1.-ZZ5*YLD(KKF,ISA)/1.5)
		  CPHT(KKF,ISA)=MAX(.001,CPHT(KKF,ISA)*ZZ2)
          HUMIN=MIN(HU(KKF,ISA),.4*PHU(KKF,IHU(KKF,ISA),ISA))
		  HU(KKF,ISA)=MAX(HUMIN,HU(KKF,ISA)-ZZ4)
          SLAI(KKF,ISA)=MAX(.01,SLAI(KKF,ISA)*(1.-X1*1.7))
          Y4=1000.*YLSD(KKF,ISA)*BN(3,KKF)
          Y5=1000.*YLSD(KKF,ISA)*BP(3,KKF)
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
          X5=MIN(YLD(KKF,ISA)*CPLV,UP1(KKF,ISA))
          YLN(KKF,ISA)=MIN(.9*(UN1(KKF,ISA)+STDN(KKF,ISA)),YLD(KKF,ISA)*CNLV)
          YLP(KKF,ISA)=MIN(.9*(UP1(KKF,ISA)+STDP(KKF,ISA)),YLD(KKF,ISA)*CPLV)
          X1=1./HE(JT1)-1.
          XZ1=YLD(KKF,ISA)*X1
          XZ2=MIN(YLSD(KKF,ISA)*X1,STD(KKF,ISA))
          X11=XZ1+XZ2
          STD(KKF,ISA)=MAX(0.,STD(KKF,ISA)-XZ2)
          XZ1=YLN(KKF,ISA)*X1
          XZ2=Y4*X1
          X10=XZ1+XZ2
          STDN(KKF,ISA)=MAX(0.,STDN(KKF,ISA)-XZ2)
          UN1(KKF,ISA)=UN1(KKF,ISA)-XZ1
          TTOT3(ISA)=TTOT3(ISA)+X10
	      JJK=KKF
	      ! ADD TRAMPLED STL AND NUTRIENT CONTENT TO LAYER 1 RESIDUE  
          CALL NCNSTD(X11,X10,LD1)
          DFOP=YLP(KKF,ISA)*X1+Y5*X1
          FOP(LD1,ISA)=FOP(LD1,ISA)+DFOP
          YLD2(KKF,ISA)=YLD2(KKF,ISA)+YLD(KKF,ISA)+YLSD(KKF,ISA)
          JD(ISA)=KKF
          SRA(KKF,ISA)=0.
          UN1(KKF,ISA)=UN1(KKF,ISA)-YLN(KKF,ISA)
          UP1(KKF,ISA)=UP1(KKF,ISA)-DFOP-YLP(KKF,ISA)
          DM(KKF,ISA)=DM(KKF,ISA)-YHE
	      IF(DM(KKF,ISA)<RW(KKF,ISA))RW(KKF,ISA)=DM(KKF,ISA)
          STL(KKF,ISA)=MAX(0.,DM(KKF,ISA)-RW(KKF,ISA))
          YLN(KKF,ISA)=YLN(KKF,ISA)+Y4
          TTOT3(ISA)=TTOT3(ISA)+YLN(KKF,ISA)
          YLP(KKF,ISA)=YLP(KKF,ISA)+Y5
          YLNF(KKF,ISA)=YLNF(KKF,ISA)+YLN(KKF,ISA)
          YLPF(KKF,ISA)=YLPF(KKF,ISA)+YLP(KKF,ISA)
          TYN(ISA)=TYN(ISA)+YLN(KKF,ISA)
          TYP(ISA)=TYP(ISA)+YLP(KKF,ISA)
          !IF(NOP>0.OR.NBSA(ISA)==ISAP)THEN
              !IF(YLD(KKF)>0..OR.YLSD(KKF)>0.)THEN
                  !Y1=1000.*YLD(KKF)
                  !Y2=1000.*YLSD(KKF)
                  !KGZ=1
                  !WRITE(KW(1),29)ISA,NBSA(ISA),IYR,MO,KDA,IOW,TIL(JT1),&
                  !CPNM(KKF),STL(KKF,ISA),STD(KKF,ISA),Y1,Y2,HUI(KKF,ISA),YLN
              !END IF 
          !END IF    
          TTOT2(ISA)=TTOT2(ISA)+UN1(KKF,ISA)+STDN(KKF,ISA)
      END DO
      DF=TTOT1(ISA)-TTOT3(ISA)-TTOT2(ISA)
      !IF(ABS(DF)>.001)WRITE(KW(1),1)IY,MO,KDA,TOT1,TOT2,TOT3,DF
      SMM(141,MO,ISA)=SMM(141,MO,ISA)+1.
      !SUM2=WPML(LD1,ISA)+WPMS(LD1,ISA)+WPMA(LD1,ISA)+WPMU(LD1,ISA)+&
      !WPOU(LD1,ISA)+WPO(LD1,ISA)+FOP(LD1,ISA)+UP1(1,ISA)
      !DF=WSAX*(SUM1-YLP-SUM2)+VAR(56,ISA)+VAR(57,ISA)
      !WRITE(KW(1),'(10E16.6)')WPML(LD1,ISA),WPMS(LD1,ISA),WPMA(LD1,ISA),&
      !WPMU(LD1,ISA),WPOU(LD1,ISA),WPO(LD1,ISA),FOP(LD1,ISA),UP1(1,ISA)
      !WRITE(KW(1),4)IY,MO,KDA,SUM1,YLP,VAR(56,ISA),VAR(57,ISA),SUM2,DF
    !4 FORMAT(1X,'$$$$$',3I4,6E16.6)     
      RETURN
    !1 FORMAT(1X,'!!!!!',3I4,4E12.5)       
    2 FORMAT(1X,2I8,1X,3I4,2X,'GNAM= ',A16,2X,'RSTK= ',F8.3,' ha/hd',2X,&
      'CPNM= ',A4,2X,'STD= ',F7.3,' t/ha',2X,'STL= ',F7.3,' t/ha',2X,&
      'DMD= ',F7.4,' t/ha',2X,'SUPLSL(P-D-U)= ',3F7.3,' t/ha',2X,&
      'SUPLSD(P-D-U)= ',3F7.3,' t/ha',2X,'HAY= ',F7.1,' kg/ha',2X,'YLD= ',&
      F7.1,' kg/ha',2X,'YLSD= ',F7.1,' kg/ha')
																	   
													   
   29 FORMAT(1X,2I8,1X,3I4,2X,'IDON=',I4,2X,A8,2X,A4,2X,'STL=',F7.3,&
      't/ha',2X,'STD=',F7.3,'t/ha',2X,'GZSL=',F7.1,' kg/ha',2X,'GZSD=',&
      F7.1,'kg/ha',2X,'HUSC = ',F6.2,2X,'YLN=',F7.4,' kg/ha')
  101 FORMAT(1X,2I8,1X,3I4,2X,A8,I4,4X,I4,F10.0,F10.2,F10.3,1F10.3,2F10.3,2F10.1,2(2X,I4),5(2X,F10.5))   
   49 FORMAT(1X,2I8,1X,3I4,2X,A8,8X,I6,'   0  ',3I4,F10.3,20X,F10.3)
	  END