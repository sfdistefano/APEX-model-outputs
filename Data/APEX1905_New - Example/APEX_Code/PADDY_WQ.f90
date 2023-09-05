      SUBROUTINE PADDY_WQ
!     APEX0806
!     THIS SUBPROGRAM COMPUTES SEDIMENT, AND NUTRIENT FROM A PADDY FIELD .
      USE PARM
      REAL:: QVOL0,SEDCON0,SEDCON1,SEDCON,NMIN0,NMIN1,NMIN,PMIN0,PMIN1,PMIN
      REAL:: NORG0,NORG1,NORG,PORG0,PORG1,PORG,SED0,SED1

      !Tentatively, mineral nutrient concentration is assumed conservative 
      !PADDY_STO(1,ISA) Water storage volume, mm
      !PADDY_STO(2,ISA) Sediment storage, t/ha
      !PADDY_STO(3,ISA) Min N storage, kg/ha
      !PADDY_STO(4,ISA) Min P storage, kg/ha
      !PADDY_STO(5,ISA) Org N storage, kg/ha
      !PADDY_STO(6,ISA) Org P storage, kg/ha

      LD1=LID(1,ISA)
      WSAX=WSA(ISA)
      QVOLX=.1*QVOL(IDO)/WSA(ISA)
      QVOL0 = PADDY_STO(1,ISA) + QVOLX !Water storage at the beginning of the day including daily rainfall and irrigation   
      PADDY_SED_NORM=AUNIF(IDG(13))*3.  !Residual sediment concentration
     
      !Initial sed conc (mg/l)
      SEDCON0 = PADDY_STO(2,ISA) / QVOL0 * 10**5 !mg/l
      IF(SEDCON0<PADDY_SED_NORM) SEDCON0 = PADDY_SED_NORM
      IF(SEDCON0<=0.) SEDCON0 = 0.01
      
      !soluble nutrient concentration
      NMIN0 = PADDY_STO(3,ISA) / QVOL0 * 100.   !Mineral N (mg/l)
      PMIN0 = PADDY_STO(5,ISA) / QVOL0 * 100.   !Mineral P (mg/l)
      
      NORG0 = PADDY_STO(4,ISA) / QVOL0 * 100.  !Organic N (mg/l)  
      PORG0 = PADDY_STO(6,ISA) / QVOL0 * 100.  !Organic P (mg/l)
         
      !Mineral N, P loss due to seepage
      DMN = MIN(PADDY_STO(3,ISA),NMIN0 * PKRZ(LD1) / 100.) !kg/ha
      DMP = MIN(PADDY_STO(4,ISA),PMIN0 * PKRZ(LD1) / 100.) !kg/ha
     
      !Sediment settling Neitsch et al. 2011
      SEDCON1 = (SEDCON0 - PADDY_SED_NORM) * EXP(-0.184 * SUB_D50(ISA)) + PADDY_SED_NORM !mg/l
      IF (SEDCON1<PADDY_SED_NORM) SEDCON1 = PADDY_SED_NORM
      SEDCON = (SEDCON0 + SEDCON1) / 2.
      DSED = MIN(PADDY_STO(2,ISA),(SEDCON0-SEDCON )*PKRZ(LD1)/ 10**5.) !sediment precipitated, t/ha 
         
      !Organic N,P average concentration after settling
      NORG = NORG0 * SEDCON / SEDCON0 !mg/l
      PORG = PORG0 * SEDCON / SEDCON0
      DNORG = MIN(PADDY_STO(5,ISA),NORG * PKRZ(LD1)/ 100.) !KG/HA
      DPORG = MIN(PADDY_STO(6,ISA),PORG * PKRZ(LD1)/ 100.) !KG/HA
      
      !add the settled organic nutrient to the top soil layer
      WNO3(1,ISA) = WNO3(1,ISA) + DMN !Add mineral N to layer 1
      WPMA(1,ISA) = WPMA(1,ISA) + DMP !Add min P to layer 1
      WON(1,ISA) = WON(1,ISA) + DNORG !kg/ha, add to layer 1
      WPO(1,ISA) = WPO(1,ISA) + DPORG !kg/ha, add to layer 1

      PADDY_STO(2,ISA) = PADDY_STO(1,ISA)*SEDCON / 10**5 !t/ha
      PADDY_STO(3,ISA) = PADDY_STO(3,ISA) - DMN  !KG/HA
      PADDY_STO(4,ISA) = PADDY_STO(4,ISA) - DMP
      PADDY_STO(5,ISA) = PADDY_STO(5,ISA) - DNORG
      PADDY_STO(6,ISA) = PADDY_STO(6,ISA) - DPORG

      !Weir discharge if any.
      IF (QVOLX>0) then
         YSD(NDRV,IDO) = MIN(PADDY_STO(2,ISA),SEDCON * QVOLX/ 10**5)*WSAX
         QN(IDO) = MIN(PADDY_STO(3,ISA),NMIN0 * QVOLX / 100.)*0.15*WSAX
         QP(IDO) = MIN(PADDY_STO(4,ISA),PMIN0 * QVOLX / 100.)*0.15*WSAX
         YN(IDO) = MIN(PADDY_STO(5,ISA),NORG * QVOLX / 100.)*WSAX
         YP(IDO) = MIN(PADDY_STO(6,ISA),PORG * QVOLX / 100.)*WSAX
         PADDY_STO(2,ISA) = PADDY_STO(1,ISA)*SEDCON / 10**5  !Update sediment balance
         PADDY_STO(3,ISA) = PADDY_STO(3,ISA) - QN(IDO)/WSAX  !Update min N balance
         PADDY_STO(4,ISA) = PADDY_STO(4,ISA) - QP(IDO)/WSAX  !Update min P balance
         PADDY_STO(5,ISA) = PADDY_STO(5,ISA) - YN(IDO)/WSAX  !Update org N balance
         PADDY_STO(6,ISA) = PADDY_STO(6,ISA) - YP(IDO)/WSAX  !Update org P balance
      else
         YSD(NDRV,IDO) = 0.
         QN(IDO) = 0.
         QP(IDO) = 0.
         YN(IDO) = 0.
         YP(IDO) = 0.
      endif
      
      RETURN
      END