      SUBROUTINE NDNITAK(DZ)
!     APEX1905
!     THIS SUBPROGRAM DEVELOPED BY ARMEN KEMANIAN ESTIMATES DAILY 
!     DENITRIFICATION AND N2O LOSSES OF SOIL NO3.
      USE PARM
      DN2O=0.
! compute water factor
      WFP=SWST(ISL,ISA)/PO(ISL,ISA)
      AIRV=MAX(0.,(PO(ISL,ISA)-SWST(ISL,ISA))/DZ)
      X1=.90+.001*CLA(ISL,ISA)
      X2=(1.0001- AIRV)/X1
      IF(X2<.8)RETURN
      H2O_F=1./(1.+X2**(-60))
! compute nitrate factor
      ONO3_C=MAX(1.E-5,1000.*WNO3(ISL,ISA)/WT(ISL,ISA)) ! g/Mg or ppm
      ONO3_F=ONO3_C/(ONO3_C+60.)
! compute respiration factor
      X3=1000.*RSPC(ISL,ISA)/WT(ISL,ISA) ! units mg/kg
      C_F=MIN(1.,X3/50.)
      D_F=ONO3_F*H2O_F*C_F
      DNITMX=PRMT(36)*32.
! units g N / Mg soil / day (32 x 10-6 kg N kg-1 soil day-1)
      WDN=.001*D_F*DNITMX*WT(ISL,ISA) ! kg N /ha
      IF(WDN>WNO3(ISL,ISA))WDN=WNO3(ISL,ISA)
! compute N2O as a fraction of WDN
      DN2O=ONO3_F*(1.-SQRT(H2O_F))*(1.-C_F**.25)*WDN
      WNO3(ISL,ISA)=WNO3(ISL,ISA)-WDN
      RETURN
      END