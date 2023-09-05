      MODULE HPERC2_INIT
      CONTAINS
      
      SUBROUTINE FREPARM(CL,OC,SA,FBD,LL,DUL,PO,alpha,n,ths,thr,OC1,IVG)
      !methods for estimating the van Genuchten soil water retention parameters for use with the fastRE soil water balance
      !CL - clay (%)
      !OC - organic carbon (%)
      !SA - sand (%)
      !FBD - bulk density (g/cm3)
      !LL - soil lower limit water content (cm3/cm3)
      !DUL - soil drained upper limit water content (cm3/cm3)
      !PO - maximum saturation water content based on porosity (cm3/cm3)
      !alpha - Van Genuchten alpha parameter (cm water-1)
      !n - Van Genuchten n parameter (unitless)
      !ths - saturation water content (cm3/cm3)
      !thr - residual water content (cm3/cm3)
      !OC1 - soil organic carbon in the top layer (%)
      !IVG - method for estimating Van Genuchten parameters based on traditional EPIC soil parameters. IVG=0 for inputting VG parameters in the soil file
      IMPLICIT NONE
      REAL,INTENT(IN)::CL,OC,SA,FBD,LL,DUL,OC1,PO
      REAL,INTENT(OUT)::alpha,n,ths,thr
      INTEGER,INTENT(IN)::IVG
      REAL::he,OM,SI,PHIE,PHIET,THmin
      INTEGER::topsoil=0,i,chngIndx
      LOGICAL::Coarse
      REAL,DIMENSION(:),ALLOCATABLE::TH,h
      REAL INDX(7)
!      !Vereecken et al., 1989, Soil Science - designed for m=1, shouldn't use this method
!      ths=0.81-0.283*BD+0.001*CL
!      thr=0.015+0.005*CL+0.014*OM
!      alpha=-2.486+0.025*SA-0.351*OM-2.617*BD-0.023*CL
!      alpha=exp(alpha)
!      n=0.053-0.009*SA-0.013*CL+0.00015*SA*SA
!      n=exp(alpha)

      SI=100.-CL-SA
      OM=1.72*OC
      IF(IVG==1)THEN
				 
        !Weynants et al., 2009, Vadose Zone Journal
        ths = 0.6355 + 0.0013*CL - 0.1631*FBD
        ths = min(ths,po)
        alpha = -4.3003 - 0.0097*CL + 0.0138*SA - 0.0992*OM
        alpha = exp(alpha)
        n = -1.0846 - 0.0236*CL - 0.0085*SA + 0.0001*SA*SA
        n = exp(n) + 1.
        thr=0.
      ELSE IF(IVG==2)THEN
        !Wosten et al, 1999, Geoderma
        IF(OC/OC1>0.5) topsoil=1 !need some definition of topsoil vs subsoil
        ths = 0.7919 + 0.001691*CL - 0.29619*FBD - 0.000001491*SI*SI + 0.0000821*OM*OM +  0.02427/CL + 0.01113/SI &
              + 0.01472*log(SI) - 0.0000733*OM*CL - 0.000619*FBD*CL -  0.001183*FBD*OM - 0.0001664*topsoil*SI
        ths = min(ths,po)
        alpha = -14.96 + 0.03135*CL + 0.0351*SI + 0.646*OM + 15.29*FBD - 0.192*topsoil - 4.671*FBD*FBD - 0.000781*CL*CL &
                - 0.00687*OM*OM + 0.0449/OM + 0.0663*log(SI) + 0.1482*log(OM) - 0.04546*FBD*SI - 0.4852*FBD*OM + 0.00673*topsoil*CL
        alpha = exp(alpha)
        n = -25.23 - 0.02195*CL + 0.0074*SI - 0.1940*OM + 45.5*FBD - 7.24*FBD*FBD + 0.0003658*CL*CL + 0.002885*OM*OM - 12.81/FBD &
            - 0.1524/SI - 0.01958*OM - 0.2876*log(SI) - 0.0709*log(OM) - 44.6*log(FBD) - 0.02264*FBD*CL + 0.0896*FBD*OM + 0.00718*topsoil*CL
        n = exp(n) + 1
        thr=0.
      ELSE IF(IVG==3)THEN
        !Jones et al., 2014, Transactions of ASABE - (merge RETC with Saxton and Rawls, 2006)
        CALL TexturalAvg(SA, SI, CL, thr, ths, he, Coarse, alpha, n)
        ths = 0.6355 + 0.0013*CL - 0.1631*FBD !Weynants et al., 2009, Vadose Zone Journal
        ths = min(ths,po)
        thr = min(thr, LL)
        phiet = -21.67*SA/100 - 27.93*CL/100 - 81.97*(ths-DUL) + 71.12*SA/100*(ths-DUL) + 8.29*CL/100*(ths-DUL) + 14.05*(SA/100)*(CL/100) + 27.16
																	
        phie = phiet + 0.02*phiet*phiet - 0.113*phiet - 0.70
        CALL Retention(Coarse,LL,DUL,ths,phie,TH,h) !Saxton and Rawls, 2006, SSSAJ
        INDX=(/0,0,1,1,0,0,0/)
        CALL RETC( DBLE(TH), DBLE(h), DBLE(INDX), ths, thr, alpha, n)  
        DEALLOCATE(h, TH)
      ELSE IF(IVG==4)THEN !this method seems unstable. consider eliminating it if it proves unstable with further testing
																											 
        !merge RETC with Rawls et al., 1982)
        ALLOCATE( h(11), TH(11) )
        CALL TexturalAvg(SA, SI, CL, thr, ths, he, Coarse, alpha, n)
        ths = 0.6355 + 0.0013*CL - 0.1631*FBD !Weynants et al., 2009, Vadose Zone Journal
        ths = min(ths,po)
        thr = min(thr, LL)
        !Rawls et al., 1982, Transactions of ASABE
        allocate(h(11))
        h = (/0.0 ,101.97 ,203.94 ,336.501 ,611.82 ,1019.7 ,2039.4 ,4078.8 ,7137.9 ,10197.0 ,15295.5/) !cm water
!        TH(1) = 0.1829 - 0.0246*OM - 0.0376*FBD + 1.89*DUL - 1.38*LL
  !      TH(2) = 0.8888 - 0.0003*SA - 0.0107*OM + 1.53*DUL - 0.81*LL
!        TH(2) = 0.4829 - 0.0035*SA - 0.0263*OM + 0.25*LL
        TH(1) = ths
        TH(2) = 0.0619 - 0.0002*SA - 0.00067*OM + 1.34*DUL - 0.51*LL
        TH(3) = 0.0319 - 0.0002*SA + 1.01*DUL - 0.06*LL
        TH(4) = 0.2391 - 0.0019*SA + 0.0210*OM + 0.72*LL
        TH(5) = 0.0136 - 0.0091*FBD + 0.66*DUL + 0.39*LL
        TH(6) = -0.0034 + 0.0022*OM + 0.52*DUL + 0.54*LL
        TH(7) = -0.0043 + 0.0026*OM + 0.36*DUL + 0.69*LL
        TH(8) = -0.0038 + 0.0026*OM + 0.24*DUL + 0.79*LL
        TH(9) = -0.0027 + 0.0024*OM + 0.16*DUL + 0.86*LL
        TH(10) = -0.0019 + 0.0022*OM + 0.11*DUL + 0.89*LL
        TH(11) = LL
        IF(Coarse==.TRUE.)chngIndx=2
        IF(Coarse==.FALSE.)THEN
          chngIndx=4
!          TH(2)=max(TH(2),DUL)
!          TH(3)=max(TH(3),DUL)
        END IF
        TH(chngIndx)=DUL
!        DO i=(chngIndx+1), 10
!          TH(i)=min(TH(i),DUL)
!        END DO
        INDX=(/0,0,1,1,0,0,0/)
        CALL RETC( DBLE(TH), DBLE(h), DBLE(INDX), ths, thr, alpha, n)  
        DEALLOCATE(h, TH)
      ELSE 
        THS=PO
        THR=0.
      END IF
      
      THmin = thr + (ths - thr)/(1 + (alpha*1000000)**n)**(1-1/n)
      
      RETURN
      END SUBROUTINE FREPARM
      
      !-------------------------------------------------------------------------------
      SUBROUTINE RETC( &
        SW_Rawls, h_Rawls,INDX &  !Input
        , SAT, wcr, &        !Input or output
        alphaVG_R, nVG_R)    !Output
      !estimates the soil water retention parameters by fitting a water retention model to values of soil water and soil tension
      !here we fit parameters selected by INDX to the Van Genuchten model assuming m=1-1/n
      !SW_Rawls - array of soil water values (cm3/cm3)
      !h_Rawls - array of soil tension values corresponding to SW_Rawls (cm water)
      !INDX - selects which parameters to fit (1) and which to keep as input (0)
      !SAT - saturation water content (cm3/cm3)
      !wcr - residual water content (cm3/cm3)
      !alphaVG_R - Van Genuchten alpha parameter (cm water -1)
      !nVG_R - Van Genuchten n parameter (unitless)
      IMPLICIT NONE
      
      INTEGER i, j, K, kin, kiter, kout, kp, kwater, method, mit, mtype
      Integer IOR, KLOG, NEXP, NIT, NOB, NP, NWC, NWC1, NW
      REAL Ksat, SAT, STOPCR, WCR
      REAL alphaVG_R, mVG_R, nVG_R
      Double Precision alphaVG, ARG, ARG1, ANGLE, DERL, EXPO, GA 
      Double Precision nVG, mVG, RPF, RPZ, RLF, RLX, RLY, RSQ 
      Double Precision SDEV, SECOEF, SSQ, SSQW, SSQ1, SSQ2, SSQW1, SSQW2
      Double Precision SUMB, SUM1, SUM2, SUM3,SUMY, SUMF, SUMYF, SUMY2
      Double Precision SUMF2, SUM, STEP
      Double Precision Temp, TMCOE, TPCOE, TSEC, TVAR,  TVALUE
      Double Precision W1, W2, W12, WA, WB, XLOG, Z 
      Double Precision A(7,7), B(14), D(7,7), DELZ(100,7), E(7), F(200) 
      Double Precision indx(7), P(7), PHI(7), Q(7), R(200), TB(14)
      Double Precision TH(14), W(200), X(200), Y(200)
      CHARACTER:: AB(14)*5
      Double Precision,DIMENSION(:) :: SW_Rawls, h_Rawls 
      DATA STOPCR/.00010/
      
      mtype = 3 !hydraulic conductivity uses Mualem's model, Retention uses VG model with m=1-1/n
      !mtype = 4 !hydraulic conductivity uses Burdine's model, Retention uses VG model with m=1-2/n
      !mtype = 5 !hydraulic conductivity uses Mualem's model, Retention uses Brook & Corey model with m=1-1/n
      !mtype = 6 !hydraulic conductivity uses Burdine's model, Retention uses Brook & Corey model with m=1-2/n 
      method = 5 !outputs diffusivity vs. water content eq. 33 of RETC.pdf
      kwater = 1 !fitting only water retention data
      kin = 1 !number of iterations to be printed
      kout = 1 !tells to print hydraulic properties
      kiter = 1 !number of iteration results to print
      mit = 50 !maximum number of iterations
      AB(8) = 'WCR' !labels of variables
      AB(9) = 'WCS'
      AB(10) = 'ALPHA'
      AB(11) = 'N'
      AB(12) = 'M'
      AB(13) = 'L'
      AB(14) = 'CONDS'
      B(8) = wcr !residual water content
      B(9) = SAT !saturation water content
      B(10) = alphaVG_R !0.059545455  !initial alpha guess - average value of alpha from RETC manual
      B(11) = nVG_R !1.258181818 !initial n guess - average value of n from RETC manual
      B(12) = 1.-1./B(11) !sets m = 1-1/n
      B(13) = .5 !L
      B(14) = 2. !saturated hydraulic conductivity
!      indx(1) = 0 !0-variable will not be fitted 1-variable will be fitted
!      indx(2) = 0
!      indx(3) = 1
!      indx(4) = 1
!      indx(5) = 0 
!      indx(6) = 0
!      indx(7) = 0
      W1 = 0
      kp = 0
    
    ! Rawls data is considered as observation data
      NOB = SIZE(h_Rawls) !number of observations(could include conductivity)
      do i = 1,NOB
        x(i) = h_Rawls(i) !tension data from rawls
        y(i) = SW_Rawls(i) !water content data
        w(i) = 1. !weighting coefficient
      end do
    
      NWC = NOB !number of retention points
      NW = NOB !number of predicted diffusivity data points

!     ----- READ INITIAL ESTIMATES -----
      IF(KWATER.EQ.1) indx(6) = 0
      IF(KWATER.EQ.1) indx(7) = 0
      IF(KWATER.EQ.2.AND.METHOD.LE.2) indx(3) = 0
      IF(KWATER.EQ.2.AND.METHOD.GE.5) indx(3) = 0

!      WRITE(KP,1008) TITLE
      IF(KP.EQ.8) CONTINUE !WRITE(7,1008) TITLE
      GO TO (5,5,1,2,3,3) MTYPE
    1 B(11) = DMAX1(1.05D0,B(11))   ! Parameter n
      B(12) = 1.-1./B(11)           ! Parameter m
      GO TO 4
    2 B(11) = DMAX1(2.05D0,B(11))
      B(12) = 1.-2./B(11)
      GO TO 4
    3 B(11) = DMAX1(0.005,B(11))
      B(12) = 1.0
    4 indx(5) = 0
    5 CONTINUE
      KLOG = 0
      IF(2 * (METHOD/2).EQ.METHOD) KLOG = 1
      IF(KWATER.EQ.1) KLOG = 0
      IF(MIT.EQ.0) GO TO 6
      IF(KWATER.NE.2) GO TO 6
      IF(METHOD.LE.2.OR.METHOD.GT.4) GO TO 6
      indx(1) = 0
      indx(2) = 0
    6 CONTINUE
      NP = 0
      IF(NOB.GT.0) CONTINUE 
      IF(NOB.EQ.0) CONTINUE
      IF(KP.EQ.8) THEN
      IF(KP.EQ.8) CONTINUE 
      IF(NOB.GT.0) CONTINUE 
      IF(NOB.EQ.0) CONTINUE 
      ENDIF
      IF(NOB.EQ.0) GO TO 14

!     ----- WRITE EXPERIMENTAL DATA -----
      WA = 0.
      IF(KWATER.EQ.2) GO TO 8
      DO 7 I = 1,NWC
      X(I) = DMAX1(X(I),1.D-5)
      IF(W(I).LT.1.D-3) W(I) = 1.0
      WA = WA + DABS(W(I) * Y(I))
      IF(KIN.EQ.0) GO TO 7
    7 CONTINUE
      WA = WA/FLOAT(NWC)
      IF(KWATER.EQ.1) GO TO 14
    8 IF(KIN.EQ.0) GO TO 9
    9 WB = 0.0
      IF(KWATER.EQ.2.AND.NWC.EQ.NOB) NWC = 0
      NWC1 = NWC + 1
      DO 10 I = NWC1,NOB
      J = I
      IF(KWATER.EQ.2) J = I-NWC
      X(J) = X(I)
      IF(METHOD.EQ.3.OR.METHOD.EQ.4) X(J) = DMAX1(X(J),1.D-5)
      Y(J) = Y(I)
      IF(KLOG .EQ. 1) Y(J) = DLOG10(Y(J))
      W(J) = W(I)
      IF(W(J) .LT. 1.D-3) W(J) = 1.0
      WB = WB + DABS(W(J) * Y(J))
      IF(KIN.EQ.0) GO TO 10
      IF(KLOG.EQ.0) CONTINUE 
      IF(KP.EQ.8) THEN
      IF(KLOG.EQ.0) CONTINUE 
      IF(KLOG.EQ.1) CONTINUE 
      ENDIF
   10 CONTINUE
      IF(KWATER.LT.2) GO TO 11
      NOB = NOB-NWC
      NWC = 0
      NWC1 = 1
   11 IF(MIT.EQ.0) GO TO 14
      IF(W1.LT.1.D-3) W1 = 1.0
      WB = WB/FLOAT(NOB-NWC)
      W2 = WA/WB
      IF(KWATER.EQ.2) W2 = 1.0
      W12 = W1 * W2
      DO 12 I = NWC1,NOB
   12 W(I) = W12 * W(I)

!      ----- INITIALIZE UNKNOWN PARAMETERS -----
   14 NP = 0
      DO 15 I = 8,14
      TB(I) = B(I) 
      IF(indx(I-7).EQ.0) GO TO 15
      NP = NP + 1
      AB(NP) = AB(I)
      B(NP) = B(I)
      TB(NP) = B(I)
      TH(NP) = B(I)
   15 TH(I) = B(I)
      GA = 0.05
      DERL = 0.002D0
      NEXP = 1 + indx(1) + indx(2) + indx(3) + indx(4) + indx(5)
      IF(KWATER.EQ.1) NOB = NWC
!
!     ----- START LEAST-SQUARES ANALYSIS -----
      CALL MODEL_retc(TH,F,X,NWC,NOB,MTYPE,METHOD,indx,IOR)
     
      IF(IOR.EQ.1) GO TO 94
      IF(MIT.EQ.0) GO TO 83
      SSQ = 0.
      DO 16 I = 1,NOB
      R(I) = W(I) * (Y(I)-F(I))
   16 SSQ = SSQ + R(I) * R(I)
      NIT = 0
!
!     ----- BEGIN OF ITERATION -----
   18 NIT = NIT + 1
      GA = 0.05 * GA
      DO 22 J = 1,NP
      TEMP = TH(J)
      TH(J) = (1.D0 + DERL) * TH(J)
      Q(J) = 0
      CALL MODEL_retc(TH,DELZ(1,J),X,NWC,NOB,MTYPE,METHOD,indx,IOR)
      DO 20 I = 1,NOB
      DELZ(I,J) = W(I) * (DELZ(I,J)-F(I))
   20 Q(J) = Q(J) + DELZ(I,J) * R(I)
      Q(J) = Q(J)/(TH(J) * DERL)
!
!     ----- STEEPEST DESCENT -----
   22 TH(J) = TEMP
      DO 28 I = 1,NP
      DO 26 J = 1,I
      SUM = 0.0
      DO 24 K = 1,NOB
   24 SUM = SUM + DELZ(K,I) * DELZ(K,J)
      D(I,J) = SUM/(TH(I) * TH(J) * DERL**2)
   26 D(J,I) = D(I,J)
   28 E(I) = DSQRT(D(I,I))
   30 DO 32 I = 1,NP
      DO 32 J = 1,NP
   32 A(I,J) = D(I,J)/(E(I) * E(J))
!
!     ----- A IS THE SCALED MOMENT MATRIX -----
      DO 34 I = 1,NP
      P(I) = Q(I)/E(I)
      PHI(I) = P(I)
   34 A(I,I) = A(I,I) + GA
      CALL MATINV(A,NP,P)
!
!     ----- P/E IS THE CORRECTION VECTOR -----
      STEP = 1.0
   36 DO 38 I = 1,NP
   38 TB(I) = P(I) * STEP/E(I) + TH(I)
      DO 40 I = 1,NP
      IF(TH(I) * TB(I))44,44,40
   40 CONTINUE
      CALL MODEL_retc(TB,F,X,NWC,NOB,MTYPE,METHOD,indx,IOR)
      SUMB = 0.0
      DO 42 I = 1,NOB
      R(I) = W(I) * (Y(I)-F(I))
   42 SUMB = SUMB + R(I) * R(I)
   44 SUM1 = 0.0
      SUM2 = 0.0
      SUM3 = 0.0
      DO 46 I = 1,NP
      SUM1 = SUM1 + P(I) * PHI(I)
      SUM2 = SUM2 + P(I) * P(I)
   46 SUM3 = SUM3 + PHI(I) * PHI(I)
      ARG = SUM1/DSQRT(SUM2 * SUM3)
      ARG1 = 0.
      IF(NP.GT.1) ARG1 = DSQRT(1.-ARG * ARG)
      ANGLE = 57.29578 * DATAN2(ARG1,ARG)
!
!     ----------
      DO 48 I = 1,NP
      IF(TH(I) * TB(I))50,50,48
   48 CONTINUE
      IF((SUMB-SSQ)/SSQ.LT.1.D-5) GO TO 56
   50 IF(ANGLE-30.D0) 52,52,54
   52 STEP = 0.5 * STEP
      GO TO 36
   54 GA = 20. * GA
      GO TO 30
!
!     ----- PRINT COEFFICIENTS AFTER EACH ITERATION -----
   56 CONTINUE
      DO 58 I = 1,14
   58 TH(I) = TB(I)
      IF(indx(1).EQ.0) GO TO 60
      IF(NIT.LE.4.OR.TH(1).GT.0.001) GO TO 60
      indx(1) = 0
      B(8) = 0.0
      GO TO 14
   60 IF(indx(6).EQ.0) GO TO 64
      EXPO = TH(NEXP)
      IF(EXPO.GT.1.D-3) GO TO 64
      IF(EXPO.LT.-1.D-3) GO TO 64
      IF(EXPO.LT.0.) GO TO 62
      B(13) = -0.2
      GO TO 14
   62 B(13) = 0.0001
      indx(6) = 0
      GO TO 14
   64 DO 66 I = 1,NP
      IF(DABS(P(I) * STEP/E(I))/(1.D-20 + DABS(TH(I)))-STOPCR) 66,66,68
  66  CONTINUE
      GO TO 70
   68 SSQ = SUMB
      IF(NIT .LE. MIT) GO TO 18

!     ----- END OF ITERATION LOOP -----
   70 CONTINUE
      if(NIT .GE. KITER)then
      alphaVG = TH(10)
      nVG = TH(11)
      mVG = TH(12)
      end if
     
      CALL MATINV(D,NP,P)
!
!     ----- WRITE CORRELATION MATRIX -----
      DO 72 I = 1,NP
   72 E(I) = DSQRT(DMAX1(D(I,I),1.D-20))
      DO 76 I = 1,NP
      DO 74 J = 1,I
   74 A(J,I) = D(J,I)/(E(I) * E(J))
   76 CONTINUE  
!
!     ----- CALCULATE R-SQUARED OF FITTED VS OBSERVED VALUES -----
   78 SUM = 0.0
      SUMY = 0.0
      SUMF = 0.0
      SUMY2 = 0.0
      SUMF2 = 0.0
      SUMYF = 0.0
      DO 80 I = 1,NOB
      SUM = SUM + W(I)
      SUMY = SUMY + Y(I) * W(I)
      SUMF = SUMF + F(I) * W(I)
      SUMY2 = SUMY2 + Y(I) **2 * W(I)
      SUMF2 = SUMF2 + F(I)**2 * W(I)
   80 SUMYF = SUMYF + Y(I) * F(I) * W(I)
      RSQ = (SUMYF-SUMY * SUMF/SUM)**2/ &
        ((SUMY2-SUMY**2/SUM) * (SUMF2-SUMF**2/SUM))
!
!     ----- CALCULATE 95% CONFIDENCE INTERVAL -----
      Z = 1./FLOAT(NOB-NP)
      SDEV = DSQRT(Z * SUMB)
      TVAR = 1.96 + Z * (2.3779+Z*(2.7135+Z*(3.187936 + 2.466666*Z**2)))
      DO 82 I = 1,NP
      SECOEF = E(I) * SDEV
      TVALUE = TH(I)/SECOEF
      TVALUE = DMIN1(TVALUE,999999.D0)
      TSEC = TVAR * SECOEF
      TMCOE = TH(I)-TSEC
82      TPCOE = TH(I) + TSEC
!
!     ----- GIVE FINAL OUTPUT -----
   83   SSQ1 = 0.0
      SSQ2 = 0.0
      SSQW1 = 0.0
      SSQW2 = 0.0
      DO 84 I = 1,NWC
      XLOG = DLOG10(DMAX1(1.D-5,X(I)))
      R(I) = Y(I)-F(I)
      SSQ1 = SSQ1 + R(I)**2
  84  SSQW1 = SSQW1 + (R(I) * W(I))**2
      IF(KWATER.EQ.1) GO TO 89
!
!     ----- WRITE CONDUCTIVITY OR DIFFUSIVITY DATA -----
      DO 88 I = NWC1,NOB
      R(I) = Y(I)-F(I)
      SSQ2 = SSQ2 + R(I)**2
      SSQW2 = SSQW2 + (R(I) * W(I))**2
      RLX = DLOG10(DMAX1(1.D-30,X(I)))
      RLY = DLOG10(DMAX1(1.D-30,Y(I)))
      RPZ = 10.**DMIN1(3.D1,Y(I))
      RLF = DLOG10(DMAX1(1.D-30,F(I)))
      RPF = 10.**DMIN1(3.D1,F(I))
   88 CONTINUE
   89 IF(MIT.EQ.0) GO TO 90
      SSQ = SSQ1 + SSQ2
      SSQW = SSQW1 + SSQW2
   90 CONTINUE
   94 CONTINUE 
      
      wcr = sngl(TH(8))
      SAT = sngl(TH(9))
      alphaVG_R = sngl(TH(10))
      nVG_R = sngl(TH(11))

      RETURN
      END SUBROUTINE RETC
      
      !-------------------------------------------------------------------------------
      SUBROUTINE MODEL_retc(B, &                          !Input
       Y, X, &                                            !Input & output
       NWC, NOB, MTYPE, METHOD, indx, &
       IOR)                                              !Output
!
!     PURPOSE: TO CALCULATE THE HYDRAULIC PROPERTIES
!
      IMPLICIT real*8 (A-H,O-Z)
      Double Precision B(14),Y(200),X(200),indx(7)
      K = 0
      IOR = 0
      DO 2 I = 8,14
      IF(indx(I-7).EQ.0) GO TO 2
      K = K + 1
      B(I) = B(K)
    2 CONTINUE
      WCR = B(8)
      WCS = B(9)
      ALPHA = B(10)
      IND = indx(1) + indx(2) + indx(3) + indx(4)
      IF(MTYPE.EQ.1.OR.MTYPE.EQ.3) B(11) = DMAX1(1.005D0,B(11))
      IF(MTYPE.EQ.3) B(12) = 1.-1./B(11)
      IF(MTYPE.EQ.2.OR.MTYPE.EQ.4) B(11) = DMAX1(2.005D0,B(11))
      IF(MTYPE.EQ.4) B(12) = 1.-2./B(11)
    4 IF(indx(4).EQ.1) B(IND) = B(11)
      RN = B(11)
      RM = B(12)
      EXPO = B(13)
      CONDS = B(14)
      RMN = RM * RN
      DLGA = DLOG10(ALPHA)
      RMT = FLOAT(MTYPE-2 * ((MTYPE-1)/2))
      IF(NOB.EQ.NWC) GO TO 12
!
!     -----CALCULATE COMPLETE BETA FUNCTION----- JZW we do not need
      IF(MTYPE.GT.2) GO TO 10
      AA = RM + RMT/RN
      BB = 1.-RMT/RN
      IF(BB.GT.0.004) GO TO 8
      IOR = 1
      GO TO 60
    8 BETA = GAMMA(AA) * GAMMA(BB)/GAMMA(RM + 1.) !JZW We should not use this
      WCL = DMAX1(2./(2. + RM),0.2D0)
      DLG1 = (3.0-RMT) * DLOG10(RN/(BETA * (RMN + RMT)))
   10 DLG2 = 3.0-RMT + EXPO + 2.0/RMN
      DLG3 = DLOG10(RMN * ALPHA * (WCS-WCR))
      DLG4 = DLOG10(CONDS)
      DLGC = -35.0
      DLGD = -35.0
!
!     ----- CALCULATE FUNCTIONAL VALUES Y(I) -----
   12 DO 54 I = 1,NOB
      IF(METHOD.EQ.3.OR.METHOD.EQ.4) GO TO 13
      IF(I.GT.NWC) GO TO 28
   13 AX = ALPHA * X(I)
      IF(AX.LT.1.D-20) GO TO 16
      EX = RN * DLOG10(AX)
      IF(MTYPE.LT.5) GO TO 14
      IF(AX.LE.1.) GO TO 16
      IF(EX.GT.10.) GO TO 20
      GO TO 22
   14 IF(EX.GT.-10.) GO TO 18
   16 RWC = 1.0
      GO TO 26
   18 IF(EX.LT.10.) GO TO 24
      EX = RM * EX
      IF(EX.LT.30.) GO TO 22
   20 RWC = 0.0
      GO TO 26
   22 RWC = AX**(-RM * RN)
      GO TO 26
   24 RWC = (1. + AX**RN)**(-RM)
   26 Y(I) = WCR + (WCS-WCR) * RWC
      IF(I.LE.NWC) GO TO 54
      GO TO 30
!
!     ----- CONDUCTIVITY DATA -----
   28 RWC = (X(I)-WCR)/(WCS-WCR)
   30 IF(RWC.GT.1.D-10) GO TO 31
      DLGC = -30
      DLGD = -30
      COND = 1.D-30
      DIF = 1.D-30
      GO TO 50
   31 IF(RWC.LT.0.999999D0) GO TO 32
      DLGC = DLG4
      COND = CONDS
      DLGD = 30.0
      DIF = 1.D30
      GO TO 50
   32 DLGW = DLOG10(RWC)
      DLGC = DLG2 * DLGW + DLG4
      DLGD = DLGC-DLG3-(RMN + 1) * DLGW/RMN
      IF(DLGC.LT.-30..OR.DLGW.LT.(-15. * RM)) GO TO 48
      IF(MTYPE.GT.4) GO TO 46
      DW = RWC**(1./RM)
      IF(MTYPE.GT.2) GO TO 42
!
!     ----- MTYPE = 1 OR 2 (VARIABLE M,N) -----
      IF(DW.GT.1.D-06) GO TO 34
      DLGC = DLGC + DLG1
      DLGD = DLGC-DLG3-(RMN + 1.) * DLGW/RMN
      GO TO 48
   34 IF(RWC-WCL) 36,36,38
   36 TERM = BINC(DW,AA,BB,BETA) !JZW we should not use this
      GO TO 44
   38 TERM = 1.-BINC(1.-DW,BB,AA,BETA)
      GO TO 44
!
!     ----- MTYPE = 3 OR 4 (RESTRICTED M,N) -----
   42 A = DMIN1(0.999999D0,DMAX1(1.D-7,1.-DW))
      TERM = 1.D0-A**RM
      IF(DW.LT.1.D-04) TERM = RM * DW * (1.-0.5 * (RM-1.) * DW)
   44 RELK = RWC**EXPO * TERM
      IF(RMT.LT.1.5) RELK = RELK * TERM
      DLGC = DLOG10(RELK) + DLG4
      DLGD = DLGC-DLG3-(RMN + 1.) * DLGW/RMN-(RN-1.) * DLOG10(1.-DW)/RN
      GO TO 48
!
!     ----- MTYPE = 5 OR 6 -----
   46 DLGD = DLG4-DLG3 + (2.0-RMT + EXPO + 1./RN) * DLGW
   48 DLGC = DMAX1(-30.D0,DLGC)
      DLGD = DMAX1(-30.D0,DLGD)
      DLGD = DMIN1(30.D0,DLGD)
      COND = 10.**DLGC
      DIF = 10.**DLGD
   50 IF(METHOD.EQ.1.OR.METHOD.EQ.3) Y(I) = COND
      IF(METHOD.EQ.2.OR.METHOD.EQ.4) Y(I) = DLGC
      IF(METHOD.EQ.5) Y(I) = DIF
      IF(METHOD.EQ.6) Y(I) = DLGD
 1000 FORMAT(I5,6D13.5)
   54 CONTINUE
   60 CONTINUE
      RETURN
      END SUBROUTINE MODEL_retc
      
      SUBROUTINE MATINV(A,NP,B)
      !TO INVERT THE MATRIX FOR PARAMETER ESTIMATION
      IMPLICIT real*8 (A-H,O-Z)
      Double Precision A(7,7),B(7),indx(7,2)
      DO 2 J = 1,7
    2 indx(J,1) = 0
      I = 0
    4 AMAX = -1.0
      DO 12 J = 1,NP
      IF(indx(J,1)) 12,6,12
    6 DO 10 K = 1,NP
      IF(indx(K,1)) 10,8,10
    8 P = ABS(A(J,K))
      IF(P.LE.AMAX) GO TO 10
      IR = J
      IC = K
      AMAX = P
   10 CONTINUE
   12 CONTINUE
      IF(AMAX) 30,30,14
   14 indx(IC,1) = IR
      IF(IR.EQ.IC) GO TO 18
      DO 16 L = 1,NP
      P = A(IR,L)
      A(IR,L) = A(IC,L)
   16 A(IC,L) = P
      P = B(IR)
      B(IR) = B(IC)
      B(IC) = P
      I = I + 1
      indx(I,2) = IC
   18 P = 1./A(IC,IC)
      A(IC,IC) = 1.0
      DO 20 L = 1,NP
   20 A(IC,L) = A(IC,L) * P
      B(IC) = B(IC) * P
      DO 24 K = 1,NP
      IF(K.EQ.IC) GO TO 24
      P = A(K,IC)
      A(K,IC) = 0.0
      DO 22 L = 1,NP
   22 A(K,L) = A(K,L)-A(IC,L) * P
      B(K) = B(K)-B(IC) * P
   24 CONTINUE
      GO TO 4
   26 IC = indx(I,2)
      IR = indx(IC,1)
      DO 28 K = 1,NP
      P = A(K,IR)
      A(K,IR) = A(K,IC)
   28 A(K,IC) = P
      I = I-1
   30 IF(I) 26,32,26
   32 RETURN
      END SUBROUTINE MATINV
      
      FUNCTION GAMMA(Z)
      !TO CALCULATE THE GAMMA FUNCTION FOR POSITIVE Z
      IMPLICIT real*8 (A-H,O-Z)
      IF(Z.LT.33.) GO TO 2
      GAMMA = 1.D36
      RETURN
    2 X = Z
      GAMMA = 1.0
      IF(X-2.0) 10,10,8
    6 IF(X-2.0) 14,14,8
    8 X = X-1.0
      GAMMA = GAMMA * X
      GO TO 6
   10 IF(X-1.0) 12,16,14
   12 GAMMA = GAMMA/X
      X = X + 1.0
   14 Y = X-1.0
      FY = 1.0-Y * (.5771017-Y * (.985854-Y * (.8764218- &
       Y * (.8328212-Y * (.5684729-Y * (.2548205-.0514993 * Y))))))
      GAMMA = GAMMA * FY
   16 RETURN
      END FUNCTION GAMMA
      
      FUNCTION BINC(X,A,B,BETA)
      !TO CALCULATE THE INCOMPLETE BETA-FUNCTION
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION T(200)
      DATA NT/10/
      NT1 = NT + 1
      T(1) = -(A + B) * X/(A + 1.0)
      DO 2 I = 2,NT,2
      Y = FLOAT(I/2)
      Y2 = FLOAT(I)
      T(I) = Y * (B-Y) * X/((A + Y2-1.0) * (A + Y2))
    2 T(I + 1) = -(A + Y) * (A + B + Y) * X/((A + Y2) * (A + Y2 + 1.0))
      BINC = 1.0
      DO 4 I = 1,NT
      K = NT1-I
    4 BINC = 1. + T(K)/BINC
      BINC = X**A * (1.-X)**B/(BINC * A * BETA)
      RETURN
      END FUNCTION BINC
      
      !-------------------------------------------------------------------------------
      SUBROUTINE TexturalAvg(SAND, SILT, CLAY,  &   !Input 
      thr, ths, he, Coarse, alpha, nvg, LL, DUL, ksat, Texture) !Output
      !calculates average soil parameters based on texture
      !SAND - sand content (%)
      !SILT - silt content (%)
      !CLAY - clay content (%)
      !thr - residual water content (cm3/cm3)
      !ths - saturation water content (cm3/cm3)
      !he - air entry tension (cm water)
      !Coarse - indicates if soil is of coarse texture
      !alpha - Van Genuchten alpha parameter (cm water -1)
      !nvg - Van Genuchten n parameter (unitless)
      !LL - soil lower limit water content (cm3/cm3)
      !DUL - soil drained upper limit water content (cm3/cm3)
      !ksat - soil saturated hydraulic conductivity (mm/hr)
      !texture - soil texture
      IMPLICIT NONE
      REAL :: SAND, SILT, CLAY, thr, ths, he
      LOGICAL :: Coarse
      REAL,OPTIONAL :: alpha, nvg, LL, DUL, ksat
      CHARACTER*12,OPTIONAL :: Texture
      
      !values from Rawls et al, 1982, Transactions of ASABE
      !alpha and n values from Rosetta Lite
      IF (SAND .GE. 85. .AND. SILT + 1.5 * CLAY .LE. 15.) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'Sand        '
        COARSE = .TRUE.
        thr=0.020
        ths=0.417
        he=7.26
        IF(PRESENT(alpha)) alpha=0.0353
        IF(PRESENT(nvg)) nvg=3.1798
        IF(PRESENT(DUL)) DUL=0.091
        IF(PRESENT(LL)) LL=0.033
        IF(PRESENT(ksat))ksat=210.0
      ELSEIF ((SAND .GE. 85. .AND. SAND .LT. 90. .AND. SILT + 1.5 * CLAY .GE. 15.) .OR. & 
             (SAND .GE. 70. .AND. SAND .LT. 85. .AND. SILT + 2.0 * CLAY .LE. 30.)) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'LoamySand   '
        COARSE = .TRUE.
        thr=0.035
        ths=0.401
        he=8.69
        IF(PRESENT(alpha)) alpha=0.0347
        IF(PRESENT(nvg)) nvg=1.7466
        IF(PRESENT(DUL)) DUL=0.125
        IF(PRESENT(LL)) LL=0.055
        IF(PRESENT(ksat))ksat=61.11
      ELSEIF ((CLAY .LE. 20. .AND. SAND .GE. 52. .AND. SILT + 2. * CLAY .GT. 30.) .OR. &
            (CLAY .LT. 7. .AND. SILT .LT. 50. .AND. SAND .GT. 43. .AND. SAND .LT. 52.)) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'SandyLoam   '
        COARSE = .TRUE.
        thr=0.041
        ths=0.412
        he=14.66
        IF(PRESENT(alpha)) alpha=0.0267
        IF(PRESENT(nvg)) nvg=1.4484
        IF(PRESENT(DUL)) DUL=0.207
        IF(PRESENT(LL)) LL=0.095
        IF(PRESENT(ksat))ksat=25.9
      ELSEIF (CLAY .GE. 7. .AND. CLAY .LE. 27. .AND. SILT .GE. 28. .AND. SILT .LT. 50. .AND. SAND .GE. 23. .AND. SAND .LT. 52.) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'Loam        '
        COARSE = .FALSE.
        thr=0.027
        ths=0.434
        he=11.15
        IF(PRESENT(alpha)) alpha=0.0111
        IF(PRESENT(nvg)) nvg=1.4737
        IF(PRESENT(DUL)) DUL=0.270
        IF(PRESENT(LL)) LL=0.117
        IF(PRESENT(ksat))ksat=13.2
      ELSEIF ((SILT .GE. 50. .AND. SILT .LE. 88. .AND. CLAY .GE. 12. .AND. CLAY .LE. 27.) .OR. &
             (SILT .GE. 50. .AND. SILT .LT. 80. .AND. CLAY .GE. 0. .AND. CLAY .LT. 12.)) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'SiltyLoam   '
        COARSE = .FALSE.
        thr=0.015
        ths=0.486
        he=20.76
        IF(PRESENT(alpha)) alpha=0.0051
        IF(PRESENT(nvg)) nvg=1.6626
        IF(PRESENT(DUL)) DUL=0.330
        IF(PRESENT(LL)) LL=0.133
        IF(PRESENT(ksat))ksat=6.8
      ELSEIF (SILT .GE. 80. .AND. CLAY .LT. 12.) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'Silt        '
        COARSE = .FALSE.
        thr=0.015
        ths=0.486
        he=20.76
        IF(PRESENT(alpha)) alpha=0.0066
        IF(PRESENT(nvg)) nvg=1.6769
        IF(PRESENT(DUL)) DUL=0.330
        IF(PRESENT(LL)) LL=0.133
        IF(PRESENT(ksat))ksat=6.8
      ELSEIF (CLAY .GE. 20. .AND. CLAY .LT. 35. .AND. SILT .GE. 0.  .AND. SILT .LT. 28. .AND. SAND .GE. 45.) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'SandClayLoam'
        COARSE = .FALSE.
        thr=0.068
        ths=0.330
        he=28.08
        IF(PRESENT(alpha)) alpha=0.0211
        IF(PRESENT(nvg)) nvg=1.3298
        IF(PRESENT(DUL)) DUL=0.255
        IF(PRESENT(LL)) LL=0.148
        IF(PRESENT(ksat))ksat=4.3
      ELSEIF (CLAY .GE. 27. .AND. CLAY .LT. 40. .AND. SAND .GE. 20. .AND. SAND .LT. 45.) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'ClayLoam    '
        COARSE = .FALSE.
        thr=0.075
        ths=0.390
        he=25.89
        IF(PRESENT(alpha)) alpha=0.0158
        IF(PRESENT(nvg)) nvg=1.4145
        IF(PRESENT(DUL)) DUL=0.318
        IF(PRESENT(LL)) LL=0.197
        IF(PRESENT(ksat))ksat=2.3
      ELSEIF (CLAY .GE. 27. .AND. CLAY .LT. 40. .AND. SAND .GE. 0. .AND. SAND .LT. 20.) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'SiltClayLoam'
        COARSE = .FALSE.
        thr=0.040
        ths=0.432
        he=32.56
        IF(PRESENT(alpha)) alpha=0.0084
        IF(PRESENT(nvg)) nvg=1.5202
        IF(PRESENT(DUL)) DUL=0.366
        IF(PRESENT(LL)) LL=0.208
        IF(PRESENT(ksat))ksat=1.5
      ELSEIF (CLAY .GE. 35. .AND. SAND .GE. 45.) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'SandyClay   '
        COARSE = .FALSE.
        thr=0.109
        ths=0.321
        he=29.17
        IF(PRESENT(alpha)) alpha=0.0334
        IF(PRESENT(nvg)) nvg=1.2067
        IF(PRESENT(DUL)) DUL=0.339
        IF(PRESENT(LL)) LL=0.239
        IF(PRESENT(ksat))ksat=1.2
      ELSEIF (CLAY .GE. 40. .AND. SILT .GE. 40.) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'SiltyClay   '
        COARSE = .FALSE.
        thr=0.056
        ths=0.423
        he=34.19
        IF(PRESENT(alpha)) alpha=0.0162
        IF(PRESENT(nvg)) nvg=1.3207
        IF(PRESENT(DUL)) DUL=0.387
        IF(PRESENT(LL)) LL=0.250
        IF(PRESENT(ksat))ksat=0.9
      ELSEIF (CLAY .GE. 40. .AND. SAND .LT. 45. .AND. SILT .LT. 40.) THEN
        IF(PRESENT(TEXTURE)) TEXTURE = 'Clay        '
        COARSE = .FALSE.
        thr=0.090
        ths=0.385
        he=37.30
        IF(PRESENT(alpha)) alpha=0.0150
        IF(PRESENT(nvg)) nvg=1.2529
        IF(PRESENT(DUL)) DUL=0.396
        IF(PRESENT(LL)) LL=0.272
        IF(PRESENT(ksat))ksat=0.6
      ELSE
        IF(PRESENT(TEXTURE)) TEXTURE = 'UNKOWN      '
        COARSE = .FALSE.
        thr=0.052
        ths=0.403
        he=22.70
        IF(PRESENT(alpha)) alpha=0.0158
        IF(PRESENT(nvg)) nvg=1.4145
        IF(PRESENT(DUL)) DUL=0.280
        IF(PRESENT(LL)) LL=0.159
        IF(PRESENT(ksat))ksat=29.8
      END IF
      
      RETURN
      END SUBROUTINE TexturalAvg
      
      !-------------------------------------------------------------------------------
      SUBROUTINE Retention(Coarse, LL, DUL, SAT, H_bub &   !Input 
      , TH, h)                                        !Output
      !calculates points along the soil water retention curve based on the soil properties
      !Coarse - if the soil is of coarse texture
      !LL - lower limit soil water content (cm3/cm3)
      !DUL - drained upper limit water content (cm3/cm3)
      !SAT - saturation water content (cm3/cm3)
      !H_bub - tension at air entry (cm water)
      !TH - water content (cm3/cm3)
      !h - tension (cm water)
      IMPLICIT NONE
      LOGICAL::Coarse
      INTEGER :: i,steps,bubSteps
      REAL :: Aparam, Bparam, tempStep, hDUL=33., H_bub, LL, DUL, SAT
      REAL,DIMENSION(:),ALLOCATABLE :: h, TH
      
      steps = CEILING( (SAT-LL)/0.005 )
      bubSteps = MAX( CEILING(steps*.05), 2 )
      ALLOCATE( h(steps+bubSteps), TH(steps+bubSteps) )
      
      !Relationships from Saxton and Rawls, 2006, SSSAJ
      IF (Coarse) hDUL=10.
      Bparam = (log(1500.) - log(hDUL)) / (log(DUL) - log(LL))
      Aparam = exp( log(hDUL) + ( Bparam * log(DUL) ) )
      
      tempStep = (SAT-LL)/steps
      DO i=1,steps
        TH(i) = LL + tempStep*(i-1)
        IF(TH(i) .LE. DUL)THEN
          h(i) = Aparam * ( TH(i) ** (- Bparam ) ) !1500 to hDUL (kPa)
        ELSE
          h(i) =hDUL-((TH(i) - DUL) * (hDUL- H_bub )/(SAT - DUL)) !hDUL to H_bub (kPa)
        END IF
      END DO
      
      !H_bub to hSat (kPa)
      tempStep = (H_bub-1)/bubSteps
      do i = steps+1, steps+bubSteps
         TH(i) = SAT
         h(i) = H_bub - tempStep*(i-steps)
      end do
      
      h = h*10.197 !convert from kPa to cm water
      
      RETURN
      END SUBROUTINE Retention
      
      END MODULE HPERC2_INIT