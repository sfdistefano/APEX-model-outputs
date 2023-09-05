      SUBROUTINE EWER(JRT)
!     APEX1905
!     THIS SUBPROGRAM ESTIMATES DAILY SOIL LOSS CAUSED BY WIND EROSION,
!     GIVEN THE AVERAGE WIND SPEED AND DIRECTION.
      USE PARM
      JRT=0
      IF(U10(IRF(ISA))<6..OR.PADDY_HWEIR(ISA)>0.)THEN
          JRT=1
          RETURN
      ELSE
          BT=PI2/4.+THW-ANG
          ALG=FL*FW/(FL*ABS(COS(BT))+FW*ABS(SIN(BT)))
          RRF=11.9*(1.-EXP(-(RRUF(ISA)/9.8)**1.3))
          RIF=ABS(SIN(BT))*1.27*RHTT(ISA)**.52
          RFB=MAX(.1,RRF+RIF)
          RFC=.77*1.002**RHTT(ISA)
          RGRF=1.-EXP(-(10./RFB)**RFC)
          X1=VAC(ISA)
          X1=MIN(10.,X1+BWN(3,JD(ISA))*RSD(LID(1,ISA),ISA))
          IF(SCLM(13)>0.)X1=MIN(X1,SCLM(13))
          F=X1/(X1+EXP(SCRP(13,1)-SCRP(13,2)*X1))
          VGF=(1.-F)*SRNG(13)+SCMN(13)
          ALG=1.-EXP(-ALG/.07)
          RTO=SWST(LID(1,ISA),ISA)/S15(LID(1,ISA),ISA)
          SMM(60,MO,ISA)=SMM(60,MO,ISA)+RTO
          NWDA(ISA)=NWDA(ISA)+1
          VAR(60,ISA)=RTO
          USTW=USTT+.5*RTO*RTO
          CALL EWNINT
          ROKF=EXP(-.047*ROK(LID(1,ISA),ISA))
          X1=PRMT(13)*TLMF(ISA)
          IF(X1>10.)THEN
              WK1=1.E-5
          ELSE
              WK1=WK(ISA)*EXP(-X1)
          END IF
          SMM(32,MO,ISA)=SMM(32,MO,ISA)+WK1
          VAR(32,ISA)=WK1
          X1=RGRF*VGF*ROKF*WK1
          IF(YW(IDO)<.001)YW(IDO)=0.   !MOVED BEFORE YW(IDO)=8640.*YW(IDO)*X1 BY LIWANG MA
          YW(IDO)=8640.*YW(IDO)*X1
          TLMF(ISA)=TLMF(ISA)+YW(IDO)
          TVGF(ISA)=TVGF(ISA)+X1
      END IF
      RETURN
      END