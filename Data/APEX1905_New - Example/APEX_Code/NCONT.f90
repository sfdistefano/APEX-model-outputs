      SUBROUTINE NCONT(I)
!     APEX1905
!     THIS SUBPROGRAM IS USED TO CHECK N & P BALANCES FOR TROUBLE-
!     SHOOTING PURPOSES.
      USE PARM
      DIMENSION ADX(16)
      DATA YLNX/0./
      WSAX=WSA(I)
      ZEK(I)=0.
      ZFK(I)=0.
      ZSK(I)=0.
      ZNO3(I)=0.
      ZNH3(I)=0.
      ZLSN(I)=0.
      ZLMN(I)=0.
      ZBMN(I)=0.
      ZHSN(I)=0.
      ZHPN(I)=0.
      ZNMU(I)=0.
      ZNOU(I)=0.
      ZLSC(I)=0.
      ZLMC(I)=0.
      ZBMC(I)=0.
      ZHSC(I)=0.
      ZHPC(I)=0.
      ZPML(I)=0.
      ZPO(I)=0.
      ZPMA(I)=0.
      ZPMU(I)=0.
      ZPOU(I)=0.
!     TRSD(I)=0.
      ZFOP(I)=0.
      ZNOS(I)=0.
      ZNOA(I)=0.
      ZPMS(I)=0.
      SW(I)=0.
      SUM=0.
      DO J=1,NBSL(I)
          ISL=LID(J,I)
          ZNO3(I)=ZNO3(I)+WNO3(ISL,I)
          ZNH3(I)=ZNH3(I)+WNH3(ISL,I)
          ZLSN(I)=ZLSN(I)+WLSN(ISL,I)
          ZLMN(I)=ZLMN(I)+WLMN(ISL,I)
          ZBMN(I)=ZBMN(I)+WBMN(ISL,I)
          ZHSN(I)=ZHSN(I)+WHSN(ISL,I)
          ZHPN(I)=ZHPN(I)+WHPN(ISL,I)
          ZNMU(I)=ZNMU(I)+WNMU(ISL,I)
          ZNOU(I)=ZNOU(I)+WNOU(ISL,I)
          ZLSC(I)=ZLSC(I)+WLSC(ISL,I)
          ZLMC(I)=ZLMC(I)+WLMC(ISL,I)
          ZBMC(I)=ZBMC(I)+WBMC(ISL,I)
          ZHSC(I)=ZHSC(I)+WHSC(ISL,I)
          ZHPC(I)=ZHPC(I)+WHPC(ISL,I)
          ZPML(I)=ZPML(I)+WPML(ISL,I)
          ZPMS(I)=ZPMS(I)+WPMS(ISL,I)
          ZPMA(I)=ZPMA(I)+WPMA(ISL,I)
          ZPMU(I)=ZPMU(I)+WPMU(ISL,I)
          ZPOU(I)=ZPOU(I)+WPOU(ISL,I)
          ZPO(I)=ZPO(I)+WPO(ISL,I)
          ZFOP(I)=ZFOP(I)+FOP(ISL,I)
          ZEK(I)=ZEK(I)+EXCK(ISL,I)
          ZFK(I)=ZFK(I)+FIXK(ISL,I)
          ZSK(I)=ZSK(I)+SOLK(ISL,I)
          ! TRSD(I)=TRSD(I)+RSD(ISL,I)
          SW(I)=SW(I)+SWST(ISL,I)
          SUM=SUM+RSPC(ISL,I)
      END DO
      SELECT CASE(IPCD)
          CASE(1) ! C
              SUM=.001*SUM*WSAX
              IF(IEXT(I)>0)THEN
                  X5=.001*VYC(IDOA(I))
                  X4=.001*VQC(IDOA(I))
                  X2=0.
                  X3=0.
              ELSE    
                  X5=.001*VYC(IDOR(I))
                  X4=.001*VQC(IDOR(I))
                  X2=.001*VYC(IDNF(I))
                  X3=.001*VQC(IDNF(I))
              END IF 
              V75=.001*VAR(75,I)
              V140=.001*VAR(140,I)
              V99=.001*VAR(99,I)
              V73=.001*VAR(73,I)
              V101=.001*VAR(101,I)
              V136=.001*VAR(136,I)
              ZOC(I)=ZLSC(I)+ZLMC(I)+ZBMC(I)+ZHSC(I)+ZHPC(I)
              FTC=.001*ZOC(I)*WSAX
              BAL=BTCZ(I)+X2-X5+X3-X4-V136-V75-SUM+V73+V99-V140-V101-FTC
              SBAL=SBAL+BAL
              IF(ABS(BAL)>.01)THEN
                  WRITE(KW(1),125)I,NBSA(I),IYR,MO,KDA,BTCZ(I),ZLSC(I),ZLMC&
                  (I),ZBMC(I),ZHSC(I),ZHPC(I),ZOC(I)
                  WRITE(KW(1),125)I,NBSA(I),IYR,MO,KDA,BTCZ(I),X2,X5,X3,X4,&
                  V136,V75,SUM,V73,V99,V140,V101,FTC,BAL,SBAL
              END IF
              BTCZ(I)=FTC
          CASE(2) ! H2O
              ESALA=SALA(ISA)
              X1=10.*BSALA(ISA)
              WSAX1=10.*WSA(ISA)
              ESW=WSAX1*SW(ISA)
              ESNO=WSAX1*SNO(ISA)
              EGWS=GWST(ISA)
              FSWL=SWLT(ISA)*WSAX1
              X3=EVRS(ISA)
              X4=VAR(4,ISA)
              X11=WSAX1*VAR(11,ISA)
              X72=VAR(72,ISA)
              X16=WSAX1*VAR(16,ISA)
              X71=VAR(71,ISA)
              X18=WSAX1*VAR(18,ISA)
              X110=WSAX1*VAR(110,ISA)
              X19=WSAX1*VAR(19,ISA)
              X86=VAR(86,ISA)
              X146=VAR(146,ISA)
              X147=VAR(147,ISA)
              X156=VAR(156,ISA)
              X157=VAR(157,ISA)
              IF(RSAE(ISA)>0.)THEN
                  X34=VARH(34,IDRO(ISA))
                  X37=VARH(37,IDRO(ISA))
              ELSE
                  X37=VAR(13,ISA)
                  IF(IEXT(ISA)>0)THEN
                      X34=VARH(34,IDOA(ISA))
                      X36=0.
                      X15=VARH(15,IDOA(ISA))
                      X17=0.
                  ELSE
                      X34=VARH(34,IDOR(ISA))
                      X36=VARH(34,IDNF(ISA))
                      X17=VARH(15,IDNF(ISA))
                      X15=VARH(15,IDOR(ISA))
                  END IF    
              END IF 
              CALL HSWBLD(X4,X37,X36,X34,X17,X15,X11,X147,X3,X16,X71,X72,&
              X18,X110,X19,X86,X146,X156,X157,BSNO(ISA),&
              ESNO,BSW(ISA),ESW,BRSV(ISA),RSV(ISA),BGWS(ISA),EGWS,&
              FSWL,BSALA(ISA),ESALA,ISA,NBSA(ISA),IY,MO,KDA,KW(1))
              BSW(ISA)=ESW
              BSNO(ISA)=ESNO
              BRSV(ISA)=RSV(ISA)
              BGWS(ISA)=EGWS
              BSALA(ISA)=ESALA
          CASE(3) ! P
              IF(IEXT(I)>0)THEN
                  X5=VAYP(IDOA(I))
                  X4=VAQP(IDOA(I))
                  X2=0.
                  X3=0.
              ELSE    
                  X5=VAYP(IDOR(I))
                  X4=VAQP(IDOR(I))
                  X2=VAYP(IDNF(I))
                  X3=VAQP(IDNF(I))
              END IF
              YLPX=WSAX*(TYP(I)-YLPX)
              ADD=0.
              TOT=0.
              DO J=1,12
                  TOT=TOT+STDP(J,I)
                  ADD=ADD+UP1(J,I)
              END DO
              V135=VAR(135,ISA)
              V51=VAR(51,ISA)
              V56=VAR(56,ISA)
              V57=VAR(57,ISA)
              V106=VAR(106,ISA)
              FTP=WSAX*(ZPML(I)+ZPMS(I)+ZPMA(I)+ZPO(I)+ZFOP(I)+ZPOU(I)+&
              ZPMU(I)+TOT+ADD)
              BAL=BTPZ(I)+X2-X5+X3-X4-V135-V51-YLPX+V56+V57+V106-V90-FTP
              SBAL=SBAL+BAL
              IF(ABS(BAL)>.001)THEN
                  WRITE(KW(1),3)I,NBSA(I),IY,MO,KDA
                  WRITE(KW(1),'(5X,12E12.5)')BTPZ(I),ZPML(I),ZPMS(I),&
                  ZPMA(I),ZPO(I),ZFOP(I),ZPMU(I),ZPOU(I),TOT,ADD,FTP
                  WRITE(KW(1),2)X2,X5,X3,X4,V135,V51,YLPX,V56,V57,V106,V90,&
                  FTP,BAL,SBAL
              END IF
              V56=0.
              V57=0.
              BTPZ(I)=FTP
              YLPX=TYP(I)
          CASE(4) ! N
              IF(IEXT(I)>0)THEN
                  X5=VAYN(IDOA(I))
                  X4=VAQN(IDOA(I))
                  X15=VSSN(IDOA(I))
                  X17=0.
                  X2=0.
                  X3=0.
              ELSE    
                  X5=VAYN(IDOR(I))
                  X4=VAQN(IDOR(I))
                  X2=VAYN(IDNF(I))
                  X3=VAQN(IDNF(I))
                  X15=VSSN(IDOR(I))
                  X17=VSSN(IDNF(I))
              END IF
              X6=(TYN(I)-YLNX)*WSAX 
              ADD=0.
              TOT=0.
              DO J=1,12
                  TOT=TOT+STDN(J,I)
                  ADD=ADD+UN1(J,I)
              END DO
              ADX(1)=ADX(1)+X5
              ADX(2)=ADX(2)+X4
              ADX(3)=ADX(3)+X15
              ADX(4)=ADX(4)+RFQN(I)
              ADX(5)=ADX(5)+VAR(134,I)
              ADX(6)=ADX(6)+VAR(42,I)
              ADX(7)=ADX(7)+X6
              ADX(8)=ADX(8)+VAR(46,I)
              ADX(9)=ADX(9)+VAR(54,I)
              ADX(10)=ADX(10)+VAR(55,I)
              ADX(11)=ADX(11)+VAR(53,I)
              ADX(12)=ADX(12)+VAR(43,I)
              ADX(13)=ADX(13)+VAR(82,I)
              ADX(14)=ADX(14)+VAR(89,I)
              ADX(15)=ADX(15)+VAR(96,I)
              ADX(16)=ADX(16)+VAR(105,I)
              ZON(I)=ZLSN(I)+ZLMN(I)+ZBMN(I)+ZHSN(I)+ZHPN(I)
              SUM=ZNO3(I)+ZNH3(I)+ZON(I)+TOT+ADD+ZNMU(I)+ZNOU(I)
              FTN=SUM*WSAX+GWSN(I)
              BAL=BTNZ(I)+X2-X5+X3-X4+X17-X15+RFQN(I)-VAR(134,I)-VAR(42,I)&
              -X6-VAR(46,I)+VAR(54,I)+VAR(55,I)+VAR(53,I)+VAR(43,I)&
              -VAR(82,I)-VAR(89,I)-VAR(96,I)+VAR(105,I)-FTN
              SBAL=SBAL+BAL
              IF(ABS(BAL)>0.)THEN
                  WRITE(KW(1),'(A,2I8,3I4)')'N BAL (kg)',I,NBSA(I),IY,MO,KDA
                  WRITE(KW(1),'(5X,12F20.6)')BTNZ(I),ZNO3(I),ZNH3(I),ZON(I),&
                  TOT,ADD,ZNMU(I),ZNOU(I),SUM,GWSN(I),FTN
                  WRITE(KW(1),1)X2,X5,X3,X4,X17,X15,RFQN(I),VAR(134,I),&
                  VAR(42,I),X6,VAR(46,I),VAR(54,I),VAR(55,I),VAR(53,I),&
                  VAR(43,I),VAR(82,I),VAR(89,I),VAR(96,I),VAR(105,I),FTN,&
                  BAL,SBAL
                  WRITE(KW(1),'(10F20.6)')ADX
              END IF
              BTNZ(I)=FTN
              YLNX=TYN(I)
      END SELECT
      RETURN
    1 FORMAT(5X,'  YI=',F20.6,2X,'  YO=',F20.6,2X,'  QI=',F20.6,2X,&
      '  QO=',F20.6,2X,'SSNI=',F20.6,2X,'SSNO=',F20.6/5X,'  RF=',F20.6,&
      2X,'YNWN=',F20.6,2X,'DNIT=',F20.6,2X,' YLN=',F20.6,2X,'AVOL=',&
      F20.6,2X,'FNO3=',F20.6/5X,'FNH3=',F20.6,2X,' FNO=',F20.6,2X,&
      'NFIX=',F20.6,2X,'BURN=',F20.6,2X,'SCRP=',F20.6,2X,'DPKN=',F20.6/&
      5X,'PSON=',F20.6,2X,' FTN=',F20.6,2X,' BAL=',F20.6,2X,'SBAL=',F20.6)  
    2 FORMAT(5X,'  YI=',E12.5,2X,'  YO=',E12.5,2X,'  QI=',E12.5,2X,&
      '  QO=',E12.5,2X,'YPWN=',E12.5/5X,'PLCH=',E12.5,2X,' YLP=',E12.5,&
      2X,'FPMN=',E12.5,2X,' FPO=',E12.5,2X,'PSOP=',E12.5/5X,'SCRP=',&
      E12.5,2X,' FTP=',E12.5,2X,' BAL=',E12.5,2X,'SBAL=',E12.5)  
    3 FORMAT(/T10,'P BAL (kg)',2X,'SA#= ',I8,1X,'ID= ',I8,2X,3I4)
  125 FORMAT(1X,'^^^^^',2I8,1X,3I4,2X,20F10.2)
      END