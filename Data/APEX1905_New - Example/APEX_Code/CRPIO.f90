      SUBROUTINE CRPIO
!     THIS SUBPROGRAM WRITES INPUT CROP PARAMETERS.
      USE PARM
      IF(KFL(1)>0)THEN
          CALL APAGE(0)
          WRITE(KW(1),11)
          WRITE(KW(1),1)(CPNM(I),I=1,LC)
          WRITE(KW(1),6)'WA',(WA(I),I=1,LC)
          WRITE(KW(1),7)'HI',(HI(I),I=1,LC)
          WRITE(KW(1),6)'TB',(TOPC(I),I=1,LC)
          WRITE(KW(1),6)'TG',(TBSC(I),I=1,LC)
          WRITE(KW(1),7)'DMLA',(DMLA(I),I=1,LC)
          WRITE(KW(1),7)'DLAI',(DLAI(I),I=1,LC)
          WRITE(KW(1),8)'LAP1',(DLAP(1,I),I=1,LC)
          WRITE(KW(1),8)'LAP2',(DLAP(2,I),I=1,LC)
          WRITE(KW(1),8)'PPL1',(PPLP(1,I),I=1,LC)
          WRITE(KW(1),8)'PPL2',(PPLP(2,I),I=1,LC)
          WRITE(KW(1),8)'FRS1',(FRST(1,I),I=1,LC)
          WRITE(KW(1),8)'FRS2',(FRST(2,I),I=1,LC)
          WRITE(KW(1),7)'RLAD',(RLAD(I),I=1,LC)
          WRITE(KW(1),7)'RBMD',(RBMD(I),I=1,LC)
          WRITE(KW(1),6)'ALT',(ALT(I),I=1,LC)
          WRITE(KW(1),7)'CAF',(CAF(I),I=1,LC)
          WRITE(KW(1),9)'GSI',(GSI(I),I=1,LC)
          WRITE(KW(1),7)'WAC2',(WAC2(2,I),I=1,LC)
          WRITE(KW(1),6)'WAVP',(WAVP(I),I=1,LC)
          WRITE(KW(1),6)'VPTH',(VPTH(I),I=1,LC)
          WRITE(KW(1),7)'VPD2',(VPD2(I),I=1,LC)
          WRITE(KW(1),6)'SDW',(SDW(I),I=1,LC)
          WRITE(KW(1),7)'HMX',(HMX(I),I=1,LC)
          WRITE(KW(1),7)'RDMX',(RDMX(I),I=1,LC)
          WRITE(KW(1),8)'RWP1',(RWPC(1,I),I=1,LC)
          WRITE(KW(1),8)'RWP2',(RWPC(2,I),I=1,LC)
          WRITE(KW(1),6)'GMHU',(GMHU(I),I=1,LC)      
          WRITE(KW(1),9)'CNY',(CNY(I),I=1,LC)
          WRITE(KW(1),9)'CPY',(CPY(I),I=1,LC)
          WRITE(KW(1),9)'CKY',(CKY(I),I=1,LC)
          WRITE(KW(1),9)'WSYF',(WSYF(I),I=1,LC)
          WRITE(KW(1),7)'PST',(PST(I),I=1,LC)
          WRITE(KW(1),7)'CSTS',(CSTS(I),I=1,LC)
          WRITE(KW(1),7)'PRYG',(PRYG(I),I=1,LC)
          WRITE(KW(1),7)'PRYF',(PRYF(I),I=1,LC)
          WRITE(KW(1),7)'WCY',(WCY(I),I=1,LC)
          WRITE(KW(1),9)'BN1',(BN(1,I),I=1,LC)
          WRITE(KW(1),9)'BN2',(BN(2,I),I=1,LC)
          WRITE(KW(1),9)'BN3',(BN(3,I),I=1,LC)
          WRITE(KW(1),9)'BP1',(BP(1,I),I=1,LC)
          WRITE(KW(1),9)'BP2',(BP(2,I),I=1,LC)
          WRITE(KW(1),9)'BP3',(BP(3,I),I=1,LC)
          WRITE(KW(1),9)'BK1',(BK(1,I),I=1,LC)
          WRITE(KW(1),9)'BK2',(BK(2,I),I=1,LC)
          WRITE(KW(1),9)'BK3',(BK(3,I),I=1,LC)
          WRITE(KW(1),8)'BW1',(BWN(1,I),I=1,LC)
          WRITE(KW(1),8)'BW2',(BWN(2,I),I=1,LC)
          WRITE(KW(1),8)'BW3',(BWN(3,I),I=1,LC)
          WRITE(KW(1),8)'STX1',(STX(1,I),I=1,LC)
          WRITE(KW(1),8)'STX2',(STX(2,I),I=1,LC)
          WRITE(KW(1),8)'BLG1',(BLG(1,I),I=1,LC)
          WRITE(KW(1),8)'BLG2',(BLG(2,I),I=1,LC)
          WRITE(KW(1),8)'FTO ',(FTO(I),I=1,LC)
          WRITE(KW(1),8)'FLT ',(FLT(I),I=1,LC)
          WRITE(KW(1),8)'TDNX',(TDNX(I),I=1,LC)
          WRITE(KW(1),8)'TDNN',(TDNN(I),I=1,LC)
          WRITE(KW(1),8)'ANQX',(ANTQX(I),I=1,LC)
          WRITE(KW(1),8)'ANQN',(ANTQN(I),I=1,LC)
          WRITE(KW(1),8)'DETN',(DETIN(I),I=1,LC)
          WRITE(KW(1),8)'RDAT',(RDEAT(I),I=1,LC)
          WRITE(KW(1),8)'EXTC',(EXTC(I),I=1,LC)
          WRITE(KW(1),10)'IDC',(IDC(I),I=1,LC)
      END IF    
      IPL=0
      DO J=1,LC
          BLG(3,J)=BLG(2,J)
          BLG(1,J)=BLG(1,J)/BLG(2,J)
          BLG(2,J)=.99
          CALL ASCRV(BLG(1,J),BLG(2,J),.5,1.,SCLM,32,KW(1))
          IF(NUPC>0)THEN
              BN(4,J)=BN(1,J)
              X1=BN(1,J)-BN(3,J)
              BN(1,J)=1.-(BN(2,J)-BN(3,J))/X1
              BN(2,J)=1.-.00001/X1
              CALL ASCRV(BN(1,J),BN(2,J),.5,1.,SCLM,33,KW(1))
              BP(4,J)=BP(1,J)
              X1=BP(1,J)-BP(3,J)
              BP(1,J)=1.-(BP(2,J)-BP(3,J))/X1
              BP(2,J)=1.-.00001/X1
              CALL ASCRV(BP(1,J),BP(2,J),.5,1.,SCLM,34,KW(1))
              BK(4,J)=BK(1,J)
              X1=BK(1,J)-BK(3,J)
              BK(1,J)=1.-(BK(2,J)-BK(3,J))/X1
              BK(2,J)=1.-.00001/X1
              CALL ASCRV(BK(1,J),BK(2,J),.5,1.,SCLM,35,KW(1))
          ELSE
              CALL NCONC(BN(1,J),BN(2,J),BN(3,J),BN(4,J),KW,MSO)
              CALL NCONC(BP(1,J),BP(2,J),BP(3,J),BP(4,J),KW,MSO)
              CALL NCONC(BK(1,J),BK(2,J),BK(3,J),BK(4,J),KW,MSO)
          END IF    
          X1=ASPLT(DLAP(1,J))*.01
          X2=ASPLT(DLAP(2,J))*.01
          CALL ASCRV(DLAP(1,J),DLAP(2,J),X1,X2,SCLM,36,KW(1))
          X1=ASPLT(FRST(1,J))
          X2=ASPLT(FRST(2,J))
          CALL ASCRV(FRST(1,J),FRST(2,J),X1,X2,SCLM,37,KW(1))
          WAC2(1,J)=.01
          X2=ASPLT(WAC2(2,J))
          WAX(J)=100.*WAC2(2,J)-WAI(J)
          WAC2(2,J)=.9
          CALL ASCRV(WAC2(1,J),WAC2(2,J),330.,X2,SCLM,38,KW(1))
          X2=ASPLT(VPD2(J))
          VPD2(J)=(1.-VPD2(J))/(X2-VPTH(J))
      END DO
      IF(KFL(1)>0)THEN
          CALL APAGE(0)
          WRITE(KW(1),11)
          WRITE(KW(1),1)(CPNM(I),I=1,LC)
          WRITE(KW(1),9)'LAP1',(DLAP(1,I),I=1,LC)
          WRITE(KW(1),9)'LAP2',(DLAP(2,I),I=1,LC)
          WRITE(KW(1),9)'PPC1',(PPCF(1,I),I=1,LC)
          WRITE(KW(1),9)'PPC2',(PPCF(2,I),I=1,LC)
          WRITE(KW(1),9)'FRS1',(FRST(1,I),I=1,LC)
          WRITE(KW(1),9)'FRS2',(FRST(2,I),I=1,LC)
          WRITE(KW(1),9)'WAC1',(WAC2(1,I),I=1,LC)
          WRITE(KW(1),9)'WAC2',(WAC2(2,I),I=1,LC)
          WRITE(KW(1),9)'BLG1',(BLG(1,I),I=1,LC)
          WRITE(KW(1),9)'BLG2',(BLG(2,I),I=1,LC)
      END IF    
      RETURN
    1 FORMAT(10X,100(6X,A4))
    6 FORMAT(T7,A4,100F10.1)
    7 FORMAT(T7,A4,100F10.2)
    8 FORMAT(T7,A4,100F10.3)
    9 FORMAT(T7,A4,100F10.4)
   10 FORMAT(T7,A4,100I10)
   11 FORMAT(//1X,'____________________CROP PARAMETERS__________________&
      __'/)
      END