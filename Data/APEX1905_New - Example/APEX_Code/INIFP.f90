      SUBROUTINE INIFP(I3,II,JJ,JRT)
!     APEX1905
!     THIS SUBPROGRAM ALLOWS INPUT TO OPERATION SCHEDULE FOR IRRIGATION,
!     FERTILIZER, OR PESTICIDE APPLICATION
      USE PARM
      J1=0
      SELECT CASE(I3)
          CASE(7)
              KP=KP+1
              NPC(KP)=JJ
              KPC(KP)=JX(2)
              JPC(KP)=JX(3)
              PSTR(II,KP,ISA)=OPV(1)
              PSTE(II,KP,ISA)=OPV(2)
              CALL PSTTBL
              LPC(II,KP,ISA)=KDP1(JX(7))
          CASE(8)
              VIRR(II,JJ,ISA)=OPV(1)
              IF(ABS(OPV(3))>1.E-5)BIR(ISA)=OPV(3)
              IF(OPV(4)>0.)EFI(ISA)=OPV(4)
              IF(OPV(5)>0.)FIRG(ISA)=OPV(5)
              IF(OPV(6)>0.)HMIN_IRR(II,JJ,ISA)=OPV(6) !PADDY minimum ponding depth that triggers irrigation JAEHAK JEONG 2016
              KI=KI+1
              NIR(KI)=JX(2)
              IIR(KI)=JX(3)
              KIR(KI)=JJ
          CASE(9)
              WFA(II,JJ,ISA)=OPV(1)
              CALL NFTBL(L)
              LFT(II,JJ,ISA)=KDF1(JX(7))      
              KF=KF+1
          CASE(15) !SET PADDY OUTLET WEIR HEIGHT 
              HWEIR(II,JJ,ISA)=OPV(1) !PADDY MODEL JAEHAK JEONG 2014
              LWEIR(II,JJ,ISA)=OPV(2) !PADDY WEIR WIDTH, M
          CASE(19)
              JGRZ(II,JJ,ISA)=JX(8)
              RSTK(II,JJ,ISA)=OPV(1)
              DO KHD=1,NHRD(IOW)
              IF(OPV(4)>0)MAXD(II,JJ,ISA)=OPV(4)	!read maximum gazing days limit in OPC
              IF(OPV(5)>0)MIND(II,JJ,ISA)=OPV(5)	!read minimum gazing days limit in OPC
              if(OPV(6)>0)GZLM(IHD,ISA)=OPV(6)		!read biomass limit in OPC
                  IF(KOW(KHD,IOW)<=0)EXIT
                  IHD=KOW(KHD,IOW)
                  IF(JX(8)==IHDM(IHD,ISA))THEN
                      IF(OPV(2)>0.)GZWX(II,JJ,ISA)=OPV(2)
                      IGZX(IHD,IOW)=ISA
                      !IGZR(IHD,ISA)=1
                  END IF    
              END DO
          CASE(20)
              JGRZ(II,JJ,ISA)=JX(8)
              RSTK(II,JJ,ISA)=0.
              DO KHD=1,NHRD(IOW)
                  IF(KOW(KHD,IOW)<=0)EXIT
                  IHD=KOW(KHD,IOW)
                  IF(JX(8)==IHDM(IHD,ISA))THEN
                      IGZR(IHD,ISA)=0
                  END IF
              END DO    
          CASE(27)
              WMUCH(II,JJ,ISA)=OPV(1)
          CASE(28)
              X1=WSA(ISA)/RSTK(II,JJ,ISA)
              X1=X1+OPV(1)
              IF(X1<1.E-5)THEN
                  GZWT(IHD,IOW)=0.
                  RSTK(II,JJ,ISA)=0.
                  GZNB(IHD,ISA)=0.
              ELSE
                  RSTK(II,JJ,ISA)=X1/WSA(ISA)
              END IF
          CASE DEFAULT
              J1=1
      END SELECT
      IF(J1==0)THEN
          TIR(II,JJ,ISA)=BIR(ISA)
          CND(II,JJ,ISA)=CN2(ISA)
          QIR(II,JJ,ISA)=EFI(ISA)
          FIRX(II,JJ,ISA)=FIRG(ISA)
          JRT=1
      ELSE    
          IF(OPV(2)<0.)THEN
              CN2(ISA)=-OPV(2)      
          ELSE
              IF(OPV(2)>0.)THEN
	              LUN(ISA)=OPV(2)
	              LUN(ISA)=LUN(ISA)+LUNS(ISA)
                  CALL HSGCN
              END IF
          END IF
          CND(II,JJ,ISA)=CN2(ISA)
          IF(ABS(OPV(3))>1.E-5)BIR(ISA)=OPV(3)
          TIR(II,JJ,ISA)=BIR(ISA)
          IF(OPV(4)>0.)EFI(ISA)=OPV(4)
          QIR(II,JJ,ISA)=EFI(ISA)
          FIRX(II,JJ,ISA)=FIRG(ISA)
          JRT=0
      END IF    
      RETURN
      END