      SUBROUTINE APRNTD
      ! APEX1905                                                                       
      ! THIS SUBPROGRAM PRINTS DAILY IN THE .OUT FILE
      USE PARM         
      WRITE(KW(1),1090)ISA,NBSA(ISA),IYR,MO,KDA,(HED(KD(K)),VAR(KD(K),&
      ISA),K=1,NKD)
      WRITE(KW(1),1201)(HEDS(KS(K)),VARS(KS(K)),K=1,NKS)
      WRITE(KW(1),1201)((HEDC(K),VARC(K,LY(IRO(ISA),J,ISA),ISA),K=1,21),&
      J=1,NCP(IRO(ISA),ISA))
!     WRITE(KW(1),1200)(LID(K,ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'Z',(Z(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'ST',(SWST(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'S15',(S15(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'FC',(FC(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'PO',(PO(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'BD',(BD(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'WT',(WT(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'SAN',(SAN(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'SIL',(SIL(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'ROK',(ROK(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'HK',(HK(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'WNO3',(WNO3(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'WNH3',(WNH3(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'WPML',(WPML(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'WPMA',(WPMA(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'WPMS',(WPMS(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'RWT',(RWT(LID(K,ISA),JJK,ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'SSO3',(SSO3(LID(K,ISA)),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'U',(U(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'SEV',(SEV(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'O',(O(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'SSF',(SSF(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'UN',(UN(LID(K,ISA)),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'UP',(UP(LID(K,ISA)),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'WPO',(WPO(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'FOP',(FOP(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1040)'RSD',(RSD(LID(K,ISA),ISA),K=1,NBSL(ISA))
!     WRITE(KW(1),1110)ISA,NBSA(ISA),IYR,MO,KDA,(STMP(LID(K,ISA),ISA),
!     K=1,NBSL(ISA)),DD,DST0(ISA)
      RETURN
 1090 FORMAT(/1X,2I8,1X,I4,2I2,2X,5(1X,A4,F8.2)/(19X,5(1X,A4,F8.2)))
 1201 FORMAT((19X,5(1X,A4,F8.2)))      
      END