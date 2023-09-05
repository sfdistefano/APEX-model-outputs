      SUBROUTINE HPERC(DZ,SATX)
!     APEX1905
!     THIS SUBPROGRAM COMPUTES PERCOLATION AND LATERAL SUBSURFACE FLOW
!     FROM A SOIL LAYER WHEN FIELD CAPACITY IS EXCEEDED.
      USE PARM
      SEP=0.
      SST(IDO)=0.
      SUP=SWST(ISL,ISA)-FC(ISL,ISA)
      IF(SUP>0.)THEN
          X1=24./(PO(ISL,ISA)-FC(ISL,ISA))
          Z1=X1*SATX
          X2=X1*HCL(ISL,ISA)
          XZ=X2+Z1
          IF(XZ<20.)THEN
              X3=SUP*(1.-EXP(-XZ))
          ELSE
              X3=SUP
          END IF
          SEP=X3/(1.+X2/Z1)
          X3=X3-SEP
          IF(ISL/=IDR(ISA))THEN      
              Z2=PRMT(90)*DZ*X2/SPLG(ISA)
              X3=X3*(1.-EXP(-Z2))
              X1=MIN(X3,.001*X3*SPLG(ISA)/RCHL(ISA))
              SST(IDO)=X1
              QRF(IDO)=X3-X1
          ELSE
              QRF(IDO)=X3
              SST(IDO)=0.
          END IF
          IF(RSAE(ISA)<1.E-10)RETURN
      END IF
      QRF(IDO)=0.
      RETURN
      END


      
      