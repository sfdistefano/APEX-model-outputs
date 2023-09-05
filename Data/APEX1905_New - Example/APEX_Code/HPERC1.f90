      SUBROUTINE HPERC1(DZ,SATX)
!     APEX1905
!     THIS SUBPROGRAM COMPUTES PERCOLATION AND LATERAL SUBSURFACE FLOW
!     FROM A SOIL LAYER BY PASSING 4 mm SLUGS THROUGH AS A FUNCTION OF 
!     HYDRAULIC CONDUCTIVITY.
      USE PARM
      SEP=0.
      SSX=0.
      SST(IDO)=0.
      QRF(IDO)=0.
      ICW=0
      AVW=SWST(ISL,ISA)-FC(ISL,ISA)
      IF(AVW>.01)THEN
          POFC=PO(ISL,ISA)-FC(ISL,ISA)
          X1=24./POFC
          DO WHILE(AVW>.01)
              X5=MIN(AVW/POFC,1.)
              X4=MAX(1.E-5,X5**PRMT(82))
              H=SATX*X4
              X2=X1*HCL(ISL,ISA)*X4
              ZZ=X1*H
              XZ=X2+ZZ
              XX=MIN(4.,AVW)
              IF(XZ<20.)THEN
                  X3=XX*(1.-EXP(-XZ))
              ELSE
                  X3=XX
              END IF
              X6=X3/(1.+X2/ZZ)
              SEP=SEP+X6
              X7=X3-X6
              SSX=SSX+X7
              AVW=AVW-4.
              ICW=ICW+1
              IF(ICW>50)EXIT
          END DO    
          IF(ISL/=IDR(ISA))THEN      
              Z2=PRMT(90)*DZ*X2/SPLG(ISA)
              SSX=SSX*(1.-EXP(-Z2))
              X1=MIN(SSX,.001*SSX*SPLG(ISA)/RCHL(ISA))
              SST(IDO)=X1
              QRF(IDO)=SSX-X1
          ELSE
              QRF(IDO)=SSX
              SST(IDO)=0.
          END IF 
          IF(RSAE(ISA)>0.)QRF(IDO)=0.
      END IF
      RETURN
      END
