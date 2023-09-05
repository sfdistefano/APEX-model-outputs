      SUBROUTINE HTR55
!     APEX1905
!     THIS SUBPROGRAM ESTIMATES PEAK RUNOFF RATES USING THE SCS TR55
!     EXTENDED METHOD.
      USE PARM
      DIMENSION PIAF(18)
      DATA PIAF/.0,.1,.2,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1./
      IF(QRB>.35)THEN
          XOX=(QRB-.35)/.05+5.
      ELSE
          XOX=QRB/.1+1.
      END IF
      INT=XOX
      INT1=INT+1
      RTO=(QRB-PIAF(INT))/(PIAF(INT1)-PIAF(INT))
      X1=LOG(TC(IDO))
      Y1=HQP(X1,CQRB,ITYP,INT)
      IF(INT>=17)THEN
          Y=EXP(Y1)*(1.0-RTO)
      ELSE
          Y2=HQP(X1,CQRB,ITYP,INT1)
          Y=Y1+(Y2-Y1)*RTO
          Y=EXP(Y)
      END IF
      RQRB(IDO)=Y*RFV(IRF(ISA))
      RETURN
      END