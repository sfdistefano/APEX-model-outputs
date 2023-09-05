      SUBROUTINE RTMCCEL(QFLO,AFLO,XMC,CLTY)
!     THIS SUB PROGRAM COMPUTES THE EXPONENT M OF FLOW AREA TO COMPUTE 
!     FLOW RATE
      USE PARM                                                     
      IF(QFLO<SCFS(1,IDN2))THEN
          XMC=XMS(1,IDN2)
      ELSE
          DO J=2,NPRC
              J1=J-1
              IF(QFLO<SCFS(J,IDN2))EXIT
          END DO
          J=MIN(J,NPRC)
          RTO=(QFLO-SCFS(J1,IDN2))/(SCFS(J,IDN2)-SCFS(J1,IDN2))
          XMC=RTO*(XMS(J,IDN2)-XMS(J1,IDN2))+XMS(J1,IDN2)
      END IF
      CLTY=XMC*QFLO/AFLO
      RETURN
      END