      SUBROUTINE GASTRANS2(CONC,DPRM,CUP,DSURF,FLUXOUT)
      ! APEX1905
      ! SUBPROGRAM SOLVES THE GAS TRANSPORT EQUATION
      USE PARM
      DIMENSION A(MSC),B(MSC),C(MSC),D(MSC),CONC(MSC),CONCTMP(MSC),&
      DPRM(MSC) ! arrays for storing tridiagonal matrix and RHS
      LOGICAL::CONVG          ! convergence flag
      DATA ALX/.5/,DTMIN/.01/ ! ALX is the explicit fraction, 1-ALX the implicit 
      C=0.
      CONVG=.FALSE.
      DTSUB=DTG 
      DO WHILE(.NOT.CONVG)
          !! set up the integration
          CONCTMP=CONC        ! set concentration to initial values
          T=0.                ! set time to top of the hour
          FLUXOUT=0.          ! zero the hourly flux
          CONVG=.TRUE.        ! will be unset when and if we observe a failure
          DO WHILE(T<DTG)
              !! These must be initialized here because we may take a
              !! short step at the end of the integration interval
              R=DTSUB/(DZDN*DZDN)
              R1=ALX*R
              R2=(1.-ALX)*R
              ! UPPER BOUNDARY CONDITION
              A(1)=-R2*DPRM(1)
              B(1)=0.0
              C(1)=AFP(1,ISA)*CONCTMP(1)+R1*(DPRM(1)*(CONCTMP(2)-CONCTMP(1))-&
              DSURF*(CONCTMP(1)-CUP))+R2*DSURF*CUP
              D(1)=AFP(1,ISA)+R2*(DSURF+DPRM(1))
              ! Calculate surface flux using the old concentration values
              F1=-DSURF*(CONCTMP(1)-CUP)/DZDN
              !      print *,'start conctmp(1), cup, F1,alx : ', conctmp(1), cup, F1, alx !XXX
              ! MAIN COMPUTATIONS
              DO ID=2,IUN
                  A(ID)=-R2*DPRM(ID)
                  B(ID)=-R2*DPRM(ID-1)
                  C(ID)=AFP(ID,ISA)*CONCTMP(ID)+R1*(DPRM(ID)*(CONCTMP(ID+1)-&
                  CONCTMP(ID))-DPRM(ID-1)*(CONCTMP(ID)-CONCTMP(ID-1)))
                  D(ID)=AFP(ID,ISA)+R2*(DPRM(ID-1)+DPRM(ID))
              END DO
              ! LOWER BOUNDARY CONDITION
              A(NBCL)=0.0
              B(NBCL)=-R2*DPRM(NBCL-1)
              C(NBCL)=AFP(NBCL,ISA)*CONCTMP(NBCL)-R1*DPRM(NBCL-1)*(CONCTMP&
              (NBCL)-CONCTMP(NBCL-1))
              D(NBCL)=AFP(NBCL,ISA)+R2*DPRM(NBCL-1)
              ! SOLVE TRIADIAGONAL SYSTEM
              CALL ATRIDIAG(B,D,A,C,NBCL)
              !DO ID=1,NBCL
                  !IF(ABS(C(ID)-CONCTMP(ID))/ABS(CONCTMP(ID)+1.E-4)>.5.AND.DTSUB>DTMIN)THEN
              IF(ANY(ABS(C-CONCTMP)/(ABS(CONCTMP)+1.0E-4)>0.5).AND.DTSUB>DTMIN)THEN
                  !! excessive change.  reduce time step and try again
                  !! (but if the timestep is already too small, then we
                  !! just live with it)
                  DTSUB=DTSUB*.5
                  CONVG=.FALSE.
                  EXIT
              ELSE 
                  CONCTMP=MAX(C,1.E-10)
              END IF
              !END DO
              !IF(CONVG==.FALSE.)EXIT
              ! Calculate surface flux using the new concentration values
              F2=-DSURF*(C(1)-CUP)/DZDN ! use C instead of conctmp, so we don't see the effect of MAX above
              ! Combine the fluxes in the same ratio used in the Crank-Nicholson scheme
              FLUXOUT=FLUXOUT+(ALX*F1+(1.-ALX)*F2)*DTSUB
              ! update the time
              T=T+DTSUB
              IF(DTG-T<DTSUB)THEN
                  DTSUB=DTG-T    ! ensure that we don't overshoot our intended finish time
              END IF 
          END DO                 ! while(t < dtg
      END DO                     ! while(.not. convg)
      CONC=CONCTMP
      RETURN 
      END
