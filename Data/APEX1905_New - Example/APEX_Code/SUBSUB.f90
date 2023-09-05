      SUBROUTINE SUBSUB(NN,GZSS,GZSL,GZSD,SPSL,SPSD,GZSLT,GZSDT,SPSLT,SPSDT,SPGZ)
!     APEX1905
!     THIS SUBPROGRAM SIMULATES ANIMAL GRAZING USING CONCEPTS FROM PHY-GROW
!     GRAZING PREFERENCE COMPONENT.
      USE PARM
      DIMENSION GZSD(3,KOMX(IOW),MNC,MSA),GZSL(3,KOMX(IOW),MNC,MSA),SPSL(3,KOMX(IOW),MNC,MSA), &
                SPSD(3,KOMX(IOW),MNC,MSA),GZSLT(3,KOMX(IOW),MNC),GZSDT(3,KOMX(IOW),MNC), &
                GZSS(MNC,MSA),SPSLT(3,KOMX(IOW),MNC),SPSDT(3,KOMX(IOW),MNC),SPGZ(KOMX(IOW)), &
                WSAV(MNC,MSA)
      GZSS=1.0
      GZSLT=0.0
      GZSDT=0.0
      SPSLT=0.0
      SPSDT=0.0
      WSAV=0.0
      WSAT=0.0
      SPGZ=0.0
      DO KKP=1,3 ! PREFERENCE LOOP														
          DO KHD=1,NHRD(IOW) ! HERD LOOP
              IHD=KOW(KHD,IOW)
              DO I=1,NN ! FORAGE LOOP													
				  DO J=1,MSA
                     DO K=1,MSA
                        IF (NBSA(J)==NPAS(K)) THEN
                        KKF=JE(I,K)
                        IF(KKF==MNC) CYCLE
                        X2=STL(KKF,K)
                        IF(IDC(KKF)==NDC(8))X2=X2*FTO(KKF)*SLAI(KKF,K)/XLAI(KKF,K)
						    WSAV(KKF,K)=WSA(K)*(X2+STD(KKF,K))
						    WSAT(KKF,J)=WSAT(KKF,J)+WSAV(KKF,K)
					    END IF
					 END DO
!                        IF(KKF==MNC) CYCLE
!                        IF (KKP==1.AND.KHD==1) GZSS(KKF,J)=WSAV(KKF,J)/WSAT(KKF,J)
                  END DO
				  DO J=1,MSA
                        KKF=JE(I,J)
                        GZSS(KKF,J)=0.25   !WSAV(KKF,J)/WSAT(KKF,J)
                     DO K=1,MSA
                     IF (NPAS(K)==NBSA(J)) THEN
                        IF(KKF==MNC) CYCLE
				        GZSLT(KKP,IHD,KKF)=GZSLT(KKP,IHD,KKF)+GZSL(KKP,IHD,KKF,J)*GZSS(KKF,J)
				        GZSDT(KKP,IHD,KKF)=GZSDT(KKP,IHD,KKF)+GZSD(KKP,IHD,KKF,J)*GZSS(KKF,J)
				        SPSLT(KKP,IHD,KKF)=SPSLT(KKP,IHD,KKF)+SPSL(KKP,IHD,KKF,J)*GZSS(KKF,J)
				        SPSDT(KKP,IHD,KKF)=SPSDT(KKP,IHD,KKF)+SPSD(KKP,IHD,KKF,J)*GZSS(KKF,J)
!                    ELSE
!                        GZSS(KKF,J)=0.25   !WSAV(KKF,J)/WSAT(KKF,J)
!				        GZSLT(KKP,IHD,KKF)=GZSLT(KKP,IHD,KKF)+GZSL(KKP,IHD,KKF,K)*GZSS(KKF,J)
!				        GZSDT(KKP,IHD,KKF)=GZSDT(KKP,IHD,KKF)+GZSD(KKP,IHD,KKF,K)*GZSS(KKF,J)
!				        SPSLT(KKP,IHD,KKF)=SPSLT(KKP,IHD,KKF)+SPSL(KKP,IHD,KKF,K)*GZSS(KKF,J)
!				        SPSDT(KKP,IHD,KKF)=SPSDT(KKP,IHD,KKF)+SPSD(KKP,IHD,KKF,K)*GZSS(KKF,J)
                        SPGZ(IHD)=SPGZ(IHD)+SPSLT(KKP,IHD,KKF)+SPSDT(KKP,IHD,KKF) 
                    ENDIF
                  END DO
                  END DO
              END DO    
          END DO
      END DO
	  RETURN
END