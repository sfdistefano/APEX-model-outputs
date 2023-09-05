C========== DEVELOPED BY LIWANG MA === FOR Single PRECISIOIN ===
C ==============================================================
      SUBROUTINE REALMATCH (NLAYRI,DLAYR,VI,NLAYRO,DS)

      IMPLICIT Real (A-H,O-Z)

      INTEGER  J,K,L,NLAYRI,NLAYRO,mxnod
      parameter (mxnod=300)

      DIMENSION    VI(mxnod),DS(nlayro),VS(mxnod),DLAYR(nlayri)

      IF (NLAYRI .LT. 1) RETURN
      IF (NLAYRO .LT. 1) RETURN

      DO L = 1, mxnod
         VS(L) = 0.0
      END DO
      TSUMI = 0.0d0

      K   = 1
      ZIL = 0.0
      ZOL = 0.0
         SUMZ = 0.0
         SUMV = 0.0

      DO 200 L = 1, NLAYRO

201         ZT   = MAX (ZOL,ZIL)
            ZB   = MIN (DS(L),DLAYR(K))
            SUMZ = SUMZ + (ZB - ZT)
            SUMV = SUMV + VI(K)*(ZB - ZT)
            IF (DS(L) .LT. DLAYR(K)) THEN
                   VS(L) = VI(K)
                IF (SUMZ .GT. 0.0) THEN
                   VS(L) = SUMV/SUMZ
                ENDIF
            TSUMI=TSUMI + sumv
            ZOL = DS(L)
            SUMZ = 0.0d0
          SUMV = 0.0d0
            GOTO 200
            ENDIF

            IF (DS(L).EQ.DLAYR(K)) THEN
                  VS(L) = VI(K)
                IF (SUMZ .GT. 0.0) THEN
                   VS(L) = SUMV/SUMZ
                ENDIF
            TSUMI=TSUMI + sumv
            SUMZ = 0.0
            SUMV = 0.0
            ZOL = DS(L)
            ZIL = DLAYR (K)
            K = K+1
         if (k.gt.nlayri) then
           goto 450
         endif
            GOTO 200
            ENDIF
         if (k.eq.nlayri) then
           goto 450
         else
            ZIL = DLAYR (K)
            K = K+1
            goto 201
         endif

200    CONTINUE

450      DO 400 J = 1, NLAYRO
         VI(J) = VS(J)
400   CONTINUE
      IF ((NLAYRO.GT.0).AND.(NLAYRO.LT.NLAYRI)) THEN
      DO 500 J = 1+NLAYRO, NLAYRI
         VI(J) = 0.0D0
500   CONTINUE
      ENDIF
      TSUMO = VI(1) * DS(1)
      DO L=2,NLAYRO
      TSUMO=TSUMO + VI(L) * (DS(L)-DS(L-1))
      END DO
      if (tsumo.eq.0.0) return

      IF (ABS(TSUMI/TSUMO-1.0).GT.0.01) THEN
      STOP 'YOU HAVE A PROBLEM WITH TRANSFORMATION BETWEEN GRIDS in REALMATCH'
      ENDIF

      RETURN

      END
