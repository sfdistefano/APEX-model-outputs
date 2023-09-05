      SUBROUTINE ASORT3(D,NX,M1)
!     SHELL SORT
      DIMENSION D(M1),NX(M1)
      M=1
      DO
          IF((2**M)>M1)EXIT
          M=M+1
      END DO
      M=M-1
      M=2**M
      DO
          K=M1-M
          DO I=1,K
              DO J=I,1,-M
                  IF(D(NX(J+M))>=D(NX(J)))EXIT
                  N1=NX(J)
                  NX(J)=NX(J+M)
                  NX(J+M)=N1
              END DO
          END DO
          M=M/2
          IF(M<=0)EXIT
      END DO
      RETURN
      END