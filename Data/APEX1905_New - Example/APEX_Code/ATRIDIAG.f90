      SUBROUTINE ATRIDIAG(B,D,A,C,N)
      ! APEX1905
      ! SUBPROGRAM ATRIDIAG USES THE THOMAS ALHORITHM
      ! D()...SOLUTION RETURNED AS C()
      ! B()...BELOW DIAGONAL ELEMENTS
      ! D()...DIAGONAL ELEMENTS
      ! A()...ABOVE DIAGONAL ELEMENTS
      ! C()...RIGHT HAND SIDE
      DIMENSION A(N),B(N),C(N),D(N)
      ! FORWARD ELIMINATION
      DO I=2,N
          R=B(I)/D(I-1)
          D(I)=D(I)-R*A(I-1)
          C(I)=C(I)-R*C(I-1)
      END DO
      ! BACK SUBSTITUTION
      C(N)=C(N)/D(N)
      DO I=2,N
          J=N-I+1
          C(J)=(C(J)-A(J)*C(J+1))/D(J)
      END DO           
      RETURN
      END