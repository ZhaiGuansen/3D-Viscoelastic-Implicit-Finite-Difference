      module Matrix_inverse
      
      
      !use INI
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
     
      contains
      
      SUBROUTINE LUM(N,A)
      INTEGER::N,I,J,K
      integer, parameter :: rk = kind ( 1.0D+00 )
      REAL(KIND=RK)::S
      REAL(KIND=RK),allocatable,DIMENSION(:,:)::A,IA
      
      !LU·Ö½â
      DO I=1,N
          DO J=1,N
              S=A(I,J)
              DO K=1,MIN(I,J)-1
                  S=S-A(I,K)*A(K,J)
              ENDDO
              IF(I<=J)THEN
                  A(I,J)=S
              ELSE
                  A(I,J)=S/A(J,J)
              ENDIF
          ENDDO
      ENDDO
      ALLOCATE(IA(N,N))
      !ÇóÄæ
      DO I=1,N
          !L
          DO J=1,I-1
              S=-A(I,J)
              DO K=J+1,I-1
                  S=S-A(I,K)*IA(K,J)
              ENDDO
              IA(I,J)=S
          ENDDO
          !U
          DO J=I,N
              S=0
              IF(I==J)S=1
              DO K=I,J-1
                  S=S-IA(I,K)*A(K,J)
              ENDDO
              IA(I,J)=S/A(J,J)
          ENDDO
          
      ENDDO
      
      DO I=1,N
          DO J=1,N
              A(I,J)=0
              DO K=MAX(I,J),N
                  IF(K==J)THEN
                      A(I,J)=IA(I,J)
                  ELSE
                      A(I,J)=A(I,J)+IA(I,K)*IA(K,J)
                  ENDIF
              ENDDO
          ENDDO
      ENDDO
      
      END
      
      
      
      end module Matrix_inverse