      module INI
      USE DIA_SM
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)   
      contains

      SUBROUTINE initial_value_B(NX,NY,NZ,DT,RHO,N,X,U,V,W,OUTPUT_B)
      ! Initial value input
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      REAL,PARAMETER::PAI=3.1415926
      INTEGER::NX,NY,NZ,N,I,POSITION
      REAL::DT,RHO
      INTEGER,DIMENSION(N,3)::X
      REAL(KIND=RK),DIMENSION(N)::U,V,W
      REAL(KIND=RK) OUTPUT_B(3*NX*NY*NZ)

      position(X0,X1,X2,X3)=X1+(X2-1)*NX+(X3-1)*NX*NY+(X0-1)*NX*NY*NZ
      
      DO I=1,N
          IF(X(I,1)<1.OR.X(I,1)>NX)CYCLE
          IF(X(I,2)<1.OR.X(I,2)>NY)CYCLE
          IF(X(I,3)<1.OR.X(I,3)>NZ)CYCLE
          OUTPUT_B(POSITION(1,X(I,1),X(I,2),X(I,3)))=OUTPUT_B(POSITION
     +     (1,X(I,1),X(I,2),X(I,3)))+U(I)
          OUTPUT_B(POSITION(2,X(I,1),X(I,2),X(I,3)))=OUTPUT_B(POSITION
     +     (2,X(I,1),X(I,2),X(I,3)))+V(I)
          OUTPUT_B(POSITION(3,X(I,1),X(I,2),X(I,3)))=OUTPUT_B(POSITION
     +     (3,X(I,1),X(I,2),X(I,3)))+W(I)
      END DO

      END
      
      
      SUBROUTINE CHECK_M(N,MN)
      ! M-file validation
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      INTEGER::N,I,FSIZE,ESIZE
      LOGICAL::IS
      INTEGER,DIMENSION(:),allocatable::MN
      CHARACTER(LEN=30)::F
      CHARACTER(LEN=10)::C
      
      ESIZE=0
      DO I=1,N
          WRITE(C,'(I6)')I
          F='.\M\M'//TRIM(ADJUSTL(C))//'.DAT'
      
          INQUIRE(FILE=F, EXIST=IS, SIZE=FSIZE)
          IF(.NOT.IS)THEN
              MN=[MN,I]
          ELSE
              IF(ESIZE==0)ESIZE=FSIZE
              IF(ESIZE/=FSIZE)MN=[MN,I]
          ENDIF
      ENDDO
      
      
      
      END
      
      
      end module INI