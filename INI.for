      module INI
      USE DIA_SM
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)   
      contains

      
      
      SUBROUTINE initial_value_A(NX,NY,NZ,N,X,A,AI,AJ)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::NX,NY,NZ,N,I,J,POSITION
      INTEGER,DIMENSION(N,3)::X
      REAL(KIND=RK),ALLOCATABLE::A(:)
      INTEGER,ALLOCATABLE::AI(:),AJ(:)

      position(X0,X1,X2,X3)=X1+(X2-1)*NX+(X3-1)*NX*NY+(X0-1)*NX*NY*NZ
      
      DO I=1,N 
          DO J=1,SIZE(A)
              IF(AI(J).EQ.POSITION(1,X(I,1),X(I,2),X(I,3))) A(J)=0
              IF(AI(J).EQ.POSITION(2,X(I,1),X(I,2),X(I,3))) A(J)=0
              IF(AI(J).EQ.POSITION(3,X(I,1),X(I,2),X(I,3))) A(J)=0
          ENDDO
          A=[A,real(1,kind=rk)]
          A=[A,real(1,kind=rk)]
          A=[A,real(1,kind=rk)]
          AI=[AI,POSITION(1,X(I,1),X(I,2),X(I,3))]
          AJ=[AJ,POSITION(1,X(I,1),X(I,2),X(I,3))]
          AI=[AI,POSITION(2,X(I,1),X(I,2),X(I,3))]
          AJ=[AJ,POSITION(2,X(I,1),X(I,2),X(I,3))]
          AI=[AI,POSITION(3,X(I,1),X(I,2),X(I,3))]
          AJ=[AJ,POSITION(3,X(I,1),X(I,2),X(I,3))]
      END DO
      CALL SORT(A,AI,AJ)
      END
      
      
      
      SUBROUTINE initial_value(NX,NY,NZ,N,X,T_START,T_END,UT,VT,WT,
     +                        A,AI,AJ,OUTPUT_X)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::NX,NY,NZ,N,T_START,T_END,I,J,K,POSITION
      INTEGER,DIMENSION(N,3)::X
      REAL(KIND=RK),DIMENSION(N,T_START:T_END)::UT,VT,WT
      REAL(KIND=RK),ALLOCATABLE::A(:)
      INTEGER,ALLOCATABLE::AI(:),AJ(:)
      REAL(KIND=RK) OUTPUT_X(T_START:T_END,3*NX*NY*NZ)

      position(X0,X1,X2,X3)=X1+(X2-1)*NX+(X3-1)*NX*NY+(X0-1)*NX*NY*NZ
      
      DO I=1,N 
          DO J=1,SIZE(A)
              IF(AI(J).EQ.POSITION(1,X(I,1),X(I,2),X(I,3))) A(J)=0
              IF(AI(J).EQ.POSITION(2,X(I,1),X(I,2),X(I,3))) A(J)=0
              IF(AI(J).EQ.POSITION(3,X(I,1),X(I,2),X(I,3))) A(J)=0
          ENDDO
          A=[A,real(1,kind=rk)]
          A=[A,real(1,kind=rk)]
          A=[A,real(1,kind=rk)]
          AI=[AI,POSITION(1,X(I,1),X(I,2),X(I,3))]
          AJ=[AJ,POSITION(1,X(I,1),X(I,2),X(I,3))]
          AI=[AI,POSITION(2,X(I,1),X(I,2),X(I,3))]
          AJ=[AJ,POSITION(2,X(I,1),X(I,2),X(I,3))]
          AI=[AI,POSITION(3,X(I,1),X(I,2),X(I,3))]
          AJ=[AJ,POSITION(3,X(I,1),X(I,2),X(I,3))]
          DO T=T_START,T_END
              OUTPUT_X(T,POSITION(1,X(I,1),X(I,2),X(I,3)))=UT(I,T)
              OUTPUT_X(T,POSITION(2,X(I,1),X(I,2),X(I,3)))=VT(I,T)
              OUTPUT_X(T,POSITION(3,X(I,1),X(I,2),X(I,3)))=WT(I,T)
          END DO
      END DO
      END
C 已知点值控制
      SUBROUTINE initial_value_B(NX,NY,NZ,DT,RHO,N,X,U,V,W,OUTPUT_B)
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
      
      
      
      
      
      
      
      
      end module INI