      module matrix_construction
      use DIA_SM
      use INI
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)  
      
      REAL::PI=3.141592653
      contains
      
      SUBROUTINE RICKER(NPTS,FP,DT,AMPL,SR)
      integer,parameter :: rk = kind ( 1.0D+00 )
      INTEGER NPTS
      REAL(KIND=RK)::SR(NPTS)
      REAL::DT,AMPL,FP

 !
 !	Subroutine RICKER is designed to Sample time, calculate the zero-phase 
 !     Ricker wavelet and Fourier transform it from domain to frequency domain.
 !

    DO I=1,NPTS
        SR(I)=0.
    ENDDO
    DO I=2,NPTS/2+1
        T=PI*FP*DT*REAL(I-1)
        SR(I)=AMPL*(1.0-2.0*T*T)*EXP(-T*T)
        SR(NPTS+2-I)=SR(I)
    ENDDO
    SR(1)=AMPL
!    CALL FORK(NPTS,SW,-1)  
      END

      SUBROUTINE MAKE1(LX,LY,LZ,VP,VS,DX,DY,DZ,DT,A1N,A1M1,A1M2,A1S1,&
   A1S2,A1S3,A1,B1N,B1M1,B1M2,B1S1,B1S2,B1S3,B1,C1N,C1M1,C1M2, C1S1,C1S2,C1S3,C1)
      ! Diagonal block matrix (7-diagonal structure)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::I,J,K,LX,LY,LZ,LINE,N,ROW
      REAL::DX,DY,DZ,DT
      real:: vp(LX,LY,LZ),vs(LX,LY,LZ)
      REAL(KIND=RK)::C0=0,X2,Y2,Z2,P,S
      INTEGER,DIMENSION(2)::A1N,B1N,C1N
      INTEGER,ALLOCATABLE,DIMENSION(:)::A1M1,A1M2,A1S1,A1S2,A1S3,B1M1,B1M2,B1S1,B1S2,B1S3,C1M1,C1M2,C1S1,C1S2,C1S3
      REAL(KIND=RK),ALLOCATABLE::A1(:),B1(:),C1(:)
      
      A1N(1)=LX*LY*LZ
      B1N(1)=LX*LY*LZ
      C1N(1)=LX*LY*LZ
      A1N(2)=7
      B1N(2)=7
      C1N(2)=7
      A1M1=[A1N(1)-LX*LY,A1N(1)-LX,A1N(1)-1,A1N(1),A1N(1)+1,A1N(1)+LX,A1N(1)+LX*LY]
      B1M1=A1M1
      C1M1=A1M1
      A1M2=[1,1,1,1,1,1,1]
      B1M2=A1M2
      C1M2=A1M2
      A1S1=[1,1,1,1,1,1,1]
      B1S1=A1S1
      C1S1=A1S1
      A1S3=[A1N(1)-LX*LY,A1N(1)-LX,A1N(1)-1,A1N(1),A1N(1)-1,A1N(1)-LX,A1N(1)-LX*LY]
      B1S3=A1S3
      C1S3=A1S3
      A1S2=[A1S2,1]
      DO I=2,7
          A1S2=[A1S2,A1S2(I-1)+A1S3(I-1)]
      ENDDO
      B1S2=A1S2
      C1S2=A1S2
      
      
      x2=DT/(dx*dx)
      y2=DT/(dy*dy)
      z2=DT/(dz*dz)

      DO LINE=1,7
          DO N=1,A1N(1)-ABS(A1M1(LINE)-A1N(1))
              IF(LINE<4)THEN
                  ROW=A1N(1)-A1M1(LINE)+N
              ELSE
                  ROW=N
              ENDIF
              
              K=(ROW-1)/(LX*LY)
              I=ROW-K*LX*LY
              K=K+1
              J=(I-1)/LX
              I=I-J*LX
              J=J+1
      IF(I==1.OR.J==1.OR.K==1.OR.I==LX.OR.J==LY.OR.K==LZ)THEN
          A1=[A1,C0]
          B1=[B1,C0]
          C1=[C1,C0]
      ELSE
          P=VP(I,J,K)**2
          S=VS(I,J,K)**2
          SELECT CASE(LINE)
          CASE(1)
              A1=[A1,S*z2]
              B1=[B1,S*z2]
              C1=[C1,P*z2]
          CASE(2)
              A1=[A1,S*y2]
              B1=[B1,P*y2]
              C1=[C1,S*y2]
          CASE(3)
              A1=[A1,P*x2]
              B1=[B1,S*x2]
              C1=[C1,S*x2]
              
          CASE(4)
              A1=[A1,-2*(P*x2+S*y2+S*z2)]
              B1=[B1,-2*(P*y2+S*x2+S*z2)]
              C1=[C1,-2*(P*z2+S*y2+S*x2)]
          CASE(5)
              A1=[A1,P*x2]
              B1=[B1,S*x2]
              C1=[C1,S*x2]  

          CASE(6)
              A1=[A1,S*y2]
              B1=[B1,P*y2]
              C1=[C1,S*y2] 

          CASE(7)
              A1=[A1,S*z2]
              B1=[B1,S*z2]
              C1=[C1,P*z2]

          END SELECT
      ENDIF
          ENDDO
      ENDDO
      END
 
      
      SUBROUTINE MAKED(LX,LY,LZ,D1N,D1M1,D1M2,D1S1,D1S2,D1S3,D1,D2N,D2M1,D2M2,D2S1,D2S2,D2S3,D2,D3N,D3M1,D3M2,D3S1,D3S2,D3S3,D3)
      ! Off-diagonal block matrix (4-diagonal structure)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::I,J,K,LX,LY,LZ,LINE,ROW,N
      REAL(KIND=RK)::C0=0,C1=1
      INTEGER,DIMENSION(2)::D1N,D2N,D3N
      INTEGER,ALLOCATABLE,DIMENSION(:)::D1M1,D1M2,D1S1,D1S2,D1S3,D2M1,D2M2,D2S1,D2S2,D2S3,D3M1,D3M2,D3S1,D3S2,D3S3
      REAL(KIND=RK),ALLOCATABLE::D1(:),D2(:),D3(:)
      
      D1N(1)=LX*LY*LZ
      D2N(1)=LX*LY*LZ
      D3N(1)=LX*LY*LZ
      D1N(2)=4
      D2N(2)=4
      D3N(2)=4
      D1M1=[D1N(1)-LX-LX*LY,D1N(1)+LX-LX*LY,D1N(1)-LX+LX*LY,D1N(1)+LX+LX*LY]
      D2M1=[D1N(1)-1-LX*LY,D1N(1)+1-LX*LY,D1N(1)-1+LX*LY,D1N(1)+1+LX*LY]
      D3M1=[D1N(1)-1-LX,D1N(1)+1-LX,D1N(1)-1+LX,D1N(1)+1+LX]
      D1M2=[1,1,1,1]
      D2M2=[1,1,1,1]
      D3M2=[1,1,1,1]
      D1S1=[1,1,1,1]
      D2S1=[1,1,1,1]
      D3S1=[1,1,1,1]
      D1S3=[D1N(1)-LX-LX*LY,D1N(1)+LX-LX*LY,D1N(1)+LX-LX*LY,D1N(1)-LX-LX*LY]
      D2S3=[D1N(1)-1-LX*LY,D1N(1)+1-LX*LY,D1N(1)+1-LX*LY,D1N(1)-1-LX*LY]
      D3S3=[D1N(1)-1-LX,D1N(1)+1-LX,D1N(1)+1-LX,D1N(1)-1-LX]
      
      D1S2=[1]
      D2S2=[1]
      D3S2=[1]
      DO I=2,4
          D1S2=[D1S2,D1S2(I-1)+D1S3(I-1)]
          D2S2=[D2S2,D2S2(I-1)+D2S3(I-1)]
          D3S2=[D3S2,D3S2(I-1)+D3S3(I-1)]
      ENDDO

      
      DO LINE=1,4
          DO N=1,D1N(1)-ABS(D1M1(LINE)-D1N(1))
              IF(LINE<3)THEN
                  ROW=D1N(1)-D1M1(LINE)+N
              ELSE
                  ROW=N
              ENDIF

              K=INT((ROW-1)/(LX*LY))
              I=ROW-K*LX*LY
              K=K+1
              J=INT((I-1)/LX)
              I=I-J*LX
              J=J+1
      IF(I==1.OR.J==1.OR.K==1.OR.I==LX.OR.J==LY.OR.K==LZ)THEN
          D1=[D1,C0]
      ELSE
          SELECT CASE(LINE)
          CASE(1)
              D1=[D1,C1]
          CASE(2)
              D1=[D1,-C1]
          CASE(3)
              D1=[D1,-C1]
          CASE(4)
              D1=[D1,C1]
          END SELECT
      ENDIF             
          ENDDO
          
          
          DO N=1,D2N(1)-ABS(D2M1(LINE)-D2N(1))
              IF(LINE<3)THEN
                  ROW=D2N(1)-D2M1(LINE)+N
              ELSE
                  ROW=N
              ENDIF
              
              K=(ROW-1)/(LX*LY)
              I=ROW-K*LX*LY
              K=K+1
              J=(I-1)/LX
              I=I-J*LX
              J=J+1
      IF(I==1.OR.J==1.OR.K==1.OR.I==LX.OR.J==LY.OR.K==LZ)THEN
          D2=[D2,C0]
      ELSE
          SELECT CASE(LINE)
          CASE(1)
              D2=[D2,C1]
          CASE(2)
              D2=[D2,-C1]
          CASE(3)
              D2=[D2,-C1]
          CASE(4)
              D2=[D2,C1]
          END SELECT
      ENDIF            
          ENDDO
          
          
          DO N=1,D3N(1)-ABS(D3M1(LINE)-D3N(1))
              IF(LINE<3)THEN
                  ROW=D3N(1)-D3M1(LINE)+N
              ELSE
                  ROW=N
              ENDIF
              
              K=(ROW-1)/(LX*LY)
              I=ROW-K*LX*LY
              K=K+1
              J=(I-1)/LX
              I=I-J*LX
              J=J+1
      IF(I==1.OR.J==1.OR.K==1.OR.I==LX.OR.J==LY.OR.K==LZ)THEN
          D3=[D3,C0]
      ELSE
          SELECT CASE(LINE)
          CASE(1)
              D3=[D3,C1]
          CASE(2)
              D3=[D3,-C1]
          CASE(3)
              D3=[D3,-C1]
          CASE(4)
              D3=[D3,C1]
          END SELECT
      ENDIF             
          ENDDO      
      ENDDO
      END
    
          
      !(I+M2)U(N+1)=(M1+2I)U(N)+(M2-I)U(N-1)
      
      SUBROUTINE getM(LX,LY,LZ,VP,VS,NP,NS,dx,dy,dz,dt,M1N,M1M1,M1M2,M1S1,M1S2,M1S3,M1,M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2,NLABEL)
      ! Assemble and reorganize block matrices to construct M matrix
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::LX,LY,LZ,N,NLABEL
      real::dx,dy,dz,dt
      REAL(KIND=RK)::a,b,K1(LX*LY*LZ),K2(LX*LY*LZ)
      REAL,DIMENSION(LX,LY,LZ)::VP,VS,NP,NS
      INTEGER,DIMENSION(2)::A1N,B1N,C1N,A2N,B2N,C2N,D1N,D2N,D3N,M1N,M2N,K2D1N,K2D2N,K2D3N
      INTEGER,ALLOCATABLE,DIMENSION(:)::A1M1,A1M2,A1S1,A1S2,A1S3,B1M1,B1M2,B1S1,B1S2,B1S3,C1M1,C1M2,C1S1,C1S2,C1S3,&
          A2M1,A2M2,A2S1,A2S2,A2S3,B2M1,B2M2,B2S1,B2S2,B2S3, C2M1,C2M2,C2S1,C2S2,C2S3,D1M1,D1M2,D1S1,D1S2,D1S3,&
          D2M1,D2M2,D2S1,D2S2,D2S3,D3M1,D3M2,D3S1,D3S2,D3S3, M1M1,M1M2,M1S1,M1S2,M1S3,M2M1,M2M2,M2S1,M2S2,M2S3,&
          K2D1M1,K2D1M2,K2D1S1,K2D1S2,K2D1S3,K2D2M1,K2D2M2, K2D2S1,K2D2S2,K2D2S3,K2D3M1,K2D3M2,K2D3S1,K2D3S2,K2D3S3
      REAL(KIND=RK),ALLOCATABLE,DIMENSION(:)::A1,B1,C1,A2,B2,C2, D1,D2,D3,M1,M2,K2D1,K2D2,K2D3

      call MAKE1(LX,LY,LZ,VP,VS,DX,DY,DZ,DT*DT,A1N,A1M1,A1M2,A1S1, A1S2,A1S3,A1,B1N,B1M1,B1M2,B1S1,B1S2,B1S3,B1,&
          C1N,C1M1,C1M2, C1S1,C1S2,C1S3,C1)
      IF(NLABEL/=0)THEN
      call MAKE1(LX,LY,LZ,NP,NS,DX,DY,DZ,DT/2,A2N,A2M1,A2M2,A2S1, A2S2,A2S3,A2,&
          B2N,B2M1,B2M2,B2S1,B2S2,B2S3,B2,C2N,C2M1,C2M2,C2S1,C2S2,C2S3,C2)
      ENDIF
      call MAKED(LX,LY,LZ,D1N,D1M1,D1M2,D1S1,D1S2,D1S3,D1,D2N,D2M1,D2M2,D2S1,D2S2,D2S3,D2,D3N,D3M1,D3M2,D3S1,D3S2,D3S3,D3)
      
      a=dx/dz
      b=dy/dz
      N=LX*LY*LZ
      K1=reshape((vp**2-vs**2)*dt*dt/(4*dx*dy),[N])
      K2=reshape((np**2-ns**2)*dt/(8*dx*dy),[N])
      K2D1N=D1N
      K2D1M1=D1M1
      K2D1M2=D1M2
      K2D1S1=D1S1
      K2D1S2=D1S2
      K2D1S3=D1S3
      K2D1=D1
      
      K2D2N=D2N
      K2D2M1=D2M1
      K2D2M2=D2M2
      K2D2S1=D2S1
      K2D2S2=D2S2
      K2D2S3=D2S3
      K2D2=D2
      
      K2D3N=D3N
      K2D3M1=D3M1
      K2D3M2=D3M2
      K2D3S1=D3S1
      K2D3S2=D3S2
      K2D3S3=D3S3
      K2D3=D3

      ! Multiply D matrix by coefficient K
      CALL DOT_ROW(D1N,D1M1,D1M2,D1S1,D1S2,D1S3,D1,A*K1,N)
      CALL DOT_ROW(D2N,D2M1,D2M2,D2S1,D2S2,D2S3,D2,B*K1,N)
      CALL DOT_ROW(D3N,D3M1,D3M2,D3S1,D3S2,D3S3,D3,K1,N)
      CALL DOT_ROW(K2D1N,K2D1M1,K2D1M2,K2D1S1,K2D1S2,K2D1S3,K2D1,A*K2,N)
      CALL DOT_ROW(K2D2N,K2D2M1,K2D2M2,K2D2S1,K2D2S2,K2D2S3,K2D2,B*K2,N)
      CALL DOT_ROW(K2D3N,K2D3M1,K2D3M2,K2D3S1,K2D3S2,K2D3S3,K2D3,K2,N)
      
      !Matrix concatenation
      !M1,D2,D3,D1,A1,B1,C1
      CALL SPLICE33(M1N,M1M1,M1M2,M1S1,M1S2,M1S3,M1,D2N,D2M1,D2M2,D2S1,D2S3,D2,&
          D3N,D3M1,D3M2,D3S1,D3S3,D3,D1N,D1M1,D1M2,D1S1,D1S3,D1,&
          A1N,A1M1,A1M2,A1S1,A1S3,A1,B1N,B1M1,B1M2,B1S1,B1S3,B1,&
          C1N,C1M1,C1M2,C1S1,C1S3,C1)
      !M2,K2D2,K2D3,K2D1,A2,B2,C2
      IF(NLABEL/=0)THEN
      CALL SPLICE33(M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2,K2D2N,K2D2M1,K2D2M2,K2D2S1,K2D2S3,K2D2,&
          K2D3N,K2D3M1,K2D3M2,K2D3S1,K2D3S3,K2D3,K2D1N,K2D1M1,K2D1M2,K2D1S1,K2D1S3,K2D1,&
          A2N,A2M1,A2M2,A2S1,A2S3,A2,B2N,B2M1,B2M2,B2S1,B2S3,B2,C2N,C2M1,C2M2,C2S1,C2S3,C2)
      
      M2=-M2
      ELSE
          M2N(1)=M1N(1)
          M2N(2)=0
      ENDIF

      END
     
      end module matrix_construction
