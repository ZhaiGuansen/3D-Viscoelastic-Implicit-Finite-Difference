      module PML
      USE DIA_SM
      use INI
      use omp_lib
      USE matrix_construction


      contains
      SUBROUTINE BASE_CUBE(NX,NY,NZ,L,VP,VS,RHO,DX,DY,DZ,DT,LABEL,PMLD)
       ! Classifies model blocks by labeling PML/PDC boundaries,
      ! outputs BASE files with LABEL tags for accelerated processing.
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::D1,D2,D3,N,NUM,LX,LY,LZ,NX,NY,NZ,L,PMLD(3)
      INTEGER::X,Y,Z,I,J,K,F
      REAL::DX,DY,DZ,DT
      REAL,DIMENSION(NX,NY,NZ)::VP,VS,RHO
      INTEGER,DIMENSION(:),allocatable::LABEL
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::CTEMP
      integer,allocatable::PDC(:,:,:)
      
      
      LX=(NX-3)/(L-3)
      LY=(NY-3)/(L-3)
      LZ=(NZ-3)/(L-3)
      NUM=LX*LY*LZ
      ALLOCATE(PDC(L-2,L-2,L-2))

      !$OMP PARALLEL PRIVATE(N,X,Y,Z,F,PDC,I,J,K,D1,D2,D3,CTEMP,FILENAME)
      !$omp do  schedule(dynamic)
      DO N=1,NUM
      Z=(N-1)/(LX*LY)
      X=N-Z*LX*LY
      Z=Z+1
      Y=(X-1)/LX
      X=X-Y*LX
      Y=Y+1
      !Set as top-left
      X=(L-3)*(X-1)+1
      Y=(L-3)*(Y-1)+1
      Z=(L-3)*(Z-1)+1
      
      F=0
      PDC=0
       !Determine PDC information (reflective interface direction: binary-encoded ZYX)
      do k=1,L-2
          do j=1,L-2
              do i=1,L-2
                  IF(VP(X+I,Y+J,Z+K)/=VP(X+I-1,Y+J,Z+K))PDC(I,J,K)=PDC(I,J,K)+1    
                  IF(VP(X+I,Y+J,Z+K)/=VP(X+I,Y+J-1,Z+K))PDC(I,J,K)=PDC(I,J,K)+2
                  IF(VP(X+I,Y+J,Z+K)/=VP(X+I,Y+J,Z+K-1))PDC(I,J,K)=PDC(I,J,K)+4
                  IF(PDC(I,J,K)==7)PDC(I,J,K)=0
              ENDDO
          ENDDO
      ENDDO
      do k=1,L-2
          do j=1,L-2
              do i=1,L-2
                  IF(F==1)EXIT
                  IF(PDC(I,J,K)>0)then
                  LABEL(N)=LABEL(N)+2
                  F=1
                  ENDIF
              ENDDO
          ENDDO
      ENDDO
      
      !Determine PML information (XYZ converted to center point)
      D1=0
      D2=0
      D3=0
      IF(PMLD(1)>NX/2.OR.PMLD(2)>NY/2.OR.PMLD(3)>NZ/2)then
          WRITE(*,*)'PML thickness setting is invalid'
          call flush()
          STOP
      endif
      !Negative value = smaller-number side
      !Positive value = larger-number side
      !Absolute value = penetration depth into PML layer
      IF(X<PMLD(1))D1=X-PMLD(1)-1
      IF(X>NX-L-PMLD(1)+2)D1=X-NX-1+L+PMLD(1)
      IF(Y<PMLD(2))D2=Y-PMLD(2)-1
      IF(Y>NY-L-PMLD(2)+2)D2=Y-NY-1+L+PMLD(2)
      IF(Z<PMLD(3))D3=Z-PMLD(3)-1
      IF(Z>NZ-L-PMLD(3)+2)D3=Z-NZ-1+L+PMLD(3)
      
      IF((ABS(D1)+ABS(D2)+ABS(D3)>0).AND.LABEL(N)<4)LABEL(N)=LABEL(N)+4

      WRITE(CTEMP,'(I6)')N
      FILENAME='./BASE/BASE'//TRIM(ADJUSTL(CTEMP))//'.DAT'
      !!$OMP CRITICAL                 
      OPEN(N+100,FILE=FILENAME,STATUS='REPLACE')
      WRITE(N+100,*)PMLD(1),PMLD(2),PMLD(3)
      WRITE(N+100,*)D1,D2,D3
      WRITE(N+100,*)VP(X+L/2,Y+L/2,Z+L/2),VS(X+L/2,Y+L/2,Z+L/2)
      WRITE(N+100,*)DX,DY
      WRITE(N+100,*)DZ,DT
      WRITE(N+100,*)RHO(X+L/2,Y+L/2,Z+L/2)
      do k=1,L-2
          do j=1,L-2
              do i=1,L-2
                  WRITE(N+100,*)PDC(I,J,K)
              ENDDO
          ENDDO
      ENDDO
      do k=1,L
          do j=1,L
              do i=1,L
                  WRITE(N+100,*)VP(X-1+I,Y-1+J,Z-1+K)
                  WRITE(N+100,*)VS(X-1+I,Y-1+J,Z-1+K)
                  WRITE(N+100,*)RHO(X-1+I,Y-1+J,Z-1+K)
              ENDDO
          ENDDO
      ENDDO
      CLOSE(N+100)
      !!$OMP END CRITICAL
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
      END
      
      
      SUBROUTINE getPM(LX,LY,LZ,I,J,K,vp,vs,np,ns,DX,DY,DZ,DT,FP,PMLD,DL,PM1I,PM1J,PM13,PM2I,PM2J,PM23)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::LX,LY,LZ,I,J,K
      INTEGER,DIMENSION(3)::PMLD,DL
      REAL::X,DX,DY,DZ,DT,FP,VP,VS,NP,NS,D,W
      REAL::X2,Y2,Z2,XY,XZ,YZ,P1,P2,S1,S2
      REAL(KIND=RK)::M,N,M1,M2,SUM1,SUM2
      REAL(KIND=RK),DIMENSION(3)::A,B2,B1
      INTEGER,ALLOCATABLE,DIMENSION(:)::PM1I,PM1J,PM2I,PM2J
      REAL(KIND=RK),ALLOCATABLE,DIMENSION(:)::PM13,PM23
      
      W=-(2*FP)**2
      !W=-(2*PI*FP)**2
      DO NUM=1,3
          X=REAL(DL(NUM))/PMLD(NUM)
          D=VP*8*(0.25*X+0.75*X**2)/PMLD(NUM)
          A(NUM)=W/(D**2-W)!A2=-A1
          B1(NUM)=D/(D**2-W)
          B2(NUM)=2*W*D/(D**2-W)**2
      ENDDO
      
      X2=DT/DX/DX
      Y2=DT/DY/DY
      Z2=DT/DZ/DZ
      XY=DT/4/DX/DY
      XZ=DT/4/DX/DZ
      YZ=DT/4/DY/DZ
      P1=VP*VP
      P2=NP*NP
      S1=VS*VS
      S2=NS*NS
      
      DO NUM=1,15
          PM1I=[PM1I,POSITION(1,I,J,K,LX,LY,LZ)]
          PM2I=[PM2I,POSITION(1,I,J,K,LX,LY,LZ)]
      ENDDO
      
      M=S1*DT*Z2
      N=S2*Z2/2
      M1=M*(-A(3)-1)
      M2=N*(-A(3)-1)+M*B2(3)
      PM1J=[PM1J,POSITION(1,I,J,K-1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I,J,K-1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(1,I,J,K+1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I,J,K+1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      SUM1=M1
      SUM2=M2
      
      M=S1*DT*Y2
      N=S2*Y2/2
      M1=M*(-A(2)-1)
      M2=N*(-A(2)-1)+M*B2(2)
      PM1J=[PM1J,POSITION(1,I,J-1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I,J-1,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(1,I,J+1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I,J+1,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      SUM1=SUM1+M1
      SUM2=SUM2+M2
      
      M=P1*DT*X2
      N=P2*X2/2
      M1=M*(-A(1)-1)
      M2=N*(-A(1)-1)+M*B2(1)
      PM1J=[PM1J,POSITION(1,I-1,J,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I-1,J,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(1,I+1,J,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I+1,J,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      SUM1=SUM1+M1
      SUM2=SUM2+M2

      PM1J=[PM1J,POSITION(1,I,J,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I,J,K,LX,LY,LZ)]
      PM13=[PM13,-2*SUM1]
      PM23=[PM23,-2*SUM2]
      
      M=(P1-S1)*DT*XY
      N=(P2-S2)*XY/2
      M1=M*(A(1)*A(2)-1)
      M2=N*(A(1)*A(2)-1)+M*(A(1)*B1(2)+A(2)*B1(1))
      PM1J=[PM1J,POSITION(2,I-1,J-1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I-1,J-1,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(2,I+1,J-1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I+1,J-1,K,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(2,I-1,J+1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I-1,J+1,K,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(2,I+1,J+1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I+1,J+1,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      
      M=(P1-S1)*DT*XZ
      N=(P2-S2)*XZ/2
      M1=M*(A(1)*A(3)-1)
      M2=N*(A(1)*A(3)-1)+M*(A(1)*B1(3)+A(3)*B1(1))
      PM1J=[PM1J,POSITION(3,I-1,J,K-1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I-1,J,K-1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(3,I+1,J,K-1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I+1,J,K-1,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(3,I-1,J,K+1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I-1,J,K+1,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(3,I+1,J,K+1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I+1,J,K+1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      
      
      DO NUM=1,15
          PM1I=[PM1I,POSITION(2,I,J,K,LX,LY,LZ)]
          PM2I=[PM2I,POSITION(2,I,J,K,LX,LY,LZ)]
      ENDDO
      
      M=S1*DT*Z2
      N=S2*Z2/2
      M1=M*(-A(3)-1)
      M2=N*(-A(3)-1)+M*B2(3)
      PM1J=[PM1J,POSITION(2,I,J,K-1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I,J,K-1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(2,I,J,K+1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I,J,K+1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      SUM1=M1
      SUM2=M2
      
      M=P1*DT*Y2
      N=P2*Y2/2
      M1=M*(-A(2)-1)
      M2=N*(-A(2)-1)+M*B2(2)
      PM1J=[PM1J,POSITION(2,I,J-1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I,J-1,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(2,I,J+1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I,J+1,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      SUM1=SUM1+M1
      SUM2=SUM2+M2
      
      M=S1*DT*X2
      N=S2*X2/2
      M1=M*(-A(1)-1)
      M2=N*(-A(1)-1)+M*B2(1)
      PM1J=[PM1J,POSITION(2,I-1,J,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I-1,J,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(2,I+1,J,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I+1,J,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      SUM1=SUM1+M1
      SUM2=SUM2+M2

      PM1J=[PM1J,POSITION(2,I,J,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I,J,K,LX,LY,LZ)]
      PM13=[PM13,-2*SUM1]
      PM23=[PM23,-2*SUM2]
      
      M=(P1-S1)*DT*XY
      N=(P2-S2)*XY/2
      M1=M*(A(1)*A(2)-1)
      M2=N*(A(1)*A(2)-1)+M*(A(1)*B1(2)+A(2)*B1(1))
      PM1J=[PM1J,POSITION(1,I-1,J-1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I-1,J-1,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(1,I+1,J-1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I+1,J-1,K,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(1,I-1,J+1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I-1,J+1,K,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(1,I+1,J+1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I+1,J+1,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      
      M=(P1-S1)*DT*YZ
      N=(P2-S2)*YZ/2
      M1=M*(A(2)*A(3)-1)
      M2=N*(A(2)*A(3)-1)+M*(A(2)*B1(3)+A(3)*B1(2))
      PM1J=[PM1J,POSITION(3,I,J-1,K-1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I,J-1,K-1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(3,I,J+1,K-1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I,J+1,K-1,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(3,I,J-1,K+1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I,J-1,K+1,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(3,I,J+1,K+1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I,J+1,K+1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      
      
      DO NUM=1,15
          PM1I=[PM1I,POSITION(3,I,J,K,LX,LY,LZ)]
          PM2I=[PM2I,POSITION(3,I,J,K,LX,LY,LZ)]
      ENDDO
      
      M=P1*DT*Z2
      N=P2*Z2/2
      M1=M*(-A(3)-1)
      M2=N*(-A(3)-1)+M*B2(3)
      PM1J=[PM1J,POSITION(3,I,J,K-1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I,J,K-1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(3,I,J,K+1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I,J,K+1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      SUM1=M1
      SUM2=M2
      
      M=S1*DT*Y2
      N=S2*Y2/2
      M1=M*(-A(2)-1)
      M2=N*(-A(2)-1)+M*B2(2)
      PM1J=[PM1J,POSITION(3,I,J-1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I,J-1,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(3,I,J+1,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I,J+1,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      SUM1=SUM1+M1
      SUM2=SUM2+M2
      
      M=S1*DT*X2
      N=S2*X2/2
      M1=M*(-A(1)-1)
      M2=N*(-A(1)-1)+M*B2(1)
      PM1J=[PM1J,POSITION(3,I-1,J,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I-1,J,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(3,I+1,J,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I+1,J,K,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      SUM1=SUM1+M1
      SUM2=SUM2+M2

      PM1J=[PM1J,POSITION(3,I,J,K,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(3,I,J,K,LX,LY,LZ)]
      PM13=[PM13,-2*SUM1]
      PM23=[PM23,-2*SUM2]
      
      M=(P1-S1)*DT*YZ
      N=(P2-S2)*YZ/2
      M1=M*(A(3)*A(2)-1)
      M2=N*(A(3)*A(2)-1)+M*(A(3)*B1(2)+A(2)*B1(3))
      PM1J=[PM1J,POSITION(2,I,J-1,K-1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I,J-1,K-1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(2,I,J-1,K+1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I,J-1,K+1,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(2,I,J+1,K-1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I,J+1,K-1,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(2,I,J+1,K+1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(2,I,J+1,K+1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      
      M=(P1-S1)*DT*XZ
      N=(P2-S2)*XZ/2
      M1=M*(A(1)*A(3)-1)
      M2=N*(A(1)*A(3)-1)+M*(A(1)*B1(3)+A(3)*B1(1))
      PM1J=[PM1J,POSITION(1,I-1,J,K-1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I-1,J,K-1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      PM1J=[PM1J,POSITION(1,I+1,J,K-1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I+1,J,K-1,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(1,I-1,J,K+1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I-1,J,K+1,LX,LY,LZ)]
      PM13=[PM13,-M1]
      PM23=[PM23,-M2]
      PM1J=[PM1J,POSITION(1,I+1,J,K+1,LX,LY,LZ)]
      PM2J=[PM2J,POSITION(1,I+1,J,K+1,LX,LY,LZ)]
      PM13=[PM13,M1]
      PM23=[PM23,M2]
      
      
      END
      
      
      
      
      
      
      
      
      
      END module PML