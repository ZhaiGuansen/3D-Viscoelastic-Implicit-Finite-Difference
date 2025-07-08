      module PML
      USE DIA_SM
      use INI
      use omp_lib
      USE matrix_construction
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)   
      contains
      

      SUBROUTINE BASE_CUBE(NX,NY,NZ,L,VP,VS,RHO,DX,DY,DZ,DT,LABEL,PMLD)
       ! Classifies model blocks by labeling PML/PDC boundaries,
      ! outputs BASE files with LABEL tags for accelerated processing.
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::D1,D2,D3,N,NUM,NX,NY,NZ,LX,LY,LZ,L,PMLD(3)
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
      !$OMP PARALLEL PRIVATE(N,X,Y,Z,F,PDC,I,J,K,
     + D1,D2,D3,CTEMP,FILENAME)
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
                  IF(VP(X+I,Y+J,Z+K)/=VP(X+I-1,Y+J,Z+K))
     + PDC(I,J,K)=PDC(I,J,K)+1    
                  IF(VP(X+I,Y+J,Z+K)/=VP(X+I,Y+J-1,Z+K))
     + PDC(I,J,K)=PDC(I,J,K)+2
                  IF(VP(X+I,Y+J,Z+K)/=VP(X+I,Y+J,Z+K-1))
     + PDC(I,J,K)=PDC(I,J,K)+4
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
          PRINT*,'PML thickness setting is invalid'
          pause
      endif
      
      IF(X<=PMLD(1)-L+2.OR.X>=NX-PMLD(1))D1=1! Entirely within PML region
      IF(Y<=PMLD(2)-L+2.OR.Y>=NY-PMLD(2))D2=1
      IF(Z<=PMLD(3)-L+2.OR.Z>=NZ-PMLD(3))D3=1
      
      IF((D1+D2+D3>0).AND.LABEL(N)<4)
     + LABEL(N)=LABEL(N)+4

      WRITE(CTEMP,'(I6)')N
      FILENAME='.\BASE\BASE'//TRIM(ADJUSTL(CTEMP))//'.DAT'
      !!$OMP CRITICAL                 
      OPEN(N,FILE=FILENAME,STATUS='REPLACE')
      WRITE(N,*)PMLD(1),PMLD(2),PMLD(3)
      WRITE(N,*)VP(X+L/2,Y+L/2,Z+L/2),VS(X+L/2,Y+L/2,Z+L/2)
      WRITE(N,*)DX,DY
      WRITE(N,*)DZ,DT
      WRITE(N,*)RHO(X+L/2,Y+L/2,Z+L/2)
      do k=1,L-2
          do j=1,L-2
              do i=1,L-2
                  WRITE(N,*)PDC(I,J,K)
              ENDDO
          ENDDO
      ENDDO
      do k=1,L
          do j=1,L
              do i=1,L
                  WRITE(N,*)VP(X-1+I,Y-1+J,Z-1+K)
                  WRITE(N,*)VS(X-1+I,Y-1+J,Z-1+K)
                  WRITE(N,*)RHO(X-1+I,Y-1+J,Z-1+K)
              ENDDO
          ENDDO
      ENDDO
      CLOSE(N)
      !!$OMP END CRITICAL
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
      END
      
      
      SUBROUTINE GET_PML(T,NX,NY,NZ,L,PMLD,IS_PML,DX,DY,DZ,DT)! Boundary indicators (xyz directions):D1/D2/D3 = -1(start)/1(end)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::X,Y,Z,NX,NY,NZ,N,LX,LY,LZ,IS_PML,L,ISZERO
      INTEGER::T,I,J,K,M,NUM,IOS,e,PMLD(3),IS_ZERO(3)
      REAL::DX,DY,DZ,DT
      INTEGER,DIMENSION(:),ALLOCATABLE::D1,D2,D3
      REAL(KIND=RK),DIMENSION(:,:,:,:),allocatable::U1,U2,U3![1:9]  - Component velocities [10:12] - Resultant vectors
      REAL(KIND=RK)::VALUE(13),COEF(13),C(7)
      CHARACTER(LEN=30)::FILENAME1,FILENAME2
      CHARACTER(LEN=10)::C1,C2
      logical :: file_exists
      integer(kind=8)  ::count1,count2,count_rate,count_max
      integer(kind=8)  ::count3
      real timespend
      
      ALLOCATE(U1(NX,3,NZ,12),U2(NX,3,NZ,12),U3(NX,3,NZ,12))
      ALLOCATE(D1(NX),D2(NY),D3(NZ))
      D1=0
      D2=0
      D3=0
      DO X=1,NX
              IF(X<PMLD(1)+1)D1(X)=PMLD(1)+1-X
              IF(NX-X<PMLD(1))D1(X)=PMLD(1)+X-NX
      ENDDO
      DO Y=1,NY
              IF(Y<PMLD(2)+1)D2(Y)=PMLD(2)+1-Y
              IF(NY-Y<PMLD(2))D2(Y)=PMLD(2)+Y-NY
      ENDDO
      DO Z=1,NZ
              IF(Z<PMLD(3)+1)D3(Z)=PMLD(3)+1-Z
              IF(NZ-Z<PMLD(3))D3(Z)=PMLD(3)+Z-NZ
      ENDDO
      
      LX=(NX-3)/(L-3)
      LY=(NY-3)/(L-3)
      LZ=(NZ-3)/(L-3)
      IS_ZERO=1
      DO Y=1,NY

          U3=0        
          IF(T<5)THEN
              CALL OUTPUT_PML(T,Y,NY,U3,PMLD,1)
              CYCLE
          ENDIF
          
          CALL READ_PML(T,Y,NY,U1,U2,PMLD,IS_ZERO)
          
          IF(IS_ZERO(1)+IS_ZERO(2)+IS_ZERO(3)==3)THEN
              CALL OUTPUT_PML(T,Y,NY,U3,PMLD,1)
              CYCLE
          ENDIF
          IS_PML=1
          ISZERO=1

          !$OMP PARALLEL PRIVATE(I,J,K,N,X,Z,NUM)
          !$omp do SCHEDULE(DYNAMIC)
          DO NUM=1,NX*NZ
              Z=1+(NUM-1)/NX
              X=NUM-(Z-1)*NX              
                  
              IF(D1(X)+D2(Y)+D3(Z)>0)THEN
                  
              I=X/(L-3)
              J=Y/(L-3)
              K=Z/(L-3)
              IF(I<1)I=1
              IF(I>LX)I=LX
              IF(J<1)J=1
              IF(J>LY)J=LY
              IF(K<1)K=1
              IF(K>LZ)K=LZ
              N=I+(J-1)*LX+(K-1)*LX*LY
              
              CALL PML_CUBE(N,U1,U2,U3,X,NX,Z,D1(X),D2(Y),D3(Z),DX,DY
     +           ,DZ,DT,ISZERO)!rewrite pml_cube
              ENDIF
          ENDDO
          !$OMP END DO
          !$OMP END PARALLEL
          !U3_output
          CALL OUTPUT_PML(T,Y,NY,U3,PMLD,ISZERO)
     
      ENDDO

      END
      
      
      SUBROUTINE READ_PML(T,Y,NY,U1,U2,PMLD,IS_ZERO)! Read columns y-1, y, y+1
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::X,Y,Z,NX,NY,NZ,T,PMLD(3)
      INTEGER::I,J,K,M,IS_ZERO(3),ISZERO(2)
      REAL(KIND=RK),DIMENSION(:,:,:,:),allocatable::U1,U2
      CHARACTER(LEN=30)::FILENAME1,FILENAME2
      CHARACTER(LEN=10)::C,C1,C2,C3
      LOGICAL::file_exists
      NX=SIZE(U1,1)
      NZ=SIZE(U1,3)
      WRITE(C1,'(I6)')T-2
      WRITE(C2,'(I6)')T-1
      WRITE(C3,'(I6)')T-3

      IF(Y==1)THEN
          U1=0
          U2=0
          WRITE(C,'(I6)')1
          FILENAME1='.\PML\PML'//TRIM(ADJUSTL(C))//'_'
     + //TRIM(ADJUSTL(C1))//'.DAT'
          FILENAME2='.\PML\PML'//TRIM(ADJUSTL(C))//'_'
     + //TRIM(ADJUSTL(C2))//'.DAT'
          OPEN(1,FILE=FILENAME1,STATUS='OLD')
          OPEN(NY+1,FILE=FILENAME2,STATUS='OLD')
          IS_ZERO(2)=0
          READ(1,*)ISZERO(1)
          READ(NY+1,*)ISZERO(2)
          IF(ISZERO(1)+ISZERO(2)==2)THEN
              IS_ZERO(2)=1
          ELSE
              IS_ZERO(2)=0
          DO Z=1,NZ
              DO X=1,NX
                  DO M=1,12
                      IF(ISZERO(1)/=1)READ(1,*)U1(X,2,Z,M)
                      IF(ISZERO(2)/=1)READ(NY+1,*)U2(X,2,Z,M)
                  ENDDO
              ENDDO
          ENDDO
          ENDIF
          CLOSE(1)
          CLOSE(NY+1)
          FILENAME1='.\PML\PML'//TRIM(ADJUSTL(C))//'_'
     + //TRIM(ADJUSTL(C3))//'.DAT'
          INQUIRE(FILE=FILENAME1,EXIST=file_exists)
          IF(file_exists)THEN
              OPEN(1,FILE=FILENAME1,STATUS='OLD')
              CLOSE(1,STATUS='DELETE')
          ENDIF
      ELSE
          IS_ZERO(1)=IS_ZERO(2)
          IS_ZERO(2)=IS_ZERO(3)
          IS_ZERO(3)=1
          U1(:,1,:,:)=U1(:,2,:,:)
          U1(:,2,:,:)=U1(:,3,:,:)
          U2(:,1,:,:)=U2(:,2,:,:)
          U2(:,2,:,:)=U2(:,3,:,:)
      ENDIF
      

      IF(Y==NY)THEN
          U1(:,3,:,:)=0
          U2(:,3,:,:)=0
      ELSE
          WRITE(C,'(I6)')Y+1
          FILENAME1='.\PML\PML'//TRIM(ADJUSTL(C))//'_'
     + //TRIM(ADJUSTL(C1))//'.DAT'
          FILENAME2='.\PML\PML'//TRIM(ADJUSTL(C))//'_'
     + //TRIM(ADJUSTL(C2))//'.DAT'
          OPEN(Y+1,FILE=FILENAME1,STATUS='OLD')
          OPEN(NY+Y+1,FILE=FILENAME2,STATUS='OLD')
          
          IS_ZERO(3)=0
          READ(Y+1,*)ISZERO(1)
          READ(NY+Y+1,*)ISZERO(2)
          IF(ISZERO(1)+ISZERO(2)==2)THEN
              IS_ZERO(3)=1
          ELSE
              IS_ZERO(3)=0
          DO Z=1,NZ
          DO X=1,NX
              DO M=1,12
                  IF(ISZERO(1)/=1)READ(Y+1,*)U1(X,3,Z,M)
                  IF(ISZERO(2)/=1)READ(NY+Y+1,*)U2(X,3,Z,M)
              ENDDO
          ENDDO
          ENDDO
          ENDIF
          

          CLOSE(Y+1)
          CLOSE(NY+Y+1)
          FILENAME1='.\PML\PML'//TRIM(ADJUSTL(C))//'_'
     + //TRIM(ADJUSTL(C3))//'.DAT'
          INQUIRE(FILE=FILENAME1,EXIST=file_exists)
          IF(file_exists)THEN
              OPEN(1,FILE=FILENAME1,STATUS='OLD')
              CLOSE(1,STATUS='DELETE')
          ENDIF
      ENDIF
      
      
      
      END
      
      
      
      SUBROUTINE OUTPUT_PML(T,Y,NY,U,PMLD,ISZERO)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::X,Y,Z,NX,NY,NZ,T,PMLD(3),ISZERO
      INTEGER::I,J,K,M
      REAL(KIND=RK),DIMENSION(:,:,:,:),allocatable::U
      REAL(KIND=RK),DIMENSION(:,:,:),allocatable::OUTPUT
      CHARACTER(LEN=30)::FILENAME1,FILENAME2
      CHARACTER(LEN=10)::C1,C2
      
      NX=SIZE(U,1)
      NZ=SIZE(U,3)
      ALLOCATE(OUTPUT(NX,NZ,3))
      
      WRITE(C1,'(I6)')Y
      WRITE(C2,'(I6)')T
      FILENAME1='.\PML\PML'//TRIM(ADJUSTL(C1))//'_'
     + //TRIM(ADJUSTL(C2))//'.DAT'
      FILENAME2='.\Y\Y'//TRIM(ADJUSTL(C2))//'_'
     + //TRIM(ADJUSTL(C1))//'.DAT'

      OPEN(Y+1,FILE=FILENAME1)
      

      IF(Y<PMLD(2)+1.OR.Y>NY-PMLD(2))THEN 
          WRITE(Y+1,*)ISZERO
          IF(ISZERO/=1)THEN
          DO Z=1,NZ
          DO X=1,NX
              DO M=1,12
                  WRITE(Y+1,*)U(X,2,Z,M)
              ENDDO
          ENDDO
          ENDDO
          ENDIF
      ELSE
          WRITE(Y+1,*)0
          OPEN(NY+Y+1,FILE=FILENAME2,STATUS='OLD')
          DO M=1,3
              DO Z=1,NZ
                  DO X=1,NX
                      READ(NY+Y+1,*)OUTPUT(X,Z,M)
                      IF(ABS(OUTPUT(X,Z,M))<1E-10)OUTPUT(X,Z,M)=0
                  ENDDO
              ENDDO
          ENDDO
          
          
          DO Z=1,NZ
          DO X=1,NX
              DO M=1,9
                  WRITE(Y+1,*)U(X,2,Z,M)
              ENDDO
              D1=0
              IF(X<PMLD(1)+1)D1=PMLD(1)+1-X
              IF(NX-X<PMLD(1))D1=PMLD(1)+X-NX
              D3=0
              IF(Z<PMLD(3)+1)D3=PMLD(3)+1-Z
              IF(NZ-Z<PMLD(3))D3=PMLD(3)+Z-NZ
              DO M=10,12
                  IF(D1+D3>0)THEN   
                      WRITE(Y+1,*)U(X,2,Z,M)
                  ELSE  
                      WRITE(Y+1,*)OUTPUT(X,Z,M-9)
                  ENDIF
              ENDDO
          ENDDO
          ENDDO
      ENDIF
      
      CLOSE(Y+1)
      CLOSE(NY+Y+1)
      END
      
      
      
      
      
      
      
      
      
      SUBROUTINE PML_CUBE(N,U1,U2,U3,X,NX,Z,D1,D2,D3,DX,DY,DZ,DT,
     + ISZERO)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::X,Z,T,I,J,K,M,D1,D2,D3,N,NUM,IOS,PMLD(3),ISZERO
      REAL::VP,VS,DX,DY,DZ,DT,PMLC(3)
      REAL(KIND=RK),DIMENSION(:,:,:,:),allocatable::U1,U2,U3
      REAL(KIND=RK)::VALUE(13),COEF(13),C(7)
      REAL(KIND=RK)::P2,S2,X2,Y2,Z2,XY,XZ,YZ,D0,D(3),K0,K2
      CHARACTER(LEN=30)::BASE
      CHARACTER(LEN=10)::C1


      WRITE(C1,'(I6)')N
      NUM=X+(Z-1)*NX
      BASE='.\BASE\BASE'//TRIM(ADJUSTL(C1))//'.DAT'
      OPEN(NUM,FILE=BASE,STATUS='OLD',IOSTAT=IOS)
      !READ(NUM,*,IOSTAT=IOS)X,Y,Z
      !READ(NUM,*,IOSTAT=IOS)NX,NY,NZ
      READ(NUM,*,IOSTAT=IOS)PMLD(1),PMLD(2),PMLD(3)
      READ(NUM,*,IOSTAT=IOS)VP,VS
      CLOSE(NUM,IOSTAT=IOS)

      
      !IF(X==6.AND.Z==7)THEN
      !    PRINT*,"-----------"
      !ENDIF
      

      P2=VP*VP
      S2=VS*VS
      X2=DX*DX/(DT*DT)
      Y2=DY*DY/(DT*DT)
      Z2=DZ*DZ/(DT*DT)
      XY=DX*DY/(DT*DT)
      XZ=DX*DZ/(DT*DT)
      YZ=DY*DZ/(DT*DT)
      D0=9*VP/4
      K0=0
      K2=2
      DO I=1,3
          PMLC(I)=4.5*VP/(PMLD(I))**3!!4.5=3/2*log(R)
      ENDDO

             
      ! Compute absorption coefficient
      
      D(1)=D1**2*PMLC(1)*DT/2
      D(2)=D2**2*PMLC(2)*DT/2
      D(3)=D3**2*PMLC(3)*DT/2
      C(2)=(P2-2*S2)/4/XY
      C(3)=(P2-2*S2)/4/XZ
      C(4)=(P2-2*S2)/4/YZ
      C(5)=S2/4/XY
      C(6)=S2/4/XZ
      C(7)=S2/4/YZ
      ! Evaluate value
      
      I=X
      J=2
      K=Z
      
      
      VALUE=0
      CALL GET_VALUE_CUBE(U2,I+1,J,K,10,VALUE(1))
      CALL GET_VALUE_CUBE(U2,I-1,J,K,10,VALUE(2))
      VALUE(3)=U2(I,J,K,10)
      VALUE(4)=U2(I,J,K,1)
      VALUE(5)=U1(I,J,K,1)
      CALL GET_VALUE_CUBE(U2,I+1,J+1,K,11,VALUE(6))
      CALL GET_VALUE_CUBE(U2,I+1,J-1,K,11,VALUE(7))
      CALL GET_VALUE_CUBE(U2,I-1,J+1,K,11,VALUE(8))
      CALL GET_VALUE_CUBE(U2,I-1,J-1,K,11,VALUE(9))
      CALL GET_VALUE_CUBE(U2,I+1,J,K+1,12,VALUE(10))
      CALL GET_VALUE_CUBE(U2,I+1,J,K-1,12,VALUE(11))
      CALL GET_VALUE_CUBE(U2,I-1,J,K+1,12,VALUE(12))
      CALL GET_VALUE_CUBE(U2,I-1,J,K-1,12,VALUE(13))
      
      C(1)=p2/x2

      coef=[C(1),C(1),-2*C(1),K2,D(1)-1,C(2),-C(2),-C(2),C(2),
     +C(3),-C(3),-C(3),C(3)]
      COEF=COEF/(1+D(1))
      U3(I,J,K,1)=dot_product(coef,value)
      
      VALUE=0
      VALUE(1)=U2(I,J+1,K,10)
      VALUE(2)=U2(I,J-1,K,10)
      VALUE(3)=U2(I,J,K,10)
      VALUE(4)=U2(I,J,K,2)
      VALUE(5)=U1(I,J,K,2)
      CALL GET_VALUE_CUBE(U2, I+1,J+1,K,11,VALUE(6))
      CALL GET_VALUE_CUBE(U2,  I+1,J-1,K,11,VALUE(7))
      CALL GET_VALUE_CUBE(U2,I-1,J+1,K,11,VALUE(8))
      CALL GET_VALUE_CUBE(U2, I-1,J-1,K,11,VALUE(9))
  
      C(1)=S2/Y2
      coef=[C(1),C(1),-2*C(1),K2,D(2)-1,
     +             C(5),-C(5),-C(5),C(5),K0,K0,K0,K0]
      COEF=COEF/(1+D(2))

      U3(I,J,K,2)=dot_product(coef,value)   
      
      VALUE=0
      CALL GET_VALUE_CUBE(U2,I,J,K+1,10,VALUE(1))
      CALL GET_VALUE_CUBE(U2,I,J,K-1,10,VALUE(2))
      VALUE(3)=U2(I,J,K,10)
      VALUE(4)=U2(I,J,K,3)
      VALUE(5)=U1(I,J,K,3)
      CALL GET_VALUE_CUBE(U2,I+1,J,K+1,12,VALUE(6))
      CALL GET_VALUE_CUBE(U2,I+1,J,K-1,12,VALUE(7))
      CALL GET_VALUE_CUBE(U2,I-1,J,K+1,12,VALUE(8))
      CALL GET_VALUE_CUBE(U2,I-1,J,K-1,12,VALUE(9))
      
      C(1)=S2/Z2
      coef=[C(1),C(1),-2*C(1),K2,D(3)-1,
     +    C(6),-C(6),-C(6),C(6),K0,K0,K0,K0]
      COEF=COEF/(1+D(3))
      U3(I,J,K,3)=dot_product(coef,value)
      
      VALUE=0
      CALL GET_VALUE_CUBE(U2,I+1,J,K,11,VALUE(1))
      CALL GET_VALUE_CUBE(U2,I-1,J,K,11,VALUE(2))
      VALUE(3)=U2(I,J,K,11)
      VALUE(4)=U2(I,J,K,4)
      VALUE(5)=U1(I,J,K,4)
      CALL GET_VALUE_CUBE(U2,I+1,J+1,K,10,VALUE(6))
      CALL GET_VALUE_CUBE(U2,I+1,J-1,K,10,VALUE(7))
      CALL GET_VALUE_CUBE(U2,I-1,J+1,K,10,VALUE(8))
      CALL GET_VALUE_CUBE(U2,I-1,J-1,K,10,VALUE(9))  

      
      C(1)=S2/X2
      coef=[C(1),C(1),-2*C(1),K2,D(1)-1,
     +C(5),-C(5),-C(5),C(5),K0,K0,K0,K0]
      COEF=COEF/(1+D(1))
      
      U3(I,J,K,4)=dot_product(coef,value)
      
      VALUE=0
      VALUE(1)=U2(I,J+1,K,11)
      VALUE(2)=U2(I,J-1,K,11)
      VALUE(3)=U2(I,J,K,11)
      VALUE(4)=U2(I,J,K,5)
      VALUE(5)=U1(I,J,K,5)
      CALL GET_VALUE_CUBE(U2,I+1,J+1,K,10,VALUE(6))
      CALL GET_VALUE_CUBE(U2,I+1,J-1,K,10,VALUE(7))
      CALL GET_VALUE_CUBE(U2,I-1,J+1,K,10,VALUE(8))
      CALL GET_VALUE_CUBE(U2,I-1,J-1,K,10,VALUE(9))
      CALL GET_VALUE_CUBE(U2, I,J+1,K+1,12,VALUE(10))
      CALL GET_VALUE_CUBE(U2,I,J-1,K+1,12,VALUE(11))
      CALL GET_VALUE_CUBE(U2,I,J+1,K-1,12,VALUE(12))
      CALL GET_VALUE_CUBE(U2,I,J-1,K-1,12,VALUE(13))
      
      C(1)=p2/Y2
      coef=[C(1),C(1),-2*C(1),K2,D(2)-1,
     +C(2),-C(2),-C(2),C(2),
     +C(4),-C(4),-C(4),C(4)]
      COEF=COEF/(1+D(2))
      
      
      U3(I,J,K,5)=dot_product(coef,value)
      
      VALUE=0
      CALL GET_VALUE_CUBE(U2,I,J,K+1,11,VALUE(1))
      CALL GET_VALUE_CUBE(U2,I,J,K-1,11,VALUE(2))
      VALUE(3)=U2(I,J,K,11)
      VALUE(4)=U2(I,J,K,6)
      VALUE(5)=U1(I,J,K,6)
      CALL GET_VALUE_CUBE(U2,I,J+1,K+1,12,VALUE(6))
      CALL GET_VALUE_CUBE(U2,I,J-1,K+1,12,VALUE(7))
      CALL GET_VALUE_CUBE(U2,I,J+1,K-1,12,VALUE(8))
      CALL GET_VALUE_CUBE(U2,I,J-1,K-1,12,VALUE(9))
     
      C(1)=S2/Z2
      coef=[C(1),C(1),-2*C(1),K2,D(3)-1,
     +             C(7),-C(7),-C(7),C(7),K0,K0,K0,K0]
      COEF=COEF/(1+D(3))
      
      
      U3(I,J,K,6)=dot_product(coef,value)
      
      VALUE=0
      CALL GET_VALUE_CUBE(U2,I+1,J,K,12,VALUE(1))
      CALL GET_VALUE_CUBE(U2,I-1,J,K,12,VALUE(2))
      VALUE(3)=U2(I,J,K,12)
      VALUE(4)=U2(I,J,K,7)
      VALUE(5)=U1(I,J,K,7)
      CALL GET_VALUE_CUBE(U2,I+1,J,K+1,10,VALUE(6))
      CALL GET_VALUE_CUBE(U2,I+1,J,K-1,10,VALUE(7))
      CALL GET_VALUE_CUBE(U2,I-1,J,K+1,10,VALUE(8))
      CALL GET_VALUE_CUBE(U2,I-1,J,K-1,10,VALUE(9))
      
      C(1)=S2/X2
          coef=[C(1),C(1),-2*C(1),K2,D(1)-1,
     +        C(6),-C(6),-C(6),C(6),K0,K0,K0,K0]
          COEF=COEF/(1+D(1))
      
      
      U3(I,J,K,7)=dot_product(coef,value)
      
      
      VALUE=0
      VALUE(1)=U2(I,J+1,K,12)
      VALUE(2)=U2(I,J-1,K,12)
      VALUE(3)=U2(I,J,K,12)
      VALUE(4)=U2(I,J,K,8)
      VALUE(5)=U1(I,J,K,8)
      CALL GET_VALUE_CUBE(U2,I,J+1,K+1,11,VALUE(6))
      CALL GET_VALUE_CUBE(U2,I,J-1,K+1,11,VALUE(7))
      CALL GET_VALUE_CUBE(U2,I,J+1,K-1,11,VALUE(8))
      CALL GET_VALUE_CUBE(U2,I,J-1,K-1,11,VALUE(9))

          C(1)=S2/Y2
          coef=[C(1),C(1),-2*C(1),K2,D(2)-1,
     +    C(7),-C(7),-C(7),C(7),K0,K0,K0,K0]
          COEF=COEF/(1+D(2))
      
      U3(I,J,K,8)=dot_product(coef,value)
      
      VALUE=0
      CALL GET_VALUE_CUBE(U2,I,J,K+1,12,VALUE(1))
      CALL GET_VALUE_CUBE(U2,I,J,K-1,12,VALUE(2))
      VALUE(3)=U2(I,J,K,12)
      VALUE(4)=U2(I,J,K,9)
      VALUE(5)=U1(I,J,K,9)
      CALL GET_VALUE_CUBE(U2,I+1,J,K+1,10,VALUE(6))
      CALL GET_VALUE_CUBE(U2,I+1,J,K-1,10,VALUE(7))
      CALL GET_VALUE_CUBE(U2,I-1,J,K+1,10,VALUE(8))
      CALL GET_VALUE_CUBE(U2,I-1,J,K-1,10,VALUE(9))
      CALL GET_VALUE_CUBE(U2,I,J+1,K+1,11,VALUE(10))
      CALL GET_VALUE_CUBE(U2,I,J+1,K-1,11,VALUE(11))
      CALL GET_VALUE_CUBE(U2,I,J-1,K+1,11,VALUE(12))
      CALL GET_VALUE_CUBE(U2,I,J-1,K-1,11,VALUE(13))

      C(1)=p2/Z2
      coef=[C(1),C(1),-2*C(1),K2,D(3)-1,
     +C(3),-C(3),-C(3),C(3),
     +C(4),-C(4),-C(4),C(4)]
      COEF=COEF/(1+D(3))
      U3(I,J,K,9)=dot_product(coef,value)
      
      DO NUM=1,9       
          IF(ABS(U3(I,J,K,NUM))<1E-10)THEN
              U3(I,J,K,NUM)=0
          ELSE
              ISZERO=0
          ENDIF
      ENDDO
      U3(I,J,K,10)=U3(I,J,K,1)+U3(I,J,K,2)+U3(I,J,K,3)
      U3(I,J,K,11)=U3(I,J,K,4)+U3(I,J,K,5)+U3(I,J,K,6)
      U3(I,J,K,12)=U3(I,J,K,7)+U3(I,J,K,8)+U3(I,J,K,9)
                  

      
      END
      


      SUBROUTINE GET_VALUE_CUBE(U,I,J,K,N,OUTPUT)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::I,J,K,N,NX,NZ
      REAL(KIND=RK),DIMENSION(:,:,:,:),allocatable::U
      REAL(KIND=RK)::OUTPUT
      
      NX=SIZE(U,1)
      NZ=SIZE(U,3)
      IF(I<1.OR.I>NX.OR.K<1.OR.K>NZ)THEN
          OUTPUT=0
      ELSE
          OUTPUT=U(I,J,K,N)
      ENDIF

      END
      
      

      
      !SUBROUTINE FROM_U0_TO_U(U,U0,I,J,K)
      !integer, parameter :: rk = kind ( 1.0D+00 )
      !REAL(KIND=RK),DIMENSION(5,5,5,12)::U
      !REAL(KIND=RK)::U0(375)
      !INTEGER::I,J,K   
      !U(I,J,K,10)=U0(I+(J-1)*5+(K-1)*25)
      !U(I,J,K,1)=U0(I+(J-1)*5+(K-1)*25)
      !U(I,J,K,11)=U0(I+(J-1)*5+(K-1)*25+125)
      !U(I,J,K,5)=U0(I+(J-1)*5+(K-1)*25+125)
      !U(I,J,K,12)=U0(I+(J-1)*5+(K-1)*25+250)
      !U(I,J,K,9)=U0(I+(J-1)*5+(K-1)*25+250)  
      !END
      
      
      
      
      
      
      
      
      
      END module PML