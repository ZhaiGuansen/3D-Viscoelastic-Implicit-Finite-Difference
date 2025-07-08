      module cube_of_five_meshes
      !USE sparse_matrix
      USE DIA_SM
      USE LAPACK95
      use matrix_construction
      use INI
      USE PML
      USE PDC
      use, intrinsic :: ieee_arithmetic

      IMPLICIT DOUBLE PRECISION (A-H,P-Z)   
      contains
      !Convert numbering to coordinates
      SUBROUTINE N_TO_X0(N,NX,NY,NZ,L,X,Y,Z)
      INTEGER::N,NX,NY,NZ,L,LX,LY,LZ,X,Y,Z
      
      
      LX=(NX-3)/(L-3)
      LY=(NY-3)/(L-3)
      LZ=(NZ-3)/(L-3)
      Z=(N-1)/(LX*LY)
      X=N-Z*LX*LY
      Z=Z+1
      Y=(X-1)/LX
      X=X-Y*LX
      Y=Y+1
      
      X=(L-3)*(X-1)+1
      Y=(L-3)*(Y-1)+1
      Z=(L-3)*(Z-1)+1
      
      END
      
      SUBROUTINE MESH_A(N,L,DX,DY,DZ,DT,VP,VS,NP,NS,LABEL)
      ! Pre-loading coefficient matrices
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::NUM,N,I,J,K,LENGTH(3),LABEL,NLABEL
      REAL:: DX,DY,DZ,DT
      REAL,DIMENSION(L,L,L)::VP,VS,NP,NS
      INTEGER,DIMENSION(2)::AN,BN,M1N,M2N,IMN
      INTEGER,ALLOCATABLE,DIMENSION(:)::AM1,M1M1,M2M1,IMM1,AM2,M1M2,M2M2
     + ,IMM2,AS1,M1S1,M2S1,IMS1,AS2,M1S2,M2S2,IMS2,AS3,M1S3,M2S3,IMS3,
     + BM1,BM2,BS1,BS2,BS3
      REAL(KIND=RK),ALLOCATABLE::A(:),B(:),M1(:),M2(:),IM(:)
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::CTEMP
      integer(kind=8)  ::count1,count2,count_rate,count_max
      real timespend
      logical :: error_flag = .false.
      INTEGER::RETRY_COUNT,MAX_RETRIES=10,E
      
      
      IF(LABEL>3)THEN
      WRITE(CTEMP,'(I6)')N
      FILENAME='.\M\M'//TRIM(ADJUSTL(CTEMP))//'.DAT'
      LENGTH=0
      OPEN(N,FILE=FILENAME)
      WRITE(N,*)LABEL
      CLOSE(N)
      RETURN
      ENDIF
      
      NLABEL=0!0 represents elastic
      DO I=2,L-1
          IF(NLABEL==1)EXIT
          DO J=2,L-1
              IF(NLABEL==1)EXIT
              DO K=2,L-1
                  IF(NP(I,J,K)/=0)THEN
                      NLABEL=1
                      EXIT
                  ENDIF
              ENDDO
          ENDDO
      ENDDO
      
                  
      
      
      
      
      !call system_clock(count1,count_rate)
      CALL DIA_ISM(3*L**3,IMN,IMM1,IMM2,IMS1,IMS2,IMS3,IM)  
      
      CALL getM(L,L,L,vp,vs,np,ns,DX,DY,DZ,DT,M1N,M1M1,M1M2,
     + M1S1,M1S2,M1S3,M1,M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2,NLABEL)

      IF(NLABEL==0)THEN
          
      A=IM
      AN=IMN
      AM1=IMM1
      AM2=IMM2
      AS1=IMS1
      AS2=IMS2
      AS3=IMS3
      
      IM=(-1)*IM
      
      M2=IM
      M2N=IMN
      M2M1=IMM1
      M2M2=IMM2
      M2S1=IMS1
      M2S2=IMS2
      M2S3=IMS3
          
      IM=(-2)*IM
      CALL DIA_ADDSM(M1N,M1M1,M1M2,M1S1,M1S2,M1S3,M1,IMN,IMM1,
     + IMM2,IMS1,IMS2,IMS3,IM)
          
      ELSE
          
      A=M2
      AN=M2N
      AM1=M2M1
      AM2=M2M2
      AS1=M2S1
      AS2=M2S2
      AS3=M2S3
          
      CALL DIA_ADDSM(AN,AM1,AM2,AS1,AS2,AS3,A,IMN,IMM1,
     + IMM2,IMS1,IMS2,IMS3,IM)
      IM=(-1)*IM
      CALL DIA_ADDSM(M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2,IMN,IMM1,
     + IMM2,IMS1,IMS2,IMS3,IM)
      IM=(-2)*IM
      CALL DIA_ADDSM(M1N,M1M1,M1M2,M1S1,M1S2,M1S3,M1,IMN,IMM1,
     + IMM2,IMS1,IMS2,IMS3,IM)
      
      ENDIF
      
      
      IF(LABEL==2.OR.LABEL==3)THEN

          CALL PDC_A(N,L,L,L,AN,AM1,AM2,AS1,AS2,AS3,A,M1N,M1M1,M1M2,
     +     M1S1,M1S2,M1S3,M1,M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2)
            
      ENDIF
      
      WRITE(CTEMP,'(I6)')N
      FILENAME='.\M\M'//TRIM(ADJUSTL(CTEMP))//'.DAT'
      OPEN(N,FILE=FILENAME)
      WRITE(N,*)LABEL
      WRITE(N,*)NLABEL
      WRITE(N,*)AN
      WRITE(N,*)AM1
      WRITE(N,*)AM2
      WRITE(N,*)AS1
      WRITE(N,*)AS3
      WRITE(N,*)A
      WRITE(N,*)M1N
      WRITE(N,*)M1M1
      WRITE(N,*)M1M2
      WRITE(N,*)M1S1
      WRITE(N,*)M1S3
      WRITE(N,*)M1
      WRITE(N,*)M2N
      WRITE(N,*)M2M1
      WRITE(N,*)M2M2
      WRITE(N,*)M2S1
      WRITE(N,*)M2S3
      WRITE(N,*)M2
      CLOSE(N)
      
      END
        
          
          
      !Special points: initial point, physical boundary, PML boundary
      SUBROUTINE MESH_B(N,U,L,T,INI_T,error_flag,e)
      ! Performs subsequent matrix operations to compute next timestep iteration solution
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::NUM,N,L,LENGTH(3),LABEL,T,R_T,INI_T,IOS,e,S,NLABEL
      REAL:: DX,DY,DZ,DT
      REAL,DIMENSION(L,L,L)::VP,VS,NP,NS
      REAL(KIND=RK)::U(3,3*L**3),R
      REAL(KIND=RK),DIMENSION(3*L**3)::B,B1,B2
      INTEGER,DIMENSION(2)::AN,M1N,M2N,IAN
      INTEGER,ALLOCATABLE,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3,IAM1,IAM2,
     + IAS1,IAS2,IAS3,M1M1,M1M2,M1S1,M1S2,M1S3,M2M1,M2M2,M2S1,M2S2,M2S3
      REAL(KIND=RK),ALLOCATABLE::A(:),IA(:),M1(:),M2(:)
      

      REAL(KIND=RK),ALLOCATABLE::ini_U(:),ini_V(:),ini_W(:)
      INTEGER,ALLOCATABLE::X(:,:) 
      REAL::RHO
      real(kind=rk),allocatable::LU(:)
      INTEGER,allocatable::LT(:),LX(:)
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::CTEMP,CT
      logical, intent(inout) :: error_flag
      
      

      IOS=0
      
      WRITE(CTEMP,'(I6)')N
      WRITE(CT,'(I6)')T
      FILENAME='.\M\M'//TRIM(ADJUSTL(CTEMP))//'.DAT'

      OPEN(N,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS)   
      READ(N,*,IOSTAT=IOS)LABEL
      
      IF(LABEL>3)THEN
          CLOSE(N)
          RETURN
      ENDIF
      READ(N,*,IOSTAT=IOS)NLABEL
      !load A,M1,M2
      READ(N,*,IOSTAT=IOS)AN
      ALLOCATE(AM1(AN(2)),AM2(AN(2)))
      READ(N,*,IOSTAT=IOS)AM1
      READ(N,*,IOSTAT=IOS)AM2
      S=0
      DO I=1,AN(2)
          S=S+AM2(I)
      ENDDO
      ALLOCATE(AS1(S),AS2(S),AS3(S))
      READ(N,*,IOSTAT=IOS)AS1
      READ(N,*,IOSTAT=IOS)AS3
      AS2(1)=1
      DO I=2,S
          AS2(I)=AS2(I-1)+AS3(I-1)
      ENDDO
      S=AS2(S)+AS3(S)-1
      ALLOCATE(A(S))
      READ(N,*,IOSTAT=IOS)A
      
      READ(N,*,IOSTAT=IOS)M1N
      ALLOCATE(M1M1(M1N(2)),M1M2(M1N(2)))
      READ(N,*,IOSTAT=IOS)M1M1
      READ(N,*,IOSTAT=IOS)M1M2
      S=0
      DO I=1,M1N(2)
          S=S+M1M2(I)
      ENDDO
      ALLOCATE(M1S1(S),M1S2(S),M1S3(S))
      READ(N,*,IOSTAT=IOS)M1S1
      READ(N,*,IOSTAT=IOS)M1S3
      M1S2(1)=1
      DO I=2,S
          M1S2(I)=M1S2(I-1)+M1S3(I-1)
      ENDDO
      S=M1S2(S)+M1S3(S)-1
      ALLOCATE(M1(S))
      READ(N,*,IOSTAT=IOS)M1
      
      
      READ(N,*,IOSTAT=IOS)M2N
      ALLOCATE(M2M1(M2N(2)),M2M2(M2N(2)))
      READ(N,*,IOSTAT=IOS)M2M1
      READ(N,*,IOSTAT=IOS)M2M2
      S=0
      DO I=1,M2N(2)
          S=S+M2M2(I)
      ENDDO
      ALLOCATE(M2S1(S),M2S2(S),M2S3(S))
      READ(N,*,IOSTAT=IOS)M2S1
      READ(N,*,IOSTAT=IOS)M2S3
      M2S2(1)=1
      DO I=2,S
          M2S2(I)=M2S2(I-1)+M2S3(I-1)
      ENDDO
      S=M2S2(S)+M2S3(S)-1
      ALLOCATE(M2(S))
      READ(N,*,IOSTAT=IOS)M2

      CLOSE(N)
      
      
      
      
      
      CALL DOTSM(M1N,M1M1,M1M2,M1S1,M1S2,M1S3,M1,U(2,:),B1,3*L**3)
      CALL DOTSM(M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2,U(1,:),B2,3*L**3)
      B=B1+B2

      ! Modify B
      ! LABEL ranges 0~7, binary representation: (PML boundary, physical boundary, initial value)
      IF(MOD(LABEL,2).EQ.1.AND.T<INI_T+1)THEN
          ! Initial value
          LABEL=LABEL/2
! Filename: inin.dat  
! File contents:  
!   ini_num - number of initial points  
!   ini_x   - coordinates of initial points  
!   u,v,w   - initial values
          
          FILENAME='.\INI\INI'//TRIM(ADJUSTL(CTEMP))//'_'
     +     //TRIM(ADJUSTL(CT))//'.DAT'
          OPEN(N,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS)

          READ(N,*,IOSTAT=IOS)NUM

      READ(N,*)DT,RHO
          ALLOCATE(X(NUM,3),ini_U(NUM),ini_V(NUM),ini_W(NUM))
          DO I=1,NUM
              READ(N,*,IOSTAT=IOS)X(I,:)
          ENDDO
          READ(N,*,IOSTAT=IOS)INI_U
          READ(N,*,IOSTAT=IOS)INI_V
          READ(N,*,IOSTAT=IOS)INI_W
          CLOSE(N)
          
          call initial_value_B(L,L,L,DT,RHO,num,X,ini_U,ini_V,ini_W,B)
          ENDIF
          
      IF(NLABEL==0)THEN
          U(3,:)=B
          
      ELSE
          

      U(3,:)=U(2,:)
      DO I=1,10
          CALL GMRES(AN,AM1,AM2,AS1,AS2,AS3,A,B,U(3,:),3*L**3,R)   
          IF(R<0)EXIT
      ENDDO
      
      !WRITE(*,*)N,'Number of Iterations:',I-1,R
      
      ENDIF

      END
      
      !!!!!!!!!
      !Here are several Subroutines for partitioned data processing:
! - Region-based data handling
! - Data merging operations
! - Input/output operations
      
      
      
      SUBROUTINE GET_X(NX,NY,NZ,OUTPUT,T)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,NX,NY,NZ,I,J,K,NUM
      REAL(KIND=RK)::VALUE   
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::CTEMP
      REAL(KIND=RK),DIMENSION(NX,NY,NZ,3)::OUTPUT
      WRITE(CTEMP,'(I6)')T
      FILENAME='.\OUTPUT\OUTPUT'//TRIM(ADJUSTL(CTEMP))//'.DAT'

      OPEN(N,FILE=FILENAME,status='old')
      READ(N,*)I,J,K
      
      DO NUM=1,3
          DO K=1,NZ
              DO J=1,NY
                  DO I=1,NX
                      READ(N,*)VALUE
                      OUTPUT(I,J,K,NUM)=VALUE
                  ENDDO
              ENDDO
          ENDDO
      ENDDO


      CLOSE(N)

      END
      
      SUBROUTINE GET_U(N,NX,NY,NZ,L,INPUT,U,LABEL)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::N,LABEL,T,NX,NY,NZ,L,X,Y,Z,I,J,K,NUM,ISZERO
      REAL(KIND=RK)::U(3*L**3),VALUE   
      REAL(KIND=RK),DIMENSION(NX,NY,NZ,3)::INPUT
      
      CALL N_TO_X0(N,NX,NY,NZ,L,X,Y,Z)

      ISZERO=1
      DO NUM=1,3
          DO K=1,L
              DO J=1,L
                  DO I=1,L
          U(I+(J-1)*L+(K-1)*L**2+(NUM-1)*L**3)=
     +     INPUT(X+I-1,Y+J-1,Z+K-1,NUM)
          IF(ABS(INPUT(X+I-1,Y+J-1,Z+K-1,NUM))>UMIN)ISZERO=0
                  ENDDO
              ENDDO
          ENDDO
      ENDDO

      IF(ISZERO.EQ.1)LABEL=LABEL-1
      IF(ISZERO.EQ.0)LABEL=1
      
      END
      
      SUBROUTINE ASSIGN_VALUE(N,U,NX,NY,NZ,L,OUTPUT)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::N,NX,NY,NZ,L,LX,LY,LZ,X,Y,Z,I,J,K,NUM,M
      REAL(KIND=RK)::U(3*L**3),OUTPUT(NX,NY,NZ,3)   
      
      CALL N_TO_X0(N,NX,NY,NZ,L,X,Y,Z)

      DO NUM=1,3
      DO K=2,L-1
          DO J=2,L-1
              DO I=2,L-1
                  M=0
                  IF(I==2.OR.I==L-1)M=M+1
                  IF(J==2.OR.J==L-1)M=M+1
                  IF(K==2.OR.K==L-1)M=M+1          
      OUTPUT(X+I-1,Y+J-1,Z+K-1,NUM)=OUTPUT(X+I-1,Y+J-1,Z+K-1,NUM)
     + +U(I+(J-1)*L+(K-1)*L**2+(NUM-1)*L**3)/(2**M)

              ENDDO
          ENDDO
      ENDDO
      ENDDO


      
      END
      
      SUBROUTINE OUTPUT_X_T(T,NX,NY,NZ,OUTPUT)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,NX,NY,NZ,I,J,K,NUM,N   
      REAL(KIND=RK),DIMENSION(:,:,:,:),allocatable::OUTPUT
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::C1,C2
      
      WRITE(C1,'(I6)')T
      
      DO J=1,NY
      WRITE(C2,'(I6)')J
      FILENAME='.\Y\Y'//TRIM(ADJUSTL(C1))//'_'
     + //TRIM(ADJUSTL(C2))//'.DAT'
      OPEN(J,FILE=FILENAME)
      ENDDO
      
      FILENAME='.\OUTPUT\OUTPUT'//TRIM(ADJUSTL(C1))//'.DAT'
      OPEN(NY+1,FILE=FILENAME)
      WRITE(NY+1,*)NX,NY,NZ
      
      DO NUM=1,3
          DO K=1,NZ
              DO J=1,NY
                  DO I=1,NX
                      WRITE(NY+1,*)OUTPUT(I,J,K,NUM)
                      WRITE(J,*)OUTPUT(I,J,K,NUM)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
     
      DO J=1,NY+1
      CLOSE(J)
      ENDDO

      END
      
      
      SUBROUTINE INPUT_Y(NX,NZ,X1,X2,Y,T,IS_PML)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,NX,NZ,I,J,K,Y,NUM,IS_PML,ISZERO(2)
      REAL(KIND=RK)::X1(NX,NZ,3),X2(NX,NZ,3),VALUE1,VALUE2  
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::C1,C2
      
      
      IF(IS_PML==0)THEN
          
      WRITE(C1,'(I6)')T-2
      WRITE(C2,'(I6)')Y
      FILENAME='.\Y\Y'//TRIM(ADJUSTL(C1))//'_'
     + //TRIM(ADJUSTL(C2))//'.DAT'
      OPEN(1,FILE=FILENAME,STATUS='OLD')
      WRITE(C1,'(I6)')T-1
      FILENAME='.\Y\Y'//TRIM(ADJUSTL(C1))//'_'
     + //TRIM(ADJUSTL(C2))//'.DAT'
      OPEN(2,FILE=FILENAME,STATUS='OLD')
     
      DO NUM=1,3
          DO K=1,NZ
              DO I=1,NX
                  READ(1,*)X1(I,K,NUM)
                  READ(2,*)X2(I,K,NUM)
              ENDDO
          ENDDO
      ENDDO
      CLOSE(1,STATUS='DELETE')!Files for time t-2 can be deleted after reading
      CLOSE(2)
      
      ELSE
          
      WRITE(C1,'(I6)')T-2
      WRITE(C2,'(I6)')Y
      FILENAME='.\PML\PML'//TRIM(ADJUSTL(C2))//'_'
     + //TRIM(ADJUSTL(C1))//'.DAT'!PML+Y_T
      OPEN(1,FILE=FILENAME,STATUS='OLD')
      WRITE(C1,'(I6)')T-1
      FILENAME='.\PML\PML'//TRIM(ADJUSTL(C2))//'_'
     + //TRIM(ADJUSTL(C1))//'.DAT'!PML+Y_T
      OPEN(2,FILE=FILENAME,STATUS='OLD')
      
      READ(1,*)ISZERO(1)
      READ(2,*)ISZERO(2)
      IF(ISZERO(1)==1)X1=0
      IF(ISZERO(2)==1)X2=0
      IF(ISZERO(1)+ISZERO(2)/=2)THEN
      DO K=1,NZ
          DO I=1,NX
              DO NUM=1,12
                  IF(ISZERO(1)/=1)THEN
                      READ(1,*)VALUE1
                      IF(NUM>9)X1(I,K,NUM-9)=VALUE1
                  ENDIF
                  IF(ISZERO(2)/=1)THEN
                      READ(2,*)VALUE2
                      IF(NUM>9)X2(I,K,NUM-9)=VALUE2
                  ENDIF
              ENDDO
          ENDDO
      ENDDO
      ENDIF
      CLOSE(1)
      CLOSE(2)
      
      WRITE(C1,'(I6)')T-2
      FILENAME='.\Y\Y'//TRIM(ADJUSTL(C1))//'_'
     + //TRIM(ADJUSTL(C2))//'.DAT'
      OPEN(1,FILE=FILENAME,STATUS='OLD')
      
      CLOSE(1,STATUS='DELETE')!Files for time t-2 can be deleted after reading
      ENDIF
      END
      
      SUBROUTINE OUTPUT_Y(NX,NZ,X,N,T)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,NX,NY,NZ,I,J,K,N
      REAL(KIND=RK)::X(NX,NZ,3)   
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::C1,C2
      
      WRITE(C1,'(I6)')T
      WRITE(C2,'(I6)')N
      FILENAME='.\Y\Y'//TRIM(ADJUSTL(C1))//'_'
     + //TRIM(ADJUSTL(C2))//'.DAT'

      OPEN(T,FILE=FILENAME)
      DO NUM=1,3
          DO K=1,NZ
              DO I=1,NX
                  WRITE(T,*)X(I,K,NUM)
              ENDDO
          ENDDO
      ENDDO
      CLOSE(T)
      
      END
      
      SUBROUTINE Y_TO_OUTPUT(T,NX,NY,NZ)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,NX,NY,NZ,I,J,K,NUM
      REAL(KIND=RK)::VALUE 
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::C1,C2
      
      WRITE(C1,'(I6)')T
      
      DO J=1,NY
      WRITE(C2,'(I6)')J
      FILENAME='.\Y\Y'//TRIM(ADJUSTL(C1))//'_'
     + //TRIM(ADJUSTL(C2))//'.DAT'
      OPEN(J,FILE=FILENAME,STATUS='OLD')
      ENDDO
      
      FILENAME='.\OUTPUT\OUTPUT'//TRIM(ADJUSTL(C1))//'.DAT'
      OPEN(NY+1,FILE=FILENAME)
      WRITE(NY+1,*)NX,NY,NZ

      DO NUM=1,3
          DO K=1,NZ
              DO J=1,NY
                  DO I=1,NX
                  READ(J,*)VALUE
                  IF(ABS(VALUE)<UMIN)VALUE=0
                  WRITE(NY+1,*)VALUE
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      DO J=1,NY+1
      CLOSE(J)
      ENDDO
      
      END
      
      
      SUBROUTINE PML_TO_OUTPUT(T,NX,NY,NZ)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,NX,NY,NZ,I,J,K,NUM,ISZERO
      REAL(KIND=RK)::VALUE 
      REAL(KIND=RK),DIMENSION(:,:,:,:),allocatable::U
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::C1,C2
      
      
      
      ALLOCATE(U(NX,NY,NZ,3))
      WRITE(C1,'(I6)')T
      !$OMP PARALLEL PRIVATE(J,C2,FILENAME,ISZERO,I,K,NUM)
      !$omp do SCHEDULE(DYNAMIC)
      DO J=1,NY
      WRITE(C2,'(I6)')J
      FILENAME='.\PML\PML'//TRIM(ADJUSTL(C2))//'_'
     + //TRIM(ADJUSTL(C1))//'.DAT'
      OPEN(J,FILE=FILENAME,STATUS='OLD')
      READ(J,*)ISZERO
      IF(ISZERO==1)THEN
          U(:,J,:,:)=0
      ELSE
      DO K=1,NZ
          DO I=1,NX
              DO NUM=1,12
                  IF(NUM<10)THEN
                      READ(J,*)VALUE
                  ELSE
                      READ(J,*)U(I,J,K,NUM-9)
                  ENDIF
              ENDDO
          ENDDO
      ENDDO
      ENDIF
      
      CLOSE(J)
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
      FILENAME='.\OUTPUT\OUTPUT'//TRIM(ADJUSTL(C1))//'.DAT'
      OPEN(NY+1,FILE=FILENAME)
      WRITE(NY+1,*)NX,NY,NZ

      DO NUM=1,3
          DO K=1,NZ
              DO J=1,NY
                  DO I=1,NX
                      IF(ieee_is_nan(U(I,J,K,NUM)))then
                          print*,'OUTPUT IS NaN'
                          PAUSE
                      ENDIF
                  WRITE(NY+1,*)U(I,J,K,NUM)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      CLOSE(NY+1)

      END
      
      
      SUBROUTINE DIVIDE_Y_Z(INPUT,OUTPUT,NX,NY,NZ,L,Z0)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::NX,NY,NZ,L,I,J,K,NUM,Z0
      REAL(KIND=RK)::INPUT(NX,NY,NZ,3),OUTPUT(NX,NY,L,3)
      
      DO I=1,NX
          DO J=1,NY
              DO K=1,L
                  DO NUM=1,3
                      OUTPUT(I,J,K,NUM)=INPUT(I,J,(L-3)*(Z0-1)+K,NUM)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      
      END

      SUBROUTINE UNION_Z_Y(OUTPUT,INPUT,NX,NY,NZ,L,Z0)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::NX,NY,NZ,L,I,J,K,NUM,Z0
      REAL(KIND=RK)::INPUT(NX,NY,L,3),OUTPUT(NX,NY,NZ,3)
      
      DO I=2,NX-1
          DO J=2,NY-1
              DO K=2,L-1
                  DO NUM=1,3
                      OUTPUT(I,J,(L-3)*(Z0-1)+K,NUM)=
     +OUTPUT(I,J,(L-3)*(Z0-1)+K,NUM)+INPUT(I,J,K,NUM)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      
      END
      
      
      
      end module cube_of_five_meshes