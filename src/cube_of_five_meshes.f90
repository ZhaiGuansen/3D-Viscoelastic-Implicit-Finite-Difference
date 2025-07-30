      module cube_of_five_meshes
      !USE sparse_matrix
      USE DIA_SM
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
      
      SUBROUTINE MESH_A(N,L,DX,DY,DZ,DT,FP,VP,VS,NP,NS,LABEL,M_MAP)
      ! Pre-loading coefficient matrices
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::N,I,J,K,LABEL,NLABEL
      INTEGER,DIMENSION(3)::PMLD,D,DL
      REAL:: DX,DY,DZ,DT,FP
      REAL,DIMENSION(L,L,L)::VP,VS,NP,NS
      INTEGER,DIMENSION(2)::AN,M1N,M2N,IMN,PM1N,PM2N
      INTEGER,ALLOCATABLE,DIMENSION(:)::AM1,M1M1,M2M1,IMM1,AM2,M1M2,M2M2 &
          ,IMM2,AS1,M1S1,M2S1,IMS1,AS2,M1S2,M2S2,IMS2,AS3,M1S3,M2S3,IMS3,&
          PM1I,PM1J,PM2I,PM2J,PM1M1,PM1M2,PM1S1,PM1S2,PM1S3,PM2M1,PM2M2,PM2S1,PM2S2,PM2S3,M_MAP
      REAL(KIND=RK),ALLOCATABLE,DIMENSION(:)::A,M1,M2,IM,PM13,PM23,PM1,PM2
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::CTEMP
      CHARACTER(LEN=30)::BASE
      CHARACTER(LEN=10)::C1
      
      
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
      
      CALL getM(L,L,L,vp,vs,np,ns,DX,DY,DZ,DT,M1N,M1M1,M1M2,&
          M1S1,M1S2,M1S3,M1,M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2,NLABEL)
      
  
      IF(LABEL>3)THEN
      NLABEL=1
      WRITE(C1,'(I6)')N
      BASE='./BASE/BASE'//TRIM(ADJUSTL(C1))//'.DAT'
      OPEN(N+100,FILE=BASE,STATUS='OLD',IOSTAT=IOS)
          READ(N+100,*,IOSTAT=IOS)PMLD(1),PMLD(2),PMLD(3)
          READ(N+100,*,IOSTAT=IOS)D(1),D(2),D(3)
      CLOSE(N+100,IOSTAT=IOS)

      DO K=2,L-1
          DO J=2,L-1
              DO I=2,L-1
                  IF(D(1)<1-I)THEN
                      DL(1)=ABS(D(1)+I-1)
                  ELSEIF(D(1)>L-I)THEN
                      DL(1)=D(1)-L+I
                  ELSE
                      DL(1)=0
                  ENDIF
                  IF(D(2)<1-J)THEN
                      DL(2)=ABS(D(2)+J-1)
                  ELSEIF(D(2)>L-J)THEN
                      DL(2)=D(2)-L+J
                  ELSE
                      DL(2)=0
                  ENDIF
                  IF(D(3)<1-K)THEN
                      DL(3)=ABS(D(3)+K-1)
                  ELSEIF(D(3)>L-K)THEN
                      DL(3)=D(3)-L+K
                  ELSE
                      DL(3)=0
                  ENDIF
                  IF(DL(1)+DL(2)+DL(3)==0)CYCLE
                  !m1=m1+pm1,m2=m2+pm2
                  
      CALL getPM(L,L,L,I,J,K,vp(I,J,K),vs(I,J,K),np(I,J,K),ns(I,J,K),&
          DX,DY,DZ,DT,FP,PMLD,DL,PM1I,PM1J,PM13,PM2I,PM2J,PM23)
       

                  
              ENDDO
          ENDDO
      ENDDO
      
      !3toDIA
          CALL SORT(PM13,PM1I,PM1J)
          CALL SORT(PM23,PM2I,PM2J)
      
          CALL DIAOF3(PM1I,PM1J,PM13,3*L**3,PM1N,PM1M1,PM1M2,PM1S1,PM1S2,PM1S3,PM1)
          CALL DIAOF3(PM2I,PM2J,PM23,3*L**3,PM2N,PM2M1,PM2M2,PM2S1,PM2S2,PM2S3,PM2)
      
          call DIA_ADDSM(M1N,M1M1,M1M2,M1S1,M1S2,M1S3,M1,PM1N,PM1M1,PM1M2,PM1S1,PM1S2,PM1S3,PM1)
          call DIA_ADDSM(M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2,PM2N,PM2M1,PM2M2,PM2S1,PM2S2,PM2S3,PM2)  

      ENDIF
      
      
      IF(M_MAP(N)==N)THEN
      
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
      CALL DIA_ADDSM(M1N,M1M1,M1M2,M1S1,M1S2,M1S3,M1,IMN,IMM1,IMM2,IMS1,IMS2,IMS3,IM)
          
      ELSE
          
      A=M2
      AN=M2N
      AM1=M2M1
      AM2=M2M2
      AS1=M2S1
      AS2=M2S2
      AS3=M2S3
          
      CALL DIA_ADDSM(AN,AM1,AM2,AS1,AS2,AS3,A,IMN,IMM1,IMM2,IMS1,IMS2,IMS3,IM)
      IM=(-1)*IM
      CALL DIA_ADDSM(M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2,IMN,IMM1,IMM2,IMS1,IMS2,IMS3,IM)
      IM=(-2)*IM
      CALL DIA_ADDSM(M1N,M1M1,M1M2,M1S1,M1S2,M1S3,M1,IMN,IMM1,IMM2,IMS1,IMS2,IMS3,IM)
      
      ENDIF
      
      IF(LABEL==2.OR.LABEL==3)THEN

          CALL PDC_A(N,L,L,L,AN,AM1,AM2,AS1,AS2,AS3,A,M1N,M1M1,M1M2,M1S1,M1S2,M1S3,M1,M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2)
            
      ENDIF
      
      

      WRITE(CTEMP,'(I6)')N
      FILENAME='./M/M'//TRIM(ADJUSTL(CTEMP))//'.DAT'
      OPEN(N+100,FILE=FILENAME)
      WRITE(N+100,*)LABEL
      WRITE(N+100,*)M_MAP(N),NLABEL
      WRITE(N+100,*)AN
      WRITE(N+100,*)AM1
      WRITE(N+100,*)AM2
      WRITE(N+100,*)AS1
      WRITE(N+100,*)AS3
      WRITE(N+100,*)A
      WRITE(N+100,*)M1N
      WRITE(N+100,*)M1M1
      WRITE(N+100,*)M1M2
      WRITE(N+100,*)M1S1
      WRITE(N+100,*)M1S3
      WRITE(N+100,*)M1
      WRITE(N+100,*)M2N
      WRITE(N+100,*)M2M1
      WRITE(N+100,*)M2M2
      WRITE(N+100,*)M2S1
      WRITE(N+100,*)M2S3
      WRITE(N+100,*)M2
      CLOSE(N+100)


      ELSE
                WRITE(CTEMP,'(I6)')N
          FILENAME='./M/M'//TRIM(ADJUSTL(CTEMP))//'.DAT'
          OPEN(N+100,FILE=FILENAME)
          WRITE(N+100,*)LABEL
          WRITE(N+100,*)M_MAP(N),NLABEL
          CLOSE(N+100)
      ENDIF
      
      
      END
        
          
          
      !Special points: initial point, physical boundary, PML boundary
      SUBROUTINE MESH_B(N,U,L,T,INI_T)
      ! Performs subsequent matrix operations to compute next timestep iteration solution
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::NUM,N,N0,L,LABEL,T,INI_T,IOS,S,NLABEL
      REAL:: DT
      REAL(KIND=RK)::U(3,3*L**3),R
      REAL(KIND=RK),DIMENSION(3*L**3)::B,B1,B2
      INTEGER,DIMENSION(2)::AN,M1N,M2N
      INTEGER,ALLOCATABLE,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3,&
          M1M1,M1M2,M1S1,M1S2,M1S3,M2M1,M2M2,M2S1,M2S2,M2S3
      REAL(KIND=RK),ALLOCATABLE::A(:),M1(:),M2(:)
      REAL(KIND=RK),ALLOCATABLE::ini_U(:),ini_V(:),ini_W(:)
      INTEGER,ALLOCATABLE::X(:,:) 
      REAL::RHO
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::CTEMP,CT
      
      

      IOS=0
      
      WRITE(CTEMP,'(I6)')N
      WRITE(CT,'(I6)')T
      FILENAME='./M/M'//TRIM(ADJUSTL(CTEMP))//'.DAT'

      OPEN(N+100,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS)   
      READ(N+100,*,IOSTAT=IOS)LABEL

      READ(N+100,*,IOSTAT=IOS)N0,NLABEL
      IF(N/=N0)THEN
          CLOSE(N+100)
          WRITE(CTEMP,'(I6)')N0
          WRITE(CT,'(I6)')T
          FILENAME='./M/M'//TRIM(ADJUSTL(CTEMP))//'.DAT'

          OPEN(N+100,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS)   
          READ(N+100,*,IOSTAT=IOS)LABEL

          READ(N+100,*,IOSTAT=IOS)N0,NLABEL
      ENDIF

      !load A,M1,M2
      READ(N+100,*,IOSTAT=IOS)AN
      ALLOCATE(AM1(AN(2)),AM2(AN(2)))
      READ(N+100,*,IOSTAT=IOS)AM1
      READ(N+100,*,IOSTAT=IOS)AM2
      S=0
      DO I=1,AN(2)
          S=S+AM2(I)
      ENDDO
      ALLOCATE(AS1(S),AS2(S),AS3(S))
      READ(N+100,*,IOSTAT=IOS)AS1
      READ(N+100,*,IOSTAT=IOS)AS3
      AS2(1)=1
      DO I=2,S
          AS2(I)=AS2(I-1)+AS3(I-1)
      ENDDO
      S=AS2(S)+AS3(S)-1
      ALLOCATE(A(S))
      READ(N+100,*,IOSTAT=IOS)A
      
      READ(N+100,*,IOSTAT=IOS)M1N
      ALLOCATE(M1M1(M1N(2)),M1M2(M1N(2)))
      READ(N+100,*,IOSTAT=IOS)M1M1
      READ(N+100,*,IOSTAT=IOS)M1M2
      S=0
      DO I=1,M1N(2)
          S=S+M1M2(I)
      ENDDO
      ALLOCATE(M1S1(S),M1S2(S),M1S3(S))
      READ(N+100,*,IOSTAT=IOS)M1S1
      READ(N+100,*,IOSTAT=IOS)M1S3
      M1S2(1)=1
      DO I=2,S
          M1S2(I)=M1S2(I-1)+M1S3(I-1)
      ENDDO
      S=M1S2(S)+M1S3(S)-1
      ALLOCATE(M1(S))
      READ(N+100,*,IOSTAT=IOS)M1

      READ(N+100,*,IOSTAT=IOS)M2N
      ALLOCATE(M2M1(M2N(2)),M2M2(M2N(2)))
      READ(N+100,*,IOSTAT=IOS)M2M1
      READ(N+100,*,IOSTAT=IOS)M2M2
      S=0
      DO I=1,M2N(2)
          S=S+M2M2(I)
      ENDDO
      ALLOCATE(M2S1(S),M2S2(S),M2S3(S))
      READ(N+100,*,IOSTAT=IOS)M2S1
      READ(N+100,*,IOSTAT=IOS)M2S3
      M2S2(1)=1
      DO I=2,S
          M2S2(I)=M2S2(I-1)+M2S3(I-1)
      ENDDO
      S=M2S2(S)+M2S3(S)-1
      ALLOCATE(M2(S))
      READ(N+100,*,IOSTAT=IOS)M2

      CLOSE(N+100)


      
      
      CALL DOTSM(M1N,M1M1,M1M2,M1S1,M1S2,M1S3,M1,U(2,:),B1,3*L**3)
      CALL DOTSM(M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2,U(1,:),B2,3*L**3)
      B=B1+B2

      ! Modify B
      ! LABEL ranges 0~7, binary representation: (PML boundary, physical boundary, initial value)
      IF(MOD(LABEL,2).EQ.1.AND.T<INI_T+1)THEN
          !Initial value
          LABEL=LABEL/2
! Filename: inin.dat  
! File contents:  
!   ini_num - number of initial points  
!   ini_x   - coordinates of initial points  
!   u,v,w   - initial values
          
          FILENAME='./INI/INI'//TRIM(ADJUSTL(CTEMP))//'_' //TRIM(ADJUSTL(CT))//'.DAT'
          OPEN(N+100,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS)

          READ(N+100,*,IOSTAT=IOS)NUM

      READ(N+100,*)DT,RHO
          ALLOCATE(X(NUM,3),ini_U(NUM),ini_V(NUM),ini_W(NUM))
          DO I=1,NUM
              READ(N+100,*,IOSTAT=IOS)X(I,:)
          ENDDO
          READ(N+100,*,IOSTAT=IOS)INI_U
          READ(N+100,*,IOSTAT=IOS)INI_V
          READ(N+100,*,IOSTAT=IOS)INI_W
          CLOSE(N+100)

          call initial_value_B(L,L,L,num,X,ini_U,ini_V,ini_W,B)
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
      
      SUBROUTINE GET_U(N,NX,NY,NZ,L,INPUT,U,LABEL)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::N,LABEL,NX,NY,NZ,L,X,Y,Z,I,J,K,NUM,ISZERO
      REAL(KIND=RK)::U(3*L**3)   
      REAL(KIND=RK),DIMENSION(NX,NY,NZ,3)::INPUT
      
      CALL N_TO_X0(N,NX,NY,NZ,L,X,Y,Z)

      ISZERO=1
      DO NUM=1,3
          DO K=1,L
              DO J=1,L
                  DO I=1,L
          U(I+(J-1)*L+(K-1)*L**2+(NUM-1)*L**3)=INPUT(X+I-1,Y+J-1,Z+K-1,NUM)
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
      INTEGER::N,NX,NY,NZ,L,X,Y,Z,I,J,K,NUM,M
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
      OUTPUT(X+I-1,Y+J-1,Z+K-1,NUM)=OUTPUT(X+I-1,Y+J-1,Z+K-1,NUM)+U(I+(J-1)*L+(K-1)*L**2+(NUM-1)*L**3)/(2**M)

              ENDDO
          ENDDO
      ENDDO
      ENDDO


      
      END
      
      SUBROUTINE OUTPUT_X_T(T,NX,NY,NZ,OUTPUT)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,NX,NY,NZ,I,J,K,NUM
      REAL(KIND=RK),DIMENSION(:,:,:,:),allocatable::OUTPUT
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::C1,C2
      
      WRITE(C1,'(I6)')T
      
      DO J=1,NY
      WRITE(C2,'(I6)')J
      FILENAME='./Y/Y'//TRIM(ADJUSTL(C1))//'_'//TRIM(ADJUSTL(C2))//'.DAT'
      OPEN(J+100,FILE=FILENAME)
      ENDDO
      
      FILENAME='./OUTPUT/OUTPUT'//TRIM(ADJUSTL(C1))//'.DAT'
      OPEN(NY+101,FILE=FILENAME)
      WRITE(NY+101,*)NX,NY,NZ
      
      DO NUM=1,3
          DO K=1,NZ
              DO J=1,NY
                  DO I=1,NX
                      WRITE(NY+101,*)OUTPUT(I,J,K,NUM)
                      WRITE(J+100,*)OUTPUT(I,J,K,NUM)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
     
      DO J=1,NY+1
      CLOSE(J+100)
      ENDDO

      END
      
      
      SUBROUTINE INPUT_Y(NX,NZ,X1,X2,Y,T)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,NX,NZ,I,K,Y,NUM
      REAL(KIND=RK)::X1(NX,NZ,3),X2(NX,NZ,3) 
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::C1,C2
          
      WRITE(C1,'(I6)')T-2
      WRITE(C2,'(I6)')Y
      FILENAME='./Y/Y'//TRIM(ADJUSTL(C1))//'_'//TRIM(ADJUSTL(C2))//'.DAT'
      OPEN(101,FILE=FILENAME,STATUS='OLD')
      WRITE(C1,'(I6)')T-1
      FILENAME='./Y/Y'//TRIM(ADJUSTL(C1))//'_'//TRIM(ADJUSTL(C2))//'.DAT'
      OPEN(102,FILE=FILENAME,STATUS='OLD')
     
      DO NUM=1,3
          DO K=1,NZ
              DO I=1,NX
                  READ(101,*)X1(I,K,NUM)
                  READ(102,*)X2(I,K,NUM)
              ENDDO
          ENDDO
      ENDDO
      CLOSE(101,STATUS='DELETE')!Files for time t-2 can be deleted after reading
      CLOSE(102)
      
      END
      
      SUBROUTINE OUTPUT_Y(NX,NZ,X,N,T)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,NX,NZ,I,K,N
      REAL(KIND=RK)::X(NX,NZ,3)   
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::C1,C2
      
      WRITE(C1,'(I6)')T
      WRITE(C2,'(I6)')N
      FILENAME='./Y/Y'//TRIM(ADJUSTL(C1))//'_'//TRIM(ADJUSTL(C2))//'.DAT'

      OPEN(T+100,FILE=FILENAME)
      DO NUM=1,3
          DO K=1,NZ
              DO I=1,NX
                  WRITE(T+100,*)X(I,K,NUM)
              ENDDO
          ENDDO
      ENDDO
      CLOSE(T+100)

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
      FILENAME='./Y/Y'//TRIM(ADJUSTL(C1))//'_' //TRIM(ADJUSTL(C2))//'.DAT'
      OPEN(J+100,FILE=FILENAME,STATUS='OLD')
      ENDDO
      
      FILENAME='./OUTPUT/OUTPUT'//TRIM(ADJUSTL(C1))//'.DAT'
      OPEN(NY+101,FILE=FILENAME)
      WRITE(NY+101,*)NX,NY,NZ

      DO NUM=1,3
          DO K=1,NZ
              DO J=1,NY
                  DO I=1,NX
                  READ(J+100,*)VALUE
                  IF(ABS(VALUE)<UMIN)VALUE=0
                  WRITE(NY+101,*)VALUE
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      DO J=1,NY+1
      CLOSE(J+100)
      ENDDO
      
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
                      OUTPUT(I,J,(L-3)*(Z0-1)+K,NUM)=OUTPUT(I,J,(L-3)*(Z0-1)+K,NUM)+INPUT(I,J,K,NUM)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      
      END
      
      
      
      end module cube_of_five_meshes