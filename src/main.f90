
!!!!Main program section
      PROGRAM main

      use omp_lib
      use matrix_construction
      USE DIA_SM
      use INI
      USE PML
      USE PDC
      USE cube_of_five_meshes
      USE make_cube
      USE PARTIAL_MAIN

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      
!!Basic variables
      INTEGER::NX=43,NY=43,NZ=43,T0,NT
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,INI_T=80,IS_INI=0
      !1: initialized, 0: uninitialized, -1: error in previous run, 2: M matrix calculated, 3:Rewrite specific M elements,-2: regenerate previous two time steps
      REAL:: XSTEP=1,YSTEP=1,ZSTEP=1
      REAL:: TSTEP=2E-4,FP=150,AMPL=1

      integer::x0(3)=[22,22,22]!Shot point
! PML boundary related variables
      INTEGER::PMLD(3)=[15,15,15]
      INTEGER::LX,N2,LY
      INTEGER::I,J,K,N,YN
      INTEGER::LABEL,L=11
      REAL(KIND=RK),DIMENSION(:,:),allocatable::U!DIMENSION(3,375)
      REAL(KIND=RK),DIMENSION(:,:,:,:),allocatable::X,X1,X2,XZ,XZ1,XZ2!DIMENSION(NX,NY,NZ,3)
!Time related variables
      integer  ::count1,count2,count_rate,count_max
      integer  ::count11,count12
      real timespend,timespend0
      logical :: error_flag = .false.
      INTEGER::RETRY_COUNT,E 

      !character(len=256) :: arg
      !call GET_COMMAND_ARGUMENT(1,arg) 
      !print *, 'Argument1: ', trim(arg)  
      !read(arg,*)T0
      !call GET_COMMAND_ARGUMENT(2,arg) 
      !print *, 'Argument2: ', trim(arg)  
      !read(arg,*)NT
      !call GET_COMMAND_ARGUMENT(3,arg) 
      !print *, 'Argument3: ', trim(arg)  
      !read(arg,*)IS_INI

      call omp_set_num_threads(18)
      T0=3
      NT=200
      IS_INI=0
      
      !Pad to integer block region
      LX=(NX-3-1)/(L-3)+1
      LY=(NY-3-1)/(L-3)+1
      LZ=(NZ-3-1)/(L-3)+1
      NX=LX*(L-3)+3
      NY=LY*(L-3)+3
      NZ=LZ*(L-3)+3

      N2=LX*LZ
      IF(T0>3.AND.IS_INI>-1)IS_INI=1
      IF(IS_INI==0.OR.IS_INI>1)THEN
      CALL main_a(NX,NY,NZ,L,INI_T,XSTEP,YSTEP,ZSTEP,TSTEP,FP,AMPL,X0,PMLD,IS_INI)
      ENDIF
      call system_clock(count1,count_rate)

      ALLOCATE(X(NX,L,NZ,3),X1(NX,L,NZ,3),X2(NX,L,NZ,3))
      ALLOCATE(XZ(NX,L,L,3),XZ1(NX,L,L,3),XZ2(NX,L,L,3))
      ALLOCATE(U(3,3*L**3))
      X=0
      X1=0
      X2=0
      T=T0
      timespend0=0
      RETRY_COUNT=0

      IF (T>3.OR.IS_INI==-1) CALL RETRIEVE_Y(NX,NY,NZ,T0-2)
      IF (IS_INI==-2)CALL RETRIEVE_Y(NX,NY,NZ,T0-1)

      DO T=T,NT

          WRITE(*,*)T 
          call flush()
          call system_clock(count11,count_rate,count_max)
          !Load initial three columns
          DO YN=1,3
              CALL INPUT_Y(NX,NZ,X1(:,YN,:,:),X2(:,YN,:,:),YN,T)
          ENDDO
          DO J=1,LY
          !Read columns 4 to L     
              DO YN=4,L
                  CALL INPUT_Y(NX,NZ,X1(:,YN,:,:),X2(:,YN,:,:),(L-3)*(J-1)+YN,T)
              ENDDO
              
              !Redivide in z-direction

              DO K=1,LZ
                  XZ=0
                  CALL DIVIDE_Y_Z(X1,XZ1,NX,L,NZ,L,K)
                  CALL DIVIDE_Y_Z(X2,XZ2,NX,L,NZ,L,K)
          !$OMP PARALLEL PRIVATE(I,U,N,RETRY_COUNT,error_flag,E,LABEL)REDUCTION(+:XZ)
          !$omp do SCHEDULE(DYNAMIC)  
              DO I=1,LX

                  U=0
                  !!RETRY_COUNT=0
                  LABEL=1
                  CALL GET_U(I,NX,L,L,L,XZ1,U(1,:),LABEL)
                  CALL GET_U(I,NX,L,L,L,XZ2,U(2,:),LABEL)

                  IF(LABEL<0)CYCLE          

                  N=I+(J-1)*LX+(K-1)*LX*LY
                  !WRITE(*,*)I
                  error_flag=.false.
                  
                  CALL MESH_B(N,U,L,T,INI_T)
                  
                  IF(error_flag)THEN
                  WRITE(*,*)'--------'
                  WRITE(*,*)e
                  WRITE(*,*)'An error occurred, retrying...'
                  call flush()
                  STOP
                  ENDIF

                  CALL ASSIGN_VALUE(I,U(3,:),NX,L,L,L,XZ)

              ENDDO
              !$OMP END DO
              !$OMP END PARALLEL
              CALL UNION_Z_Y(X,XZ,NX,L,NZ,L,K)
              
              ENDDO
              

              !Output x columns except last three
              DO YN=1,L-3
                  CALL OUTPUT_Y(NX,NZ,X(:,YN,:,:),(L-3)*(J-1)+YN,T)
              ENDDO
              !Retain last three columns for next cycle
              DO YN=1,3
                  X1(:,YN,:,:)=X1(:,L-3+YN,:,:)  
                  X2(:,YN,:,:)=X2(:,L-3+YN,:,:) 
                  X(:,YN,:,:)=X(:,L-3+YN,:,:) 
              ENDDO
              DO YN=4,L
                  X(:,YN,:,:)=0
              ENDDO


          ENDDO
          !Output final three columns
          DO YN=1,3
              CALL OUTPUT_Y(NX,NZ,X(:,YN,:,:),(L-3)*LY+YN,T)
          ENDDO
          
          CALL Y_TO_OUTPUT(T,NX,NY,NZ)
          
         call system_clock(count12,count_rate,count_max)
          timespend=(count12-count11)/real(count_rate)
          WRITE(*,*) 'Internal time consumed (in seconds): ',timespend
          call flush()
         



      ENDDO
          
          
      call system_clock(count2,count_rate,count_max)
      timespend=(count2-count1)/real(count_rate)/real(NT-T0+1)
      WRITE(*,*) 'Average time consumed (in seconds): ',timespend 
      call flush()
      WRITE(*,*)'finish'
      call flush()
      STOP
      END PROGRAM main

      
      