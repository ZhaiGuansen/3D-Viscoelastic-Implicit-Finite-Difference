
!!!!主程序部分
      PROGRAM main
      USE LAPACK95
      use matrix_construction
      USE DIA_SM
      use INI
      USE PML
      USE PDC
      USE cube_of_five_meshes
      USE make_cube
      use omp_lib
      USE PARTIAL_MAIN
c	

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      
C基础变量
      INTEGER::NX=131,NY=123,NZ=67,T0,NT
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,ini_num=1,INI_T=160,IS_INI=0
      !1表示已经初始化,0表示未初始化,-1表示先前运行出现报错,2:M矩阵已计算-2:前两时刻重新生成
      REAL:: XSTEP=1,YSTEP=1,ZSTEP=0.5
      REAL:: TSTEP=1E-4,FP=150,AMPL=1

C初值相关变量
      integer::x0(3)=[30,62,34]!炮点
C PML边界相关变量
      INTEGER::PMLD(3)=[20,20,20]
      INTEGER::LX,N2,LY,PX,PY,PZ
      INTEGER::I,J,K,NUM,N,ISJ,IS_PML,YN
      INTEGER::LABEL,L=11
      REAL(KIND=RK),DIMENSION(:,:),allocatable::U!DIMENSION(3,375)
      REAL(KIND=RK),DIMENSION(:,:,:,:),allocatable::X,X1,X2,XZ,XZ1,XZ2!DIMENSION(NX,NY,NZ,3)
C时间计算相关变量
      integer(kind=8)  ::count1,count2,count_rate,count_max
      integer(kind=8)  ::count11,count12,count13
      real timespend,timespend0
      logical :: error_flag = .false.
      INTEGER::RETRY_COUNT,MAX_RETRIES=10,E
      integer :: arg_count, arg_num, iostat  
      character(len=256) :: arg
      


      call omp_set_num_threads(18)
      

      call GET_COMMAND_ARGUMENT(1,arg) 
      print *, 'Argument1: ', trim(arg)  
      read(arg,*)T0
      call GET_COMMAND_ARGUMENT(2,arg) 
      print *, 'Argument2: ', trim(arg)  
      read(arg,*)NT
      call GET_COMMAND_ARGUMENT(3,arg) 
      print *, 'Argument3: ', trim(arg)  
      read(arg,*)IS_INI
      
      !T0=29
      !NT=300
      !IS_INI=-1
      
      
      LX=(NX-3)/(L-3)
      LY=(NY-3)/(L-3)
      LZ=(NZ-3)/(L-3)
      N2=LX*LZ
      IF(T0>3.AND.IS_INI>-1)IS_INI=1
      IF(IS_INI==0.OR.IS_INI==2)THEN
      CALL main_a(NX,NY,NZ,L,INI_T,T0,NT,XSTEP,YSTEP,ZSTEP,TSTEP,
     + FP,AMPL,X0,PMLD,IS_INI)
      ENDIF


      !!T=T0
      !!RETRY_COUNT=0
      !!DO WHILE(T<NT+1)
      !!    CALL main_b(NX,NY,NZ,INI_T,T,NT,RETRY_COUNT)
      !!ENDDO
      
      call system_clock(count1,count_rate)


      ALLOCATE(X(NX,L,NZ,3),X1(NX,L,NZ,3),X2(NX,L,NZ,3))
      ALLOCATE(XZ(NX,L,L,3),XZ1(NX,L,L,3),XZ2(NX,L,L,3))
      ALLOCATE(U(3,3*L**3))
      X=0
      X1=0
      X2=0
      T=T0
      PX=PMLD(1)/(L-3)
      PY=PMLD(2)/(L-3)
      PZ=PMLD(3)/(L-3)
      timespend0=0
      IS_PML=0
      IF(PX<1)PX=1
      IF(PY<1)PY=1
      IF(PZ<1)PZ=1
      RETRY_COUNT=0


      IF (T>3.OR.IS_INI==-1) CALL RETRIEVE_Y(NX,NY,NZ,T0-2)
      IF (IS_INI==-2)CALL RETRIEVE_Y(NX,NY,NZ,T0-1)

      DO T=T,NT
          !PRINT*,'------------'
          PRINT*,T   
          
          call system_clock(count11,count_rate,count_max)

          DO YN=1,3
              CALL INPUT_Y(NX,NZ,X1(:,YN,:,:),X2(:,YN,:,:),
     +         (L-3)*(PY-1)+YN,T,IS_PML)
          ENDDO
          DO J=PY,LY+1-PY
  
              DO YN=4,L
                  CALL INPUT_Y(NX,NZ,X1(:,YN,:,:),X2(:,YN,:,:),
     +              (L-3)*(J-1)+YN,T,IS_PML)
              ENDDO


              DO K=PZ,LZ+1-PZ
                  XZ=0
                  CALL DIVIDE_Y_Z(X1,XZ1,NX,L,NZ,L,K)
                  CALL DIVIDE_Y_Z(X2,XZ2,NX,L,NZ,L,K)
          !$OMP PARALLEL PRIVATE(I,U,N,RETRY_COUNT,error_flag,E,LABEL)
     +       REDUCTION(+:XZ)
          !$omp do SCHEDULE(DYNAMIC)  
              DO I=PX,LX+1-PX

                  U=0

                  LABEL=1
                  CALL GET_U(I,NX,L,L,L,XZ1,U(1,:),LABEL)
                  CALL GET_U(I,NX,L,L,L,XZ2,U(2,:),LABEL)

                  IF(LABEL<0)CYCLE          

                  N=I+(J-1)*LX+(K-1)*LX*LY
                  !PRINT*,I
                  error_flag=.false.
                  
                  CALL MESH_B(N,U,L,T,INI_T,error_flag,e)
                  
                  IF(error_flag)THEN
                  print*,'--------'
                  print*,e
                  PRINT*,'An error occurred, retrying...'
                  pause
                  ENDIF

                  !!PRINT*,N
                  !!DO NUM=1,375
                  !!    IF(ABS(U(3,NUM))>1E-5)PRINT*,(NUM,U(3,NUM))
                  !!ENDDO
                  !!PRINT*,'----------'
                  CALL ASSIGN_VALUE(I,U(3,:),NX,L,L,L,XZ)

              ENDDO
              !$OMP END DO
              !$OMP END PARALLEL
              CALL UNION_Z_Y(X,XZ,NX,L,NZ,L,K)
              
              ENDDO
              


              DO YN=1,L-3
                  CALL OUTPUT_Y(NX,NZ,X(:,YN,:,:),(L-3)*(J-1)+YN,T)
              ENDDO

              DO YN=1,3
                  X1(:,YN,:,:)=X1(:,L-3+YN,:,:)  
                  X2(:,YN,:,:)=X2(:,L-3+YN,:,:) 
                  X(:,YN,:,:)=X(:,L-3+YN,:,:) 
              ENDDO
              DO YN=4,L
                  X(:,YN,:,:)=0
              ENDDO


          ENDDO

          DO YN=1,3
              CALL OUTPUT_Y(NX,NZ,X(:,YN,:,:),(L-3)*(LY+1-PY)+YN,T)
          ENDDO
         call system_clock(count12,count_rate,count_max)
          timespend=(count12-count11)/real(count_rate)
          write(*,*) '内部耗费的时间为（单位秒）: ',timespend          
         
          CALL GET_PML(T,NX,NY,NZ,L,PMLD,IS_PML,XSTEP,YSTEP,ZSTEP,TSTEP)
          call system_clock(count13,count_rate,count_max)

          CALL PML_TO_OUTPUT(T,NX,NY,NZ)

          timespend=(count13-count12)/real(count_rate)
          write(*,*) 'pml耗费的时间为（单位秒）: ',timespend 


      ENDDO
          
          
      call system_clock(count2,count_rate,count_max)
      timespend=(count2-count1)/real(count_rate)/real(NT-T0+1)
      write(*,*) '平均耗费的时间为（单位秒）: ',timespend 

      write(*,*)'finish'
      pause
      END PROGRAM main

      
      