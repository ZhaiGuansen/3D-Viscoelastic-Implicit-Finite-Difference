      module partial_main
      USE LAPACK95
      use matrix_construction
      USE DIA_SM
      use INI
      USE PML
      USE PDC
      USE cube_of_five_meshes
      USE make_cube
      use omp_lib
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)   
      contains
      
      
      SUBROUTINE main_a(NX,NY,NZ,L,INI_T,T0,NT,XSTEP,YSTEP,ZSTEP,TSTEP,
     + FP,AMPL,X0,PMLD,IS_INI)
C基础变量
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER,intent(in)::NX,NY,NZ,T0,NT,INI_T,x0(3),PMLD(3)
      INTEGER::LX,N2,LY,L
      INTEGER::I,J,K,T,NUM,N,IS_INI,P
      REAL,intent(in):: XSTEP,YSTEP,ZSTEP,TSTEP,FP,AMPL
      INTEGER,DIMENSION(:),allocatable::LABEL
      REAL,DIMENSION(:,:,:),allocatable::VP,VS,NP,NS,RHO!DIMENSION(NX,NY,NZ)
C时间计算相关变量
      integer(kind=8)  ::count1,count2,count_rate,count_max
      real timespend
C物理边界点集
      !INTEGER::PBN
      !INTEGER,DIMENSION(PBN,3)::PBX


      call omp_set_num_threads(18)
      
      LX=(NX-3)/(L-3)
      N2=LX*((NZ-3)/(L-3))
      LY=(NY-3)/(L-3)

      ALLOCATE(LABEL(N2*LY))

C测试赋值
      ALLOCATE(VP(NX,NY,NZ),VS(NX,NY,NZ),NP(NX,NY,NZ),
     + NS(NX,NY,NZ),RHO(NX,NY,NZ))

      
      VP=3500
      VS=2100
      NP=26
      NS=22
      !NP=0
      !NS=0
      RHO=2400
      !煤层
      DO I=1,NX
          DO J=1,NY
              DO K=29,38
                  VP(I,J,K)=2200
                  VS(I,J,K)=1200
                  NP(I,J,K)=23
                  NS(I,J,K)=18
                  RHO(I,J,K)=1400
              END DO
          END DO
      END DO
      
      !DO I=70,NX
      !    DO J=1,NY
      !        DO K=32,41 
      !            VP(I,J,K)=2200
      !            VS(I,J,K)=1200
      !            NP(I,J,K)=23
      !            NS(I,J,K)=18
      !            RHO(I,J,K)=1400
      !        END DO
      !    END DO
      !END DO


      !DO I=31,80
      !    DO J=16,21
      !        DO K=16,21   
      !            VP(I,J,K)=2200
      !            VS(I,J,K)=1200
      !            NP(I,J,K)=23
      !            NS(I,J,K)=18
      !            RHO(I,J,K)=1400
      !        END DO
      !    END DO
      !END DO
      
      
      
      
      
      
      
      
      LABEL=0
      CALL SYSTEM('folder.bat')
      call system("DEL/q " // '.\PML')!pml文件提前存在会影响判断,此处为调用Windows命令进行文件夹删除,可以手动删除
      call system("DEL/q " // '.\INI')
      !call system("cls")
      !pml边界信息，标签初始化
      CALL BASE_CUBE(NX,NY,NZ,L,VP,VS,RHO,XSTEP,YSTEP,ZSTEP,TSTEP,LABEL,
     + PMLD)
      PRINT*,'BASE_CUBE FINISHED'
      
 ! 对炮点周围值初始化,output1和2.dat文件作为初始化
      CALL MAKE_INI_CUBE(NX,NY,NZ,L,INI_T,X0(1),X0(2),X0(3),FP,TSTEP,
     + RHO(X0(1),X0(2),X0(3)),VP(X0(1),X0(2),X0(3)),
     + NP(X0(1),X0(2),X0(3)),AMPL,LABEL)
      PRINT*,'INI_CUBE FINISHED'

      IF(IS_INI==0)THEN 
CCC
     
      
      call system_clock(count1,count_rate)
      N=0
      P=0
      !$OMP PARALLEL PRIVATE(NUM,I,J,K) 
      !$omp do
      DO NUM=1,N2*LY
          CALL N_TO_X0(NUM,NX,NY,NZ,L,I,J,K)
     
          CALL MESH_A(NUM,L,XSTEP,YSTEP,ZSTEP,TSTEP,VP(I:I+L-1,J:J+L-1,
     +K:K+L-1),VS(I:I+L-1,J:J+L-1,K:K+L-1),NP(I:I+L-1,J:J+L-1,K:K+L-1),
     +NS(I:I+L-1,J:J+L-1,K:K+L-1),LABEL(NUM))    
          !$OMP CRITICAL
          N=N+1
          IF(N*100/(N2*LY)>P)THEN
              P=N*100/(N2*LY)
              WRITE(*,*)P,'%'
          ENDIF
          !$OMP END CRITICAL
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
      call system_clock(count2,count_rate,count_max)
      write(*,*) '计数速率count_rate    : ', count_rate
      write(*,*) '最大时间戳count_max   : ', count_max
      timespend=(count2-count1)/real(count_rate)
      write(*,*) '构建系数矩阵耗费的时间为（单位秒）: ',timespend 
cccc
      ENDIF



      DEALLOCATE(VP,VS,NP,NS,RHO)
      end
      
      END module partial_main