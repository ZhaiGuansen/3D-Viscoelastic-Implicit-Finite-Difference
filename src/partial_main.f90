      module partial_main

      use omp_lib
      use matrix_construction
      USE DIA_SM
      use INI
      USE PML
      USE PDC
      USE cube_of_five_meshes
      USE make_cube

      IMPLICIT DOUBLE PRECISION (A-H,P-Z)   
      contains
      
      
      SUBROUTINE main_a(NX,NY,NZ,L,INI_T,XSTEP,YSTEP,ZSTEP,TSTEP,FP,AMPL,X0,PMLD,IS_INI)
!!Basic variables
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER,intent(in)::NX,NY,NZ,INI_T,x0(3),PMLD(3)
      INTEGER::LX,N2,LY,L
      INTEGER::I,J,K,NUM,N,IS_INI,P
      REAL,intent(in):: XSTEP,YSTEP,ZSTEP,TSTEP,FP,AMPL
      INTEGER,DIMENSION(:),allocatable::LABEL,M_MAP,MN
      REAL,DIMENSION(:,:,:),allocatable::VP,VS,NP,NS,RHO!DIMENSION(NX,NY,NZ)
!Time related variables
      integer ::count1,count2,count_rate,count_max
      real timespend
      character(len=10) :: os



      call omp_set_num_threads(18)
      
      LX=(NX-3)/(L-3)
      N2=LX*((NZ-3)/(L-3))
      LY=(NY-3)/(L-3)

      ALLOCATE(LABEL(N2*LY),M_MAP(N2*LY))

! Test value
      ALLOCATE(VP(NX,NY,NZ),VS(NX,NY,NZ),NP(NX,NY,NZ),NS(NX,NY,NZ),RHO(NX,NY,NZ))

      
      VP=3500
      VS=2100
      NP=SQRT(3500**2/200/2/PI/FP)
      NS=SQRT(2100**2/100/2/PI/FP)
      !NP=0
      !NS=0
      RHO=2400
      
      !Coal seam
      !DO I=1,NX
      !    DO J=1,NY
      !        DO K=26,35
      !            VP(I,J,K)=2200
      !            VS(I,J,K)=1200
      !            NP(I,J,K)=SQRT(2200**2/40/2/PI/FP)
      !            NS(I,J,K)=SQRT(1200**2/20/2/PI/FP)
      !            RHO(I,J,K)=1400
      !        END DO
      !    END DO
      !END DO
      

      
      LABEL=0
      os = get_os()
      call cleanup_dir('./INI', os)
      call cleanup_dir('./PML', os)

      call create_dir('./PML', os)
      call create_dir('./INI', os)
      call create_dir('./BASE', os)
      call create_dir('./M', os)
      call create_dir('./OUTPUT', os)
      call create_dir('./Y', os)

      !PML boundary information, label initialization
      CALL BASE_CUBE(NX,NY,NZ,L,VP,VS,RHO,XSTEP,YSTEP,ZSTEP,TSTEP,LABEL,PMLD)

      CALL MAKE_MAP(LX,LY,(NZ-3)/(L-3),L,VP,LABEL,M_MAP)

      WRITE(*,*)'BASE_CUBE FINISHED'
      call flush()
      
 !Initialize values around shot point, output1.dat and output2.dat files serve as initialization
      CALL MAKE_INI_CUBE(NX,NY,NZ,L,INI_T,X0(1),X0(2),X0(3),FP,TSTEP,RHO(X0(1),X0(2),X0(3)),VP(X0(1),X0(2),X0(3)),NP(X0(1),X0(2),X0(3)),AMPL,LABEL)
      WRITE(*,*)'INI_CUBE FINISHED'
      call flush()

      IF(IS_INI/=2)THEN 
     
      
      call system_clock(count1,count_rate)
      
      
      IF(IS_INI==3)THEN!!!Rewrite specific M elements
      
      CALL CHECK_M(N2*LY,MN)
      N=N2*LY-SIZE(MN)
      P=0
      !$OMP PARALLEL PRIVATE(NUM,I,J,K) 
      !$omp do schedule(dynamic)
      DO NUM=1,SIZE(MN)
          CALL N_TO_X0(MN(NUM),NX,NY,NZ,L,I,J,K)
     
          CALL MESH_A(MN(NUM),L,XSTEP,YSTEP,ZSTEP,TSTEP,FP,VP(I:I+L-1,J:J+L-1,K:K+L-1),VS(I:I+L-1,J:J+L-1,K:K+L-1),&
              NP(I:I+L-1,J:J+L-1,K:K+L-1),NS(I:I+L-1,J:J+L-1,K:K+L-1),LABEL(MN(NUM)),M_MAP)  
          !$OMP CRITICAL
          N=N+1
          IF(N*100/(N2*LY)>P)THEN
              P=N*100/(N2*LY)
              WRITE(*,*)P,'%'
              call flush()
          ENDIF
          !$OMP END CRITICAL
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
      ELSE
      N=0
      P=0
      !$OMP PARALLEL PRIVATE(NUM,I,J,K) 
      !$omp do schedule(dynamic)
      DO NUM=1,N2*LY
          CALL N_TO_X0(NUM,NX,NY,NZ,L,I,J,K)
     
          CALL MESH_A(NUM,L,XSTEP,YSTEP,ZSTEP,TSTEP,FP,VP(I:I+L-1,J:J+L-1,K:K+L-1),VS(I:I+L-1,J:J+L-1,K:K+L-1),&
              NP(I:I+L-1,J:J+L-1,K:K+L-1),NS(I:I+L-1,J:J+L-1,K:K+L-1),LABEL(NUM),M_MAP)  
          !WRITE(*,*)NUM
          !$OMP CRITICAL
          N=N+1
          IF(N*100/(N2*LY)>P)THEN
              P=N*100/(N2*LY)
              WRITE(*,*)P,'%'
              call flush()
          ENDIF
          !$OMP END CRITICAL
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
      ENDIF
      
      call system_clock(count2,count_rate,count_max)
      WRITE(*,*) 'count_rate    : ', count_rate
      WRITE(*,*) 'count_max   : ', count_max
      timespend=(count2-count1)/real(count_rate)
      WRITE(*,*) 'Time consumed for constructing coefficient matrix (in seconds): ',timespend 
      call flush()
      ENDIF



      DEALLOCATE(VP,VS,NP,NS,RHO)
      end

      END module partial_main