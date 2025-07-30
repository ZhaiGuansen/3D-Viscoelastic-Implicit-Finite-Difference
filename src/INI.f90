      module INI
      USE DIA_SM
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)   
      contains
      
      INTEGER FUNCTION POSITION(X0,X1,X2,X3,NX,NY,NZ)
          INTEGER :: X0, X1, X2, X3, NX, NY, NZ
          POSITION = X1 + (X2-1)*NX + (X3-1)*NX*NY + (X0-1)*NX*NY*NZ
      END FUNCTION POSITION


      SUBROUTINE initial_value_B(NX,NY,NZ,N,X,U,V,W,OUTPUT_B)
      ! Initial value input
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      REAL,PARAMETER::PAI=3.1415926
      INTEGER::NX,NY,NZ,N,I
      INTEGER,DIMENSION(N,3)::X
      REAL(KIND=RK),DIMENSION(N)::U,V,W
      REAL(KIND=RK) OUTPUT_B(3*NX*NY*NZ)
      
      DO I=1,N
          IF(X(I,1)<1.OR.X(I,1)>NX)CYCLE
          IF(X(I,2)<1.OR.X(I,2)>NY)CYCLE
          IF(X(I,3)<1.OR.X(I,3)>NZ)CYCLE
          OUTPUT_B(POSITION(1,X(I,1),X(I,2),X(I,3),NX,NY,NZ)) = OUTPUT_B(POSITION(1,X(I,1),X(I,2),X(I,3),NX,NY,NZ)) + U(I)
          OUTPUT_B(POSITION(2,X(I,1),X(I,2),X(I,3),NX,NY,NZ)) = OUTPUT_B(POSITION(2,X(I,1),X(I,2),X(I,3),NX,NY,NZ)) + V(I)
          OUTPUT_B(POSITION(3,X(I,1),X(I,2),X(I,3),NX,NY,NZ)) = OUTPUT_B(POSITION(3,X(I,1),X(I,2),X(I,3),NX,NY,NZ)) + W(I)
      END DO

      END
      
      
      SUBROUTINE CHECK_M(N,OUTPUT)
      ! M-file validation
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      INTEGER::N,I,FSIZE,ESIZE,IDX
      LOGICAL::IS
      INTEGER,DIMENSION(:),allocatable::MN,OUTPUT
      CHARACTER(LEN=30)::F
      CHARACTER(LEN=10)::C
      
      ALLOCATE(MN(N))
      ESIZE=0
      IDX=0
      DO I=1,N
          WRITE(C,'(I6)')I
          F='./M/M'//TRIM(ADJUSTL(C))//'.DAT'
      
          INQUIRE(FILE=F, EXIST=IS, SIZE=FSIZE)
          IF(.NOT.IS)THEN
              IDX=IDX+1
              MN(IDX)=I
          ELSE
              IF(ESIZE==0)ESIZE=FSIZE
              IF(ESIZE>FSIZE)THEN
                  IDX=IDX+1
                  MN(IDX)=I
              ENDIF
          ENDIF
      ENDDO
      OUTPUT=MN(1:IDX)
      DEALLOCATE(MN)
      END
      
      
      SUBROUTINE MAKE_MAP(LX,LY,LZ,L,VP,LABEL,M_MAP)
      INTEGER::LX,LY,LZ,L,I,J,K,N,NUM,X0,Y0,Z0,IDX
      REAL,DIMENSION(:,:,:),allocatable::VP
      INTEGER,DIMENSION(:),allocatable::LABEL,M_MAP,M0
      
      ALLOCATE(M0(LX*LY*LZ))
      M_NAP=0
      IDX=0
      DO K=1,LZ
          DO J=1,LY
              DO I=1,LX
                  NUM=I+(J-1)*LX+(K-1)*LX*LY
                  IF(LABEL(NUM)>3)THEN
                      M_MAP(NUM)=NUM
                  ELSEIF(IDX==0)THEN
                      IDX=IDX+1
                      M0(IDX)=NUM
                      M_MAP(NUM)=NUM
                  ELSE
                      DO N=1,IDX
                          Z0=(M0(N)-1)/(LX*LY)
                          X0=M0(N)-Z0*LX*LY
                          Z0=Z0+1
                          Y0=(X0-1)/LX
                          X0=X0-Y0*LX
                          Y0=Y0+1
                          
                          
      IF(COMPARE_V(L,VP((L-3)*(I-1)+1:(L-3)*(I-1)+L,(L-3)*(J-1)+1:(L-3)*(J-1)+L,(L-3)*(K-1)+1:(L-3)*(K-1)+L), &
    VP((L-3)*(X0-1)+1:(L-3)*(X0-1)+L,(L-3)*(Y0-1)+1:(L-3)*(Y0-1)+L,(L-3)*(Z0-1)+1:(L-3)*(Z0-1)+L)).eqv..TRUE.)THEN
                              M_MAP(NUM)=M0(N)
                              EXIT
                          ENDIF
                      ENDDO
                      IF(M_MAP(NUM)==0)THEN
                          IDX=IDX+1
                          M0(IDX)=NUM
                          M_MAP(NUM)=NUM
                      ENDIF
                  ENDIF
              ENDDO
          ENDDO
      ENDDO
      DEALLOCATE(M0)
      
      END
      
      
      LOGICAL FUNCTION COMPARE_V(L,V1,V2) RESULT(is_equal)
      INTEGER, INTENT(IN) :: L
      REAL,DIMENSION(L,L,L)::V1, V2
      INTEGER :: i, j, k

      is_equal = .TRUE.
      DO k = 1, L
        DO j = 1, L
            DO i = 1, L
                IF (ABS(V1(i,j,k) - V2(i,j,k)) > 1E-6) THEN
                    is_equal = .FALSE.
                    RETURN
                END IF
            END DO
        END DO
      END DO
      END FUNCTION COMPARE_V



function get_os() result(os)
      ! This function returns the operating system type as a string.
      ! It can be used to determine the appropriate commands for file operations.
    implicit none
    character(len=10) :: os
    character(len=32) :: env
    integer :: length

    call get_environment_variable("OS", env, length)
    if (index(env, "Windows") > 0) then
        os = "Windows"
    else
        os = "Linux/Mac"
    end if
end function get_os


subroutine cleanup_dir(dir_path,os)
!Delete folder and its contents
    implicit none
    character(len=*), intent(in) :: dir_path
    character(len=100) :: command
    character(len=10) :: os

    if (trim(os) == "Windows") then
        command = 'rmdir /S /Q "' // trim(dir_path) // '"'
    else
        command = 'rm -rf "' // trim(dir_path) // '"'
    end if

    call system(command)
end subroutine cleanup_dir

subroutine create_dir(dir_path,os)
!Create a directory if it does not exist
    implicit none
    character(len=*), intent(in) :: dir_path
    character(len=100) :: command
    character(len=10) :: os

    if (trim(os) == "Windows") then
        command = 'mkdir "' // trim(dir_path) // '" 2> nul || ver > nul'
    else
        command = 'mkdir -p "' // trim(dir_path) // '"'
    end if

    call system(command)
end subroutine create_dir
      
      
      end module INI