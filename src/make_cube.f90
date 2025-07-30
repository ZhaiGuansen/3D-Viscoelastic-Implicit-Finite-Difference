      module make_cube

      use matrix_construction
      USE DIA_SM
      use INI
      USE PML
      USE PDC
      USE cube_of_five_meshes
      
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)   
      contains
      
      SUBROUTINE MAKE_INI_CUBE(NX,NY,NZ,L,NT,X,Y,Z,FP,TSTEP,RHO,VP,NP,AMPL,LABEL)
      ! Build initialization metadata
      integer, parameter :: rk = kind ( 1.0D+00 )
      REAL,PARAMETER::PAI=3.1415926
      INTEGER::T,NX,NY,NZ,NT,X,Y,Z,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
      INTEGER::N,DX,DY,DZ,LX,LY,LZ,X2,Y2,Z2,L
      INTEGER::LABEL(((NX-3)/2)*((NY-3)/2)*((NZ-3)/2))
      REAL::TSTEP,FP,AMPL,RHO,VP,NP,VP1,NP1
      REAL(KIND=RK),DIMENSION(:,:,:,:),allocatable::OUTPUT
      REAL(KIND=RK)::U(NT),V(NT),W(NT),SR(NT)
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::CTEMP,CT
      
      U=0
      V=0
      W=0
      
      VP1=SQRT(VP**2-NP**4)
      NP1=NP**2/2
      call ricker(NT,FP,TSTEP,AMPL,SR)
      DO T=1,NT/2
          !WRITE(*,*)SR(NT-NT/2+T)
          U(T)=SR(NT-NT/2+T)
          V(T)=SR(NT-NT/2+T)
          W(T)=SR(NT-NT/2+T)
      END DO
      DO T=NT/2+1,NT
          !WRITE(*,*)SR(T-NT/2)
          U(T)=SR(T-NT/2)
          V(T)=SR(T-NT/2)
          W(T)=SR(T-NT/2)
      END DO  
            
      ALLOCATE(OUTPUT(NX,NY,NZ,3))
      OUTPUT=0
      OUTPUT(X,Y,Z,1)=U(1)
      OUTPUT(X,Y,Z,2)=V(1)
      OUTPUT(X,Y,Z,3)=W(1)
      CALL OUTPUT_X_T(1,NX,NY,NZ,OUTPUT)
      OUTPUT(X,Y,Z,1)=(U(2)+U(1)*EXP(-NP1*1E-4)*SIN(VP1*1E-4))
      OUTPUT(X,Y,Z,2)=(V(2)+U(1)*EXP(-NP1*1E-4)*SIN(VP1*1E-4))
      OUTPUT(X,Y,Z,3)=(W(2)+U(1)*EXP(-NP1*1E-4)*SIN(VP1*1E-4))
      CALL OUTPUT_X_T(2,NX,NY,NZ,OUTPUT)
      DEALLOCATE(OUTPUT)
      
      
      XMIN=(X+L-7)/(L-3)
      XMAX=(X+L-4)/(L-3)
      YMIN=(Y+L-7)/(L-3)
      YMAX=(Y+L-4)/(L-3)
      ZMIN=(Z+L-7)/(L-3)
      ZMAX=(Z+L-4)/(L-3)
      LX=(NX-3)/(L-3)
      LY=(NY-3)/(L-3)
      LZ=(NZ-3)/(L-3)
      
      
      
      
      DO X2=XMIN,XMAX
          IF(X2>LX)EXIT
          DO Y2=YMIN,YMAX
              IF(Y2>LY)EXIT
              DO Z2=ZMIN,ZMAX
                  IF(Z2>LZ)EXIT
      N=X2+(Y2-1)*LX+(Z2-1)*LX*LY
      IF(MOD(LABEL(N),2).EQ.0)LABEL(N)=LABEL(N)+1
      DX=X-(X2-1)*(L-3)
      DY=Y-(Y2-1)*(L-3)
      DZ=Z-(Z2-1)*(L-3)
      
      WRITE(CTEMP,'(I6)')N
      
      DO T=3,NT
          WRITE(CT,'(I6)')T
      FILENAME='./INI/INI'//TRIM(ADJUSTL(CTEMP))//'_' //TRIM(ADJUSTL(CT))//'.DAT'
      OPEN(N+100,FILE=FILENAME)
      WRITE(N+100,*)8
      WRITE(N+100,*)TSTEP,RHO
      WRITE(N+100,*)DX,DY,DZ
      WRITE(N+100,*)DX-1,DY,DZ
      WRITE(N+100,*)DX,DY-1,DZ
      WRITE(N+100,*)DX-1,DY-1,DZ
      WRITE(N+100,*)DX,DY,DZ-1
      WRITE(N+100,*)DX-1,DY,DZ-1
      WRITE(N+100,*)DX,DY-1,DZ-1
      WRITE(N+100,*)DX-1,DY-1,DZ-1
      WRITE(N+100,*)U(T),-U(T),U(T),-U(T),U(T),-U(T),U(T),-U(T)
      WRITE(N+100,*)V(T),V(T),-V(T),-V(T),V(T),V(T),-V(T),-V(T)
      WRITE(N+100,*)W(T),W(T),W(T),W(T),-W(T),-W(T),-W(T),-W(T)
      CLOSE(N+100)
      ENDDO
              ENDDO
          ENDDO
      ENDDO


      END
      
      SUBROUTINE RETRIEVE_Y(NX,NY,NZ,T)
      ! Reconstruct y-file using output file data
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::T,NX,NY,NZ,I,J,K,NUM
      REAL(KIND=RK)::VALUE
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::C1,C2
      
      WRITE(C1,'(I6)')T
      
      FILENAME='./OUTPUT/OUTPUT'//TRIM(ADJUSTL(C1))//'.DAT'
      OPEN(NY+101,FILE=FILENAME)
      READ(NY+101,*)I,J,K
      
      DO J=1,NY
      WRITE(C2,'(I6)')J
      FILENAME='./Y/Y'//TRIM(ADJUSTL(C1))//'_' //TRIM(ADJUSTL(C2))//'.DAT'
      OPEN(J+100,FILE=FILENAME)
      ENDDO
      
      DO NUM=1,3
          DO K=1,NZ
              DO J=1,NY
                  DO I=1,NX
                      READ(NY+101,*)VALUE
                      WRITE(J+100,*)VALUE
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      
      DO J=1,NY+1
      CLOSE(J+100)
      ENDDO

      END
      
      
      
      
      end module make_cube