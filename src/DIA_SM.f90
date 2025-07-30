      module DIA_SM
      
      !use INI
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      REAL::UMIN=1E-10
      INTEGER::NUM0=0
     
      contains


! Sparse matrix addition (A is overwritten by result)
      subroutine DIA_ADDSM(AN,AM1,AM2,AS1,AS2,AS3,A,BN,BM1, BM2,BS1,BS2,BS3,B)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )

      real(kind=rk),allocatable,DIMENSION(:)::A,B,X,Y
      integer,allocatable,DIMENSION(:)::AM1,AM2,BM1,BM2,M1,M2,AS1,AS2,AS3,BS1,BS2,BS3,S1,S2,S3,NS,NL
      INTEGER::AN(2),BN(2),I,J,K,N,II,JJ,MAX_SIZE,idx(3)!real size of m,s,x
      
      IF(AN(1)/=BN(1))THEN
          WRITE(*,*)"ADD ERROR"
          call flush()
      ENDIF
      
      ALLOCATE(Y(AN(1)))
      MAX_SIZE=MIN(AN(2) + BN(2),2*AN(1)-1)
      ALLOCATE(M1(MAX_SIZE),M2(MAX_SIZE))
      MAX_SIZE=MIN(MAX_SIZE*(MAXVAL(AM2)+MAXVAL(BM2)),SIZE(AS1)+SIZE(AS2))
      ALLOCATE(S1(MAX_SIZE),S3(MAX_SIZE))
      MAX_SIZE=SIZE(A)+SIZE(B)
      ALLOCATE(X(MAX_SIZE))
      
      IDX=0
      I=1
      J=1
      K=1
      NI=1
      NJ=1
      DO WHILE(I<=AN(2))
          IF(J>BN(2))THEN
              IDX(1)=IDX(1)+1
              M1(IDX(1))=AM1(I)
              M2(IDX(1))=AM2(I)
             DO II=1,AM2(I)
                 IDX(2)=IDX(2)+1
                 S1(IDX(2))=AS1(NI)
                 S3(IDX(2))=AS3(NI)
                 DO JJ=1,AS3(NI)
                     IDX(3)=IDX(3)+1
                     X(IDX(3))=A(AS2(NI)+JJ-1)
                 ENDDO
                 NI=NI+1
             ENDDO
             I=I+1
          ELSE
          
         IF(AM1(I)<BM1(J))THEN
             IDX(1)=IDX(1)+1
              M1(IDX(1))=AM1(I)
              M2(IDX(1))=AM2(I)
             DO II=1,AM2(I)
                 IDX(2)=IDX(2)+1
                 S1(IDX(2))=AS1(NI)
                 S3(IDX(2))=AS3(NI)
                 DO JJ=1,AS3(NI)
                     IDX(3)=IDX(3)+1
                     X(IDX(3))=A(AS2(NI)+JJ-1)
                 ENDDO
                 NI=NI+1
             ENDDO
             I=I+1
         ELSEIF(AM1(I)>BM1(J))THEN
             IDX(1)=IDX(1)+1
             M1(IDX(1))=BM1(J)
              M2(IDX(1))=BM2(J)
             DO II=1,BM2(J)
                 IDX(2)=IDX(2)+1
                 S1(IDX(2))=BS1(NJ)
                 S3(IDX(2))=BS3(NJ)
                 DO JJ=1,BS3(NJ)
                     IDX(3)=IDX(3)+1
                     X(IDX(3))=B(BS2(NJ)+JJ-1)
                 ENDDO
                 NJ=NJ+1
             ENDDO
             J=J+1
         ELSE
              IDX(1)=IDX(1)+1
              M1(IDX(1))=AM1(I) 
             Y=0
             DO II=1,AM2(I)
                 DO JJ=1,AS3(NI)
                     Y(AS1(NI)+JJ-1)=Y(AS1(NI)+JJ-1)+A(AS2(NI)+JJ-1)
                 ENDDO
                 NI=NI+1
             ENDDO
             DO II=1,BM2(J)
                 DO JJ=1,BS3(NJ)
                     Y(BS1(NJ)+JJ-1)=Y(BS1(NJ)+JJ-1)+B(BS2(NJ)+JJ-1)
                 ENDDO
                 NJ=NJ+1
             ENDDO
             
             
             CALL DLINE(Y,AN(1)-ABS(AM1(I)-AN(1)),N,NS,NL)
             M2(IDX(1))=N

            S1(IDX(2)+1:IDX(2)+SIZE(NS))=NS
            S3(IDX(2)+1:IDX(2)+SIZE(NS))=NL
            IDX(2)=IDX(2)+SIZE(NS)
             
             DO II=1,N
                 DO JJ=1,NL(II)
                     IDX(3)=IDX(3)+1
                     X(IDX(3))=Y(NS(II)+JJ-1)
                 ENDDO
             ENDDO  
             I=I+1
             J=J+1
             IF(SIZE(NS)>0)THEN
                  DEALLOCATE(NS,NL)
                  ALLOCATE(NS(0),NL(0))
             ENDIF
         ENDIF
      ENDIF
         K=K+1
      ENDDO
      DO WHILE(J<=BN(2))
          IDX(1)=IDX(1)+1
          M1(IDX(1))=BM1(J)
          M2(IDX(1))=BM2(J)
          DO II=1,BM2(J)
              IDX(2)=IDX(2)+1
              S1(IDX(2))=BS1(NJ)
              S3(IDX(2))=BS3(NJ)
              DO JJ=1,BS3(NJ)
                  IDX(3)=IDX(3)+1
                  X(IDX(3))=B(BS2(NJ)+JJ-1)
              ENDDO
              NJ=NJ+1
          ENDDO
          J=J+1
          K=K+1
      ENDDO
      ALLOCATE(S2(IDX(2)))
      S2(1)=1
      DO I=2,IDX(2)
          S2(I)=S2(I-1)+S3(I-1)
      ENDDO
      AN(2)=K-1
      AM1=M1(1:IDX(1))
      AM2=M2(1:IDX(1))
      AS1=S1(1:IDX(2))
      AS2=S2
      AS3=S3(1:IDX(2))
      A=X(1:IDX(3))
      DEALLOCATE(M1,M2,S1,S2,S3,X)
      END
      
      ! Check for zero-value entries
      subroutine DIA_SORT(AN,AM1,AM2,AS1,AS2,AS3,A)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )

      real(kind=rk),allocatable,DIMENSION(:)::A,X,Y
      integer,allocatable,DIMENSION(:)::AM1,AM2,M1,M2,AS1,AS2,AS3,S1,S2,S3,NS,NL
      INTEGER::AN(2),I,J,K,N,NN,IDX(3)!real size of m,s,x

      ALLOCATE(X(SIZE(A)))
      ALLOCATE(M1(SIZE(AM1)),M2(SIZE(AM1)))
      ALLOCATE(S1(SIZE(A)),S3(SIZE(A)))
      IDX=0
      ALLOCATE(Y(AN(1)))
      NN=1
      DO I=1,AN(2)
          Y=0


          DO J=1,AM2(I)
              DO K=1,AS3(NN)
                  Y(AS1(NN)+K-1)=Y(AS1(NN)+K-1)+A(AS2(NN)+K-1)
              ENDDO
              NN=NN+1
          ENDDO
             
          CALL DLINE(Y,AN(1)-ABS(AM1(I)-AN(1)),N,NS,NL)

          IF(N==0)THEN
              CYCLE
          ELSE
          IDX(1)=IDX(1)+1
              M1(IDX(1))=AM1(I)
          ENDIF
          M2(IDX(1))=N
          S1(IDX(2)+1:IDX(2)+SIZE(NS))=NS
          S3(IDX(2)+1:IDX(2)+SIZE(NL))=NL
          IDX(2)=IDX(2)+SIZE(NS)

          DO J=1,N
              DO K=1,NL(J)
                  IDX(3)=IDX(3)+1
                  X(IDX(3))=Y(NS(J)+K-1)
              ENDDO
          ENDDO  

          DEALLOCATE(NS,NL)
          ALLOCATE(NS(0),NL(0))
      ENDDO

      ALLOCATE(S2(IDX(2)))
      S2(1)=1
      DO I=2,IDX(2)
          S2(I)=S2(I-1)+S3(I-1)
      ENDDO
      AN(2)=IDX(1)
      AM1=M1(1:IDX(1))
      AM2=M2(1:IDX(1))
      AS1=S1(1:IDX(2))
      AS2=S2
      AS3=S3(1:IDX(2))
      A=X(1:IDX(3))
      DEALLOCATE(M1,M2,S1,S2,S3,X)


      END
      
      SUBROUTINE DLINE(Y,LENGTH,N,M,M1)
      ! Retrieve positions and values of non-zero elements in array
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::N,I,J,LENGTH,NUM,IDX
      real(kind=rk),allocatable::Y(:)
      integer,allocatable::M(:),M1(:),NS(:),NL(:)

IF(LENGTH<=0)THEN
        N=0
        RETURN
ENDIF
        ALLOCATE(NS(LENGTH),NL(LENGTH))
        IDX=0
      I=1
      DO WHILE(I<=LENGTH)
          IF(ABS(Y(I))<UMIN)THEN
              I=I+1
          ELSE
                IDX=IDX+1
              NS(IDX)=I
              NUM=0
              J=1
              DO WHILE(NUM<=NUM0)
                  I=I+1
                  IF(I>LENGTH)EXIT
                  J=J+1
                  IF(ABS(Y(I))<UMIN)THEN
                      NUM=NUM+1
                  ELSE
                      NUM=0
                  ENDIF
              ENDDO
              NL(IDX)=J-NUM
          ENDIF
      ENDDO
        IF(IDX==0)THEN
            N=0
            DEALLOCATE(NS,NL)
            RETURN
        ENDIF

      N=IDX
      M=NS(1:IDX)
      M1=NL(1:IDX)
      DEALLOCATE(NS,NL)
      END
! Create identity matrix
      SUBROUTINE DIA_ISM(N,AN,AM1,AM2,AS1,AS2,AS3,A)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::N,AN(2)
      real(kind=rk),allocatable::A(:)
      integer,allocatable,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3
      AN(1)=N
      AN(2)=1
      AM1=[N]
      AM2=[1]
      AS1=[1]
      AS2=[1]
      AS3=[N]
      ALLOCATE(A(N))
      A=1
      END
! Sparse matrix-vector multiplication
      SUBROUTINE DOTSM(AN,AM1,AM2,AS1,AS2,AS3,A,B,OUTPUT,LENGTH)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::AN(2),LENGTH,I,J,K,NS,AI,AJ
      REAL(KIND=RK),DIMENSION(LENGTH)::B,OUTPUT
      real(kind=rk),allocatable::A(:)
      integer,allocatable,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3

      IF(AN(1)/=LENGTH)THEN
          WRITE(*,*)'DOTSM ERROR'
          call flush()
      ENDIF
      J=1
      OUTPUT=0
      NS=1
      DO I=1,AN(2)
          DO J=1,AM2(I)
              DO K=1,AS3(NS)
                  IF(AM1(I)<=AN(1))THEN
                      AI=AN(1)+AS1(NS)+K-1-AM1(I)
                      AJ=AS1(NS)+K-1
                  ELSE
                      AI=AS1(NS)+K-1
                      AJ=AM1(I)+AS1(NS)+K-1-AN(1)
                  ENDIF
                  OUTPUT(AI)=OUTPUT(AI)+A(AS2(NS)+K-1)*B(AJ) 
              ENDDO
              NS=NS+1
          ENDDO
      ENDDO
      END  
      
      
      SUBROUTINE GMRES(AN,AM1,AM2,AS1,AS2,AS3,A,B,OUTPUT,LENGTH,LABEL)
      ! Solve large linear system using GMRES algorithm
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::M,AN(2),LENGTH,I,J
      INTEGER,ALLOCATABLE,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3
      REAL(KIND=RK),ALLOCATABLE::A(:)
      REAL(KIND=RK),ALLOCATABLE,DIMENSION(:,:)::Q,H,BM,YM
      REAL(KIND=RK)::NM,LABEL
      REAL(KIND=RK),DIMENSION(LENGTH)::B,OUTPUT,X,R
      
      
      
      M=6
      ALLOCATE(Q(M,LENGTH),H(M+1,M),YM(M,1),BM(M+1,1))
      Q=0
      H=0
      BM=0

      CALL DOTSM(AN,AM1,AM2,AS1,AS2,AS3,A,OUTPUT,X,LENGTH)
      R=B-X

      
      CALL NRM2(R,NM,LENGTH)
      IF(NM<UMIN)THEN
          LABEL=-1
          RETURN
      ENDIF
      LABEL=NM
      Q(1,:)=R/NM
      BM(1,1)=NM
      
      DO J=1,M
          CALL DOTSM(AN,AM1,AM2,AS1,AS2,AS3,A,Q(J,:),X,LENGTH)
          DO I=1,J
              H(I,J)=DOT_PRODUCT(X,Q(I,:))
              X=X-H(I,J)*Q(I,:)
          ENDDO
          CALL NRM2(X,NM,LENGTH)
          H(J+1,J)=NM
          IF(ABS(NM)<UMIN)THEN
              !WRITE(*,*)'M IS TOO LARGE'
              M=J
              Q=Q(1:M,:)
              H=H(1:M+1,1:M)
              YM=YM(1:M,:)
              BM=BM(1:M+1,:)
              EXIT
          ENDIF
          IF(J<M)Q(J+1,:)=X/H(J+1,J)
      ENDDO

      CALL LEAST_SQUARE(M,H,BM,YM) 
      DO I=1,M
          OUTPUT=OUTPUT+YM(I,1)*Q(I,:)
      ENDDO
          
      END
      
      
      SUBROUTINE LEAST_SQUARE(M,A,B,X)
      ! Least squares method, matrix size: (m+1, m)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::M,I,J,IPIV(M)
      REAL(KIND=RK),allocatable,DIMENSION(:,:)::A,AT,IA,C,B,X,BB
      
      ALLOCATE(AT(M,M+1),IA(M,M),C(M,M),BB(M,1))
      
      DO I=1,M+1
          DO J=1,M
              AT(J,I)=A(I,J)
          ENDDO
      ENDDO
      CALL MDOT(M,M+1,M,AT,A,IA)
      call dgetrf(M, M, IA, M, ipiv, info)
      if (info /= 0) then
        print *, 'Matrix is singular'
        STOP
      end if
      call dgetri(M, IA, M, ipiv, C,M*M, info)
      if (info /= 0) then
        print *, 'Inverse not computed'
        STOP
      end if

      
      CALL MDOT(M,M+1,1,AT,B,BB)
      CALL MDOT(M,M,1,IA,BB,X)
      
      
      END
      
      
      SUBROUTINE MDOT(M1,M2,M3,A,B,OUTPUT)
      ! Conventional matrix multiplication
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::M1,M2,M3,I,J,K
      REAL(KIND=RK),allocatable,DIMENSION(:,:)::A,B,OUTPUT
      
      IF(SIZE(OUTPUT)==0)THEN
          ALLOCATE(OUTPUT(M1,M3))
      ELSEIF(SIZE(OUTPUT,1)/=M1.OR.SIZE(OUTPUT,2)/=M3)THEN
          WRITE(*,*)'SIZE WRONG'
          call flush()
          STOP
      ENDIF
      
          
      OUTPUT=0
      DO I=1,M1
          DO J=1,M3
              DO K=1,M2
                  OUTPUT(I,J)=OUTPUT(I,J)+A(I,K)*B(K,J)
              ENDDO
          ENDDO
      ENDDO
      
      END
      ! Calculate 2-norm of matrix A
      SUBROUTINE NRM2(A,X,LENGTH)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::I,LENGTH
      REAL(KIND=RK)::X
      REAL(KIND=RK),DIMENSION(LENGTH)::A
      
      X=0
      DO I=1,LENGTH
          X=X+A(I)**2
      ENDDO
      X=SQRT(X)
      END
      
      
      SUBROUTINE DOT_ROW(AN,AM1,AM2,AS1,AS2,AS3,A,B,LENGTH)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::AN(2),LENGTH,I,J,AI
      REAL(KIND=RK),DIMENSION(LENGTH)::B
      real(kind=rk),allocatable::A(:)
      integer,allocatable,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3
      
      IF(AN(1)/=LENGTH)THEN
          WRITE(*,*)'DOTROW ERROR'
          call flush()
      ENDIF
      J=1
      NS=1
      DO I=1,AN(2)
          DO J=1,AM2(I)
              DO K=1,AS3(NS)
                  IF(AM1(I)<=AN(1))THEN
                      AI=AN(1)+AS1(NS)+K-1-AM1(I)
                  ELSE
                      AI=AS1(NS)+K-1
                  ENDIF
                  A(AS2(NS)+K-1)=A(AS2(NS)+K-1)*B(AI) 
              ENDDO
              NS=NS+1
          ENDDO
      ENDDO
      END  


      

     !! ! Symmetric matrix operation
     !! !!A4
     !! !!A2 A5
     !! !!A1 A3 A6
      ! Assemble 3*3 block matrix
      SUBROUTINE SPLICE33(N,M1,M2,S1,S2,S3,A,A1N,A1M1,A1M2,A1S1,A1S3,A1,&
          A2N,A2M1,A2M2,A2S1,A2S3,A2,A3N,A3M1,A3M2,A3S1,A3S3,A3,&
          A4N,A4M1,A4M2,A4S1,A4S3,A4,A5N,A5M1,A5M2,A5S1,A5S3,A5,A6N,A6M1,A6M2,A6S1,A6S3,A6)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::I,J,MK(6),SK(6),AK(6),MM
      !REAL(KIND=RK)::C0=0
      INTEGER,DIMENSION(2)::A1N,A2N,A3N,A4N,A5N,A6N,N
      INTEGER,ALLOCATABLE,DIMENSION(:)::M1,M2,S1,S2,S3,A1M1,A1M2,A1S1,A1S3,A2M1,A2M2,A2S1,A2S3,A3M1,A3M2,A3S1,A3S3,&
          A4M1,A4M2,A4S1,A4S3,A5M1,A5M2,A5S1,A5S3,A6M1,A6M2,A6S1,A6S3
      REAL(KIND=RK),ALLOCATABLE,DIMENSION(:)::A,A1,A2,A3,A4,A5,A6
      
      N(1)=3*A1N(1)
      DO I=1,A1N(2)
          M1=[M1,A1M1(I)]
          M1=[M1,A1M1(I)+4*A1N(1)]
      ENDDO
      DO I=1,A2N(2)
          M1=[M1,A2M1(I)+A1N(1)]
          M1=[M1,A2M1(I)+3*A1N(1)]
      ENDDO
      DO I=1,A3N(2)
          M1=[M1,A3M1(I)+A1N(1)]
          M1=[M1,A3M1(I)+3*A1N(1)]
      ENDDO
      DO I=1,A4N(2)
          M1=[M1,A4M1(I)+2*A1N(1)]
      ENDDO
      DO I=1,A5N(2)
          M1=[M1,A5M1(I)+2*A1N(1)]
      ENDDO
      DO I=1,A6N(2)
          M1=[M1,A6M1(I)+2*A1N(1)]
      ENDDO
      CALL sort_LINE(M1)
      N(2)=SIZE(M1)
      MK=1
      SK=1
      AK=1
      DO I=1,N(2)
          IF(M1(I)==3*A1N(1))THEN
              MK(1)=1
              MK(2)=1
              MK(3)=1
              SK(1)=1
              SK(2)=1
              SK(3)=1
              AK(1)=1
              AK(2)=1
              AK(3)=1
          ENDIF
          
          IF(M1(I)<=A1N(1))THEN
              IF(MK(1)<=A1N(2))THEN
              M2=[M2,A1M2(MK(1))]
              DO J=1,A1M2(MK(1))
                  S1=[S1,A1S1(SK(1))]
                  A=[A,A1(AK(1):AK(1)+A1S3(SK(1))-1)]
                  AK(1)=AK(1)+A1S3(SK(1))
                  S3=[S3,A1S3(SK(1))]
                  SK(1)=SK(1)+1
              ENDDO  
              MK(1)=MK(1)+1
            ENDIF
          
          ELSEIF(M1(I)<=2*A1N(1))THEN
              MM=0
              IF(MK(2)<=A2N(2))THEN
              IF(M1(I)==A1N(1)+A2M1(MK(2)))THEN
                  MM=MM+A2M2(MK(2))
                  DO J=1,A2M2(MK(2))
                      S1=[S1,A2S1(SK(2))]
                      A=[A,A2(AK(2):AK(2)+A2S3(SK(2))-1)]
                      AK(2)=AK(2)+A2S3(SK(2))
                      S3=[S3,A2S3(SK(2))]
                      SK(2)=SK(2)+1
                  ENDDO  
                  MK(2)=MK(2)+1
              ENDIF
              ENDIF
              
              IF(MK(1)<=A1N(2))THEN
              IF(M1(I)==A1M1(MK(1)))THEN
                  MM=MM+A1M2(MK(1))
                  DO J=1,A1M2(MK(1))
                      S1=[S1,A1S1(SK(1))+M1(I)-A1N(1)]
                      A=[A,A1(AK(1):AK(1)+A1S3(SK(1))-1)]
                      AK(1)=AK(1)+A1S3(SK(1))
                      S3=[S3,A1S3(SK(1))]
                      SK(1)=SK(1)+1
                  ENDDO  
                  MK(1)=MK(1)+1
              ENDIF
              ENDIF
              
              IF(MK(3)<=A3N(2))THEN
              IF(M1(I)==A1N(1)+A3M1(MK(3)))THEN
                  MM=MM+A3M2(MK(3))
                  DO J=1,A3M2(MK(3))
                      S1=[S1,A3S1(SK(3))+A1N(1)]
                      A=[A,A3(AK(3):AK(3)+A3S3(SK(3))-1)]
                      AK(3)=AK(3)+A3S3(SK(3))
                      S3=[S3,A3S3(SK(3))]
                      SK(3)=SK(3)+1
                  ENDDO  
                  MK(3)=MK(3)+1
              ENDIF
              ENDIF
              M2=[M2,MM]


          ELSEIF(M1(I)<=3*A1N(1))THEN
              MM=0
              IF(MK(4)<=A4N(2))THEN
              IF(M1(I)==2*A1N(1)+A4M1(MK(4)))THEN
                  MM=MM+A4M2(MK(4))
                  DO J=1,A4M2(MK(4))
                      S1=[S1,A4S1(SK(4))]
                      A=[A,A4(AK(4):AK(4)+A4S3(SK(4))-1)]
                      AK(4)=AK(4)+A4S3(SK(4))
                      S3=[S3,A4S3(SK(4))]
                      SK(4)=SK(4)+1
                  ENDDO  
                  MK(4)=MK(4)+1
              ENDIF
              ENDIF
              
              IF(MK(2)<=A2N(2))THEN
              IF(M1(I)==A1N(1)+A2M1(MK(2)))THEN
                  MM=MM+A2M2(MK(2))
                  DO J=1,A2M2(MK(2))
                      S1=[S1,A2S1(SK(2))+M1(I)-2*A1N(1)]
                      A=[A,A2(AK(2):AK(2)+A2S3(SK(2))-1)]
                      AK(2)=AK(2)+A2S3(SK(2))
                      S3=[S3,A2S3(SK(2))]
                      SK(2)=SK(2)+1
                  ENDDO  
                  MK(2)=MK(2)+1
              ENDIF
              ENDIF
              
              IF(MK(5)<=A5N(2))THEN
              IF(M1(I)==2*A1N(1)+A5M1(MK(5)))THEN
                  MM=MM+A5M2(MK(5))
                  DO J=1,A5M2(MK(5))
                      S1=[S1,A5S1(SK(5))+A1N(1)]
                      A=[A,A5(AK(5):AK(5)+A5S3(SK(5))-1)]
                      AK(5)=AK(5)+A5S3(SK(5))
                      S3=[S3,A5S3(SK(5))]
                      SK(5)=SK(5)+1
                  ENDDO  
                  MK(5)=MK(5)+1
              ENDIF
              ENDIF
              IF(MK(3)<=A3N(2))THEN
              IF(M1(I)==A1N(1)+A3M1(MK(3)))THEN
                  MM=MM+A3M2(MK(3))
                  DO J=1,A3M2(MK(3))
                      S1=[S1,A3S1(SK(3))+M1(I)-A1N(1)]
                      A=[A,A3(AK(3):AK(3)+A3S3(SK(3))-1)]
                      AK(3)=AK(3)+A3S3(SK(3))
                      S3=[S3,A3S3(SK(3))]
                      SK(3)=SK(3)+1
                  ENDDO  
                  MK(3)=MK(3)+1
              ENDIF
              ENDIF
              IF(MK(6)<=A6N(2))THEN
              IF(M1(I)==2*A1N(1)+A6M1(MK(6)))THEN
                  MM=MM+A6M2(MK(6))
                  DO J=1,A6M2(MK(6))
                      S1=[S1,A6S1(SK(6))+2*A1N(1)]
                      A=[A,A6(AK(6):AK(6)+A6S3(SK(6))-1)]
                      AK(6)=AK(6)+A6S3(SK(6))
                      S3=[S3,A6S3(SK(6))]
                      SK(6)=SK(6)+1
                  ENDDO  
                  MK(6)=MK(6)+1
              ENDIF
              ENDIF
              
              M2=[M2,MM]
              
          ELSEIF(M1(I)<=4*A1N(1))THEN
              MM=0
              IF(MK(4)<=A4N(2))THEN
              IF(M1(I)==2*A1N(1)+A4M1(MK(4)))THEN
                  MM=MM+A4M2(MK(4))
                  DO J=1,A4M2(MK(4))
                      S1=[S1,A4S1(SK(4))]
                      A=[A,A4(AK(4):AK(4)+A4S3(SK(4))-1)]
                      AK(4)=AK(4)+A4S3(SK(4))
                      S3=[S3,A4S3(SK(4))]
                      SK(4)=SK(4)+1
                  ENDDO  
                  MK(4)=MK(4)+1
              ENDIF
              ENDIF
              
              IF(MK(2)<=A2N(2))THEN
              IF(M1(I)==3*A1N(1)+A2M1(MK(2)))THEN
                  MM=MM+A2M2(MK(2))
                  DO J=1,A2M2(MK(2))
                      S1=[S1,A2S1(SK(2))+4*A1N(1)-M1(I)]
                      A=[A,A2(AK(2):AK(2)+A2S3(SK(2))-1)]
                      AK(2)=AK(2)+A2S3(SK(2))
                      S3=[S3,A2S3(SK(2))]
                      SK(2)=SK(2)+1
                  ENDDO  
                  MK(2)=MK(2)+1
              ENDIF
              ENDIF
              
              IF(MK(5)<=A5N(2))THEN
              IF(M1(I)==2*A1N(1)+A5M1(MK(5)))THEN
                  MM=MM+A5M2(MK(5))
                  DO J=1,A5M2(MK(5))
                      S1=[S1,A5S1(SK(5))+A1N(1)]
                      A=[A,A5(AK(5):AK(5)+A5S3(SK(5))-1)]
                      AK(5)=AK(5)+A5S3(SK(5))
                      S3=[S3,A5S3(SK(5))]
                      SK(5)=SK(5)+1
                  ENDDO  
                  MK(5)=MK(5)+1
              ENDIF
              ENDIF
              IF(MK(3)<=A3N(2))THEN
              IF(M1(I)==3*A1N(1)+A3M1(MK(3)))THEN
                  MM=MM+A3M2(MK(3))
                  DO J=1,A3M2(MK(3))
                      S1=[S1,A3S1(SK(3))+5*A1N(1)-M1(I)]
                      A=[A,A3(AK(3):AK(3)+A3S3(SK(3))-1)]
                      AK(3)=AK(3)+A3S3(SK(3))
                      S3=[S3,A3S3(SK(3))]
                      SK(3)=SK(3)+1
                  ENDDO  
                  MK(3)=MK(3)+1
              ENDIF
              ENDIF
              IF(MK(6)<=A6N(2))THEN
              IF(M1(I)==2*A1N(1)+A6M1(MK(6)))THEN
                  MM=MM+A6M2(MK(6))
                  DO J=1,A6M2(MK(6))
                      S1=[S1,A6S1(SK(6))+2*A1N(1)]
                      A=[A,A6(AK(6):AK(6)+A6S3(SK(6))-1)]
                      AK(6)=AK(6)+A6S3(SK(6))
                      S3=[S3,A6S3(SK(6))]
                      SK(6)=SK(6)+1
                  ENDDO  
                  MK(6)=MK(6)+1
              ENDIF
              ENDIF
              
              M2=[M2,MM]
          ELSEIF(M1(I)<=5*A1N(1))THEN
              MM=0
              IF(MK(2)<=A2N(2))THEN
              IF(M1(I)==3*A1N(1)+A2M1(MK(2)))THEN
                  MM=MM+A2M2(MK(2))
                  DO J=1,A2M2(MK(2))
                      S1=[S1,A2S1(SK(2))]
                      A=[A,A2(AK(2):AK(2)+A2S3(SK(2))-1)]
                      AK(2)=AK(2)+A2S3(SK(2))
                      S3=[S3,A2S3(SK(2))]
                      SK(2)=SK(2)+1
                  ENDDO  
                  MK(2)=MK(2)+1
              ENDIF
              ENDIF
              
              IF(MK(1)<=A1N(2))THEN
              IF(M1(I)==4*A1N(1)+A1M1(MK(1)))THEN
                  MM=MM+A1M2(MK(1))
                  DO J=1,A1M2(MK(1))
                      S1=[S1,A1S1(SK(1))+5*A1N(1)-M1(I)]
                      A=[A,A1(AK(1):AK(1)+A1S3(SK(1))-1)]
                      AK(1)=AK(1)+A1S3(SK(1))
                      S3=[S3,A1S3(SK(1))]
                      SK(1)=SK(1)+1
                  ENDDO  
                  MK(1)=MK(1)+1
              ENDIF
              ENDIF
              
              IF(MK(3)<=A3N(2))THEN
              IF(M1(I)==3*A1N(1)+A3M1(MK(3)))THEN
                  MM=MM+A3M2(MK(3))
                  DO J=1,A3M2(MK(3))
                      S1=[S1,A3S1(SK(3))+A1N(1)]
                      A=[A,A3(AK(3):AK(3)+A3S3(SK(3))-1)]
                      AK(3)=AK(3)+A3S3(SK(3))
                      S3=[S3,A3S3(SK(3))]
                      SK(3)=SK(3)+1
                  ENDDO  
                  MK(3)=MK(3)+1
              ENDIF
              ENDIF
              M2=[M2,MM]
          ELSE
              IF(MK(1)<=A1N(2))THEN
              M2=[M2,A1M2(MK(1))]
              DO J=1,A1M2(MK(1))
                  S1=[S1,A1S1(SK(1))]
                  A=[A,A1(AK(1):AK(1)+A1S3(SK(1))-1)]
                  AK(1)=AK(1)+A1S3(SK(1))
                  S3=[S3,A1S3(SK(1))]
                  SK(1)=SK(1)+1
              ENDDO  
              MK(1)=MK(1)+1
            ENDIF
 
          ENDIF
      ENDDO
      S2=[S2,1]
      DO I=2,SIZE(S3)
          S2=[S2,S2(I-1)+S3(I-1)]
      ENDDO
      
      CALL DIA_SORT(N,M1,M2,S1,S2,S3,A)    
      
      END
      
      !3toDIA
      SUBROUTINE DIAOF3(AI,AJ,A3,N,AN,AM1,AM2,AS1,AS2,AS3,A)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::AN(2),I,J,K,N,S,MM
      real(kind=rk),allocatable::A(:),A3(:)
      integer,allocatable,DIMENSION(:)::AM1,AM2,M,AI,AJ,AS1,AS2,AS3
      

      AN(1)=N
      AN(2)=0
      ! Find diagonals with elements 
      ALLOCATE(M(2*N-1))
      M=0
      S=SIZE(A3)
      DO I=1,S
          K=N+AJ(I)-AI(I)
          IF(M(K)==0)THEN
              M(K)=1
              AN(2)=AN(2)+1
          ENDIF
      ENDDO
      
      I=1
      K=1
      DO WHILE(I<=AN(2).AND.K<2*N)
         IF(M(K)==1)THEN
             MM=0
             DO J=1,SIZE(A3)
                 IF(K==N-AI(J)+AJ(J))THEN
                     MM=MM+1
                     AS1=[AS1,MIN(AI(J),AJ(J))]
                     AS3=[AS3,1]
                     A=[A,A3(J)]
                     AS2=[AS2,SIZE(A)]
                 ENDIF
             ENDDO
             AM2=[AM2,MM]
             AM1=[AM1,K]
             I=I+1
         ENDIF
         K=K+1
      ENDDO
      
      CALL DIA_SORT(AN,AM1,AM2,AS1,AS2,AS3,A)


      END
      
      
      subroutine sort_LINE(a)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer::I,J,K,P,IDX
      INTEGER,ALLOCATABLE::A(:),B(:)

      ALLOCATE(B(SIZE(A)))
      IDX=0
      DO I=1,SIZE(A)
          K=I
          DO J=I+1,SIZE(A)
              IF(A(J)<A(K))K=J
          ENDDO
          IF(K/=I)THEN
              P=A(I)
              A(I)=A(K)
              A(K)=P
          ENDIF
          IF(I==1)THEN
            IDX=1
            B(IDX)=A(1)
          ELSEIF(A(I-1)/=A(I))THEN
              IDX=IDX+1
              B(IDX)=A(I)
          ENDIF
      ENDDO
      A=B(1:IDX)
      end

      ! Sort sparse matrices by position and merge identical terms 
      subroutine sort(a,ai,aj)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::I,J,length,lengthb,R
      real(kind=rk),allocatable::A(:),b(:)
      integer,allocatable::AI(:),AJ(:),bi(:),bj(:)
      
      LENGTH=SIZE(A)
      ALLOCATE(BI(LENGTH),BJ(LENGTH),B(LENGTH))
      IF(LENGTH.NE.SIZE(AI))WRITE(*,*)'LENGTH REEOE'
      IF(LENGTH.NE.SIZE(AJ))WRITE(*,*)'LENGTH REEOE'
      call flush()
      lengthb=0
      do i=1,length
          do j=i+1,length
              call compare(ai(i),aj(i),ai(j),aj(j),r)
              if (r.eq.1)then
                  call change(a,ai,aj,i,j)
              endif
          enddo
          if (lengthb.eq.0)then
            lengthb=lengthb+1
              b(LENGTHB)=a(i)
              bi(LENGTHB)=ai(i)
              bj(LENGTHB)=aj(i)

          else
              call compare(ai(i),aj(i),bi(lengthb),bj(lengthb),r)
              if(r.eq.0)then
                  b(lengthb)=b(lengthb)+a(i)
              elseif(r.eq.1)then
                lengthb=lengthb+1
                  b(lengthb)=a(i)
                  bi(lengthb)=ai(i)
                  bj(lengthb)=aj(i)
              else
                  WRITE(*,*)'sort error'
                  call flush()
              endif
          endif
          
      end do
! Remove zero entries
      length=0
      A=0
      AI=0
      AJ=0
      DO I=1,LENGTHB
          IF(ABS(B(I))>UMIN)THEN
                LENGTH=LENGTH+1
                a(LENGTH)=b(i)
                ai(LENGTH)=bi(i)
                aj(LENGTH)=bj(i)
          endif
      enddo
      B=A(1:LENGTH)
      BI=AI(1:LENGTH)
      BJ=AJ(1:LENGTH)
      A=B
      AI=BI
      AJ=BJ
      DEALLOCATE(BI,BJ,B)
      end
      
      !a>b,r=1;a<b,r=-1;a=b,r=0
      subroutine compare(ai,aj,bi,bj,r)
      integer::AI,AJ,bi,bj
      integer::r
      if (ai>bi)then
          r=1
      elseif(ai<bi)then
          r=-1
      else
          if(aj>bj)then
              r=1
          elseif(aJ<bj)then
              r=-1
          else 
              r=0
          endif
      endif
      end
      
      subroutine change(a,ai,aj,x,y)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::x,y,t
      real(kind=rk)::at
      real(kind=rk),allocatable::A(:)
      integer,allocatable::AI(:),AJ(:)
      
      if (x>size(a).or.y>size(a))then
          WRITE(*,*)error
          call flush()
      endif
      
      t=ai(x)
      ai(x)=ai(y)
      ai(y)=t
      
      t=aj(x)
      aj(x)=aj(y)
      aj(y)=t
      
      at=a(x)
      a(x)=a(y)
      a(y)=at
      
      end
      
      end module DIA_SM