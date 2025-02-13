      module DIA_SM
      USE Matrix_inverse
      
      !use INI
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      REAL::UMIN=1E-10
      INTEGER::NUM0=0
     
      contains


!稀疏矩阵加法,A被结果覆盖
      subroutine DIA_ADDSM(AN,AM1,AM2,AS1,AS2,AS3,A,BN,BM1,
     + BM2,BS1,BS2,BS3,B)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )

      real(kind=rk),allocatable,DIMENSION(:)::A,B,X,Y
      integer,allocatable,DIMENSION(:)::AM1,AM2,BM1,BM2,M1,M2,
     + AS1,AS2,AS3,BS1,BS2,BS3,S1,S2,S3,NS,NL
      INTEGER::AN(2),BN(2),I,J,K,N,NN,II,JJ
      
      IF(AN(1)/=BN(1))THEN
          PRINT*,"ADD ERROR"
      ENDIF
      
      ALLOCATE(Y(AN(1)))
      I=1
      J=1
      K=1
      NI=1
      NJ=1
      DO WHILE(I<=AN(2))
          IF(J>BN(2))THEN
              M1=[M1,AM1(I)]
              M2=[M2,AM2(I)]
             DO II=1,AM2(I)
                 S1=[S1,AS1(NI)]
                 S3=[S3,AS3(NI)]
                 DO JJ=1,AS3(NI)
                     X=[X,A(AS2(NI)+JJ-1)]
                 ENDDO
                 NI=NI+1
             ENDDO
             I=I+1
          ELSE
          
         IF(AM1(I)<BM1(J))THEN
             M1=[M1,AM1(I)]
              M2=[M2,AM2(I)]
             DO II=1,AM2(I)
                 S1=[S1,AS1(NI)]
                 S3=[S3,AS3(NI)]
                 DO JJ=1,AS3(NI)
                     X=[X,A(AS2(NI)+JJ-1)]
                 ENDDO
                 NI=NI+1
             ENDDO
             I=I+1
         ELSEIF(AM1(I)>BM1(J))THEN
             M1=[M1,BM1(J)]
              M2=[M2,BM2(J)]
             DO II=1,BM2(J)
                 S1=[S1,BS1(NJ)]
                 S3=[S3,BS3(NJ)]
                 DO JJ=1,BS3(NJ)
                     X=[X,B(BS2(NJ)+JJ-1)]
                 ENDDO
                 NJ=NJ+1
             ENDDO
             J=J+1
         ELSE
             M1=[M1,AM1(I)]     !AM1(I)=BM1(J),用Y接收整个斜线非0，再重新分组
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
             M2=[M2,N]
             S1=[S1,NS]
             S3=[S3,NL]
             DO II=1,N
                 DO JJ=1,NL(II)
                     X=[X,Y(NS(II)+JJ-1)]
                 ENDDO
             ENDDO  
             I=I+1
             J=J+1
             DEALLOCATE(NS,NL)
             ALLOCATE(NS(0),NL(0))
         ENDIF
      ENDIF
         K=K+1
      ENDDO
      DO WHILE(J<=BN(2))
          M1=[M1,BM1(J)]
          M2=[M2,BM2(J)]
          DO II=1,BM2(J)
              S1=[S1,BS1(NJ)]
              S3=[S3,BS3(NJ)]
              DO JJ=1,BS3(NJ)
                  X=[X,B(BS2(NJ)+JJ-1)]
              ENDDO
              NJ=NJ+1
          ENDDO
          J=J+1
          K=K+1
      ENDDO
      S2=[S2,1]
      DO I=2,SIZE(S3)
          S2=[S2,S2(I-1)+S3(I-1)]
      ENDDO
      AN(2)=K-1
      AM1=M1
      AM2=M2
      AS1=S1
      AS2=S2
      AS3=S3
      A=X
      END
      
      !检查值为0项
      subroutine DIA_SORT(AN,AM1,AM2,AS1,AS2,AS3,A)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )

      real(kind=rk),allocatable,DIMENSION(:)::A,B,X,Y
      integer,allocatable,DIMENSION(:)::AM1,AM2,M1,M2,
     + AS1,AS2,AS3,S1,S2,S3,NS,NL
      INTEGER::AN(2),I,J,K,N,NN,II,JJ
      

      
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
              M1=[M1,AM1(I)]
          ENDIF
          M2=[M2,N]
          S1=[S1,NS]
          S3=[S3,NL]
          DO J=1,N
              DO K=1,NL(J)
                  X=[X,Y(NS(J)+K-1)]
              ENDDO
          ENDDO  
          
          DEALLOCATE(NS,NL)
          ALLOCATE(NS(0),NL(0))
      ENDDO

      S2=[S2,1]
      DO I=2,SIZE(S3)
          S2=[S2,S2(I-1)+S3(I-1)]
      ENDDO
      AN(2)=SIZE(M1)
      AM1=M1
      AM2=M2
      AS1=S1
      AS2=S2
      AS3=S3
      A=X
      END
      
      SUBROUTINE DLINE(Y,LENGTH,N,M,M1)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::N,I,J,LENGTH,NUM
      real(kind=rk),allocatable::Y(:)
      integer,allocatable::M(:),M1(:)
      I=1
      DO WHILE(I<=LENGTH)
          IF(ABS(Y(I))<UMIN)THEN
              I=I+1
          ELSE
              M=[M,I]
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
              M1=[M1,J-NUM]
          ENDIF
      ENDDO
      N=SIZE(M)   
      END
!创建单位矩阵
      SUBROUTINE DIA_ISM(N,AN,AM1,AM2,AS1,AS2,AS3,A)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::N,I,AN(2)
      REAL(KIND=RK)::X=1
      real(kind=rk),allocatable::A(:)
      integer,allocatable,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3
      AN(1)=N
      AN(2)=1
      AM1=[AM1,N]
      AM2=[AM2,1]
      AS1=[AS1,1]
      AS2=[AS2,1]
      AS3=[AS3,N]
      DO I=1,N
          A=[A,X]
      ENDDO
      END
!稀疏矩阵与向量乘法
      SUBROUTINE DOTSM(AN,AM1,AM2,AS1,AS2,AS3,A,B,OUTPUT,LENGTH)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::N,AN(2),LENGTH,I,J,K,NS
      REAL(KIND=RK),DIMENSION(LENGTH)::B,OUTPUT
      real(kind=rk),allocatable::A(:)
      integer,allocatable,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3

      IF(AN(1)/=LENGTH)THEN
          PRINT*,'DOTSM ERROR'
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
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::N,M,AN(2),LENGTH,I,J,K,NS
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
      
      !OPEN(1,FILE='TEST2.TXT')
      !DO I=1,SIZE(B)
      !WRITE(1,*)B(I),X(I),R(I)
      !ENDDO
      !CLOSE(1)
      
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
              !PRINT*,'M IS TOO LARGE'
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
      
      
      SUBROUTINE LEAST_SQUARE(M,A,B,X)!系数矩阵大小（m+1，m）
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
      !对IA求逆
      !CALL LUM(M,IA)
      call dgetrf(M, M, IA, M, ipiv, info)
      if (info /= 0) then
        print *, 'Matrix is singular'
        PAUSE
      end if
      call dgetri(M, IA, M, ipiv, C,M*M, info)
      if (info /= 0) then
        print *, 'Inverse not computed'
        PAUSE
      end if

      
      CALL MDOT(M,M+1,1,AT,B,BB)
      CALL MDOT(M,M,1,IA,BB,X)
      
      
      END
      
      
      SUBROUTINE MDOT(M1,M2,M3,A,B,OUTPUT)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::M1,M2,M3,I,J,K
      REAL(KIND=RK),allocatable,DIMENSION(:,:)::A,B,OUTPUT
      
      IF(SIZE(OUTPUT)==0)THEN
          ALLOCATE(OUTPUT(M1,M3))
      ELSEIF(SIZE(OUTPUT,1)/=M1.OR.SIZE(OUTPUT,2)/=M3)THEN
          PRINT*,'SIZE WRONG'
          PAUSE
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
      !取A的2范数
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
      INTEGER::N,AN(2),LENGTH,I,J,AI
      REAL(KIND=RK),DIMENSION(LENGTH)::B
      real(kind=rk),allocatable::A(:)
      integer,allocatable,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3
      
      IF(AN(1)/=LENGTH)THEN
          PRINT*,'DOTROW ERROR'
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


      

     !! !对称
     !! !!A4
     !! !!A2 A5
     !! !!A1 A3 A6
      !一般化拼接（含判断）
      SUBROUTINE SPLICE33(N,M1,M2,S1,S2,S3,A,
     + A1N,A1M1,A1M2,A1S1,A1S3,A1,
     + A2N,A2M1,A2M2,A2S1,A2S3,A2,
     + A3N,A3M1,A3M2,A3S1,A3S3,A3,
     + A4N,A4M1,A4M2,A4S1,A4S3,A4,
     + A5N,A5M1,A5M2,A5S1,A5S3,A5,
     + A6N,A6M1,A6M2,A6S1,A6S3,A6)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::I,J,K,MK(6),SK(6),AK(6),MM
      REAL(KIND=RK)::C0=0
      INTEGER,DIMENSION(2)::A1N,A2N,A3N,A4N,A5N,A6N,N
      INTEGER,ALLOCATABLE,DIMENSION(:)::M1,M2,S1,S2,S3,A1M1,A1M2,A1S1,
     + A1S3,A2M1,A2M2,A2S1,A2S3,A3M1,A3M2,A3S1,A3S3,A4M1,
     + A4M2,A4S1,A4S3,A5M1,A5M2,A5S1,A5S3,A6M1,A6M2,A6S1,
     + A6S3
      REAL(KIND=RK),ALLOCATABLE,DIMENSION(:)::A,A1,A2,A3,A4,A5,A6
      
      N(1)=3*A1N(1)
      !改为一般化
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
      CALL sort_LINE(M1)!补上去重复
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
      
      
      !A从小到大排序
      SUBROUTINE FINDL(A,X,IS)
      INTEGER,ALLOCATABLE::A(:)
      INTEGER::I,X,IS
      IS=0
      DO I=1,SIZE(A)
          IF(X==A(I))THEN
              IS=I
              EXIT
          ELSEIF(X<A(I))THEN
              EXIT
          ENDIF
      ENDDO
      END
      
      !DIA转三元格式
      SUBROUTINE DIA23(N,M1,M2,S1,S2,S3,A,AI,AJ,A3)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::N(2),I,J,K,S,NS
      real(kind=rk),allocatable::A(:),A3(:)
      integer,allocatable,DIMENSION(:)::M1,M2,S1,S2,S3,AI,AJ
      
      S=SIZE(A)
      ALLOCATE(AI(S),AJ(S),A3(S))
      
      
      NS=1
      DO I=1,N(2)
          DO J=1,M2(I)
              DO K=1,S3(NS)
                  IF(M1(I)<=N(1))THEN
                      AI=[AI,S1(NS)+K-1+N(1)-M1(I)]
                      AJ=[AJ,S1(NS)+K-1]
                  ELSE
                      AI=[AI,S1(NS)+K-1]
                      AJ=[AJ,S1(NS)+K-1-N(1)+M1(I)]
                  ENDIF
                  A3=[A3,A(S2(NS)+K-1)]
              ENDDO
              NS=NS+1
          ENDDO
      ENDDO
      
      
      CALL sort(A3,AI,AJ)
      END
      !三元格式转DIA
      SUBROUTINE DIAOF3(AI,AJ,A3,N,AN,AM1,AM2,AS1,AS2,AS3,A)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::AN(2),I,J,K,N,S,MM
      real(kind=rk),allocatable::A(:),A3(:)
      integer,allocatable,DIMENSION(:)::AM1,AM2,M,AI,AJ,AS1,AS2,AS3
      

      AN(1)=N
      AN(2)=0
      !找到有元素的斜线
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
      
      SUBROUTINE LUSM(AN,AM1,AM2,AS1,AS2,AS3,A,L,LI,LJ,U,UI,UJ)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::AN(2),I,J,K,N,LL,LN,UN
      real(kind=rk),allocatable::A(:),L(:),U(:)
      REAL(KIND=RK)::SUM
      integer,allocatable,DIMENSION(:)::AM1,AM2,AM3,AS1,AS2,AS3,
     + LI,LJ,UI,UJ
      
      ALLOCATE(AM3(SIZE(AM2)))
      AM3(1)=1
      DO I=2,SIZE(AM2)
          AM3(I)=AM3(I-1)+AM2(I-1)
      ENDDO
      
      LN=1
      DO I=1,AN(1)
          DO J=1,AN(1)
              K=1
              DO WHILE(AM1(K)<AN(1)-I+J)
                  K=K+1
                  IF(K>AN(2))EXIT
              ENDDO
              SUM=0
              IF(K<=AN(2))THEN
              IF(AM1(K)==AN(1)-I+J)THEN
              DO II=1,AM2(K)
                  IF(AS1(AM3(K)+II-1)<=MIN(I,J).AND.AS1(AM3(K)+II-1)
     +            +AS3(AM3(K)+II-1)>MIN(I,J))THEN
                      SUM=A(AS2(AM3(K)+II-1)+MIN(I,J)-AS1(AM3(K)+II-1))
                  ENDIF
              ENDDO
              ENDIF
              ENDIF
              UN=1
              !找l对应值
              DO LN=LN,SIZE(L)
                  IF(LI(LN).EQ.I)THEN
                      !找u对应值
                      DO UN=UN,SIZE(U)
                          IF(UI(UN).EQ.LJ(LN))THEN
                              IF(UJ(UN).EQ.J)THEN
                                  SUM=SUM-L(LN)*U(UN)
                              ELSEIF(UJ(UN)>J)THEN
                                  EXIT
                              ENDIF
                          ELSEIF(UI(UN)>LJ(LN))THEN
                              EXIT
                          ENDIF
                      ENDDO
                      UN=UN-1
                  ELSEIF(LI(LN)>I)THEN
                      EXIT
                  ENDIF
              ENDDO
              
              LN=LN-I+1
              IF(LN<1)LN=1

              
              !赋值
              IF(ABS(SUM)>UMIN)THEN
                IF(I>J)THEN            
                  DO UN=UN,SIZE(U)
                      IF((UI(UN).EQ.J).AND.(UJ(UN).EQ.J))THEN
                          L=[L,SUM/U(UN)]
                          LI=[LI,I]
                          LJ=[LJ,J]
                         EXIT
                      ELSEIF(UI(UN)>J)THEN
                          EXIT
                      ENDIF
                  ENDDO
                ELSE
                  U=[U,SUM]
                  UI=[UI,I]
                  UJ=[UJ,J]
                ENDIF  
              ENDIF 
          ENDDO
      ENDDO

      
      END
      
      SUBROUTINE IL_DIA23(N,M1,M2,S1,S2,S3,A,AI,AJ,A3)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::N(2),I,J,K,NS
      real(kind=rk)::C
      real(kind=rk),allocatable::A(:),A3(:)
      integer,allocatable,DIMENSION(:)::M1,M2,S1,S2,S3,AI,AJ

      NS=1
      C=1
      DO I=1,N(2)
          DO J=1,M2(I)
              DO K=1,S3(NS)
                  IF(ABS(A(S2(NS)+K-1))<1E-3)CYCLE
                  IF(M1(I)<N(1))THEN
                      AI=[AI,S1(NS)+K-1+N(1)-M1(I)]
                      AJ=[AJ,S1(NS)+K-1]
                      A3=[A3,-A(S2(NS)+K-1)]
                  ELSEIF(M1(I)==N(1))THEN
                      AI=[AI,S1(NS)+K-1]
                      AJ=[AJ,S1(NS)+K-1]
                      A3=[A3,C]
                  ENDIF
              ENDDO
              NS=NS+1
          ENDDO
      ENDDO
      CALL sort(A3,AI,AJ)
      END
      SUBROUTINE IU_DIA23(N,M1,M2,S1,S2,S3,A,AI,AJ,A3)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::N(2),I,J,K,NS
      real(kind=rk),allocatable::A(:),A3(:)
      integer,allocatable,DIMENSION(:)::M1,M2,S1,S2,S3,AI,AJ

      NS=1
      DO I=1,N(2)
          DO J=1,M2(I)
              DO K=1,S3(NS)
                  IF(ABS(A(S2(NS)+K-1))<1E-3)CYCLE
                  IF(M1(I)>N(1))THEN
                      AI=[AI,S1(NS)+K-1]
                      AJ=[AJ,S1(NS)+K-1-N(1)+M1(I)]
                      A3=[A3,-A(S2(NS)+K-1)]
                  ELSEIF(M1(I)==N(1))THEN
                      AI=[AI,S1(NS)+K-1]
                      AJ=[AJ,S1(NS)+K-1]
                      A3=[A3,1/A(S2(NS)+K-1)]
                  ENDIF
              ENDDO
              NS=NS+1
          ENDDO
      ENDDO
      CALL sort(A3,AI,AJ)
      END
      !求伪逆
      SUBROUTINE DIA_ILU(AN,AM1,AM2,AS1,AS2,AS3,A,BN,BM1,BM2,
     + BS1,BS2,BS3,B)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::i,j,K,N,N1,N2,AN(2),BN(2)
      real(kind=rk),allocatable,DIMENSION(:)::A,B,L,U,A3,IM
      REAL(KIND=RK)::SUM
      integer,allocatable,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3,BM1,BM2,BS1,
     + BS2,BS3,AI,AJ,LI,LJ,UI,UJ,IMI,IMJ
      integer(kind=8)  ::count1,count2,count_rate,count_max
      real timespend
      logical :: error_flag = .false.
      INTEGER::RETRY_COUNT,MAX_RETRIES=10,E
      
      call system_clock(count1,count_rate)
      
      NS=1
      DO I=1,AN(2)
          IF(AM1(I)>=AN(1))EXIT
          NS=NS+AM2(I)
      ENDDO
      
      BM1=AM1(1:I-1)
      BM2=AM2(1:I-1)
      BS1=AS1(1:NS-1)
      BS2=AS2(1:NS-1)
      BS3=AS3(1:NS-1)
      B=A(1:AS2(NS)-1)
      CALL IL_DIA23([AN(1),I-1],BM1,BM2,BS1,BS2,BS3,B,LI,LJ,L)
      
      BM1=AM1(I:AN(2))
      BM2=AM2(I:AN(2))
      B=A(AS2(NS):SIZE(A))      
      BS1=AS1(NS:SIZE(AS1))
      BS2=AS2(NS:SIZE(AS2))-AS2(NS)+1
      BS3=AS3(NS:SIZE(AS3))
      CALL IU_DIA23([AN(1),AN(2)-I+1],BM1,BM2,BS1,BS2,BS3,B,UI,UJ,U)
      
      
      CALL ISM(IM,IMI,IMJ,AN(1))
      CALL ADDSM(L,LI,LJ,IM,IMI,IMJ)

      call system_clock(count2,count_rate,count_max)
      timespend=(count2-count1)/real(count_rate)
      write(*,*) '分解: ',timespend 
      
      call system_clock(count1,count_rate)
      CALL SMDOTSM2(AN(1),U,UI,UJ,L,LI,LJ,A3,AI,AJ)
      
      call system_clock(count2,count_rate,count_max)
      timespend=(count2-count1)/real(count_rate)
      write(*,*) '乘法: ',timespend 
      !!CALL DIA23(AN,AM1,AM2,AS1,AS2,AS3,A,BI,BJ,B3)
      !!CALL SMDOTSM(375,A3,AI,AJ,B3,BI,BJ,C3,CI,CJ)
      call system_clock(count1,count_rate)
      CALL DIAOF3(AI,AJ,A3,AN(1),BN,BM1,BM2,BS1,BS2,BS3,B)
      
      call system_clock(count2,count_rate,count_max)
      timespend=(count2-count1)/real(count_rate)
      write(*,*) '格式: ',timespend 


      
      END
      
      
      
      !利用lu分解对稀疏矩阵a求逆，方阵，维度n
      SUBROUTINE ILUSM(AN,AM1,AM2,AS1,AS2,AS3,A)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::i,j,N,N1,N2,AN(2)
      real(kind=rk),allocatable,DIMENSION(:)::A,L,U,IL,IU,A3,IM,B3,C3
      REAL(KIND=RK)::SUM
      integer,allocatable,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3,AI,AJ,LI,LJ,
     +    UI,UJ,ILI,ILJ,IUI,IUJ,IMI,IMJ,BI,BJ,CI,CJ

      CALL LUSM(AN,AM1,AM2,AS1,AS2,AS3,A,L,LI,LJ,U,UI,UJ)
      
      !为il补充对角元素
      CALL ISM(IM,IMI,IMJ,N)
      CALL ADDSM(L,LI,LJ,IM,IMI,IMJ)

      N=AN(1)
      N1=1
      !对u求逆iu
      DO I=1,N
          DO J=I,N
              N2=1
              IF(I.EQ.J)THEN 
                  SUM=1 
              ELSE
                  SUM=0
              ENDIF
              DO N1=N1,SIZE(IU)
                  IF(IUI(N1).EQ.I)THEN
                      DO N2=N2,SIZE(U)
                          IF(UI(N2).EQ.IUJ(N1))THEN
                              IF(UJ(N2).EQ.J)THEN
                                  SUM=SUM-IU(N1)*U(N2)
                              ELSEIF(UJ(N2)>J)THEN
                                  EXIT
                              ENDIF
                          ELSEIF(UI(N2)>IUJ(N1))THEN
                              EXIT
                          ENDIF
                      ENDDO
                      N2=N2-1
                  ELSEIF(IUI(N1)>I)THEN
                      EXIT
                  ENDIF
              ENDDO
              
              N1=N1-J
              IF(N1<1)N1=1
                          
              IF(ABS(SUM)>UMIN)THEN
                  DO N2=N2,SIZE(U)
                      IF((UI(N2).EQ.J).AND.(UJ(N2).EQ.J))THEN
                          IU=[IU,SUM/U(N2)]
                          IUI=[IUI,I]
                          IUJ=[IUJ,J]
                          EXIT
                      ENDIF
                  ENDDO
              ENDIF
          ENDDO
      ENDDO
      DEALLOCATE(U,UI,UJ)
      !对l求逆il
      DO I=N,2,-1
          DO J=I-1,1,-1
              SUM=0
              N2=SIZE(L)
              DO N2=N2,1,-1
                  IF((LI(N2).EQ.I).AND.(LJ(N2).EQ.J))THEN
                      SUM=SUM-L(N2)
                      EXIT
                  ELSEIF((LJ(N2)<J).AND.(LI(N2)<=I))THEN
                      EXIT
                  ENDIF
              ENDDO
              N2=N2+1
              
              DO N1=1,SIZE(IL)
                  IF(ILI(N1).EQ.I)THEN
                      DO N2=N2,1,-1
                          IF(LI(N2).EQ.ILJ(N1))THEN
                              IF(LJ(N2).EQ.J)THEN
                                  SUM=SUM-IL(N1)*L(N2)
                              ELSEIF(LJ(N2)<J)THEN
                                  EXIT
                              ENDIF
                          ELSEIF(LI(N2)<ILJ(N1))THEN
                              EXIT
                          ENDIF
                      ENDDO
                      N2=N2+1
                  ELSEIF(LI(N1)>I)THEN
                      EXIT
                  ENDIF
              ENDDO
                          
              IF(ABS(SUM)>UMIN)THEN
                  IL=[SUM,IL]
                  ILI=[I,ILI]
                  ILJ=[J,ILJ]
              ENDIF
          ENDDO
      ENDDO
      IF(SIZE(L).NE.0)DEALLOCATE(L,LI,LJ)
      
      !为il补充对角元素
      CALL ISM(IM,IMI,IMJ,N)
      CALL ADDSM(IL,ILI,ILJ,IM,IMI,IMJ)

      CALL SMDOTSM(N,IU,IUI,IUJ,IL,ILI,ILJ,A3,AI,AJ)
      
      !!CALL DIA23(AN,AM1,AM2,AS1,AS2,AS3,A,BI,BJ,B3)
      !!CALL SMDOTSM(375,A3,AI,AJ,B3,BI,BJ,C3,CI,CJ)
      
      DEALLOCATE(IM,IMI,IMJ,AM1,AM2,AS1,AS2,AS3,A)
      ALLOCATE(AM1(0),AM2(0),AS1(0),AS2(0),AS3(0),A(0))
      CALL DIAOF3(AI,AJ,A3,N,AN,AM1,AM2,AS1,AS2,AS3,A)

      
      END
      !稀疏矩阵与稀疏矩阵点乘，ab均为方阵，输出c
      SUBROUTINE SMDOTSM(N,A,AI,AJ,B,BI,BJ,C,CI,CJ)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::i,j,N,N1,N2
      real(kind=rk),allocatable::A(:),B(:),C(:)
      REAL(KIND=RK)::SUM
      integer,allocatable::AI(:),AJ(:),BI(:),BJ(:),CI(:),CJ(:)
      
      N1=1
      DO I=1,N
          DO J=1,N
              SUM=0
              N2=1
              DO N1=N1,SIZE(A)
                  IF(AI(N1).EQ.I)THEN
                      DO N2=N2,SIZE(B)
                          IF(BI(N2).EQ.AJ(N1))THEN
                              IF(BJ(N2).EQ.J)THEN
                                  SUM=SUM+A(N1)*B(N2)
                              ELSEIF(BJ(N2)>J)THEN
                                  EXIT
                              ENDIF
                          ELSEIF(BI(N2)>AJ(N1))THEN
                              EXIT
                          ENDIF
                      ENDDO
                      IF(N2>1)N2=N2-1
                  ELSEIF(AI(N1)>I)THEN
                      EXIT
                  ENDIF
              ENDDO
              
              N1=N1-N
              IF(N1<1)N1=1
              
              IF(ABS(SUM)>UMIN)THEN
                  C=[C,SUM]
                  CI=[CI,I]
                  CJ=[CJ,J]
              ENDIF
              
          ENDDO
      ENDDO
      end
      
      SUBROUTINE SMDOTSM2(N,A,AI,AJ,B,BI,BJ,C,CI,CJ)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::i,j,N,N1,N2
      real(kind=rk),allocatable::A(:),B(:),C(:)
      REAL(KIND=RK)::SUM,S
      integer,allocatable::AI(:),AJ(:),BI(:),BJ(:),CI(:),CJ(:)
      
      N1=SIZE(A)
      N2=SIZE(B)
      DO I=1,N1
          DO J=1,N2
              IF(AJ(I)==BI(J))THEN
                  S=A(I)*B(J)
                  IF(ABS(S)>1E-3)THEN
                  CI=[CI,AI(I)]
                  CJ=[CJ,BJ(J)]
                  C=[C,A(I)*B(J)]
                  ENDIF
              ELSEIF(AJ(I)<BI(J))THEN
                  EXIT
              ENDIF
          ENDDO
      ENDDO
      
      CALL sort(C,CI,CJ)
      
      end
      
      
      
      !稀疏矩阵加法,A被结果覆盖
      subroutine ADDSM(A,AI,AJ,B,BI,BJ)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )

      real(kind=rk),allocatable::A(:),b(:)
      integer,allocatable::AI(:),AJ(:),bi(:),bj(:)
      
      A=[A,B]
      AI=[AI,BI]
      AJ=[AJ,BJ]
      CALL sort(A,AI,AJ)
      
      end
      
!创建单位矩阵
      SUBROUTINE ISM(A,AI,AJ,N)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER::N,I
      REAL(KIND=RK)::X=1
      real(kind=rk),allocatable::A(:)
      integer,allocatable::AI(:),AJ(:)
      DO I=1,N
          AI=[AI,I]
          AJ=[AJ,I]
          A=[A,X]
      ENDDO
      END
      
      
      subroutine sort_LINE(a)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer::I,J,K,P
      INTEGER,ALLOCATABLE::A(:),B(:)
      
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
              B=[B,A(1)]
          ELSEIF(A(I-1)/=A(I))THEN
              B=[B,A(I)]
          ENDIF
      ENDDO
      A=B
      end
      
      
      !稀疏矩阵按位置排序并合并相同项
      subroutine sort(a,ai,aj)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::I,J,length,lengthb,R
      real(kind=rk),allocatable::A(:),b(:)
      integer,allocatable::AI(:),AJ(:),bi(:),bj(:)
      
      LENGTH=SIZE(A)
      IF(LENGTH.NE.SIZE(AI))PRINT*,'LENGTH REEOE'
      IF(LENGTH.NE.SIZE(AJ))PRINT*,'LENGTH REEOE'
      lengthb=0
      do i=1,length
          do j=i+1,length
              call compare(ai(i),aj(i),ai(j),aj(j),r)
              if (r.eq.1)then
                  call change(a,ai,aj,i,j)
              endif
          enddo
          if (lengthb.eq.0)then
              b=[b,a(i)]
              bi=[bi,ai(i)]
              bj=[bj,aj(i)]
              lengthb=lengthb+1
          else
              call compare(ai(i),aj(i),bi(lengthb),bj(lengthb),r)
              if(r.eq.0)then
                  b(lengthb)=b(lengthb)+a(i)
              elseif(r.eq.1)then
                  b=[b,a(i)]
                  bi=[bi,ai(i)]
                  bj=[bj,aj(i)]
                  lengthb=lengthb+1
              else
                  print*,'sort error'
              endif
          endif
          
      end do
!去零项
      length=0
      A=[0]
      AI=[O]
      AJ=[0]      
      DO I=1,LENGTHB
          IF(ABS(B(I))>UMIN)THEN
              if(length.eq.0)then
                  a=[b(i)]
                  ai=[bi(i)]
                  aj=[bj(i)]
                  LENGTH=1
              else
                  a=[a,b(i)]
                  ai=[ai,bi(i)]
                  aj=[aj,bj(i)]
              endif
          endif
      enddo
      end
      
! 矩阵位置比较，a>b,r=1;a<b,r=-1;a=b,r=0
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
      
      
      
!交换x，y处
      subroutine change(a,ai,aj,x,y)
      integer, parameter :: rk = kind ( 1.0D+00 )
      integer::x,y,t
      real(kind=rk)::at
      real(kind=rk),allocatable::A(:)
      integer,allocatable::AI(:),AJ(:)
      
      if (x>size(a).or.y>size(a))then
          print*,error
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