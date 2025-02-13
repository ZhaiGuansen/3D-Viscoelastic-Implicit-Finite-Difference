      module PDC          !physical boundary condition
      USE DIA_SM
      use INI
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)   
      contains
      !1为分界上方，2为分界下方
      SUBROUTINE PDC_A(N,NX,NY,NZ,AN,AM1,AM2,AS1,AS2,AS3,A,M1N,M1M1,
     + M1M2,M1S1,M1S2,M1S3,M1,M2N,M2M1,M2M2,M2S1,M2S2,M2S3,M2)
      integer, parameter :: rk = kind ( 1.0D+00 )
      INTEGER,DIMENSION(2)::AN,M1N,M2N
      INTEGER,ALLOCATABLE,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3,M1M1,
     + M1M2,M1S1,M1S2,M1S3,M2M1,M2M2,M2S1,M2S2,M2S3
      REAL(KIND=RK),ALLOCATABLE::A(:),M1(:),M2(:)
      INTEGER::NX,NY,NZ,X,Y,Z,I,J,K
      INTEGER::POSITION,PDC(NX-2,NY-2,NZ-2)
      REAL::DX,DY,DZ,DT
      REAL,DIMENSION(NX,NY,NZ)::VP,VS,RHO
      CHARACTER(LEN=30)::FILENAME
      CHARACTER(LEN=10)::CTEMP
      !POSITION(N0,X0,Y0,Z0)=X0+(Y0-1)*NX+(Z0-1)*NX*NY+(N0-1)*NX*NY*NZ
      
      WRITE(CTEMP,'(I6)')N
      FILENAME='.\BASE\BASE'//TRIM(ADJUSTL(CTEMP))//'.DAT'
      OPEN(N,FILE=FILENAME)
      READ(N,*)X,Y,Z
      READ(N,*)X,Y
      READ(N,*)DX,DY
      READ(N,*)DZ,DT
      READ(N,*)Z
      do k=1,NZ-2
          do j=1,NY-2
              do i=1,NX-2
                  READ(N,*)PDC(I,J,K)
              ENDDO
          ENDDO
      ENDDO
      do k=1,NZ
          do j=1,NY
              do i=1,NX
                  READ(N,*)VP(I,J,K)
                  READ(N,*)VS(I,J,K)
                  READ(N,*)RHO(I,J,K)
              ENDDO
          ENDDO
      ENDDO
      CLOSE(N)
      !确定pdc信息,反射界面方向：二进制编码zyx
      !三方向均不均匀时不再采用该方法
      do k=1,NZ-2
          do j=1,NY-2
              do i=1,NX-2
                  IF(PDC(I,J,K)/=0)THEN
                      CALL PDC_0(NX,NY,NZ,I+1,J+1,K+1,AN,AM1,AM2,AS1,
     +                  AS2,AS3,A,1.)
                      CALL PDC_0(NX,NY,NZ,I+1,J+1,K+1,M2N,M2M1,M2M2,
     +                 M2S1,M2S2,M2S3,M2,-1.)
                  ENDIF
                  SELECT CASE(PDC(I,J,K))
                  CASE(1)
                      CALL PDC_1(NX,NY,NZ,I+1,J+1,K+1,M1N,M1M1,M1M2,
     + M1S1,M1S2,M1S3,M1,DX,DY,DZ,DT,VP,VS,RHO,1)
                      
                  CASE(2)
                      CALL PDC_1(NX,NY,NZ,I+1,J+1,K+1,M1N,M1M1,M1M2,
     + M1S1,M1S2,M1S3,M1,DX,DY,DZ,DT,VP,VS,RHO,2)
                      
                  CASE(4)
                      CALL PDC_1(NX,NY,NZ,I+1,J+1,K+1,M1N,M1M1,M1M2,
     + M1S1,M1S2,M1S3,M1,DX,DY,DZ,DT,VP,VS,RHO,3)
                  ENDSELECT
                  
              ENDDO
          ENDDO
      ENDDO
      

      
      
      
      END
      
      SUBROUTINE PDC_0(NX,NY,NZ,X,Y,Z,AN,AM1,AM2,AS1,AS2,AS3,A,C)
      INTEGER,PARAMETER::RK=KIND(1.0D+00)
      INTEGER,DIMENSION(2)::AN
      INTEGER,ALLOCATABLE,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3
      REAL(KIND=RK),ALLOCATABLE::A(:)
      REAL::C
      INTEGER::NX,NY,NZ,X,Y,Z,I,J,K,ROW1,ROW2,ROW3,NS
      INTEGER::POSITION
      POSITION(N0,X0,Y0,Z0)=X0+(Y0-1)*NX+(Z0-1)*NX*NY+(N0-1)*NX*NY*NZ
      
      ROW1=POSITION(1,X,Y,Z)
      ROW2=POSITION(2,X,Y,Z)
      ROW3=POSITION(3,X,Y,Z)
      !这三行除对角线为c外其余全为0
      
      NS=1
      DO I=1,AN(2)
          DO J=1,AM2(I)
              IF(AM1(I)<AN(1))THEN
                  K=ROW1-AS1(NS)-AN(1)+AM1(I)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
                  K=ROW2-AS1(NS)-AN(1)+AM1(I)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
                  K=ROW3-AS1(NS)-AN(1)+AM1(I)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
              ELSEIF(AM1(I)>AN(1))THEN
                  K=ROW1-AS1(NS)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
                  K=ROW2-AS1(NS)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
                  K=ROW3-AS1(NS)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
              ELSE
                  K=ROW1-AS1(NS)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=C
                  K=ROW2-AS1(NS)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=C
                  K=ROW3-AS1(NS)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=C
              ENDIF
              NS=NS+1
          ENDDO
      ENDDO   
      END

      SUBROUTINE PDC_1(NX,NY,NZ,X,Y,Z,AN,AM1,AM2,AS1,AS2,AS3,A,
     + DX,DY,DZ,DT,VP,VS,RHO,D)
      INTEGER,PARAMETER::RK=KIND(1.0D+00)
      INTEGER,DIMENSION(2)::AN,BN
      INTEGER,ALLOCATABLE,DIMENSION(:)::AM1,AM2,AS1,AS2,AS3,
     + BM1,BM2,BS1,BS2,BS3,BI,BJ
      REAL(KIND=RK),ALLOCATABLE::A(:),B(:),B3(:)
      INTEGER::NX,NY,NZ,X,Y,Z,I,J,K,LA,ROW1,ROW2,ROW3,D,NS
      INTEGER::POSITION
      REAL::DX,DY,DZ,DT,L1,L2,M1,M2,P1,P2,X2,Y2,Z2,T2,PP,C(4)
      REAL,DIMENSION(NX,NY,NZ)::VP,VS,RHO
      POSITION(N0,X0,Y0,Z0)=X0+(Y0-1)*NX+(Z0-1)*NX*NY+(N0-1)*NX*NY*NZ
      
      SELECT CASE(D)
      CASE(1)
          P2=RHO(X-1,Y,Z)
          L2=(VP(X-1,Y,Z)**2-2*VS(X-1,Y,Z)**2)*P2
          M2=VS(X-1,Y,Z)**2*P2
      CASE(2)
          P2=RHO(X,Y-1,Z)
          L2=(VP(X,Y-1,Z)**2-2*VS(X,Y-1,Z)**2)*P2
          M2=VS(X,Y-1,Z)**2*P2
      CASE(3)
          P2=RHO(X,Y,Z-1)
          L2=(VP(X,Y,Z-1)**2-2*VS(X,Y,Z-1)**2)*P2
          M2=VS(X,Y,Z-1)**2*P2
      ENDSELECT

      
      P1=RHO(X,Y,Z)
      L1=(VP(X,Y,Z)**2-2*VS(X,Y,Z)**2)*P1
      M1=VS(X,Y,Z)**2*P1
      T2=DT**2
      X2=T2/DX**2
      Y2=T2/DY**2
      Z2=T2/DZ**2
      PP=P1+P2
      ROW1=POSITION(1,X,Y,Z)
      ROW2=POSITION(2,X,Y,Z)
      ROW3=POSITION(3,X,Y,Z)
!先将所在行归零，再利用稀疏矩阵加法修改矩阵信息
      
      NS=1
      DO I=1,AN(2)
          DO J=1,AM2(I)
              IF(AM1(I)<AN(1))THEN
                  K=ROW1-AS1(NS)-AN(1)+AM1(I)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
                  K=ROW2-AS1(NS)-AN(1)+AM1(I)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
                  K=ROW3-AS1(NS)-AN(1)+AM1(I)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
              ELSE
                  K=ROW1-AS1(NS)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
                  K=ROW2-AS1(NS)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
                  K=ROW3-AS1(NS)+1
                  IF(K>0.AND.K<=AS3(NS))A(AS2(NS)+K-1)=0
              ENDIF
              NS=NS+1
          ENDDO
      ENDDO  

      
      
      
      ALLOCATE(B3(57),BI(57),BJ(57))
      DO I=1,19
          BI(I)=ROW1
      ENDDO
      DO I=20,38
          BI(I)=ROW2
      ENDDO
      DO I=39,57
          BI(I)=ROW3
      ENDDO
      !临时计算数据
      C(1)=(L1+L2+2*M1+2*M2)/PP*X2
      C(2)=(M1+M2)/PP*Y2
      C(3)=(M1+M2)/PP*Z2
      
      BJ(1)=ROW1
      B3(1)=2-2*(C(1)+C(2)+C(3))
      BJ(2)=ROW1+1
      B3(2)=C(1)
      BJ(3)=ROW1-1
      B3(3)=C(1)      
      BJ(4)=ROW1+NX
      B3(4)=C(2)
      BJ(5)=ROW1-NX
      B3(5)=C(2)    
      BJ(6)=ROW1+NX*NY
      B3(6)=C(3)   
      BJ(7)=ROW1-NX*NY
      B3(7)=C(3)
      SELECT CASE(D)
      CASE(1)
          B3(2)=2*(L1+2*M1)*X2/PP
          B3(3)=2*(L2+2*M2)*X2/PP
      CASE(2)
          B3(4)=2*M1*Y2/PP
          B3(5)=2*M2*Y2/PP
      CASE(3)
          B3(6)=2*M1*Z2/PP 
          B3(7)=2*M2*Z2/PP 
      ENDSELECT
      
      C(1)=(L1+M1)*T2/PP/4/DX
      C(2)=(L2+M2)*T2/PP/4/DX
      C(3)=(L1-L2)/PP
      C(4)=(M1-M2)/PP
      
      BJ(8)=POSITION(2,X+1,Y+1,Z)
      BJ(9)=POSITION(2,X+1,Y-1,Z)
      BJ(10)=POSITION(2,X-1,Y+1,Z)
      BJ(11)=POSITION(2,X-1,Y-1,Z)
      BJ(12)=POSITION(3,X+1,Y,Z+1)
      BJ(13)=POSITION(3,X+1,Y,Z-1)
      BJ(14)=POSITION(3,X-1,Y,Z+1)
      BJ(15)=POSITION(3,X-1,Y,Z-1)
      B3(8)=2*C(1)/DY
      B3(9)=-B3(8)
      B3(10)=-2*C(2)/DY
      B3(11)=-B3(10)
      B3(12)=2*C(1)/DZ
      B3(13)=-B3(12)
      B3(14)=-2*C(2)/DZ
      B3(15)=-B3(14)
      
      BJ(16)=POSITION(2,X+1,Y,Z)
      BJ(17)=POSITION(2,X-1,Y,Z)
      BJ(18)=POSITION(3,X+1,Y,Z)
      BJ(19)=POSITION(3,X-1,Y,Z)

      B3(16)=0
      B3(17)=0
      B3(18)=0
      B3(19)=0
      
      SELECT CASE(D)
      CASE(1)
          BJ(16)=POSITION(2,X,Y+1,Z)
          BJ(17)=POSITION(2,X,Y-1,Z)
          BJ(18)=POSITION(3,X,Y,Z+1)
          BJ(19)=POSITION(3,X,Y,Z-1)
          C(4)=(2*(C(2)-C(1))+C(3)*X2)
          B3(16)=C(4)/DY
          B3(17)=-B3(16)
          B3(18)=C(4)/DZ
          B3(19)=-B3(18)
      CASE(2)
          B3(12)=(C(1)+C(2))/DZ
          B3(13)=-B3(12)
          B3(14)=-B3(12)
          B3(15)=B3(12)

          B3(16)=2*(C(2)-C(1))/DY+C(4)*Y2/DX
          B3(17)=-B3(16)
      CASE(3)
          B3(8)=(C(1)+C(2))/DY
          B3(9)=-B3(8)
          B3(10)=-B3(8)
          B3(11)=B3(8)
          
          B3(18)=2*(C(2)-C(1))/DZ+C(4)*Z2/DX
          B3(19)=-B3(18)
      ENDSELECT

      
      
      !Y
      
      C(1)=(M1+M2)/PP*X2
      C(2)=(L1+L2+2*M1+2*M2)/PP*Y2
      C(3)=(M1+M2)/PP*Z2
      
      BJ(20)=ROW2
      B3(20)=2-2*(C(1)+C(2)+C(3))
      BJ(21)=ROW2+1
      B3(21)=C(1)
      BJ(22)=ROW2-1
      B3(22)=C(1)      
      BJ(23)=ROW2+NX
      B3(23)=C(2)
      BJ(24)=ROW2-NX
      B3(24)=C(2)    
      BJ(25)=ROW2+NX*NY
      B3(25)=C(3)   
      BJ(26)=ROW2-NX*NY
      B3(26)=C(3)
      SELECT CASE(D)
      CASE(1)
          B3(21)=2*M1*X2/PP
          B3(22)=2*M2*X2/PP
      CASE(2)
          B3(23)=2*(L1+2*M1)*Y2/PP
          B3(24)=2*(L2+2*M2)*Y2/PP
      CASE(3)
          B3(25)=2*M1*Z2/PP 
          B3(26)=2*M2*Z2/PP 
      ENDSELECT
      
      C(1)=(L1+M1)*T2/PP/4/DY
      C(2)=(L2+M2)*T2/PP/4/DY
      C(3)=(L1-L2)/PP
      C(4)=(M1-M2)/PP
      
      BJ(27)=POSITION(1,X+1,Y+1,Z)
      BJ(28)=POSITION(1,X-1,Y+1,Z)
      BJ(29)=POSITION(1,X+1,Y-1,Z)
      BJ(30)=POSITION(1,X-1,Y-1,Z)
      BJ(31)=POSITION(3,X,Y+1,Z+1)
      BJ(32)=POSITION(3,X,Y+1,Z-1)
      BJ(33)=POSITION(3,X,Y-1,Z+1)
      BJ(34)=POSITION(3,X,Y-1,Z-1)
      B3(27)=2*C(1)/DX
      B3(28)=-B3(27)
      B3(29)=-2*C(2)/DX
      B3(30)=-B3(29)
      B3(31)=2*C(1)/DZ
      B3(32)=-B3(31)
      B3(33)=-2*C(2)/DZ
      B3(34)=-B3(33)
      
      BJ(35)=POSITION(1,X,Y+1,Z)
      BJ(36)=POSITION(1,X,Y-1,Z)
      BJ(37)=POSITION(3,X,Y+1,Z)
      BJ(38)=POSITION(3,X,Y-1,Z)

      B3(35)=0
      B3(36)=0
      B3(37)=0
      B3(38)=0
      
      SELECT CASE(D)
      CASE(1)
          B3(31)=(C(1)+C(2))/DZ
          B3(32)=-B3(31)
          B3(33)=-B3(31)
          B3(34)=B3(31)

          B3(35)=2*(C(2)-C(1))/DX+C(4)*X2/DY
          B3(36)=-B3(35)
      CASE(2)
          BJ(35)=POSITION(1,X+1,Y,Z)
          BJ(36)=POSITION(1,X-1,Y,Z)
          BJ(37)=POSITION(3,X,Y,Z+1)
          BJ(38)=POSITION(3,X,Y,Z-1)
          C(4)=(2*(C(2)-C(1))+C(3)*Y2)
          B3(35)=C(4)/DX
          B3(36)=-B3(35)
          B3(37)=C(4)/DZ
          B3(38)=-B3(37)
      CASE(3)
          B3(27)=(C(1)+C(2))/DX
          B3(28)=-B3(27)
          B3(29)=-B3(27)
          B3(30)=B3(27)
          
          B3(37)=2*(C(2)-C(1))/DZ+C(4)*Z2/DY
          B3(38)=-B3(37)
      ENDSELECT
      
      

      
      
      !Z
      
      C(1)=(M1+M2)/PP*X2
      C(2)=(M1+M2)/PP*Y2
      C(3)=(L1+L2+2*M1+2*M2)/PP*Z2
      
      BJ(39)=ROW3
      B3(39)=2-2*(C(1)+C(2)+C(3))
      BJ(40)=ROW3+1
      B3(40)=C(1)
      BJ(41)=ROW3-1
      B3(41)=C(1)      
      BJ(42)=ROW3+NX
      B3(42)=C(2)
      BJ(43)=ROW3-NX
      B3(43)=C(2)    
      BJ(44)=ROW3+NX*NY
      B3(44)=C(3)   
      BJ(45)=ROW3-NX*NY
      B3(45)=C(3)
      SELECT CASE(D)
      CASE(1)
          B3(40)=2*M1*X2/PP
          B3(41)=2*M2*X2/PP
      CASE(2)
          B3(42)=2*M1*Y2/PP
          B3(43)=2*M2*Y2/PP
      CASE(3)
          B3(44)=2*(L1+2*M1)*Z2/PP 
          B3(45)=2*(L2+2*M2)*Z2/PP 
      ENDSELECT
      
      C(1)=(L1+M1)*T2/PP/4/DZ
      C(2)=(L2+M2)*T2/PP/4/DZ
      C(3)=(L1-L2)/PP
      C(4)=(M1-M2)/PP
      
      BJ(46)=POSITION(1,X+1,Y,Z+1)
      BJ(47)=POSITION(1,X-1,Y,Z+1)
      BJ(48)=POSITION(1,X+1,Y,Z-1)
      BJ(49)=POSITION(1,X-1,Y,Z-1)
      BJ(50)=POSITION(2,X,Y+1,Z+1)
      BJ(51)=POSITION(2,X,Y-1,Z+1)
      BJ(52)=POSITION(2,X,Y+1,Z-1)
      BJ(53)=POSITION(2,X,Y-1,Z-1)
      B3(46)=2*C(1)/DX
      B3(47)=-B3(46)
      B3(48)=-2*C(2)/DX
      B3(49)=-B3(48)
      B3(50)=2*C(1)/DY
      B3(51)=-B3(50)
      B3(52)=-2*C(2)/DY
      B3(53)=-B3(52)
      
      BJ(54)=POSITION(1,X,Y,Z+1)
      BJ(55)=POSITION(1,X,Y,Z-1)
      BJ(56)=POSITION(2,X,Y,Z+1)
      BJ(57)=POSITION(2,X,Y,Z-1)

      B3(54)=0
      B3(55)=0
      B3(56)=0
      B3(57)=0
      
      SELECT CASE(D)
      CASE(1)
          B3(50)=(C(1)+C(2))/DY
          B3(51)=-B3(50)
          B3(52)=-B3(50)
          B3(53)=B3(50)

          B3(54)=2*(C(2)-C(1))/DX+C(4)*X2/DZ
          B3(55)=-B3(54)
      CASE(2)
          B3(46)=(C(1)+C(2))/DX
          B3(47)=-B3(46)
          B3(48)=-B3(46)
          B3(49)=B3(46)
          
          B3(56)=2*(C(2)-C(1))/DY+C(4)*Y2/DZ
          B3(57)=-B3(56)
      CASE(3)
          BJ(54)=POSITION(1,X+1,Y,Z)
          BJ(55)=POSITION(1,X-1,Y,Z)
          BJ(56)=POSITION(2,X,Y+1,Z)
          BJ(57)=POSITION(2,X,Y-1,Z)
          C(4)=(2*(C(2)-C(1))+C(3)*Z2)
          B3(54)=C(4)/DX
          B3(55)=-B3(54)
          B3(56)=C(4)/DY
          B3(57)=-B3(56)
          
      ENDSELECT
      
      CALL DIAOF3(BI,BJ,B3,AN(1),BN,BM1,BM2,BS1,BS2,BS3,B)
      CALL DIA_ADDSM(AN,AM1,AM2,AS1,AS2,AS3,A,BN,BM1,
     + BM2,BS1,BS2,BS3,B)

      END
      
      
      
      
      
      
      
      END MODULE PDC