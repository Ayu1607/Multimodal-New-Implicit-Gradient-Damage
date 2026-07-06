      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)     
!     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      PARAMETER (ZERO= 0.D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,SIX=6.0D0,
     + NEWTON=10,TOLER=1.D-6, NGP = 8, NEL = 75601 ,HALF=0.5D0
     1, FOURTH=.25D0)   


      PARAMETER (YSTRESS_MIN = 2.D0, MAXITER = 200, TOL = 1.D-3, 
     +       PHI_MAX = 0.99D0)
      ! KEEPING MINIMUM YIELDSTRESS SO HIGH THAT PLASTIC PHASE DOES NOT EXIST.  
!     
      DIMENSION STRESS(NTENS),STRAN(NTENS),DSTRAN(NTENS),STRESS0(NTENS),
     +  DDSDDE(NTENS,NTENS),EELAS(NTENS),PROPS(*),EPLAS(NTENS),
     +  FLOW(NTENS),STATEV(NSTATV),TIME(2) 
!     
      DOUBLE PRECISION P_STR(6),P(6,6),Z(6,6), ESTRESS(6),M(6,6),M_INV(6,6),
     4 IDENT(6,6),DS_DDLAMBDA(6),F_DL_T1_1(6,6),F_DL_T1(6,6),PSTRAIN(6),
     5 F_DL_T2(6),DS_DDL_T1(6,6),DEQS_DDL_T2(6,6),NF(6),DEQS_DDL_T3(6),
     6 MAT(6,6),TV1(6),V2(6),C2(6),AIN(6,6),C(6,6),DDS_T3_1(6,6),T_INV(6,6),
     7 ESTRAIN(6),TSTRAIN(6),TSTRESS(6),CF(6,6),SE(6,6),COMPL(6,6),
     8 LSTRESS(6),LCF(6,6),LDSTRAIN(6),LSTRAIN(6),T(6,6),LSTRESSOLD(6),
     9 SE_FF1(1,1),SE_FF2(1,1),SE_IFF1(1,1),SE_IFF2(1,1),SE_IFF3(1,1)

      DOUBLE PRECISION    
     1 E1,E2,E3,G12,G13,G23,V12,V13,V23,V21,V31,V32,IFF1,IFF2,IFF3,
     2 XT,XC,YT,YC,VF12,EF1,EQPSTRAIN, F_DLAMBDA,
     3 S21,FFT,FFC,MFT,MFC,D_FT,D_FC,D_MT,D_MC,DMLC,FF10,FF20,IFF10,     
     5 IFF20,IFF30,ZERO,ONE,LC,THETA,E1O,XTF,A,FN0,C1,FS1

      DOUBLE PRECISION SIG_IFFV(6),SIG_IFF3V(6),NG_IFF(6),NG_IFF3(6),
     1  DEPSP_IFF(6), DEPSP_IFF3(6), DWP_IFF, DWP_IFF3, DWP

      DOUBLE PRECISION
     1 CF_INV(6,6), H_TILDE, NF_MAT(6,1), NF_MAT_TRANS(1,6),
     2 XI_TILDE, XI_TILDE_SQUARE(1,1), NF_NFTRANS(6,6), D_EP_INV(6,6) ,
     3 DPSTRAIN(6), ESTRESS_MAT(6,1), C_FF(6,6), C_IFF1(6,6), C_IFF3(6,6)

      DOUBLE PRECISION 
     1 I_PLUS_DLAMBDA_DE_P(6,6), I_PLUS_DLAMBDA_DE_P_INV(6,6), 
     2 I_PLUS_DLAMBDA_DE_P_INV_SIG(6,1), ESTRESS_NPLUS1(6,1),
     3 I_PLUS_DLAMBDA_DE_P_INV_SIG_TR(1,6), DSIG_BY_DDLAMBDA(6,1),
     4 DF_BY_DDLAMBDA_T1(1,1), DF_BY_DDLAMBDA_T1_SCALAR, P_SIG_NPLUS1(6,1),
     5 P_SIG_NPLUS1_TR(1,6), DP_BY_DDLAMBDA_T4(6,1), DP_BY_DDLAMBDA(1,1),
     6 DSIGY_BY_DDLAMBDA, PSIE_PLUS_F, PSIE_PLUS_M, PSIE_PLUS_S, INC_PCD_M, 
     7 INC_PCD_S, PHI_F, PHI_M, PHI_S, PCD_S, PCD_M, H_PCD_M, H_PCD_S,
     8 H_PLUS_F, H_PLUS_M, H_PLUS_S, YSTRESS_MIN, DEELAS(1:6), CF_NFMAT(6,1),
     9 X_B, X_A, NFMAT_CF_NFMAT, GD_F, GD_M, GD_S, G1, G2, PHI1, PHI2, PHI3,
     3 GFFIB,GFMAT, ESTRESS_NPLUS1_TR(1,6)

      DOUBLE PRECISION
     1 GC_S, GC_C, GC_T, MDOT, B1, B2, I_23_5, PHI_MAX, 
     2 INC_PCD_1, INC_PCD_2, INC_PCD_3, XSIG1(6), XSIG2(6), XSIG3(6),
     3 GC_FT, GC_FC, FF1, FF2, DFVOLD, DMVOLD, DSVOLD, ETA, DT,
     4 DFV, DMV, DSV, D_F, D_M, D_S, STRAN_TRANS(1,6), ESTRAN_TRANS(1,6)

      DOUBLE PRECISION SIG_FF(6,1),NF_MAT_FF(6,1),XI_TILDE_SQUARE_FF(1,1),
     1  SIG_IFF(6,1),NF_MAT_IFF(6,1),XI_TILDE_SQUARE_IFF(1,1),
     2  SIG_IFF3(6,1),NF_MAT_IFF3(6,1),XI_TILDE_SQUARE_IFF3(1,1),
     3  STRAN_MAT(6,1),ESTRAN_MAT(6,1),PD_FF(1,1),PD_IFF(1,1),PD_IFFS(1,1),
     4  IFF1_HIST,IFF2_HIST,IFF3_HIST,F1_HIST,FF2_HIST,FE_EFF_HIST,
     5  IFF1N,IFF2N,IFF3N,FF1N,FF2N,
     6  PHI_C_FF1,PHI_C_FF2,PHI_C_IFF1,PHI_C_IFF2,PHI_C_IFF3

      DOUBLE PRECISION AM(6,6),AT(1,6),JM(6,6),B(6,6),NG(6),F_INV(6,6),
     1 FM(6,6),NG_J_S,NG_J,ESTRESS_A(6),ESTRESS_B(6),ESTRESS_C(6),
     2 HM(6,6),NG_MAT(6,1),DDSDDE_1(6,6),DDSDDE_2(1,1)


      DOUBLE PRECISION ESTRESS_TR(6), SIGA(6), SIGB(6), SIGC(6)
      DOUBLE PRECISION YS_TR, fA, YS_B, fC
      DOUBLE PRECISION gA, gB, gC, gN, dg
      DOUBLE PRECISION EQPSN_A, EQPSN_B, EQPSN_C
      DOUBLE PRECISION DFDEP, MARG
      DOUBLE PRECISION DSIG_DG(6), NF_BAR(6)
      DOUBLE PRECISION qBTv(6)
      DOUBLE PRECISION DEQPS_DG, DF_DG
      DOUBLE PRECISION HNG(6), NBARH(6), DEN, LAMBDA
      INTEGER IT, MAXIT
      DOUBLE PRECISION TOL_F, TOL_G 

      INTEGER MAXITER, ITER,ITER1! 

      COMMON/KUSER/USRVAR(NEL,20,NGP),TIME_VAR,ITER_VAR          
      
!---------------------READ PROPERTIES-----------------------------------
      E1 = PROPS(1)           !YOUNG'S MODULUS IN DIRECTION 1 (L)
      E2 = PROPS(2)           !YOUNG'S MODULUS IN DIRECTION 2 (T) 
      E3=E2               
      G12 = PROPS(3)          !SHEAR MODULUS IN 12 PLANE
      G13=G12                 !SHEAR MODULUS IN 13 PLANE
      V12=PROPS(4)            !POISSON RATIO IN 12
      V23=PROPS(5)            !POISSON RATIO IN 23
      V13=V12                 !POISSON RATIO IN 13 
      EF1=PROPS(6)            !MODULUS OF FIBER PARALLEL TO FIBER
      VF12 = PROPS(7)         !POISSON RATIO OF FIBER
      XT = PROPS(8)           !TENSILE STRENGTH PARALLEL TO FIBER
      XC = PROPS(9)           !COMPRESSIVE STRENGTH PARALLEL TO FIBER
      YT = PROPS(10)          !TENSILE STRENGTH PERPENDICULAR TO FIBER
      YC = PROPS(11)          !COMPRESSIVE STRENGTH PERPENDICULAR TO FIBER
      S21 = PROPS(12)         !IN PLANE SHEAR STRENGTH
      MAT = PROPS(13)         !MATERIAL TYPE FOR INCLINATION PARAMETERS
      GC_T = PROPS(14)        !IFF1 GC
      GC_C = PROPS(15)        !IFF2 GC
      GC_S = PROPS(16)        !IFF3 GC
      GC_FT = PROPS(17)        !IFF2 GC
      GC_FC = PROPS(18)        !IFF3 GC
      G23 = E2/2/(1.+V23)     !SHEAR MODULUS IN 23 PLANE
      XTF=XT*EF1/E1           !EFFECTIVE TENSILE STRENGTH OF FIBER
C
      V21=(E2/E1)*V12
      V31=(E3/E1)*V13
      V32=(E3/E2)*V23 
!------------ INITIALIZATION OF MATRICES-------------------------------- 
      DO K1=1,6
          DO K2=1,6
              CF(K1,K2)=0.D0
              C_FF(K1,K2)=0.D0
              C_IFF1(K1,K2)=0.D0
              C_IFF3(K1,K2)=0.D0
              B(K1,K2)=0.D0
              AM(K1,K2)=0.D0 
              IDENT(K1,K2) = 0.0D0
              JM(K1,K2)=0.0D0    
          ENDDO
          IDENT(K1,K1) = 1.0D0
          AT(1,K1) = 0.0D0
      ENDDO 


      M = IDENT
      M_INV = IDENT    
      DO I=1,NTENS  ! INITIALIZATION OF JACOBIAN MATRIX
          DO J=1,NTENS
              DDSDDE(I,J) = ZERO
          ENDDO
      ENDDO
!------------ INITIALISATION OF STATE VARIABLES-------------------------
      IF (TIME(2) .EQ. ZERO) THEN
        STATEV = ZERO
        STATEV(35:39) = ONE
C         STATEV(26) = ONE
      END IF

      EELAS(1:6)    = STATEV(1:6)
      EPLAS(1:6)    = STATEV(7:12)
      ESTRESS(1:6)  = STATEV(13:18)
      EQPSN     = STATEV(19)   ! Equivalent Plastic Strain
      H_FFT= STATEV(20)
      H_FFC= STATEV(21) 
      H_IFF1 = STATEV(22) 
      H_IFF2 = STATEV(23) 
      H_IFF3 = STATEV(24) 
      D_FTVOLD = STATEV(25) 
      D_FCVOLD = STATEV(26) 
      D_MTVOLD = STATEV(27)
      D_MCVOLD = STATEV(28)
      D_SVOLD = STATEV(29)
      D_FT = STATEV(25) 
      D_FC = STATEV(26) 
      D_MT = STATEV(27)
      D_MC = STATEV(28)
      D_S = STATEV(29)
C       PHI_C_FF1 = STATEV(32)
C       PHI_C_FF2 = STATEV(33)
C       PHI_C_IFF1 = STATEV(34)
      FF10 = STATEV(35)
      FF20 = STATEV(36)
      IFF10 = STATEV(37)  
      IFF20 = STATEV(38) 
      IFF30 = STATEV(39) 
      FF1N = STATEV(40)
      FF2N= STATEV(41) 
      IFF1N = STATEV(42) 
      IFF2N = STATEV(43) 
      IFF3N = STATEV(44)  
      FE_EFFN = STATEV(45)
      DAM_COUNTER = STATEV(46)
      PD_IFF1 = STATEV(52)
      PD_IFF3 = STATEV(53)
      MAXIM = ZERO
      LC = CELENT  
!-------------------ELASTIC STIFFNESS MATRIX----------------------------
      DELTA=1.D0/(1.D0-V12*V21-V23*V32-V13*V31-2.D0*V21*V32*V13)
      CF(1,1) = E1*(1.D0-V23*V32)*DELTA
      CF(1,2) = E2*(V12+V32*V13)*DELTA
      CF(1,3) = E1*(V31+V21*V32)*DELTA
      CF(2,1) = CF(1,2)   
      CF(2,2) = E2*(1.D0-V13*V31)*DELTA     
      CF(2,3) = E2*(V32+V12*V31)*DELTA      
      CF(3,1) = CF(1,3)    
      CF(3,2) = CF(2,3)    
      CF(3,3) = E3*(1.D0-V12*V21)*DELTA   
      CF(4,4) = G12   
      CF(5,5) = G13   
      CF(6,6) = G23
      DDSDDE = CF
!---------------DEFINE THE COMPLIANCE MATRIX, COMPL---------------------
      COMPL = 0.D0 ! INITIALIZE AS ZEROS
      COMPL(1,1) = 1.D0/E1
      COMPL(1,2) = -V12/E1
      COMPL(1,3) = -V13/E1
      COMPL(2,1) = COMPL(1,2)
      COMPL(2,2) = 1.D0/E2     
      COMPL(2,3) = -V23/E2      
      COMPL(3,1) = COMPL(1,3)    
      COMPL(3,2) = COMPL(2,3)    
      COMPL(3,3) = 1.D0/E3   
      COMPL(4,4) = 1.D0/G12   
      COMPL(5,5) = 1.D0/G13
      COMPL(6,6) = 1.D0/G23
!-----------STIFFNESS MATRIX CORRESPONDING TO DIFFERENT MODES-----------
      C_FF(1,1) = CF(1,1)
      C_IFF3(4,4) = CF(4,4)
      C_IFF3(5,5) = CF(5,5)
      C_IFF1       = CF
      C_IFF1(1,1)  = 0.D0
      C_IFF1(4,4)  = 0.D0
      C_IFF1(5,5)  = 0.D0 
!------------------------PLASTICITY PARAMETERS--------------------------   
      NU_P = 0.32D0
      A = 1.50D0
      ALPHA = 0.35D0
      BETA = 355.0274D0
      ETA = 0.0002D0
      XK = 1.0D-6   

      JM(1,1) = 0.5D0
      JM(2,2) = 0.5D0
      JM(3,3) = 0.5D0
      JM(4,4) = 0.25D0
      JM(5,5) = 0.25D0
      JM(6,6) = 0.25D0

      B(2,2) = 1.0D0
      B(2,3) = -NU_P
      B(3,2) = -NU_P
      B(3,3) = 1.0D0
      B(4,4) = 2.0D0*(1.0D0 + NU_P)
      B(5,5) = 2.0D0*(1.0D0 + NU_P)
      B(6,6) = 2.0D0*(1.0D0 + NU_P)

! INITIALISE PARAMETERS      
      NG(:) = 0.0d0
      FF1 = 0.0D0
      FF2 = 0.0D0
      IFF1 = 0.0D0  
      IFF2 = 0.0D0 
      IFF3 = 0.0D0

 !     Phase-field int pt value  
C         IF (ITER.EQ.ZERO) THEN           !  from 2nd (phase-field) layer at the end of last increment 
C          PHI1=USRVAR(NOEL,11,NPT)
C          PHI2=USRVAR(NOEL,12,NPT)
C          PHI3=USRVAR(NOEL,13,NPT)
C          PHI4=USRVAR(NOEL,14,NPT)
C          PHI5=USRVAR(NOEL,15,NPT)
C         ELSE                             !  from 1st (displacement) layer value
         PHI1=USRVAR(NOEL,1,NPT)
         PHI2=USRVAR(NOEL,2,NPT)
         PHI3=USRVAR(NOEL,3,NPT)
         PHI4=USRVAR(NOEL,4,NPT)
         PHI5=USRVAR(NOEL,5,NPT) 
C         ENDIF 


     

!-------------- CALCULATE PREDICTOR STRESS AND ELASTIC STRAIN ----------
      DO K1=1,NTENS
        EELAS(K1) = EELAS(K1) + DSTRAN(K1)       ! TRIAL ELASTIC STRAIN    
        DO K2=1,NTENS
          ESTRESS(K2)= ESTRESS(K2)+DDSDDE(K2,K1)*DSTRAN(K1) ! TRIAL VALUE
        ENDDO                
      ENDDO 
      
      DO I = 1,6
        ESTRESS_MAT(I,1) = ESTRESS(I)
      END DO  


C=======================================================================
C  V O G L E R   P L A S T I C I T Y   (11,22,33,12,13,23)
C  Non-associated: A != B
C  Mixed NR + bisection safeguarded; consistent DDSDDE (Van der Meer)
C=======================================================================
      MAXIT = 200
      TOL_F = 1.D-4
      TOL_G = 1.D-11

C-----------------------------------------------------------------------
C Trial stress saved
C-----------------------------------------------------------------------
      ESTRESS_TR = ESTRESS

C-----------------------------------------------------------------------
C Trial yield check
C-----------------------------------------------------------------------
      CALL CALC_YS_FULL(EQPSN,ESTRESS_TR,YS_TR,AM,AT,DFDEP)
C        YS_TR = -1.0d0
      IF ((YS_TR .LE. 0.D0) .OR. (DAM_COUNTER .GE. 1.D0)) THEN
C       IF (YS_TR .LE. 0.D0) THEN
        DGAMMA = 0.D0
        ESTRESS = ESTRESS_TR
        DDSDDE = CF

      ELSE
C-----------------------------------------------------------------------
C 1) Bracketing for DGAMMA: find gA,gB such that f(gA)>0, f(gB)<0
C    start gA = 0, gB small then expand.
C-----------------------------------------------------------------------
      gA = 0.D0
      SIGA = ESTRESS_TR
      EQPSN_A  = EQPSN
      fA   = YS_TR

      gB = 1.D-10
      DO IT = 1, 60
        F_INV = IDENT + gB*MATMUL(CF,B)
        CALL INVERSE(F_INV,FM,6)
        SIGB = MATMUL(FM,ESTRESS_TR)

        NG = MATMUL(B,SIGB)
        CALL TV_MAT_V(NG,JM,NG,NG_J_S)
        NG_J = SQRT(MAX(NG_J_S,1.D-30))
        EQPSN_B = EQPSN + gB*NG_J

        CALL CALC_YS_FULL(EQPSN_B,SIGB,YS_B,AM,AT,DFDEP)

        IF (YS_B .LT. 0.D0) EXIT
        gB = gB*2.D0
      END DO

      IF (YS_B .GE. 0.D0) THEN
C       Could not bracket -> fall back to very small step elastic-like
        DGAMMA = 0.D0
        ESTRESS = ESTRESS_TR
        DDSDDE = CF
        GO TO 9999
      END IF

C-----------------------------------------------------------------------
C 2) Safeguarded Newton iterations
C-----------------------------------------------------------------------
      gC = 0.5D0*(gA+gB)

      DO IT = 1, MAXIT

C       Evaluate at gC
        F_INV = IDENT + gC*MATMUL(CF,B)
        CALL INVERSE(F_INV,FM,6)
        SIGC = MATMUL(FM,ESTRESS_TR)

        NG = MATMUL(B,SIGC)
        CALL TV_MAT_V(NG,JM,NG,NG_J_S)
        NG_J = SQRT(MAX(NG_J_S,1.D-30))
        EQPSN_C  = EQPSN + gC*NG_J

        CALL CALC_YS_FULL(EQPSN_C,SIGC,fC,AM,AT,DFDEP)

C       Converged?
        IF (ABS(fC) .LT. TOL_F) EXIT

C------------------------------------------------------------
C       Compute DF/DGAMMA for Newton step (consistent chain rule)
C       dSigma/dG = - FM * CF * B * sigma
C------------------------------------------------------------
        DSIG_DG = MATMUL( MATMUL( -FM, MATMUL(CF,B) ), SIGC )

C       n_f = A*sigma + a
        DO I=1,6
          NF(I) = AT(1,I)
          DO J=1,6
            NF(I) = NF(I) + AM(I,J)*SIGC(J)
          END DO
        END DO

C       TV1 = (JM*ng)/||ng||  (6-vector)
        TV1 = MATMUL(JM,NG)
        DO I=1,6
          TV1(I) = TV1(I)/NG_J
        END DO

C       qBTv = B^T * TV1
        DO I=1,6
          qBTv(I) = 0.D0
          DO J=1,6
            qBTv(I) = qBTv(I) + B(J,I)*TV1(J)
          END DO
        END DO

C       dEqps/dG = ||ng|| + (dEqps/dSigma)·(dSigma/dG)
C       with dEqps/dSigma = -gC * (B^T*(JM*ng/||ng||))
        DEQPS_DG = NG_J
        DO I=1,6
          DEQPS_DG = DEQPS_DG - gC*qBTv(I)*DSIG_DG(I)
        END DO

C       DF/DG = n_f·dSigma/dG + df/dep * dEqps/dG
        DF_DG = 0.D0
        DO I=1,6
          DF_DG = DF_DG + NF(I)*DSIG_DG(I)
        END DO
        DF_DG = DF_DG + DFDEP*DEQPS_DG

C       Newton proposal
        IF (ABS(DF_DG) .LT. 1.D-18) THEN
          gN = 0.5D0*(gA+gB)
        ELSE
          gN = gC - fC/DF_DG
        END IF

C       Safeguard: if outside bracket -> bisection
        IF ( (gN .LE. gA) .OR. (gN .GE. gB) ) THEN
          gN = 0.5D0*(gA+gB)
        END IF

C       Update bracket using sign
        IF (fA*fC .LT. 0.D0) THEN
          gB = gC
          YS_B = fC
        ELSE
          gA = gC
          fA = fC
        END IF

C       Update iterate
        dg = ABS(gN - gC)
        gC = gN

        IF (dg .LT. TOL_G) EXIT
      END DO

C-----------------------------------------------------------------------
C Accept solution
C-----------------------------------------------------------------------
      DGAMMA  = gC
      ESTRESS = SIGC
      EQPSN   = EQPSN_C

C-----------------------------------------------------------------------
C 3) Consistent tangent (Van der Meer style)
C    H = FM*CF
C    nbar = nf - df/dep * DGAMMA * (B^T*(JM*ng/||ng||))
C    lambda = - df/dep * ||ng||_J
C    DDSDDE = H - (H*ng ⊗ (nbar^T*H)) / (lambda + nbar^T*H*ng)
C-----------------------------------------------------------------------
      HM = MATMUL(FM,CF)

C     nbar
      DO I=1,6
        NF_BAR(I) = NF(I) - DFDEP*DGAMMA*qBTv(I)
      END DO
      LAMBDA = -DFDEP*NG_J

C     H*ng
      DO I=1,6
        HNG(I) = 0.D0
        DO J=1,6
          HNG(I) = HNG(I) + HM(I,J)*NG(J)
        END DO
      END DO

C     nbar^T * H
      DO J=1,6
        NBARH(J) = 0.D0
        DO I=1,6
          NBARH(J) = NBARH(J) + NF_BAR(I)*HM(I,J)
        END DO
      END DO

C     denom
      DEN = LAMBDA
      DO I=1,6
        DEN = DEN + NF_BAR(I)*HNG(I)
      END DO
      IF (ABS(DEN) .LT. 1.D-18) DEN = SIGN(1.D-18,DEN)

C     DDSDDE
      DO I=1,6
        DO J=1,6
          DDSDDE(I,J) = HM(I,J) - (HNG(I)*NBARH(J))/DEN
        END DO
      END DO

9999  CONTINUE
      ENDIF

      EPLAS(:) = EPLAS(:) + DGAMMA*NG(:)
      EELAS(:) = STRAN(:) + DSTRAN(:) - EPLAS(:)
C       STRESS = ESTRESS


      CALL CUNTZE(
     1     ESTRESS,S21,XT,XC,YT,YC,FF1,FF2,IFF1,IFF2,IFF3,FE_EFF,MDOT)
      
      IFF1_HIST=MAX(IFF1, IFF1N)
      IFF2_HIST=MAX(IFF2, IFF2N)
      IFF3_HIST=MAX(IFF3, IFF3N)
      FF1_HIST =MAX(FF1,  FF1N)
      FF2_HIST =MAX(FF2,  FF2N)
      FE_EFF_HIST = MAX(FE_EFF, FE_EFFN)

      IF (FE_EFF.GT.FE_EFFN .AND. FE_EFF.GT.1.D0) THEN  ! ADDED DAMCOUNTER CONDITION
        X_IND = 1.D0 ! BECOMES 1 IF ACTIVATED.
        DAM_COUNTER = DAM_COUNTER + 1.D0
        FN0 = MAX(FF1,FF2,IFF1,IFF2,IFF3)
        IF(FN0 .EQ. FF1)   FF10 = MIN(FN0,FF10)
        IF(FN0 .EQ. FF2)   FF20 = MIN(FN0,FF20)
        IF(FN0 .EQ. IFF1)  IFF10 = MIN(FN0,IFF10)
        IF(FN0 .EQ. IFF2)  IFF20 = MIN(FN0,IFF20)
        IF(FN0 .EQ. IFF3)  IFF30 = MIN(FN0,IFF30)
      END IF
      ! DAMAGE DRIVING 
      ESTRAN_TRANS(1,:) = EELAS(:)  
      ESTRAN_MAT(:,1) = EELAS(:) 
      SE_FF1 = 0.5D0 * MATMUL(MATMUL(ESTRAN_TRANS,C_FF),TRANSPOSE(
     1       ESTRAN_TRANS))
      SE_FF2 = 0.5D0 * MATMUL(MATMUL(ESTRAN_TRANS,C_FF),TRANSPOSE(
     1       ESTRAN_TRANS))
      SE_IFF1 = 0.5D0 * MATMUL(MATMUL(ESTRAN_TRANS,C_IFF1),TRANSPOSE(
     1       ESTRAN_TRANS))
      SE_IFF2 = 0.5D0 * MATMUL(MATMUL(ESTRAN_TRANS,C_IFF1),TRANSPOSE(
     1       ESTRAN_TRANS))
      SE_IFF3 = 0.5D0 * MATMUL(MATMUL(ESTRAN_TRANS,C_IFF3),TRANSPOSE(
     1       ESTRAN_TRANS)) 

! ---- compute mode stresses from mode stiffness matrices
      SIG_IFFV  = MATMUL(C_IFF1, EELAS)   ! 6-vector
      SIG_IFF3V = MATMUL(C_IFF3, EELAS)   ! 6-vector

! ---- build projected flow directions (simple projectors)
      NG_IFF  = 0.D0
      NG_IFF3 = 0.D0

! IFF1: (22,33,23) -> indices 2,3,6
      NG_IFF(2) = NG(2)     
      NG_IFF(3) = NG(3)
      NG_IFF(6) = NG(6)

! IFF3: (12,13) -> indices 4,5
      NG_IFF3(4) = NG(4)
      NG_IFF3(5) = NG(5)

! ---- plastic strain increments in each mode
      DEPSP_IFF(:)  = DGAMMA * NG_IFF(:)
      DEPSP_IFF3(:) = DGAMMA * NG_IFF3(:)

! ---- mode plastic work contributions
      DWP_IFF = SIG_IFFV(1)*DEPSP_IFF(1) + SIG_IFFV(2)*DEPSP_IFF(2)  +
     1     SIG_IFFV(3)*DEPSP_IFF(3)  + SIG_IFFV(4)*DEPSP_IFF(4)  +
     2     SIG_IFFV(5)*DEPSP_IFF(5)  + SIG_IFFV(6)*DEPSP_IFF(6)

      DWP_IFF3 = SIG_IFF3V(1)*DEPSP_IFF3(1) + SIG_IFF3V(2)*DEPSP_IFF3(2) +
     1     SIG_IFF3V(3)*DEPSP_IFF3(3) + SIG_IFF3V(4)*DEPSP_IFF3(4) +
     2     SIG_IFF3V(5)*DEPSP_IFF3(5) + SIG_IFF3V(6)*DEPSP_IFF3(6)


      DWP = DWP_IFF + DWP_IFF3
      DWP = MAX(0.D0, DWP)   ! ENFORCE DISSIPATION 
      DWP_IFF  = MAX(0.D0, DWP_IFF)
      DWP_IFF3 = MAX(0.D0, DWP_IFF3)

      PD_IFF1 = PD_IFF1 + DWP_IFF
      PD_IFF3 = PD_IFF3 + DWP_IFF3


    
      H_FFT  = MAX(H_FFT,SE_FF1(1,1))
      H_FFC  = MAX(H_FFC,SE_FF2(1,1))
      H_IFF1 = MAX(H_IFF1,(SE_IFF1(1,1) + PD_IFF1))
      H_IFF2 = MAX(H_IFF2,(SE_IFF2(1,1) + PD_IFF1))
      H_IFF3 = MAX(H_IFF3,(SE_IFF3(1,1) + PD_IFF3))
      
C       H_FFT  = MAX(H_FFT,SE_FF1(1,1))
C       H_FFC  = MAX(H_FFC,SE_FF2(1,1))
C       H_IFF1 = MAX(H_IFF1,SE_IFF1(1,1))
C       H_IFF2 = MAX(H_IFF2,SE_IFF2(1,1))
C       H_IFF3 = MAX(H_IFF3,SE_IFF3(1,1))

C       PRINT*,"DAMAGE DRIVING:"
C       PRINT*,H_FFT,H_FFC,H_IFF1,H_IFF2,H_IFF3

! CALCULATE THRESHOLD ENERGIES IF DAMAGE HAS JUST BEGUN
      IF ((FF1 .LT. FF10) .AND. (DAM_COUNTER .LE. ONE))  THEN
        PHI_C_FF1 = H_FFT 
      ELSE
        PHI_C_FF1 = STATEV(30)         
      ENDIF

      IF ((FF2 .LT. FF20) .AND. (DAM_COUNTER .LE. ONE))  THEN
        PHI_C_FF2 = H_FFC 
      ELSE
        PHI_C_FF2 = STATEV(31)         
      ENDIF

      IF ((IFF1 .LT. IFF10) .AND. (DAM_COUNTER .LE. ONE))  THEN
        PHI_C_IFF1 = H_IFF1 
      ELSE
        PHI_C_IFF1 = STATEV(32)         
      ENDIF

      IF ((IFF2 .LT. IFF20) .AND. (DAM_COUNTER .LE. ONE))  THEN
        PHI_C_IFF2 = H_IFF2
      ELSE
        PHI_C_IFF2 = STATEV(33)         
      ENDIF

      IF ((IFF3 .LT. IFF30) .AND. (DAM_COUNTER .LE. ONE))  THEN
        PHI_C_IFF3 = H_IFF3 
      ELSE
        PHI_C_IFF3 = STATEV(34)         
      ENDIF

      PHI_S_IFF1 = MAX( (GC_T /LC - PHI_C_IFF1), 0.001D0 )
      PHI_S_IFF2 = MAX( (GC_C /LC - PHI_C_IFF2), 0.001D0 )
      PHI_S_IFF3 = MAX( (GC_S /LC - PHI_C_IFF3), 0.001D0 )
      PHI_S_FF1  = MAX( (GC_FT/LC - PHI_C_FF1 ), 0.001D0 )
      PHI_S_FF2  = MAX( (GC_FC/LC - PHI_C_FF2 ), 0.001D0 )

      IF(DAM_COUNTER .GE. ONE) THEN
        IF(FF1 .GT. FF10) THEN
          CALL DAMAGE(H_FFT,PHI_C_FF1,PHI_S_FF1,D_FT)
          D_FC = ZERO
C           PRINT*, "H_FFT",D_FT,D_FC
        ELSE IF(FF2 .GT. FF20) THEN
          CALL DAMAGE(H_FFC,PHI_C_FF2,PHI_S_FF2,D_FC)
          D_FT = ZERO
C          PRINT*,"H_FFC",D_FT,D_FC
        END IF

        IF(IFF1 .GT. IFF10) THEN
          CALL DAMAGE(H_IFF1,PHI_C_IFF1,PHI_S_IFF1,D_MT)
          D_MC = ZERO
        ELSE IF(IFF2 .GT. IFF20) THEN
          CALL DAMAGE(H_IFF2,PHI_C_IFF2,PHI_S_IFF2,D_MC)
          D_MT = ZERO
        END IF

        IF(IFF3 .GT. IFF30) THEN
          CALL DAMAGE(H_IFF3,PHI_C_IFF3,PHI_S_IFF3,D_S)
        END IF

C         IF (NOEL .EQ. 1.0D0) PRINT*,D_MC,H_IFF2,PHI_C_IFF2,PHI_S_IFF2,GC_C,2*LC
        D_FTV = DTIME/(ETA + DTIME)*D_FT + ETA/(ETA + DTIME)*D_FTVOLD
        D_FCV = DTIME/(ETA + DTIME)*D_FC + ETA/(ETA + DTIME)*D_FCVOLD
        D_MTV = DTIME/(ETA + DTIME)*D_MT + ETA/(ETA + DTIME)*D_MTVOLD
        D_MCV = DTIME/(ETA + DTIME)*D_MC + ETA/(ETA + DTIME)*D_MCVOLD
        D_SV  = DTIME/(ETA + DTIME)*D_S  + ETA/(ETA + DTIME)*D_SVOLD

!     DEGRADATRION FUNCTIONS DEFINED BELOW
C         GD1 = (1.D0 - MIN(D_FTV,0.9999D0))**2 
C         GD2 = (1.D0 - MIN(D_FCV,0.9999D0))**2 
C         GD3 = (1.D0 - MIN(D_MTV,0.9999D0))**2 
C         GD4 = (1.D0 - MIN(D_MCV,0.9999D0))**2 
C         GD5 = (1.D0 - MIN(D_SV, 0.9999D0))**2  
!     DEGRADATRION FUNCTIONS DEFINED BELOW
        GD1 = (1.D0 - MIN(PHI1,1.0D0))**2 + XK
        GD2 = (1.D0 - MIN(PHI2,1.0D0))**2 + XK
        GD3 = (1.D0 - MIN(PHI3,1.0D0))**2 + XK
        GD4 = (1.D0 - MIN(PHI4,1.0D0))**2 + XK
        GD5 = (1.D0 - MIN(PHI5,1.0D0))**2 + XK


C         STRESS(1) = GD1*GD2*ESTRESS(1)
C         STRESS(2) = GD3*GD5*GD4*ESTRESS(2)
C         STRESS(3) = GD3*GD5*GD4*ESTRESS(3)
C         STRESS(4) = GD5*ESTRESS(4)
C         STRESS(5) = GD5*ESTRESS(5)
C         STRESS(6) = GD3*GD5*GD4*ESTRESS(6)
!     DEGRADATRION FUNCTIONS DEFINED BELOW
C       GD1 = (1.D0 - MIN(PHI1,0.99D0))**2 + XK
C       GD2 = (1.D0 - MIN(PHI2,0.99D0))**2 + XK
C       GD3 = (1.D0 - MIN(PHI3,0.99D0))**2 + XK
C       GD4 = (1.D0 - MIN(PHI4,0.99D0))**2 + XK
C       GD5 = (1.D0 - MIN(PHI5,0.99D0))**2 + XK

        DDSDDE(1,1) = GD1*GD2*DDSDDE(1,1)
        DDSDDE(2,2) = GD3*GD4*DDSDDE(2,2)
        DDSDDE(3,3) = GD3*GD4*DDSDDE(3,3)
        DDSDDE(4,4) = GD5*DDSDDE(4,4)
        DDSDDE(5,5) = GD5*DDSDDE(5,5)
        DDSDDE(6,6) = GD3*GD4*DDSDDE(6,6)
        DDSDDE(1,2) = GD3*GD4*DDSDDE(1,2)
        DDSDDE(2,1) = GD3*GD4*DDSDDE(2,1)
        DDSDDE(1,3) = GD3*GD4*DDSDDE(1,3)
        DDSDDE(3,1) = GD3*GD4*DDSDDE(3,1)
        DDSDDE(3,2) = GD3*GD4*DDSDDE(3,2)
        DDSDDE(2,3) = GD3*GD4*DDSDDE(2,3)

C         ESTRESS_TR = MATMUL(CF(:,:),EELAS(:)) 

        CF(1,1) = GD1*GD2*CF(1,1)
        CF(2,2) = GD3*GD4*CF(2,2)
        CF(3,3) = GD3*GD4*CF(3,3)
        CF(4,4) = GD5*CF(4,4)
        CF(5,5) = GD5*CF(5,5)
        CF(6,6) = GD3*GD4*CF(6,6)
        CF(1,2) = GD3*GD4*CF(1,2)
        CF(2,1) = GD3*GD4*CF(2,1)
        CF(1,3) = GD3*GD4*CF(1,3)
        CF(3,1) = GD3*GD4*CF(3,1)
        CF(3,2) = GD3*GD4*CF(3,2)
        CF(2,3) = GD3*GD4*CF(2,3)

        STRESS(:) = MATMUL(CF(:,:),EELAS(:))
C         PRINT*,GD1,GD2,GD3,GD4,GD5,STRESS(2),ESTRESS_TR(2)

      ELSE
        STRESS = ESTRESS
      END IF

      SSE = ZERO
      DO K1 = 1, NTENS
         SSE= SSE+ ESTRESS(K1)*EELAS(K1)/TWO  !Update specific elastic strain energy
      ENDDO
!------------- STORE STRAINS IN STATE VARIABLE ARRAY ----------------------
      DO K1=1,NTENS ! FROM 1 TO 18
        STATEV(K1)= EELAS(K1)
        STATEV(K1+NTENS)= EPLAS(K1)
        STATEV(K1+2*NTENS)= ESTRESS(K1)
      ENDDO 
      
      STATEV(19) = EQPSN
      STATEV(20) = H_FFT
      STATEV(21) = H_FFC
      STATEV(22) = H_IFF1
      STATEV(23) = H_IFF2
      STATEV(24) = H_IFF3
      STATEV(25) = D_FTV
      STATEV(26) = D_FCV
      STATEV(27) = D_MTV
      STATEV(28) = D_MCV      
      STATEV(29) = D_SV
      STATEV(30) = PHI_C_FF1
      STATEV(31) = PHI_C_FF2
      STATEV(32) = PHI_C_IFF1
      STATEV(33) = PHI_C_IFF2
      STATEV(34) = PHI_C_IFF3
      STATEV(35) = FF10
      STATEV(36) = FF20
      STATEV(37) = IFF10
      STATEV(38) = IFF20
      STATEV(39) = IFF30
      STATEV(40) = FF1
      STATEV(41) = FF2
      STATEV(42) = IFF1
      STATEV(43) = IFF2
      STATEV(44) = IFF3
      STATEV(45) = FE_EFF
      STATEV(46) = DAM_COUNTER
      STATEV(47) = PHI1
      STATEV(48) = PHI2
      STATEV(49) = PHI3
      STATEV(50) = PHI4
      STATEV(51) = PHI5 
      STATEV(52) = PD_IFF1  
      STATEV(53) = PD_IFF3    
      STATEV(54) = LC         
      USRVAR(NOEL, 6, NPT) = D_FTV
      USRVAR(NOEL, 7, NPT) = D_FCV
      USRVAR(NOEL, 8, NPT) = D_MTV
      USRVAR(NOEL, 9, NPT) = D_MCV
      USRVAR(NOEL, 10, NPT)= D_SV
      USRVAR(NOEL, 11, NPT)= LC
      RETURN
      END



!     ==================================================================
!     User subroutine UEL for phase-field element
!     ==================================================================
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     &     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     &     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     &     PERIOD)

      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*3 JTYPE

C === PARAMETERS ===
      PARAMETER (NGP=8,ZERO=0.D0, ONE=1.D0, NINPT=8, NEL=75601)
      COMMON/KUSER/USRVAR(NEL,20,NGP),TIME_VAR,ITER_VAR

C === INPUT ARRAYS ===
      DIMENSION RHS(MLVARX,1), AMATRX(NDOFEL,NDOFEL), SVARS(NSVARS),
     1          ENERGY(8), PROPS(NPROPS), COORDS(3,8),
     2          U(NDOFEL), DU(MLVARX,1), V(NDOFEL), A(NDOFEL),
     3          TIME(2), PARAMS(3), JDLTYP(MDLOAD,*), ADLMAG(MDLOAD,*),
     4          DDLMAG(MDLOAD,*), PREDEF(2,NPREDF,NNODE),
     5          LFLAGS(*), JPROPS(*)

C === LOCAL VARIABLES ===
      DOUBLE PRECISION XLC, DETJ, WEIGHT, XI(3), AINTW(8)
      DOUBLE PRECISION N(8), DNDXI(8,3), DNDX(8,3), JAC(3,3), JINV(3,3)
      DOUBLE PRECISION BPHI(3,8), N_MAT(8,1)

      DOUBLE PRECISION D_FFT_NODAL(8), D_FFC_NODAL(8), D_IFF1_NODAL(8)
      DOUBLE PRECISION D_IFF2_NODAL(8), D_IFF3_NODAL(8)
      DOUBLE PRECISION D_FFT_GP, D_FFC_GP, D_IFF1_GP, D_IFF2_GP
      DOUBLE PRECISION D_IFF3_GP, GP_LC

      DOUBLE PRECISION PHI1_GP, PHI2_GP, PHI3_GP, PHI4_GP, PHI5_GP

      DOUBLE PRECISION K_FFT(8,8), R_FFT(8), M(3), MDYADM(3,3)
      DOUBLE PRECISION K_FFC(8,8), R_FFC(8), STRUCT_TENSOR_FIBRE(3,3)
      DOUBLE PRECISION K_IFF1(8,8), R_IFF1(8), STRUCT_TENSOR(3,3)
      DOUBLE PRECISION K_IFF2(8,8), R_IFF2(8), IDEN(3,3)
      DOUBLE PRECISION K_IFF3(8,8), R_IFF3(8), SDV(5)

      INTEGER I, J, GP, K
      DOUBLE PRECISION THETA, BETA, PI, PLY_ANGLE, COSTHETA, SINTHETA

C === INITIALIZATION
C       XLC   = PROPS(1)
      THETA = PROPS(2)

      AMATRX(:,:) = ZERO
      RHS(:,:)    = ZERO
      ENERGY(:)   = ZERO

      BETA = 25.D0
      PI   = 4.0D0*ATAN(1.0D0)
      PLY_ANGLE = THETA*PI/180.D0

      COSTHETA = COS(PLY_ANGLE)
      SINTHETA = SIN(PLY_ANGLE)
      M = (/ COSTHETA, SINTHETA, 0.D0 /)

      CALL ONEM(IDEN)
      CALL DYADIC(M, M, 3, MDYADM)

      STRUCT_TENSOR       = IDEN + BETA*MDYADM
      STRUCT_TENSOR_FIBRE = IDEN + BETA*(IDEN - MDYADM)

      K_FFT(:,:)  = ZERO
      K_FFC(:,:)  = ZERO
      K_IFF1(:,:) = ZERO
      K_IFF2(:,:) = ZERO
      K_IFF3(:,:) = ZERO

      R_FFT(:)  = ZERO
      R_FFC(:)  = ZERO
      R_IFF1(:) = ZERO
      R_IFF2(:) = ZERO
      R_IFF3(:) = ZERO

      DO I=1,NGP
         AINTW(I) = ONE
      END DO

C === EXTRACT NODAL DOFs FOR THE FIVE PHASE FIELDS
      DO I = 1, 8
         D_FFT_NODAL(I)  = U(5*(I-1) + 1)
         D_FFC_NODAL(I)  = U(5*(I-1) + 2)
         D_IFF1_NODAL(I) = U(5*(I-1) + 3)
         D_IFF2_NODAL(I) = U(5*(I-1) + 4)
         D_IFF3_NODAL(I) = U(5*(I-1) + 5)
      END DO

C === LOOP OVER GAUSS POINTS
      DO GP = 1, NGP

C ---- Load current element state variables at this GP
         DO I = 1, 5
            SDV(I) = SVARS(5*(GP-1) + I)
         END DO

C ---- Shape functions / Jacobian
         CALL GAUSSPOINT(GP, XI)
         CALL SHAPEFUN(N, DNDXI, XI)
         CALL JACOBIAN(COORDS, DNDXI, JAC, DETJ)

         IF (DETJ .LT. ZERO) THEN
            WRITE(7,*) 'Negative Jacobian', DETJ
            CALL XIT
         END IF

         CALL INVERSE3X3(JAC, JINV)

C ---- Compute dNdx, BPHI, N_MAT
         DO K = 1, 8
            DO I = 1, 3
               DNDX(K,I) = 0.D0
               DO J = 1, 3
                  DNDX(K,I) = DNDX(K,I) + DNDXI(K,J)*JINV(J,I)
               END DO
            END DO
         END DO

         DO K = 1, 8
            BPHI(1,K) = DNDX(K,1)
            BPHI(2,K) = DNDX(K,2)
            BPHI(3,K) = DNDX(K,3)
            N_MAT(K,1)= N(K)
         END DO

         WEIGHT = DETJ * AINTW(GP)

C ---- Receive local damage/driving from UMAT: slots 6:10
         D_FFT_GP  = USRVAR(JELEM-NEL, 6, GP)
         D_FFC_GP  = USRVAR(JELEM-NEL, 7, GP)
         D_IFF1_GP = USRVAR(JELEM-NEL, 8, GP)
         D_IFF2_GP = USRVAR(JELEM-NEL, 9, GP)
         D_IFF3_GP = USRVAR(JELEM-NEL,10, GP)
         GP_LC = USRVAR(JELEM-NEL,11, GP)
         XLC = 2.0D0 * GP_LC

C ---- Interpolate current nonlocal phase fields at this GP from nodal U
         PHI1_GP = 0.D0
         PHI2_GP = 0.D0
         PHI3_GP = 0.D0
         PHI4_GP = 0.D0
         PHI5_GP = 0.D0

         DO I = 1, 8
            PHI1_GP = PHI1_GP + N(I)*D_FFT_NODAL(I)
            PHI2_GP = PHI2_GP + N(I)*D_FFC_NODAL(I)
            PHI3_GP = PHI3_GP + N(I)*D_IFF1_NODAL(I)
            PHI4_GP = PHI4_GP + N(I)*D_IFF2_NODAL(I)
            PHI5_GP = PHI5_GP + N(I)*D_IFF3_NODAL(I)
         END DO

C ---- Store GP phase fields in SVARS
         SDV(1) = PHI1_GP
         SDV(2) = PHI2_GP
         SDV(3) = PHI3_GP
         SDV(4) = PHI4_GP
         SDV(5) = PHI5_GP

C ---- Assemble residual and tangent for each field
         CALL ASSEMBLE_FIELD(D_FFT_NODAL, D_FFT_GP, BPHI, N_MAT, XLC,
     1        WEIGHT, K_FFT, R_FFT, STRUCT_TENSOR_FIBRE)

         CALL ASSEMBLE_FIELD(D_FFC_NODAL, D_FFC_GP, BPHI, N_MAT, XLC,
     1        WEIGHT, K_FFC, R_FFC, STRUCT_TENSOR_FIBRE)

         CALL ASSEMBLE_FIELD(D_IFF1_NODAL, D_IFF1_GP, BPHI, N_MAT, XLC,
     1        WEIGHT, K_IFF1, R_IFF1, STRUCT_TENSOR)

         CALL ASSEMBLE_FIELD(D_IFF2_NODAL, D_IFF2_GP, BPHI, N_MAT, XLC,
     1        WEIGHT, K_IFF2, R_IFF2, STRUCT_TENSOR)

         CALL ASSEMBLE_FIELD(D_IFF3_NODAL, D_IFF3_GP, BPHI, N_MAT, XLC,
     1        WEIGHT, K_IFF3, R_IFF3, STRUCT_TENSOR)

C ---- Transfer nonlocal phase fields back to UMAT: slots 1:5
         USRVAR(JELEM-NEL,1,GP) = MIN(0.999D0, MAX(0.D0, PHI1_GP))
         USRVAR(JELEM-NEL,2,GP) = MIN(0.999D0, MAX(0.D0, PHI2_GP))
         USRVAR(JELEM-NEL,3,GP) = MIN(0.999D0, MAX(0.D0, PHI3_GP))
         USRVAR(JELEM-NEL,4,GP) = MIN(0.999D0, MAX(0.D0, PHI4_GP))
         USRVAR(JELEM-NEL,5,GP) = MIN(0.999D0, MAX(0.D0, PHI5_GP))

C ---- Store updated GP SVARS
         DO I = 1, 5
            SVARS(5*(GP-1) + I) = SDV(I)
         END DO

      END DO

C === FINAL ASSEMBLY INTO ELEMENT RHS / STIFFNESS
      DO I = 1, 8
         RHS(5*I-4,1) = RHS(5*I-4,1) + R_FFT(I)
         RHS(5*I-3,1) = RHS(5*I-3,1) + R_FFC(I)
         RHS(5*I-2,1) = RHS(5*I-2,1) + R_IFF1(I)
         RHS(5*I-1,1) = RHS(5*I-1,1) + R_IFF2(I)
         RHS(5*I  ,1) = RHS(5*I  ,1) + R_IFF3(I)

         DO J = 1, 8
            AMATRX(5*I-4,5*J-4) = AMATRX(5*I-4,5*J-4) + K_FFT(I,J)
            AMATRX(5*I-3,5*J-3) = AMATRX(5*I-3,5*J-3) + K_FFC(I,J)
            AMATRX(5*I-2,5*J-2) = AMATRX(5*I-2,5*J-2) + K_IFF1(I,J)
            AMATRX(5*I-1,5*J-1) = AMATRX(5*I-1,5*J-1) + K_IFF2(I,J)
            AMATRX(5*I  ,5*J  ) = AMATRX(5*I  ,5*J  ) + K_IFF3(I,J)
         END DO
      END DO

      PNEWDT = ONE
      RETURN
      END

!     ==================================================================      
!     Global Assembly
!     ==================================================================
      SUBROUTINE ASSEMBLE_FIELD(dnodal,d_gp,BPHI,N_MAT,XLC,WEIGHT,K,R,A)
      DOUBLE PRECISION dnodal(8), d_gp, BPHI(3,8), N_MAT(8,1)
      DOUBLE PRECISION XLC, WEIGHT, K(8,8), R(8)
      DOUBLE PRECISION Kloc(8,8), A(3,3)

      ZETA = 1.0D-7 
      Kloc = (XLC**2 * MATMUL(TRANSPOSE(BPHI),MATMUL(A,BPHI))
     1     + MATMUL(N_MAT, TRANSPOSE(N_MAT))/((1-d_gp)+ZETA)) * WEIGHT
C       Kloc =  (XLC**2 * MATMUL(TRANSPOSE(BPHI),BPHI)
C      1     + MATMUL(N_MAT, TRANSPOSE(N_MAT))/((1-d_gp)+ZETA)) * WEIGHT
C       Kloc = (XLC**2 * MATMUL(TRANSPOSE(BPHI),MATMUL(A,BPHI))
C      1     + MATMUL(N_MAT, TRANSPOSE(N_MAT))) * WEIGHT
      K = K + Kloc

      R = R - MATMUL(Kloc, dnodal) 
      R = R + N_MAT(:,1) * d_gp/((1-d_gp)+ZETA) * WEIGHT
C       R = R + N_MAT(:,1) * d_gp * WEIGHT

      RETURN
      END

!     ==================================================================      
!     GAUSS POINT
!     ==================================================================       
      SUBROUTINE GAUSSPOINT(GP, XI)
      IMPLICIT NONE

      INTEGER GP
      DOUBLE PRECISION XI(3)
      DOUBLE PRECISION SQRT3

      PARAMETER (SQRT3 = 0.5773502691896257D0)

      SELECT CASE (GP)
      CASE (1)
        XI(1) = -SQRT3
        XI(2) = -SQRT3
        XI(3) = -SQRT3
      CASE (2)
        XI(1) =  SQRT3
        XI(2) = -SQRT3
        XI(3) = -SQRT3
      CASE (3)
        XI(1) =  SQRT3
        XI(2) =  SQRT3
        XI(3) = -SQRT3
      CASE (4)
        XI(1) = -SQRT3
        XI(2) =  SQRT3
        XI(3) = -SQRT3
      CASE (5)
        XI(1) = -SQRT3
        XI(2) = -SQRT3
        XI(3) =  SQRT3
      CASE (6)
        XI(1) =  SQRT3
        XI(2) = -SQRT3
        XI(3) =  SQRT3
      CASE (7)
        XI(1) =  SQRT3
        XI(2) =  SQRT3
        XI(3) =  SQRT3
      CASE (8)
        XI(1) = -SQRT3
        XI(2) =  SQRT3
        XI(3) =  SQRT3
      END SELECT

      RETURN
      END

!     ==================================================================      
!     Shape function
!     ==================================================================       
      SUBROUTINE SHAPEFUN(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      DOUBLE PRECISION AN(8),dNdxi(8,3)
      DOUBLE PRECISION XI(3)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0,EIGHT=8.D0)

!     Values of shape functions as a function of local coord.
      AN(1) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(2) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(3) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(4) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(5) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(6) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(7) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE+XI(3))
      AN(8) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE+XI(3))
      
!     Derivatives of shape functions respect to local coordinates
      DO I=1,8
        DO J=1,3
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(1,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(1,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(2,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(2,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(2,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(3,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(3,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(3,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(4,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(4,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(4,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
      dNdxi(5,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(5,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(5,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(6,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(6,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(6,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(7,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(7,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(7,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(8,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(8,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(8,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
     
      RETURN
      END

!     ==================================================================      
!     JACOBIAN
!     ==================================================================       
      SUBROUTINE JACOBIAN(COORDS, dNdxi, JAC, DETJ)
      IMPLICIT NONE

      DOUBLE PRECISION COORDS(3,8)         ! Physical coordinates of 8 nodes
      DOUBLE PRECISION dNdxi(8,3)          ! Shape function derivatives w.r.t ξ,η,ζ
      DOUBLE PRECISION JAC(3,3), DETJ
      INTEGER I, J, K

      ! Initialize Jacobian matrix
      DO I = 1, 3
        DO J = 1, 3
          JAC(I,J) = 0.D0
          DO K = 1, 8
            JAC(I,J) = JAC(I,J) + COORDS(I,K) * dNdxi(K,J)
          END DO
        END DO
      END DO

      ! Compute determinant of Jacobian
      DETJ = JAC(1,1)*(JAC(2,2)*JAC(3,3) - JAC(2,3)*JAC(3,2)) - 
     1        JAC(1,2)*(JAC(2,1)*JAC(3,3) - JAC(2,3)*JAC(3,1)) + 
     2       JAC(1,3)*(JAC(2,1)*JAC(3,2) - JAC(2,2)*JAC(3,1))

      RETURN
      END

!     ==================================================================      
!     INVERSE3X3
!     ==================================================================       
      SUBROUTINE INVERSE3X3(JAC, JINV)
      IMPLICIT NONE

      DOUBLE PRECISION JAC(3,3), JINV(3,3)
      DOUBLE PRECISION DET

      ! Compute determinant of JAC
      DET = JAC(1,1)*(JAC(2,2)*JAC(3,3) - JAC(2,3)*JAC(3,2)) - 
     1       JAC(1,2)*(JAC(2,1)*JAC(3,3) - JAC(2,3)*JAC(3,1)) + 
     2       JAC(1,3)*(JAC(2,1)*JAC(3,2) - JAC(2,2)*JAC(3,1))

      IF (ABS(DET) .LT. 1.0D-12) THEN
         WRITE(*,*) 'Singular Jacobian — determinant too small'
         CALL XIT
      END IF

      ! Inverse using adjugate formula
      JINV(1,1) =  (JAC(2,2)*JAC(3,3) - JAC(2,3)*JAC(3,2)) / DET
      JINV(1,2) = -(JAC(1,2)*JAC(3,3) - JAC(1,3)*JAC(3,2)) / DET
      JINV(1,3) =  (JAC(1,2)*JAC(2,3) - JAC(1,3)*JAC(2,2)) / DET

      JINV(2,1) = -(JAC(2,1)*JAC(3,3) - JAC(2,3)*JAC(3,1)) / DET
      JINV(2,2) =  (JAC(1,1)*JAC(3,3) - JAC(1,3)*JAC(3,1)) / DET
      JINV(2,3) = -(JAC(1,1)*JAC(2,3) - JAC(1,3)*JAC(2,1)) / DET

      JINV(3,1) =  (JAC(2,1)*JAC(3,2) - JAC(2,2)*JAC(3,1)) / DET
      JINV(3,2) = -(JAC(1,1)*JAC(3,2) - JAC(1,2)*JAC(3,1)) / DET
      JINV(3,3) =  (JAC(1,1)*JAC(2,2) - JAC(1,2)*JAC(2,1)) / DET

      RETURN
      END

C******************************************************************************
C CUNTZE FAILRE CRITERIA******************************************************
C******************************************************************************
      SUBROUTINE CUNTZE(S,S21,XT,XC,YT,YC,FF1,FF2,IFF1,IFF2,IFF3,
     1       FE_EFF,MDOT)   

      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DOUBLE PRECISION
     1 S21,XT,XC,YT,YC,S(6),
     2 FF1,FF2,IFF1,IFF2,IFF3,FE_EFF,
     3 B1,B2, MDOT,I_23_5,TSTRAIN(6),STRESS(6)    
      PARAMETER (ZERO=0.D0, ONE=1.D0)
C       estress was defined in abaqus notation - 11, 22, 33, 12, 23, 13
C       tot_strain also
      SIG11=S(1)
      SIG22=S(2)
      SIG33=S(3)
      SIG12=S(4)
      SIG13=S(5)
      SIG23=S(6) 

C       TSN11=TSTRAIN(1)
C       TSN22=TSTRAIN(2)
C       TSN33=TSTRAIN(3)
C       TSN12=TSTRAIN(4)
C       TSN13=TSTRAIN(5)
C       TSN23=TSTRAIN(6)
      
      B1 = 1.27D0
      B2 = 0.62D0
      MDOT = 2.2D0
      VF = 1.0D0
      ! INITIALIZE ALL AS ZERO
C       FF1=0.D0; FF2=0.D0
      
      IF(SIG11 .GE. ZERO) THEN
          FF1 = VF*SIG11/XT
      ELSE
          FF2 = ABS(VF*SIG11)/XC
      END IF
C       Wrong expressions. Check paper Petersen 2016 
      IF(SIG22 .GE. ZERO) THEN
        IFF1 = (SIG22 + SIG33 + SQRT((SIG22 - SIG33)**2 
     1  + 4. * SIG23**2)) / (2.0d0*YT)
      ELSE 
        IFF2 = ABS((B1 * SQRT((SIG22 - SIG33)**2 + 4.D0 * SIG23**2) 
     1  + (B1 - 1.D0) * (SIG22 + SIG33)) /YC)
      END IF

        I_23_5 = 2.*SIG22*SIG12**2 + 2.*SIG33*SIG13**2 
     1  + 4.*SIG23*SIG13*SIG12
C         IFF3 = SQRT((SQRT(B2**2 * I_23_5**2 + 4. * S21**2 *
C      1  (SIG13**2 +SIG12**2)**2) + B2*I_23_5)/(2.*S21**3))
        X1 = B2**2 * I_23_5**2 + 4. *S21**2 *(SIG13**2 +SIG12**2)**2
        X2 = B2*I_23_5
        IF (ABS(X2).LT.1.D-8) X2 = 0.D0
        X3 = (SQRT(X1) + X2) / (2.D0 * S21**3)
        IF (ABS(X3).LT.1.D-8) X3 = 0.D0
        IF (X3.LT.0.D0) X3 = 0.D0
        IFF3 = SQRT(X3)

   
        FE_EFF = (FF1**MDOT+FF2**MDOT+IFF1**MDOT+IFF2**MDOT+IFF3**MDOT)
     1       **(1/MDOT)
C         FE_EFF = IFF1**MDOT + IFF2**MDOT + IFF3**MDOT
        IF(FE_EFF.NE.FE_EFF) THEN
          PRINT*, "INSIDE THE CUNTZE FUNCTION! THIS MEANS SOMETHING HAS GONE WRONG"
          PRINT*, "FE_EFF IS NOT DEFINED"
          PRINT*, "IFF1, IFF2, IFF3",  IFF1, IFF2, IFF3
C           PRINT*, "MDOT, B1, B2 = ", MDOT, B1, B2
          PRINT*, "S = ", S
          PRINT*, "I_23_5 = ", I_23_5
          PRINT*, "X1, X2, X3 = ", X1, X2, X3
          CALL XIT
        ENDIF

        

      RETURN
      END SUBROUTINE CUNTZE
C********************************************************************     
!     
C******************************************************************************
C CALCULATION OF DAMAGE PARAMETER**********************************************
C******************************************************************************
C       SUBROUTINE DAMAGE(H,PHI_C,PHI_S,D)

C       IMPLICIT DOUBLE PRECISION (A-H,O-Z)  
 
C       DOLD = D
C       D = (H - PHI_C)/(H - PHI_C + PHI_S)

C       IF(D .GT. 0.99D0) then
C         D = 0.99D0
C       ENDIF
      
C       IF(D .LT. 0.0D0) then
C         PRINT*,'NEGATIVE D',D,H,PHI_C,PHI_S
C         D = DOLD
C         CALL XIT
C       ENDIF

C       IF(DOLD .GT. D) then
C         PRINT*,'LOW D',D,H,PHI_C,PHI_S
C         D = DOLD
C       ENDIF

C       RETURN
C       END 

      SUBROUTINE DAMAGE(H,PHI_C,PHI_S,D)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION H,PHI_C,PHI_S,D,DOLD
      DOUBLE PRECISION F_POS,DTRIAL,DEN
      PARAMETER (ZERO=0.D0, ONE=1.D0)

      DOLD = D

C     Macaulay bracket: <H - PHI_C>
      F_POS = 0.5D0 * ( (H - PHI_C) + ABS(H - PHI_C) )

C     If below threshold, no damage growth
      IF (F_POS .LE. ZERO) THEN
         D = DOLD
         RETURN
      END IF

C     Algebraic damage law
      DEN = F_POS + PHI_S

      IF (DEN .LE. ZERO) THEN
         PRINT*, 'BAD DAMAGE DENOMINATOR', H, PHI_C, PHI_S, F_POS, DEN
         D = DOLD
         CALL XIT
      END IF

      DTRIAL = F_POS / DEN

C     Bound damage
      DTRIAL = MIN(MAX(0.0D0,DTRIAL),0.999D0)

C     Irreversibility: d_{n+1} >= d_n
      D = MAX(DOLD,DTRIAL)

      RETURN
      END

C      
C**********************************************************       
      SUBROUTINE INVERSE(AA,C,N)
!============================================================
! INVERSE MATRIX
! METHOD: BASED ON DOOLITTLE LU FACTORIZATION FOR AX=B
! ALEX G. DECEMBER 2009
!-----------------------------------------------------------
! INPUT ...
! A(N,N) - ARRAY OF COEFFICIENTS FOR MATRIX A
! N      - DIMENSION
! OUTPUT ...
! C(N,N) - INVERSE MATRIX OF A
! COMMENTS ...
! THE ORIGINAL MATRIX A(N,N) WILL BE DESTROYED 
! DURING THE CALCULATION
!===========================================================
        IMPLICIT NONE 
        INTEGER N
        DOUBLE PRECISION A(N,N), C(N,N), AA(N,N)
        DOUBLE PRECISION L(N,N), U(N,N), B(N), D(N), X(N)
        DOUBLE PRECISION COEFF
        INTEGER I, J, K

! STEP 0: INITIALIZATION FOR MATRICES L AND U AND B
! FORTRAN 90/95 ALOOWS SUCH OPERATIONS ON MATRICES
        L=0.0
        U=0.0
        B=0.0
        A = AA
! STEP 1: FORWARD ELIMINATION
        DO K=1, N-1
            DO I=K+1,N
              COEFF=A(I,K)/A(K,K)
              L(I,K) = COEFF
                DO J=K+1,N
                   A(I,J) = A(I,J)-COEFF*A(K,J)
                END DO
             END DO
        END DO

! STEP 2: PREPARE L AND U MATRICES 
! L MATRIX IS A MATRIX OF THE ELIMINATION COEFFICIENT
! + THE DIAGONAL ELEMENTS ARE 1.0
       DO I=1,N
         L(I,I) = 1.0
       END DO
! U MATRIX IS THE UPPER TRIANGULAR PART OF A
        DO J=1,N
          DO I=1,J
              U(I,J) = A(I,J)
          END DO
        END DO

! STEP 3: COMPUTE COLUMNS OF THE INVERSE MATRIX C
          DO K=1,N
              B(K)=1.0
              D(1) = B(1)
! STEP 3A: SOLVE LD=B USING THE FORWARD SUBSTITUTION
                  DO I=2,N
                      D(I)=B(I)
                      DO J=1,I-1
                          D(I) = D(I) - L(I,J)*D(J)
                      END DO
                  END DO
! STEP 3B: SOLVE UX=D USING THE BACK SUBSTITUTION
                  X(N)=D(N)/U(N,N)
                  DO I = N-1,1,-1
                       X(I) = D(I)
                       DO J=N,I+1,-1
                          X(I)=X(I)-U(I,J)*X(J)
                       END DO
                       X(I) = X(I)/U(I,I)
                  END DO
! STEP 3C: FILL THE SOLUTIONS X(N) INTO COLUMN K OF C
                  DO I=1,N
                      C(I,K) = X(I)
                  END DO
                  B(K)=0.0
          END DO
       RETURN   
       END SUBROUTINE INVERSE
C************************************************************
      SUBROUTINE TV_MAT_V(TV1,MAT,V2,C1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      
      DOUBLE PRECISION C1,TV1(6),V2(6),MAT(6,6),F(1,1),MV1(1,6),MV2(6,1)

        DO I = 1,6
          MV1(1,I) = TV1(I)
          MV2(I,1) = V2(I)
        ENDDO
      
        F = MATMUL(MATMUL(MV1,MAT),MV2)
        C1 = F(1,1)

      RETURN
      END SUBROUTINE TV_MAT_V
C**********************************************************************
      SUBROUTINE onem(A)

C     THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C     3 BY 3 MATRIX [A]
C**********************************************************************

        REAL*8 A(3,3)
        DATA ZERO/0.D0/
        DATA ONE/1.D0/

      DO 1 I=1,3
        DO 1 J=1,3
          IF (I .EQ. J) THEN
              A(I,J) = 1.0
            ELSE
              A(I,J) = 0.0
            ENDIF
1       CONTINUE

      RETURN
      END

      subroutine dyadic(vector1,vector2, vlen, dyadicprod)
                  
        integer  vlen, i, j
        real*8 vector1(vlen),vector2(vlen)
        real*8 dyadicprod(vlen,vlen)

        do i = 1, vlen
              do j = 1, vlen
              dyadicprod(i,j) = vector1(i) * vector2(j)
              end do
        end do

        return
      end subroutine dyadic 
C=======================================================================
C  CALC_YS_FULL  (CORRECTED)
C  --------------------------------
C  Uses PIECEWISE-LINEAR hardening curves provided by you (E_* vs Y_*)
C  instead of any fitted equation.
C
C  INPUT:
C    EQPSN   : equivalent plastic strain
C    SIG(6)  : stress (11,22,33,12,13,23)
C
C  OUTPUT:
C    F       : yield function value  f = 0.5*s^T*A*s + a·s - 1
C    AM(6,6) : matrix A(EQPSN,SIG sign-branch)
C    AT(1,6) : vector a(EQPSN,SIG sign-branch)
C    DFDEP   : df/d(EQPSN) for consistent tangent
C
C  Notes:
C   - Compression data Y_UC, Y_BC are NEGATIVE (as provided). This is OK:
C     the Vogler-form uses these signed values in ALPHA3 and the squares.
C   - Outside tabulated EQPS range: LINEAR EXTRAPOLATION using end segments.
C=======================================================================

      SUBROUTINE CALC_YS_FULL(EQPSN,SIG,F,AM,AT,DFDEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION EQPSN,SIG(6),F,DFDEP
      DOUBLE PRECISION AM(6,6),AT(1,6)
      DOUBLE PRECISION dAM(6,6),dAT(1,6)

      INTEGER NBT,NUT,NBC,NUC,NLS,NTS
      PARAMETER (NBT=10, NUT=10, NBC=13, NUC=11, NLS=12, NTS=12)

C       DOUBLE PRECISION E_BT(NBT),Y_BT(NBT)
C       DOUBLE PRECISION E_UT(NUT),Y_UT(NUT)
C       DOUBLE PRECISION E_BC(NBC),Y_BC(NBC)
C       DOUBLE PRECISION E_UC(NUC),Y_UC(NUC)
C       DOUBLE PRECISION E_LS(NLS),Y_LS(NLS)
C       DOUBLE PRECISION E_TS(NTS),Y_TS(NTS)

C C----- your tabulated data (as-is)
C       DATA E_BT / 0.D0, 0.001118D0, 0.002294D0, 0.003765D0, 0.005294D0,
C      &            0.006588D0, 0.007941D0, 0.009471D0, 0.011D0, 0.012765D0 /
C       DATA Y_BT / 30.14354D0, 55.98086D0, 70.33493D0, 78.94737D0, 84.689D0,
C      &            88.99522D0, 91.86603D0, 96.17225D0, 97.60766D0, 101.9139D0 /

C       DATA E_UT / 0.D0, 0.000941D0, 0.002176D0, 0.003765D0, 0.005529D0,
C      &            0.007353D0, 0.009118D0, 0.011235D0, 0.012882D0, 0.014294D0 /
C       DATA Y_UT / 33.1579D0, 70.33493D0, 93.30144D0, 107.6555D0, 116.2679D0,
C      &            123.445D0, 129.1866D0, 134.9282D0, 137.799D0, 139.799D0 /

C       DATA E_BC / 0.D0, 0.000471D0, 0.001471D0, 0.003235D0, 0.005882D0,
C      &            0.009176D0, 0.012471D0, 0.015941D0, 0.019294D0, 0.022706D0,
C      &            0.026118D0, 0.029706D0, 0.032059D0 /
C       DATA Y_BC / -236.842D0, -308.612D0, -381.818D0, -449.282D0, -496.651D0,
C      &            -526.794D0, -549.761D0, -564.115D0, -578.469D0, -589.952D0,
C      &            -598.565D0, -605.742D0, -608.612D0 /

C       DATA E_UC / 0.D0, 0.000882D0, 0.002235D0, 0.004235D0, 0.006294D0,
C      &            0.008941D0, 0.011471D0, 0.014D0, 0.016588D0, 0.018706D0,
C      &            0.020824D0 /

C C       DATA Y_UC / -106.22D0, -153.12D0, -186.603D0, -216.746D0, -229.665D0,
C C      &            -241.148D0, -248.325D0, -255.502D0, -261.244D0, -261.244D0,
C C      &            -261.244D0 /

C       DATA Y_UC / -105.902D0, -152.673D0, -184.939D0, -214.176D0, -226.271D0, 
C      &            -236.884D0, -243.217D0, -249.514D0, -254.376D0, -258.806D0,
C      &            -261.424D0 /

C       DATA E_LS / 0.D0, 0.001007D0, 0.002937D0, 0.007301D0, 0.012755D0,
C      &            0.018294D0, 0.025762D0, 0.03449D0, 0.042042D0, 0.047832D0,
C      &            0.054042D0, 0.059413D0 /
C       DATA Y_LS / 30.04418D0, 42.76878D0, 65.39028D0, 84.12371D0, 94.19735D0,
C      &            100.5596D0, 106.215D0, 110.6333D0, 113.9912D0, 116.2887D0,
C      &            118.2327D0, 119.8233D0 /

C       DATA E_TS / 0.D0, 0.001343D0, 0.003273D0, 0.006545D0, 0.009986D0,
C      &            0.017538D0, 0.024084D0, 0.03172D0, 0.037678D0, 0.04414D0,
C      &            0.050098D0, 0.05958D0 /
C       DATA Y_TS / 30.04418D0, 50.01473D0, 65.03682D0, 77.40795D0, 85.18409D0,
C      &            94.72754D0, 99.85272D0, 103.9175D0, 106.7452D0, 109.0427D0,
C      &            111.1635D0, 113.8144D0 /

C----- locals
      DOUBLE PRECISION YBT,YUT,YBC,YUC,YLS,YTS
      DOUBLE PRECISION dYBT,dYUT,dYBC,dYUC,dYLS,dYTS
      DOUBLE PRECISION ALPHA1,ALPHA2,ALPHA32,ALPHA3
      DOUBLE PRECISION dA1,dA2,dA32,dA3
      DOUBLE PRECISION DENOM,NUMER,dDEN,dNUM
      DOUBLE PRECISION FS1,FS1d,YSTRESS,SIGSUM
      INTEGER I,J

C---------------- reset matrices
      DO I=1,6
        AT(1,I)=0.D0
        dAT(1,I)=0.D0
        DO J=1,6
          AM(I,J)=0.D0
          dAM(I,J)=0.D0
        END DO
      END DO

C       EQPSN_U = 1.346D0 * EQPSN
C       EQPSN_S = 2.D0 * EQPSN

C C---------------- tabulated Y(EQPSN) + slopes dY/dEQPSN
C       CALL PWL1D(EQPSN, E_BT, Y_BT, NBT, YBT, dYBT)
C       CALL PWL1D(EQPSN_U, E_UT, Y_UT, NUT, YUT, dYUT)
C       CALL PWL1D(EQPSN, E_BC, Y_BC, NBC, YBC, dYBC)
C       CALL PWL1D(EQPSN_U, E_UC, Y_UC, NUC, YUC, dYUC)
C       CALL PWL1D(EQPSN_S, E_LS, Y_LS, NLS, YLS, dYLS)
C       CALL PWL1D(EQPSN_S, E_TS, Y_TS, NTS, YTS, dYTS)

      CALL CALC_HARDENING_VOCE_LINEAR(EQPSN,
     1     YBT,dYBT,YUT,dYUT,YBC,dYBC,YUC,dYUC,YLS,dYLS,YTS,dYTS)

C---------------- alpha1, alpha2 based on transverse/longitudinal shear curves
      ALPHA1 = 1.D0/(YTS*YTS)
      dA1    = -2.D0*dYTS/(YTS**3)

      ALPHA2 = 1.D0/(YLS*YLS)
      dA2    = -2.D0*dYLS/(YLS**3)

C---------------- tension/compression switch
C   Keep your original logic: use SIG(2)+SIG(3) as branch indicator.
      SIGSUM = SIG(2) + SIG(3)

      IF (SIGSUM .GE. 0.D0) THEN
C       Tension branch uses YUT, YBT
        DENOM = YUT**2 - 2.D0*YBT*YUT
        NUMER = 1.D0 - (YUT/(2.D0*YBT)) - ALPHA1*(YUT**2/4.D0)
        ALPHA32 = NUMER/DENOM

        dNUM = -( (dYUT/(2.D0*YBT)) - (YUT*dYBT)/(2.D0*YBT*YBT) )
     &          - dA1*(YUT**2/4.D0) - ALPHA1*(YUT*dYUT/2.D0)
        dDEN = 2.D0*YUT*dYUT - 2.D0*(dYBT*YUT + YBT*dYUT)
        dA32 = (dNUM*DENOM - NUMER*dDEN)/(DENOM*DENOM)

        ALPHA3 = 1.D0/(2.D0*YBT) - 2.D0*ALPHA32*YBT
        dA3    = -(dYBT/(2.D0*YBT*YBT))
     &           - 2.D0*(dA32*YBT + ALPHA32*dYBT)

      ELSE
C       Compression branch uses YUC, YBC (NEGATIVE values allowed)
        DENOM = YUC**2 - 2.D0*YBC*YUC
        NUMER = 1.D0 - (YUC/(2.D0*YBC)) - ALPHA1*(YUC**2/4.D0)
        ALPHA32 = NUMER/DENOM

        dNUM = -( (dYUC/(2.D0*YBC)) - (YUC*dYBC)/(2.D0*YBC*YBC) )
     &          - dA1*(YUC**2/4.D0) - ALPHA1*(YUC*dYUC/2.D0)
        dDEN = 2.D0*YUC*dYUC - 2.D0*(dYBC*YUC + YBC*dYUC)
        dA32 = (dNUM*DENOM - NUMER*dDEN)/(DENOM*DENOM)

        ALPHA3 = 1.D0/(2.D0*YBC) - 2.D0*ALPHA32*YBC
        dA3    = -(dYBC/(2.D0*YBC*YBC))
     &           - 2.D0*(dA32*YBC + ALPHA32*dYBC)
      END IF

C---------------- A matrix: (22-33 coupling) + shears
      AM(2,2) = 0.5D0*ALPHA1 + 2.D0*ALPHA32
      AM(2,3) =-0.5D0*ALPHA1 + 2.D0*ALPHA32
      AM(3,2) = AM(2,3)
      AM(3,3) = 0.5D0*ALPHA1 + 2.D0*ALPHA32

C shear terms: 12,13 use ALPHA2 ; 23 uses ALPHA1 (as in your original)
      AM(4,4) = 2.D0*ALPHA2
      AM(5,5) = 2.D0*ALPHA2
      AM(6,6) = 2.D0*ALPHA1

C---------------- a vector
      AT(1,2) = ALPHA3
      AT(1,3) = ALPHA3

C---------------- yield value
      CALL TV_MAT_V(SIG,AM,SIG,FS1)
      YSTRESS = AT(1,2)*SIG(2) + AT(1,3)*SIG(3)
      F = 0.5D0*FS1 + YSTRESS - 1.D0

C---------------- df/dep
      dAM(2,2) = 0.5D0*dA1 + 2.D0*dA32
      dAM(2,3) =-0.5D0*dA1 + 2.D0*dA32
      dAM(3,2) = dAM(2,3)
      dAM(3,3) = 0.5D0*dA1 + 2.D0*dA32

      dAM(4,4) = 2.D0*dA2
      dAM(5,5) = 2.D0*dA2
      dAM(6,6) = 2.D0*dA1

      dAT(1,2) = dA3
      dAT(1,3) = dA3

      CALL TV_MAT_V(SIG,dAM,SIG,FS1d)
      DFDEP = 0.5D0*FS1d + dAT(1,2)*SIG(2) + dAT(1,3)*SIG(3)

      RETURN
      END SUBROUTINE CALC_YS_FULL


C=======================================================================
C  PWL1D: Piecewise-linear interpolation + slope
C  - Extrapolates linearly outside the table using end segments
C=======================================================================
      SUBROUTINE PWL1D(X, XA, YA, N, Y, DYDX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N, I
      DOUBLE PRECISION X, XA(N), YA(N), Y, DYDX
      DOUBLE PRECISION X0,X1,Y0,Y1,T

C----- guard N>=2
      IF (N .LT. 2) THEN
        Y    = YA(1)
        DYDX = 0.D0
        RETURN
      END IF

C----- left extrapolation
      IF (X .LE. XA(1)) THEN
        X0 = XA(1); X1 = XA(2)
        Y0 = YA(1); Y1 = YA(2)
        DYDX = (Y1-Y0)/(X1-X0)
        Y    = Y0 + DYDX*(X - X0)
        RETURN
      END IF

C----- right extrapolation
      IF (X .GE. XA(N)) THEN
        X0 = XA(N-1); X1 = XA(N)
        Y0 = YA(N-1); Y1 = YA(N)
        DYDX = (Y1-Y0)/(X1-X0)
        Y    = Y1 + DYDX*(X - X1)
        RETURN
      END IF

C----- find interval
      DO I=1,N-1
        IF (X .GE. XA(I) .AND. X .LE. XA(I+1)) THEN
          X0 = XA(I);   X1 = XA(I+1)
          Y0 = YA(I);   Y1 = YA(I+1)
          DYDX = (Y1-Y0)/(X1-X0)
          T = (X - X0)/(X1 - X0)
          Y = (1.D0-T)*Y0 + T*Y1
          RETURN
        END IF
      END DO

C----- fallback (should not happen)
      Y    = YA(N)
      DYDX = 0.D0
      RETURN
      END SUBROUTINE PWL1D

      SUBROUTINE CALC_HARDENING_VOCE_LINEAR(EQPSN,
     1     YBT,dYBT,YUT,dYUT,YBC,dYBC,YUC,dYUC,YLS,dYLS,YTS,dYTS)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION EQPSN, EQPSN_U, EQPSN_S
      DOUBLE PRECISION YBT,dYBT,YUT,dYUT,YBC,dYBC,YUC,dYUC,YLS,dYLS,YTS,dYTS
      DOUBLE PRECISION TAIL

C     ------------------------------------------------------------
C     Mapping used already in your model
C     ------------------------------------------------------------
      EQPSN_U = 1.346D0 * EQPSN
      EQPSN_S = 2.D0   * EQPSN

C     ============================================================
C     BT
C     Y_BT = 30.1624698244 + 55.7209902798*(1-exp(-555.4995779317*EQPS))
C    &      + 1869.8089624996*max(0,EQPS-0.0042565399)
C     ============================================================
      TAIL = MAX(0.D0, EQPSN - 0.0042565399D0)
      YBT  = 30.1624698244D0
     &     + 55.7209902798D0*(1.D0-EXP(-555.4995779317D0*EQPSN))
     &     + 1869.8089624996D0*TAIL

      dYBT = 55.7209902798D0*555.4995779317D0
     &     * EXP(-555.4995779317D0*EQPSN)

      IF (EQPSN .GT. 0.0042565399D0) THEN
         dYBT = dYBT + 1869.8089624996D0
      END IF

C     ============================================================
C     UT
C     ============================================================
      TAIL = MAX(0.D0, EQPSN_U - 0.0038352261D0)
      YUT  = 33.4948076825D0
     &     + 82.5188770914D0*(1.D0-EXP(-606.5605258063D0*EQPSN_U))
     &     + 2397.0100057674D0*TAIL

      dYUT = 82.5188770914D0*606.5605258063D0
     &     * EXP(-606.5605258063D0*EQPSN_U)

      IF (EQPSN_U .GT. 0.0038352261D0) THEN
         dYUT = dYUT + 2397.0100057674D0
      END IF

C     chain rule because YUT depends on EQPSN_U = 1.346*EQPSN
      dYUT = 1.346D0 * dYUT

C     ============================================================
C     BC
C     ============================================================
      TAIL = MAX(0.D0, EQPSN - 0.0058819884D0)
      YBC  = -246.5699096531D0
     &     - 278.8960931346D0*(1.D0-EXP(-428.9581028583D0*EQPSN))
     &     - 3479.6198624682D0*TAIL

      dYBC = -278.8960931346D0*428.9581028583D0
     &     * EXP(-428.9581028583D0*EQPSN)

      IF (EQPSN .GT. 0.0058819884D0) THEN
         dYBC = dYBC - 3479.6198624682D0
      END IF

C     ============================================================
C     UC
C     ============================================================
      TAIL = MAX(0.D0, EQPSN_U - 0.0045355935D0)
      YUC  = -107.1683887260D0
     &     - 122.2572702495D0*(1.D0-EXP(-480.5773669390D0*EQPSN_U))
     &     - 2040.9896514827D0*TAIL

      dYUC = -122.2572702495D0*480.5773669390D0
     &     * EXP(-480.5773669390D0*EQPSN_U)

      IF (EQPSN_U .GT. 0.0045355935D0) THEN
         dYUC = dYUC - 2040.9896514827D0
      END IF

C     chain rule
      dYUC = 1.346D0 * dYUC

C     ============================================================
C     LS
C     ============================================================
      TAIL = MAX(0.D0, EQPSN_S - 0.0121057874D0)
      YLS  = 30.1176886553D0
     &     + 69.2338015144D0*(1.D0-EXP(-219.2624831020D0*EQPSN_S))
     &     + 459.5754088487D0*TAIL

      dYLS = 69.2338015144D0*219.2624831020D0
     &     * EXP(-219.2624831020D0*EQPSN_S)

      IF (EQPSN_S .GT. 0.0121057874D0) THEN
         dYLS = dYLS + 459.5754088487D0
      END IF

C     chain rule
      dYLS = 2.D0 * dYLS

C     ============================================================
C     TS
C     ============================================================
      TAIL = MAX(0.D0, EQPSN_S - 0.0099859873D0)
      YTS  = 31.7034344390D0
     &     + 61.2719745511D0*(1.D0-EXP(-227.1910881684D0*EQPSN_S))
     &     + 453.8647523823D0*TAIL

      dYTS = 61.2719745511D0*227.1910881684D0
     &     * EXP(-227.1910881684D0*EQPSN_S)

      IF (EQPSN_S .GT. 0.0099859873D0) THEN
         dYTS = dYTS + 453.8647523823D0
      END IF

C     chain rule
      dYTS = 2.D0 * dYTS

      RETURN
      END