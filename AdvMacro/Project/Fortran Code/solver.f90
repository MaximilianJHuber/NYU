MODULE solver
  USE modeldef
  USE params
  USE stat_helper
  USE glob
  USE output
  
  IMPLICIT NONE
  
CONTAINS
  !----- CODE TO COMPUTE FIXED POINT IN S(X,Y,Z)
  SUBROUTINE solveModel(m,p,e)
    USE modeldef
    USE params
    IMPLICIT NONE
    
    TYPE (model), TARGET, INTENT(inout) :: m
    TYPE (ExogenousParameters), INTENT(inout)  :: p
    ! 0 for convergence 1 no convergence in S 2 no convergence in H
    INTEGER, INTENT(inout) :: e  
    INTEGER i, j, x, y, z
    REAL(wp) tol 
    
    !--- Define any temporaty arrays used during computation here
    REAL(wp), ALLOCATABLE :: ES_xyz(:,:,:)
    REAL(wp), ALLOCATABLE :: EBB_xz(:,:)
    REAL(wp), ALLOCATABLE :: PmCmB_xyz(:,:,:)
    REAL(wp), ALLOCATABLE :: S0_xy(:,:) 
    REAL(wp), ALLOCATABLE :: H0p_xy(:,:)
    REAL(wp), ALLOCATABLE :: U_x(:)
    REAL(wp), ALLOCATABLE :: Jtmp1_y(:)
    REAL(wp), ALLOCATABLE :: Jtmp2_y(:)
    REAL(wp), ALLOCATABLE :: J_y(:)
    REAL(wp), ALLOCATABLE :: v_y(:)
    REAL(wp), ALLOCATABLE :: q_y_v_y(:)
    REAL(wp), ALLOCATABLE :: Htmp1_xy(:,:)
    REAL(wp), ALLOCATABLE :: Htmp2_xy(:,:)
    REAL(wp), ALLOCATABLE :: Htmp3_xy(:,:)
    REAL(wp), ALLOCATABLE :: H0tmp_xy(:,:)

    REAL(wp) theta, L, Jc, V, Matches, q_theta 
    
    !-------- DEFINING POINTERS FOR THE MODELS ---------
    TYPE (model), TARGET 		:: m_tmp 
    TYPE (model), POINTER 		:: m_old , m_new , m_swap
    
    !--- Initialize the model with the current parameter values
    CALL Init(m_tmp,p) 
    
    !--- Allocate the temporary arrays
    ALLOCATE ( ES_xyz(p%nx,p%ny,p%nz) )
    ALLOCATE ( EBB_xz(p%nx,p%nz) )
    ALLOCATE ( PmCmB_xyz(p%nx,p%ny,p%nz) )
    ALLOCATE ( S0_xy(p%nx,p%ny) ) 
    ALLOCATE ( H0p_xy(p%nx,p%ny) )
    ALLOCATE ( U_x(p%nx) )
    ALLOCATE ( Jtmp1_y(p%ny) )
    ALLOCATE ( Jtmp2_y(p%ny) )
    ALLOCATE ( J_y(p%ny) )
    ALLOCATE ( v_y(p%ny) )
    ALLOCATE ( q_y_v_y(p%ny) )
    ALLOCATE ( Htmp1_xy(p%nx,p%ny) )
    ALLOCATE ( Htmp2_xy(p%nx,p%ny) )
    ALLOCATE ( Htmp3_xy(p%nx,p%ny) )
    ALLOCATE ( H0tmp_xy(p%nx,p%ny) )
    
    !-------- PREPARATION BEFORE THE LOOP ---------
    m_new=>m_tmp
    m_old=>m

    ! iterate on S(x,y,z)
    ! calculate p(x,y,z) - b(x,z)
    CALL PmCmBxyz( m_old%p_xyz, m_old%b_xz, PmCmB_xyz )
    
    e = 1  ! initialize flag to non-convergence

    DO i = 1,p%MaxIter
       
       ! calculate E[max(S(x,y,z2),0)|z]
       CALL ESxyz( m_old%S_xyz, m_old%Q_zz, ES_xyz )
       
       ! update iteration of S(x,y,z)
       m_new%S_xyz = PmCmB_xyz + ((1.0-p%delta)/(1.0+p%r)) * ES_xyz 
    
       tol = MAXVAL( ABS( m_new%S_xyz - m_old%S_xyz ) )
       
       IF ( tol .LE. p%tol ) THEN
          m_swap => m_new
          m_new => m_old
          m_old => m_swap
          
          EXIT
       END IF
       m_swap => m_new
       m_new => m_old
       m_old => m_swap
       
    END DO

    IF ( i .LT. p%MaxIter ) THEN
       e = 0 ! S(x,y,z,) converged 
       IF ( GLOBAL_DISPLAY_SWITCH ) THEN
         WRITE(*,*) "S(x,y,z) converged in "
         WRITE(*,*) i
         WRITE(*,*) "itereations, Max(abs(S1(x,y,z)-S0(x,y,z))) is:"
         WRITE(*,*) tol
       END IF
    ELSE
       e = 1
    END IF

    ! If requested, calculate B(x,z)
    IF ( GLOBAL_DISPLAY_SWITCH ) THEN
       DO i = 1,p%MaxIter
          CALL EBBxz(m_old%BB_xz, m_old%Q_zz, EBB_xz)
          m_new%BB_xz = m_old%b_xz + (1.0/(1.0+p%r)) * EBB_xz
          tol = MAXVAL( ABS( m_new%BB_xz - m_old%BB_xz ) )
          IF ( tol .LE. p%tol ) THEN
             m_swap => m_new
             m_new => m_old
             m_old => m_swap
             
             EXIT
          END IF
          m_swap => m_new
          m_new => m_old
          m_old => m_swap
          
       END DO
    END IF
    
    IF ( e .EQ. 0 ) THEN
       ! only iterate on H0(x,y,z0) if S(x,y,z) converged
       
       ! Use S(x,y) at exp(z)=1
       S0_xy = m_new%S_xyz(:,:, (p%nz+1)/2 )
       
       DO i=1,p%MaxIter
          ! measure of workers who are employed after delta and z shock
          CALL Hxytp( S0_xy, m_old%H0_xy, H0p_xy, p%delta )
          ! calculae distribution of unemployed workers
          CALL Uxtp( m_old%ell_x, H0p_xy, p%wy, U_x )
          
          !--- Calculate effective worker Search
          CALL Lt( U_x, H0p_xy, p%s0, p%s1, p%wx, p%wy, L )
          
          !--- Calculate optimal vacancy creation by fimrs (y-j types)
          CALL J1y( U_x, S0_xy, p%wx, Jtmp1_y )
          CALL J2y( H0p_xy, S0_xy, p%wx, p%wy, Jtmp2_y )

          J_y = (p%s0/L) * Jtmp1_y + (p%s1/L) * Jtmp2_y
          
          CALL IntJ( J_y, p%eta, p%gamma0, p%gamma1, p%y, p%wy, Jc ) 
       
          CALL computeTheta(Jc, L, p%MatchingTechnology, &
               p%alpha, p%omega, p%eta, theta )

          V = theta*L ! by definition of theta

          CALL computeMatching(L, V, p%MatchingTechnology, &
               p%alpha, p%omega, Matches)

          ! create optimal vacancy creation v_y(y)
          CALL vy(J_y, Matches, V, p%gamma0, p%gamma1, p%eta, p%y, v_y)
          
          !--- Calculate update to H(x,y,z)

          q_theta = Matches/V

          ! since q_theta always appears multiplied by v_y it is efficient to do this up front
          q_y_v_y = q_theta * v_y
          
          CALL H1xyt( S0_xy, U_x, q_y_v_y,    Htmp1_xy )
          CALL H2xyt( S0_xy, H0p_xy, q_y_v_y, p%wy, Htmp2_xy )
          CALL H3xyt( S0_xy, H0p_xy, q_y_v_y, p%wy, Htmp3_xy )
          
          m_new%H0_xy = MAX(0.0,  H0p_xy + (p%s0/L)*Htmp1_xy &
               + (p%s1/L)*Htmp2_xy - (p%s1/L)*Htmp3_xy)
          
          tol = MAXVAL( ABS( m_new%H0_xy - m_old%H0_xy ) )

          IF (tol .LE. p%tol_H ) THEN
             m_swap => m_new
             m_new => m_old
             m_old => m_swap
             
             EXIT
          END IF
          m_swap => m_new
          m_new => m_old
          m_old => m_swap
          
       END DO
       IF ( i .LT. p%MaxIter ) THEN
          e = 0 ! S(x,y,z,) and H(x,y,0) both converged 
          IF ( GLOBAL_DISPLAY_SWITCH ) THEN
            WRITE(*,*) "H(x,y,0) converged in "
            WRITE(*,*) i
            WRITE(*,*) "itereations,  Max(abs(H1(x,y,0)-H0(x,y,0))) is:"
            WRITE(*,*) tol
          END IF
       ELSE
          e = 2
       END IF
       
       IF ( GLOBAL_DISPLAY_SWITCH ) THEN
          CALL FirmshareNewMatch( U_x, H0p_xy, q_y_v_y, S0_xy, p%s0, p%s1, p%wx, p%wy )
       END IF

    END IF

    !--- clean up temporary arrays and temporary models
    DEALLOCATE ( ES_xyz )
    DEALLOCATE ( EBB_xz )
    DEALLOCATE ( PmCmB_xyz )
    DEALLOCATE ( S0_xy )
    DEALLOCATE ( H0p_xy )
    DEALLOCATE ( U_x )
    DEALLOCATE ( Jtmp1_y )
    DEALLOCATE ( Jtmp2_y )
    DEALLOCATE ( J_y )
    DEALLOCATE ( v_y )
    DEALLOCATE ( q_y_v_y )
    DEALLOCATE ( Htmp1_xy )
    DEALLOCATE ( Htmp2_xy )
    DEALLOCATE ( Htmp3_xy )
    DEALLOCATE ( H0tmp_xy )

    CALL FreeModel(m_tmp)
  END SUBROUTINE solveModel
  
  

  SUBROUTINE ESxyz(S,Q,ES)
    !ES = E[S(x,y,z')|z]
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: S(:,:,:)
    REAL(wp), INTENT(IN) :: Q(:,:)
    REAL(wp), INTENT(INOUT) :: ES(:,:,:)
    INTEGER nx, ny, nz
    INTEGER x,y,z,z2
    
    nx = SIZE(S,1)
    ny = SIZE(S,2)
    nz = SIZE(S,3)
    
    ES = 0.0

    DO x=1,nx
       DO y=1,ny
          DO z=1,nz
             DO z2=1,nz
                ES(x,y,z) = ES(x,y,z) + MAX(0.0, S(x,y,z2)) * Q(z,z2)
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE ESxyz

  SUBROUTINE EBBxz(BB,Q,EBB)
    !ES = E[S(x,y,z')|z]
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: BB(:,:)
    REAL(wp), INTENT(IN) :: Q(:,:)
    REAL(wp), INTENT(INOUT) :: EBB(:,:)
    INTEGER nx, nz
    INTEGER x,z,z2
    
    nx = SIZE(BB,1)
    nz = SIZE(Q,1)
    
    EBB = 0.0

    DO x=1,nx
       DO z=1,nz
          DO z2=1,nz
             EBB(x,z) = EBB(x,z) + MAX(0.0, BB(x,z2)) * Q(z,z2)
          END DO
       END DO
    END DO
  END SUBROUTINE EBBxz
  
  SUBROUTINE ESxyBBxH0( S_xy, B_x, H0_xy )
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: S_xy(:,:)
    REAL(wp), INTENT(IN) :: B_x(:)
    REAL(wp), INTENT(IN) :: H0_xy(:,:)

    REAL(wp), ALLOCATABLE :: EGain(:)
    REAL(wp), ALLOCATABLE :: MaxGain(:)
    
    INTEGER nx, ny, x, y
    REAL(wp) Htmp

    nx = SIZE(S_xy,1)
    ny = SIZE(S_xy,2)

    ALLOCATE( EGain(nx) )
    ALLOCATE( MaxGain(nx) )

    EGain = 0.0
    MaxGain = 0.0 
    Htmp = 0.0

    DO x=1,nx
       Htmp = 0.0
       DO y=1,ny
          Htmp = Htmp + H0_xy(x,y)
          EGain(x) = h0_xy(x,y) * ( (S_xy(x,y)+B_x(x)) / B_x(x) )
       END DO
       EGain(x) = EGain(x) / Htmp
       MaxGain(x) = (MAXVAL(S_xy(x,:))+B_x(x)) / B_x(x)
    END DO

    WRITE(*,*) "saving to ./output/EGain_x.dat"
    CALL writetocsvvec( "./output/EGain_x.dat", EGain )
    
    WRITE(*,*) "saving to ./output/MaxGain_x.dat"
    CALL writetocsvvec( "./output/MaxGain_x.dat", MaxGain )
    
    DEALLOCATE( EGain )
    DEALLOCATE( MaxGain )
  END SUBROUTINE ESxyBBxH0

  SUBROUTINE EpxybxH0( p_xyz, b_xz, H0_xy )
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: p_xyz(:,:,:)
    REAL(wp), INTENT(IN) :: b_xz(:,:)
    REAL(wp), INTENT(IN) :: H0_xy(:,:)

    REAL(wp), ALLOCATABLE :: EGain(:)
    REAL(wp), ALLOCATABLE :: MaxGain(:,:)
    
    INTEGER nx, ny, nz, x, y, z
    REAL(wp) Htmp

    nx = SIZE(p_xyz,1)
    ny = SIZE(p_xyz,2)
    nz = SIZE(p_xyz,3)

    ALLOCATE( EGain(nx) )
    ALLOCATE( MaxGain(nx,nz) )

    EGain = 0.0
    MaxGain = 0.0 
    Htmp = 0.0

    DO x=1,nx
       Htmp = 0.0
       DO y=1,ny
          Htmp = Htmp + H0_xy(x,y)
          EGain(x) = h0_xy(x,y) * ( p_xyz(x,y,(nz+1)/2) / b_xz(x,(nz+1)/2) ) 
       END DO
       EGain(x) = EGain(x) / Htmp
       DO z=1,nz
          MaxGain(x,z) = ( b_xz(x,z) / MAXVAL(p_xyz(x,:,z)) )
       END DO
    END DO
    
    WRITE(*,*) "saving to ./output/Epxybx_x.dat"
    CALL writetocsvvec( "./output/Epxybx_x.dat", EGain )
    
    WRITE(*,*) "saving to ./output/Maxpxyzbxz_xz.dat"
    CALL writetocsv( "./output/Maxpxyzbxz_xz.dat", MaxGain )
    
    DEALLOCATE( EGain )
    DEALLOCATE( MaxGain )
  END SUBROUTINE EpxybxH0

  SUBROUTINE FirmShareNewMatch( U_x, H0p_xy, q_y_v_y, S_xy, s0, s1, wx, wy )
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: U_x(:)
    REAL(wp), INTENT(IN) :: H0p_xy(:,:)
    REAL(wp), INTENT(IN) :: q_y_v_y(:)
    REAL(wp), INTENT(IN) :: S_xy(:,:)
    REAL(wp), INTENT(IN) :: s0
    REAL(wp), INTENT(IN) :: s1
    REAL(wp), INTENT(IN) :: wx(:)
    REAL(wp), INTENT(IN) :: wy(:)
    INTEGER nx, ny, x, y, y2
    REAL(wp) share, hires_U, hires_E, share_E
    
    nx = SIZE(H0p_xy,1)
    ny = SIZE(H0p_xy,2)

    hires_U = 0.0
    hires_E = 0.0
    share_E = 0.0
    
    ! count hires from unemployment
    DO x=1,nx
       DO y=1,ny
          IF (S_xy(x,y) .GT. 0.0 ) THEN
             hires_U = hires_U + U_x(x) * s0 * q_y_v_y(y)
          END IF
       END DO
    END DO

    ! count hires from employment and accumulate surplus share
    DO x=1,nx
       DO y=1,ny
          DO y2=1,ny
             IF ( (S_xy(x,y) .GT. S_xy(x,y2)) &
                  .AND. (S_xy(x,y2) .GE. 0.0 ) ) THEN
                hires_E = hires_E + H0p_xy(x,y2) * s1 * q_y_v_y(y) * wy(y2)
                share_E = share_E + H0p_xy(x,y2) * s1 * q_y_v_y(y) * ( (S_xy(x,y) - S_xy(x,y2)) / S_xy(x,y) ) * wy(y2)
             END IF
          END DO
       END DO
    END DO
    ! the share of surplus to hte firm is all the surplus from hires out of unemployment and a share from j2j hires
    WRITE(*,*) "Hires from unemployment:         ", hires_U 
    WRITE(*,*) "Hires from other firms:          ", hires_E 
    WRITE(*,*) "Average Surplus from other firms:", share_E / hires_E 

    WRITE(*,*) "Average surplus share at new matches going to firm:", ( hires_U + share_E ) / ( hires_U + hires_E )
  END SUBROUTINE FirmShareNewMatch

  SUBROUTINE PmCmBxyz(p, b, PmCmB)
    ! PmB(x,y,z) = p(x,y,z) - b(x,z)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: p(:,:,:)
    REAL(wp), INTENT(IN) :: b(:,:)
    REAL(wp), INTENT(INOUT) :: PmCmB(:,:,:)
    INTEGER nx, ny, nz
    INTEGER x,y,z
    
    nx = SIZE(p,1)
    ny = SIZE(p,2)
    nz = SIZE(p,3)

    DO x=1,nx
       DO y=1,ny
          DO z=1,nz
             PmCmB(x,y,z) = p(x,y,z) - b(x,z) 
          END DO
       END DO
    END DO
  END SUBROUTINE PmCmBxyz

  ! calculate distribution of employed workers after delta and z shocks
  SUBROUTINE Hxytp(Sxy, H0xy, H0p_xy, delta)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: Sxy(:,:)
    REAL(wp), INTENT(IN) :: H0xy(:,:)
    REAL(wp), INTENT(INOUT) :: H0p_xy(:,:)
    REAL(wp), INTENT(IN) :: delta

    INTEGER nx, ny, x, y
    
    nx = SIZE(H0xy,1)
    ny = SIZE(H0xy,2)

    H0p_xy = 0.0

    DO x=1,nx
       DO y=1,ny
          IF ( Sxy(x,y) .GE. 0.0 ) THEN
             H0p_xy(x,y) = (1.0-delta) * H0xy(x,y) 
          END IF
       END DO
    END DO
  END SUBROUTINE Hxytp

  ! calculate distribution of unemployed workers after delta and z shocks
  SUBROUTINE Uxtp( ell_x, H0p_xy, wy, U_x )
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: ell_x(:)
    REAL(wp), INTENT(IN) :: H0p_xy(:,:)
    REAL(wp), INTENT(IN) :: wy(:)
    REAL(wp), INTENT(INOUT) :: U_x(:)
    
    INTEGER nx, x, ny, y
    REAL(wp) tmp
    
    nx = SIZE(H0p_xy, 1)
    ny = SIZE(H0p_xy, 2)
    
    DO x=1,nx
       tmp = 0.0
       DO y=1,ny
          tmp = tmp + H0p_xy(x,y)*wy(y)
       END DO
       U_x(x) = MAX( 0.0, ell_x(x) - tmp )
    END DO
  END SUBROUTINE Uxtp

  SUBROUTINE Lt( U_x, H0p_xy, s0, s1, wx, wy, L )
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: U_x(:)
    REAL(wp), INTENT(IN) :: H0p_xy(:,:)
    REAL(wp), INTENT(IN) :: s0
    REAL(wp), INTENT(IN) :: s1
    REAL(wp), INTENT(IN) :: wx(:)
    REAL(wp), INTENT(IN) :: wy(:)
    REAL(wp), INTENT(INOUT) :: L

    INTEGER nx, ny, x, y
    REAL(wp) tmp1, tmp2
    
    nx = SIZE(H0p_xy, 1)
    ny = SIZE(H0p_xy, 2)
    
    tmp1 = 0.0
    tmp2 = 0.0

    DO x=1,nx
       tmp1 = tmp1 + U_x(x)*wx(x)
    END DO

    DO x=1,nx
       DO y=1,ny
          tmp2 = tmp2 + H0p_xy(x,y)*wx(x)*wy(y)
       END DO
    END DO
    
    L = s0*tmp1 + s1*tmp2
  END SUBROUTINE Lt

  
  ! calculate the distribution of newly unemployed workers after delta and z 
  SUBROUTINE UNxtp( H_xy, Hp_xy, wy, UN_x )
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: H_xy(:,:)
    REAL(wp), INTENT(in) :: Hp_xy(:,:)
    REAL(wp), INTENT(in) :: wy(:)
    REAL(wp), INTENT(inout) :: UN_x(:)

    INTEGER nx, x, ny, y

    nx = SIZE(H_xy,1)
    ny = SIZE(H_xy,2)

    UN_x = 0.0

    DO x=1,nx
       DO y=1,ny
          UN_x(x) = UN_x(x) + ( H_xy(x,y) - Hp_xy(x,y) ) * wy(y)
       END DO
    END DO
  END SUBROUTINE UNxtp

  ! calculate the distribution of newly unemployed workers after delta and z 
  SUBROUTINE Hxtp( Hp_xy, wy, H_x )
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: Hp_xy(:,:)
    REAL(wp), INTENT(in) :: wy(:)
    REAL(wp), INTENT(inout) :: H_x(:)
    
    INTEGER nx, ny, x, y
    
    nx = SIZE(Hp_xy,1)
    ny = SIZE(Hp_xy,2)

    H_x = 0.0

    DO x=1,nx
       DO y=1,ny
          H_x(x) = H_x(x) +  Hp_xy(x,y) * wy(y)
       END DO
    END DO
  END SUBROUTINE Hxtp

  ! calculate the correlation between worker-firm types in matches
  SUBROUTINE sortxy( x_x, y_y, H_xy, corr_xy )
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: x_x(:)
    REAL(wp), INTENT(in) :: y_y(:)
    REAL(wp), INTENT(in) ::  H_xy(:,:)
    REAL(wp), INTENT(inout) :: corr_xy
    REAL(wp) mean_x, mean_y, sd_x, sd_y, H
    INTEGER x,y,nx,ny

    nx = SIZE(H_xy,1)
    ny = SIZE(H_xy,2)

    mean_x = 0.0
    mean_y = 0.0
    sd_x   = 0.0
    sd_y   = 0.0
    H = SUM(H_xy)
    corr_xy = 0.0

    DO x=1,nx
       DO y=1,ny
          mean_x = mean_x + H_xy(x,y) * x_x(x)
          mean_y = mean_y +  H_xy(x,y) * y_y(y)
       END DO
    END DO
    mean_x = mean_x / H
    mean_y = mean_y / H
    
    DO x=1,nx
       DO y=1,ny
          sd_x = sd_x + H_xy(x,y) * (( x_x(x) - mean_x )**2)
          sd_y = sd_y + H_xy(x,y) * (( y_y(y) - mean_y )**2) 
          corr_xy = corr_xy + H_xy(x,y)*(x_x(x)-mean_x)*(y_y(y)-mean_y)
       END DO
    END DO
    sd_x = SQRT( sd_x / H )
    sd_y = SQRT( sd_y / H )
    corr_xy = corr_xy / H  ! this is the covariance at this point

    corr_xy = corr_xy / ( sd_x * sd_y )
  END SUBROUTINE SORTXY


  ! calculate output 
  SUBROUTINE GDP( H_xy, P_xy, wx, wy, gdp_pc )
    IMPLICIT NONE
    REAL(wp), INTENT(in) ::  H_xy(:,:)
    REAL(wp), INTENT(in) ::  P_xy(:,:)
    REAL(wp), INTENT(in) ::  wx(:)
    REAL(wp), INTENT(in) ::  wy(:)
    
    REAL(wp) gdp_pc
    
    INTEGER x,y,nx,ny

    gdp_pc = 0.0

    nx = SIZE(H_xy,1)
    ny = SIZE(H_xy,2)
    
    DO x=1,nx
       DO y=1,ny
          gdp_pc = gdp_pc + H_xy(x,y) * &
               ( P_xy(x,y) ) *wx(x)*wy(y)
       END DO
    END DO
    gdp_pc = gdp_pc ! This version uses level of gdp, not productivity
  END SUBROUTINE GDP

  ! calculate output per y-type firm 
  SUBROUTINE GDP_y( H_xy, P_xy, wx, gdpy )
    IMPLICIT NONE
    REAL(wp), INTENT(in) ::  H_xy(:,:)
    REAL(wp), INTENT(in) ::  P_xy(:,:)
    REAL(wp), INTENT(in) ::  wx(:)
    REAL(wp), intent(out) ::  gdpy(:)
    
    INTEGER x,y,nx,ny

    gdpy = 0.0

    nx = SIZE(H_xy,1)
    ny = SIZE(H_xy,2)
    
    DO y=1,ny
       DO x=1,nx
          gdpy(y) = gdpy(y) + H_xy(x,y)*P_xy(x,y)*wx(x)
       END DO
    END DO
  END SUBROUTINE GDP_y

  ! calculate employment per y-type firm 
  SUBROUTINE E_y( H_xy, wx, Ey )
    IMPLICIT NONE
    REAL(wp), INTENT(in) ::  H_xy(:,:)
    REAL(wp), INTENT(in) ::  wx(:)
    REAL(wp), intent(out) ::  Ey(:)
    
    INTEGER x,y,nx,ny

    Ey = 0.0

    nx = SIZE(H_xy,1)
    ny = SIZE(H_xy,2)
    
    DO y=1,ny
       DO x=1,nx
          Ey(y) = Ey(y) + H_xy(x,y)*wx(x)
       END DO
    END DO
  END SUBROUTINE E_y



 
  ! calculate standard deviation of value added per worker by y-type firms 
  ! (large firm assumption)
  SUBROUTINE sdp( H_xy, P_xy, wx, wy, epsilon, sd )
    IMPLICIT NONE
    REAL(wp), INTENT(in) ::    H_xy(:,:)
    REAL(wp), INTENT(in) ::    P_xy(:,:)
    REAL(wp), INTENT(in) ::    wx(:)
    REAL(wp), INTENT(in) ::    wy(:)
    REAL(wp), INTENT(in) ::    epsilon
    REAL(wp), INTENT(inout) :: sd
    
    INTEGER x,y,nx,ny
    REAL(wp)  mean, H

    mean = 0.0
    sd = 0.0
    H = quad_xy( H_xy, wx, wy )

    nx = SIZE(H_xy,1)
    ny = SIZE(H_xy,2)
    
    DO x=1,nx
       DO y=1,ny
          IF ( H_xy(x,y) .GT. 0.0 ) THEN
             mean = mean + H_xy(x,y) &
                  * LOG( max(epsilon, P_xy(x,y)) )*wx(x)*wy(y)
          END IF
       END DO
    END DO
    mean = mean / H

    DO x=1,nx
       DO y=1,ny
          IF ( H_xy(x,y) .GT. 0.0 ) THEN
             sd = sd + H_xy(x,y)* &
                  ( ( LOG( max(epsilon, P_xy(x,y) ) ) - mean )**2 ) &
                  *wx(x)*wy(y)
          END IF
       END DO
    END DO
    sd = sd / H
    sd = SQRT( sd )
  END SUBROUTINE SDP

  SUBROUTINE sdlogx( xx, H_xy, wx, wy, epsilon, sdlogxx )
    IMPLICIT NONE
    REAL(wp), INTENT(in) ::    xx(:)
    REAL(wp), INTENT(in) ::    H_xy(:,:)
    REAL(wp), INTENT(in) ::    wx(:)
    REAL(wp), INTENT(in) ::    wy(:)
    REAL(wp), INTENT(in) ::    epsilon
    REAL(wp), INTENT(inout) :: sdlogxx

    INTEGER x,y,nx,ny
    REAL(wp) mean_xx, H 

    mean_xx  = 0.0
    sdlogxx = 0.0

    nx = SIZE(H_xy,1)
    ny = SIZE(H_xy,2)

    H = quad_xy(H_xy, wx, wy)

    DO x=1,nx
       DO y=1,ny
          mean_xx = mean_xx + LOG(xx(x)+epsilon) * H_xy(x,y)*wx(x)*wy(y)
       END DO
    END DO
    mean_xx = mean_xx / H
    
    DO x=1,nx
       DO y=1,ny
          sdlogxx = sdlogxx + ((LOG(xx(x)+epsilon)-mean_xx)**2) &
               * H_xy(x,y)*wx(x)*wy(y)
       END DO
    END DO
    
    sdlogxx = sdlogxx / H
    sdlogxx = SQRT( sdlogxx )
  END SUBROUTINE sdlogx

   
  SUBROUTINE J1y( U_x, S0_xy, wx, Jtmp1 )
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: U_x(:)
    REAL(wp), INTENT(IN) :: S0_xy(:,:)
    REAL(wp), INTENT(IN) :: wx(:)
    REAL(wp), INTENT(INOUT) :: Jtmp1(:)
    
    INTEGER x, y, nx, ny
    
    nx = SIZE(S0_xy,1)
    ny = SIZE(S0_xy,2)

    Jtmp1 = 0.0
    DO y=1,ny
       DO x=1,nx
          Jtmp1(y) = Jtmp1(y) + MAX(0.0, S0_xy(x,y))*U_x(x)*wx(x)
       END DO
    END DO
  END SUBROUTINE J1y

  SUBROUTINE J2y( H0p_xy, S0_xy, wx, wy, Jtmp2 )
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: H0p_xy(:,:)
    REAL(wp), INTENT(IN) :: S0_xy(:,:)
    REAL(wp), INTENT(IN) :: wx(:)
    REAL(wp), INTENT(IN) :: wy(:)
    REAL(wp), INTENT(INOUT) :: Jtmp2(:)
    
    INTEGER x, y, y2, nx, ny
    
    nx = SIZE(H0p_xy,1)
    ny = SIZE(H0p_xy,2)

    Jtmp2 = 0.0

    DO y=1,ny
       DO x=1,nx
          DO y2=1,ny
             Jtmp2(y) = Jtmp2(y) + &
                  MAX( 0.0, (S0_xy(x,y) - S0_xy(x,y2)) ) &
                  * H0p_xy(x,y2)*wx(x)*wy(y2)
          END DO
       END DO
    END DO
  END SUBROUTINE J2y
  
  SUBROUTINE IntJ( J_y, eta, gamma0, gamma1, y_y, wy, Jc )
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: J_y(:)
    REAL(wp), INTENT(IN) :: eta
    REAL(wp), INTENT(IN) :: gamma0
    REAL(wp), INTENT(IN) :: gamma1
    REAL(wp), INTENT(IN) :: y_y(:)
    REAL(wp), INTENT(IN) :: wy(:)
    REAL(wp), INTENT(INOUT) :: Jc
    
    INTEGER y, ny
    REAL(wp) etainv
    
    ny = SIZE(J_y,1)

    etainv = 1.0/eta
    
    Jc = 0.0

	if (GLOBAL_HET_COSTS) then
	    DO y=1,ny
    	   Jc = Jc + wy(y) * ( ( J_y(y) / (gamma0*(y_y(y)**gamma1) ) ) ** etainv )
    	END DO
    else
	    DO y=1,ny
    	   Jc = Jc + wy(y) * ( ( J_y(y) / (gamma0 ) ) ** etainv )
    	END DO
	end if    
    
  END SUBROUTINE IntJ

  SUBROUTINE computeTheta(Jc, L, MT, alpha, omega, eta, theta )
    IMPLICIT NONE
    REAL(wp), INTENT(IN)    :: Jc
    REAL(wp), INTENT(IN)    :: L
    INTEGER,  INTENT(IN)    :: MT       ! type of matching technology
    REAL(wp), INTENT(IN)    :: alpha
    REAL(wp), INTENT(IN)    :: omega
    REAL(wp), INTENT(IN)    :: eta
    REAL(wp), INTENT(INOUT) :: theta

    IF ( MT .EQ. 1 ) THEN
       ! Cobb-Douglas Matching
       theta = ( alpha**(1.0/(eta+omega)) ) * ( (Jc/L)**(eta/(eta+omega)) )
    ELSE IF ( MT .EQ. 2 ) THEN
       ! Stevens CES Matching
       theta = 0.0 
       ! NOTE: I STILL NEED TO CODE THIS.  IT WILL REQUIRE A NUMERICAL SOLUTION UNLESS eta=1
    ELSE
       ! DenHaan CES Matching
       theta = 0.0
       ! NOTE: I STILL NEED TO CODE THIS.  IT WILL REQUIRE A NUMERICAL SOLUTION
    END IF
  END SUBROUTINE computeTheta
  
  SUBROUTINE computeMatching(L, V, MT, alpha, omega, M)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: L
    REAL(wp), INTENT(IN) :: V
    INTEGER, INTENT(IN) :: MT       ! type of matching technology
    REAL(wp), INTENT(IN) :: alpha
    REAL(wp), INTENT(IN) :: omega
    REAL(wp), INTENT(INOUT) :: M

    IF ( MT .EQ. 1 ) THEN
       ! Cobb-Douglas Matching
       M = alpha * (L**omega) * (V**(1.0-omega))
       M = MIN( M, L, V )
    ELSE IF ( MT .EQ. 2 ) THEN
       ! Stevens CES Matching
       M = alpha * omega * L * V / ( omega * L + V )
    ELSE 
       ! DenHaan CES Matching
       M = alpha*L*V / ( (L**omega + V**omega)**(1.0/omega) )
    END IF
  END SUBROUTINE computeMatching
    
  SUBROUTINE vy(J_y, M, V, gamma0, gamma1, eta, y_y, v_y)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: J_y(:)
    REAL(wp), INTENT(IN) :: M
    REAL(wp), INTENT(IN) :: V
    REAL(wp), INTENT(IN) :: eta
    REAL(wp), INTENT(IN) :: gamma0
    REAL(wp), INTENT(IN) :: gamma1
    REAL(wp), INTENT(IN) :: y_y(:)
    REAL(wp), INTENT(INOUT) :: v_y(:)

    REAL(wp) inveta
    
    INTEGER y, ny
    ny = SIZE(v_y,1)

    inveta = 1.0/eta
    v_y = 0.0
    if (GLOBAL_HET_COSTS) then
	    DO y=1,ny
    	   v_y(y) = ( (M/V)  * ( J_y(y) / ( gamma0*(y_y(y)**gamma1) ) ) ) ** inveta
    	END DO
    else
	    DO y=1,ny
    	   v_y(y) = ( (M/V)  * ( J_y(y) / gamma0 ) ) ** inveta
    	END DO
	end if    
  END SUBROUTINE vy

  SUBROUTINE H1xyt( S_xy, U_x, q_y_v_y, Htmp1_xy )
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: S_xy(:,:)
    REAL(wp), INTENT(IN) :: U_x(:)
    REAL(wp), INTENT(IN) :: q_y_v_y(:)
    REAL(wp), INTENT(INOUT) :: Htmp1_xy(:,:)
    
    INTEGER x, y, nx, ny

    nx = SIZE(Htmp1_xy,1)
    ny = SIZE(Htmp1_xy,2)
    
    Htmp1_xy = 0.0
    
    DO x=1,nx
       DO y=1,ny
          IF (S_xy(x,y) .GT. 0.0 ) THEN
             Htmp1_xy(x,y) = U_x(x) * q_y_v_y(y)
          END IF
       END DO
    END DO
  END SUBROUTINE H1xyt

  SUBROUTINE H2xyt( S_xy, H0p_xy, q_y_v_y, wy, Htmp2_xy )
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: S_xy(:,:)
    REAL(wp), INTENT(IN) :: H0p_xy(:,:)
    REAL(wp), INTENT(IN) :: q_y_v_y(:)
    REAL(wp), INTENT(IN) :: wy(:)
    REAL(wp), INTENT(INOUT) :: Htmp2_xy(:,:)
    
    INTEGER x, y, y2, nx, ny
    
    nx = SIZE(H0p_xy,1)
    ny = SIZE(H0p_xy,2)
    
    Htmp2_xy = 0.0
    
    DO x=1,nx
       DO y=1,ny
          DO y2=1,ny
             IF ( (S_xy(x,y) .GT. S_xy(x,y2)) &
                  .AND. (S_xy(x,y2) .GE. 0.0 ) ) THEN
                Htmp2_xy(x,y) = Htmp2_xy(x,y) &
                     + H0p_xy(x,y2) * q_y_v_y(y) * wy(y2)
             END IF
          END DO
       END DO
    END DO
  END SUBROUTINE H2xyt

  SUBROUTINE H3xyt( S_xy, H0p_xy, q_y_v_y, wy, Htmp3_xy )
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: S_xy(:,:)
    REAL(wp), INTENT(IN) :: H0p_xy(:,:)
    REAL(wp), INTENT(IN) :: q_y_v_y(:)
    REAL(wp), INTENT(IN) :: wy(:)
    REAL(wp), INTENT(INOUT) :: Htmp3_xy(:,:)
    
    INTEGER x, y, y2, nx, ny
    
    nx = SIZE(H0p_xy,1)
    ny = SIZE(H0p_xy,2)
    
    Htmp3_xy = 0.0
    
    DO x=1,nx
       DO y=1,ny
          DO y2=1,ny
             IF ( (S_xy(x,y2) .GT. S_xy(x,y)) &
                  .AND. (S_xy(x,y) .GE. 0.0 ) ) THEN
                Htmp3_xy(x,y) = Htmp3_xy(x,y) &
                     + H0p_xy(x,y) * q_y_v_y(y2) * wy(y2)
             END IF
          END DO
       END DO
    END DO
  END SUBROUTINE H3xyt

  ! compute numerical integration in one dimension
  FUNCTION quad_x( F, wx  )
    IMPLICIT NONE
    REAL(wp) quad_x
    REAL(wp), INTENT(IN) :: F(:)
    REAL(wp), INTENT(IN) :: wx(:)

    INTEGER x, nx

    nx = SIZE(F, 1)

    quad_x = 0.0

    DO x=1,nx
       quad_x = quad_x + F(x) * wx(x)
    END DO
  END FUNCTION quad_x

  ! compute numerical integration in 2 dimensions
  FUNCTION quad_xy( F, wx, wy )
    IMPLICIT NONE
    REAL(wp) quad_xy
    REAL(wp), INTENT(IN) :: F(:,:)
    REAL(wp), INTENT(IN) :: wx(:)
    REAL(wp), INTENT(IN) :: wy(:)

    INTEGER x, y, nx, ny 
    
    nx = SIZE(F, 1)
    ny = SIZE(F, 2)

    quad_xy = 0.0
    
    DO x=1,nx
       DO y=1,ny
          quad_xy = quad_xy + F(x,y) * wx(x) * wy(y)
       END DO
    END DO
  END FUNCTION quad_xy

END MODULE solver
