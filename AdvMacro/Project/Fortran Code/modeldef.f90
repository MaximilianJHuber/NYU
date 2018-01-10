MODULE modeldef
  USE params
  USE stat_helper
  USE glob
  USE output
  
  IMPLICIT NONE
  
  !This represents an equilibrium solution
  TYPE Model
     REAL(wp), ALLOCATABLE :: S_xyz(:,:,:)     ! Surplus function S(x,y,z)
     REAL(wp), ALLOCATABLE :: BB_xz(:,:)       ! value of unemployment
     REAL(wp), ALLOCATABLE :: Q_zz(:,:)        ! Markov transition matrix
     REAL(wp), ALLOCATABLE :: Q_zz_CDF(:,:)    ! CDF for simulation
     REAL(wp), ALLOCATABLE :: p_xyz(:,:,:)     ! Market production
     REAL(wp), ALLOCATABLE :: b_xz(:,:)        ! Home production
     REAL(wp), ALLOCATABLE :: H0_xy(:,:)       ! dist in non-stochastic SS
     REAL(wp), ALLOCATABLE :: ell_x(:)         ! distn of workers
  END TYPE Model
  
  INTERFACE Init
     MODULE PROCEDURE InitModel
  END INTERFACE
  
CONTAINS
  
  SUBROUTINE InitModel(m,p)
    IMPLICIT NONE
    TYPE (Model), INTENT(inout) :: m
    TYPE (ExogenousParameters), INTENT(in) :: p

    INTEGER i,j,x,y,z
    REAL(wp) tmp, tmp_p_xy
    
    ALLOCATE( m%S_xyz ( p%nx, p%ny, p%nz) )
    ALLOCATE( m%BB_xz ( p%nx, p%nz ) )
    ALLOCATE( m%Q_zz ( p%nz, p%nz) )
    ALLOCATE( m%Q_zz_CDF ( p%nz, p%nz) )
    ALLOCATE( m%p_xyz( p%nx, p%ny, p%nz) )
    ALLOCATE( m%b_xz( p%nx, p%nz ) )
    ALLOCATE( m%H0_xy( p%nx, p%ny ) )
    ALLOCATE( m%ell_x( p%nx ) )


    ! Initialize objects
    
    ! use Beta distributed worker types
    m%ell_x = ( p%x ** (p%shape1 - 1.0) ) * & 
         ( (1.0-p%x ) ** (p%shape2 - 1.0) )
    m%ell_x = real(p%nx, wp) *  m%ell_x / SUM(m%ell_x)
    
    ! Production Function (value added)
    DO x=1,p%nx
       DO y=1,p%ny
          DO z=1,p%nz
             ! NOTE: we have merged the operating costs into p_xyz which is now directly
             !       value added
             ! 
             m%p_xyz(x,y,z) = p%f0 * EXP(p%sigma*p%z(z)) & 
                  * ( p%f1 + p%f2*p%x(x) + p%f3*p%y(y) &
                  + p%f4*(p%x(x)**2) + p%f5*(p%y(y)**2) &
                  + p%f6 * p%x(x) * p%y(y) ) * p%dt
          END DO
       END DO
    END DO
        
    ! Home production
    DO x=1,p%nx
       tmp_p_xy = maxval( m%p_xyz(x,:,(p%nz+1)/2) )
       m%b_xz(x,:) = 0.7 * tmp_p_xy
    END DO

    ! Initialize S_xyz to non-stochastic steady state value
    DO x=1,p%nx
       DO y=1,p%ny
          DO z=1,p%nz
             m%S_xyz(x,y,z)=&
                  ((1.0+p%r)/(p%r+p%delta)) * &
                  (m%p_xyz(x,y,z)-m%b_xz(x,z))
          END DO
       END DO
    END DO

    ! Initialize H0_xyij as workers x-i unifomly distributed over jobs
    DO x=1,p%nx
       m%H0_xy(x,:) = m%ell_x(x)
    END DO
    
    ! Create Markov transition matrix Q_zz
    DO i=1,p%nz
       DO j=1,p%nz
          ! Gaussian Copula:
          m%Q_zz(i,j) =  EXP(-(1.0/(2.0*(1.0-p%rho**2))) &
               *(p%z(i)**2+p%z(j)**2-2.0*p%rho*p%z(i)*p%z(j)))
!!$          ! t-copula:
!!$          m%Q_zz(i,j) = (1.0 + (p%z(i)**2+p%z(j)**2-2.0*p%rho*p%z(i)*p%z(j)) &
!!$               / (p%nu*(1.0-p%rho**2)) )**(-(p%nu+2.0)/2.0)
       END DO
    END DO
    DO i=1,p%nz
       tmp = sum( m%Q_zz(i,:) )
       DO j=1,p%nz
          m%Q_zz(i,j) = m%Q_zz(i,j) / tmp
       END DO
    END DO

   ! Create conditional CDF of Markov Transiton matrix for simulation
    
    DO i=1,p%nz
       m%Q_zz_CDF(i,1) = m%Q_zz(i,1)
       DO j=2,p%nz
          m%Q_zz_CDF(i,j) = m%Q_zz_CDF(i,j-1) + m%Q_zz(i,j)
       END DO
    END DO
    ! normalize to make sure upper bound is 1 numerically
    DO i=1,p%nz
       DO j=1,p%nz
          m%Q_zz_CDF(i,:) = m%Q_zz_CDF(i,:) /  m%Q_zz_CDF(i,p%nz)
       END DO
    END DO
    
    
  END SUBROUTINE InitModel

  SUBROUTINE FreeModel(m)
    IMPLICIT NONE
    TYPE (Model), INTENT(inout) :: m
    
    DEALLOCATE( m%S_xyz )
    DEALLOCATE( m%BB_xz )
    DEALLOCATE( m%Q_zz )
    DEALLOCATE( m%Q_zz_CDF )
    DEALLOCATE( m%p_xyz )
    DEALLOCATE( m%b_xz )
    DEALLOCATE( m%H0_xy )
    DEALLOCATE( m%ell_x )
    
  END SUBROUTINE FreeModel

  !! Put any other needed functions in here ...


END MODULE modeldef
