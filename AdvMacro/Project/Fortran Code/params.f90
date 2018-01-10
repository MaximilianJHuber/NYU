MODULE params
  ! This module defines the parameters that will be used in solving the model
  
  USE glob
  USE array_helper
  use stat_helper
  
  IMPLICIT NONE
  
  PUBLIC
  TYPE ExogenousParameters
     ! Calibrated parameters for Low Skilled workers

     ! grid points (not estimated)
     INTEGER nx           ! grid points for human capital
     INTEGER ny           ! grid points for firm productivity   
     INTEGER nz           ! grid points for aggregate shock
 
     !! Distribution of worker types, Beta distributed
     !NOTE: these have dimenstion equal to the number of education groups
     REAL(wp) shape1
     REAL(wp) shape2

     REAL(wp) gamma0       ! scale of search costs
     REAL(wp) gamma1       ! power term in scale of search costs (function of firm type)
     REAL(wp) eta          ! elasticity of vacancy costs with respect to 
                           ! vacancy creation (minus 1)
     
     ! paramaters for simulation (not estimated)
     REAL(wp) dt           ! length of a period in years

     REAL(wp) epsilon      ! distance to stay away from 0 and 1 for z shocks

     ! discounting and separation
     REAL(wp) r                ! interest rate  
     REAL(wp) delta            ! separation rate
     
     ! home production
     REAL(wp) b0               ! fixed amount of home production
     REAL(wp) b1               ! variable amount of home production (quadratic)
     REAL(wp) b2               ! variable amount of home production (linear)
     REAL(wp) b3               
     REAL(wp) b4               
     REAL(wp) b5               

     ! value added
     REAL(wp) f0               ! scale of all flow values  
     REAL(wp) f1               ! loading on aggregate shock
     REAL(wp) f2               ! loading on xy component
     REAL(wp) f3               ! loading on y^2
     REAL(wp) f4               ! loading on x^2
     REAL(wp) f5               ! loading on y
     REAL(wp) f6               ! loading on x
    
     ! Stochastic process for aggregate shock
     REAL(wp) rho              ! parameter for gaussian Copula
     REAL(wp) sigma            ! standard deviation of the  distributin for 
     REAL(wp) nu               ! degrees of freedom when using the t-copula

     ! Meeting function 
     INTEGER  MatchingTechnology    ! 1: Cobb-Douglas, 2: Stevens CES, 
                                    ! 3: DenHaan CES  
     REAL(wp) alpha        ! in [0,1]  ! Scales the number of matches
     REAL(wp) omega        ! in [0,1]  
     ! 1: alpha*L^omega*V^(1-omega), 2: alpha*omega*L*V / (omega*L + V) 
     ! 3: alpha*L*V/ ((L^omega + V^omega)^(1/omega)) 

     REAL(wp) s0         ! unemployed search intensity
     REAL(wp) s1         ! employed search intensity
     
     !! create vectors for x,y,z
     REAL(wp), ALLOCATABLE ::  x(:)
     REAL(wp), ALLOCATABLE ::  y(:)
     REAL(wp), ALLOCATABLE ::  z(:)
     REAL(wp), ALLOCATABLE ::  z_grid(:)  ! i don't know if this is needed
     REAL(wp), ALLOCATABLE ::  wx(:) ! weights for numerical integration
     REAL(wp), ALLOCATABLE ::  wy(:) ! weights for numerical integration

     !! For simulation NOT ESTIMATED
     ! for simulation (not estimated)
     INTEGER bigT      ! number of weeks to simulate, including burn-in
     INTEGER sample    ! number of quarters to base estimation on
     INTEGER resample  ! number of replications to average out 
                       ! simulation error in calculation of moments

     INTEGER quarters  ! number of quarters for estimation (excuding burnin)


     ! numerical methods constants  (not estimated)
     INTEGER MaxIter   ! maximum number of iterations
     REAL(wp) tol      ! numerical tolerance for convergence
     REAL(wp) tol_H    ! numerical tolerance for convergence

     ! HP filtering parameter
     REAL(wp) hp_lambda 

  END TYPE ExogenousParameters
  
  INTERFACE Init
     MODULE PROCEDURE InitExogenousParameters
  END INTERFACE
CONTAINS
  SUBROUTINE InitExogenousParameters(p)
    IMPLICIT NONE

    integer i

    TYPE(ExogenousParameters), INTENT(INOUT) :: p

    ! grid points for ESTIMATION
    p%nx        =  21            ! grid points for human capital
    p%ny        =  21            ! grid points for firm productivity   
    p%nz        =  51            ! grid points for aggregate shock

!!$    ! grid points for creating SMOOTHED figures 2, 3 and 4
!!$    p%nx        =  151            ! grid points for human capital
!!$    p%ny        =  151            ! grid points for firm productivity   
!!$    p%nz        =  151            ! grid points for aggregate shock

    p%shape1    = 1.0           ! shape of worker distribution 1
    p%shape2    = 1.0           ! shpae of worker distribution 2

    p%gamma0    = 1.00
    p%gamma1    = 0.00
    p%eta       = 1.00          ! elasticity of vacancy costs with respect to vacancy creation (minus 1)
    
    ! paramaters for simulation (not estimated)
    p%dt        = 1.0/52.0         ! length of a period in years
    
    
    p%epsilon   = 0.001         ! distance to stay away from 0 and 1 on the unit interval.
    
    ! discounting and separation
    p%r         = (1.05) ** p%dt - 1.0  ! interest rate  (not estimated)
    p%delta     = (1.40) ** p%dt - 1.0  ! separation rate
    
    ! home production
    p%b0        = 0.40            ! fixed amount of home production
    p%b1        = 0.00            ! variable amount of home production indexed to aggregate productivity
    p%b2        = 0.00
    p%b3        = 0.00
    p%b4        = 0.00
    p%b5        = 0.00

    ! value added
    p%f0     =    (1.0 / 0.398401314845053) / 0.412335720378949
    ! scale of GDP (just divide by realized mean VA to get to 1)
    ! I rescaled this to match the mean for the data through 2012
    p%f1     =    1.0             ! loading on aggregate shock
    p%f2     =    3.0             ! loading on xy part of production
    p%f3     =    1.0             ! loading on y^2
    p%f4     =    0.0             ! loading on x^2 
    p%f5     =    0.0             ! loading on y
    p%f6     =    0.0             ! loading on x 

    ! Stochastic process for aggregate shock
    p%rho       = 0.996                ! parameter for gaussian Copula
    ! p%rho       = 0.999                ! parameter for gaussian Copula
    p%sigma     = 0.067                ! standard deviation of the G distributin for aggregate shocks
    p%nu        = 3.5                  ! estiamted from quarterly gdp series as a good starting value

    ! Meeting function 
    ! NOTE: ONLY COBB-DOUGLAS IS CODES FOR THE TIME BEING.  RETURN TO THIS
    p%MatchingTechnology = 1     ! 1: Cobb-Douglas, 2: Stevens CES, 3: DenHaan CES  (not estimated)
    p%alpha  = 1.0       ! in [0,1]  ! Scales the number of matches
    p%omega  = 0.5       ! in [0,1]  ! 1: alpha*L^omega*V^(1-omega), 2: alpha*omega*L*V / (omega*L + V) 3: alpha*L*V/ ((L^omega + V^omega)^(1/omega)) 
    
    p%s0  = 1.00        ! unemployed search intensity
    p%s1  = 0.10        ! employed search intensity
    
    
    !! create vectors of x,y,z
    ALLOCATE ( p%x (p%nx) )
    ALLOCATE ( p%y (p%ny) )
    ALLOCATE ( p%z (p%nz) )
    ALLOCATE ( p%z_grid (p%nz) )
    ALLOCATE ( p%wx (p%nx) )
    ALLOCATE ( p%wy (p%ny) )

    CALL linspace( p%x,  p%epsilon, 1.0-p%epsilon ) 
    CALL linspace( p%y,  p%epsilon, 1.0-p%epsilon ) 
    CALL linspace( p%z_grid, p%epsilon, 1.0-p%epsilon ) 

    ! For now, just use equally spaced discrete grid
    p%wx = 1.0 / real( p%nx, wp )
    p%wy = 1.0 / real( p%ny, wp )

    ! used to find position in array 
    ! when drawing random uniform numbers
    do i=1,p%nz
       p%z(i) = norminv( p%z_grid(i) )
    end do
    
    !!! for testing, use 60 years !!!
    ! p%quarters = 60  ! 30 years for testing
    p%quarters = 2400  ! 600 years for estimation
    
    p%sample =  (p%quarters/4) / p%dt        ! base moments on ## year simulaiton

    !!! for testing use 10 years of burnin
    ! p%bigT = int(10 / p%dt) + p%sample ! 100 years of burn in
    p%bigT = int(100 / p%dt) + p%sample ! 100 years of burn in

    p%resample = 1 ! number of replications to average out simulation error in calculation of moments
        
    ! numerical methods constants  (not estimated)
    p%MaxIter = 2000
    p%tol     = 10e-8
    p%tol_H   = 10e-5   ! Note this is only the tolerance for teh initial condition of the 
                        ! simulation.
    
    ! HP filtering parameter
    p%hp_lambda = 100000.0
    !p%hp_lambda = 1600.0

END SUBROUTINE InitExogenousParameters


SUBROUTINE FreeExogenousParameters(p)
  IMPLICIT NONE
  TYPE(ExogenousParameters), INTENT(INOUT) :: p
  
  DEALLOCATE ( p%x  )
  DEALLOCATE ( p%y  )
  DEALLOCATE ( p%z  )
  DEALLOCATE ( p%z_grid )
  DEALLOCATE ( p%wx )
  DEALLOCATE ( p%wy )
  
END SUBROUTINE FreeExogenousParameters

end MODULE params

