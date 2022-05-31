!------------------------------------------------------------------------------
FUNCTION getSourceField( model, n ) RESULT(Hx)
!------------------------------------------------------------------------------
! Calculate magnetic field strength, Hs, from a line current source
! Reference: R.Gupta, Field Calculations and Computations
!
! Author: Ye Yang
! Date  : 11.11.2021
!
!------------------------------------------------------------------------------
    USE DefUtils
    IMPLICIT NONE

    TYPE(Model_t) :: model
    INTEGER :: n
!------------------------------------------------------------------------------
    INTEGER :: i
    INTEGER, PARAMETER :: nwire = 2140, nmax = 5 

    REAL(KIND=dp) :: x, y
    REAL(KIND=dp) :: Hx, Hy, Hr, Ht, a, r
    REAL(KIND=dp) :: tt, pp
    REAL(KIND=dp) :: xs(nwire), ys(nwire)
    REAL(KIND=dp), PARAMETER :: curr = 264.9

    LOGICAL :: FirstVisit = .TRUE.
    LOGICAL :: Found
    
    SAVE FirstVisit, xs, ys
!------------------------------------------------------------------------------

    ! do something only the first time the function is visited
    IF (FirstVisit) THEN
      OPEN(1,file='CoilWirePosition.txt',status='old',action='read')

      PRINT *, 'READING WIRE LOCATION.'
      PRINT *, 'X[m]', 'Y[m]'

      DO i = 1, nwire
        READ(1,*) xs(i), ys(i)
        !PRINT *, i, xs(i), ys(i)
      END DO

      PRINT *, 'FININSHED LOADING WIRE LOCATION.'
      FirstVisit = .FALSE.
    END IF

    ! get the coordinates of the point
    x = model % Nodes % x(n)
    y = model % Nodes % y(n)

    PRINT *, SIZE(xs)
    PRINT *, x, y, xs(3), ys(3)

    ! loop for the calculation of magnetic field
    Hr = 0.0d0
    Ht = 0.0d0

    DO i = 1, nwire
      r = SQRT(  x**2 +  y**2 )
      a = SQRT( xs(i)**2 + ys(i)**2 )
      tt= DATAN2( y, x)
      pp= DATAN2(ys(i),xs(i))

      IF (r.LT.a) THEN
        DO n = 0, nmax
          Hr = Hr + (r/a)**(n-1) * SIN(n*(pp-tt))
          Ht = Ht + (r/a)**(n-1) * COS(n*(pp-tt))
        END DO
      END IF

      IF (r.GT.a) THEN
        DO n = 0, nmax
          Hr = Hr + (a/r)**(n+1) * SIN(n*(pp-tt))
          Ht = Ht + (a/r)**(n+1) * COS(n*(pp-tt))
        END DO
      END IF

    END DO

    Hr = Hr * curr / (2.0*PI*a)
    Ht = Ht * curr / (2.0*PI*a)

    ! coordinate transformation for magnetic field
    Hx = Hr*COS(tt) - Ht*SIN(tt)
    Hy = Hr*SIN(tt) + Ht*COS(tt) 

!------------------------------------------------------------------------------
END FUNCTION getSourceField
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!SUBROUTINE GetFieldFromLineCurrent(xs, ys, x, y, curr, nmax, Hx, Hy )
!------------------------------------------------------------------------------
! Calculate magnetic field strength, Hs, from a line current source
! Reference: R.Gupta, Field Calculations and Computations
!
! xs, ys: source point of line current
! x , y : observe point
! curr  : operation current
! nmax  : maximum order of field calculation related to the precision
! Hx, Hy: magnetic field strength
!
! Author: Ye Yang
! Date  : 11.11.2021
!------------------------------------------------------------------------------
!    IMPLICIT NONE
!
!    REAL(KIND=dp) :: x, y, xs, ys, Hx, Hy, Hr, Ht, a, r, theta, phi, curr
!    INTEGER :: nmax, n
!------------------------------------------------------------------------------
!    Hr = 0.0d0
!    Ht = 0.0d0
!    Hx = 0.0d0
!    Hy = 0.0d0

    ! coordinate transformation
!    r    = SQRT(  x**2 +  y**2 )
!    a    = SQRT( xs**2 + ys**2 )
!    theta= DATAN2( y, x)
!    phi  = DATAN2(ys,xs)

!    IF (r.EQ.a) THEN
!      EXIT
!    END IF

!    IF (r.LT.a) THEN
!      DO n = 0, nmax
!        Hr += (r/a)**(n-1) * SIN(n*(phi-theta))
!        Ht += (r/a)**(n-1) * COS(n*(phi-theta))
!      END DO
!      Hr *= curr / (2.0*PI*a)
!      Ht *= curr / (2.0*PI*a)
!    END IF

!    IF (r.GT.a) THEN
!      DO n = 0, nmax
!        Hr += (a/r)**(n+1) * SIN(n*(phi-theta))
!        Ht += (a/r)**(n+1) * COS(n*(phi-theta))
!      END DO
!      Hr *= curr / (2.0*PI*a)
!      Ht *= curr / (2.0*PI*a)
!    END IF

    ! coordinate transformation for magnetic field
!    Hx = Hr*COS(theta) - Ht*SIN(theta)
!    Hy = Hr*SIN(theta) + Ht*COS(theta) 

!------------------------------------------------------------------------------
!END SUBROUTINE GetFieldFromLineCurrent
!------------------------------------------------------------------------------




