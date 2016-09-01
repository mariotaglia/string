!***********************************************************************************************************************************
!
!                                                        L I N T E R P O L
!
!  Subroutine:   LINTERPOL
!
!  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         October 29, 2013
!
!  Language:     Fortran-90
!
!  Version:      1.00a
!
!  Description:  Piecewise liner interpolation.
!
!  Notes:        Piecewise linear interpolation.  Given input arrays XX (independent variable) and YY (dependent variable),
!                both of dimension NN, this routine finds, by linear interpolation, the value of Y(X).  Array XX must be in
!                ascending order.
!
!                The flag IERR is returned as -1 if X is below the low end of XX (an error), +1 if X is above the high end
!                of XX (also an error), or 0 if there was no error.
!***********************************************************************************************************************************

      FUNCTION LINTERPOL (NN, XX, YY, X, IERR) RESULT (Y)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NN                                                     ! dimension of XX and YY arrays
      DOUBLE PRECISION, DIMENSION(NN), INTENT(IN) :: XX, YY                         ! indep and dep variable arrays
      DOUBLE PRECISION, INTENT(IN) :: X                                             ! interpolate at X
      INTEGER, INTENT(OUT) :: IERR                                                  ! returned error code
      DOUBLE PRECISION :: Y                                                         ! interpolated value Y(X)
      INTEGER :: I                                                                  ! loop counter


      IF (X .LT. XX(1)) THEN                                                        ! if below low end of XX (error)..
         Y = YY(1)                                                                  !  set Y = first YY value
         IERR = -1                                                                  !  return error code -1
      ELSE IF (X .GT. XX(NN)) THEN                                                  ! if above high end of XX (error)..
         Y = YY(NN)                                                                 !  set Y = last YY value
         IERR = +1                                                                  !  return error code +1
      ELSE                                                                          ! if OK
         DO I = 2, NN                                                               ! loop to find first XX > X
            IF (XX(I) .GT. X) EXIT
         END DO
         Y = (YY(I)-YY(I-1))/(XX(I)-XX(I-1))*(X-XX(I-1))+YY(I-1)                    ! interpolate
         IERR = 0                                                                   ! set error code to 0 (no error)
      END IF

      RETURN

      END FUNCTION LINTERPOL

