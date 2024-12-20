        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 21 22:25:03 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SPLINE_DERIV__genmod
          INTERFACE 
            SUBROUTINE SPLINE_DERIV(X,Y,Y2,Y1,N)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: X(N)
              REAL(KIND=8), INTENT(IN) :: Y(N)
              REAL(KIND=8), INTENT(IN) :: Y2(N)
              REAL(KIND=8), INTENT(OUT) :: Y1(N)
            END SUBROUTINE SPLINE_DERIV
          END INTERFACE 
        END MODULE SPLINE_DERIV__genmod
