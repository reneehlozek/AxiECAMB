        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 21 22:25:03 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SPLINE_INTEGRATE__genmod
          INTERFACE 
            SUBROUTINE SPLINE_INTEGRATE(X,Y,Y2,YINT,N)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: X(N)
              REAL(KIND=8), INTENT(IN) :: Y(N)
              REAL(KIND=8), INTENT(IN) :: Y2(N)
              REAL(KIND=8), INTENT(OUT) :: YINT(N)
            END SUBROUTINE SPLINE_INTEGRATE
          END INTERFACE 
        END MODULE SPLINE_INTEGRATE__genmod
