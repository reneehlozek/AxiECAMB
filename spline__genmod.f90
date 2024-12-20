        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 21 22:25:03 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SPLINE__genmod
          INTERFACE 
            SUBROUTINE SPLINE(X,Y,N,D11,D1N,D2)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: X(N)
              REAL(KIND=8), INTENT(IN) :: Y(N)
              REAL(KIND=8), INTENT(IN) :: D11
              REAL(KIND=8), INTENT(IN) :: D1N
              REAL(KIND=8), INTENT(OUT) :: D2(N)
            END SUBROUTINE SPLINE
          END INTERFACE 
        END MODULE SPLINE__genmod
