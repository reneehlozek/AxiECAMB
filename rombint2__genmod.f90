        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 21 22:25:03 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ROMBINT2__genmod
          INTERFACE 
            FUNCTION ROMBINT2(F,A,B,TOL,MAXIT,MINSTEPS)
              REAL(KIND=8) :: F
              EXTERNAL F
              REAL(KIND=8), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B
              REAL(KIND=8), INTENT(IN) :: TOL
              INTEGER(KIND=4), INTENT(IN) :: MAXIT
              INTEGER(KIND=4), INTENT(IN) :: MINSTEPS
              REAL(KIND=8) :: ROMBINT2
            END FUNCTION ROMBINT2
          END INTERFACE 
        END MODULE ROMBINT2__genmod
