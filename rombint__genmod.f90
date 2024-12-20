        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 21 22:25:03 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ROMBINT__genmod
          INTERFACE 
            FUNCTION ROMBINT(F,A,B,TOL)
              REAL(KIND=8) :: F
              EXTERNAL F
              REAL(KIND=8), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8) :: ROMBINT
            END FUNCTION ROMBINT
          END INTERFACE 
        END MODULE ROMBINT__genmod
