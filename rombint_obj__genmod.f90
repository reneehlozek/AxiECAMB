        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 21 22:25:03 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ROMBINT_OBJ__genmod
          INTERFACE 
            FUNCTION ROMBINT_OBJ(OBJ,F,A,B,TOL,MAXIT)
              REAL(KIND=4) :: OBJ
              REAL(KIND=8) :: F
              EXTERNAL F
              REAL(KIND=8), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B
              REAL(KIND=8), INTENT(IN) :: TOL
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: MAXIT
              REAL(KIND=8) :: ROMBINT_OBJ
            END FUNCTION ROMBINT_OBJ
          END INTERFACE 
        END MODULE ROMBINT_OBJ__genmod
