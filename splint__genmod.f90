        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 21 22:25:03 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SPLINT__genmod
          INTERFACE 
            SUBROUTINE SPLINT(Y,Z,N)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: Y(N)
              REAL(KIND=8), INTENT(OUT) :: Z
            END SUBROUTINE SPLINT
          END INTERFACE 
        END MODULE SPLINT__genmod
