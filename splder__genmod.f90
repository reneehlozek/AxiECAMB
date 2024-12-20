        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 21 22:25:03 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SPLDER__genmod
          INTERFACE 
            SUBROUTINE SPLDER(Y,DY,N,G)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: Y(N)
              REAL(KIND=8), INTENT(OUT) :: DY(N)
              REAL(KIND=8), INTENT(IN) :: G(N)
            END SUBROUTINE SPLDER
          END INTERFACE 
        END MODULE SPLDER__genmod
