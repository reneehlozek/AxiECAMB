        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 21 22:25:03 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DVERK__genmod
          INTERFACE 
            SUBROUTINE DVERK(EV,N,FCN,X,Y,XEND,TOL,IND,C,NW,W)
              INTEGER(KIND=4) :: NW
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: EV
              EXTERNAL FCN
              REAL(KIND=8) :: X
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: XEND
              REAL(KIND=8) :: TOL
              INTEGER(KIND=4) :: IND
              REAL(KIND=8) :: C(*)
              REAL(KIND=8) :: W(NW,9)
            END SUBROUTINE DVERK
          END INTERFACE 
        END MODULE DVERK__genmod
