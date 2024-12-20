        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec 21 22:25:03 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SPLINE_OUT__genmod
          INTERFACE 
            SUBROUTINE SPLINE_OUT(XARR,YARR,YARR_BUFF,N,X,Y)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: XARR(N)
              REAL(KIND=8) :: YARR(N)
              REAL(KIND=8) :: YARR_BUFF(N)
              REAL(KIND=8) :: X
              REAL(KIND=8) :: Y
            END SUBROUTINE SPLINE_OUT
          END INTERFACE 
        END MODULE SPLINE_OUT__genmod
