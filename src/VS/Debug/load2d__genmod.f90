        !COMPILER-GENERATED INTERFACE MODULE: Thu Aug 25 15:09:19 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LOAD2D__genmod
          INTERFACE 
            SUBROUTINE LOAD2D(RLOAD,NDOF,COORD,NUMBEROFINTEGRATIONPOINTS&
     &,ILOADCON,LOADVALUE,LOADTYPE)
              USE MODELEMENTEVALUATION
              INTEGER(KIND=4), INTENT(IN) :: NUMBEROFINTEGRATIONPOINTS
              REAL(KIND=8), INTENT(INOUT) :: RLOAD(COUNTERS%NODTOT*     &
     &NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: NDOF(COUNTERS%NODTOT)
              REAL(KIND=8), INTENT(IN) :: COORD(COUNTERS%NODTOT,NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: ILOADCON(                  &
     &ELEMENTBOUNDARYNODES_XI)
              REAL(KIND=8), INTENT(IN) :: LOADVALUE(                    &
     &ELEMENTBOUNDARYNODES_XI,NVECTOR)
              INTEGER(KIND=4), INTENT(IN) :: LOADTYPE
            END SUBROUTINE LOAD2D
          END INTERFACE 
        END MODULE LOAD2D__genmod