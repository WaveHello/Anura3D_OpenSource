        !COMPILER-GENERATED INTERFACE MODULE: Sat Sep  2 13:56:20 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ASSERT__genmod
          INTERFACE 
            SUBROUTINE ASSERT(CONDITION,MESSAGE)
              LOGICAL(KIND=4), INTENT(IN) :: CONDITION
              CHARACTER(*), INTENT(IN) :: MESSAGE
            END SUBROUTINE ASSERT
          END INTERFACE 
        END MODULE ASSERT__genmod
