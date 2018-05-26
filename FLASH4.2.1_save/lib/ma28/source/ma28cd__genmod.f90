        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb  5 15:03:42 2015
        MODULE MA28CD__genmod
          INTERFACE 
            SUBROUTINE MA28CD(N,A,LICN,ICN,IKEEP,RHS,W,MTYPE)
              INTEGER(KIND=4) :: LICN
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LICN)
              INTEGER(KIND=4) :: ICN(LICN)
              INTEGER(KIND=4) :: IKEEP(N,5)
              REAL(KIND=8) :: RHS(N)
              REAL(KIND=8) :: W(N)
              INTEGER(KIND=4) :: MTYPE
            END SUBROUTINE MA28CD
          END INTERFACE 
        END MODULE MA28CD__genmod
