        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb  5 15:03:42 2015
        MODULE MC21A__genmod
          INTERFACE 
            SUBROUTINE MC21A(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
              INTEGER(KIND=4) :: LICN
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ICN(LICN)
              INTEGER(KIND=4) :: IP(N)
              INTEGER(KIND=4) :: LENR(N)
              INTEGER(KIND=4) :: IPERM(N)
              INTEGER(KIND=4) :: NUMNZ
              INTEGER(KIND=4) :: IW(N,4)
            END SUBROUTINE MC21A
          END INTERFACE 
        END MODULE MC21A__genmod
