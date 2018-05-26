        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb  5 15:03:42 2015
        MODULE MC23AD__genmod
          INTERFACE 
            SUBROUTINE MC23AD(N,ICN,A,LICN,LENR,IDISP,IP,IQ,LENOFF,IW,  &
     &IW1)
              INTEGER(KIND=4) :: LICN
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ICN(LICN)
              REAL(KIND=8) :: A(LICN)
              INTEGER(KIND=4) :: LENR(N)
              INTEGER(KIND=4) :: IDISP(2)
              INTEGER(KIND=4) :: IP(N)
              INTEGER(KIND=4) :: IQ(N)
              INTEGER(KIND=4) :: LENOFF(N)
              INTEGER(KIND=4) :: IW(N,5)
              INTEGER(KIND=4) :: IW1(N,2)
            END SUBROUTINE MC23AD
          END INTERFACE 
        END MODULE MC23AD__genmod
