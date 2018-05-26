        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb  5 15:03:42 2015
        MODULE MC22AD__genmod
          INTERFACE 
            SUBROUTINE MC22AD(N,ICN,A,NZ,LENROW,IP,IQ,IW,IW1)
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ICN(NZ)
              REAL(KIND=8) :: A(NZ)
              INTEGER(KIND=4) :: LENROW(N)
              INTEGER(KIND=4) :: IP(N)
              INTEGER(KIND=4) :: IQ(N)
              INTEGER(KIND=4) :: IW(N,2)
              INTEGER(KIND=4) :: IW1(NZ)
            END SUBROUTINE MC22AD
          END INTERFACE 
        END MODULE MC22AD__genmod
