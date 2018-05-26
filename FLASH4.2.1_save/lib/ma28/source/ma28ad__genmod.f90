        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb  5 15:03:42 2015
        MODULE MA28AD__genmod
          INTERFACE 
            SUBROUTINE MA28AD(N,NZ,A,LICN,IRN,LIRN,ICN,U,IKEEP,IW,W,    &
     &IFLAG)
              INTEGER(KIND=4) :: LIRN
              INTEGER(KIND=4) :: LICN
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NZ
              REAL(KIND=8) :: A(LICN)
              INTEGER(KIND=4) :: IRN(LIRN)
              INTEGER(KIND=4) :: ICN(LICN)
              REAL(KIND=8) :: U
              INTEGER(KIND=4) :: IKEEP(N,5)
              INTEGER(KIND=4) :: IW(N,8)
              REAL(KIND=8) :: W(N)
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE MA28AD
          END INTERFACE 
        END MODULE MA28AD__genmod
