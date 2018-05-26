        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb  5 15:03:42 2015
        MODULE MA28BD__genmod
          INTERFACE 
            SUBROUTINE MA28BD(N,NZ,A,LICN,IVECT,JVECT,ICN,IKEEP,IW,W,   &
     &IFLAG)
              INTEGER(KIND=4) :: LICN
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LICN)
              INTEGER(KIND=4) :: IVECT(NZ)
              INTEGER(KIND=4) :: JVECT(NZ)
              INTEGER(KIND=4) :: ICN(LICN)
              INTEGER(KIND=4) :: IKEEP(N,5)
              INTEGER(KIND=4) :: IW(N,5)
              REAL(KIND=8) :: W(N)
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE MA28BD
          END INTERFACE 
        END MODULE MA28BD__genmod
