        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb  5 15:03:42 2015
        MODULE MA28ID__genmod
          INTERFACE 
            SUBROUTINE MA28ID(N,NZ,AORG,IRNORG,ICNORG,LICN,A,ICN,IKEEP, &
     &RHS,X,R,W,MTYPE,PREC,IFLAG)
              INTEGER(KIND=4) :: LICN
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: AORG(NZ)
              INTEGER(KIND=4) :: IRNORG(NZ)
              INTEGER(KIND=4) :: ICNORG(NZ)
              REAL(KIND=8) :: A(LICN)
              INTEGER(KIND=4) :: ICN(LICN)
              INTEGER(KIND=4) :: IKEEP(N,5)
              REAL(KIND=8) :: RHS(N)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: R(N)
              REAL(KIND=8) :: W(N)
              INTEGER(KIND=4) :: MTYPE
              REAL(KIND=8) :: PREC
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE MA28ID
          END INTERFACE 
        END MODULE MA28ID__genmod
