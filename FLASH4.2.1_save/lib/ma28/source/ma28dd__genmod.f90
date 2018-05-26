        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb  5 15:03:42 2015
        MODULE MA28DD__genmod
          INTERFACE 
            SUBROUTINE MA28DD(N,A,LICN,IVECT,JVECT,NZ,ICN,LENR,LENRL,   &
     &LENOFF,IP,IQ,IW1,IW,W1,IFLAG)
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: LICN
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LICN)
              INTEGER(KIND=4) :: IVECT(NZ)
              INTEGER(KIND=4) :: JVECT(NZ)
              INTEGER(KIND=4) :: ICN(LICN)
              INTEGER(KIND=4) :: LENR(N)
              INTEGER(KIND=4) :: LENRL(N)
              INTEGER(KIND=4) :: LENOFF(N)
              INTEGER(KIND=4) :: IP(N)
              INTEGER(KIND=4) :: IQ(N)
              INTEGER(KIND=4) :: IW1(N,3)
              INTEGER(KIND=4) :: IW(N,2)
              REAL(KIND=8) :: W1
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE MA28DD
          END INTERFACE 
        END MODULE MA28DD__genmod
