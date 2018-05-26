        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb  5 15:03:42 2015
        MODULE MA30CD__genmod
          INTERFACE 
            SUBROUTINE MA30CD(N,ICN,A,LICN,LENR,LENRL,LENOFF,IDISP,IP,IQ&
     &,X,W,MTYPE)
              INTEGER(KIND=4) :: LICN
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ICN(LICN)
              REAL(KIND=8) :: A(LICN)
              INTEGER(KIND=4) :: LENR(N)
              INTEGER(KIND=4) :: LENRL(N)
              INTEGER(KIND=4) :: LENOFF(N)
              INTEGER(KIND=4) :: IDISP(2)
              INTEGER(KIND=4) :: IP(N)
              INTEGER(KIND=4) :: IQ(N)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: W(N)
              INTEGER(KIND=4) :: MTYPE
            END SUBROUTINE MA30CD
          END INTERFACE 
        END MODULE MA30CD__genmod
