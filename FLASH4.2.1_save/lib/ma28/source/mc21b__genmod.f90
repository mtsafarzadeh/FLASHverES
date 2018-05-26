        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb  5 15:03:42 2015
        MODULE MC21B__genmod
          INTERFACE 
            SUBROUTINE MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,  &
     &OUT)
              INTEGER(KIND=4) :: LICN
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ICN(LICN)
              INTEGER(KIND=4) :: IP(N)
              INTEGER(KIND=4) :: LENR(N)
              INTEGER(KIND=4) :: IPERM(N)
              INTEGER(KIND=4) :: NUMNZ
              INTEGER(KIND=4) :: PR(N)
              INTEGER(KIND=4) :: ARP(N)
              INTEGER(KIND=4) :: CV(N)
              INTEGER(KIND=4) :: OUT(N)
            END SUBROUTINE MC21B
          END INTERFACE 
        END MODULE MC21B__genmod
