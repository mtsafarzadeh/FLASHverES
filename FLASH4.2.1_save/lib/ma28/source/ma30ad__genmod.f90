        !COMPILER-GENERATED INTERFACE MODULE: Thu Feb  5 15:03:42 2015
        MODULE MA30AD__genmod
          INTERFACE 
            SUBROUTINE MA30AD(NN,ICN,A,LICN,LENR,LENRL,IDISP,IP,IQ,IRN, &
     &LIRN,LENC,IFIRST,LASTR,NEXTR,LASTC,NEXTC,IPTR,IPC,U,IFLAG)
              INTEGER(KIND=4) :: LIRN
              INTEGER(KIND=4) :: LICN
              INTEGER(KIND=4) :: NN
              INTEGER(KIND=4) :: ICN(LICN)
              REAL(KIND=8) :: A(LICN)
              INTEGER(KIND=4) :: LENR(NN)
              INTEGER(KIND=4) :: LENRL(NN)
              INTEGER(KIND=4) :: IDISP(2)
              INTEGER(KIND=4) :: IP(NN)
              INTEGER(KIND=4) :: IQ(NN)
              INTEGER(KIND=4) :: IRN(LIRN)
              INTEGER(KIND=4) :: LENC(NN)
              INTEGER(KIND=4) :: IFIRST(NN)
              INTEGER(KIND=4) :: LASTR(NN)
              INTEGER(KIND=4) :: NEXTR(NN)
              INTEGER(KIND=4) :: LASTC(NN)
              INTEGER(KIND=4) :: NEXTC(NN)
              INTEGER(KIND=4) :: IPTR(NN)
              INTEGER(KIND=4) :: IPC(NN)
              REAL(KIND=8) :: U
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE MA30AD
          END INTERFACE 
        END MODULE MA30AD__genmod
