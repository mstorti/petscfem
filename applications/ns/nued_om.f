C=============================================================== NUEVOSD
C**** FUNDAMI
C
      INTEGER*4      IMICR,IMICO
C
      COMMON/FUNDAMI/IMICR,IMICO
C
      INTEGER*4       NPROPM,NHISTM
      INTEGER*4       NNUPM,NNUPO
C
      COMMON/FUNDAMIA/NPROPM,NHISTM
      COMMON/FUNDAMIB/NNUPM,NNUPO
C--------------------------------------------------------------- NUEVOSD
C**** INTRD0AS
C
      INTEGER*4       NSSTEPS
C
      COMMON/INTRD0AS/NSSTEPS
C--------------------------------------------------------------- NUEVOSD
C**** OUTPUT & OUTPUTT
C
c     INTEGER*4      NCKGLOS
C
c     COMMON/OUTPUSA/NCKGLOS
C--------------------------------------------------------------- NUEVOSD
C**** OUTPUT & OUTPUTT
C
      INTEGER*4       ICOSTAS,ITERMES,ISTAGGS
      INTEGER*4       ITERMGS,ITERMPS,ITERMDS
C
      COMMON/COSTAGGS/ICOSTAS,ITERMES,ISTAGGS
      COMMON/COSTAGHS/ITERMGS,ITERMPS,ITERMDS
C--------------------------------------------------------------- NUEVOSD
C**** ELEMENT
C
c     INTEGER*4      NPOINC,NELEMC
C
c     COMMON/COSTAEL/NPOINC,NELEMC
C--------------------------------------------------------------- NUEVOSD
C**** ELEMENN
C
c     REAL*8         COUFAC,CENKEL
      INTEGER*4      IITERCS,INTERCS
      INTEGER*4      NITERCS,NFCOUCS,ICONVCS
c     INTEGER*4      NFPCH
C
c     COMMON/COSTAE1S/COUFAC,CENKEL
      COMMON/COSTAE2S/IITERCS,INTERCS
      COMMON/COSTAE3S/NITERCS,NFCOUCS,ICONVCS
c     COMMON/COSTAE4S/NFPCH
C--------------------------------------------------------------- NUEVOSD
C**** NAMECSA
C
c     CHARACTER*80   CCOA,CCOB
C
c     COMMON/NAMECSA/CCOA,CCOB
C
C**** LOGUNCS
C
C     To add units, see:                        ! to be confirmed
C
c     INTEGER*4      LUDATCS,LURESCS
C
c     COMMON/LOGUNCS/LUDATCS,LURESCS
C=============================================================== NUEVOSD
