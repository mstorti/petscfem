C=============================================================== NUEVOSC
C**** INTRD0A
C
      INTEGER*4      ISSTEP,NSSTEP,ISSTEPT,NSSTEPT
C
      COMMON/INTRD0A/ISSTEP,NSSTEP,ISSTEPT,NSSTEPT
C--------------------------------------------------------------- NUEVOSC
C**** OUTPUT & OUTPUTT
C
      INTEGER*4      NCKGLO
C
      COMMON/OUTPUTA/NCKGLO
C--------------------------------------------------------------- NUEVOSC
C**** OUTPUT & OUTPUTT
C
      INTEGER*4      ICOSTA,ITERME,ISTAGG,INICOU
      INTEGER*4      ITERMG,ITERMP,ITERMD
C
      COMMON/COSTAGG/ICOSTA,ITERME,ISTAGG,INICOU
      COMMON/COSTAGH/ITERMG,ITERMP,ITERMD
C--------------------------------------------------------------- NUEVOSC
C**** ELEMENT
C
      INTEGER*4      NPOINC,NELEMC
C
      COMMON/COSTAEL/NPOINC,NELEMC
C--------------------------------------------------------------- NUEVOSC
C**** ELEMENN
C
      REAL*8         COUFAC,CENKEL
      INTEGER*4      IITERC,INTERC
      INTEGER*4      NITERC,NFCOUC,ICONVC
      INTEGER*4      NFPCH
C
      COMMON/COSTAE1/COUFAC,CENKEL
      COMMON/COSTAE2/IITERC,INTERC
      COMMON/COSTAE3/NITERC,NFCOUC,ICONVC
      COMMON/COSTAE4/NFPCH
C--------------------------------------------------------------- NUEVOSC
C**** NAMECA
C
      CHARACTER*80  CCOA,CCOB
C
      COMMON/NAMECA/CCOA,CCOB
C
      INTEGER*4      INUS1C,ILUS1C
      COMMON/USEFULC/INUS1C,ILUS1C
C
C**** LOGUNC
C
C     To add units, see:                        ! to be confirmed
C
      INTEGER*4     LUDATC,LURESC
C
      COMMON/LOGUNC/LUDATC,LURESC
C=============================================================== NUEVOSC
