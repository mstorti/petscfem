C=============================================================== NUEVOSF
C**** INTRD0A
C
      INTEGER*4       ISSTEPF,NSSTEPF
C
      COMMON/INTRD0AF/ISSTEPF,NSSTEPF
C--------------------------------------------------------------- NUEVOSF
C**** OUTPUT & OUTPUTT
C
      INTEGER*4       NCKGLOF
C
      COMMON/OUTPUTAF/NCKGLOF
C--------------------------------------------------------------- NUEVOSF
C**** OUTPUT & OUTPUTT
C
      INTEGER*4       ITERMEF,ISTAGF,ITERMPF
C
      COMMON/COSTAFAF/ITERMEF,ISTAGF,ITERMPF
C--------------------------------------------------------------- NUEVOSF
C**** ELEMENT
C
c     INTEGER*4      NPOINC,NELEMC
C
c     COMMON/COSTAEL/NPOINC,NELEMC
C--------------------------------------------------------------- NUEVOSF
C**** ELEMENN
C
c     REAL*8         COUFAC,CENKEL
c     INTEGER*4      IITERC
c     INTEGER*4      NITERC,NFCOUC,ICONVC
c     INTEGER*4      NFPCH
C
c     COMMON/COSTAE1/COUFAC,CENKEL
c     COMMON/COSTAE2/IITERC
c     COMMON/COSTAE3/NITERC,NFCOUC,ICONVC
c     COMMON/COSTAE4/NFPCH
C--------------------------------------------------------------- NUEVOSF
C**** IMPROVED COUPLING ALGORITHM & VISCOUS HEAT SOURCE
C
      REAL*8         GRATEFF,TREFEFF,GRAVYFF,GVECTFF(3),FPARATT(20),
     .               VISC1F,VISC2F,DILA1F,DILA2F
C
      COMMON/COUPMAT/GRATEFF,TREFEFF,GRAVYFF,GVECTFF,FPARATT,
     .               VISC1F,VISC2F,DILA1F,DILA2F
C
      INTEGER*4      KCONST,KVILAT
C
      COMMON/VISCOHE/KCONST,KVILAT
C--------------------------------------------------------------- NUEVOSF
C**** CHECKING OPERATIONS FOR THERMALLY-COUPLED FLOWS
C
      REAL*8          VELMAC,TEMINC,TEMAXC
C
      COMMON/CHECKOP1/VELMAC,TEMINC,TEMAXC
C
      INTEGER*4       ICHECKOP
C
      COMMON/CHECKOP2/ICHECKOP
C--------------------------------------------------------------- NUEVOSF
C**** NAMECA
C
c     CHARACTER*80  CCOA,CCOB
C
c     COMMON/NAMECA/CCOA,CCOB
C
c     INTEGER*4      INUS1C,ILUS1C
c     COMMON/USEFULC/INUS1C,ILUS1C
C
C**** LOGUNC
C
C     To add units, see:                        ! to be confirmed
C
c     INTEGER*4     LUDATC,LURESC
C
c     COMMON/LOGUNC/LUDATC,LURESC
C=============================================================== NUEVOSF
