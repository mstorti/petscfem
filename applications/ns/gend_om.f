C***********************************************************************
C
C**** GENERAL DIMENSIONS OF PROGRAM VULCAN (mechanical variables)
C
C     IF ANY CHANGE IS PRODUCED, COMPILE:
C
C       adddat.f addelm.f addpri.f addsol.f addwor.f
C       datbas.f memodt.f (comp-m)
C
C     IF "segmentation fault" is produced just starting VULCAN, check
C     maximum dimensions
C
C     Command "max" works for every MMACHIM except for MMACHIM=5 & 8
C     (see also para_om.f)
C
C***********************************************************************
C
C     PROBLEM: crankshaft
C
C***********************************************************************
C
C**** WRITE MAXIMUM DIMENSIONS (0=write, 1=write & stop, 2=write more)
C
      PARAMETER(
     .   MWRITEM=2 )
C
C**** GENERAL DIMENSIONS (change if necessary; see input data)
C
      PARAMETER(
     .   MPOIN=1500,                       ! including MPOIC
     .   MELEM=1000,
     .   MDIME=3,
     .   MNODE=20,
     .   MGAUS=27,
     .   MGRUP=20,
     .   MMATS=20,
     .   MFUNC=2 )
C
      PARAMETER(
     .   MDYNA=1,            ! 0:static; 1:dynamic
     .   MSMOM=1,            ! smoothing (0:no; 1:yes)
     .   MRENU=0,            ! renumbering (0:no; 1:yes)
     .   MHOUR=0,            ! hourglass control (0:no; 1:yes)
     .   MPORE=0,            ! pore water pressure prob. (0:no; 2:yes)
     .   MSKEW=10 )          ! skew system (0:no; n:yes)
C
      PARAMETER(
     .   MACTI=1 )           ! 0: no active elements; 1: active elements
C
C**** COUPLING (thermal-mechanical) DIMENSIONS
C     (change if necessary; see input data)
C
      PARAMETER(
     .   MTERMEM=1,  !-1: pure mech. 0:no mech. coupling; 1:mech. coupl.
     .   MFPCHM=6,   ! 2*number of phase-changes (see setdatt.f)
     .   MITERCM=0 ) ! 0:standard stag. scheme; 1:improved stag. scheme
C
C**** COUPLING (thermal-microstructural) DIMENSIONS
C     (change if necessary; see input data)
C
      PARAMETER(
     .   MMICRM=1,       ! 0:no microstructural coupling; 1:micr. coupl.
     .   MMICOM=0 )      ! 0:weak micr. coupling; 1:full micr. coupl.
C
C**** OTHER DIMENSIONS (generally, do not change)
C
C     Only change:
C
C     MPROP:  mechanical & microstructural proerties
C     MPROPM: mechanical properties
C     MPROP-MPROPM: microstructural properties
C
C     In a mechanical analyisis: MPROP=MPROPM
C
C     Recommended values:
C
C     MPROP=1000 & MPROPM=500 (mechanical & microstructural analysis)
C     MPROP=500 (mechanical analysis only)
C
      PARAMETER(
     .   MSTR1=2*MDIME,        ! nstr1=2*ndime; see conset
     .   MPROP=500+MMICRM*500, ! see setdat & setdats (1-500; 501-1000)
     .   MDOFC=MDIME,          ! always for solids elements; see conset
     .   MHLOD=5,              ! see setdat
     .   MSUBF=50,             ! see setdat
     .   MPREL=18,             ! see setdat
     .   MHIST=36,             ! 3D symm. const. tensor; see conset
     .   MNUIN=23 )            ! comp. of int. var. to print; see setdat
C
C**** DEGREES OF FREEDOMS (do not change)
C
      PARAMETER(
     .   MEVAB=MNODE*MDOFC,
     .   MTOTV=MPOIN*MDOFC )
C
C***********************************************************************
C
C**** ADDITIONAL PARAMETERS OF PROGRAM VULCAN (thermal problem)
C
C     it should be read in the input data !!!  (see check0.f) ctm
C
C***********************************************************************
C
C**** ADDITIONAL PARAMETERS
C
C               MDISKDM=0: database is out of core (datbas.f)
C               MDISKDM=1: database is in core (datbas.f)
C
C               MFURESM=0: a future restart will not be made (datrst.f)
C               MFURESM>0: a future restart will be made (datrst.f)
C
C               MMACHIM=1: CONVEX
C               MMACHIM=2: SILICON GRAPHICS
C               MMACHIM=3: VAX COMPUTER (not implemented)
C               MMACHIM=4: SUN
C               MMACHIM=5: PERSONAL COMPUTER
C               MMACHIM=6: SILICON GRAPHICS POWER CHALLENGE
C               MMACHIM=7: HEWLETT PACKARD
C
C     THEORETICAL DEFAULTS: MDISKDM=0, MMURESM=1, MMACHIM=1
C     REAL DEFAULTS:        MDISKDM=1, MFURESM=0, MMACHIM=2
C
C***********************************************************************
C
      PARAMETER(
     .   MDISKDM=1,
     .   MFURESM=1,
     .   MMACHIM=8 )
C
      PARAMETER(
c    .   MCHA1=max((1+(MMACHIM/6)*6-MMACHIM),0),     ! 1 for MMACHI=6
c    .   MCHA2=max((1+(MMACHIM/10)*10-MMACHIM),0),   ! 1 for MMACHI=10
c    .   MCHAA=max(MCHA1,MCHA2),
c    .   MCHAL=max(4,8*MCHAA) )
     .   MCHAL=4 )                                   ! PC & linux
C
C***********************************************************************
C
C**** FULL OR PARTIAL MEMORY (MMEMOM=1: full; MMEMOM=0: partial)
C
C     MMEMOM concerns to: EPMTX (elastic const. tensor)
C                         RMAT1 (anisotropic const. tensor)
C                         STRA0 (initial strains)
C                         STRS0 (initial stresses)
C                         TEMPC (temperature)
C
C     (change if necessary; see input data [now in setdat])
C
C***********************************************************************
C
      PARAMETER(
     .   MMEMOM=0 )
C
C***********************************************************************
C
C**** ADDITIONAL MEMORY PARAMETER 1
C
C     MMEMO1M concerns to: COORD(NDIME,NNODE) as a global array
C
C               NMEMO1M=0: COORD is an elemental array (ELCOD)
C               NMEMO1M=1: COORD is a global array
C
C     MMEMA1M=1-MMEMO1M (do not change)
C
C     THEORETICAL DEFAULT: NMEMO1M=0
C
C
C     MMEMO2M concerns to: shape functions, cartesian derivatives,
C                          etc. computed initially or every time
C                          when needed
C
C               MMEMO2M=0: SHAPE, CARTD, etc. computed initially
C               MMEMO2M=1: SHAPE, CARTD, etc. computed every time
C
C     THEORETICAL DEFAULT: MMEMO2M=0
C
C
C     MMEMO3M concerns to: ---
C
C
C     MMEMO4M concerns to: --- 
C
C
C     MMEMO5M concerns to: ELDIS in/out of ELVAR
C
C               MMEMO5M=0: ELDIS in ELVAR
C               MMEMO5M=1: ELDIS out of ELVAR
C
C     THEORETICAL DEFAULT: MMEMO5M=0
C
C
C     MMEMO6M concerns to: elemental assembly process in other/same
C                          NELEM loop as the evaluation of each
C                          contribution of the jacobian matrix
C
C               MMEMO6M=0: elemental assembly process in other loop
C               MMEMO6M=1: elemental assembly process in the same loop
C
C     THEORETICAL DEFAULT: MMEMO6M=0
C
C
C     MMEMO7M concerns to: the jacobian matrix is (not) evaluated in
C                          the NELEM loop in the solver routines. 
C
C               MMEMO7M=0: jacobian not evaluated in NELEM solver loop
C               MMEMO7M=1: jacobian evaluated in NELEM solver loop. For
C                         this case, NMEMO6 must be equals 1
C
C     THEORETICAL DEFAULT: MMEMO7M=0
C
C
C     MMEMO8M concerns to: the mass matrix is computed in other/same
C                          NELEM loop as the stiffness matrix
C
C               MMEMO8M=0: mass matrix computed in other loop (MMEMO6M
C                          must be equals 0)
C               MMEMO8M=1: mass matrix computed in the same (stiffness)
C                          loop
C
C     THEORETICAL DEFAULT: NMEMO8M=0
C
C
C     MMEMO9M concerns to: the second component of DISIT
C
C               MMEMO9M=0: second component of DISIT not used
C               MMEMO9M=1: second component of DISIT used
C
C     THEORETICAL DEFAULT: MMEMO9M=0
C
C
C     MMEMO10M concerns to: ---
C
C
C     MMEMO11M concerns to: normal gap & normal pressure vectors
C                           considered in mechanical problem are
C                           necessary for printing (average values for
C                           contact elements ITYPE=4 or ITYPE=32)
C
C     => MMEMO11M is not used
C
C
C     (change if necessary; see input data [now in setdat])
C
C***********************************************************************
C
      PARAMETER(
     .   MMEMO1M=1,
     .   MMEMO2M=1,
     .   MMEMO5M=1,
     .   MMEMO6MA=1,
     .   MMEMO7M=0,
     .   MMEMO8M=1,
     .   MMEMO9M=0 )
C
C***********************************************************************
C
C     SOLVER VARIABLES
C
C     1) DEFAULTS OR CHOSEN VARIABLES
C     2) SOLVER TO BE USED
C     3) SYMMETRIC OR UNSYMMETRIC CASE
C
C     (change if necessary; see input data)
C
C***********************************************************************
C
C**** 1) DEFAULTS (MDEF1=1 & MDEF2=0) OR
C        CHOSEN VARIABLES (MDEF1=0 & MDEF2=1)
C
C        (change if necessary; see input data)
C
C***********************************************************************
C
      PARAMETER(
     .   MDEF1=0,
     .   MDEF2=1 )
C
C***********************************************************************
C
C**** 2) SOLVER TO BE USED:
C
C        MSOL1=1 & MSOL2=0 & MSOL3=0: SKYLINE SOLVER
C        MSOL1=0 & MSOL2=1 & MSOL3=0: FRONTAL SOLVER
C        MSOL1=0 & MSOL2=0 & MSOL3=1: PCG SOLVER
C
C        (change if necessary; see input data)
C
C***********************************************************************
C
      PARAMETER(
     .   MSOL1=0,
     .   MSOL2=1,
     .   MSOL3=0 )
C
C***********************************************************************
C
C**** 3) SYMMETRIC CASE (MUNS1=1 & MUNS2=0) (skyline, frontal & pcg)
C        OR
C        UNSYMMETRIC CASE (MUNS1=0 & MUNS2=1) (skyline & frontal)
C
C        (change if necessary; see input data)
C
C***********************************************************************
C
      PARAMETER(
     .   MUNS1=0,
     .   MUNS2=1 )
C
C***********************************************************************
C
C     (from here onwards change only for chosen variables of the 
C      chosen solver)
C
C***********************************************************************
C
C**** SKYLINE SOLVER
C
C     VARIABLES:
C     MEQNS: number of equations
C     MLAST: bandwidth
C
      PARAMETER(
     .   MEQNS1=MDEF1*MTOTV,                  ! default
     .   MEQNS2=MDEF2*1 )                     ! chosen variable
C
      PARAMETER(
c    .   MEQNS=max(MEQNS1,MEQNS2) )
     .   MEQNS=(MEQNS1+MEQNS2) )              ! PC & linux
C
      PARAMETER(
     .   MLAST1=MDEF1*MEQNS*(MEQNS+1)/2,      ! default
     .   MLAST2=MDEF2*1 )                     ! chosen variable
C
      PARAMETER(
c    .   MLAST=max(MLAST1,MLAST2) )
     .   MLAST=(MLAST1+MLAST2) )              ! PC & linux
C
      PARAMETER(         ! mfite=0: direct solver
     .   MFITE=0 )       ! mfite=1: semidirect solver
C
C**** FRONTAL SOLVER
C
C     VARIABLES:
C     MFRON: frontwidth
C     MBUFA: buffer size
C
      PARAMETER(
     .   MFRON1=MDEF1*MTOTV,                  ! default
     .   MFRON2=MDEF2*1800 )                  ! chosen variable
C
      PARAMETER(
c    .   MFRON=max(MFRON1,MFRON2) )
     .   MFRON=(MFRON1+MFRON2) )              ! PC & linux
C
      PARAMETER(
     .   MBUFA1=MDEF1*(MFRON+1),              ! default
     .   MBUFA2=MDEF2*10 )                    ! chosen variable
C
      PARAMETER(
c    .   MBUFA=max(MBUFA1,MBUFA2) )
     .   MBUFA=(MBUFA1+MBUFA2) )              ! PC & linux
C
      PARAMETER(
     .   MSTIF1=MUNS1*MFRON*(MFRON+1)/2,      ! symmetric case
     .   MSTIF2=MUNS2*MFRON*MFRON  )          ! unsymmetric case 
C
      PARAMETER(
c    .   MSTIF=max(MSTIF1,MSTIF2) )
     .   MSTIF=(MSTIF1+MSTIF2) )              ! PC & linux
C
C**** PCG SOLVER
C
      PARAMETER(
     .   MWIDT1=0,                            ! default
     .   MWIDT2=MDEF2 )                       ! chosen variable
C
      PARAMETER(
c    .   MWIDT=max(MWIDT1,MWIDT2) )
     .   MWIDT=(MWIDT1+MWIDT2) )              ! PC & linux
C
      PARAMETER(
     .   MSIZE1=MDOFC,
     .   MSIZE2=MDOFC*MDOFC*MWIDT )
C
      PARAMETER(
c    .   MSIZE=max(MSIZE1,MSIZE2) )
     .   MSIZE=(MSIZE1+MSIZE2) )              ! PC & linux
C
C***********************************************************************
C
C**** ADDITIONAL MEMORY PARAMETERS (do not change)
C
C***********************************************************************
C
      PARAMETER(
c    .   MMEMO7MA=max(MMEMO7M,0),
c    .   MMEMO6M=max(MMEMO7MA,MMEMO6MA) )
     .   MMEMO6M=1 )                          ! PC & linux
C
      PARAMETER(
     .   MMEMA1M=1-MMEMO1M,
     .   MMEMA2M=1-MMEMO2M,
     .   MMEMA5M=1-MMEMO5M,
     .   MMEMA6M=1-MMEMO6M,
     .   MMEMA7M=1-MMEMO7M,
     .   MMEMA9M=1-MMEMO9M )
C
C***********************************************************************
C
C**** MATRICES DIMENSIONS (do not change)
C
C***********************************************************************
C
      PARAMETER(
     .   MEVAC=MEVAB,
     .   MKOVA1=MUNS1*MEVAB*(MEVAB+1)/2,      ! symmetric case
     .   MKOVA2=MUNS2*MEVAB*MEVAB )           ! unsymmetric case
C
      PARAMETER(
c    .   MKOVA=max(MKOVA1,MKOVA2) )
     .   MKOVA=(MKOVA1+MKOVA2) )              ! PC & linux
C
      PARAMETER(
     .   MKOND=MKOVA )
C
      PARAMETER(
     .   MKOST1=MUNS1*MSTR1*(MSTR1+1)/2,      ! symmetric case
     .   MKOST2=MUNS2*MSTR1*MSTR1 )           ! unsymmetric case
C
      PARAMETER(
c    .   MKOST=max(MKOST1,MKOST2) )
     .   MKOST=(MKOST1+MKOST2) )              ! PC & linux
C
