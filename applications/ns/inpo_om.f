C=============================================================== INP_OUT
C**** LISTENA:
C**** LISTENB:
C**** LISTEND:
        INTEGER*4 MAXWP
        PARAMETER (MAXWP=50)
C
        INTEGER*4      NNWOR,NNPAR
        REAL*8         PARAM(MAXWP),DPARA(MAXWP)
        CHARACTER*5    WORDS(MAXWP),DWORD(MAXWP)
C
        COMMON/LISTENA/NNWOR,NNPAR
        COMMON/LISTENB/PARAM,DPARA
        COMMON/LISTEND/WORDS,DWORD
C--------------------------------------------------------------- INP_OUT
C**** PLOTERA:
C**** PLOTERB:
        INTEGER*4   MMCUR,MSPLO          ! Max curves & Size of MPLOT
        PARAMETER   (MMCUR=40,MSPLO=240) ! MSPLO=MMCUR*6
C
        INTEGER*4     NCOLD,NCURV,NPONT(MMCUR),MPLOT(MMCUR,2,3)
        CHARACTER*8   FORMA
C
        INTEGER*4     IPRCO
C
        COMMON/PLOTERA/NCOLD,NCURV,NPONT,MPLOT
        COMMON/PLOTERB/FORMA
        COMMON/PLOTERC/IPRCO
C--------------------------------------------------------------- INP_OUT
C**** PRIOUTA:
        INTEGER*4      KFEMV,KPRI0,KPRI1,KPRI2,KPRI3,KPRI4,KPRI5,KPRI6,
     .                             KPRI7,KPRI8,KPRI9,
     .                 KPRI10,KPRI11,KPRI12,KPRI13,KPRI15,KPRI16,KPRI17,
     .                 KPRI18,KPRI19,KPRI20
C
        COMMON/PRIOUTA/KFEMV,KPRI0,KPRI1,KPRI2,KPRI3,KPRI4,KPRI5,KPRI6,
     .                             KPRI7,KPRI8,KPRI9,
     .                 KPRI10,KPRI11,KPRI12,KPRI13,KPRI15,KPRI16,KPRI17,
     .                 KPRI18,KPRI19,KPRI20
C=============================================================== INP_OUT
