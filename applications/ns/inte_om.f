C=============================================================== INTERVL
C**** CONSTNA
      REAL*8         DITER,DTIME,GRAVY,GVECT,STIFI,
     .               TALFA,TBETA,TGAMA,TDELT,TOLER,XTIME,WLUMP
C
      COMMON/CONSTNA/DITER,DTIME,GRAVY,GVECT(3),STIFI,
     .               TALFA,TBETA,TGAMA,TDELT,TOLER,XTIME,WLUMP
C--------------------------------------------------------------- INTERVL
C**** CURRENA
C**** CURRENB
      INTEGER*4      ITIME,ISTEP,IITER,KDTIM,KRESL,KSTIF,KUNLD,
     .               NCHEK
      REAL*8         ARCLN,FACTO,FACPR,PITER,PVALU,STICU,
     .               STIIN,TFACT,TTIME
C
      COMMON/CURRENA/ITIME,ISTEP,IITER,KDTIM,KRESL,KSTIF,KUNLD,
     .               NCHEK
      COMMON/CURRENB/ARCLN,FACTO,FACPR,PITER,PVALU,STICU,
     .               STIIN,TFACT,TTIME
C--------------------------------------------------------------- INTERVL
C**** INTERVA
      INTEGER*4      IMPOS,KALGO,KARCL,KINTE,KOPTI,KCONV,KSAVE,
     .               LACCE,LAUTO,LINES,MITER,NALGO,NBACK,NCDIS,
     .               NSTEP,NOUTP
C
      COMMON/INTERVA/IMPOS,KALGO,KARCL,KINTE,KOPTI,KCONV,KSAVE,
     .               LACCE,LAUTO,LINES,MITER,NALGO,NBACK,NCDIS,
     .               NSTEP,NOUTP(50)
C--------------------------------------------------------------- INTERVL
C**** INTERVA
      INTEGER*4       MKONT,NFORZ,MSUBP
      REAL*8          TOPLA,ALFAP,EXCTP
C
      COMMON/RETURNPA/MKONT,NFORZ,MSUBP
      COMMON/RETURNPB/TOPLA,ALFAP,EXCTP
C=============================================================== INTERVL
