C=============================================================== PROBLEM
C**** DATBASA
C
      INTEGER*4      IDATP,IDATC,
     .               LENRC,NLENC,NLENP,NRECC,NRECG,NRECP,NWORP
C
      COMMON/DATBASA/IDATP(12,5),IDATC(12),
     .               LENRC,NLENC,NLENP,NRECC,NRECG,NRECP,NWORP
C--------------------------------------------------------------- PROBLEM
C**** DIMEMNA
C
      INTEGER*4      NDIME,NELEM,NFUNC,NGRUP,NHLOD,NHIST,NMATS,
     .               NPOIN,NPREL,NPROP,NSTR1,NTOTV,NSUBF,NTOTG,NPOIC
C
      COMMON/DIMENNA/NDIME,NELEM,NFUNC,NGRUP,NHLOD,NHIST,NMATS,
     .               NPOIN,NPREL,NPROP,NSTR1,NTOTV,NSUBF,NTOTG,NPOIC
C--------------------------------------------------------------- PROBLEM
C**** ELDATAA
C
      INTEGER*4      IDATA,IPREV,ISTAT,IMATX
C
      COMMON/ELDATAA/IDATA(15),IPREV(3),ISTAT(4),IMATX(7)
C--------------------------------------------------------------- PROBLEM
C**** ELEMNTA
C
      INTEGER*4      NDOFC,NDOFN,NEVAB,NEVAC,NGAUS,NKOND,NKOVA,NMOVA,
     .               NNODE,NDATA,NPREV,NSTAT,NMATX,NNODC
C
      COMMON/ELEMNTA/NDOFC,NDOFN,NEVAB,NEVAC,NGAUS,NKOND,NKOVA,NMOVA,
     .               NNODE,NDATA,NPREV,NSTAT,NMATX,NNODC
C--------------------------------------------------------------- PROBLEM
C**** HOURGLA
C**** HOURGLB
C
      INTEGER*4     NHOUR,KELAS
      REAL*8        HPARA
C
      COMMON/HOURGLA/NHOUR,KELAS
      COMMON/HOURGLB/HPARA
C--------------------------------------------------------------- PROBLEM
C**** PROBLMA
C
      INTEGER*4      KDYNA,KPORE,KPOST,KPROB,KSGAU,KSMUS,KTEMP,
     .               LARGE,NSKEW,LSKEW,NDISR,NDISO
C
      COMMON/PROBLMA/KDYNA,KPORE,KPOST,KPROB,KSGAU,KSMUS,KTEMP,
     .               LARGE,NSKEW,LSKEW,NDISR,NDISO
C--------------------------------------------------------------- PROBLEM
C**** SOLVERA
C**** SOLVERB
C
      INTEGER*4      KRENU,KSOLV,KSYMM,NWIDT,MITCG,NBUFA,NPRIR
      REAL*8         TOLCG,TOLC1
C
      COMMON/SOLVERA/KRENU,KSOLV,KSYMM,NWIDT,MITCG,NBUFA,NPRIR
      COMMON/SOLVERB/TOLCG,TOLC1
C--------------------------------------------------------------- PROBLEM
C**** TITLESA
C
      CHARACTER*8    TITLE,SUBTI
C
      COMMON/TITLESA/TITLE(8),SUBTI(8)
C--------------------------------------------------------------- PROBLEM
C**** ACTIVE ELEMENTS
C
      INTEGER*4     NACTI
      COMMON/ACTIV1/NACTI
C--------------------------------------------------------------- PROBLEM
C**** NAMEA
C
C    To add names, see: frodes.f, fronts.f & skydes.f
C
      CHARACTER*80 CA,CB,CC,CD,CE,CF,CG,CH,CI,CJ,CK,CL,
     .             CM,CN,CO,CP,CQ,CR,CS,
     .             CA1,CB1,CC1,CD1,CE1,CF1,CG1,CH1,CI1,CJ1,
     .             C1M,C2M,C3M,C4M,C5M,C6M,C7M,C8M,C9M,C10M
C
      COMMON/NAMEA/CA,CB,CC,CD,CE,CF,CG,CH,CI,CJ,CK,CL,
     .             CM,CN,CO,CP,CQ,CR,CS,
     .             CA1,CB1,CC1,CD1,CE1,CF1,CG1,CH1,CI1,CJ1,
     .             C1M,C2M,C3M,C4M,C5M,C6M,C7M,C8M,C9M,C10M
C
      INTEGER*4      INUS1,ILUS1
      COMMON/USEFULM/INUS1,ILUS1
C
C**** LOGUN
C
C     To add units, see: chek01.f, froapp.f, frodes.f, froeli.f,
C                        fronts.f, frotap.f, inpcek.f, jacbou.f,
C                        jacobs.f, pcgite.f, renum0.f,
C                        skydes.f, skygra.f, skyite.f,
C                        solver.f, update.f
C
      INTEGER*4    LUDTS,LUSOL,LUFRO,LUFRH,LUDAT,LUPRI,LURES,
     .             LUSO2,LUFR2,LUPOS,LURST,LUBFG,LUPIP,LUPAN,
     .             LUGEO,LUSET,LUMAT,LUINI,LULOA,LUFIX,LUIN1,
     .             LUTUN,LUCON,LUACT,LUFAN,
     .             LUCU1,LUCU2,LUCU3,LUCU4,LUCU5,LUCU6,LUCU7,
     .             LUCU8,LUCU9,LUC10
C
      COMMON/LOGUN/LUDTS,LUSOL,LUFRO,LUFRH,LUDAT,LUPRI,LURES,
     .             LUSO2,LUFR2,LUPOS,LURST,LUBFG,LUPIP,LUPAN,
     .             LUGEO,LUSET,LUMAT,LUINI,LULOA,LUFIX,LUIN1,
     .             LUTUN,LUCON,LUACT,LUFAN,
     .             LUCU1,LUCU2,LUCU3,LUCU4,LUCU5,LUCU6,LUCU7,
     .             LUCU8,LUCU9,LUC10
C--------------------------------------------------------------- PROBLEM
C**** AUGMENTED
C
      INTEGER*4        IAUGM,ICONC
C
      COMMON/AUGMENTED/IAUGM,ICONC
C=============================================================== PROBLEM
