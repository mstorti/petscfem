C=============================================================== AUXLIAR
C**** INDWORA
      INTEGER*4      IAQUA,IASEM,IFIXY,IFORC,IFORD,
     .               IFORW,IGSMO,IIPCW,ILOSU,IQUAS,
     .               ISETM,ISOLV,ISTAR,ISTIF,ILOSY
C
      COMMON/INDWORA/IAQUA(10),IASEM(10),IFIXY(20),IFORC(50),IFORD(50),
     .               IFORW(10),IGSMO(50),IIPCW(10),ILOSU(50),IQUAS(10),
     .               ISETM(50),ISOLV(20),ISTAR(30),ISTIF(50),ILOSY(10)
C--------------------------------------------------------------- AUXLIAR
C**** MEMORYA
      INTEGER*4      IADPR,LPRIN,IADW1,LWOR1,IADSO,LSOLV,IRELE,
     .               LBYTS,IADDB,LDABA,IREL1,LBYTB,NBLIM,LSTIF
C
      COMMON/MEMORYA/IADPR,LPRIN,IADW1,LWOR1,IADSO,LSOLV,IRELE,
     .               LBYTS,IADDB,LDABA,IREL1,LBYTB,NBLIM,LSTIF
C--------------------------------------------------------------- AUXLIAR
C**** NEWINTA
      INTEGER*4      NEWBO,NEWLO,NEWST,NEWFU,NEWAC,
     .               IOFIX,IOLOA,            IOACT,IOINI
C
      COMMON/NEWINTA/NEWBO,NEWLO,NEWST,NEWFU,NEWAC,
     .               IOFIX,IOLOA,            IOACT,IOINI
C--------------------------------------------------------------- AUXLIAR
C**** PROCESA
      INTEGER*4      IELEM,LMATS,NDIML,NNODL,NRULE,
     .               NGAUL,NTYPE,NCRIT,NKOST,NSTRE,NSTRS,NNODS,
     .               NQUTR
C
      COMMON/PROCESA/IELEM,LMATS,NDIML,NNODL,NRULE,
     .               NGAUL,NTYPE,NCRIT,NKOST,NSTRE,NSTRS,NNODS,
     .               NQUTR
C--------------------------------------------------------------- AUXLIAR
C**** RESIDUA
      REAL*8         AALPH,GZERO,GCURN,AACOE,BBCOE
C
      COMMON/RESIDUA/AALPH,GZERO,GCURN,AACOE,BBCOE
C--------------------------------------------------------------- AUXLIAR
C**** RSTARTA
C**** RSTARTB
      INTEGER*4      INITI,IREST,ISAVE,ISKIP,KTIME,KSTEP,KTSTE,
     .               NWPOS,KPPCG,KPPCN
      REAL*8         TIMST,TLIMT
C
      COMMON/RSTARTA/INITI,IREST,ISAVE,ISKIP,KTIME,KSTEP,KTSTE,
     .               NWPOS,KPPCG,KPPCN
      COMMON/RSTARTB/TIMST,TLIMT
C--------------------------------------------------------------- AUXLIAR
C**** RUNTIMA
      REAL*8         CPUIN,CPUDA,CPUST,CPUSF,CPUAS,CPUSO,CPURE,
     .               CPURS,CPUOU,CPURN
C
      COMMON/RUNTIMA/CPUIN,CPUDA,CPUST,CPUSF,CPUAS,CPUSO,CPURE,
     .               CPURS,CPUOU,CPURN
C--------------------------------------------------------------- AUXLIAR
C**** SMIPTI
      INTEGER*4     IPLAS,IPLAO,IPLAN,IPLAM
      INTEGER*4     NNUIN,NNUNO
      INTEGER*4     KPLA1,KPLA2,KPLA3,KPLA4,KPLA5,KPLA6,KPLA7
      INTEGER*4     NNUPCM,NNUPTM
C
      COMMON/SMIPT1/IPLAS(50),IPLAO(50),IPLAN(50),IPLAM(50)
      COMMON/SMIPT2/NNUIN,NNUNO
      COMMON/SMIPT3/KPLA1,KPLA2,KPLA3,KPLA4,KPLA5,KPLA6,KPLA7
      COMMON/SMIPMICM/NNUPCM,NNUPTM
C
C--------------------------------------------------------------- AUXLIAR
C**** SMIVTI
      INTEGER*4     IVLAS,IVLAO,IVLAN,IVLAM
      INTEGER*4     KVLA1,KVLA2,KVLA3,KVLA4,KVLA5,KVLA6,KVLA7
C
      COMMON/SMIVT1/IVLAS(50),IVLAO(50),IVLAN(50),IVLAM(50)
      COMMON/SMIVT3/KVLA1,KVLA2,KVLA3,KVLA4,KVLA5,KVLA6,KVLA7
C
C--------------------------------------------------------------- AUXLIAR
C**** PROPER1
      INTEGER*4      NYOUN,NPOIS,NALPH,NCCER,
     .               IVIFL,NVISC,NEXPO,ISOTT,IKINE,
     .               NCCOE,NCCOB,NCCOQ,NKCOE,NKCOB,NKCOQ,
     .               NCCOEE,
     .              NPLATM,IDAMG,NDAMA,NFRAC,IPORO,NPCOE,NPCO1,NPCO2,
     .               NZABA,NZABB,NZABC,ISHRI,IFREN
C
      COMMON/PROPER1/NYOUN,NPOIS,NALPH,NCCER,
     .               IVIFL,NVISC,NEXPO,ISOTT,IKINE,
     .               NCCOE,NCCOB,NCCOQ,NKCOE,NKCOB,NKCOQ,
     .               NCCOEE,
     .              NPLATM,IDAMG,NDAMA,NFRAC,IPORO,NPCOE,NPCO1,NPCO2,
     .               NZABA,NZABB,NZABC,ISHRI,IFREN
C
C**** PROPER2
      INTEGER*4      NYOUNF,NYOUNA,NPOISF,NPOISA,NALPHF,NALPHA,
     .               NCCERF,NCCERA,NVISCF,NVISCA,NEXPOF,NEXPOA,
     .               NCCOEF,NCCOEA,NCCOBF,NCCOBA,NCCOQF,NCCOQA,
     .               NKCOEF,NKCOEA,NKCOBF,NKCOBA,NKCOQF,NKCOQA,
     .               NDAMAF,NDAMAA,NFRACF,NFRACA,NPCOEF,NPCOEA,
     .               NZABAF,NZABAA,NZABBF,NZABBA,NZABCF,NZABCA
C
      COMMON/PROPER2/NYOUNF,NYOUNA,NPOISF,NPOISA,NALPHF,NALPHA,
     .               NCCERF,NCCERA,NVISCF,NVISCA,NEXPOF,NEXPOA,
     .               NCCOEF,NCCOEA,NCCOBF,NCCOBA,NCCOQF,NCCOQA,
     .               NKCOEF,NKCOEA,NKCOBF,NKCOBA,NKCOQF,NKCOQA,
     .               NDAMAF,NDAMAA,NFRACF,NFRACA,NPCOEF,NPCOEA,
     .               NZABAF,NZABAA,NZABBF,NZABBA,NZABCF,NZABCA
C
C**** PROPER3
      REAL*8         VYOUN,VPOIS,VALPH,VCCER,
     .               VVISC,VEXPO,
     .               VCCOE,VCCOB,VCCOQ,VKCOE,VKCOB,VKCOQ,
     .               VCCOEE,
     .               TEREF,TEMPS,TEMPL,VPLATM,
     .               VDAMA,VFRAC,VPCOE,
     .               VPCO1,VPCO2,
     .               VZABA,VZABB,VZABC,VALPU
C
      COMMON/PROPER3/VYOUN(20,2),VPOIS(20,2),VALPH(20,2),VCCER(20,2),
     .               VVISC(20,2),VEXPO(20,2),
     .               VCCOE(20,2),VCCOB(20,2),VCCOQ(20,2),
     .               VCCOEE(20,4),
     .               VKCOE(20,2),VKCOB(20,2),VKCOQ(20,2),
     .               TEREF,TEMPS,TEMPL,VPLATM(5,20),
     .               VDAMA(20,2),VFRAC(20,2),VPCOE(20,2),
     .               VPCO1(20,2),VPCO2(20,2),
     .               VZABA(20,2),VZABB(20,2),VZABC(20,2),VALPU(20,2)
C
C**** PROPER4
      REAL*8         VYOUNF,VYOUNA,VPOISF,VPOISA,VALPHF,VALPHA,
     .               VCCERF,VCCERA,VVISCF,VVISCA,VEXPOF,VEXPOA,
     .               VCCOEF,VCCOEA,VCCOBF,VCCOBA,VCCOQF,VCCOQA,
     .               VKCOEF,VKCOEA,VKCOBF,VKCOBA,VKCOQF,VKCOQA,
     .               EXPAN,EXPAU,
     .               VDAMAF,VDAMAA,VFRACF,VFRACA,VPCOEF,VPCOEA,
     .               VZABAF,VZABAA,VZABBF,VZABBA,VZABCF,VZABCA,
     .               VALPUF,VALPUA
C
      COMMON/PROPER4/VYOUNF(20,2),VYOUNA(20,2),
     .              VPOISF(20,2),VPOISA(20,2),VALPHF(20,2),VALPHA(20,2),
     .              VCCERF(20,2),VCCERA(20,2),VVISCF(20,2),VVISCA(20,2),
     .              VEXPOF(20,2),VEXPOA(20,2),
     .              VCCOEF(20,2),VCCOEA(20,2),VCCOBF(20,2),VCCOBA(20,2),
     .              VCCOQF(20,2),VCCOQA(20,2),
     .              VKCOEF(20,2),VKCOEA(20,2),VKCOBF(20,2),VKCOBA(20,2),
     .              VKCOQF(20,2),VKCOQA(20,2),
     .              EXPAN(5,20,3),EXPAU(5,20),
     .              VDAMAF(20,2),VDAMAA(20,2),VFRACF(20,2),VFRACA(20,2),
     .              VPCOEF(20,2),VPCOEA(20,2),
     .              VZABAF(20,2),VZABAA(20,2),VZABBF(20,2),VZABBA(20,2),
     .              VZABCF(20,2),VZABCA(20,2),VALPUF(20,2),VALPUA(20,2)
C
C**** PROPER5
      INTEGER*4      NCOE1,NCOE2,NCOE3,NCOE4,NCOE5,NCOE6,
     .               NCOE7,NCOE8,NCOE9
C
      COMMON/PROPER5/NCOE1,NCOE2,NCOE3,NCOE4,NCOE5,NCOE6,
     .               NCOE7,NCOE8,NCOE9
C
C**** PROPER6
      REAL*8         VCOE1,VCOE2,VCOE3,
     .               VCOE4,VCOE5,VCOE6,
     .               VCOE7,VCOE8,VCOE9
C
      COMMON/PROPER6/VCOE1(20,2),VCOE2(20,2),VCOE3(20,2),
     .               VCOE4(20,2),VCOE5(20,2),VCOE6(20,2),
     .               VCOE7(20,2),VCOE8(20,2),VCOE9(20,2)
C
C**** PROPER7
      REAL*8         GAMMAM,GAMMAP
C
      COMMON/PROPER7/GAMMAM,GAMMAP
C
C--------------------------------------------------------------- AUXLIAR
C**** FATIGUE PROBLEMS
      INTEGER*4       IFATI,IFATM
C
      COMMON/FATIGUE1/IFATI,IFATM
C
      REAL*8          VFATI
C
      COMMON/FATIGUE2/VFATI(20)
C--------------------------------------------------------------- AUXLIAR
C**** COSOLTRU
      INTEGER*4      NCETA
      COMMON/TRUCHIM/NCETA
C=============================================================== AUXLIAR
