      subroutine initf(ndim,nel,npg,nen,
     *     nhist_a,nstr1_a,nkost_a,nprop_a,
     *     nstre_a,nstrs_a,nfpch_a,props,
     *     E,dnu,ccero,dt,kdyna_a,inwt,ndof)
      implicit real*8(a-h,o-z)
c
c**** additional parameters
c
      include 'addi_om.f'
c
c
c**** coupling variables
c
      include 'nuec_om.f'
c
c**** mechanical variables
c
      include 'prob_om.f'
      include 'inte_om.f'
      include 'auxl_om.f'
      include 'inpo_om.f'
c
c --------------- COMMON ON FRIN30 ---------------------
      COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
      COMMON/JACOBSA/IEROR,KEROR

c --------------- COMMON ON STIF30 ---------------------
      COMMON/REHMAT/KFLAG,IFLAG
c
      dimension props(nprop_a)
c
c     Variables with `_a' are used to define variables in common
      ndime = ndim
      nnodl = nel
      ngaul = npg
      nstr1 = nstr1_a
      nkost = nkost_a
      nhist = nhist_a
      nevab = nen
      nprop = nprop_a
      nstre = nstre_a
      nstrs = nstrs_a
      nfpch = nfpch_a
      kdyna = kdyna_a

c =========================      =========================      =============
c     ADDI_OM.F
      nmemom = 0

c->   C=============================================================== NUEVOSC
c->   C**** INTRD0A
c->   C
c->         INTEGER*4      ISSTEP,NSSTEP,ISSTEPT,NSSTEPT
c->   C
c->         COMMON/INTRD0A/ISSTEP,NSSTEP,ISSTEPT,NSSTEPT
c->   C--------------------------------------------------------------- NUEVOSC
c->   C**** OUTPUT & OUTPUTT
c->   C
c->         INTEGER*4      NCKGLO
c->   C
c->         COMMON/OUTPUTA/NCKGLO
c->   C--------------------------------------------------------------- NUEVOSC
c->   C**** OUTPUT & OUTPUTT
c->   C
c->         INTEGER*4      ICOSTA,ITERME,ISTAGG,INICOU
c->         INTEGER*4      ITERMG,ITERMP,ITERMD
c->   C
c->         COMMON/COSTAGG/ICOSTA,ITERME,ISTAGG,INICOU
c->         COMMON/COSTAGH/ITERMG,ITERMP,ITERMD

      iterme = 0

c->   C--------------------------------------------------------------- NUEVOSC
c->   C**** ELEMENT
c->   C
c->         INTEGER*4      NPOINC,NELEMC
c->   C
c->         COMMON/COSTAEL/NPOINC,NELEMC
c->   C--------------------------------------------------------------- NUEVOSC
c->   C**** ELEMENN
c->   C
c->         REAL*8         COUFAC,CENKEL
c->         INTEGER*4      IITERC,INTERC
c->         INTEGER*4      NITERC,NFCOUC,ICONVC
c->         INTEGER*4      NFPCH
c->   C
c->         COMMON/COSTAE1/COUFAC,CENKEL
c->         COMMON/COSTAE2/IITERC,INTERC
c->         COMMON/COSTAE3/NITERC,NFCOUC,ICONVC
c->         COMMON/COSTAE4/NFPCH
c->   C--------------------------------------------------------------- NUEVOSC
c->   C**** NAMECA
c->   C
c->         CHARACTER*80  CCOA,CCOB
c->   C
c->         COMMON/NAMECA/CCOA,CCOB
c->   C
c->         INTEGER*4      INUS1C,ILUS1C
c->         COMMON/USEFULC/INUS1C,ILUS1C
c->   C
c->   C**** LOGUNC
c->   C
c->   C     To add units, see:                        ! to be confirmed
c->   C
c->         INTEGER*4     LUDATC,LURESC
c->   C
c->         COMMON/LOGUNC/LUDATC,LURESC
C=============================================================== NUEVOSC



c->   C=============================================================== AUXLIAR
c->   C**** INDWORA
c->         INTEGER*4      IAQUA,IASEM,IFIXY,IFORC,IFORD,
c->        .               IFORW,IGSMO,IIPCW,ILOSU,IQUAS,
c->        .               ISETM,ISOLV,ISTAR,ISTIF,ILOSY
c->   C
c->         COMMON/INDWORA/IAQUA(10),IASEM(10),IFIXY(20),IFORC(50),IFORD(50),
c->        .               IFORW(10),IGSMO(50),IIPCW(10),ILOSU(50),IQUAS(10),
c->        .               ISETM(50),ISOLV(20),ISTAR(30),ISTIF(50),ILOSY(10)
c->   C--------------------------------------------------------------- AUXLIAR
c->   C**** MEMORYA
c->         INTEGER*4      IADPR,LPRIN,IADW1,LWOR1,IADSO,LSOLV,IRELE,
c->        .               LBYTS,IADDB,LDABA,IREL1,LBYTB,NBLIM,LSTIF
c->   C
c->         COMMON/MEMORYA/IADPR,LPRIN,IADW1,LWOR1,IADSO,LSOLV,IRELE,
c->        .               LBYTS,IADDB,LDABA,IREL1,LBYTB,NBLIM,LSTIF
c->   C--------------------------------------------------------------- AUXLIAR
c->   C**** NEWINTA
c->         INTEGER*4      NEWBO,NEWLO,NEWST,NEWFU,NEWAC,
c->        .               IOFIX,IOLOA,            IOACT,IOINI
c->   C
c->         COMMON/NEWINTA/NEWBO,NEWLO,NEWST,NEWFU,NEWAC,
c->        .               IOFIX,IOLOA,            IOACT,IOINI
c->   C--------------------------------------------------------------- AUXLIAR
c->   C**** PROCESA
c->         INTEGER*4      IELEM,LMATS,NDIML,NNODL,NRULE,
c->        .               NGAUL,NTYPE,NCRIT,NKOST,NSTRE,NSTRS,NNODS,
c->        .               NQUTR
c->   C
c->         COMMON/PROCESA/IELEM,LMATS,NDIML,NNODL,NRULE,
c->        .               NGAUL,NTYPE,NCRIT,NKOST,NSTRE,NSTRS,NNODS,
c->        .               NQUTR

      ntype = 2                 !  plain strain
      ncrit = int(props(36))
      nnods = nnodl

c->   C--------------------------------------------------------------- AUXLIAR
c->   C**** RESIDUA
c->         REAL*8         AALPH,GZERO,GCURN,AACOE,BBCOE
c->   C
c->         COMMON/RESIDUA/AALPH,GZERO,GCURN,AACOE,BBCOE
c->   C--------------------------------------------------------------- AUXLIAR
c->   C**** RSTARTA
c->   C**** RSTARTB
c->         INTEGER*4      INITI,IREST,ISAVE,ISKIP,KTIME,KSTEP,KTSTE,
c->        .               NWPOS,KPPCG,KPPCN
c->         REAL*8         TIMST,TLIMT
c->   C
c->         COMMON/RSTARTA/INITI,IREST,ISAVE,ISKIP,KTIME,KSTEP,KTSTE,
c->        .               NWPOS,KPPCG,KPPCN
c->         COMMON/RSTARTB/TIMST,TLIMT
c->   C--------------------------------------------------------------- AUXLIAR
c->   C**** RUNTIMA
c->         REAL*8         CPUIN,CPUDA,CPUST,CPUSF,CPUAS,CPUSO,CPURE,
c->        .               CPURS,CPUOU,CPURN
c->   C
c->         COMMON/RUNTIMA/CPUIN,CPUDA,CPUST,CPUSF,CPUAS,CPUSO,CPURE,
c->        .               CPURS,CPUOU,CPURN
c->   C--------------------------------------------------------------- AUXLIAR
c->   C**** SMIPTI
c->         INTEGER*4     IPLAS,IPLAO,IPLAN,IPLAM
c->         INTEGER*4     NNUIN,NNUNO
c->         INTEGER*4     KPLA1,KPLA2,KPLA3,KPLA4,KPLA5,KPLA6,KPLA7
c->         INTEGER*4     NNUPCM,NNUPTM
c->   C
c->         COMMON/SMIPT1/IPLAS(50),IPLAO(50),IPLAN(50),IPLAM(50)
c->         COMMON/SMIPT2/NNUIN,NNUNO
c->         COMMON/SMIPT3/KPLA1,KPLA2,KPLA3,KPLA4,KPLA5,KPLA6,KPLA7
c->         COMMON/SMIPMICM/NNUPCM,NNUPTM
c->   C
      
      kpla1 = 1
      kpla2 = 0
      kpla3 = 0
      kpla4 = 0
      kpla5 = 0
      kpla6 = 0
      kpla7 = 0

      NBASE=1                                     ! eff. plast. strain
      IF(KPLA1.EQ.1) NBASE=NBASE+1                ! isotropic hardening
      IF(KPLA2.EQ.1) NBASE=NBASE+NSTR1            ! kinematic hardening
C
      NNDAM=0
c     IF(KPLA3.EQ.1) NNDAM=1                      ! damage (standard)
      IF(KPLA3.EQ.1) NNDAM=2                      ! damage (concrete)
C
      IF(KPLA4.EQ.1) NBASE=NBASE+2                ! porosity+total hard.
C
      NSHRI=0
      IF(KPLA5.EQ.1) NSHRI=1                      ! shrinkage
C
      NFATI=0
      IF(KPLA6.EQ.1) NFATI=8                      ! fatigue
C
      IF(KPLA7.EQ.1) NBASE=NBASE+1                ! dot e^p
C
      NMECH=0
C
C**** PLASTIC MODELS
C
      IPLAS(1)=1                       ! STRAP(NSTR1)
      IPLAS(2)=IPLAS(1)+NSTR1          ! DMTEP(NKOST)
      IPLAS(3)=IPLAS(2)+NKOST          ! EBASE(NBASE)
      IPLAS(4)=IPLAS(3)+NBASE          ! DBASE(NNDAM)
      IPLAS(5)=IPLAS(4)+NNDAM          ! SHRIN
      IPLAS(6)=IPLAS(5)+NSHRI          ! COUTD
      IPLAS(7)=IPLAS(6)+NMECH          ! FATIG(NFATI)
      NPLAS   =IPLAS(7)+NFATI
      IF(NPLAS.GT.NHIST)
     . CALL RUNEND('ERROR: NPLAS GT NHIST')

c->   C--------------------------------------------------------------- AUXLIAR
c->   C**** SMIVTI
c->         INTEGER*4     IVLAS,IVLAO,IVLAN,IVLAM
c->         INTEGER*4     KVLA1,KVLA2,KVLA3,KVLA4,KVLA5,KVLA6,KVLA7
c->   C
c->         COMMON/SMIVT1/IVLAS(50),IVLAO(50),IVLAN(50),IVLAM(50)
c->         COMMON/SMIVT3/KVLA1,KVLA2,KVLA3,KVLA4,KVLA5,KVLA6,KVLA7
c->   C
c->   C--------------------------------------------------------------- AUXLIAR
c->   C**** PROPER1
c->         INTEGER*4      NYOUN,NPOIS,NALPH,NCCER,
c->        .               IVIFL,NVISC,NEXPO,ISOTT,IKINE,
c->        .               NCCOE,NCCOB,NCCOQ,NKCOE,NKCOB,NKCOQ,
c->        .               NCCOEE,
c->        .              NPLATM,IDAMG,NDAMA,NFRAC,IPORO,NPCOE,NPCO1,NPCO2,
c->        .               NZABA,NZABB,NZABC,ISHRI,IFREN
c->   C
c->         COMMON/PROPER1/NYOUN,NPOIS,NALPH,NCCER,
c->        .               IVIFL,NVISC,NEXPO,ISOTT,IKINE,
c->        .               NCCOE,NCCOB,NCCOQ,NKCOE,NKCOB,NKCOQ,
c->        .               NCCOEE,
c->        .              NPLATM,IDAMG,NDAMA,NFRAC,IPORO,NPCOE,NPCO1,NPCO2,
c->        .               NZABA,NZABB,NZABC,ISHRI,IFREN
c->   C
c->   C**** PROPER2
c->         INTEGER*4      NYOUNF,NYOUNA,NPOISF,NPOISA,NALPHF,NALPHA,
c->        .               NCCERF,NCCERA,NVISCF,NVISCA,NEXPOF,NEXPOA,
c->        .               NCCOEF,NCCOEA,NCCOBF,NCCOBA,NCCOQF,NCCOQA,
c->        .               NKCOEF,NKCOEA,NKCOBF,NKCOBA,NKCOQF,NKCOQA,
c->        .               NDAMAF,NDAMAA,NFRACF,NFRACA,NPCOEF,NPCOEA,
c->        .               NZABAF,NZABAA,NZABBF,NZABBA,NZABCF,NZABCA
c->   C
c->         COMMON/PROPER2/NYOUNF,NYOUNA,NPOISF,NPOISA,NALPHF,NALPHA,
c->        .               NCCERF,NCCERA,NVISCF,NVISCA,NEXPOF,NEXPOA,
c->        .               NCCOEF,NCCOEA,NCCOBF,NCCOBA,NCCOQF,NCCOQA,
c->        .               NKCOEF,NKCOEA,NKCOBF,NKCOBA,NKCOQF,NKCOQA,
c->        .               NDAMAF,NDAMAA,NFRACF,NFRACA,NPCOEF,NPCOEA,
c->        .               NZABAF,NZABAA,NZABBF,NZABBA,NZABCF,NZABCA
c->   C
c->   C**** PROPER3
c->         REAL*8         VYOUN,VPOIS,VALPH,VCCER,
c->        .               VVISC,VEXPO,
c->        .               VCCOE,VCCOB,VCCOQ,VKCOE,VKCOB,VKCOQ,
c->        .               VCCOEE,
c->        .               TEREF,TEMPS,TEMPL,VPLATM,
c->        .               VDAMA,VFRAC,VPCOE,
c->        .               VPCO1,VPCO2,
c->        .               VZABA,VZABB,VZABC,VALPU

      vyoun(1,1) = E
      vpois(1,1) = dnu
      vccer(1,1) = ccero

c->   C
c->         COMMON/PROPER3/VYOUN(20,2),VPOIS(20,2),VALPH(20,2),VCCER(20,2),
c->        .               VVISC(20,2),VEXPO(20,2),
c->        .               VCCOE(20,2),VCCOB(20,2),VCCOQ(20,2),
c->        .               VCCOEE(20,4),
c->        .               VKCOE(20,2),VKCOB(20,2),VKCOQ(20,2),
c->        .               TEREF,TEMPS,TEMPL,VPLATM(5,20),
c->        .               VDAMA(20,2),VFRAC(20,2),VPCOE(20,2),
c->        .               VPCO1(20,2),VPCO2(20,2),
c->        .               VZABA(20,2),VZABB(20,2),VZABC(20,2),VALPU(20,2)
c->   C
c->   C**** PROPER4
c->         REAL*8         VYOUNF,VYOUNA,VPOISF,VPOISA,VALPHF,VALPHA,
c->        .               VCCERF,VCCERA,VVISCF,VVISCA,VEXPOF,VEXPOA,
c->        .               VCCOEF,VCCOEA,VCCOBF,VCCOBA,VCCOQF,VCCOQA,
c->        .               VKCOEF,VKCOEA,VKCOBF,VKCOBA,VKCOQF,VKCOQA,
c->        .               EXPAN,EXPAU,
c->        .               VDAMAF,VDAMAA,VFRACF,VFRACA,VPCOEF,VPCOEA,
c->        .               VZABAF,VZABAA,VZABBF,VZABBA,VZABCF,VZABCA,
c->        .               VALPUF,VALPUA
c->   C
c->         COMMON/PROPER4/VYOUNF(20,2),VYOUNA(20,2),
c->        .              VPOISF(20,2),VPOISA(20,2),VALPHF(20,2),VALPHA(20,2),
c->        .              VCCERF(20,2),VCCERA(20,2),VVISCF(20,2),VVISCA(20,2),
c->        .              VEXPOF(20,2),VEXPOA(20,2),
c->        .              VCCOEF(20,2),VCCOEA(20,2),VCCOBF(20,2),VCCOBA(20,2),
c->        .              VCCOQF(20,2),VCCOQA(20,2),
c->        .              VKCOEF(20,2),VKCOEA(20,2),VKCOBF(20,2),VKCOBA(20,2),
c->        .              VKCOQF(20,2),VKCOQA(20,2),
c->        .              EXPAN(5,20,3),EXPAU(5,20),
c->        .              VDAMAF(20,2),VDAMAA(20,2),VFRACF(20,2),VFRACA(20,2),
c->        .              VPCOEF(20,2),VPCOEA(20,2),
c->        .              VZABAF(20,2),VZABAA(20,2),VZABBF(20,2),VZABBA(20,2),
c->        .              VZABCF(20,2),VZABCA(20,2),VALPUF(20,2),VALPUA(20,2)
c->   C
c->   C**** PROPER5
c->         INTEGER*4      NCOE1,NCOE2,NCOE3,NCOE4,NCOE5,NCOE6,
c->        .               NCOE7,NCOE8,NCOE9
c->   C
c->         COMMON/PROPER5/NCOE1,NCOE2,NCOE3,NCOE4,NCOE5,NCOE6,
c->        .               NCOE7,NCOE8,NCOE9
c->   C
c->   C**** PROPER6
c->         REAL*8         VCOE1,VCOE2,VCOE3,
c->        .               VCOE4,VCOE5,VCOE6,
c->        .               VCOE7,VCOE8,VCOE9
c->   C
c->         COMMON/PROPER6/VCOE1(20,2),VCOE2(20,2),VCOE3(20,2),
c->        .               VCOE4(20,2),VCOE5(20,2),VCOE6(20,2),
c->        .               VCOE7(20,2),VCOE8(20,2),VCOE9(20,2)
c->   C
c->   C**** PROPER7
c->         REAL*8         GAMMAM,GAMMAP
c->   C
c->         COMMON/PROPER7/GAMMAM,GAMMAP
c->   C
c->   C--------------------------------------------------------------- AUXLIAR
c->   C**** FATIGUE PROBLEMS
c->         INTEGER*4       IFATI,IFATM
c->   C
c->         COMMON/FATIGUE1/IFATI,IFATM
c->   C
c->         REAL*8          VFATI
c->   C
c->         COMMON/FATIGUE2/VFATI(20)
c->   C--------------------------------------------------------------- AUXLIAR
c->   C**** COSOLTRU
c->         INTEGER*4      NCETA
c->         COMMON/TRUCHIM/NCETA
C=============================================================== AUXLIAR

c$$$C
c$$$        COMMON/LISTENA/NNWOR,NNPAR
c$$$        COMMON/LISTENB/PARAM,DPARA
c$$$        COMMON/LISTEND/WORDS,DWORD
c$$$C--------------------------------------------------------------- INP_OUT
c$$$C**** PLOTERA:
c$$$C**** PLOTERB:
c$$$        INTEGER*4   MMCUR,MSPLO          ! Max curves & Size of MPLOT
c$$$        PARAMETER   (MMCUR=40,MSPLO=240) ! MSPLO=MMCUR*6
c$$$C
c$$$        INTEGER*4     NCOLD,NCURV,NPONT(MMCUR),MPLOT(MMCUR,2,3)
c$$$        CHARACTER*8   FORMA
c$$$C
c$$$        INTEGER*4     IPRCO
c$$$C
c$$$        COMMON/PLOTERA/NCOLD,NCURV,NPONT,MPLOT
c$$$        COMMON/PLOTERB/FORMA
c$$$        COMMON/PLOTERC/IPRCO
c$$$C--------------------------------------------------------------- INP_OUT
c$$$C**** PRIOUTA:
c$$$        INTEGER*4      KFEMV,KPRI0,KPRI1,KPRI2,KPRI3,KPRI4,KPRI5,KPRI6,
c$$$     .                             KPRI7,KPRI8,KPRI9,
c$$$     .                 KPRI10,KPRI11,KPRI12,KPRI13,KPRI15,KPRI16,KPRI17,
c$$$     .                 KPRI18,KPRI19,KPRI20
c$$$C
c$$$        COMMON/PRIOUTA/KFEMV,KPRI0,KPRI1,KPRI2,KPRI3,KPRI4,KPRI5,KPRI6,
c$$$     .                             KPRI7,KPRI8,KPRI9,
c$$$     .                 KPRI10,KPRI11,KPRI12,KPRI13,KPRI15,KPRI16,KPRI17,
c$$$     .                 KPRI18,KPRI19,KPRI20
c$$$C=============================================================== INP_OUT

c$$$C=============================================================== INTERVL
c$$$C**** CONSTNA
c$$$      REAL*8         DITER,DTIME,GRAVY,GVECT,STIFI,
c$$$     .               TALFA,TBETA,TGAMA,TDELT,TOLER,XTIME,WLUMP
c$$$C
c$$$      COMMON/CONSTNA/DITER,DTIME,GRAVY,GVECT(3),STIFI,
c$$$     .               TALFA,TBETA,TGAMA,TDELT,TOLER,XTIME,WLUMP

      dtime = dt
      iiter = inwt + 1

c$$$C--------------------------------------------------------------- INTERVL
c$$$C**** CURRENA
c$$$C**** CURRENB
c$$$      INTEGER*4      ITIME,ISTEP,IITER,KDTIM,KRESL,KSTIF,KUNLD,
c$$$     .               NCHEK
c$$$      REAL*8         ARCLN,FACTO,FACPR,PITER,PVALU,STICU,
c$$$     .               STIIN,TFACT,TTIME
c$$$C
c$$$      COMMON/CURRENA/ITIME,ISTEP,IITER,KDTIM,KRESL,KSTIF,KUNLD,
c$$$     .               NCHEK
c$$$      COMMON/CURRENB/ARCLN,FACTO,FACPR,PITER,PVALU,STICU,
c$$$     .               STIIN,TFACT,TTIME
c$$$C--------------------------------------------------------------- INTERVL
c$$$C**** INTERVA
c$$$      INTEGER*4      IMPOS,KALGO,KARCL,KINTE,KOPTI,KCONV,KSAVE,
c$$$     .               LACCE,LAUTO,LINES,MITER,NALGO,NBACK,NCDIS,
c$$$     .               NSTEP,NOUTP

      nalgo = 1

c$$$C
c$$$      COMMON/INTERVA/IMPOS,KALGO,KARCL,KINTE,KOPTI,KCONV,KSAVE,
c$$$     .               LACCE,LAUTO,LINES,MITER,NALGO,NBACK,NCDIS,
c$$$     .               NSTEP,NOUTP(50)
c$$$C--------------------------------------------------------------- INTERVL
c$$$C**** INTERVA
c$$$      INTEGER*4       MKONT,NFORZ,MSUBP
c$$$      REAL*8          TOPLA,ALFAP,EXCTP
c$$$C
c$$$      COMMON/RETURNPA/MKONT,NFORZ,MSUBP
c$$$      COMMON/RETURNPB/TOPLA,ALFAP,EXCTP
c$$$C=============================================================== INTERVL

c->   C=============================================================== PROBLEM
c->   C**** DATBASA
c->   C
c->         INTEGER*4      IDATP,IDATC,
c->        .               LENRC,NLENC,NLENP,NRECC,NRECG,NRECP,NWORP
c->   C
c->         COMMON/DATBASA/IDATP(12,5),IDATC(12),
c->        .               LENRC,NLENC,NLENP,NRECC,NRECG,NRECP,NWORP
c->   C--------------------------------------------------------------- PROBLEM
c->   C**** DIMEMNA
c->   C
c->         INTEGER*4      NDIME,NELEM,NFUNC,NGRUP,NHLOD,NHIST,NMATS,
c->        .               NPOIN,NPREL,NPROP,NSTR1,NTOTV,NSUBF,NTOTG,NPOIC
c->   C
c->         COMMON/DIMENNA/NDIME,NELEM,NFUNC,NGRUP,NHLOD,NHIST,NMATS,
c->        .               NPOIN,NPREL,NPROP,NSTR1,NTOTV,NSUBF,NTOTG,NPOIC
c->   C--------------------------------------------------------------- PROBLEM
c->   C**** ELDATAA
c->   C
c->         INTEGER*4      IDATA,IPREV,ISTAT,IMATX
c->   C
c->         COMMON/ELDATAA/IDATA(15),IPREV(3),ISTAT(4),IMATX(7)
c->   C--------------------------------------------------------------- PROBLEM
c->   C**** ELEMNTA
c->   C
c->         INTEGER*4      NDOFC,NDOFN,NEVAB,NEVAC,NGAUS,NKOND,NKOVA,NMOVA,
c->        .               NNODE,NDATA,NPREV,NSTAT,NMATX,NNODC
c->   C
c->         COMMON/ELEMNTA/NDOFC,NDOFN,NEVAB,NEVAC,NGAUS,NKOND,NKOVA,NMOVA,
c->        .               NNODE,NDATA,NPREV,NSTAT,NMATX,NNODC

      ndofc = ndof
      ndofn = ndof
      ngaus = ngaul
      nkova = nevab * nevab 
      nnode = nnodl

c->   C--------------------------------------------------------------- PROBLEM
c->   C**** HOURGLA
c->   C**** HOURGLB
c->   C
c->         INTEGER*4     NHOUR,KELAS
c->         REAL*8        HPARA
c->   C
c->         COMMON/HOURGLA/NHOUR,KELAS
c->         COMMON/HOURGLB/HPARA
c->   C--------------------------------------------------------------- PROBLEM
c->   C**** PROBLMA
c->   C
c->         INTEGER*4      KDYNA,KPORE,KPOST,KPROB,KSGAU,KSMUS,KTEMP,
c->        .               LARGE,NSKEW,LSKEW,NDISR,NDISO
c->   C
c->         COMMON/PROBLMA/KDYNA,KPORE,KPOST,KPROB,KSGAU,KSMUS,KTEMP,
c->        .               LARGE,NSKEW,LSKEW,NDISR,NDISO

      large = 0
      kprob = 1

c->   C--------------------------------------------------------------- PROBLEM
c->   C**** SOLVERA
c->   C**** SOLVERB
c->   C
c->         INTEGER*4      KRENU,KSOLV,KSYMM,NWIDT,MITCG,NBUFA,NPRIR
c->         REAL*8         TOLCG,TOLC1
c->   C
c->         COMMON/SOLVERA/KRENU,KSOLV,KSYMM,NWIDT,MITCG,NBUFA,NPRIR
c->         COMMON/SOLVERB/TOLCG,TOLC1

      ksymm = 0

c->   C--------------------------------------------------------------- PROBLEM
c->   C**** TITLESA
c->   C
c->         CHARACTER*8    TITLE,SUBTI
c->   C
c->         COMMON/TITLESA/TITLE(8),SUBTI(8)
c->   C--------------------------------------------------------------- PROBLEM
c->   C**** ACTIVE ELEMENTS
c->   C
c->         INTEGER*4     NACTI
c->         COMMON/ACTIV1/NACTI
c->   C--------------------------------------------------------------- PROBLEM
c->   C**** NAMEA
c->   C
c->   C    To add names, see: frodes.f, fronts.f & skydes.f
c->   C
c->         CHARACTER*80 CA,CB,CC,CD,CE,CF,CG,CH,CI,CJ,CK,CL,
c->        .             CM,CN,CO,CP,CQ,CR,CS,
c->        .             CA1,CB1,CC1,CD1,CE1,CF1,CG1,CH1,CI1,CJ1,
c->        .             C1M,C2M,C3M,C4M,C5M,C6M,C7M,C8M,C9M,C10M
c->   C
c->         COMMON/NAMEA/CA,CB,CC,CD,CE,CF,CG,CH,CI,CJ,CK,CL,
c->        .             CM,CN,CO,CP,CQ,CR,CS,
c->        .             CA1,CB1,CC1,CD1,CE1,CF1,CG1,CH1,CI1,CJ1,
c->        .             C1M,C2M,C3M,C4M,C5M,C6M,C7M,C8M,C9M,C10M
c->   C
c->         INTEGER*4      INUS1,ILUS1
c->         COMMON/USEFULM/INUS1,ILUS1
c->   C
c->   C**** LOGUN
c->   C
c->   C     To add units, see: chek01.f, froapp.f, frodes.f, froeli.f,
c->   C                        fronts.f, frotap.f, inpcek.f, jacbou.f,
c->   C                        jacobs.f, pcgite.f, renum0.f,
c->   C                        skydes.f, skygra.f, skyite.f,
c->   C                        solver.f, update.f
c->   C
c->         INTEGER*4    LUDTS,LUSOL,LUFRO,LUFRH,LUDAT,LUPRI,LURES,
c->        .             LUSO2,LUFR2,LUPOS,LURST,LUBFG,LUPIP,LUPAN,
c->        .             LUGEO,LUSET,LUMAT,LUINI,LULOA,LUFIX,LUIN1,
c->        .             LUTUN,LUCON,LUACT,LUFAN,
c->        .             LUCU1,LUCU2,LUCU3,LUCU4,LUCU5,LUCU6,LUCU7,
c->        .             LUCU8,LUCU9,LUC10
c->   C
c->         COMMON/LOGUN/LUDTS,LUSOL,LUFRO,LUFRH,LUDAT,LUPRI,LURES,
c->        .             LUSO2,LUFR2,LUPOS,LURST,LUBFG,LUPIP,LUPAN,
c->        .             LUGEO,LUSET,LUMAT,LUINI,LULOA,LUFIX,LUIN1,
c->        .             LUTUN,LUCON,LUACT,LUFAN,
c->        .             LUCU1,LUCU2,LUCU3,LUCU4,LUCU5,LUCU6,LUCU7,
c->        .             LUCU8,LUCU9,LUC10
c->   C--------------------------------------------------------------- PROBLEM
c->   C**** AUGMENTED
c->   C
c->         INTEGER*4        IAUGM,ICONC
c->   C
c->         COMMON/AUGMENTED/IAUGM,ICONC
c->   C=============================================================== PROBLEM


c --------------- COMMON IN FRIN30 ---------------------
c->         COMMON/TUNING30/PPART,PPARB,PPARI,ESTAB
c->         COMMON/JACOBSA/IEROR,KEROR
      ppart = 0
      pparb = 0
      ppari = 0
      estab = 1.e7
      ieror = 0
      keror = 0

c      COMMON/REHMAT/KFLAG,IFLAG
      kflag = 0
      iflag = 0

      return
c
      end
