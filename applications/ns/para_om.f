C***********************************************************************
C
C***** MECHANICAL PARAMETERS
C
C***********************************************************************
C
C     DIFFERENCES WITH para_omt.f >> !*
C
C***********************************************************************
C
      INCLUDE 'gend_om.f'
C
      PARAMETER(
c    . MNODG=max(MNODE,MGAUS) )               ! SG
     . MNODG=(MNODE+MGAUS) )                  ! PC & linux
C
C**** ADDELM
C
      PARAMETER(
c    . MDATA1AM=max(MMEMA1M,MMEMA2M,MMEMOM,(MSMOM*MMEMA2M),MHOUR) ) ! SG
     . MDATA1AM=0 )                           ! PC & linux
C
      PARAMETER(
     . MDATA1=MDATA1AM,                       ! ELCOD(NDIME,NNODE)
     . MDATA2=MDATA1+MDIME*MNODE*MMEMA1M,     ! CARTD(NDIME,NNODE,NGAUS)
     . MDATA3=MDATA2+
     .             MDIME*MNODE*MGAUS*MMEMA2M, ! DVOLU(NGAUS)
     . MDATA4=MDATA3+MGAUS*MMEMA2M,           ! GPCOD(NDIME,NGAUS)
     . MDATA5=MDATA4+MDIME*MGAUS*MMEMA2M,     ! SHAPE(NNODE,NGAUS)
     . MDATA6=MDATA5+MNODE*MGAUS*MMEMA2M,     ! EPMTX(NKOST,NGAUS)
     . MDATA7=MDATA6+MKOST*MGAUS*MMEMOM,      ! RMAT1(NDIME,NDIME)
     . MDATA8=MDATA7+MDIME*MDIME*MMEMOM,      ! EMASS(NNODE,NNODE)
     . MDATA9=MDATA8+
     .             MNODE*MNODE*MSMOM*MMEMA2M, ! STIFH(NEVAB,NEVAB)
     . MDATA10=MDATA9+MEVAB*MEVAB*MHOUR,      !*BSBAR(NSTR1,NEVAB,NGAUS)
     . MDATA11=MDATA10+
     .             MSTR1*MEVAB*MGAUS*MMEMA2M, ! VNORI(NDIME,NGAUS)
     . MDATA  =MDATA11+MDIME*MGAUS*MMEMA2M )
C
      PARAMETER(
     . MPREV1=1*MMEMOM,                       ! STRA0(NSTR1,NGAUS)
     . MPREV2=MPREV1+MSTR1*MGAUS*MMEMOM,      ! STRS0(NSTR1,NGAUS)
     . MPREV3=MPREV2+MSTR1*MGAUS*MMEMOM,      ! TEMPC(4)
     . MPREV=MPREV3+4*MMEMOM )
C
      PARAMETER(
     . MSTAT1=1,                              ! ELDIS(NDOFC,NNODE)
     . MSTAT2=MSTAT1+MDOFC*MNODE*MMEMA5M,     ! EHIST(NHIST,NGAUS)
     . MSTAT3=MSTAT2+MHIST*MGAUS,             ! STRAN(NSTR1,NGAUS) 
     . MSTAT4=MSTAT3+MSTR1*MGAUS,             ! STRSG(NSTR1,NGAUS)
     . MSTAT=MSTAT4+MSTR1*MGAUS )
C
      PARAMETER(
c    . MPOREA=max(MPORE-1,0) )                ! SG
     . MPOREA=0 )                             ! PC & linux
C
      PARAMETER(
c    . MMATXA=max(MMEMA7M,MPOREA) )           ! SG
     . MMATXA=1 )                             ! PC & linux
C
      PARAMETER(
     . MMATX1=MMATXA,                         ! CSTIF(NEVAC,NEVAC)
     . MMATX2=MMATX1+MEVAC*MEVAC*MMEMA7M,     ! KSTIF(NKOVA)
     . MMATX3=MMATX2+MKOVA*MMEMA6M,           ! MSTIF(NKOVA)
     . MMATX4=MMATX3+MKOVA*MDYNA*MMEMA6M,     ! PSTIF(NKOND)
     . MMATX5=MMATX4+MKOND*MPOREA,            ! QSTIF(NKOND)
     . MMATX6=MMATX5+MKOND*MPOREA,            ! HSTIF(NEVAB*NNODE)
     . MMATX=MMATX6+MEVAB*MNODE*MPOREA )
C
C**** ADDPRI
C
      PARAMETER(
c    . MTER12M=max((MTERMEM+1),0) )           ! SG
     . MTER12M=2 )                            ! PC & linux
      PARAMETER(
c    . MTER2M=max(MTERMEM,0) )                ! SG
     . MTER2M=1 )                             ! PC & linux
      PARAMETER(
     . MTER1M=MTER12M/(1+MTER2M) )
C
      PARAMETER(                              ! simplified
c    . MITERCTM=max(MTERMEM,0) )              ! SG
     . MITERCTM=1 )                           ! PC & linux
C
      PARAMETER(
c    . MPOREB=max(MPORE,1) )                  ! SG
     . MPOREB=1 )                             ! PC & linux
      PARAMETER(
     . MPOREC=MPORE/MPOREB )
C
      PARAMETER(
c    . MSKEA=max(MSKEW,1) )                   ! SG
     . MSKEA=1 )                              ! PC & linux
      PARAMETER(
c    . MSKE1=max((MSKEW/MSKEA),0) )           ! SG
     . MSKE1=1 )                              ! PC & linux
C
      PARAMETER(
     . mprin1= 1,                             ! LNODS(NNODE,NELEM)
     . mprin2=mprin1+(MNODE*MELEM*MCHAL+4)/8, ! MATNO(NELEM)
     . mprin3=mprin2+(MELEM*MCHAL+4)/8,       ! PROEL(NPREL,NGRUP)
     . mprin4=mprin3+MPREL*MGRUP,             ! PROPS(NPROP,NMATS)
     . mprin5=mprin4+MPROP*MMATS,             ! COORD(NDIME,NPOIN)
     . mprin6=mprin5+MDIME*MPOIN*MMEMO1M,     ! HTLOD(NHLOD,NSUBF,NFUNC)
     . mprin7=mprin6+MHLOD*MSUBF*MFUNC,       ! IFFIX(NTOTV,2)
     . mprin8=mprin7+(MTOTV*MCHAL+4)/8*
     .               (1+MPOREB/2)*(1+MACTI),  ! PRESC(NTOTV,2)
     . mprin9=mprin8+MTOTV*2,                 ! RLOAD(NTOTV)
     . mprin10=mprin9+MTOTV,                  ! RLOAH(NTOTV,NFUNC)
     . mprin11=mprin10+MTOTV*MFUNC,           ! FICTO(NFUNC)
     . mprin12=mprin11+MFUNC )                ! TFICT(NFUNC)
      PARAMETER(
     . mprin13=mprin12+MFUNC,                 ! DISIT(NTOTV,2)
     . mprin14=mprin13+MTOTV*(1+MMEMO9M),     ! DISPR(NTOTV,NDISR=3)
     . mprin15=mprin14+MTOTV*(1+2*MDYNA),     ! DISTO(NTOTV,NDISO=5)
     . mprin16=mprin15+MTOTV*(1+4*MDYNA),     ! HEADS(NPOIN,4)
     . mprin17=mprin16+MPOIN*4*MPOREC,        ! REFOR(NTOTV,2)
     . mprin18=mprin17+MTOTV*2,               ! TLOAD(NTOTV,2)
     . mprin19=mprin18+MTOTV*2,               ! LPNTN(NPOIN)
     . mprin20=mprin19+
     .               (MPOIN*MCHAL+4)/8*MRENU, ! ELDAT(NDATA)
     . mprin21=mprin20+MDATA,                 ! ELPRE(NPREV)
     . mprin22=mprin21+MPREV,                 ! ELVAR(NSTAT)
     . mprin23=mprin22+MSTAT,                 ! ELMAT(NMATX)
     . mprin24=mprin23+MMATX )                ! TEMPN(NPOIN,2)
      PARAMETER(
     . mprin25=mprin24+MPOIN*2*MTER1M,        ! DTEMN(NPOIN)
     . mprin26=mprin25+MPOIN*MTER1M,          ! INFRI(NPOIN)
     . mprin27=mprin26+
     .               (MPOIN*MCHAL+4)/8*MSKE1, ! COFRI(NSKEW,NDIME,NDIME)
     . mprin28=mprin27+MSKEW*MDIME*MDIME,     ! PWORK(NPOIN,2)
     . mprin29=mprin28+
     .             MPOIN*(1+MITERCTM)*MTER2M, ! PREAS(NPOIN)
     . mprin30=mprin29+MPOIN*MTER1M,          ! TGAPS(NPOIN)
     . mprin31=mprin30+MPOIN*MTER1M,          ! VNORM(NTOTV)
     . mprin32=mprin31+MTOTV,                 ! FPCHA(NFPCH,NPOIN)
     . mprin33=mprin32+MFPCHM*MPOIN*MTER1M,   ! LACTI(NELEM)
     . MPRIN  =mprin33+MACTI*MELEM )
C
C**** RSSETP
C
      PARAMETER(
     .       i11 =MDATA,                      ! ELDAT
     .       i21 =MPREV,                      ! ELPRE
     .       i31 =MSTAT,                      ! ELVAR
     .       i41 =MPREV,                      ! ELPRE
     .       i51 =MSTAT,                      ! ELVAR
     .       i61 =MEVAC*MEVAC*MMEMA7M,        ! CSTIF
     .       i71 =MKOVA*MMEMA6M,              ! ESTIF
     .       i81 =MKOVA*MDYNA*MMEMA6M,        ! WSTIF
     .       i91 =MKOND*MPOREA,               ! PSTIF
     .       i101=MKOND*MPOREA,               ! QSTIF
     .       i111=MEVAB*MNODE*MPOREA,         ! HSTIF
     .       i121=MTOTV )                     ! DISTO
C
      PARAMETER(
     .       i23 =MDATA,
     .       i33 =(i23+MPREV),
     .       i43 =(i33+MSTAT),
     .       i53 =(i43+MPREV),
     .       i63 =(i53+MSTAT),
     .       i73 =(i63+MEVAC*MEVAC*MMEMA7M),
     .       i83 =(i73+MDYNA*i71),
     .       i93 =(i83+i81),
     .       i103=(i93+i91),
     .       i113=(i103+i101),
     .       i123=(i113+i111) )
C
C**** DATBAS (BUFFER)
C
C     DISTO is always IN-CORE for MBUFFER
C
      PARAMETER(
c    . MWORP =max(i23,i33,i43,i53,i63,i73,    ! SG
c    .            i83,i93,i103,i113,i123),
     . MWORP =i123,                           ! PC & linux
c    . index1=max(MDATA,MPREV,MSTAT,(MEVAC*MEVAC)*MMEMA7M,  ! SG
c    .            MKOVA*MMEMA6M,MKOND*MPOREA,(MEVAB*MNODE),MTOTV),
     . index1=   (MDATA+MPREV+MSTAT+(MEVAC*MEVAC)*MMEMA7M+  ! PC & linux
     .            MKOVA*MMEMA6M+MKOND*MPOREA+(MEVAB*MNODE)+MTOTV),
     . index2=   (MDATA+MPREV+MSTAT+MEVAC*MEVAC*MMEMA7M+
     .            MKOVA*MMEMA6M+MKOND*MPOREA+MEVAB*MNODE+MTOTV),
     . MBUFFER=MTOTV+(MELEM*MWORP+INDEX2+1+INDEX1)*MDISKDM )
C
C**** ADDWOR
C
      PARAMETER(
     . ISTA1 =1,                              ! COORD(NDIME,NPOIN)
     . ISTA2 =ISTA1+MDIME*MPOIN*MMEMA1M,      ! GPCOD(NDIME,NGAUS)
     . ISTA3 =ISTA2+MDIME*MGAUS,              ! CARTD(NDIME,NNODE,NGAUS)
     . ISTA4 =ISTA3+
     .             MDIME*MNODE*MGAUS*MMEMO2M, ! DVOLU(NGAUS)
     . ISTA5 =ISTA4+MGAUS*MMEMO2M,            ! GPCOD(NDIME,NGAUS)
     . ISTA6 =ISTA5+MDIME*MGAUS*MMEMO2M,      ! SHAPE(NNODE,NGAUS)
     . ISTA7 =ISTA6+MNODE*MGAUS*MMEMO2M,      ! BSBAR(NSTR1,NEVAB,NGAUS)
     . ISTA8 =ISTA7+
     .             MSTR1*MEVAB*MGAUS*MMEMO2M, ! DERIV(NDIME,NNODE,NGAUS)
     . ISTA9 =ISTA8+
     .             MDIME*MNODE*MGAUS*MMEMO2M, ! POSGP(NDIME,NGAUS)
     . ISTA10=ISTA9 +MDIME*MGAUS*MMEMO2M,     ! WEIGP(NGAUS)
     . ISTA11=ISTA10+MGAUS*MMEMO2M,           ! XJACM(NDIME,NDIME)
     . ISTA12=ISTA11+MDIME*MDIME*MMEMO2M,     ! RMAT2(NSTR1,NSTR1)
     . ISTA13=ISTA12+MSTR1*MSTR1*MMEMO2M,     ! CMEAN(3,NNODE)
     . ISTA14=ISTA13+3*MNODE*MMEMO2M,         ! SHAPR(NNODE,NGAUS)
     . ISTA15=ISTA14+MNODE*MGAUS*MMEMO2M,     ! WSTIR(NNODE,NNODE)
     . ISTA16=ISTA15+MNODE*MNODE*MMEMO2M )    ! VNORL(NDIME,NNODG)
      PARAMETER(
     . ISTA17=ISTA16+MDIME*MNODG*MMEMO2M,     ! ELCO1(NDIME,NNODE)
     . ISTA18=ISTA17+MDIME*MNODE*MMEMO2M,     ! EMASS(NNODE,NNODE)
     . ISTA19=ISTA18+MNODE*MNODE*MMEMO2M,     ! ELDI1(NDOFN,NNODE)
     . ISTA20=ISTA19+MDOFC*MNODE*MMEMO2M,     ! CARTS(NDIME,NNODE,NGAUS)
     . ISTA21=ISTA20+
     .             MDIME*MNODE*MGAUS*MMEMO2M, ! SHAPS(NNODE,NGAUS)
     . ISTA22=ISTA21+MNODE*MGAUS*MMEMO2M,     ! GPCOS(NDIME,NGAUS)
     . LIND  =ISTA22+MDIME*MGAUS*MMEMO2M )
C
      PARAMETER(
     . LREN =MRENU*MPOIN*(10*MNODE+3)/2 )
C
      PARAMETER(
     . IQUA1=1,                               ! VVECT(NTOTV)
     . IQUA2=IQUA1+MTOTV,                     ! WVECT(NTOTV)
     . LQUA =IQUA2+MTOTV )
C
      PARAMETER(
     . IFIX1=1,                               ! LEQNS(NEVAC)
     . IFIX2=IFIX1+(MEVAC*MCHAL+4)/8,         ! LNUEQ(NTOTV)
     . IFIX3=IFIX2+(MTOTV*MCHAL+4)/8,         ! LOCEL(NEVAC)
     . IFIX4=IFIX3+(MEVAC*MCHAL+4)/8,         ! LPONT(NTOTV)
     . IFIX5=IFIX4+(MTOTV*MCHAL+4)/8,         ! NACVA(NTOTV)
     . IFIX6=IFIX5+(MTOTV*MCHAL+4)/8,         ! NDEST(NEVAC)
     . IFIX7=IFIX6+(MEVAC*MCHAL+4)/8,         ! NDFRO(NELEM)
     . IFIX8=IFIX7+(MELEM*MCHAL+4)/8,         ! PRESC(NDOFC)
     . LFIX =IFIX8+MDOFC )
C
      PARAMETER(
     . IFOT1=1,                               ! BMSIG(NEVAB)
     . IFOT2=IFOT1+MEVAB,                     ! BMATX(NSTR1,NEVAB)
     . IFOT3=IFOT2+MSTR1*MEVAB,               ! DESIG(NSTR1) 
     . IFOT4=IFOT3+MSTR1,                     ! DMATX(NSTR1,NSTR1)
     . IFOT5=IFOT4+MSTR1*MSTR1,               ! DSTRA(NSTR1)
     . IFOT6=IFOT5+MSTR1,                     ! PRESG(NSTR1)
     . IFOT7=IFOT6+MSTR1,                     ! SGTOT(NSTR1)
     . IFOT8=IFOT7+MSTR1,                     ! SIGMA(NSTR1)
     . IFOT9=IFOT8+MSTR1,                     ! TSTRA(NSTR1)
     . IFOT10=IFOT9+MSTR1,                    ! XJACM(NDIME,NDIME)
     . IFOT11=IFOT10+MDIME*MDIME,             ! TENOD(NNODE)
     . IFOT12=IFOT11+MNODE,                   ! DTENO(NNODE)
     . IFOT13=IFOT12+MNODE,                   ! CMEAN(3,NNODE) not used
     . IFOT14=IFOT13+3*MNODE,                 ! PWOEL(NNODE)
     . IFOT15=IFOT14+MNODE,                   ! DISPL(NEVAB)
     . IFOT16=IFOT15+MEVAB )                  ! PREAL(NNODE)
      PARAMETER(
     . IFOT17=IFOT16+MNODE,                   ! TGAPL(NNODE)
     . IFOT18=IFOT17+MNODE,                   ! VNORL(NDIME,NNODG)
     . IFOT19=IFOT18+MDIME*MNODG,             ! ELCO1(NDIME,NNODE)
     . IFOT20=IFOT19+MDIME*MNODE,             ! TENOI(NNODE)
     . IFOT21=IFOT20+MNODE,                   ! CARTD(NDIME,NNODE,NGAUS)
     . IFOT22=IFOT21+MDIME*MNODE*MGAUS,       ! DVOLU(NGAUS)
     . IFOT23=IFOT22+MGAUS,                   ! GPCOD(NDIME,NGAUS)
     . IFOT24=IFOT23+MDIME*MGAUS,             ! SHAPE(NNODE,NGAUS)
     . IFOT25=IFOT24+MNODE*MGAUS,             ! BSBAR(NSTR1,NEVAB,NGAUS)
     . IFOT26=IFOT25+MSTR1*MEVAB*MGAUS,       ! DERIV(NDIME,NNODE,NGAUS)
     . IFOT27=IFOT26+MDIME*MNODE*MGAUS,       ! POSGP(NDIME,NGAUS)
     . IFOT28=IFOT27+MDIME*MGAUS,             ! WEIGP(NGAUS)
     . IFOT29=IFOT28+MGAUS,                   ! XJACM(NDIME,NDIME)
     . IFOT30=IFOT29+MDIME*MDIME,             ! RMAT2(NSTR1,NSTR1)
     . IFOT31=IFOT30+MSTR1*MSTR1,             ! CMEAN(3,NNODE)
     . IFOT32=IFOT31+3*MNODE )                ! SHAPR(NNODE,NGAUS)
      PARAMETER(
     . IFOT33=IFOT32+MNODE*MGAUS,             ! WSTIR(NNODE,NNODE)
     . IFOT34=IFOT33+MNODE*MNODE,             ! VNORL(NDIME,NNODG)
     . IFOT35=IFOT34+MDIME*MNODG,             ! EMASS(NNODE,NNODE)
     . IFOT36=IFOT35+MNODE*MNODE,             ! FPCHL(NFPCH,NNODE)
     . IFOT37=IFOT36+MFPCHM*MNODE*MTER1M,     ! ELDI1(NDOFN,NNODE)
     . IFOT38=IFOT37+MDOFC*MNODE,             ! DISIL(NEVAB)
     . IFOT39=IFOT38+MEVAB,                   ! CARTS(NDIME,NNODE,NGAUS)
     . IFOT40=IFOT39+MDIME*MNODE*MGAUS,       ! SHAPS(NNODE,NGAUS)
     . IFOT41=IFOT40+MNODE*MGAUS,             ! GPCOS(NDIME,NGAUS)
     . LFOT  =IFOT41+MDIME*MGAUS )
C
      PARAMETER(
     . IFOR1=1,                               ! ELELM(NEVAB)
     . IFOR2=IFOR1+MEVAB,                     ! ACELM(NEVAB)
     . IFOR3=IFOR2+MEVAB,                     ! VEELM(NEVAB)
     . IFOR4=IFOR3+MEVAB,                     ! ESTIF(NKOVA)
     . IFOR5=IFOR4+MKOVA,                     ! DSTIF(NKOVA)
     . IFOR6 =IFOR5 +MKOVA,                   ! CARTD(NDIME,NNODE,NGAUS)
     . IFOR7 =IFOR6 +MDIME*MNODE*MGAUS,       ! DVOLU(NGAUS)
     . IFOR8 =IFOR7 +MGAUS,                   ! GPCOD(NDIME,NGAUS)
     . IFOR9 =IFOR8 +MDIME*MGAUS,             ! SHAPE(NNODE,NGAUS)
     . IFOR10=IFOR9 +MNODE*MGAUS,             ! BSBAR(NSTR1,NEVAB,NGAUS)
     . IFOR11=IFOR10+MSTR1*MEVAB*MGAUS,       ! DERIV(NDIME,NNODE,NGAUS)
     . IFOR12=IFOR11+MDIME*MNODE*MGAUS,       ! POSGP(NDIME,NGAUS)
     . IFOR13=IFOR12+MDIME*MGAUS,             ! WEIGP(NGAUS)
     . IFOR14=IFOR13+MGAUS,                   ! XJACM(NDIME,NDIME)
     . IFOR15=IFOR14+MDIME*MDIME,             ! RMAT2(NSTR1,NSTR1)
     . IFOR16=IFOR15+MSTR1*MSTR1 )            ! CMEAN(3,NNODE)
      PARAMETER(
     . IFOR17=IFOR16+3*MNODE,                 ! SHAPR(NNODE,NGAUS)
     . IFOR18=IFOR17+MNODE*MGAUS,             ! WSTIR(NNODE,NNODE)
     . IFOR19=IFOR18+MNODE*MNODE,             ! VNORL(NDIME,NNODG)
     . IFOR20=IFOR19+MDIME*MNODG,             ! ELCO1(NDIME,NNODE)
     . IFOR21=IFOR20+MDIME*MNODE,             ! EMASS(NNODE,NNODE)
     . IFOR22=IFOR21+MNODE*MNODE,             ! ELDI1(NDOFN,NNODE)
     . IFOR23=IFOR22+MDOFC*MNODE,             ! CARTS(NDIME,NNODE,NGAUS)
     . IFOR24=IFOR23+MDIME*MNODE*MGAUS,       ! SHAPS(NNODE,NGAUS)
     . IFOR25=IFOR24+MNODE*MGAUS,             ! GPCOS(NDIME,NGAUS)
     . LFOR  =(IFOR25+MDIME*MGAUS)*MDYNA )
C
      PARAMETER(
c    . MGASS=max(    3,MGAUS) )               ! SG
     . MGASS=   (    3+MGAUS) )               ! PC & linux
C
      PARAMETER(
     . ILOS1=1,                               ! ALOAD(NEVAB)
     . ILOS2=ILOS1+MEVAB,                     ! DERIV(NDIME,NNODE)
     . ILOS3=ILOS2+MDIME*MNODE,               ! ELEDG(NDIME,NNODE)
     . ILOS4=ILOS3+MDIME*MNODE,               ! GPCOD(NDIME)
     . ILOS5=ILOS4+MDIME,                     ! GVECT(NDIME)
     . ILOS6=ILOS5+MDIME,                     ! NOPRS(NNODE)
     . ILOS7=ILOS6+MNODE,                     ! POSGP(NDIME,NGAUS)
     . ILOS8=ILOS7+MDIME*MGASS,               ! PRESS(NNODE,NDIME)
     . ILOS9=ILOS8+MNODE*MDIME,               ! PVECT(NDIME)
     . ILOS10=ILOS9+MDIME,                    ! SHAPE(NNODE)
     . ILOS11=ILOS10+MNODE,                   ! VAREA(NDIME)
     . ILOS12=ILOS11+MDIME,                   ! WEIGP(NGAUS)
     . ILOS13=ILOS12+MGASS,                   ! XJACM(NDIME*NDIME)
     . ILOS14=ILOS13+MDIME*MDIME,             ! CARTD(NDIME,NNODE)
     . ILOS15=ILOS14+MDIME*MDIME,             ! PRESC(NDOFC)
     . ILOS16=ILOS15+MDOFC )                  ! ELCO1(NDIME,NNODE)
      PARAMETER(
     . ILOS17=ILOS16+MDIME*MNODE,             ! CARTD(NDIME,NNODE,NGAUS)
     . ILOS18=ILOS17+MDIME*MNODE*MGASS,       ! DVOLU(NGAUS)
     . ILOS19=ILOS18+MGASS,                   ! GPCOD(NDIME,NGAUS)
     . ILOS20=ILOS19+MDIME*MGASS,             ! SHAPE(NNODE,NGAUS)
     . ILOS21=ILOS20+MNODE*MGASS,             ! BSBAR(NSTR1,NEVAB,NGAUS)
     . ILOS22=ILOS21+MSTR1*MEVAB*MGASS,       ! DERIV(NDIME,NNODE,NGAUS)
     . ILOS23=ILOS22+MDIME*MNODE*MGASS,       ! POSGP(NDIME,NGAUS)
     . ILOS24=ILOS23+MDIME*MGASS,             ! WEIGP(NGAUS)
     . ILOS25=ILOS24+MGASS,                   ! XJACM(NDIME,NDIME)
     . ILOS26=ILOS25+MDIME*MDIME,             ! RMAT2(NSTR1,NSTR1)
     . ILOS27=ILOS26+MSTR1*MSTR1,             ! CMEAN(3,NNODE)
     . ILOS28=ILOS27+3*MNODE,                 ! SHAPR(NNODE,NGAUS)
     . ILOS29=ILOS28+MNODE*MGASS,             ! WSTIR(NNODE,NNODE)
     . ILOS30=ILOS29+MNODE*MNODE,             ! VNORL(NDIME,NNODG)
     . ILOS31=ILOS30+MDIME*MNODG,             ! EMASS(NNODE,NNODE)
     . ILOS32=ILOS31+MNODE*MNODE )            ! ELDI1(NDOFN,NNODE)
      PARAMETER(
     . ILOS33=ILOS32+MDOFC*MNODE,             ! CARTS(NDIME,NNODE,NGAUS)
     . ILOS34=ILOS33+MDIME*MNODE*MGASS,       ! SHAPS(NNODE,NGAUS)
     . ILOS35=ILOS34+MNODE*MGASS,             ! GPCOS(NDIME,NGAUS)
     . ILOS36=ILOS35+MDIME*MGASS,             ! ISKPO(NSKEW)
     . LLOS  =ILOS36+(MSKEW*MCHAL+4)/8*MSKE1 )
C
      PARAMETER(
     . IGSM1 =1,                              ! SMSTS(NSTR1,NPOIN)
     . IGSM2 =IGSM1+MSTR1*MPOIN*MSMOM,        ! SMSTN(NSTR1,NPOIN)
     . IGSM3 =IGSM2+MSTR1*MPOIN*MSMOM,        ! ACCPN(NPOIN)
     . IGSM4 =IGSM3+MPOIN*MSMOM,              ! STREB(NSTR1,NNODE)
     . IGSM5 =IGSM4+MSTR1*MNODE*MSMOM,        ! STREA(NNODE,NSTR1)
     . IGSM6 =IGSM5+MNODE*MSTR1*MSMOM,        ! SFISB(NSTR1,NNODE)
     . IGSM7 =IGSM6+MSTR1*MNODE*MSMOM,        ! SFISA(NNODE,NSTR1)
     . IGSM8 =IGSM7+MNODE*MSTR1*MSMOM,        ! SMSTP(NNUIN,NPOIN)
     . IGSM9 =IGSM8+MNUIN*MPOIN*MSMOM,        ! SFIPA(NNUIN,NNODE)
     . IGSM10=IGSM9+MNUIN*MNODE*MSMOM,        ! SFIPB(NNODE,NNUIN)
     . IGSM11=IGSM10+MNODE*MNUIN*MSMOM,       ! CARTD(NDIME,NNODE,NGAUS)
     . IGSM12=IGSM11+MDIME*MNODE*MGAUS,       ! DVOLU(NGAUS)
     . IGSM13=IGSM12+MGAUS,                   ! GPCOD(NDIME,NGAUS)
     . IGSM14=IGSM13+MDIME*MGAUS,             ! SHAPE(NNODE,NGAUS)
     . IGSM15=IGSM14+MNODE*MGAUS,             ! BSBAR(NSTR1,NEVAB,NGAUS)
     . IGSM16=IGSM15+MSTR1*MEVAB*MGAUS )      ! DERIV(NDIME,NNODE,NGAUS)
      PARAMETER(
     . IGSM17=IGSM16+MDIME*MNODE*MGAUS,       ! POSGP(NDIME,NGAUS)
     . IGSM18=IGSM17+MDIME*MGAUS,             ! WEIGP(NGAUS)
     . IGSM19=IGSM18+MGAUS,                   ! XJACM(NDIME,NDIME)
     . IGSM20=IGSM19+MDIME*MDIME,             ! RMAT2(NSTR1,NSTR1)
     . IGSM21=IGSM20+MSTR1*MSTR1,             ! CMEAN(3,NNODE)
     . IGSM22=IGSM21+3*MNODE,                 ! SHAPR(NNODE,NGAUS)
     . IGSM23=IGSM22+MNODE*MGAUS,             ! WSTIR(NNODE,NNODE)
     . IGSM24=IGSM23+MNODE*MNODE,             ! VNORL(NDIME,NNODG)
     . IGSM25=IGSM24+MDIME*MNODG,             ! ELCO1(NDIME,NNODE)
     . IGSM26=IGSM25+MDIME*MNODE,             ! EMASS(NNODE,NNODE)
     . IGSM27=IGSM26+MNODE*MNODE,             ! ELDI1(NDOFN,NNODE)
     . IGSM28=IGSM27+MDOFC*MNODE,             ! CARTS(NDIME,NNODE,NGAUS)
     . IGSM29=IGSM28+MDIME*MNODE*MGAUS,       ! SHAPS(NNODE,NGAUS)
     . IGSM30=IGSM29+MNODE*MGAUS,             ! GPCOS(NDIME,NGAUS)
     . IGSM31=IGSM30+MNODE*MGAUS,             ! XJACI(NDIME,NDIME)
     . LGSM  =IGSM31+MDIME*MDIME )
C
      PARAMETER(
     . ISET1=1,                               ! DERIV(NDIME,NNODE,NGAUS)
     . ISET2=ISET1+MDIME*MNODE*MGAUS,         ! POSGP(NDIME,NGAUS)
     . ISET3=ISET2+MDIME*MGAUS,               ! WEIGP(NGAUS)
     . ISET4=ISET3+MGAUS,                     ! XJACM(NDIME,NDIME)
     . ISET5=ISET4+MDIME*MDIME,               ! RMAT2(NSTR1,NSTR1)
     . ISET6=ISET5+MSTR1*MSTR1,               ! CMEAN(3,NNODE)
     . ISET7=ISET6+3*MNODE,                   ! SHAPR(NNODE,NGAUS) B-bar
     . ISET8=ISET7+MNODE*MGAUS,               ! WSTIR(NNODE,NNODE) B-bar
     . ISET9=ISET8+MNODE*MNODE,               ! VNORL(NDIME,NNODG)
     . ISET10=ISET9+MDIME*MNODG,              ! ELCO1(NDIME,NNODE)
     . ISET11=ISET10+MDIME*MNODE,             ! CARTD(NDIME,NNODE,NGAUS)
     . ISET12=ISET11+MDIME*MNODE*MGAUS,       ! DVOLU(NGAUS)
     . ISET13=ISET12+MGAUS,                   ! GPCOD(NDIME,NGAUS)
     . ISET14=ISET13+MDIME*MGAUS,             ! SHAPE(NNODE,NGAUS)
     . ISET15=ISET14+MNODE*MGAUS,             ! BSBAR(NSTR1,NEVAB,NGAUS)
     . ISET16=ISET15+MSTR1*MEVAB*MGAUS )      ! EMASS(NNODE,NNODE)
      PARAMETER(
     . ISET17=ISET16+MNODE*MNODE,             ! ELDI1(NDOFN,NNODE)
     . ISET18=ISET17+MDOFC*MNODE,             ! CARTS(NDIME,NNODE,NGAUS)
     . ISET19=ISET18+MDIME*MNODE*MGAUS,       ! SHAPS(NNODE,NGAUS)
     . ISET20=ISET19+MNODE*MGAUS,             ! GPCOS(NDIME,NGAUS)
     . LSET  =ISET20+MDIME*MGAUS )
C
      PARAMETER(
     . ISTI1=1,                               ! BMATX(NSTR1,NEVAB)
     . ISTI2=ISTI1+MSTR1*MEVAB,               ! DMATX(NSTR1,NSTR1)
     . ISTI3=ISTI2+MSTR1*MSTR1,               ! SIGMA(NSTR1)
     . ISTI4=ISTI3+MSTR1,                     ! XJACM(NDIME,NDIME)
     . ISTI5=ISTI4+MDIME*MDIME,               ! TENOD(NNODE)
     . ISTI6=ISTI5+MNODE,                     ! DTENO(NNODE)
     . ISTI7=ISTI6+MNODE,                     ! DISPL(NEVAB)
     . ISTI8=ISTI7+MEVAB,                     ! VNORL(NDIME,NNODG)
     . ISTI9=ISTI8+MDIME*MNODG,               ! ELCO1(NDIME,NNODE)
     . ISTI10=ISTI9 +MDIME*MNODE,             ! CARTD(NDIME,NNODE,NGAUS)
     . ISTI11=ISTI10+MDIME*MNODE*MGAUS,       ! DVOLU(NGAUS)
     . ISTI12=ISTI11+MGAUS,                   ! GPCOD(NDIME,NGAUS)
     . ISTI13=ISTI12+MDIME*MGAUS,             ! SHAPE(NNODE,NGAUS)
     . ISTI14=ISTI13+MNODE*MGAUS,             ! BSBAR(NSTR1,NEVAB,NGAUS)
     . ISTI15=ISTI14+MSTR1*MEVAB*MGAUS,       ! DERIV(NDIME,NNODE,NGAUS)
     . ISTI16=ISTI15+MDIME*MNODE*MGAUS )      ! POSGP(NDIME,NGAUS)
      PARAMETER(
     . ISTI17=ISTI16+MDIME*MGAUS,             ! WEIGP(NGAUS)
     . ISTI18=ISTI17+MGAUS,                   ! XJACM(NDIME,NDIME)
     . ISTI19=ISTI18+MDIME*MDIME,             ! RMAT2(NSTR1,NSTR1)
     . ISTI20=ISTI19+MSTR1*MSTR1,             ! CMEAN(3,NNODE)
     . ISTI21=ISTI20+3*MNODE,                 ! SHAPR(NNODE,NGAUS)
     . ISTI22=ISTI21+MNODE*MGAUS,             ! WSTIR(NNODE,NNODE)
     . ISTI23=ISTI22+MNODE*MNODE,             ! VNORL(NDIME,NNODG)
     . ISTI24=ISTI23+MDIME*MNODG,             ! EMASS(NNODE,NNODE)
     . ISTI25=ISTI24+MNODE*MNODE,             ! ESTII(NKOVA)
     . ISTI26=ISTI25+MKOVA,                   ! WSTII(NKOVA)
     . ISTI27=ISTI26+MKOVA*MDYNA,             ! ELDI1(NDOFN,NNODE)
     . ISTI28=ISTI27+MDOFC*MNODE,             ! CSTI1(NEVAB,NEVAB)
     . ISTI29=ISTI28+MEVAB*MEVAB*MMEMO7M,     ! CARTS(NDIME,NNODE,NGAUS)
     . ISTI30=ISTI29+MDIME*MNODE*MGAUS,       ! SHAPS(NNODE,NGAUS)
     . ISTI31=ISTI30+MNODE*MGAUS,             ! GPCOS(NDIME,NGAUS)
     . ISTI32=ISTI31+MDIME*MGAUS )            ! AUXMA(NEVAB,NEVAB)
      PARAMETER(
     . ISTI33=ISTI32+MEVAB*MEVAB*MSKE1,       ! TRAMA(NEVAB,NEVAB)
     . LSTI  =ISTI33+MEVAB*MEVAB*MSKE1 )
C
      PARAMETER(
     . ILOSY1=1,                              ! AUXMA(NEVAB,NEVAB)
     . ILOSY2=ILOSY1+MEVAB*MEVAB,             ! TRAMA(NEVAB,NEVAB)
     . LSKE  =(ILOSY2+MEVAB*MEVAB)*MSKE1 )
C
C**** ADDSOL
C
C**** 1) SKYLINE SOLVER
C
      PARAMETER(
     . ISOL1S=1,                              ! GSTDI(NEQNS)
     . ISOL2S=ISOL1S+MEQNS,                   ! GSTLO(NLAST*IUNSY)
     . ISOL3S=ISOL2S+MLAST*MUNS2,             ! GSTUP(NLAST)
     . ISOL4S=ISOL3S+MLAST,                   ! CSTIF(NEVAC,NEVAC)
     . ISOL5S=ISOL4S+MEVAC*MEVAC,             ! ELOAD(NEQNS)
     . ISOL6S=ISOL5S+MEQNS,                   ! LNUEQ(NTOTV)
     . ISOL7S=ISOL6S+(MTOTV*MCHAL+4)/8,       ! LPONT(NTOTV)
     . ISOL8S=ISOL7S+(MTOTV*MCHAL+4)/8,       ! DISIM(NEVAC)
     . ISOL9S=ISOL8S+MEVAC,                   ! FOREL(NEVAC)
     . ISOL10S=ISOL9S+MEVAC,                  ! ALOAD(NEQNS*IFITE)
     . ISOL11S=ISOL10S+MEQNS*MFITE,           ! DELTA(NEQNS*IFITE)
     . ISOL12S=ISOL11S+MEQNS*MFITE,           ! DISIC(NEQNS*IFITE)
     . ISOL13S=ISOL12S+MEQNS*MFITE,           ! LOCAL(NEVAC*IFITE)
     . LSOLS  =ISOL13S+(MEVAC*MFITE*MCHAL+4)/8 )
C
      PARAMETER(
     . LSOL1=MSOL1*LSOLS )
C
C**** 2) FRONTAL SOLVER
C
      PARAMETER(
     . ISOL1F=1,                              ! EQRHS(NBUFA)
     . ISOL2F=ISOL1F+MBUFA,                   ! EQUAT(NFRON,NBUFA)
     . ISOL3F=ISOL2F+MFRON*MBUFA,             ! GLOAD(NFRON)
     . ISOL4F=ISOL3F+MFRON,                   ! GSTIF(NSTIF)
     . ISOL5F=ISOL4F+MSTIF,                   ! LOCEL(NEVAC)
     . ISOL6F=ISOL5F+(MEVAC*MCHAL+4)/8,       ! NACVA(NFRON)
     . ISOL7F=ISOL6F+(MFRON*MCHAL+4)/8,       ! NAMEV(NBUFA)
     . ISOL8F=ISOL7F+(MBUFA*MCHAL+4)/8,       ! NDEST(NEVAC)
     . ISOL9F=ISOL8F+(MEVAC*MCHAL+4)/8,       ! NPIVO(NBUFA)
     . ISOL10F=ISOL9F+(MBUFA*MCHAL+4)/8,      ! VECRV(NFRON)
     . ISOL11F=ISOL10F+MFRON,                 ! CSTIF(NEVAC,NEVAC)
     . ISOL12F=ISOL11F+MEVAC*MEVAC,           ! EQCOL(IUNSY*NFRON*NBUFA)
     . LSOLF  =ISOL12F+MUNS2*MFRON*MBUFA )
C
      PARAMETER(
     . LSOL2=MSOL2*LSOLF )
C
C**** 3) PCG SOLVER
C
      PARAMETER(
     . ISOL1P=1,                              ! GSTDI(NPOIN*NSIZE)
     . ISOL2P=ISOL1P+MPOIN*MSIZE,             ! CSTIF(NEVAC,NEVAC)
     . ISOL3P=ISOL2P+MEVAC*MEVAC,             ! DISIM(NEVAC)
     . ISOL4P=ISOL3P+MEVAC,                   ! FOREL(NEVAC)
     . ISOL5P=ISOL4P+MEVAC,                   ! ALOAD(NTOTV)
     . ISOL6P=ISOL5P+MTOTV,                   ! DELTA(NTOTV)
     . LSOLP =ISOL6P+MTOTV )
C
      PARAMETER(
     . LSOL3=MSOL3*LSOLP )
C
      PARAMETER(
c    . LSOL=max(LSOL1,LSOL2,LSOL3) )          ! SG
     . LSOL=(LSOL1+LSOL2+LSOL3) )             ! PC & linux
C
C**** SCRATCH VECTOR (WORK1)
C
      PARAMETER(
c    . MWORK1=max(LIND,LREN,LQUA,LFIX,LFOT,LFOR,   ! SG
c    .            LLOS,LGSM,LSET,LSTI,LSKE,LSOL) )
     . MWORK1=(LIND+LREN+LQUA+LFIX+LFOT+LFOR+      ! PC & linux
     .         LLOS+LGSM+LSET+LSTI+LSKE+LSOL) )
C
