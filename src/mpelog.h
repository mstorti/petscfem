//__INSERT_LICENSE__
// $Id: mpelog.h,v 1.3.2.1 2004/09/25 23:20:04 mstorti Exp $

extern int mpe_initialized;

#define MPE_LOG_INIT2(event) extern int start_##event,end_##event

#define MPE_START(event) MPE_Log_event(start_##event,0,"start-" #event)
#define MPE_END(event) MPE_Log_event(end_##event,0,"end-" #event)

MPE_LOG_INIT2(comp);
MPE_LOG_INIT2(assmbly);
MPE_LOG_INIT2(upl);
MPE_LOG_INIT2(passm);
MPE_LOG_INIT2(passmb);
MPE_LOG_INIT2(aux);
MPE_LOG_INIT2(aux1);
MPE_LOG_INIT2(aux2);
MPE_LOG_INIT2(aux3);

void mpe_initialize();

