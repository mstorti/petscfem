//__INSERT_LICENSE__
// $Id: mpelog.cpp,v 1.3.2.1 2004/09/25 23:20:04 mstorti Exp $

#include <mpi.h>
#include <mpe.h>
#include <src/mpelog.h>

#define MPE_LOG_INIT3(event) int start_##event,end_##event

MPE_LOG_INIT3(comp);
MPE_LOG_INIT3(assmbly);
MPE_LOG_INIT3(upl);
MPE_LOG_INIT3(passm);
MPE_LOG_INIT3(passmb);
MPE_LOG_INIT3(aux);
MPE_LOG_INIT3(aux1);
MPE_LOG_INIT3(aux2);
MPE_LOG_INIT3(aux3);

extern int MY_RANK;
int mpe_initialized = 0;

#define MPE_LOG_INIT(event,color)				\
start_##event = MPE_Log_get_event_number();			\
end_##event = MPE_Log_get_event_number();			\
if (!MY_RANK) {							\
   MPE_Describe_state(start_##event,end_##event,#event,#color);	\
}

void mpe_initialize() {
  if (!mpe_initialized) {
    mpe_initialized = 1;
    MPE_LOG_INIT(comp,green);
    MPE_LOG_INIT(assmbly,red);
    MPE_LOG_INIT(upl,blue);
    MPE_LOG_INIT(passm,yellow);
    MPE_LOG_INIT(passmb,brown);
    MPE_LOG_INIT(aux,magenta);
    MPE_LOG_INIT(aux1,cyan);
    MPE_LOG_INIT(aux2,turquoise);
    MPE_LOG_INIT(aux3,aquamarine);
  }
}

