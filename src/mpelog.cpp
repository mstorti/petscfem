//__INSERT_LICENSE__
// $Id: mpelog.cpp,v 1.1 2004/07/29 22:37:49 mstorti Exp $

#include <src/mpelog.h>

int mpe_initialize() {
  if (!mpe_initialized) {
    mpe_initialized = 1;
    start_comp = MPE_Log_get_event_number();
    end_comp = MPE_Log_get_event_number();
    start_assmbly = MPE_Log_get_event_number();
    end_assmbly = MPE_Log_get_event_number();
    start_assm = MPE_Log_get_event_number();
    end_assm = MPE_Log_get_event_number();
    if (!myrank) {
      MPE_Describe_state(start_comp,end_comp,"comp","green:gray");
      MPE_Describe_state(start_assmbly,end_assmbly,"assmbly","red:white");
      MPE_Describe_state(start_assm,end_assm,"assm","red:white");
    }
  }
}
