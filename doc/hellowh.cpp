#include <src/hook.h>
#include <src/dlhook.h>

class hello_world_hook {
public:
  void init(Mesh &mesh_a,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name) {
    printf("Hello world! I'm in the \"init\" hook\n");
  }
  void time_step_pre(double time,int step) {
    printf("Hello world! I'm in the \"time_step_pre\" hook\n");
  }
  void time_step_post(double time,int step,
		      const vector<double> &gather_values) {
    printf("Hello world! I'm in the \"time_step_post\" hook\n");
  }
  void close() {
    printf("Hello world! I'm in the \"close\" hook\n");
  }
};

DL_GENERIC_HOOK(hello_world_hook);
