#define _GNU_SOURCE
#include <sched.h>

int findmycpu_()
{
  int cpu;
  cpu = sched_getcpu();
  return cpu;
}
