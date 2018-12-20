#include <stdio.h>
#include <malloc.h>

/* #include <invent.h> */

#include <sys/resource.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>

/* #include <sys/swap.h> */

#include "gen.h"

double DynamicMemoryUsage(double *used, double *freed, int kilo)
{
  struct mallinfo mi;
  double alloc;
  
  /*
  mi = mallinfo();
  alloc = (double)mi.arena*GEN_MB_PER_B;
  if (kilo) {
  */
                                 /* report in kilobytes */
  /*
    *used  = (double)(mi.uordblks+mi.usmblks)*GEN_KB_PER_B;
    *freed = (double)(mi.fordblks+mi.fsmblks)*GEN_KB_PER_B;
  }
  else {
  */
                                 /* report in megabytes */
  /*
    *used  = (double)(mi.uordblks+mi.usmblks)*GEN_MB_PER_B;
    *freed = (double)(mi.fordblks+mi.fsmblks)*GEN_MB_PER_B;
  }
  */
  
  return alloc;
}

double LogicalSwapSpace(double *rss, double *swap_p, double *swap_v)
{
  
  /*
  struct rlimit rlp;
  inventory_t *inv;
  off_t arg;
  double core, heap, logical, swap, swap_l;

  while (inv = getinvent())
    if (inv->inv_class == INV_MEMORY)
      if (inv->inv_type == INV_MAIN_MB) {
        core = (double)(inv->inv_state);
	break;
      }
  endinvent();

  getrlimit(RLIMIT_DATA, &rlp);
  heap = (double)(rlp.rlim_cur)*GEN_MB_PER_B;
  getrlimit(RLIMIT_RSS,  &rlp);
  *rss  = (double)(rlp.rlim_cur)*GEN_MB_PER_B;
  swapctl(SC_GETSWAPTOT,  &arg);
  *swap_p = (double)arg*GEN_B_PER_BLOCK*GEN_MB_PER_B;
  swapctl(SC_GETSWAPVIRT, &arg);
  *swap_v = (double)arg*GEN_B_PER_BLOCK*GEN_MB_PER_B;
  swapctl(SC_GETLSWAPTOT, &arg);
  swap_l = (double)arg*GEN_B_PER_BLOCK*GEN_MB_PER_B;
  swap = *swap_p + *swap_v;

  logical = *rss + swap;
  */

  double logical;
  

  return logical;
}

char *MemoryStatus(void)
{
  double alloc, freed, logical, rss, swap_p, swap_v, used;
  static char status[256];

  alloc = DynamicMemoryUsage(&used, &freed, 0);
  logical = LogicalSwapSpace(&rss, &swap_p, &swap_v);

  sprintf(status, " arena: %5.1fMB = %5.1fMB used + %5.1fMB free\n limit: %5.1fMB = %5.1fMB core + %5.1fMB pswap + %5.1fMB vswap",
	 alloc, used, freed, logical, rss, swap_p, swap_v);

  /* beg */
  printf("Ignore the returned value from MemoryStatus()\n");
  /* end */

  return &status[0];
}
