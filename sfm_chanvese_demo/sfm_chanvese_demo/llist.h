#ifndef __LLIST_H
#define __LLIST_H

#include <stdlib.h>
#include <mex.h>

/* pt = point */
/* ll = linked list */

struct pt{
    long x;
    long y;
    long z;
    long idx;
    struct pt *prev;
    struct pt *next;
};

typedef struct pt PT;

struct ll{
  PT *head;
  PT *curr;
  long length;
};

typedef struct ll LL;

void ll_push(LL *list, PT *add);
void ll_pushnew(LL *list, long x, long y, long z, long idx);
void ll_remcurr_free(LL *list);
PT *ll_remcurr(LL* list);
void ll_destroy(LL *list);
PT *pt_create(long x, long y, long z, long idx);
LL *ll_create();
void ll_pop_free(LL *list);
PT *ll_pop(LL *list);
void ll_init(LL *list);
void ll_step(LL *list);



#endif
