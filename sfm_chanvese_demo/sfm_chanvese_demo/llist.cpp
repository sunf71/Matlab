#include "llist.h"

PT *pt_create(long x, long y, long z, long idx){
  PT* newpt = (PT*)mxMalloc(sizeof(PT));
  if(newpt == NULL) return NULL;

  newpt->x = x;
  newpt->y = y;
  newpt->z = z;
  newpt->idx = idx;
  newpt->prev = NULL;
  newpt->next = NULL;
  return newpt;
}

LL *ll_create(){
  LL *newll = (LL*)mxMalloc(sizeof(LL));
  if(newll == NULL) return NULL;
  newll->head = NULL;
  newll->curr = NULL;
  newll->length = 0;
  return newll;
}

void ll_push(LL *list, PT *add){
  if(add == NULL || list == NULL) return;
  add->next = list->head;
  add->prev = NULL;
  if(add->next != NULL){
    add->next->prev = add;
  }
  list->head = add;
  list->length++;
}

void ll_pushnew(LL *list, long x, long y, long z, long idx){
  if(list == NULL) return;
  PT* add = pt_create(x,y,z,idx);
  if(add == NULL) return;
  ll_push(list,add);
}

void ll_destroy(LL *list){
  if(list==NULL) return;
  while(list->head != NULL){
    ll_pop_free(list);
  }
  mxFree(list);
}

void ll_remcurr_free(LL *list){
  PT* p = ll_remcurr(list);
  if(p != NULL) mxFree(p);
}

PT *ll_remcurr(LL *list){
  if(list == NULL) return NULL;
  PT* out = list->curr;
  if(out == list->head){
    return ll_pop(list);
  }
  else
  {
    if(out != NULL)
    {
      if(out->next != NULL){
        out->next->prev = out->prev;
      }
      if(out->prev != NULL){
        out->prev->next = out->next;
      }
      list->curr = out->next;
      list->length--;
    }
    return out;
  }
}

void ll_pop_free(LL *list){
  PT* p = ll_pop(list);
  if(p != NULL) mxFree(p);
}

PT *ll_pop(LL *list){
  if(list == NULL) return NULL;
  PT *out = list->head;
  if(out != NULL){
    list->head = out->next;
    if(list->curr == out) list->curr = list->head;
    if(list->head != NULL) list->head->prev = NULL;
    list->length--;
  }
  return out;
}

void ll_init(LL *list){
  if(list == NULL) return;
  list->curr = list->head;
}

void ll_step(LL *list){
  if(list == NULL) return;
  if(list->curr != NULL){
    list->curr = list->curr->next;
  }
}











