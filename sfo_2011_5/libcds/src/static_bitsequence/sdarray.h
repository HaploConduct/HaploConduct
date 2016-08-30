
#ifndef SDARRAY_H
#define SDARRAY_H

#include <basics.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>

typedef struct {
  int n,m;
  int size;
  uchar *buf;
  uint *lp;
  uint *sl;
  ushort *ss;
  uint ss_len, sl_len;
  uint *p;
} selectd2;

typedef struct {
  int n,m,d;
  int size;
  uchar *hi;
  uint *low;
  selectd2 *sd0,*sd1;
	uint hi_len, low_len;
	//uint lasti, lasts;
} selects3;

int selects3_construct(selects3 *select, int n, uint *buf);
int selects3_select(selects3 *select, int i);
int selects3_rank(selects3 *select, int i);
int selects3_selectnext(selects3 *select, int i);

void make___selecttbl(void);
int __setbit(uint *B, int i,int x);
int selectd2_save(selectd2 * s, FILE * fp);
int selects3_save(selects3 * s, FILE * fp);

int selectd2_load(selectd2 * s, FILE * fp);
int selects3_load(selects3 * s, FILE * fp);

void selectd2_free(selectd2 * s);
void selects3_free(selects3 * s);


#endif

