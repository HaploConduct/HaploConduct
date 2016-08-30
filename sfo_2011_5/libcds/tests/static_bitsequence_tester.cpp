/* static_bitsequence_tester.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_bitsequence_tester
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <iostream>
#include <algorithm>
#include <sstream>
#include <basics.h>
#include <static_bitsequence.h>

using namespace std;

/* Time meassuring */
double ticks= (double)sysconf(_SC_CLK_TCK);
struct tms t1,t2;

void start_clock() {
  times (&t1);
}


double stop_clock() {
  times (&t2);
  return (t2.tms_utime-t1.tms_utime)/ticks;
}
/* end Time meassuring */

uint NQUERIES=10000000;
uint SEED=47;

void load(char *fname, uint ** text, uint * n) {
  FILE * fp = fopen(fname,"r");
  if(fp==NULL) {
    cout << "could not open " << fname << endl;
    return;
  }
  if(fread(n,sizeof(uint),1,fp)!=1) {
    cout << "Error reading file " << fname << endl;
    return;
  }
  *text = new uint[uint_len(*n,1)];

  if(fread(*text,sizeof(uint),uint_len(*n,1),fp)!=uint_len(*n,1)) {
    cout << "Error reading file " << fname << endl;
    return;
  }
}

void test_bitsequence(uint * bitseq, uint len, static_bitsequence * bs) {
  uint ones = 0;
	uint last_one = 0;
  bool error = false;
  for(uint i=0;i<len && !error;i++) {
    //cout << "i="<<i<< endl;
    if(i>0) {
      if(i%max((uint)1,(bs->length()/100))==0) { cout << "."; cout.flush(); }
      if(i%max((uint)1,(bs->length()/10))==0) { cout << endl; cout.flush(); }
    }
    if(bitget(bitseq,i)) {
			for(uint k=last_one; !error && k<i;k++) {
				if(bs->select_next1(k)!=i) {
					uint ans= bs->select_next1(k);
					cout << "Error select_next1" << endl;
					cout << " got: (k=" << k << ") " << ans << " expected: " << i << endl;
					cout << " rank(" << k << ")=" << bs->rank1(k) << " access(" << k << ")=" << bs->access(k) << endl;
					cout << " rank(" << ans << ")=" << bs->rank1(ans) << " access(" << ans << ")=" << bs->access(ans) << endl;
					error = true;
				}
			}
			last_one = i;
			ones++;
		}
    if(bs->access(i) != (bitget(bitseq,i)!=0)) {
      cout << "Access error for position " << i << endl;
      cout << " got: " << bs->access(i) << " expected: " << (bitget(bitseq,i)!=0) << endl;
      error = true;
    }
    if(bs->rank1(i) != ones) {
      cout << "Rank1 error for position " << i << endl;
      cout << " got: " << bs->rank1(i) << " expected: " << ones << endl;
      error = true;
    } 
    if(bitget(bitseq,i) && bs->select1(ones) != i) {
      cout << "Select1 error for position " << i << " ones:" << ones << endl;
      cout << " got: " << bs->select1(ones) << " expected: " << i << endl;
      error = true;
    }
    if(bs->rank0(i) != i+1-ones) {
      cout << "Rank0 error for position " << i << endl;
      cout << " got: " << bs->rank0(i) << " expected: " << ones << endl;
      error = true;
    }
    if(!bitget(bitseq,i) && bs->select0(i+1-ones) != i) {
      cout << "Select0 error for position " << i << endl;
      cout << " got: " << bs->select0(i+1-ones) << " expected: " << i << endl;
      error = true;
    }
  }
	cout << "." << endl;
}

void speed_access(static_bitsequence * ss, uint * bitseq, uint n) {
  uint acc=0;
  srand(SEED);

  start_clock();
  for(uint i=0;i<NQUERIES;i++) {
    uint pos = rand()%n;
    acc += ss->access(pos);
  }
  double t = stop_clock();
  cout << " * Time for " << NQUERIES << " accesses: " << t << " secs" << endl;
  cout << " * Time per access: " << 1000*t/NQUERIES << " msecs" << endl;
  cout << " - Check sum: " << acc << endl;
}

void speed_rank0(static_bitsequence * ss, uint * bitseq, uint n) {
  uint acc=0;
  srand(SEED);

  start_clock();
  for(uint i=0;i<NQUERIES;i++) {
    uint pos = rand()%n;
    acc += ss->rank0(pos);
  }
  double t = stop_clock();
  cout << " * Time for " << NQUERIES << " rank0s: " << t << " secs" << endl;
  cout << " * Time per rank0: " << 1000*t/NQUERIES << " msecs" << endl;
  cout << " - Check sum: " << acc << endl;
}

void speed_rank1(static_bitsequence * ss, uint * bitseq, uint n) {
  uint acc=0;
  srand(SEED);

  start_clock();
  for(uint i=0;i<NQUERIES;i++) {
    uint pos = rand()%n;
    acc += ss->rank1(pos);
  }
  double t = stop_clock();
  cout << " * Time for " << NQUERIES << " rank1s: " << t << " secs" << endl;
  cout << " * Time per rank1: " << 1000*t/NQUERIES << " msecs" << endl;
  cout << " - Check sum: " << acc << endl;
}

void speed_select0(static_bitsequence * ss, uint * bitseq, uint n) {
  uint acc=0;
  uint ones=ss->rank0(n-1);
  srand(SEED);

  start_clock();
  for(uint i=0;i<NQUERIES;i++) {
    uint pos = rand()%ones+1;
    acc += ss->select0(pos);
  }
  double t = stop_clock();
  cout << " * Time for " << NQUERIES << " select0s: " << t << " secs" << endl;
  cout << " * Time per select0: " << 1000*t/NQUERIES << " msecs" << endl;
  cout << " - Check sum: " << acc << endl;
}

void speed_select1(static_bitsequence * ss, uint * bitseq, uint n) {
  uint acc=0;
  uint ones=ss->rank1(n-1);
  srand(SEED);

  start_clock();
  for(uint i=0;i<NQUERIES;i++) {
    uint pos = rand()%ones+1;
    acc += ss->select1(pos);
  }
  double t = stop_clock();
  cout << " * Time for " << NQUERIES << " select1s: " << t << " secs" << endl;
  cout << " * Time per select1: " << 1000*t/NQUERIES << " msecs" << endl;
  cout << " - Check sum: " << acc << endl;
}

void speed_selectnext1(static_bitsequence * ss, uint * bitseq, uint n) {
  uint acc=0;
  srand(SEED);

  start_clock();
  for(uint i=0;i<NQUERIES;i++) {
    uint pos = rand()%n;
    acc += ss->select_next1(pos);
  }
  double t = stop_clock();
  cout << " * Time for " << NQUERIES << " select_next1s: " << t << " secs" << endl;
  cout << " * Time per select_next1: " << 1000*t/NQUERIES << " msecs" << endl;
  cout << " - Check sum: " << acc << endl;
}
