/* static_sequence_tester.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_sequence_tester
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

#include "static_sequence_tester.h"

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

uint NQUERIES=100000;
uint SEED=47;

void test_static_sequence(uint * symbols, uint n, static_sequence * ss) {
  cout << "Size: " << ss->size() << endl;
  uint max_v=0;
  for(uint i=0;i<n;i++) 
    max_v = max(max_v,symbols[i]);
  uint * occ = new uint[max_v+1];
  for(uint i=0;i<=max_v;i++)
    occ[i] = 0;
  bool error = false;
  for(uint i=0;i<n && !error;i++) {
    if(i!=0 && i%max((uint)1,(n-1)/100)==0) { cout << "."; cout.flush(); }
    if(i!=0 && i%max((uint)1,(n-1)/10)==0) cout << endl;
    occ[symbols[i]]++;
    uint a = /*symbols[i]; /*/ss->access(i);
    uint r = /*occ[symbols[i]];/*/ss->rank(symbols[i],i);
    uint s = /*i; /*/ss->select(symbols[i],occ[symbols[i]]);
    uint rM1 = (i==0)?0:ss->rank(symbols[i],i-1);
    if(r!=occ[symbols[i]]) {
      cout << "Error in rank for symbol " << symbols[i] << " at position " << i << endl;
      cout << "value: " << r << endl;
      cout << "Expected: " << occ[symbols[i]] << endl;
      error = true;
    }
    if(s!=i) {
      cout << "Error in select for symbol " << symbols[i] << " at position " << occ[symbols[i]] << endl;
      cout << "value: " << s << endl;
      cout << "Expected: " << i << endl;
      error = true;
    }
    if(a!=symbols[i]) {
      cout << "Error in access at position " << i << endl;
      cout << "value: " << a << endl;
      cout << "Expected: " << symbols[i] << endl;
      error = true;
    }
    if(rM1!=occ[symbols[i]]-1) {
      cout << "Error in rankM1 for symbol " << symbols[i] << " at position " << i-1 << endl;
      cout << "value: " << rM1 << endl;
      cout << "Expected: " << occ[symbols[i]]-1 << endl;
      error = true;
    }
  }
  if(!error)
    cout << "Test OK! It works :)" << endl;
  delete [] occ;  
}

void load(char *fname, uint ** text, uint * n) {
  struct stat text_info;
  if(stat(fname,&text_info)<0) {
    cout << "could not stat: " << fname << endl;
    return;
  }
  
  *n= (uint)text_info.st_size/4;
  *text = new uint[*n+1];
  FILE * fp = fopen(fname,"r");
  if(fp==NULL) {
    cout << "could not open " << fname << endl;
    return;
  }

  cout << "File: " << fname << endl;
  cout << "Length: " << *n << endl;

  uint max_symbol = 0;
  for(uint i=0;i<*n;i++) {
    uint c;
    uint read=fread(&c,sizeof(uint),1,fp);
    //assert(read==1);
    (*text)[i] = 1+(uint)c;
    c += read;
    max_symbol = max(max_symbol,(*text)[i]);
  }
  max_symbol++;
  fclose(fp);

  /*static_sequence * ss = ssb->build(*text,*n);

  char * fname2 = new char[10+string(fname).length()];
  sprintf(fname2,"%s.wt",fname);
  fp = fopen(fname2,"w");
  ss->save(fp);
  fclose(fp);
  delete ss;
  fp = fopen(fname2,"r");
  ss = static_sequence::load(fp);
  fclose(fp);
  delete [] fname2;
  return ss;*/
}

static_sequence * savetest(char * bname, static_sequence * ss) {
	char * fname = new char[10+string(bname).length()];
	sprintf(fname,"%s.ss",bname);
	FILE * fp = fopen(fname,"w");
	cout << "Saving structure ... "; cout.flush();
	ss->save(fp);
	cout << "done" << endl; cout.flush();
	fclose(fp);
	cout << "Deleting structure ... "; cout.flush();
	delete ss;
	cout << "done" << endl; cout.flush();
	fp = fopen(fname,"r");
	cout << "Loading structure ... "; cout.flush();
	ss = static_sequence::load(fp);
	cout << "done" << endl; cout.flush();
	fclose(fp);
	if(ss==NULL) cout << "Error loading static_sequence" << endl;
	//cout << ss << endl;
	delete [] fname;
	return ss;
}

void speed_rank(static_sequence * ss, uint * text, uint n) {
  uint max_symbol = 0;
  for(uint i=0;i<n;i++) {
    max_symbol = max(max_symbol,text[i]);
  }
  max_symbol++;
  uint *occ = new uint[max_symbol];
  for(uint i=0;i<max_symbol;i++)
    occ[i] = 0;
  for(uint i=0;i<n;i++)
    occ[text[i]]++;

  uint * valid_symbols = new uint[max_symbol]; uint c=0;
  for(uint i=0;i<max_symbol;i++) {
    if(occ[i]>0)
      valid_symbols[c++]=i;
  }

  uint acc=0;
  srand(SEED);

  start_clock();
  for(uint i=0;i<NQUERIES;i++) {
    uint symb = rand()%c;
    uint pos = rand()%n;
    acc += ss->rank(valid_symbols[symb],pos);
  }
  double t = stop_clock();
  cout << " * Time for " << NQUERIES << " ranks: " << t << " secs" << endl;
  cout << " * Time per rank: " << 1000*t/NQUERIES << " msecs" << endl;
  cout << " - Check sum: " << acc << endl;
  delete [] valid_symbols;
  delete [] occ;
}

void speed_select(static_sequence * ss, uint * text, uint n) {
  uint max_symbol = 0;
  for(uint i=0;i<n;i++) {
    max_symbol = max(max_symbol,text[i]);
  }
  max_symbol++;
  uint *occ = new uint[max_symbol];
  for(uint i=0;i<max_symbol;i++)
    occ[i] = 0;
  for(uint i=0;i<n;i++)
    occ[text[i]]++;

  uint * valid_symbols = new uint[max_symbol]; uint c=0;
  for(uint i=0;i<max_symbol;i++) {
    if(occ[i]>0)
      valid_symbols[c++]=i;
  }

  uint acc=0;
  srand(SEED);

  start_clock();
  for(uint i=0;i<NQUERIES;i++) {
    uint symb = rand()%c;
    uint pos = rand()%occ[valid_symbols[symb]]+1;
    acc += ss->select(valid_symbols[symb],pos);
  }
  double t = stop_clock();
  cout << " * Time for " << NQUERIES << " selects: " << t << " secs" << endl;
  cout << " * Time per select: " << 1000*t/NQUERIES << " msecs" << endl;
  cout << " - Check sum: " << acc << endl;
  delete [] valid_symbols;
  delete [] occ;
}

void speed_access(static_sequence * ss, uint * text, uint n) {
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
