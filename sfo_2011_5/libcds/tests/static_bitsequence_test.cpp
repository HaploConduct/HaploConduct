/* static_bitsequence_test.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_bitsequence_test
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
#include <sstream>
#include <basics.h>
#include <static_bitsequence.h>
#include "static_bitsequence_tester.h"

using namespace std;

int main(int argc, char ** argv) {
	if(argc!=5) {
		cout << "usage: " << argv[0] << " <bitmap_file> <b|r|s> <sample_rate> <t|s>" << endl;
		return 0;
	}
	FILE * fp = fopen(argv[1],"r");
	if(fp==NULL) {
		cout << "Error opening " << argv[1] << endl;
		return -1;
	}
	uint *bitseq, len;//, ones;
	uint l=fread(&len, sizeof(uint), 1, fp);
	//l += fread(&ones,sizeof(uint),1,fp);
	bitseq = new uint[uint_len(len,1)];
	l+=fread(bitseq, sizeof(uint), uint_len(len,1), fp);
	fclose(fp);

  uint sample_rate;
  stringstream ss(argv[3]);
  ss >> sample_rate;
  
	static_bitsequence * bs;
  
  if(string(argv[2])==string("r")) bs = new static_bitsequence_rrr02(bitseq,len,sample_rate);
  if(string(argv[2])==string("s")) bs = new static_bitsequence_sdarray(bitseq,len);
  else bs = new static_bitsequence_brw32(bitseq,len,sample_rate);
	
  cout << "Size: " << bs->size() << endl;
  cout << "bpb = " << bs->size()*8./len << endl;
  
	/*for(uint kk=0;kk<30;kk++)
		cout << bs->access(kk);
	cout << endl;*/

	/*for(uint kk=0;kk<20;kk++) {
		bs->select_next1(kk);
	}*/

	if(string(argv[4])==string("t"))
	  test_bitsequence(bitseq,len,bs);
  cout << "******************************************" << endl;
  speed_access(bs, bitseq, len);
  cout << "******************************************" << endl;
  speed_rank0(bs, bitseq, len);
  cout << "******************************************" << endl;
  speed_rank1(bs, bitseq, len);
  cout << "******************************************" << endl;
  speed_select0(bs, bitseq, len);
  cout << "******************************************" << endl;
  speed_select1(bs, bitseq, len);
  cout << "******************************************" << endl;
  speed_selectnext1(bs, bitseq, len);
}
