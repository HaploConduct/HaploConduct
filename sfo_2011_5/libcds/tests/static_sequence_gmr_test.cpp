/* static_sequence_gmr_test.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_sequence_gmr_test
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

#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <basics.h>
#include <static_bitsequence.h>
#include <static_bitsequence_builder.h>
#include <alphabet_mapper.h>
#include <static_sequence.h>
#include <wt_coder.h>
#include <static_sequence_builder.h>
#include "static_sequence_tester.h"

int main(int argc, char ** argv) {
  if(argc!=8) {
    cout << "Usage: " << argv[0] << " <file> <b|r> <p|w|c> <sampling> <chunk_length> <perm_samp> <t|s>" << endl;
    return 0;
  }
  stringstream ss;
  ss << argv[4];
  uint samp;
  ss >> samp;
  stringstream ss2;
  ss2 << argv[5];
  uint chunk_length;
  ss2 >> chunk_length;
	stringstream ss3;
	ss3 << argv[6];
	uint perm_samp;
	ss3 >> perm_samp;
  
  uint * text;
  uint n;
  load(argv[1],&text,&n);
    
  static_bitsequence_builder * bmb;
  if(string(argv[2])==string("b"))
    bmb = new static_bitsequence_builder_brw32(samp);
  else
    bmb = new static_bitsequence_builder_rrr02(samp);
  
  static_sequence_builder * ssb;
  
	if(string(argv[3])==string("w")) {
    alphabet_mapper * am = new alphabet_mapper_cont(text,n,bmb);
		ssb = new static_sequence_builder_wvtree_noptrs(bmb,am);
	}	else if(string(argv[3])==string("p")) {
    alphabet_mapper * am = new alphabet_mapper_none();
		wt_coder * wc = new wt_coder_huff(text,n,am);
		ssb = new static_sequence_builder_wvtree(wc,bmb,am);
  } else {
    ssb = new static_sequence_builder_gmr_chunk(bmb, new static_permutation_builder_mrrr(perm_samp,bmb));
  }
  
  static_sequence * sseq = new static_sequence_gmr(text,n,chunk_length,bmb,ssb);

	delete bmb;
	delete ssb;

  sseq = savetest(argv[1], sseq);
  if(string(argv[7])==string("t"))
    test_static_sequence(text,n,sseq);
  else 
    cout << "Size: " << sseq->size() << endl;
  cout << "*************************************" << endl;
  speed_access(sseq,text,n);
  cout << "*************************************" << endl;
  speed_rank(sseq,text,n);
  cout << "*************************************" << endl;
  speed_select(sseq,text,n);
  
  delete sseq;
  delete [] text;
}

