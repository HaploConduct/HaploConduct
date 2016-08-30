/* static_sequence.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_sequence definition
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
 
#include <static_sequence.h>

static_sequence::static_sequence() {}
static_sequence::~static_sequence() {}
uint static_sequence::length() { return len; }

uint static_sequence::count(uint s) {
  return rank(s,len-1);
}

static_sequence * static_sequence::load(FILE * fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  fseek(fp,-sizeof(uint),SEEK_CUR);
  switch(rd) {
    case WVTREE_HDR: return static_sequence_wvtree::load(fp);
    case GMR_CHUNK_HDR: return static_sequence_gmr_chunk::load(fp);
    case GMR_HDR: return static_sequence_gmr::load(fp);
    case WVTREE_NOPTRS_HDR: return static_sequence_wvtree_noptrs::load(fp);
		case BS_HDR: return static_sequence_bs::load(fp);
  }
  return NULL;
}

uint static_sequence::select_next(uint c, uint i) { 
	return select(c,rank(c,i)+1);
}

bool static_sequence::test(uint * seq, uint n) {
	uint sigma = 0;
	for(uint i=0;i<n;i++)
		sigma = max(sigma,seq[i]);
	uint * occ = new uint[sigma+1];
	for(uint i=0;i<=sigma;i++)
		occ[i] = 0;
	for(uint i=0;i<n;i++) {
		occ[seq[i]]++;
		if(rank(seq[i],i)!=occ[seq[i]]) {
			cout << "rank failed!" << endl;
			cout << "rank("<<seq[i]<<","<<i<<")="<<rank(seq[i],i)<<endl;
			cout << "expected result: " << occ[seq[i]] << endl;
			delete [] occ;
			return false;
		}
		if(i>0 && rank(seq[i],i-1)!=occ[seq[i]]-1) {
			cout << "rank-1 failed!" << endl;
			delete [] occ;
			return false;
		}
		if(select(seq[i],occ[seq[i]])!=i) {
			cout << "select failed!" << endl;
			cout << "select(" << seq[i] << "," << occ[seq[i]] << ")="<<select(seq[i],occ[seq[i]]) << endl;
			cout << "i=" << i << "  rank(" << seq[i] << ",i)=" << rank(seq[i],i) << endl;
			delete [] occ;
			return false;
		}
		if(access(i)!=seq[i]) {
			cout << "access failed!" << endl;
			delete [] occ;
			return false;
		}
	}
	delete [] occ;
	return true;
}

