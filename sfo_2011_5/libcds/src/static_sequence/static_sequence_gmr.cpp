/* static_sequence_gmr.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * GMR
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

#include <static_sequence_gmr.h>

static_sequence_gmr::static_sequence_gmr(uint * sequence, uint n, uint chunk_length, static_bitsequence_builder * bmb, static_sequence_builder * ssb) {
  len = n;
	if(len%chunk_length) len+=chunk_length-len%chunk_length;
	uint * new_seq = new uint[len];
  sigma = 0;
	for(uint i=0;i<n;i++){
		new_seq[i] = sequence[i]+1;
    sigma = max(sigma,new_seq[i]);
	}
  sigma++;
	for(uint i=n;i<len;i++)
		new_seq[i] = sigma;
	if(len!=n) sigma++;
  this->chunk_length = chunk_length;
  build(new_seq,bmb,ssb);
	delete [] new_seq;
}

static_sequence_gmr::static_sequence_gmr() {
}

static_sequence_gmr::~static_sequence_gmr() {
  delete B;
  for (uint i=0;i<len/chunk_length;i++)
    delete chunk[i];
  delete [] chunk;
}


void static_sequence_gmr::build(uint * sequence, static_bitsequence_builder * bmb, static_sequence_builder * ssb) {
  uint num_chunks = len/chunk_length;
  chunk = new static_sequence*[num_chunks];
  assert(chunk!=NULL);
  for (uint i=0;i<num_chunks;i++) {
    chunk[i] = ssb->build(sequence+i*chunk_length, chunk_length);
    //cout << "1." << i << endl; cout.flush();
    assert(chunk[i]!=NULL);
  }
  uint * ones = get_ones(sequence);
  uint *B_bitmap = new uint[(2+len+(unsigned long long)num_chunks*sigma)/W+1];
    assert(B_bitmap!=NULL);
  for (uint i=0;i<(2+len+(unsigned long long)num_chunks*sigma)/W+1;i++)
    B_bitmap[i] = 0;
  uint pos=0;
  for (unsigned long long i=0;i<(unsigned long long)num_chunks*sigma;i++) {
    for (uint j=0;j<ones[i];j++) {
      bitset(B_bitmap, pos);
      pos++;
    }
    pos++;
  }
  pos++;
  //cout << "5  pos=" << pos << endl; cout.flush();
  B = bmb->build(B_bitmap, pos);
  //cout << "6" << endl; cout.flush();
  delete [] B_bitmap;
  delete [] ones;
}


uint * static_sequence_gmr::get_ones(uint * sequence) {
  uint * ones = new uint[(unsigned long long)(len/chunk_length)*sigma];
  assert(ones!=NULL);
  for (uint i=0;i<(unsigned long long)(len/chunk_length)*sigma;i++) ones[i] = 0;
  for (uint i=0;i<len;i++) {
    uint whichChunk = (uint)(((unsigned long long)sequence[i]*len+i)/chunk_length);
    ones[whichChunk]++;
  }
  return ones;
}


uint static_sequence_gmr::rank(uint c, uint j) {
  c++;
  uint i = j/chunk_length;
  uint bp = (c)*(len/chunk_length);
  uint rank_pos = B->select0(bp);
  uint prev = rank_pos-bp+1;
  uint sum = B->rank1(B->select0(bp+i)) - prev;
  uint cr = chunk[i]->rank(c,j-i*chunk_length);
	/*if(c==0) {
		cout << "c=" << c << " j=" << j << endl;
		cout << "i=" << i << endl;
		cout << "bp=" << bp << endl;
		cout << "rank_pos=" << rank_pos << endl;
		cout << "prev=" << prev << endl;
		cout << "sum=" << sum << endl;
		cout << "cr=" << cr << endl;
	}*/
  return sum + cr;
}


uint static_sequence_gmr::select(uint c, uint j) {
   c++;
  uint rank_pos = B->select0(c*(len/chunk_length));
  uint prev = B->rank1(rank_pos);
  uint sel = prev+j;
  uint block = (B->select1(sel));
  uint i = block-sel+1;
  uint desp = B->rank1(B->select0((i)))-prev;
  if (desp+1==0) desp=0;
  uint rchunk = i%(len/chunk_length);
	/*if(j==90) {
		cout << "------------------------------" << endl;
		cout << "c=" << c << "  j=" << j << endl;
		cout << "chunk_length=" << chunk_length << endl;
		cout << "rank_pos=" << rank_pos << endl;
		cout << "prev=" << prev << endl;
		cout << "sel=" << sel << endl;
		cout << "block=" << block << endl;
		cout << "i=" << i << endl;
		cout << "desp=" << desp << endl;
		cout << "rchunk=" << rchunk << endl;
		cout << "j-desp=" << j-desp << endl;
	}*/
  return (rchunk*chunk_length)+chunk[rchunk]->select(c, j-desp);
}


uint static_sequence_gmr::access(uint j) {
  return chunk[j/chunk_length]->access(j%chunk_length)-1;
}


uint static_sequence_gmr::size() {
  uint s = 0;
  for (uint i=0;i<len/chunk_length;i++)
    s += sizeof(void*)+chunk[i]->size();
  return s+B->size()+sizeof(static_sequence_gmr);
}

uint static_sequence_gmr::save(FILE *fp) {
  uint wr = GMR_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  wr += fwrite(&len,sizeof(uint),1,fp);
  wr += fwrite(&sigma,sizeof(uint),1,fp);
  wr += fwrite(&chunk_length,sizeof(uint),1,fp);
  if(wr!=4) return 1;
  if(B->save(fp)) return 1;
  for(uint i=0;i<len/chunk_length;i++)
    if(chunk[i]->save(fp)) return 1;
  return 0;
}

static_sequence_gmr * static_sequence_gmr::load(FILE *fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd!=GMR_HDR) return NULL;
  static_sequence_gmr * ret = new static_sequence_gmr();
  rd = fread(&ret->len,sizeof(uint),1,fp);
  rd += fread(&ret->sigma,sizeof(uint),1,fp);
  rd += fread(&ret->chunk_length,sizeof(uint),1,fp);
  if(rd!=3) {
    delete ret;
    return NULL;
  }
  ret->B = static_bitsequence::load(fp);
  if(ret->B==NULL) {
    delete ret;
    return NULL;
  }
  ret->chunk = new static_sequence*[ret->len/ret->chunk_length];
  for(uint i=0;i<ret->len/ret->chunk_length;i++) {
    ret->chunk[i] = static_sequence::load(fp);
    if(ret->chunk[i]==NULL) {
      delete ret;
      return NULL;
    }
  }
  return ret;
}
