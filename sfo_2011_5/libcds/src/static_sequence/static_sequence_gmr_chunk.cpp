/* static_sequence_gmr_chunk.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * gmr_chunk
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
 
#include "static_sequence_gmr_chunk.h"

static_sequence_gmr_chunk::static_sequence_gmr_chunk(uint * sequence, uint chunk_length, static_bitsequence_builder *bmb, static_permutation_builder *pmb) {
  sigma = 0;
  for(uint i=0;i<chunk_length;i++) {
    sigma = max(sigma,sequence[i]);
  }sigma++;
  uint * X_bitmap = new uint[uint_len(1+chunk_length+sigma,1)];
  assert(X_bitmap!=NULL);
  for(uint i=0;i<uint_len(1+chunk_length+sigma,1);i++) X_bitmap[i]=0;
  uint pi_blen = bits(chunk_length-1);
  uint * pi = new uint[uint_len(pi_blen,chunk_length)];
  assert(pi!=NULL);
  for(uint i=0;i<uint_len(pi_blen,chunk_length);i++) pi[i] = 0;
  uint X_pos = 0;
  uint * counter = new uint[sigma+2];
  for(uint c=0;c<=sigma+1;c++) counter[c]=0;
  for(uint i=0;i<chunk_length;i++) counter[sequence[i]+1]++;

  for(uint c=0;c<sigma;c++) {
    X_pos++;
    for(uint i=0;i<counter[c+1];i++) {
      bitset(X_bitmap, X_pos);
      X_pos++;
    }
    counter[c+1]+=counter[c];
  }
  X_pos++;
  for(uint i=0;i<chunk_length;i++) {
    set_field(pi, pi_blen,counter[sequence[i]], i);
    counter[sequence[i]]++;
  }
  //cout << "pi_blen=" << pi_blen << endl;
  this->X = bmb->build(X_bitmap,X_pos); //new BitRankW32Int(X_bitmap, X_pos, true,20);
  assert(X!=NULL);
  delete [] X_bitmap;
  //cout << "a" << endl; cout.flush();
  this->permutation = pmb->build(pi,chunk_length); //createPerm(pi, chunk_length, t);
  //cout << "a" << endl; cout.flush();
  assert(permutation!=NULL);
  this->sigma = sigma;
  this->len = chunk_length;
	delete [] counter;
}

static_sequence_gmr_chunk::static_sequence_gmr_chunk() {
}

static_sequence_gmr_chunk::~static_sequence_gmr_chunk() {
	delete X;
  delete permutation;
}


uint static_sequence_gmr_chunk::access(uint j) {
  uint invPerm = permutation->rev_pi(j); //inversePerm(permutation, j);
  //cout << "invPerm=" << invPerm << endl;
  uint rank_pos = X->select1(invPerm+1);
  //cout << "rank_pos=" << rank_pos << endl;
  uint ret = rank_pos - X->rank1(rank_pos);// - 1;
  //cout << "ret = " << ret << endl;
  return ret;
}


uint static_sequence_gmr_chunk::select(uint i, uint j) {
  uint pos = X->select0(i+1) + j - i -1;
  /*cout << "pos=" << pos << endl;
  cout << "pos'=" << X->rank1(X->select0(i+1)+j) << endl;
  cout << "perm_pos=" << permutation->pi(pos) << endl;*/
  return permutation->pi(pos); //getelemPerm(permutation, pos);
}


uint static_sequence_gmr_chunk::rank(uint i, uint j) {
  uint ini = X->select0(i+1)-i;
  uint ini_o = ini;
  uint fin = X->select0(i+2);
	if(fin<i+2) return 0;
	fin = fin-(i+2);
	if(fin<ini) return 0;
	if(permutation->pi(ini) > j) return 0;
	if(permutation->pi(ini) == j) return 1;
	if(ini==fin) return 1;
  while(ini < fin-1) {
		uint med = (ini+fin)/2;
    uint elem = permutation->pi(med); //getelemPerm(permutation, med);
    if(elem >= j) fin = med;
    else ini = med;
  }
	while(fin>ini_o && permutation->pi(fin)>j) fin--;
  return fin-ini_o+1;
}


uint static_sequence_gmr_chunk::size() {
  return sizeof(static_sequence_gmr_chunk)+permutation->size()+X->size();
}

uint static_sequence_gmr_chunk::save(FILE *fp) {
  uint wr = GMR_CHUNK_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  wr += fwrite(&len,sizeof(uint),1,fp);
  wr += fwrite(&sigma,sizeof(uint),1,fp);
  if(wr!=3) return 1;
  if(X->save(fp)) return 1;
  if(permutation->save(fp)) return 1;
  return 0;
}

static_sequence_gmr_chunk * static_sequence_gmr_chunk::load(FILE *fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd!=GMR_CHUNK_HDR) return NULL;
  static_sequence_gmr_chunk * ret = new static_sequence_gmr_chunk();
  rd = fread(&ret->len,sizeof(uint),1,fp);
  rd += fread(&ret->sigma,sizeof(uint),1,fp);
  ret->X = static_bitsequence::load(fp);
  ret->permutation = static_permutation::load(fp);
  if(rd!=2 || ret->X==NULL || ret->permutation==NULL) {
		/*cout << "rd=" << rd << endl;
		cout << "X =" << ret->X << endl;
		cout << "P =" << ret->permutation << endl;*/
    delete ret;
    return NULL;
  }
  return ret;
}
