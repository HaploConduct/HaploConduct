/* static_bitsequence_naive.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * Naive Bitsequence - don't use, only for testing
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

#include "static_bitsequence_naive.h"

static_bitsequence_naive::static_bitsequence_naive(uint * bitseq, uint len) {
	this->len = len;
	this->bitseq = new uint[len/W+(len%W>0)];
  for(uint i=0;i<len/W+(len%W>0);i++)
    this->bitseq[i] = bitseq[i];
	uint ret = 0;
	for(uint k=0;k<len;k++)
		if(bitget(bitseq,k))
			ret++;
	this->ones = ret;
}

static_bitsequence_naive::~static_bitsequence_naive() {
	delete [] bitseq;
}

uint static_bitsequence_naive::rank1(uint i) {
  if(i>=len) return ones;
	uint ret = 0;
	for(uint k=0;k<=i;k++)
		if(bitget(bitseq,k))
			ret++;
	return ret;
}

uint static_bitsequence_naive::select1(uint i) {
	if(i==0) return (uint)-1;
	if(i>ones) return len;
	uint cnt = 0;
	for(uint k=0;k<len;k++) {
		if(bitget(bitseq,k))
			cnt++;
		if(cnt==i)
			return k;
	}
	return len;
}

bool static_bitsequence_naive::access(uint i) {
	return bitget(bitseq,i)!=0;
}

uint static_bitsequence_naive::size() {
	return BW*uint_len(len,1);
}

int static_bitsequence_naive::save(FILE * fp) { return -1; }
