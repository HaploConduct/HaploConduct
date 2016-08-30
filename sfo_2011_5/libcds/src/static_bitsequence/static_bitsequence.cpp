/* static_bitsequence.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_bitsequence definition
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

#include "static_bitsequence.h"

uint static_bitsequence::rank0(uint i) {
	return i+1-rank1(i);
}

uint static_bitsequence::rank1(uint i) {
  if(i>=len) return (uint)-1;
  if(ones==0) return 0;
	if(ones==len) return i+1;
  uint ini = 1;
  uint fin = ones;
  while(ini<fin) {
    uint pos = (ini+fin)/2;
    uint bp = select1(pos);
    if(bp==i) return pos;
    if(bp<i)
      ini = pos+1;
    else
      fin = pos-1;
  }
	if(select1(ini)>i) return ini-1;
	return ini;
}

uint static_bitsequence::select0(uint i) {
  if(i>len-ones) return -1;
  if(i==0) return -1;
	if(ones==0) return i-1;
  uint ini = 0;
  uint fin = len-1;
  while(ini<fin) {
    uint pos = (ini+fin)/2;
    uint br = rank0(pos);
    if(br<i)
      ini = pos+1;
    else
      fin = pos;
  }
	return ini;
}

uint static_bitsequence::select1(uint i) {
  if(i>ones) return -1;
  if(i==0) return -1;
	if(ones==len) return i-1;
  uint ini = 0;
  uint fin = len-1;
  while(ini<fin) {
    uint pos = (ini+fin)/2;
    uint br = rank1(pos);
    if(br<i)
      ini = pos+1;
    else
      fin = pos;
  }
	return ini;
}

uint static_bitsequence::select_next1(uint i) {
	return select1(rank1(i)+1);
}

uint static_bitsequence::select_next0(uint i) {
	return select0(rank0(i)+1);
}

bool static_bitsequence::access(uint i) {
  return (rank1(i)-(i!=0?rank1(i-1):0))>0;
}

uint static_bitsequence::length() {
	return len;
}

uint static_bitsequence::count_one() {
	return ones;
}

uint static_bitsequence::count_zero() {
	return len-ones;
}

static_bitsequence * static_bitsequence::load(FILE * fp) {
  uint r;
  if(fread(&r,sizeof(uint),1,fp)!=1) return NULL;
  fseek(fp,-1*sizeof(uint),SEEK_CUR);
  switch(r) {
    case RRR02_HDR: return static_bitsequence_rrr02::load(fp);
    case BRW32_HDR: return static_bitsequence_brw32::load(fp);
    case RRR02_LIGHT_HDR: return static_bitsequence_rrr02_light::load(fp);
    case SDARRAY_HDR: return static_bitsequence_sdarray::load(fp);
  }
  return NULL;
}
