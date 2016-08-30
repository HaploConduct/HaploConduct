/* wt_coder_binary.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * wt_coder_binary definition
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
 
#include <wt_coder_binary.h>

wt_coder_binary::wt_coder_binary(uint * seq, uint n, alphabet_mapper * am) {
	uint max_v = 0;
  for(uint i=0;i<n;i++)
    max_v = max(am->map(seq[i]),max_v);
  h=bits(max_v);
}

wt_coder_binary::wt_coder_binary(uchar * seq, uint n, alphabet_mapper * am) {
	uint max_v = 0;
  for(uint i=0;i<n;i++)
      max_v = max(am->map((uint)seq[i]),max_v);
  h=bits(max_v);
}

wt_coder_binary::wt_coder_binary() {}

wt_coder_binary::~wt_coder_binary() {}

bool wt_coder_binary::is_set(uint symbol, uint l) {
	if((1<<(h-l-1))&symbol) return true;
	return false;
}

bool wt_coder_binary::done(uint symbol, uint l) {
	if(l==h) return true;
	return false;
}

uint wt_coder_binary::size() {
  return sizeof(wt_coder_binary);
}

uint wt_coder_binary::save(FILE *fp) {
  uint wr = WT_CODER_BINARY_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  wr += fwrite(&h,sizeof(uint),1,fp);
  return wr-2;
}

wt_coder_binary * wt_coder_binary::load(FILE *fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd!=WT_CODER_BINARY_HDR) return NULL;
  wt_coder_binary * ret = new wt_coder_binary();
  if(fread(&ret->h,sizeof(uint),1,fp)!=1) {
    delete ret;
    return NULL;
  }
  return ret;
}
