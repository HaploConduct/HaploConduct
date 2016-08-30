/* wt_coder_huff.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * wt_coder_huff definition
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
 
#include <wt_coder_huff.h>

wt_coder_huff::wt_coder_huff(uint * symbs, uint n, alphabet_mapper * am) {
  for(uint i=0;i<n;i++)
    symbs[i] = am->map(symbs[i]);
	hc = new huffman_codes(symbs, n);
  buffer = new uint[hc->max_length()/W+1]; 
  s_len = 0; last_symbol = (uint)-1;
  for(uint i=0;i<n;i++)
    symbs[i] = am->unmap(symbs[i]);
}

wt_coder_huff::wt_coder_huff(uchar * symbs, uint n, alphabet_mapper * am) {
  for(uint i=0;i<n;i++)
    symbs[i] = (uchar)am->map((uint)symbs[i]);
	hc = new huffman_codes(symbs, n);
  buffer = new uint[hc->max_length()/W+1]; 
  s_len = 0; last_symbol = (uint)-1;
  for(uint i=0;i<n;i++)
    symbs[i] = (uchar)am->unmap((uint)symbs[i]);
}

wt_coder_huff::wt_coder_huff() {}

wt_coder_huff::~wt_coder_huff() {
  delete hc;
  delete [] buffer;
}

bool wt_coder_huff::is_set(uint symbol, uint l) {
	if(symbol!=last_symbol) {
    s_len = (uint)hc->encode(symbol, buffer, (ulong)0);
    last_symbol = symbol;
  }
  return bitget(buffer,l);
}

bool wt_coder_huff::done(uint symbol, uint l) {
  if(symbol!=last_symbol) {
    s_len = (uint)hc->encode(symbol, buffer, (ulong)0);
    last_symbol = symbol;
  }
  return l==s_len;
}

uint wt_coder_huff::size() {
  return 2*sizeof(uint)+sizeof(wt_coder_huff)+hc->size()+(hc->max_length()/W+1)*sizeof(uint);
}

uint wt_coder_huff::save(FILE * fp) {
  uint wr = WT_CODER_HUFF_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  if(wr!=1) return 1;
  if(hc->save(fp)) return 1;
  //if(am->save(fp)) return 1;
  return 0;
}

wt_coder_huff * wt_coder_huff::load(FILE *fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd!=WT_CODER_HUFF_HDR) return NULL;
  wt_coder_huff * ret = new wt_coder_huff();
  ret->hc = huffman_codes::load(fp);
  ret->buffer = new uint[ret->hc->max_length()/W+1]; 
  ret->s_len = 0; ret->last_symbol = (uint)-1;
  return ret;
}
