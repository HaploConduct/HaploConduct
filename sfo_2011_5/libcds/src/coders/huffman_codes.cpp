/* huffman_codes.cpp
   Copyright (C) 2008, Francisco Claude, all rights reserved.

   Wrapper for huff written by Gonzalo Navarro

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#include <huffman_codes.h>

huffman_codes::huffman_codes(uint * symb, uint n) {
  uint max_v = 0;
  for(uint i=0;i<n;i++)
    max_v = max(max_v,symb[i]);
  uint * occ = new uint[max_v+1];
  for(uint i=0;i<max_v+1;i++)
    occ[i] = 0;
  for(uint i=0;i<n;i++)
    occ[symb[i]]++;
  huff_table = createHuff(occ, max_v);
  delete [] occ;
}

huffman_codes::huffman_codes(uchar * symb, uint n) {
  uchar max_v = 0;
  for(uint i=0;i<n;i++)
    max_v = max(max_v,symb[i]);
  uint * occ = new uint[max_v+1];
  for(uint i=0;i<(uint)max_v+1;i++)
    occ[i] = 0;
  for(uint i=0;i<n;i++)
    occ[symb[i]]++;
  huff_table = createHuff(occ, max_v);
  delete [] occ;
}

huffman_codes::huffman_codes() {
}

huffman_codes::~huffman_codes() {
  freeHuff(huff_table);
}

uint huffman_codes::max_length() {
  return huff_table.depth;
}

uint huffman_codes::size() {
  return sizeof(huffman_codes)+sizeHuff(huff_table);
}

ulong huffman_codes::encode(uint symb, uint * stream, ulong pos) {
  return encodeHuff(huff_table, symb, stream, pos);
}

ulong huffman_codes::decode(uint * symb, uint * stream, ulong pos) {
  return decodeHuff(huff_table, symb, stream, pos);
}

uint huffman_codes::save(FILE *fp) {
  saveHuff(huff_table,fp);
  return 0;
}

huffman_codes * huffman_codes::load(FILE * fp) {
  huffman_codes * ret = new huffman_codes();
  ret->huff_table = loadHuff(fp,1);
  return ret;
}
