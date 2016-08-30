/* alphabet_mapper_cont.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * alphabet_mapper_cont definition
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
 
#include <alphabet_mapper_cont.h>

alphabet_mapper_cont::alphabet_mapper_cont(uint * seq, uint n, static_bitsequence_builder *bmb) { 
  uint max_v = 0;
  for(uint i=0;i<n;i++)
    max_v = max(max_v,seq[i]);
  max_v++;
  uint blen = uint_len(max_v,1);
  uint * bmap = new uint[blen];
  for(uint i=0;i<blen;i++)
    bmap[i] = 0;
  for(uint i=0;i<n;i++)
    bitset(bmap,seq[i]);
  m = bmb->build(bmap,max_v);
  delete [] bmap;
}

alphabet_mapper_cont::alphabet_mapper_cont(unsigned char * seq, uint n, static_bitsequence_builder *bmb) { 
  uint max_v = 0;
  for(uint i=0;i<n;i++)
    max_v = max(max_v, (uint)seq[i]);
  max_v++;
  uint blen = uint_len(max_v,1);
  uint * bmap = new uint[blen];
  for(uint i=0;i<blen;i++)
    bmap[i] = 0;
  for(uint i=0;i<n;i++)
    bitset(bmap, (uint)seq[i]);
  m = bmb->build(bmap,max_v);
  delete [] bmap;
}

alphabet_mapper_cont::alphabet_mapper_cont() {
}

alphabet_mapper_cont::~alphabet_mapper_cont() {
  delete m;
}

uint alphabet_mapper_cont::map(uint s) {
  return m->rank1(s);
}

uint alphabet_mapper_cont::unmap(uint s) {
  return m->select1(s);
}

uint alphabet_mapper_cont::size() { 
  return sizeof(alphabet_mapper_cont)+m->size(); 
}

uint alphabet_mapper_cont::save(FILE *fp) {
  uint wr = ALPHABET_MAPPER_CONT_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  if(wr!=1) return 1;
  if(m->save(fp)) return 1;
  return 0;
}

alphabet_mapper_cont * alphabet_mapper_cont::load(FILE * fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd!=ALPHABET_MAPPER_CONT_HDR) return NULL;
  alphabet_mapper_cont * ret = new alphabet_mapper_cont();
  ret->m = static_bitsequence::load(fp);
  if(ret->m==NULL) {
    delete ret;
    return NULL;
  }
  return ret;
}
