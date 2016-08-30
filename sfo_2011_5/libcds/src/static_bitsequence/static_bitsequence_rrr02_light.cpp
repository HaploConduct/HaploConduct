/* static_bitsequence_rrr02_light.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_bitsequence_rrr02_light definition
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

#include <static_bitsequence_rrr02_light.h>

#define VARS_NEEDED uint C_len = len/BLOCK_SIZE_LIGHT + (len%BLOCK_SIZE_LIGHT!=0);\
uint C_field_bits = bits(BLOCK_SIZE_LIGHT);\
uint O_len = uint_len(1,O_bits_len);\
uint C_sampling_len = C_len/sample_rate+2;\
uint C_sampling_field_bits = bits(ones);\
uint O_pos_len = C_len/sample_rate+1;\
uint O_pos_field_bits = bits(O_bits_len);


table_offset * static_bitsequence_rrr02_light::E = NULL;

static_bitsequence_rrr02_light::static_bitsequence_rrr02_light() {
  ones=0;
  len=0;
  if(E==NULL) E = new table_offset(BLOCK_SIZE_LIGHT);
  E->use();
  C = NULL;
  O = NULL;
  C_sampling = NULL;
  O_pos = NULL;
  sample_rate = DEFAULT_SAMPLING_LIGHT;
  O_bits_len = 0;
}

static_bitsequence_rrr02_light::static_bitsequence_rrr02_light(uint * bitseq, uint len, uint sample_rate) {
  ones = 0;
  this->len = len;
  if(E==NULL) E = new table_offset(BLOCK_SIZE_LIGHT);
  E->use();
  // Table C
  uint C_len = len/BLOCK_SIZE_LIGHT + (len%BLOCK_SIZE_LIGHT!=0);
  uint C_field_bits = bits(BLOCK_SIZE_LIGHT);
  C = new uint[uint_len(C_len,C_field_bits)];
  for(uint i=0;i<uint_len(C_len,C_field_bits);i++)
    C[i] = 0;
  O_bits_len = 0;
  for(uint i=0;i<C_len;i++) {
    uint value = popcount(get_var_field(bitseq,i*BLOCK_SIZE_LIGHT,min((uint)len-1,(i+1)*BLOCK_SIZE_LIGHT-1)));
    assert(value<=BLOCK_SIZE_LIGHT);
    set_field(C,C_field_bits,i,value);
    ones += value;
    O_bits_len += E->get_log2binomial(BLOCK_SIZE_LIGHT,value);
  }
  // Table O
  uint O_len = uint_len(1,O_bits_len);
  O = new uint[O_len];
  for(uint i=0;i<O_len;i++)
    O[i] = 0;
  uint O_pos = 0;
  for(uint i=0;i<C_len;i++) {
    uint value = (ushort)get_var_field(bitseq,i*BLOCK_SIZE_LIGHT,min((uint)len-1,(i+1)*BLOCK_SIZE_LIGHT-1));
    set_var_field(O,O_pos,O_pos+E->get_log2binomial(BLOCK_SIZE_LIGHT,popcount(value))-1,E->compute_offset((ushort)value));
    O_pos += E->get_log2binomial(BLOCK_SIZE_LIGHT,popcount(value));
  }
  C_sampling = NULL;
  this->O_pos = NULL;
  
  create_sampling(sample_rate);
}

void static_bitsequence_rrr02_light::create_sampling(uint sample_rate) {
  this->sample_rate = sample_rate;
/*  for(uint i=0;i<C_len;i++) {
    O_bits_len += E->get_log2binomial(BLOCK_SIZE_LIGHT,get_field(C,C_field_bits,i));
  }*/
  // Sampling for C
  uint C_len = len/BLOCK_SIZE_LIGHT + (len%BLOCK_SIZE_LIGHT!=0);
  uint C_field_bits = bits(BLOCK_SIZE_LIGHT);
  uint C_sampling_len = C_len/sample_rate+2;
  uint C_sampling_field_bits = bits(ones);
  if(C_sampling!=NULL) delete [] C_sampling;
  C_sampling = new uint[max((uint)1,uint_len(C_sampling_len,C_sampling_field_bits))];
  for(uint i=0;i<max((uint)1,uint_len(C_sampling_len,C_sampling_field_bits));i++)
    C_sampling[i] = 0;
  uint sum = 0;
  for(uint i=0;i<C_len;i++) {
    if(i%sample_rate==0)
      set_field(C_sampling,C_sampling_field_bits,i/sample_rate,sum);
    sum += get_field(C,C_field_bits,i);
  }
  for(uint i=(C_len-1)/sample_rate+1;i<C_sampling_len;i++)
    set_field(C_sampling,C_sampling_field_bits,i,sum);
  // Sampling for O (table S) (Code separated from previous construction for readability)
  uint O_pos_len = C_len/sample_rate+1;
  uint O_pos_field_bits = bits(O_bits_len);
  if(O_pos!=NULL) delete [] O_pos;
  O_pos = new uint[uint_len(O_pos_len,O_pos_field_bits)];
  for(uint i=0;i<uint_len(O_pos_len,O_pos_field_bits);i++)
    O_pos[i] = 0;
  uint pos = 0;
  for(uint i=0;i<C_len;i++) {
    if(i%sample_rate==0)
      set_field(O_pos,O_pos_field_bits,i/sample_rate,pos);
    pos += E->get_log2binomial(BLOCK_SIZE_LIGHT,get_field(C,C_field_bits,i));
  }
}

bool static_bitsequence_rrr02_light::access(uint i) {
  uint C_field_bits = bits(BLOCK_SIZE_LIGHT);
  uint O_pos_field_bits = bits(O_bits_len);
  uint nearest_sampled_value = i/BLOCK_SIZE_LIGHT/sample_rate;
  uint pos_O = get_field(O_pos,O_pos_field_bits,nearest_sampled_value);
  uint pos = i/BLOCK_SIZE_LIGHT;
  for(uint k=nearest_sampled_value*sample_rate;k<pos;k++) {
    uint aux = get_field(C,C_field_bits,k);
    pos_O += E->get_log2binomial(BLOCK_SIZE_LIGHT,aux);
  }
  uint c = get_field(C,C_field_bits,pos);
  return ((1<<(i%BLOCK_SIZE_LIGHT))&E->short_bitmap(c,get_var_field(O,pos_O,pos_O+E->get_log2binomial(BLOCK_SIZE_LIGHT,c)-1)))!=0;
}

uint static_bitsequence_rrr02_light::rank0(uint i) {
  if(i+1==0) return 0;
  return 1+i-rank1(i);
}

uint static_bitsequence_rrr02_light::rank1(uint i) {
  uint C_field_bits = bits(BLOCK_SIZE_LIGHT);
  uint C_sampling_field_bits = bits(ones);
  uint O_pos_field_bits = bits(O_bits_len);
  if(i+1==0) return 0;
  uint nearest_sampled_value = i/BLOCK_SIZE_LIGHT/sample_rate;
  uint sum = get_field(C_sampling,C_sampling_field_bits,nearest_sampled_value);
  uint pos_O = get_field(O_pos,O_pos_field_bits,nearest_sampled_value);
  uint pos = i/BLOCK_SIZE_LIGHT;
  uint k=nearest_sampled_value*sample_rate;
  if(k%2==1 && k<pos) {
    uint aux = get_field(C,C_field_bits,k);
    sum += aux;
    pos_O += E->get_log2binomial(BLOCK_SIZE_LIGHT,aux);
    k++;
  }
  uchar * a = (uchar *)C;
  uint mask = 0x0F;
  a += k/2;
  while(k<(uint)max(0,(int)pos-1)) {
    assert(((*a)&mask)==get_field(C,C_field_bits,k));
    assert((*a)/16==get_field(C,C_field_bits,k+1));
    sum += ((*a)&mask)+(*a)/16;
    pos_O += E->get_log2binomial(BLOCK_SIZE_LIGHT,((*a)&mask))+E->get_log2binomial(BLOCK_SIZE_LIGHT,((*a)/16));
    a++;
    k+=2;
  }
  if(k<pos) {
    uint aux = get_field(C,C_field_bits,k);
    sum += aux;
    pos_O += E->get_log2binomial(BLOCK_SIZE_LIGHT,aux);
    k++;
  }
  uint c = get_field(C,C_field_bits,pos);
  sum += popcount(((2<<(i%BLOCK_SIZE_LIGHT))-1) & E->short_bitmap(c,get_var_field(O,pos_O,pos_O+E->get_log2binomial(BLOCK_SIZE_LIGHT,c)-1)));
  return sum;
}

uint static_bitsequence_rrr02_light::select0(uint i) {
  uint C_len = len/BLOCK_SIZE_LIGHT + (len%BLOCK_SIZE_LIGHT!=0);
  uint C_field_bits = bits(BLOCK_SIZE_LIGHT);
  uint C_sampling_len = C_len/sample_rate+2;
  uint C_sampling_field_bits = bits(ones);
  uint O_pos_field_bits = bits(O_bits_len);
  if(i==0) return -1;
  if(i>len-ones) return len;
  // Search over partial sums
  uint start=0;
  uint end=C_sampling_len-1;
  uint med, acc=0, pos;
  while(start<end-1) {
    med = (start+end)/2;
    acc = med*sample_rate*BLOCK_SIZE_LIGHT-get_field(C_sampling,C_sampling_field_bits,med);
    if(acc<i) {
      if(med==start) break;
      start=med;
    }
    else  {
      if(end==0) break;
      end = med-1;
    }
  }
  acc = get_field(C_sampling,C_sampling_field_bits,start);
  while(start<C_len-1 && acc+sample_rate*BLOCK_SIZE_LIGHT==get_field(C_sampling,C_sampling_field_bits,start+1)) {
    start++;
    acc +=sample_rate*BLOCK_SIZE_LIGHT;
  }
  acc = start*sample_rate*BLOCK_SIZE_LIGHT-acc;
  pos = (start)*sample_rate;
  uint pos_O = get_field(O_pos,O_pos_field_bits,start);
  // Sequential search over C
  uint s = 0;
  for(;pos<C_len;pos++) {
    s = get_field(C,C_field_bits,pos);
    if(acc+BLOCK_SIZE_LIGHT-s>=i) break;
    pos_O += E->get_log2binomial(BLOCK_SIZE_LIGHT,s);
    acc += BLOCK_SIZE_LIGHT-s;
  }
  pos = (pos)*BLOCK_SIZE_LIGHT;
  // Search inside the block
  
  while(acc<i) {
    uint new_posO = pos_O+E->get_log2binomial(BLOCK_SIZE_LIGHT,s);
    uint block = E->short_bitmap(s,get_var_field(O,pos_O,new_posO-1));
    pos_O = new_posO;
    new_posO = 0;
    while(acc<i && new_posO<BLOCK_SIZE_LIGHT) {
      pos++;new_posO++;
      acc += (((block&1)==0)?1:0);
      block = block/2;
    }
  }
  pos--;
  assert(acc==i);
  assert(rank0(pos)==i);
  assert(!access(pos));
  return pos;
}

uint static_bitsequence_rrr02_light::select1(uint i) {
  uint C_len = len/BLOCK_SIZE_LIGHT + (len%BLOCK_SIZE_LIGHT!=0);
  uint C_field_bits = bits(BLOCK_SIZE_LIGHT);
  uint C_sampling_len = C_len/sample_rate+2;
  uint C_sampling_field_bits = bits(ones);
  uint O_pos_field_bits = bits(O_bits_len);
  if(i==0) return -1;
  if(i>ones) return len;
  // Search over partial sums
  uint start=0;
  uint end=C_sampling_len-1;
  uint med, acc=0, pos;
  while(start<end-1) {
    med = (start+end)/2;
    acc = get_field(C_sampling,C_sampling_field_bits,med);
    if(acc<i) {
      if(med==start) break;
      start=med;
    }
    else  {
      if(end==0) break;
      end = med-1;
    }
  }
  acc = get_field(C_sampling,C_sampling_field_bits,start);
  while(start<C_len-1 && acc==get_field(C_sampling,C_sampling_field_bits,start+1)) start++;
  pos = (start)*sample_rate;
  uint pos_O = get_field(O_pos,O_pos_field_bits,start);
  acc = get_field(C_sampling,C_sampling_field_bits,start);
  // Sequential search over C
  uint s = 0;
  for(;pos<C_len;pos++) {
    s = get_field(C,C_field_bits,pos);
    if(acc+s>=i) break;
    pos_O += E->get_log2binomial(BLOCK_SIZE_LIGHT,s);
    acc += s;
  }
  pos = (pos)*BLOCK_SIZE_LIGHT;
  //cout << "pos=" << pos << endl;
  // Search inside the block
  while(acc<i) {
    uint new_posO = pos_O+E->get_log2binomial(BLOCK_SIZE_LIGHT,s);
    uint block = E->short_bitmap(s,get_var_field(O,pos_O,new_posO-1));
    pos_O = new_posO;
    new_posO = 0;
    while(acc<i && new_posO<BLOCK_SIZE_LIGHT) {
      pos++;new_posO++;
      acc += (((block&1)!=0)?1:0);
      block = block/2;
    }
    //cout << "i=" << i << " len=" << len << " ones=" << ones << " pos=" << pos << " acc=" << acc << " rank=" << rank1(pos) << endl;
  }
  pos--;
  assert(acc==i);
  assert(rank1(pos)==i);
  assert(access(pos));
  return pos;
}

uint static_bitsequence_rrr02_light::size() {
  VARS_NEEDED
  /*cout << "RRR02 SIZE: " << endl;
  cout << "Default: " << 9*sizeof(uint)+sizeof(uint*)*4 << endl;
  cout << "Cs:      " << uint_len(C_len,C_field_bits)*sizeof(uint) << endl;
  cout << "Os:      " << O_len*sizeof(uint) << endl;
  cout << "CSamp:   " << uint_len(C_sampling_len,C_sampling_field_bits)*sizeof(uint) << endl;
  cout << "OSamp:   " << uint_len(O_pos_len,O_pos_field_bits)*sizeof(uint) << endl;
  cout << "E:       " << E->size() << endl;*/
  uint sum = sizeof(uint)*8;//sizeof(static_bitsequence_rrr02_light);
  sum += uint_len(C_len,C_field_bits)*sizeof(uint);
  sum += O_len*sizeof(uint);
  sum += uint_len(C_sampling_len,C_sampling_field_bits)*sizeof(uint);
  sum += uint_len(O_pos_len,O_pos_field_bits)*sizeof(uint);
  //sum += E->size();
  return sum;
}

static_bitsequence_rrr02_light::~static_bitsequence_rrr02_light() {
  if(C!=NULL) delete [] C;
  if(O!=NULL) delete [] O;
  if(C_sampling!=NULL) delete [] C_sampling;
  if(O_pos!=NULL) delete [] O_pos;
  E = E->unuse();
}

int static_bitsequence_rrr02_light::save(FILE * fp) {
  uint C_len = len/BLOCK_SIZE_LIGHT + (len%BLOCK_SIZE_LIGHT!=0);
  uint C_field_bits = bits(BLOCK_SIZE_LIGHT);
  uint O_len = uint_len(1,O_bits_len);
  uint wr = RRR02_LIGHT_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  wr += fwrite(&len,sizeof(uint),1,fp);
  wr += fwrite(&ones,sizeof(uint),1,fp);
  wr += fwrite(&O_bits_len,sizeof(uint),1,fp);
  wr += fwrite(&sample_rate,sizeof(uint),1,fp);
  if(wr!=5) return -1;
  wr = fwrite(C,sizeof(uint),uint_len(C_len,C_field_bits),fp);
  if(wr!=uint_len(C_len,C_field_bits)) return -1;
  wr = fwrite(O,sizeof(uint),O_len,fp);
  if(wr!=O_len) return -1;
  return 0;
}

static_bitsequence_rrr02_light * static_bitsequence_rrr02_light::load(FILE * fp) {
  static_bitsequence_rrr02_light * ret = new static_bitsequence_rrr02_light();
  uint rd = 0, type;
  rd += fread(&type,sizeof(uint),1,fp);
  rd += fread(&ret->len,sizeof(uint),1,fp);
  rd += fread(&ret->ones,sizeof(uint),1,fp);
  rd += fread(&ret->O_bits_len,sizeof(uint),1,fp);
  rd += fread(&ret->sample_rate,sizeof(uint),1,fp);
  uint C_len = ret->len/BLOCK_SIZE_LIGHT + (ret->len%BLOCK_SIZE_LIGHT!=0);
  uint C_field_bits = bits(BLOCK_SIZE_LIGHT);
  uint O_len = uint_len(1,ret->O_bits_len);
  if(rd!=5 || type!=RRR02_LIGHT_HDR) {
    delete ret;
    return NULL;
  }
  ret->C = new uint[uint_len(C_len,C_field_bits)];
  rd = fread(ret->C,sizeof(uint),uint_len(C_len,C_field_bits),fp);
  if(rd!=uint_len(C_len,C_field_bits)) {
    ret->C=NULL;
    delete ret;
    return NULL;
  }
  ret->O = new uint[O_len];
  rd = fread(ret->O,sizeof(uint),O_len,fp);
  if(rd!=O_len) {
    ret->O=NULL;
    delete ret;
    return NULL;
  }
  ret->create_sampling(ret->sample_rate);
  return ret;
}
