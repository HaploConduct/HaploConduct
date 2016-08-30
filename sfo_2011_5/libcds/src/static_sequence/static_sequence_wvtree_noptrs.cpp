/* static_sequence_wvtree_noptrs.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_sequence_wvtree_noptrs definition
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
 
#include <static_sequence_wvtree_noptrs.h>

static_sequence_wvtree_noptrs::static_sequence_wvtree_noptrs(uint * symbols, uint n, static_bitsequence_builder * bmb, alphabet_mapper * am, bool deleteSymbols) {
  this->n=n;
  this->am=am;
	am->use();
  for(uint i=0;i<n;i++)
    symbols[i] = am->map(symbols[i]);
  max_v=max_value(symbols,n);
  height=bits(max_v);
  uint *occurrences=new uint[max_v+1];
  for(uint i=0;i<=max_v;i++) occurrences[i]=0;
  for(uint i=0;i<n;i++)
    occurrences[symbols[i]]++;
  uint to_add=0;
  for(uint i=0;i<max_v;i++)
    if(occurrences[i]==0) to_add++;
  uint * new_symb = new uint[n+to_add];
  for(uint i=0;i<n;i++)
    new_symb[i] = symbols[i];

  if (deleteSymbols)
  {
      delete [] symbols;
      symbols = 0;
  }

  to_add = 0;
  for(uint i=0;i<max_v;i++)
    if(occurrences[i]==0) {
      occurrences[i]++;
      new_symb[n+to_add]=i;
      to_add++;
    }
  uint new_n = n+to_add;
  for(uint i=1;i<=max_v;i++)
    occurrences[i] += occurrences[i-1];
  uint *oc = new uint[(new_n+1)/W+1];
  for(uint i=0;i<(new_n+1)/W+1;i++)
    oc[i] = 0;
  for(uint i=0;i<=max_v;i++)
    bitset(oc,occurrences[i]-1);
  bitset(oc,new_n);
  occ = bmb->build(oc,new_n+1);
  delete [] occurrences;
  this->n = new_n;
  uint ** _bm=new uint*[height];
  for(uint i=0;i<height;i++) {
    _bm[i] = new uint[new_n/W+1];
    for(uint j=0;j<new_n/W+1;j++)
      _bm[i][j]=0;
  }
  build_level(_bm,new_symb,0,new_n,0);
  bitstring = new static_bitsequence*[height];
  for(uint i=0;i<height;i++) {
  	bitstring[i] = bmb->build(_bm[i],new_n);
    delete [] _bm[i];
  }
  delete [] _bm;
  
  if (!deleteSymbols)
      for(uint i=0;i<n;i++)
          symbols[i] = am->unmap(symbols[i]);

// delete [] new_symb; // already deleted in build_level()!
  delete [] oc;
}

// symbols is an array of elements of "width" bits
static_sequence_wvtree_noptrs::static_sequence_wvtree_noptrs(uint * symbols, uint n, unsigned width, static_bitsequence_builder * bmb, alphabet_mapper * am, bool deleteSymbols) {
  this->n=n;
  this->am=am;
	am->use();
  for(uint i=0;i<n;i++)
      set_field(symbols, width, i, am->map(get_field(symbols, width, i)));
  max_v=max_value(symbols, width, n);
  height=bits(max_v);
  uint *occurrences=new uint[max_v+1];
  for(uint i=0;i<=max_v;i++) occurrences[i]=0;
  for(uint i=0;i<n;i++)
      occurrences[get_field(symbols, width, i)]++;
  uint to_add=0;
  for(uint i=0;i<max_v;i++)
    if(occurrences[i]==0) to_add++;
  uint * new_symb = new uint[((n+to_add)*width)/W + 1];
  for(uint i=0;i<n;i++)
      set_field(new_symb, width, i, get_field(symbols, width, i));

  if (deleteSymbols)
  {
      delete [] symbols;
      symbols = 0;
  }

  to_add = 0;
  for(uint i=0;i<max_v;i++)
    if(occurrences[i]==0) {
      occurrences[i]++;
      set_field(new_symb, width, n+to_add, i);
      to_add++;
    }
  uint new_n = n+to_add;
  for(uint i=1;i<=max_v;i++)
    occurrences[i] += occurrences[i-1];
  uint *oc = new uint[(new_n+1)/W+1];
  for(uint i=0;i<(new_n+1)/W+1;i++)
    oc[i] = 0;
  for(uint i=0;i<=max_v;i++)
    bitset(oc,occurrences[i]-1);
  bitset(oc,new_n);
  occ = bmb->build(oc,new_n+1);
  delete [] occurrences;
  this->n = new_n;
  uint ** _bm=new uint*[height];
  for(uint i=0;i<height;i++) {
    _bm[i] = new uint[new_n/W+1];
    for(uint j=0;j<new_n/W+1;j++)
      _bm[i][j]=0;
  }
  build_level(_bm,new_symb,width,0,new_n,0);
  bitstring = new static_bitsequence*[height];
  for(uint i=0;i<height;i++) {
  	bitstring[i] = bmb->build(_bm[i],new_n);
    delete [] _bm[i];
  }
  delete [] _bm;
  
  if (!deleteSymbols)
      for(uint i=0;i<n;i++)
          set_field(symbols, width, i, am->unmap(get_field(symbols, width, i)));

// delete [] new_symb; // already deleted in build_level()!
  delete [] oc;
}

static_sequence_wvtree_noptrs::static_sequence_wvtree_noptrs() {
}

static_sequence_wvtree_noptrs::~static_sequence_wvtree_noptrs() {
  for(uint i=0;i<height;i++)
    delete bitstring[i];
  delete [] bitstring;
  delete occ;
 	am->unuse();
}

uint static_sequence_wvtree_noptrs::save(FILE *fp) {
  uint wr = WVTREE_NOPTRS_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  wr += fwrite(&n,sizeof(uint),1,fp);
  wr += fwrite(&max_v,sizeof(uint),1,fp);
  wr += fwrite(&height,sizeof(uint),1,fp);
  if(wr!=4) return 1;
  if(am->save(fp)) return 1;
  for(uint i=0;i<height;i++)
    if(bitstring[i]->save(fp)) return 1;
	if(occ->save(fp)) return 1;
  return 0;
}

static_sequence_wvtree_noptrs * static_sequence_wvtree_noptrs::load(FILE *fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd!=WVTREE_NOPTRS_HDR) return NULL;
  static_sequence_wvtree_noptrs * ret = new static_sequence_wvtree_noptrs();
  rd = fread(&ret->n,sizeof(uint),1,fp);
  rd += fread(&ret->max_v,sizeof(uint),1,fp);
  rd += fread(&ret->height,sizeof(uint),1,fp);
  if(rd!=3) {
    delete ret;
    return NULL;
  }
  ret->am = alphabet_mapper::load(fp);
  if(ret->am==NULL) {
    delete ret;
    return NULL;
  }
	ret->am->use();
  ret->bitstring = new static_bitsequence*[ret->height];
  for(uint i=0;i<ret->height;i++) {
    ret->bitstring[i] = static_bitsequence::load(fp);
    if(ret->bitstring[i]==NULL){
      delete ret;
      return NULL;
    }
  }
	ret->occ = static_bitsequence::load(fp);
	if(ret->occ==NULL) {
		delete ret;
		return NULL;
	}
  return ret;
}

uint static_sequence_wvtree_noptrs::access(uint pos) {
  uint level=0;
  uint ret=0;
  uint start=0;
  uint end=n-1;
  while(level<height) {
    assert(pos>=start && pos<=end);
    if(bitstring[level]->access(pos)) {
      ret=set(ret,level);
      pos=bitstring[level]->rank1(pos-1)-bitstring[level]->rank1(start-1);
      start=(bitstring[level]->rank1(end)-bitstring[level]->rank1(start-1));
      start=end-start+1;
      pos+=start;
    }
    else {
      pos=pos-start-(bitstring[level]->rank1(pos)-bitstring[level]->rank1(start-1));
      end=end-start-(bitstring[level]->rank1(end)-bitstring[level]->rank1(start-1));
      end+=start;
      pos+=start;
    }
    level++;
  }
  return am->unmap(ret);
}

uint static_sequence_wvtree_noptrs::rank(uint symbol, uint pos) {
  symbol = am->map(symbol);
  uint level=0;
  uint start=0;
  uint end=n-1;
  uint count=0;
  while(level<height) {
    if(is_set(symbol,level)) {
      pos=bitstring[level]->rank1(pos)-bitstring[level]->rank1(start-1)-1;
      count=pos+1;
      start=(bitstring[level]->rank1(end)-bitstring[level]->rank1(start-1));
      start=end-start+1;
      pos+=start;
    }
    else {
      pos=pos-start+bitstring[level]->rank1(start-1)-bitstring[level]->rank1(pos);
      count=pos+1;
      end=end-start-(bitstring[level]->rank1(end)-bitstring[level]->rank1(start-1));
      end+=start;
      pos+=start;
    }
    level++;
    if(count==0) return 0;
  }
  return count;
}

vector<int> static_sequence_wvtree_noptrs::access(uint i, uint j, uint min, uint max)
{
    vector<int> resultSet;
//    cout << "height = " << height << endl;
    access(resultSet, i, j, am->map(min), am->map(max), 0, 0, 0, n-1);
    return resultSet;
}

void static_sequence_wvtree_noptrs::access(vector<int> &result, uint i, uint j, uint min, uint max, uint l, uint pivot, uint start, uint end)
{
    uint symbol = pivot | (1 << (height-l-1));
    //std::cout << "At l = " << l << ", [" << i << ", " << j  << "], [" << min << ", " << max << "], [" << start << ", " << end << "], symbol = " << symbol << std::endl;

    if (l == height)
    {
        if (i <= j && pivot >= min && pivot <= max && start <= end)
            result.push_back(am->unmap((int)pivot));
        return;
    }

    if (j < i || max < min || end < start)
        return;

    if (min < symbol)
    {
        // Recurse left
        uint newi = i + bitstring[l]->rank1(start-1) - bitstring[l]->rank1(i-1);
        uint newend = end - (bitstring[l]->rank1(end) - bitstring[l]->rank1(start-1));
        uint newj = j + bitstring[l]->rank1(start-1) - bitstring[l]->rank1(j) + 1;

        uint newmax = max < symbol - 1 ? max : symbol - 1;
        if (newj > start)
            access(result, newi, newj-1, min, newmax, l+1, pivot, start, newend);
    }

    if (max >= symbol)
    {
        // Recurse right
        uint newstart = (bitstring[l]->rank1(end)-bitstring[l]->rank1(start-1));
        newstart = end - newstart + 1;
        uint newi = bitstring[l]->rank1(i-1)-bitstring[l]->rank1(start-1) + newstart;
        uint newj = bitstring[l]->rank1(j)-bitstring[l]->rank1(start-1) + newstart;

        uint newmin = min > symbol ? min : symbol;
        if (newj > newstart)
            access(result, newi, newj-1, newmin, max, l+1, symbol, newstart, end);
    }
}


vector<int> static_sequence_wvtree_noptrs::accessAll(uint i, uint j)
{
    vector<int> resultSet;
    if (j < i)
        return resultSet;

    resultSet.reserve(j-i+1);
    accessAll(resultSet, i, j, 0, 0, 0, n-1);
    return resultSet;
}

void static_sequence_wvtree_noptrs::accessAll(vector<int> &result, uint i, uint j, uint l, uint pivot, uint start, uint end)
{
    uint symbol = pivot | (1 << (height-l-1));
//    std::cout << "At l = " << l << ", [" << i << ", " << j  << "], [" << start << ", " << end << "], symbol = " << symbol << std::endl;

    if (l == height)
    {
        if (i <= j && start <= end)
            result.push_back(am->unmap((int)pivot));
        return;
    }

    if (j < i || end < start)
        return;

    {
        // Recurse left
        uint newi = i + bitstring[l]->rank1(start-1) - bitstring[l]->rank1(i-1);
        uint newend = end - (bitstring[l]->rank1(end) - bitstring[l]->rank1(start-1));
        uint newj = j + bitstring[l]->rank1(start-1) - bitstring[l]->rank1(j) + 1;

        if (newj > start)
            accessAll(result, newi, newj-1, l+1, pivot, start, newend);
    }

    {
        // Recurse right
        uint newstart = (bitstring[l]->rank1(end)-bitstring[l]->rank1(start-1));
        newstart = end - newstart + 1;
        uint newi = bitstring[l]->rank1(i-1)-bitstring[l]->rank1(start-1) + newstart;
        uint newj = bitstring[l]->rank1(j)-bitstring[l]->rank1(start-1) + newstart;

        if (newj > newstart)
            accessAll(result, newi, newj-1, l+1, symbol, newstart, end);
    }
}


uint static_sequence_wvtree_noptrs::count(uint i, uint j, uint min, uint max)
{
    return count(i, j, am->map(min), am->map(max), 0, 0, 0, n-1);
}

uint static_sequence_wvtree_noptrs::count(uint i, uint j, uint min, uint max, uint l, uint pivot, uint start, uint end)
{
    uint symbol = pivot | (1 << (height-l-1));
    //std::cout << "At l = " << l << ", [" << i << ", " << j  << "], [" << min << ", " << max << "], [" << start << ", " << end << "], symbol = " << symbol << std::endl;

    if (l == height)
    {
        if (i <= j && pivot >= min && pivot <= max && start <= end)
            return 1;
        return 0;
    }

    if (j < i || max < min || end < start)
        return 0;

    uint result = 0;
    if (min < symbol)
    {
        // Recurse left
        uint newi = i + bitstring[l]->rank1(start-1) - bitstring[l]->rank1(i-1);
        uint newend = end - (bitstring[l]->rank1(end) - bitstring[l]->rank1(start-1));
        uint newj = j + bitstring[l]->rank1(start-1) - bitstring[l]->rank1(j) + 1;

        uint newmax = max < symbol - 1 ? max : symbol - 1;
        if (newj > start)
            result += count(newi, newj-1, min, newmax, l+1, pivot, start, newend);
    }

    if (max >= symbol)
    {
        // Recurse right
        uint newstart = (bitstring[l]->rank1(end)-bitstring[l]->rank1(start-1));
        newstart = end - newstart + 1;
        uint newi = bitstring[l]->rank1(i-1)-bitstring[l]->rank1(start-1) + newstart;
        uint newj = bitstring[l]->rank1(j)-bitstring[l]->rank1(start-1) + newstart;

        uint newmin = min > symbol ? min : symbol;
        if (newj > newstart)
            result += count(newi, newj-1, newmin, max, l+1, symbol, newstart, end);
    }
    return result;
}



inline uint get_start(uint symbol, uint mask) {
  return symbol&mask;
}

inline uint get_end(uint symbol, uint mask) {
  return get_start(symbol,mask)+!mask+1;
}

uint static_sequence_wvtree_noptrs::select(uint symbol, uint j) {
  symbol = am->map(symbol);
  uint mask = (1<<height)-2;
  uint sum=2;
  uint level = height-1;
  uint pos=j;
  while(true) {
    uint start = get_start(symbol,mask);
    uint end = min(max_v+1,start+sum);
    start = (start==0)?0:(occ->select1(start)+1);
    end = occ->select1(end+1)-1;
    if(is_set(symbol,level)) {
      uint ones_start = bitstring[level]->rank1(start-1);
      pos = bitstring[level]->select1(ones_start+pos)-start+1;
    }
    else {
      uint ones_start = bitstring[level]->rank1(start-1);
      pos = bitstring[level]->select0(start-ones_start+pos)-start+1;
    }
    mask <<=1;
    sum <<=1;
    if(level==0) break;
    level--;
  }
  return pos-1;
}

uint static_sequence_wvtree_noptrs::size() {
  uint ptrs = sizeof(static_sequence_wvtree_noptrs)+height*sizeof(static_sequence*);
  uint bytesBitstrings = 0;
  for(uint i=0;i<height;i++)
    bytesBitstrings += bitstring[i]->size();
  return bytesBitstrings+occ->size()+ptrs;
}

void static_sequence_wvtree_noptrs::build_level(uint **bm, uint *symbols, uint level, uint length, uint offset) {
  if(level==height)
  {
      delete [] symbols;
      return;
  }
  uint cleft=0;
  for(uint i=0;i<length;i++)
    if(!is_set(symbols[i],level))
      cleft++;
  uint cright=length-cleft;
  uint *left=new uint[cleft], *right=new uint[cright];
  cleft=cright=0;
  for(uint i=0;i<length;i++)
  if(!is_set(symbols[i],level)) {
    left[cleft++]=symbols[i];
    bitclean(bm[level],offset+i);
  }
  else {
    right[cright++]=symbols[i];
    bitset(bm[level],offset+i);
  }
  
  delete [] symbols;
  symbols = 0;
  
  build_level(bm,left,level+1,cleft,offset);
  left = 0; // Gets deleted in recursion.
  build_level(bm,right,level+1,cright,offset+cleft);
  right = 0; // Gets deleted in recursion.
  //delete [] left;
  //delete [] right;
}

// symbols is an array of elements of "width" bits.
void static_sequence_wvtree_noptrs::build_level(uint **bm, uint *symbols, unsigned width, uint level, uint length, uint offset) {
    if(level==height)
    {
        delete [] symbols;
        return;
    }
    uint cleft=0;
    for(uint i=0;i<length;i++)
        if(!is_set(get_field(symbols, width, i),level))
            cleft++;
    uint cright=length-cleft;
    uint *left=new uint[(cleft*width)/W + 1], 
        *right=new uint[(cright*width)/W + 1];
    cleft=cright=0;
    for(uint i=0;i<length;i++)
        if(!is_set(get_field(symbols,width,i),level)) {
            set_field(left,width,cleft++,get_field(symbols, width,i));
            bitclean(bm[level],offset+i);
        }
        else {
            set_field(right,width,cright++,get_field(symbols,width,i));
            bitset(bm[level],offset+i);
        }
  
    delete [] symbols;
    symbols = 0;
  
    build_level(bm,left,width,level+1,cleft,offset);
    left = 0; // Gets deleted in recursion.
    build_level(bm,right,width,level+1,cright,offset+cleft);
    right = 0; // Gets deleted in recursion.
    //delete [] left;
    //delete [] right;
}

uint static_sequence_wvtree_noptrs::max_value(uint *symbols, uint n) {
  uint max_v = 0;
  for(uint i=0;i<n;i++)
    max_v = max(symbols[i],max_v);
  return max_v;
}

uint static_sequence_wvtree_noptrs::max_value(uint *symbols, unsigned width, uint n) {
  uint max_v = 0;
  for(uint i=0;i<n;i++)
      max_v = max(get_field(symbols, width, i),max_v);
  return max_v;
}

uint static_sequence_wvtree_noptrs::bits(uint val) {
  uint ret = 0;
  while(val!=0) {
    ret++;
    val >>= 1;
  }
  return ret;
}

bool static_sequence_wvtree_noptrs::is_set(uint val, uint ind) {
  assert(ind<height);
  return (val & (1<<(height-ind-1)))!=0;
}


uint static_sequence_wvtree_noptrs::set(uint val, uint ind) {
  assert(ind<=height);
  return val | (1<<(height-ind-1));
}
