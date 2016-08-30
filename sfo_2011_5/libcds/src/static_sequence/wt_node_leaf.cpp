/* wt_node_leaf.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * wt_node_leaf
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

#include <wt_node_leaf.h>

wt_node_leaf::wt_node_leaf(uint symbol, uint count) {
	this->symbol = symbol;
	this->count = count;
}

wt_node_leaf::wt_node_leaf() {}

wt_node_leaf::~wt_node_leaf() {}

uint wt_node_leaf::rank(uint symbol, uint pos, uint l, wt_coder * c) {
	if(symbol!=this->symbol) return 0;
	pos++;
	return pos;
}

uint wt_node_leaf::rankLessThan(uint &symbol, uint pos) {
//    std::cout <<"this-symbol: " << (uchar)this->symbol << ", symbol = " << (uchar)symbol << ", pos = " << pos << std::endl;
    if (pos == (uint)-1 || symbol < this->symbol)
        return -1;
    symbol = this->symbol;
    pos++;
    return pos;
}

uint wt_node_leaf::select(uint symbol, uint pos, uint l, wt_coder * c) {
	if(symbol!=this->symbol) return (uint)-1;
	if(pos==0 || pos>count) return (uint)-1;
	return pos;
}

uint wt_node_leaf::access(uint pos) {
//   std::cout <<"this-symbol: " << (uchar)this->symbol << ", pos = " << pos << std::endl;

	return symbol;
}

uint wt_node_leaf::access(uint pos, uint &rank) {
    rank = pos+1;
    return symbol;
}

void wt_node_leaf::access(vector<int> &result, uint i, uint j, uint min, uint max, uint l, uint pivot)
{
//    std::cout << "At l = " << l << ", [" << i << ", " << j  << "], [" << min << ", " << max << "], symbol = " << symbol << std::endl;
    
    if (i <= j && symbol >= min && symbol <= max)
        result.push_back((int)symbol);
}

void wt_node_leaf::access(vector<int> &result, uint i, uint j)
{
//    std::cout << "At l = " << l << ", [" << i << ", " << j  << "], [" << min << ", " << max << "], symbol = " << symbol << std::endl;
    
    if (i <= j)
        result.push_back((int)symbol);
}

uint wt_node_leaf::access(uint i, uint j, uint min, uint max, uint l, uint pivot)
{
//    std::cout << "At l = " << l << ", [" << i << ", " << j  << "], [" << min << ", " << max << "], symbol = " << symbol << std::endl;
    
    if (i <= j && symbol >= min && symbol <= max)
        return 1;
    return 0;
}

uint wt_node_leaf::size() {
	return sizeof(wt_node_leaf);
}

uint wt_node_leaf::save(FILE *fp) {
  uint wr = WT_NODE_LEAF_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  wr += fwrite(&count,sizeof(uint),1,fp);
  wr += fwrite(&symbol,sizeof(uint),1,fp);
  return wr-3;
}

wt_node_leaf * wt_node_leaf::load(FILE *fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd!=WT_NODE_LEAF_HDR) return NULL;
  wt_node_leaf * ret = new wt_node_leaf();
  rd = fread(&(ret->count),sizeof(uint),1,fp);
  rd += fread(&(ret->symbol),sizeof(uint),1,fp);
  if(rd!=2) { delete ret; return NULL; }
  return ret;
}
