/* static_sequence_wvtree.cpp
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_sequence_wvtree definition
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
 
#include <static_sequence_wvtree.h>

static_sequence_wvtree::static_sequence_wvtree(uint * symbols, uint n, wt_coder * c, static_bitsequence_builder * bmb, alphabet_mapper * am) {
    this->n = n;
  for(uint i=0;i<n;i++) 
    symbols[i] = am->map(symbols[i]);
  this->am = am;
	am->use();
  this->c=c;
	c->use();
	root = new wt_node_internal(symbols, n, 0, c, bmb);
  for(uint i=0;i<n;i++) 
    symbols[i] = am->unmap(symbols[i]);  
}

static_sequence_wvtree::static_sequence_wvtree(uchar * symbols, uint n, wt_coder * c, static_bitsequence_builder * bmb, alphabet_mapper * am) {
    this->n = n;
  for(uint i=0;i<n;i++) 
    symbols[i] = (uchar)am->map((uint)symbols[i]);
  this->am = am;
	am->use();
  this->c=c;
	c->use();
        uint *done = new uint[n/W+1];
        for (uint i = 0; i < n/W+1; i++)
            done[i] = 0;
	root = new wt_node_internal(symbols, n, 0, c, bmb, 0, done);
        delete [] done;
        delete [] symbols;
        symbols = 0; // Already deleted!
//  for(uint i=0;i<n;i++) 
//    symbols[i] = (uchar)am->unmap((uint)symbols[i]);  
}

static_sequence_wvtree::static_sequence_wvtree() {}

static_sequence_wvtree::~static_sequence_wvtree() {
	delete root;
	am->unuse();
  c->unuse(); 
}

uint static_sequence_wvtree::rank(uint symbol, uint pos) {
	return root->rank(am->map(symbol), pos, 0, c);
}

uint static_sequence_wvtree::rankLessThan(uint &symbol, uint pos) {
    uint s = am->map(symbol);
//    std::cout << "lessthan..." << std::endl;
    uint r = root->rankLessThan(s, pos);
    symbol = am->unmap(s);
    return r;
}


uint static_sequence_wvtree::count(uint s) {
  return root->rank(am->map(s), len-1, 0, c);
}

uint static_sequence_wvtree::select(uint symbol, uint pos) {
	uint ret = root->select(am->map(symbol), pos, 0, c);
	if(ret==((uint)-1)) return (uint)-1;
	return ret-1;
}

uint static_sequence_wvtree::access(uint pos) {
	return am->unmap(root->access(pos));
}

vector<int> static_sequence_wvtree::access(uint i, uint j, uint min, uint max)
{
    vector<int> resultSet;
    root->access(resultSet, i, j, am->map(min), am->map(max), c->depth()-1, 0);
    for (vector<int>::iterator it = resultSet.begin(); it != resultSet.end(); ++it)
        *it = am->unmap(*it);
    return resultSet;
}

vector<int> static_sequence_wvtree::accessAll(uint i, uint j)
{
    vector<int> resultSet;
    if (j < i)
        return resultSet;

    // resultSet.reserve(j-i+1); // avoid reallocation
    root->access(resultSet, i, j);
    for (vector<int>::iterator it = resultSet.begin(); it != resultSet.end(); ++it)
        *it = am->unmap(*it);
    return resultSet;
}

uint static_sequence_wvtree::count(uint i, uint j, uint min, uint max)
{
    return root->access(i, j, am->map(min), am->map(max), c->depth()-1, 0);
}


uint static_sequence_wvtree::size() {
	/*cout << "WT: " << root->size() << endl;
	cout << "Coder: " << c->size() << endl;
	cout << "AM: " << am->size() << endl;*/
	return sizeof(static_sequence_wvtree)+sizeof(uint)+root->size()+am->size()+c->size();
}

uint static_sequence_wvtree::save(FILE * fp) { 
  uint wr = WVTREE_HDR;
  wr = fwrite(&wr,sizeof(uint),1,fp);
  if(wr!=1) return 1;
  wr = fwrite(&n,sizeof(uint),1,fp);
  if(wr!=1) return 1;
  if(c->save(fp)) return 1;
  if(am->save(fp)) return 1;
  if(root->save(fp)) return 1;
  return 0;
}

static_sequence_wvtree * static_sequence_wvtree::load(FILE *fp) {
  uint rd;
  if(fread(&rd,sizeof(uint),1,fp)!=1) return NULL;
  if(rd!=WVTREE_HDR) return NULL;
  static_sequence_wvtree * ret = new static_sequence_wvtree();
  if(fread(&ret->n,sizeof(uint),1,fp)!=1) return NULL;
  ret->c = wt_coder::load(fp);
	ret->c->use();
  ret->am = alphabet_mapper::load(fp);
	ret->am->use();
  ret->root = wt_node::load(fp);
  return ret;
}
