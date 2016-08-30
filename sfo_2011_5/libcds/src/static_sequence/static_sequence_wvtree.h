/* static_sequence_wvtree.h
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
 
#ifndef STATIC_SEQUENCE_WVTREE_H
#define STATIC_SEQUENCE_WVTREE_H

#include <iostream>
#include <cassert>
#include <basics.h>
#include <static_bitsequence.h>
#include <static_bitsequence_builder.h>
#include <wt_node_internal.h>
#include <wt_coder_binary.h>
#include <alphabet_mapper.h>
#include <static_sequence.h>

using namespace std;

/** Wavelet tree implementation using pointers. 
 * 
 *  @author Francisco Claude
 */
class static_sequence_wvtree : public static_sequence {
  public:

    /** Builds a Wavelet Tree for the string
     * pointed by symbols assuming its length
     * equals n */
    static_sequence_wvtree(uint * symbols, uint n, wt_coder * coder, static_bitsequence_builder * bmb, alphabet_mapper * am);

    static_sequence_wvtree(uchar * symbols, uint n, wt_coder * coder, static_bitsequence_builder * bmb, alphabet_mapper * am);

    virtual ~static_sequence_wvtree();

    virtual uint rank(uint symbol, uint pos);
    virtual uint rankLessThan(uint &symbol, uint pos);

    virtual uint select(uint symbol, uint i);

    virtual uint access(uint pos);
    virtual uint access(uint pos, uint &rank)
    {
      return am->unmap(root->access(pos, rank));
    }
    
    // Returns all elements from interval [i, j] such that 
    // their value is in [min, max].
    virtual vector<int> access(uint i, uint j, uint min, uint max);
    virtual vector<int> accessAll(uint i, uint j);
    virtual uint count(uint i, uint j, uint min, uint max);

    virtual uint count(uint s);

    virtual uint size();
    
    virtual uint save(FILE * fp);
    static static_sequence_wvtree * load(FILE *fp);

  protected:

    static_sequence_wvtree();
    
    wt_node * root;
    wt_coder * c;
    alphabet_mapper * am;
    //bitmap_builder * bmb;
    
    /** Length of the string. */
    uint n;

    /** Height of the Wavelet Tree. */
    uint max_v;

    /** Flag for testing for correcteness. */
    bool test;

		
};
#endif /* _STATIC_SEQUENCE_WVTREE_H */
