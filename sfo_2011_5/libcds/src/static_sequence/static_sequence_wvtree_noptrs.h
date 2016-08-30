/* static_sequence_wvtree_noptrs.h
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

#ifndef _STATIC_SEQUENCE_WVTREE_NOPTRS_H
#define _STATIC_SEQUENCE_WVTREE_NOPTRS_H

#include <iostream>
#include <cassert>
#include <basics.h>
#include <static_bitsequence.h>
#include <static_bitsequence_builder.h>
#include <static_sequence.h>
#include <alphabet_mapper.h>

using namespace std;

class static_sequence_wvtree_noptrs : public static_sequence {
  public:

    /** Builds a Wavelet Tree for the string
     * pointed by symbols assuming its length
     * equals n and uses bmb to build the bitsequence */
    static_sequence_wvtree_noptrs(uint * symbols, uint n, static_bitsequence_builder * bmb, alphabet_mapper * am, bool deleteSymbols = false);

    // symbols is an array of elements of "width" bits.
    static_sequence_wvtree_noptrs(uint * symbols, uint n, unsigned width, static_bitsequence_builder * bmb, alphabet_mapper * am, bool deleteSymbols = false);

    /** Destroys the Wavelet Tree */
    virtual ~static_sequence_wvtree_noptrs();

    virtual uint rank(uint symbol, uint pos);
    virtual uint select(uint symbol, uint i);
    virtual uint access(uint pos);
    virtual uint size();

    virtual vector<int> access(uint i, uint j, uint min, uint max);
    virtual vector<int> accessAll(uint i, uint j);
    virtual uint count(uint i, uint j, uint min, uint max);
    
    virtual uint save(FILE *fp);
    static static_sequence_wvtree_noptrs * load(FILE *fp);

  protected:
    void access(vector<int> &result, uint i, uint j, uint min, uint max, uint l, uint pivot, uint start, uint end);
    void accessAll(vector<int> &result, uint i, uint j, uint l, uint pivot, uint start, uint end);
    uint count(uint i, uint j, uint min, uint max, uint l, uint pivot, uint start, uint end);

    static_sequence_wvtree_noptrs();
    
    alphabet_mapper * am;
    /** Only one bit-string for the Wavelet Tree. */
    static_bitsequence **bitstring, *occ;

    /** Length of the string. */
    uint n;

    /** Height of the Wavelet Tree. */
    uint height,max_v;

    /** Obtains the maximum value from the string
     * symbols of length n */
    uint max_value(uint * symbols, uint n);
    uint max_value(uint * symbols, unsigned width, uint n);

    /** How many bits are needed to represent val */
    uint bits(uint val);

    /** Returns true if val has its ind-th bit set
     * to one. */
    bool is_set(uint val, uint ind);

    /** Sets the ind-th bit in val */
    uint set(uint val, uint ind);

    /** Recursive function for building the Wavelet Tree. */
    void build_level(uint **bm, uint *symbols, uint level, uint length, uint offset);
    void build_level(uint **bm, uint *symbols, unsigned width, uint level, uint length, uint offset);
};
#endif
