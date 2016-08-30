/*****************************************************************************
 * Copyright (C) 2005, Rodrigo Gonzalez, all rights reserved.                *
 *                                                                           *
 * New RANK, SELECT, SELECT-NEXT and SPARSE RANK implementations.            *
 *                                                                           *
 * This library is free software; you can redistribute it and/or             *
 * modify it under the terms of the GNU Lesser General Public                *
 * License as published by the Free Software Foundation; either              *
 * version 2.1 of the License, or (at your option) any later version.        *
 *                                                                           *
 * This library is distributed in the hope that it will be useful,           * 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         *
 * Lesser General Public License for more details.                           *
 *                                                                           *
 * You should have received a copy of the GNU Lesser General Public          *
 * License along with this library; if not, write to the Free Software       *
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA *
 ****************************************************************************/

#ifndef _BITSTREAM_H_
#define _BITSTREAM_H_

#include "BlockArray.h"
#include "Tools.h"

class BitRank {
private:
    // Check word length
#if W == 32
    static const unsigned wordShift = 5;
    static const unsigned superFactor = 8; // 256 bit blocks
#else
    static const unsigned wordShift = 6;
    static const unsigned superFactor = 4; // 256 bit blocks
#endif
    static const unsigned s = W*superFactor;
    
    ulong *data; //here is the bit-array
    bool owner;
    ulong n,integers;
    BlockArray *Rs; //superblock array
    uchar *Rb; //block array
    ulong BuildRankSub(ulong,  ulong); //internal use of BuildRank
    void BuildRank(); //crea indice para rank
public:
    BitRank(ulong *, ulong, bool);
    BitRank(FILE *);
    ~BitRank(); //destructor    
    void Save(FILE *);
    ulong rank(ulong i); //Rank from 0 to n-1
    ulong select(ulong x); // gives the position of the x:th 1.
    ulong select0(ulong x); // gives the position of the x:th 0.

    bool IsBitSet(ulong i);
    ulong NumberOfBits();
};

#endif
