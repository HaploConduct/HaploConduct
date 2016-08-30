/* static_sequence.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_sequence definition
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

#ifndef _STATIC_SEQUENCE_H
#define _STATIC_SEQUENCE_H


#include <basics.h>
#include <iostream>
#include <vector>

#define WVTREE_HDR 2
#define GMR_CHUNK_HDR 3
#define GMR_HDR 4
#define WVTREE_NOPTRS_HDR 5
#define BS_HDR 6

using namespace std;

/** Base class for static sequences, contains many abstract functions, so this can't
 *  be instantiated. 
 * 
 *  @author Francisco Claude
 */
class static_sequence {
  
public:
  static_sequence();
  virtual ~static_sequence();

  /** Returns the number of occurrences of c until position i */
  virtual uint rank(uint c, uint i)=0;
  virtual uint rankLessThan(uint &i, uint j)
  {
      //assert(0); // Implemented only in static_sequence_wvtree
      return -1;
  }
  
  /** Returns the position of the i-th c 
   * @return (uint)-1 if i=0, len if i exceeds the number of cs */
  virtual uint select(uint c, uint i)=0;
	virtual uint select_next(uint c, uint i);

  /** Returns the i-th element */
  virtual uint access(uint i)=0;
  virtual uint access(uint i, uint &rank)
  {
      //assert(0); // Implemented only in static_sequence_wvtree
      return -1;
  }

  // Returns all elements from interval [i, j] such that 
  // their value is in [min, max].
  virtual vector<int> access(uint i, uint j, uint min, uint max)
  {
      //assert(0); // Implemented only in static_sequence_wvtree
      return vector<int>();
  }

  // Returns all elements from interval [i, j] 
  virtual vector<int> accessAll(uint i, uint j)
  {
      //assert(0); // Implemented only in static_sequence_wvtree
      return vector<int>();
  }

  // Counts the number of elements in interval [i,j] such that
  // their values are in [min,max]
  virtual uint count(uint i, uint j, uint min, uint max)
  {
      //assert(0); // Implemented only in static_sequence_wvtree
      return 0;
  }

  /** Returns the length of the sequence */
  virtual uint length();

  /** Returns how many cs are in the sequence */
  virtual uint count(uint c);

  /** Returns the size of the structure in bytes */
  virtual uint size()=0;

  /** Stores the bitmap given a file pointer, return 0 in case of success */
  virtual uint save(FILE * fp)=0;

	virtual bool test(uint * seq, uint n);
  
  /** Reads a bitmap determining the type */
  static static_sequence * load(FILE * fp);
  
protected:
  /** Length of the bitstring */
  uint len;
  
};

#include <static_sequence_wvtree.h>
#include <static_sequence_gmr_chunk.h>
#include <static_sequence_wvtree_noptrs.h>
#include <static_sequence_gmr.h>
#include <static_sequence_bs.h>

#endif  /* _STATIC_SEQUENCE_H */
