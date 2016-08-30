/* static_bitsequence_rrr02.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * RRR02 Bitsequence - 
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


#ifndef _STATIC_BITSEQUENCE_RRR02_H
#define	_STATIC_BITSEQUENCE_RRR02_H

#define BLOCK_SIZE 15
#define DEFAULT_SAMPLING 32

#include <static_bitsequence.h>
#include <table_offset.h>
#include <cassert>
#include <iostream>

using namespace std;

/** Implementation of Raman, Raman and Rao's [1] proposal for rank/select capable
 *  data structures, it achieves space nH_0, O(sample_rate) time for rank and O(log len) 
 *  for select. The practial implementation is based on [2]
 * 
 *  [1] R. Raman, V. Raman and S. Rao. Succinct indexable dictionaries with applications 
 *     to encoding $k$-ary trees and multisets. SODA02.
 *  [2] F. Claude and G. Navarro. Practical Rank/Select over Arbitrary Sequences. SPIRE08.
 *
 *  @author Francisco Claude
 */
class static_bitsequence_rrr02: public static_bitsequence {
public:
  static_bitsequence_rrr02(uint * bitseq, uint len, uint sample_rate=DEFAULT_SAMPLING);
  virtual ~static_bitsequence_rrr02();
  
  /** Returns the number of zeros until position i */
  virtual uint rank0(uint i);
  
  /** Returns the number of ones until position i */
  virtual uint rank1(uint i);
  
  /** Returns the position of the i-th zero 
   * @return (uint)-1 if i=0, len if i>num_zeros or the position */
  virtual uint select0(uint i);
  
  /** Returns the position of the i-th one 
   * @return (uint)-1 if i=0, len if i>num_ones or the position */
  virtual uint select1(uint i);
  
  /** Returns the i-th bit */
  virtual bool access(uint i);
  
  /** Returns the size of the structure in bytes */
  virtual uint size();
  
  /** Stores the bitmap given a file pointer, return 0 in case of success */
	virtual int save(FILE * fp);
  
  /** Reads the bitmap from a file pointer, returns NULL in case of error */
	static static_bitsequence_rrr02 * load(FILE * fp);
  
  /** Creates a new sampling for the queries */
	void create_sampling(uint sampling_rate);
  
  /** Frees the space required by the table E, which is static and global
   *  to all instances.
   */
  static void delete_E() {
    delete E;
  }

  
protected:
  static_bitsequence_rrr02();
	/** Classes and offsets */
  uint *C, *O;
	/** Length of C and O (in uints) */
	uint C_len, O_len;
	/** Bits required per field for C and in total for O */
	uint C_field_bits, O_bits_len;
	/** C and O samplings */
	uint *C_sampling, *O_pos;
	/** Length of the samplings */
	uint C_sampling_len,O_pos_len;
	/** Lenght in bits per field */
	uint C_sampling_field_bits,O_pos_field_bits;
	/** Sample rate */
	uint sample_rate;

	static table_offset * E;
};

#endif	/* _STATIC_BITSEQUENCE_RRR02_H */
