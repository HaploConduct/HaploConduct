/* static_bitsequence.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * static_bitsequence definition
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

#ifndef _STATIC_BITSEQUENCE_H
#define	_STATIC_BITSEQUENCE_H

#define RRR02_HDR 2
#define BRW32_HDR 3
#define RRR02_LIGHT_HDR 4
#define SDARRAY_HDR 5

#include <basics.h>
#include <iostream>


using namespace std;

/** Base class for static bitsequences, contains many abstract functions, so this can't
 *  be instantiated. It includes base implementations for rank0, select0 and select1 based
 *  on rank0.
 * 
 *  @author Francisco Claude
 */
class static_bitsequence {
  
public:
  virtual ~static_bitsequence() {};

	/** Returns the number of zeros until position i */
  virtual uint rank0(uint i);

	/** Returns the position of the i-th zero 
	 * @return (uint)-1 if i=0, len if i>num_zeros or the position */
  virtual uint select0(uint i);

	/** Returns the number of ones until position i */
  virtual uint rank1(uint i);

	/** Returns the position of the i-th one 
	 * @return (uint)-1 if i=0, len if i>num_ones or the position */
  virtual uint select1(uint i);

	virtual uint select_next1(uint i);
	virtual uint select_next0(uint i);

	/** Returns the i-th bit */
  virtual bool access(uint i);

	/** Returns the length in bits of the bitmap */
  virtual uint length();

	/** Returns how many ones are in the bitstring */
  virtual uint count_one();

	/** Returns how many zeros are in the bitstring */
  virtual uint count_zero();

	/** Returns the size of the structure in bytes */
  virtual uint size()=0;

  /** Stores the bitmap given a file pointer, return 0 in case of success */
	virtual int save(FILE * fp)=0;
  
  /** Reads a bitmap determining the type */
  static static_bitsequence * load(FILE * fp);
  
protected:
	/** Length of the bitstring */
  uint len;
	/** Number of ones in the bitstring */
	uint ones;
  
};

#include <static_bitsequence_rrr02.h>
#include <static_bitsequence_rrr02_light.h>
#include <static_bitsequence_naive.h>
#include <static_bitsequence_brw32.h>
#include <static_bitsequence_sdarray.h>

#endif	/* _STATIC_BITSEQUENCE_H */
