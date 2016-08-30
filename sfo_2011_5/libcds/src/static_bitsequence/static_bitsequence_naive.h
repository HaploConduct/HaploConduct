/* static_bitsequence_naive.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * Naive Bitsequence - don't use, only for testing
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


#ifndef _STATIC_BITSEQUENCE_NAIVE_H
#define	_STATIC_BITSEQUENCE_NAIVE_H

#include <static_bitsequence.h>

/** Class used for testing, should not be used with long bitmaps
 *  @author Francisco Claude
 */
class static_bitsequence_naive: public static_bitsequence {
public:
  /** Builds a naive bitsequence, receives the bitmap and the length
   *  in bits
   */
  static_bitsequence_naive(uint * bitseq, uint len);
  
  virtual ~static_bitsequence_naive();
  
  /** Returns the number of ones until position i */
  virtual uint rank1(uint i);
  
  /** Returns the position of the i-th one 
   * @return (uint)-1 if i=0, len if i>num_ones or the position */
  virtual uint select1(uint i);
  
  /** Returns the i-th bit */
  virtual bool access(uint i);
  
  /** Returns the size of the structure in bytes */
  virtual uint size();
  
  /** - Not implemented - */
	virtual int save(FILE * fp);
  
  
protected:
  uint * bitseq;
};

#endif	/* _STATIC_BITSEQUENCE_NAIVE_H */

