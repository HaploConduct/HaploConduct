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

#ifndef _STATIC_SEQUENCE_BS_H
#define _STATIC_SEQUENCE_BS_H


#include <basics.h>
#include <static_sequence.h>
#include <static_bitsequence.h>

/** static_sequence represented using one bitmap per symbol, doesn't support efficient access
 * 
 *  @author Francisco Claude
 */
class static_sequence_bs : public static_sequence {
  
public:
  static_sequence_bs(uint * seq, uint n, alphabet_mapper * am, static_bitsequence_builder * bmb);
  virtual ~static_sequence_bs();

  virtual uint rank(uint c, uint i);
  
  virtual uint select(uint c, uint i);
	virtual uint select_next(uint c, uint i);

  virtual uint access(uint i);

  virtual uint size();

  virtual uint save(FILE * fp);

  /** Reads a bitmap determining the type */
  static static_sequence_bs * load(FILE * fp);
  
protected:
  uint sigma;
	static_bitsequence ** bitmaps;
	alphabet_mapper * am;

	static_sequence_bs();
  
};


#endif  /* _STATIC_SEQUENCE_BS_H */

