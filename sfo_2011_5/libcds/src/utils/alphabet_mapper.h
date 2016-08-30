/* alphabet_mapper.h
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

#ifndef _ALPHABET_MAPPER_H
#define _ALPHABET_MAPPER_H

#include <basics.h>
#include <iostream>

#define ALPHABET_MAPPER_NONE_HDR 2
#define ALPHABET_MAPPER_CONT_HDR 3

using namespace std;

/** Base class for alphabet mappers
 * 
 *  @author Francisco Claude
 */
class alphabet_mapper {
  public:
		alphabet_mapper();
    virtual ~alphabet_mapper() {}
    /** Maps the symbol */
    virtual uint map(uint s)=0;
    /** Unmaps the symbol */
    virtual uint unmap(uint s)=0;
    /** Returns the size of the mapper */
    virtual uint size()=0;
    /** Saves the mapper to a file */
    virtual uint save(FILE *fp)=0;
    /** Loads the mapper from a file */
    static alphabet_mapper * load(FILE * fp);
		virtual void use();
		virtual void unuse();
	protected:
		uint user_count;
};

#include <alphabet_mapper_none.h>
#include <alphabet_mapper_cont.h>

#endif /* _ALPHABET_MAPPER_H */
