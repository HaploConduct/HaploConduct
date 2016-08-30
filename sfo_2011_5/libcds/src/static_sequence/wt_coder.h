/* wt_coder.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * wt_coder definition
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
 
#ifndef wt_coder_h
#define wt_coder_h

#include <basics.h>
#include <iostream>

using namespace std;

#define WT_CODER_HUFF_HDR 2
#define WT_CODER_BINARY_HDR 3

/** Coder that defines the shape of a wavelet tree 
 * 
 *  @author Francisco Claude
 */
class wt_coder {
	public:
		wt_coder();
		virtual void use();
		virtual void unuse();
    virtual ~wt_coder() {}; 
    /** Tells if at level l the symbol is represented by a one or a zero */
		virtual bool is_set(uint symbol, uint l)=0;
    /** Tells if the path of symbol becomes unique at level l */
		virtual bool done(uint symbol, uint l)=0;
    /** Returns the size of the coder */
    virtual uint size()=0;
    /** Returns the depth of the tree */
    virtual uint depth() {
        return -1; // Implemented in wt_coder_binary
    }
    /** Saves the coder to a file, returns 0 in case of success */
    virtual uint save(FILE *fp)=0;
    /** Loads a coder from a file, returns NULL in case of error */
    static wt_coder * load(FILE *fp);
	protected:
		uint user_count;
};

#include <wt_coder_huff.h>
#include <wt_coder_binary.h>

#endif
