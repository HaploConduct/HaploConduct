/* wt_coder_huff.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * wt_coder_huff definition
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
 
#ifndef wt_coder_huff_h
#define wt_coder_huff_h

#include <basics.h>
#include <wt_coder.h>
#include <huffman_codes.h>
#include <alphabet_mapper.h>

/** Uses huffman codes to determine the shape of the wavelet tree  
 * 
 *  @author Francisco Claude
 */
class wt_coder_huff: public wt_coder {
	public:
    /** Buils a wt_coder_huff using the sequence of length n and the alphabet_mapper
     *  to determine the huffman codes */
		wt_coder_huff(uint *symbs, uint n, alphabet_mapper * am);
		wt_coder_huff(uchar *symbs, uint n, alphabet_mapper * am);
		virtual ~wt_coder_huff();
		virtual bool is_set(uint symbol, uint l);
		virtual bool done(uint symbol, uint l);
    virtual uint size();
    virtual uint save(FILE *fp);
    static wt_coder_huff * load(FILE *fp);
    //uint * get_buffer(uint symbol, uint *n);

	protected:
    wt_coder_huff();
		huffman_codes * hc;
    uint * buffer;
    uint last_symbol, s_len;
};

#endif
