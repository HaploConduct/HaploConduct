/* huffman_codes.h
   Copyright (C) 2008, Francisco Claude, all rights reserved.

   Wrapper for huff written by Gonzalo Navarro

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#ifndef HUFFMAN_CODES_H
#define HUFFMAN_CODES_H

#include <basics.h>
#include <huff.h>

/** Wrapper for the canonical huffman implementation of Gonzalo Navarro. 
 * 
 *  @author Francisco Claude
 */
class huffman_codes {

  public:
    /** Creates the codes for the sequence seq of length n */
    huffman_codes(uint * seq, uint n);
    huffman_codes(uchar * seq, uint n);
    ~huffman_codes();
    
    /** Encodes symb into stream at bit-position pos, return the ending position (bits) */
    ulong encode(uint symb, uint * stream, ulong pos);
    
    /** decodes into symb from stream at bit-position pos, returns the new position */
    ulong decode(uint * symb, uint * stream, ulong pos);
    
    /** Returns the maximum length of a code */
    uint max_length();
    
    /** Returns the size of the table */
    uint size();
    
    /** Saves the coder to a file */
    uint save(FILE *fp);
    
    /** Loads a coder from a file */
    static huffman_codes * load(FILE *fp);
    
  protected:
    huffman_codes();
    THuff huff_table;
};

#endif
