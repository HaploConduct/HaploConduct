/* huff.h
   Copyright (C) 2008, Gonzalo Navarro, all rights reserved.

   Canonical Huffman

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

#ifndef HUFFINCLUDED
#define HUFFINCLUDED

#include <basics.h>

typedef struct
   { uint max,lim;   // maximum symbol (0..max), same excluding zero freqs
     uint depth; // max symbol length
     union
       { uint *spos; // symbol positions after sorting by decr freq (enc)
	 uint *symb; // symbols sorted by freq (dec)
       } s;
     uint *num;  // first pos of each length (dec), number of each length (enc)
     uint *fst;  // first code (numeric) of each length (dec)
     ulong total; // total length to achieve, in bits
   } THuff;


/** Creates Huffman encoder given symbols 0..lim with frequencies 
 *  freq[i], ready for compression 
 * 
 *  @author Gonzalo Navarro
 */
THuff createHuff (uint *freq, uint lim);

/** Encodes symb using H, over stream[ptr...lim] (ptr and lim are
 *  bit positions of stream). Returns the new ptr. 
 * 
 *  @author Gonzalo Navarro
 */
ulong encodeHuff (THuff H, uint symb, uint *stream, ulong ptr);

/** Decodes *symb using H, over stream[ptr...lim] (ptr and lim are
 *  bit positions of stream). Returns the new ptr. 
 * 
 *  @author Gonzalo Navarro
 */
ulong decodeHuff (THuff H, uint *symb, uint *stream, ulong ptr);

/** Writes H in file f 
 * 
 *  @author Gonzalo Navarro
 */
void saveHuff (THuff H, FILE *f);

/** Size of H written on file 
 * 
 *  @author Gonzalo Navarro
 */
uint sizeHuff (THuff H);

/** Frees H 
 * 
 *  @author Gonzalo Navarro
 */	
void freeHuff (THuff H);

/** Loads H from file f, prepared for encoding or decoding depending
 *  on enc 
 * 
 *  @author Gonzalo Navarro
 */
THuff loadHuff (FILE *f, int enc);

#endif
