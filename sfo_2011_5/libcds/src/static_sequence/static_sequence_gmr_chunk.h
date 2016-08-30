/* static_sequence_gmr_chunk.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * gmr_chunk
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

#ifndef _STATIC_SEQUENCE_GMR_CHUNK_H
#define _STATIC_SEQUENCE_GMR_CHUNK_H

#include <basics.h>
#include <static_sequence.h>
#include <static_bitsequence.h>
#include <static_bitsequence_builder.h>
#include <static_permutation.h>
#include <static_permutation_builder.h>
#include <cassert>
#include <iostream>

using namespace std;

/** Implementation of the Chunk of Golynski et al's rank/select
 * data structure [1].
 *
 * [1] A. Golynski and I. Munro and S. Rao. 
 * Rank/select operations on large alphabets: a tool for text indexing.
 * SODA 06.
 *
 * @author Francisco Claude
 */
class static_sequence_gmr_chunk: public static_sequence {
  public:
    /** Builds the structures needed for the chunk */
    static_sequence_gmr_chunk(uint * sequence, uint chunk_length, static_bitsequence_builder *bmb, static_permutation_builder *pmb);

    /** Destroy the chunk */
    ~static_sequence_gmr_chunk();

    virtual uint access(uint j);
    virtual uint select(uint i, uint j);
    virtual uint rank(uint i, uint j);
    virtual uint size();
    virtual uint save(FILE *fp);
    static static_sequence_gmr_chunk * load(FILE *fp);

  protected:
    /** Bitmap */
    static_bitsequence * X;
    /** Permutation */
    static_permutation * permutation;
    /** Size of the alphabet */
    uint sigma;
    /** Length of the chunk */
    //uint chunk_length;
    static_sequence_gmr_chunk();
};
#endif
