/* static_sequence_builder_gmr_chunk.h
 * Copyright (C) 2008, Francisco Claude, all rights reserved.
 *
 * Sequence builder gmr chunk
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

#ifndef _STATIC_SEQUENCE_BUILDER_GMR_CHUNK_H
#define _STATIC_SEQUENCE_BUILDER_GMR_CHUNK_H

#include <basics.h>
#include <static_sequence_gmr_chunk.h>
#include <static_sequence_builder.h>
#include <static_permutation_builder.h>

/** gmr chunk builder
 *  @author Francisco Claude
 */
class static_sequence_builder_gmr_chunk : public static_sequence_builder {
  public:
    static_sequence_builder_gmr_chunk(static_bitsequence_builder *bmb, static_permutation_builder *pmb);
    virtual ~static_sequence_builder_gmr_chunk() {}
    virtual static_sequence * build(uint * seq, uint len);
    
  protected:
    static_bitsequence_builder *bmb;
    static_permutation_builder *pmb;
};

#endif
