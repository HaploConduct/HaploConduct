/* static_sequence_builder_gmr_chunk.cpp
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

#include <static_sequence_builder_gmr_chunk.h>

static_sequence_builder_gmr_chunk::static_sequence_builder_gmr_chunk(static_bitsequence_builder *bmb, static_permutation_builder *pmb) {
  this->bmb = bmb;
  this->pmb = pmb;
}

static_sequence * static_sequence_builder_gmr_chunk::build(uint * seq, uint len) {
  return new static_sequence_gmr_chunk(seq,len,bmb,pmb);
}
